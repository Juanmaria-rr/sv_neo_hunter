"""
neosv_trace.py — stage-2 backend: NeoSV-Trace.

    ┌──────────────────────────────────────────────────────────────────┐
    │  THIS FILE IS THE ONLY PLACE THE PIPELINE TOUCHES NeoSV-Trace.   │
    │  The vendored copy lives in vendor/neosv_trace/ (see             │
    │  vendor/VENDOR.md for provenance, licence and patch policy).     │
    │  No other module imports it. Replacing the generator means       │
    │  adding a sibling of this file — nothing else changes.           │
    └──────────────────────────────────────────────────────────────────┘

WHAT NeoSV-Trace DOES
---------------------
Annotates each breakend against Ensembl transcripts, reconstructs the junction
nucleotide sequence, checks frame, translates, and emits every 8-11mer spanning
the junction via a sliding window. Peptides are generated PER TRANSCRIPT, which
is what makes isoform-level expression meaningful at stage 7.

WHY THE CLI IS NOT USED
-----------------------
Its `main()` offers two modes, neither of which fits a recurrence test:

    --anno-only   writes <prefix>.anno.txt and stops — no peptides at all
    (default)     writes <prefix>.all_neopeptides.txt and then immediately calls
                  netMHCpan/MHCflurry, which needs a licensed binary and an HLA
                  file the identity test deliberately does not use

The peptides are written *before* the MHC call, so running the CLI and letting it
fail would "work" — while leaving a non-zero exit code and empty `.net.in` /
`.net.out` files that read as a failed run to anyone who opens the directory
later. This backend calls the same library functions in the same order as
`main()` and stops after `write_all_neopeptides`.

It imports the tool rather than reimplementing it, deliberately: the reference
catalogue was built with this code, so an identity test against that catalogue is
only valid with the same peptide-generation logic.

THREE UPSTREAM CONSTRAINTS HANDLED HERE
---------------------------------------
1. Only `.vcf` / `.bedpe` are accepted, checked by file EXTENSION — a `.vcf.gz`
   is rejected outright, so a gzipped input is decompressed to a temporary file.
2. The `prefix` column of both outputs is derived internally as
   `basename(file).split('.')[0]`, so that temporary file must be named after the
   sample or every row carries a meaningless label.
3. The tool applies NO admission logic of its own: handed a raw VCF it loads PON
   and INFERRED records too (observed: 176 SVs from a 350-record file whose PASS
   subset is 108). It must be given the stage-1 admitted VCF.
"""
from __future__ import annotations

import gzip
import os
import shutil
import sys

import pandas as pd

from .. import criteria
from .base import GenerationResult, GeneratorError, PeptideGenerator

#: The vendored copy, resolved relative to this file: <repo>/vendor.
#: Tried before any configured path and before an installed package, so a run is
#: reproducible from the repository alone.
VENDOR_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "vendor"))


class NeoSVTraceGenerator(PeptideGenerator):
    """NeoSV-Trace peptide generation, stopping before the MHC step."""

    name = "neosv_trace"

    def __init__(self, neosv_path: str | None = None):
        #: Search order: explicit config path -> vendored copy -> installed package.
        self.search_paths = [p for p in (neosv_path, VENDOR_DIR) if p]

    # -- import ------------------------------------------------------------

    def _import(self):
        for path in self.search_paths:
            if os.path.isdir(os.path.join(path, "neosv_trace")) and path not in sys.path:
                sys.path.insert(0, path)
        try:
            import neosv_trace
            from neosv_trace.input import vcf_load, bedpe_load, ensembl_load, get_window_range
            from neosv_trace.sv_utils import (sv_pattern_infer_vcf,
                                              sv_pattern_infer_bedpe, remove_duplicate)
            from neosv_trace.annotation_utils import sv_to_sveffect
            from neosv_trace.fusion_utils import sv_to_svfusion
            from neosv_trace.sequence_utils import set_nt_seq, set_aa_seq, generate_neoepitopes
            from neosv_trace.output import write_annot, write_all_neopeptides
        except ImportError as error:
            raise GeneratorError(
                f"NeoSV-Trace is not importable ({error}).\n"
                f"Searched: {self.search_paths}\n"
                f"The repository ships a vendored copy at {VENDOR_DIR}/neosv_trace; "
                f"if it is missing, restore it (see vendor/VENDOR.md) or set "
                f"`resources.neosv_path` in the config.") from error
        return dict(
            module=neosv_trace, vcf_load=vcf_load, bedpe_load=bedpe_load,
            ensembl_load=ensembl_load, get_window_range=get_window_range,
            sv_pattern_infer_vcf=sv_pattern_infer_vcf,
            sv_pattern_infer_bedpe=sv_pattern_infer_bedpe,
            remove_duplicate=remove_duplicate, sv_to_sveffect=sv_to_sveffect,
            sv_to_svfusion=sv_to_svfusion, set_nt_seq=set_nt_seq,
            set_aa_seq=set_aa_seq, generate_neoepitopes=generate_neoepitopes,
            write_annot=write_annot, write_all_neopeptides=write_all_neopeptides)

    def _version(self) -> str:
        """Vendored commit, read from VENDOR.md so the manifest carries it."""
        vendor_doc = os.path.join(VENDOR_DIR, "VENDOR.md")
        if os.path.exists(vendor_doc):
            with open(vendor_doc) as handle:
                for line in handle:
                    if line.startswith("| Commit "):
                        return line.split("`")[1] if "`" in line else "unknown"
        return "unknown"

    # -- input preparation -------------------------------------------------

    @staticmethod
    def _plain_vcf(vcf_path: str, sample: str, work_dir: str) -> tuple[str, bool]:
        """Decompress if needed, naming the copy after the sample (constraint 2)."""
        if not vcf_path.endswith(".gz"):
            return vcf_path, False
        os.makedirs(work_dir, exist_ok=True)
        plain = os.path.join(work_dir, f"{sample}.vcf")
        with gzip.open(vcf_path, "rt") as src, open(plain, "w") as dst:
            shutil.copyfileobj(src, dst)
        return plain, True

    # -- the stage ---------------------------------------------------------

    def generate(self, admitted_vcf: str, out_dir: str, prefix: str,
                 release: int = 115, cache_dir: str | None = None,
                 complete_transcript: bool = False, keep_temp: bool = False,
                 **_) -> GenerationResult:
        api = self._import()
        os.makedirs(out_dir, exist_ok=True)
        work_dir = os.path.join(out_dir, "_work")
        vcf_file, is_temp = self._plain_vcf(admitted_vcf, prefix, work_dir)

        try:
            ensembl = api["ensembl_load"](str(release), None, None, cache_dir)

            if vcf_file.endswith(".vcf"):
                svs = [api["sv_pattern_infer_vcf"](v) for v in api["vcf_load"](vcf_file)]
            elif vcf_file.endswith(".bedpe"):
                svs = [api["sv_pattern_infer_bedpe"](b) for b in api["bedpe_load"](vcf_file)]
            else:
                raise GeneratorError(f"input must be .vcf or .bedpe, got {vcf_file}")

            n_loaded = len(svs)
            svs = api["remove_duplicate"](svs)
            row_prefix = os.path.basename(vcf_file).split(".")[0]

            # --- SV annotation ---
            effects = [unit for sv in svs
                       for unit in api["sv_to_sveffect"](sv, ensembl, complete_transcript)]
            anno_path = os.path.join(out_dir, f"{prefix}.anno.txt")
            api["write_annot"](anno_path, effects, prefix=row_prefix)

            # --- peptide generation ---
            window = api["get_window_range"](
                f"{min(criteria.PEPTIDE_LENGTHS)}-{max(criteria.PEPTIDE_LENGTHS)}")
            fusions = [f for f in (api["sv_to_svfusion"](sv, ensembl) for sv in svs)
                       if not f.is_empty()]
            for fusion in fusions:
                fusion.nt_sequence = api["set_nt_seq"](fusion)
                fusion.aa_sequence = api["set_aa_seq"](fusion)
                fusion.neoepitopes = api["generate_neoepitopes"](fusion, window)
                fusion.neoepitope_meta = {}

            peptide_path = os.path.join(out_dir, f"{prefix}.all_neopeptides.txt")
            api["write_all_neopeptides"](peptide_path, fusions, prefix=row_prefix)

            # STOP. main() would now run MHC prediction; the recurrence test is
            # HLA-independent by design (criteria.RUN_MHC_PREDICTION).
            peptides = pd.read_csv(peptide_path, sep="\t", dtype=str)
            annotation = pd.read_csv(anno_path, sep="\t", dtype=str)
            self.validate_output(peptides, annotation, source="NeoSV-Trace")
            for note in self.warn_optional(peptides):
                print(f"    NOTE: {note}")

            return GenerationResult(
                peptides=peptides, annotation=annotation, generator=self.name,
                version=self._version(),
                stats={"svs_loaded": n_loaded, "svs_after_dedup": len(svs),
                       "fusions_non_empty": len(fusions),
                       "mhc_prediction": "skipped by design"},
                outputs={"peptides": peptide_path, "annotation": anno_path})
        finally:
            if is_temp and not keep_temp and os.path.exists(vcf_file):
                os.remove(vcf_file)
