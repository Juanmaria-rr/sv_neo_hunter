"""
neosv.py — stage-2 backend: NeoSV (published, MIT) plus this repository's own
candidate-peptide output.

    ┌──────────────────────────────────────────────────────────────────────┐
    │  Peptide sequences come from NeoSV, unmodified except for a single   │
    │  documented frame patch (vendor/patches/001-…). Everything this       │
    │  pipeline additionally needs — SV identifiers, junction offsets, and  │
    │  the pre-MHC candidate table — is in _neosv_extensions.py and is our  │
    │  code, not the tool's.                                                │
    └──────────────────────────────────────────────────────────────────────┘

NeoSV:  Shi, Y., Jing, B. & Xi, R. Comprehensive analysis of neoantigens derived
        from structural variation across whole genomes from 2528 tumors.
        Genome Biology 24, 169 (2023). MIT licence — see vendor/LICENSE.NeoSV.

WHY THE VENDORED COPY CARRIES A PATCH
-------------------------------------
pyensembl changed its coding-sequence accessor at v2.3.13: the stop codon is no
longer excluded from `coding_sequence_position_ranges`. NeoSV compensates for the
old behaviour by subtracting 3, which with any later pyensembl shifts the 3' side
of a fusion one amino acid out of frame. Measured on one sample, the unpatched
tool emits 170 peptides against 64 for the patched one, with **none in common** —
the frame-shifted sequences do not match the wild-type protein and survive the
`set(mut) − set(wt)` filter as spurious neopeptides.

The patch removes the subtraction. It is three characters, documented in
`vendor/patches/001-pyensembl-stop-codon-frame.md`, and it is the only
modification to upstream code in this repository.

Credit: the frame issue was identified in NeoSV-Trace
(github.com/winterga/NeoSV-Trace) by Greyson Wintergerst.

WHY THE CLI IS NOT USED
-----------------------
NeoSV writes only post-netMHCpan binders. A recurrence test asks whether a
peptide *exists* in a sample, which is HLA-independent; filtering by predicted
binding would discard true matches and would require a licensed predictor for a
question that does not need one. This backend calls the library functions and
writes the candidate table itself.
"""
from __future__ import annotations

import gzip
import os
import shutil
import sys

import pandas as pd

from .. import criteria
from . import _neosv_extensions as ext
from .base import GenerationResult, GeneratorError, PeptideGenerator

#: <repo>/vendor, resolved relative to this file so a clone is self-contained.
VENDOR_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..", "vendor"))


class NeoSVGenerator(PeptideGenerator):
    """Published NeoSV + this repository's candidate-peptide output."""

    name = "neosv"
    #: Package directory name inside vendor/, and the module it exposes.
    package = "neosv"

    def __init__(self, neosv_path: str | None = None, **_):
        self.search_paths = [p for p in (neosv_path, VENDOR_DIR) if p]

    def _import(self):
        for path in self.search_paths:
            if os.path.isdir(os.path.join(path, self.package)) and path not in sys.path:
                sys.path.insert(0, path)
        try:
            import importlib
            module = importlib.import_module(self.package)
            return {name: importlib.import_module(f"{self.package}.{name}")
                    for name in ("input", "sv_utils", "fusion_utils",
                                 "sequence_utils")} | {"root": module}
        except ImportError as error:
            raise GeneratorError(
                f"NeoSV is not importable ({error}). Searched: {self.search_paths}. "
                f"The repository ships a vendored copy at {VENDOR_DIR}/{self.package}; "
                f"restore it (see vendor/VENDOR.md) or set `resources.neosv_path`."
            ) from error

    @staticmethod
    def _plain_vcf(vcf_path: str, sample: str, work_dir: str) -> tuple[str, bool]:
        """NeoSV accepts only .vcf/.bedpe, checked by extension."""
        if not vcf_path.endswith(".gz"):
            return vcf_path, False
        os.makedirs(work_dir, exist_ok=True)
        plain = os.path.join(work_dir, f"{sample}.vcf")
        with gzip.open(vcf_path, "rt") as src, open(plain, "w") as dst:
            shutil.copyfileobj(src, dst)
        return plain, True

    def _version(self) -> str:
        patches = os.path.join(VENDOR_DIR, "patches")
        applied = sorted(f for f in os.listdir(patches)
                         if f.endswith(".patch")) if os.path.isdir(patches) else []
        return f"NeoSV (MIT) + {len(applied)} patch(es): {', '.join(applied) or 'none'}"

    def generate(self, admitted_vcf: str, out_dir: str, prefix: str,
                 release: int = 115, cache_dir: str | None = None,
                 complete_transcript: bool = False, keep_temp: bool = False,
                 **_) -> GenerationResult:
        api = self._import()
        os.makedirs(out_dir, exist_ok=True)
        work_dir = os.path.join(out_dir, "_work")
        vcf_file, is_temp = self._plain_vcf(admitted_vcf, prefix, work_dir)

        try:
            ensembl = api["input"].ensembl_load(str(release), None, None, cache_dir)

            if vcf_file.endswith(".vcf"):
                svs = [api["sv_utils"].sv_pattern_infer_vcf(v)
                       for v in api["input"].vcf_load(vcf_file)]
            elif vcf_file.endswith(".bedpe"):
                svs = [api["sv_utils"].sv_pattern_infer_bedpe(b)
                       for b in api["input"].bedpe_load(vcf_file)]
            else:
                raise GeneratorError(f"input must be .vcf or .bedpe, got {vcf_file}")

            n_loaded = len(svs)
            svs = api["sv_utils"].remove_duplicate(svs)

            # SV identifiers are recovered from the VCF rather than threaded
            # through the library, keeping the upstream diff to the frame patch.
            sv_id_map = ext.build_sv_id_map(admitted_vcf)

            window = api["input"].get_window_range(
                f"{min(criteria.PEPTIDE_LENGTHS)}-{max(criteria.PEPTIDE_LENGTHS)}")
            sequence_utils = api["sequence_utils"]

            fusions = []
            for sv in svs:
                fusion = api["fusion_utils"].sv_to_svfusion(sv, ensembl)
                if fusion.is_empty():
                    continue
                fusion.sv = sv
                fusion.nt_sequence = sequence_utils.set_nt_seq(fusion)
                fusion.aa_sequence = sequence_utils.set_aa_seq(fusion)
                # Peptides come from upstream, untouched.
                fusion.neoepitopes = sequence_utils.generate_neoepitopes(fusion, window)
                fusions.append(fusion)

            peptide_path = os.path.join(out_dir, f"{prefix}.all_neopeptides.txt")
            anno_path = os.path.join(out_dir, f"{prefix}.anno.txt")
            n_rows = ext.write_all_neopeptides(peptide_path, fusions, sv_id_map, prefix)
            ext.write_annotation(anno_path, fusions, sv_id_map, prefix)

            peptides = pd.read_csv(peptide_path, sep="\t", dtype=str)
            annotation = pd.read_csv(anno_path, sep="\t", dtype=str)
            self.validate_output(peptides, annotation, source="NeoSV + extensions")
            for note in self.warn_optional(peptides):
                print(f"    NOTE: {note}")

            unmatched = int((peptides["sv_id"] == "NA").sum()) if n_rows else 0
            if unmatched:
                print(f"    NOTE: {unmatched} peptide rows could not be matched to a "
                      f"VCF ID; they will not deduplicate to an event")

            return GenerationResult(
                peptides=peptides, annotation=annotation, generator=self.name,
                version=self._version(),
                stats={"svs_loaded": n_loaded, "svs_after_dedup": len(svs),
                       "fusions_non_empty": len(fusions),
                       "peptide_rows": n_rows,
                       "rows_without_sv_id": unmatched,
                       "mhc_prediction": "not part of this backend"},
                outputs={"peptides": peptide_path, "annotation": anno_path})
        finally:
            if is_temp and not keep_temp and os.path.exists(vcf_file):
                os.remove(vcf_file)
