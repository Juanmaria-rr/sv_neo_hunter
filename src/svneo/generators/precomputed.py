"""
precomputed.py — stage-2 backend: reuse an existing generator output.

Two uses, both routine:

1. **Reproduce a run without the tool.** Stage 2 is the slow, dependency-heavy
   step; everything after it is fast and ours. Re-running the analysis with new
   thresholds should not require regenerating peptides, and must not risk a
   different tool version silently changing them.
2. **Bring in peptides from elsewhere** — another pipeline, a collaborator, an
   older run — as long as the tables satisfy the contract in `base.py`.

The backend performs no generation and no filtering: it loads, validates against
the contract, and records where the files came from so the provenance survives
into the run manifest.
"""
from __future__ import annotations

import os

import pandas as pd

from .base import GenerationResult, GeneratorError, PeptideGenerator


class PrecomputedGenerator(PeptideGenerator):
    """Load `<prefix>.all_neopeptides.txt` + `<prefix>.anno.txt` from disk."""

    name = "precomputed"

    def __init__(self, source_prefix: str | None = None):
        #: Path prefix of the stored output. May also be passed per call.
        self.source_prefix = source_prefix

    @staticmethod
    def _find(prefix: str) -> tuple[str, str]:
        """Locate the pair, tolerating either naming convention."""
        peptide = next((p for p in (f"{prefix}.all_neopeptides.txt",
                                    f"{prefix}_all_neopeptides.txt")
                        if os.path.exists(p)), None)
        anno = next((p for p in (f"{prefix}.anno.txt", f"{prefix}_anno.txt")
                     if os.path.exists(p)), None)
        if not peptide or not anno:
            raise GeneratorError(
                f"no stored generator output for prefix '{prefix}'. Expected "
                f"'{prefix}.all_neopeptides.txt' and '{prefix}.anno.txt'.")
        return peptide, anno

    def generate(self, admitted_vcf: str, out_dir: str, prefix: str,
                 source_prefix: str | None = None, **_) -> GenerationResult:
        source = source_prefix or self.source_prefix
        if not source:
            raise GeneratorError(
                "the `precomputed` backend needs a source prefix: set "
                "`resources.peptides_from.<sample>` in the config.")

        peptide_path, anno_path = self._find(source)
        peptides = pd.read_csv(peptide_path, sep="\t", dtype=str)
        annotation = pd.read_csv(anno_path, sep="\t", dtype=str)
        self.validate_output(peptides, annotation, source=peptide_path)
        for note in self.warn_optional(peptides):
            print(f"    NOTE: {note}")

        # The admitted VCF is NOT re-applied to a stored table: the stored output
        # was generated from whatever call set it was generated from. Recorded so
        # a mismatch is visible rather than assumed away.
        return GenerationResult(
            peptides=peptides, annotation=annotation, generator=self.name,
            version=f"loaded from {os.path.basename(source)}",
            stats={"loaded_from": source,
                   "note": "peptides were NOT regenerated; the admitted VCF of "
                           "this run was not applied to them"},
            outputs={"peptides": peptide_path, "annotation": anno_path})
