"""
base.py — the contract every peptide generator must satisfy (stage 2).

WHY THIS IS A PLUGGABLE STAGE
-----------------------------
Stage 2 is the only stage that depends on an external tool, and it is the one
most likely to be replaced: by a newer NeoSV-Trace, by a different SV-to-peptide
tool, or by a previously computed output being reused. Every other stage —
admission, the cross, QC, event dedup, confidence, RNA evidence, the null model —
is ours and is generator-agnostic.

Isolating it behind this contract means:
  * the rest of the pipeline never imports the tool, so swapping it touches one
    directory;
  * a run can be reproduced from a stored output without the tool installed;
  * the boundary with third-party code is visible in the directory tree
    (`generators/` + `vendor/`), not buried in a function call.

THE CONTRACT
------------
A generator receives an ADMITTED VCF (stage 1 output — never a raw call set) and
returns two tables:

`peptides`  one row per (SV x transcript x window position). Required columns:
    sv_id            the sample's own SV identifier -> the event dedup key
    neopeptide       the peptide string             (criteria.PEPTIDE_COLUMN)
    spans_junction   bool; does the peptide cross the breakpoint?
    gene1            gene at breakend 1             -> gene concordance
    transcript_id1   transcript at breakend 1       -> isoform-level expression
  Recommended: chrom1, pos1, chrom2, pos2, gene2, transcript_id2, svtype,
    frame_effect, pep_length.

`annotation`  one row per annotated SV. Required: sv_id, chrom1, pos1, gene1,
    transcript_id1, chrom2, pos2, gene2, transcript_id2, svtype.

A generator MUST NOT apply MHC-binding prediction, expression filters, or any
admission logic of its own: those are separate, later, and configurable stages.
Generating peptides is the only job.

`validate_output()` enforces the required columns, because a silently missing
column has already cost this project a real error — reading the wrong peptide
column collapsed a count from 20,710 to 4.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

#: Columns without which downstream stages cannot run.
REQUIRED_PEPTIDE_COLUMNS = ("sv_id", "neopeptide", "gene1", "transcript_id1")
REQUIRED_ANNOTATION_COLUMNS = ("sv_id", "chrom1", "pos1", "gene1", "chrom2", "pos2")

#: Optional but used when present; their absence degrades a feature, not the run.
OPTIONAL_PEPTIDE_COLUMNS = ("spans_junction", "svtype", "frame_effect", "pep_length",
                            "gene2", "transcript_id2", "chrom1", "pos1",
                            "chrom2", "pos2")


class GeneratorError(RuntimeError):
    """Raised when peptides can be neither generated nor loaded."""


@dataclass
class GenerationResult:
    """What a generator returns, plus the provenance to reproduce it."""
    peptides: pd.DataFrame
    annotation: pd.DataFrame
    generator: str                       # backend name, e.g. "neosv_trace"
    version: str = "unknown"             # tool version / vendored commit
    stats: dict = field(default_factory=dict)   # counts the backend reports
    outputs: dict = field(default_factory=dict) # paths it wrote

    def manifest(self) -> dict:
        """Provenance block, echoed into the run's summary.json so a peptide set
        can always be traced to the generator and version that produced it."""
        return {"generator": self.generator, "version": self.version,
                "stats": self.stats, "outputs": self.outputs,
                "n_peptide_rows": len(self.peptides),
                "n_unique_peptides": int(self.peptides["neopeptide"].nunique())
                                     if len(self.peptides) else 0}


class PeptideGenerator:
    """Interface. Subclass, set `name`, implement `generate()`."""

    name = "base"

    def generate(self, admitted_vcf: str, out_dir: str, prefix: str,
                 **kwargs) -> GenerationResult:
        raise NotImplementedError

    # -- shared helpers ---------------------------------------------------

    @staticmethod
    def validate_output(peptides: pd.DataFrame, annotation: pd.DataFrame,
                        source: str = "") -> None:
        """Fail loudly on a missing required column.

        A generator that silently omits a column produces a run that looks
        successful and is wrong — the failure mode this check exists to prevent.
        """
        for frame, required, what in ((peptides, REQUIRED_PEPTIDE_COLUMNS, "peptides"),
                                      (annotation, REQUIRED_ANNOTATION_COLUMNS,
                                       "annotation")):
            missing = [c for c in required if c not in frame.columns]
            if missing:
                raise GeneratorError(
                    f"{what} table from {source or 'the generator'} is missing "
                    f"required column(s) {missing}; it has {list(frame.columns)}. "
                    f"See generators/base.py for the contract.")

    @staticmethod
    def warn_optional(peptides: pd.DataFrame) -> list[str]:
        """Report absent optional columns and what is lost, rather than silently
        degrading. Missing `spans_junction`, for instance, disables the filter
        that removes peptides not crossing the breakpoint."""
        notes = []
        if "spans_junction" not in peptides.columns:
            notes.append("no `spans_junction` column: cannot drop peptides that "
                         "do not cross the breakpoint (criteria.REQUIRE_SPANS_JUNCTION)")
        if "svtype" not in peptides.columns:
            notes.append("no `svtype` column: SV-type concordance not evaluated")
        return notes
