"""
config.py — the run definition, and the generalisation boundary of this repo.

WHAT MAKES THE ANALYSIS GENERAL
-------------------------------
The original question was "do reference-cohort SV-neoantigens recur in the sample set?".
Nothing in the method is specific to that cohort or to the parental line. Two things are:

    a REFERENCE  — a catalogue of peptides to look for (any cohort, any caller)
    a set of SAMPLES — anything with an SV VCF, optionally RNA and DNA

Everything else (admission, peptide generation, the cross, QC, event dedup,
confidence, RNA evidence, the null, attribution) is identical whatever those two
are. So they, and only they, live in the config; the code never names a cohort
or a cell line.

SAMPLES FORM A LINEAGE, NOT A LIST
----------------------------------
`parent:` makes the sample set a tree. It drives two things:
  1. ORDER — a parent is always analysed before its children, so a hit in the
     parent is attributable to the parent rather than to whatever was done to
     derive the child.
  2. ATTRIBUTION — the earliest sample carrying a breakpoint explains it.
For unrelated samples (a patient cohort, say) leave `parent: null` everywhere:
the lineage collapses to a flat list and attribution is skipped.

VCF KIND IS A PROPERTY OF THE SAMPLE, NOT OF THE FILE NAME
----------------------------------------------------------
`vcf_kind: germline` means "this VCF describes the background genome shared with
descendants" and gets the population-frequency filters; `somatic` means "acquired
relative to the parent" and is assumed already panel-filtered by the caller.
Getting this wrong silently re-analyses a parental genome as if it were acquired
— an observed failure mode where every pair's germline VCF was the *same*
parental genome (86.2% of breakends at identical positions).
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import yaml


@dataclass
class Sample:
    """One analysable unit: an SV call set plus whatever modalities exist."""
    name: str
    vcf: str
    vcf_kind: str = "somatic"          # germline | somatic
    parent: str | None = None          # lineage; None = root / unrelated
    rna_bam: str | None = None
    isofox_dir: str | None = None
    isofox_prefix: str | None = None
    dna_bam: str | None = None
    notes: str = ""

    @property
    def has_rna(self) -> bool:
        return bool(self.rna_bam and os.path.exists(self.rna_bam))

    @property
    def has_expression(self) -> bool:
        return bool(self.isofox_dir and os.path.isdir(self.isofox_dir))

    def modalities(self) -> dict:
        """What can actually be run for this sample. A missing modality yields
        `NA` downstream — never `fail`. Missing data is not negative evidence."""
        return {"sv": os.path.exists(self.vcf), "rna_bam": self.has_rna,
                "expression": self.has_expression,
                "dna_bam": bool(self.dna_bam and os.path.exists(self.dna_bam))}


@dataclass
class Reference:
    """The query set: peptides to look for, and what to match them on."""
    peptides: str                       # TSV/CSV, one row per reference peptide
    peptide_column: str = "neoantigen"
    gene_column: str | None = "gene1"   # None disables gene concordance
    svtype_column: str | None = "svtype"
    pool: str | None = None             # optional larger set WITH coordinates
    pool_peptide_column: str = "neoantigen"
    proteome: str | None = None         # for the self-peptide test
    name: str = "reference"


@dataclass
class Branch:
    """A sensitivity branch: the same sample judged under a different threshold,
    so that any difference in the result is attributable to that threshold.

    Branches differ in stage-1 admission and, for the panel count only, in
    whether that count votes on privacy at stage 6. The panel needs both levers
    because it reaches a candidate twice — once as an admission rule and once
    inside `is_private` — and relaxing either alone leaves the other filtering.
    Nothing else about a branch may differ, or the comparison stops being
    attributable.
    """
    name: str
    pon_max: int | None = None          # None = this pipeline's PON threshold off
    gnomad_max_af: float | None = None  # None = no population filter

    #: Admit records the CALLER rejected as panel hits (somatic `FILTER=PON`).
    #: Without this, `pon_max=None` is a no-op on somatic call sets: those
    #: records are dropped for not being PASS before `pon_max` is consulted.
    admit_panel_filtered: bool = False

    #: Whether the panel count votes on `is_private` (stage 6). False keeps it as
    #: a reported column, leaving privacy to population frequency alone.
    pon_in_privacy: bool = True


@dataclass
class RunConfig:
    run_name: str
    reference: Reference
    samples: list[Sample]
    branches: list[Branch] = field(default_factory=lambda: [Branch("default")])
    genome_build: str = "GRCh38"
    ensembl_release: int = 115
    resources: dict = field(default_factory=dict)
    criteria_overrides: dict = field(default_factory=dict)
    out_dir: str = "results"
    #: Stage-2 backend name; see src/svneo/generators/. The only pluggable stage.
    peptide_generator: str = "neosv_trace"

    # -- lineage helpers ---------------------------------------------------
    def by_name(self, name: str) -> Sample | None:
        return next((s for s in self.samples if s.name == name), None)

    def ordered_samples(self) -> list[Sample]:
        """Parents before children (topological). A sample whose parent is not in
        the config is treated as a root, with a warning at load time."""
        done, out = set(), []
        pending = list(self.samples)
        while pending:
            progressed = False
            for sample in list(pending):
                if sample.parent is None or sample.parent in done \
                        or self.by_name(sample.parent) is None:
                    out.append(sample)
                    done.add(sample.name)
                    pending.remove(sample)
                    progressed = True
            if not progressed:                    # a cycle: emit deterministically
                out.extend(sorted(pending, key=lambda s: s.name))
                break
        return out

    def ancestors(self, name: str) -> list[str]:
        """Chain from the sample up to the root, used for attribution."""
        chain, seen = [], set()
        node = self.by_name(name)
        while node and node.parent and node.parent not in seen:
            seen.add(node.parent)
            chain.append(node.parent)
            node = self.by_name(node.parent)
        return chain


class ConfigError(Exception):
    """A configuration the user can fix, as opposed to a bug in the pipeline.

    Raised so `run.py` can print the message alone. A traceback in front of
    "/path/to/reference_peptides.tsv does not exist" tells the reader nothing
    they can act on and suggests the tool is broken rather than unconfigured.
    """


def load(path: str) -> RunConfig:
    """Load and validate a run config. Fails loudly: a silently mis-specified
    input is far more expensive than a refused run."""
    with open(path) as handle:
        raw = yaml.safe_load(handle)

    for key in ("run_name", "reference", "samples"):
        if key not in raw:
            raise ConfigError(f"config is missing required key '{key}'")

    reference = Reference(**raw["reference"])
    if not os.path.exists(reference.peptides):
        # Every missing path at once, not the first one. A new user's first
        # action is to copy the template and run --dry-run; reporting one
        # placeholder per invocation turns that into a guessing game.
        missing = [f"reference.peptides: {reference.peptides}"]
        for entry in raw.get("samples", []):
            for key in ("vcf", "rna_bam", "dna_bam", "isofox_dir"):
                value = entry.get(key)
                if value and not os.path.exists(value):
                    missing.append(f"samples[{entry.get('name', '?')}].{key}: {value}")
        raise ConfigError(
            "these configured paths do not exist:\n  "
            + "\n  ".join(missing)
            + "\n\nIf this is a freshly copied config/template.yaml, replace the "
              "/path/to/... placeholders with real paths.")

    samples = [Sample(**s) for s in raw["samples"]]
    if not samples:
        raise ConfigError("config lists no samples")

    names = [s.name for s in samples]
    if len(names) != len(set(names)):
        raise ConfigError("sample names must be unique")
    for sample in samples:
        if sample.vcf_kind not in ("germline", "somatic"):
            raise ValueError(f"{sample.name}: vcf_kind must be 'germline' or "
                             f"'somatic', got '{sample.vcf_kind}'")
        if not os.path.exists(sample.vcf):
            raise FileNotFoundError(f"{sample.name}: VCF not found: {sample.vcf}")
        if sample.parent and sample.parent not in names:
            print(f"[config] WARNING: {sample.name}.parent='{sample.parent}' is "
                  f"not a configured sample; treating {sample.name} as a root")

    branches = [Branch(**b) for b in raw.get("branches", [{"name": "default"}])]

    # A germline VCF shared by several samples is the single most consequential
    # input mistake in this design, so warn about it at load time.
    germline = [s for s in samples if s.vcf_kind == "germline"]
    if len(germline) > 1:
        print(f"[config] WARNING: {len(germline)} samples declare vcf_kind="
              f"'germline' ({', '.join(s.name for s in germline)}). If these are "
              f"a derived lineage, they very likely share ONE parental genome; "
              f"analysing it repeatedly inflates any 'shared across samples' claim.")

    return RunConfig(
        run_name=raw["run_name"], reference=reference, samples=samples,
        branches=branches, genome_build=raw.get("genome", {}).get("build", "GRCh38"),
        ensembl_release=raw.get("genome", {}).get("ensembl_release", 115),
        resources=raw.get("resources", {}),
        criteria_overrides=raw.get("criteria", {}),
        out_dir=raw.get("out_dir", "results"),
        peptide_generator=raw.get("peptide_generator", "neosv_trace"),
    )


def describe(config: RunConfig) -> str:
    """Human-readable summary printed at the start of every run, so the run's
    scope is visible in the log rather than inferred from the config file."""
    lines = [f"run: {config.run_name}",
             f"genome: {config.genome_build} / Ensembl {config.ensembl_release}",
             f"reference: {config.reference.name} ({config.reference.peptides})",
             f"branches: {', '.join(b.name for b in config.branches)}",
             f"peptide generator: {config.peptide_generator}",
             f"samples ({len(config.samples)}), in execution order:"]
    for sample in config.ordered_samples():
        mods = sample.modalities()
        available = ",".join(k for k, v in mods.items() if v) or "none"
        parent = f" <- {sample.parent}" if sample.parent else ""
        lines.append(f"  {sample.name}{parent}  [{sample.vcf_kind}]  {available}")
    return "\n".join(lines)
