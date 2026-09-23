#!/usr/bin/env python
"""
build_sample_neoantigens.py — one sample's own neoantigen catalogue.

A DIFFERENT QUESTION FROM THE RECURRENCE TEST
---------------------------------------------
The pipeline proper asks whether a sample's neopeptides reappear in a reference
cohort, so its outputs are *matches*. This script asks what a sample carries,
whether or not anyone else carries it, so its input is the candidate universe
(`<sample>.all_neopeptides.txt`) and the cohort never enters.

The two differ by an order of magnitude, and the candidate table is the raw
generator output: it has had NO sequence QC. A peptide that occurs in the normal
proteome is not a neoantigen, and nothing upstream has checked. So the checks
that the cross applies to matched peptides are applied here to every candidate,
reusing `svneo.cross.sequence_qc` rather than reimplementing them.

LINEAGE: WHERE A PEPTIDE COMES FROM IS PART OF THE ANSWER
---------------------------------------------------------
A derived clone is analysed from its somatic VCF alone, because its germline
genome IS the parent's. "Neoantigens in clone X" therefore has two readings:
what the clone ACQUIRED (its somatic events), or what its genome CONTAINS (those
plus everything inherited from the parent). The second is the larger and usually
the one that matters for what a cell could present; the first is what is
attributable to whatever was done to make the clone.

Both are emitted, and `origin` says which is which, so neither reading is
imposed. Counts are reported per origin as well as in total, because a total
that mixes them answers neither question cleanly.

`pon_count` is carried as a COLUMN, never applied as a filter. A junction common
in unrelated normals is not a neoantigen, but that is a judgement to make while
looking at the number, not one to lose silently upstream.

Usage
-----
    python tools/build_sample_neoantigens.py \\
        --run-dir  <results>/<clone>          \\
        --parent-run-dir <results>/<parental> \\
        --out-dir  <results>/<clone>/neoantigens
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
from svneo import criteria, cross      # noqa: E402

#: Junction-level evidence to carry onto every peptide of that junction. Taken
#: from the stages 6-8 pass over EVERY admitted junction, not from the credible
#: subset, because a candidate catalogue must cover junctions that never matched.
JUNCTION_COLUMNS = [
    "sample_sv_id", "sv_hc", "pon_count", "pon_fraction", "gnomad_af_used",
    "gnomad_af_popmax", "pass_pon", "pass_gnomad", "is_private", "privacy_note",
    "vf_bp1", "vf_bp2", "qual_bp1", "qual_bp2", "segmapq_bp1", "segmapq_bp2",
    "expressed", "min_side_TPM", "rna_tier", "junction_reads", "test",
    "test_reason",
    # Isofox's own splice calls, independent of this pipeline's read counting.
    # They are here because they catch a confound nothing else does: an SV whose
    # breakpoints sit beside a splice site in a highly expressed gene inherits
    # that site's spliced reads as apparent junction-crossing evidence. The tell
    # is a huge `junction_reads` on tiny `vf_*`, with `nearest_alt_sj_bp` small
    # and `alt_sj_frags` large. Compare the junction_reads/vf ratio against the
    # other events of the same run before believing such a candidate.
    "nearest_alt_sj_bp", "alt_sj_frags", "alt_sj_type", "alt_sj_within_window",
    "retained_intron", "isofox_fusion", "isofox_fusion_support",
    "gene1_TPM", "gene2_TPM",
]


def self_peptide_flags(peptides: pd.Series, proteome: str | None) -> pd.Series:
    """Is each peptide present in the normal proteome?

    `svneo.cross.sequence_qc` asks this with `peptide in proteome`, which is
    correct and fine for the few thousand peptides a cross produces. A candidate
    catalogue is a different scale — tens of thousands of peptides against a
    ~19 MB string is roughly an hour of CPU, since every lookup rescans the whole
    proteome.

    Indexing inverts that: one pass builds the set of proteome substrings at each
    peptide length actually present, after which every lookup is a hash probe.
    Same answer, seconds instead of an hour. The index is built per length so
    memory stays proportional to the proteome, not to its square.
    """
    if proteome is None:
        return pd.Series(pd.NA, index=peptides.index)
    lengths = sorted({len(p) for p in peptides.dropna().unique()})
    flags = pd.Series(False, index=peptides.index)
    for length in lengths:
        wanted = {p for p in peptides.dropna().unique() if len(p) == length}
        if not wanted:
            continue
        # Only substrings that could be a peptide are kept: the NUL separators
        # that load_proteome inserts mark protein boundaries, and a window
        # spanning one is not a real sequence.
        present = {proteome[i:i + length]
                   for i in range(len(proteome) - length + 1)}
        hits = wanted & present
        flags |= peptides.isin(hits)
        print(f"    self-test {length}mer: {len(hits):,} of {len(wanted):,} "
              "are proteome sequences")
    return flags


def read_binding(path: pathlib.Path) -> dict[str, dict[str, dict]]:
    """peptide -> {allele: {margin, cut, allele}} for rows passing every cut.

    The binder rule and the margin are the same ones the patient-side analysis
    uses, imported rather than restated so a peptide cannot be a binder in one
    table and not the other.
    """
    from run_netmhcpan import BINDER_THRESHOLDS
    binding: dict[str, dict[str, dict]] = {}
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    for column in BINDER_THRESHOLDS:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    passes = frame[list(BINDER_THRESHOLDS)].le(pd.Series(BINDER_THRESHOLDS)).all(axis=1)
    kept = frame[passes]
    for row in kept.itertuples(index=False):
        values = {c: getattr(row, c) for c in BINDER_THRESHOLDS}
        distances = {c: (t - values[c]) / t for c, t in BINDER_THRESHOLDS.items()}
        cut = min(distances, key=distances.get)
        allele = str(row.allele).replace("*", "").upper()
        binding.setdefault(row.peptide, {})[allele] = {
            "margin": distances[cut], "cut": cut, "allele": allele}
    print(f"  binding: {len(kept):,} of {len(frame):,} peptide-allele rows pass, "
          f"{len(binding):,} peptides bind >= 1 allele")
    return binding


def read_candidates(run_dir: pathlib.Path, origin: str) -> pd.DataFrame:
    """The generator's candidate peptides for one run, tagged with its origin."""
    hits = sorted(run_dir.glob("*.all_neopeptides.txt"))
    if not hits:
        raise SystemExit(f"no *.all_neopeptides.txt in {run_dir}")
    table = pd.read_csv(hits[0], sep="\t", dtype=str, low_memory=False)
    table["origin"] = origin
    table["origin_run"] = run_dir.name
    return table


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-dir", type=pathlib.Path, required=True,
                        help="the sample's own run directory (somatic events)")
    parser.add_argument("--parent-run-dir", type=pathlib.Path, action="append",
                        default=[], metavar="DIR",
                        help="a parental run whose events this sample inherits; "
                             "repeatable, innermost ancestor first. Omit for a "
                             "somatic-only catalogue.")
    parser.add_argument("--proteome", type=pathlib.Path,
                        help="Ensembl pep.all.fa.gz for the self-peptide test. "
                             "Without it `is_self` is NA — which is NOT the same "
                             "as tested and found non-self.")
    parser.add_argument("--predictions", type=pathlib.Path,
                        help="netMHCpan table for this sample's OWN class I "
                             "genotype, from tools/run_netmhcpan.py. Adds the "
                             "presentability columns; without it they are absent "
                             "rather than blank, because an unrun prediction and "
                             "a negative prediction are different facts.")
    parser.add_argument("--origin", choices=["somatic", "parental", "both", "all"],
                        default="all",
                        help="restrict the catalogue to peptides of this origin. "
                             "`somatic` is what this clone ACQUIRED — the events "
                             "attributable to whatever was done to make it. The "
                             "default keeps everything and leaves the filtering "
                             "to the `peptide_origin` column.")
    parser.add_argument("--out-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()

    frames = [read_candidates(args.run_dir, "somatic")]
    for parent in args.parent_run_dir:
        frames.append(read_candidates(parent, "parental"))
    table = pd.concat(frames, ignore_index=True)
    print(f"  candidate rows: {len(table):,} "
          f"({', '.join(f'{k} {v:,}' for k, v in table['origin'].value_counts().items())})")

    peptide_column = criteria.PEPTIDE_COLUMN
    table = table[table[peptide_column].notna()].copy()

    # A peptide can arise in BOTH the parent and the clone. That is not a
    # duplicate to drop: it says the peptide is inherited AND regenerated, so the
    # origins are collapsed into one label per peptide rather than one row being
    # discarded arbitrarily.
    origins = (table.groupby(peptide_column)["origin"]
                    .agg(lambda s: "both" if s.nunique() > 1 else s.iloc[0]))
    table["peptide_origin"] = table[peptide_column].map(origins)

    # -- sequence QC, the checks the raw candidate table has never had ---------
    proteome = None
    if args.proteome:
        print(f"  loading proteome {args.proteome.name}")
        proteome = cross.load_proteome(str(args.proteome))

    table["low_complexity"] = table[peptide_column].map(criteria.is_low_complexity)
    table["is_self"] = self_peptide_flags(table[peptide_column], proteome)

    # -- junction-level evidence ---------------------------------------------
    evidence = []
    for run_dir in [args.run_dir, *args.parent_run_dir]:
        path = run_dir / "stage6_8_all_junctions.tsv"
        if not path.exists():
            print(f"  WARNING: no stage6_8_all_junctions.tsv in {run_dir} — "
                  "confidence, privacy and RNA columns will be blank for it")
            continue
        frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        keep = [c for c in JUNCTION_COLUMNS if c in frame.columns]
        evidence.append(frame[keep].drop_duplicates("sample_sv_id"))
    if evidence:
        merged = pd.concat(evidence, ignore_index=True).drop_duplicates("sample_sv_id")
        table["sample_sv_id"] = table["sv_id"].astype(str)
        table = table.merge(merged, on="sample_sv_id", how="left")

    # -- presentability by this sample's own MHC ------------------------------
    # For a cell line this is the autologous question and there is no panel
    # branch to confuse it with: the alleles ARE the line's genotype, so a binder
    # here is "this cell could present this peptide", not "somebody could".
    presentable = None
    if args.predictions:
        binding = read_binding(args.predictions)
        table["n_alleles_binding"] = table[peptide_column].map(
            lambda p: len(binding.get(p, ())))
        table["binding_alleles"] = table[peptide_column].map(
            lambda p: ";".join(sorted(binding.get(p, ()))))
        table["presentable"] = table["n_alleles_binding"] > 0
        best = table[peptide_column].map(
            lambda p: max(binding.get(p, {}).values(),
                          key=lambda v: v["margin"], default=None))
        table["presentation_margin"] = [b["margin"] if b else "" for b in best]
        table["presentation_limiting_cut"] = [b["cut"] if b else "" for b in best]
        table["presentation_best_allele"] = [b["allele"] if b else "" for b in best]
        table["presentation_robustness"] = [
            "" if not b else
            "flippable" if b["margin"] <= 0.10 else
            "marginal" if b["margin"] <= 0.25 else
            "solid" if b["margin"] <= 0.50 else "robust" for b in best]
        presentable = table.drop_duplicates(peptide_column)

    if args.origin != "all":
        before = table[peptide_column].nunique()
        table = table[table["peptide_origin"] == args.origin]
        print(f"  origin={args.origin}: {table[peptide_column].nunique():,} "
              f"peptides of {before:,}")
        if table.empty:
            raise SystemExit(f"no peptides with origin={args.origin}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out_dir / "candidate_neoantigens.tsv"
    table.to_csv(out, sep="\t", index=False)

    peptides = table.drop_duplicates(peptide_column)
    (args.out_dir / "candidate_peptides.txt").write_text(
        "\n".join(sorted(peptides[peptide_column])) + "\n")

    # -- funnel, per origin and in total --------------------------------------
    def flag(frame, column, value="True"):
        return frame[column].astype(str).str.strip().str.lower() == value.lower()

    # THE SEQUENCE CUTS CHAIN; THE EVIDENCE CUTS DO NOT.
    # Sequence QC is cumulative: each row is a subset of the one above it. The
    # evidence rows are NOT — they are marginal counts over the same QC'd set,
    # because SV confidence, transcription and junction reads are independent
    # properties and chaining them would hide which one actually bites. Printing
    # them in one column would make the funnel appear to widen (more peptides
    # transcribed than on a high-confidence SV), so they are labelled and the
    # conjunctions are given their own rows, which are the ones to quote.
    steps = []
    surviving = peptides
    steps.append(("candidate_peptides", surviving, "cumulative"))
    if proteome is not None:
        surviving = surviving[~flag(surviving, "is_self")]
        steps.append(("non_self", surviving, "cumulative"))
    surviving = surviving[~flag(surviving, "low_complexity")]
    steps.append(("high_complexity", surviving, "cumulative"))

    qc_passed = surviving
    has = lambda c: c in qc_passed.columns                       # noqa: E731
    hc = flag(qc_passed, "sv_hc") if has("sv_hc") else None
    expressed = flag(qc_passed, "expressed") if has("expressed") else None
    rna = (qc_passed["rna_tier"].isin(["STRONG", "SUGGESTIVE", "WEAK"])
           if has("rna_tier") else None)
    rna_strong = (qc_passed["rna_tier"].isin(["STRONG", "SUGGESTIVE"])
                  if has("rna_tier") else None)

    for name, mask in (("on_high_confidence_sv", hc), ("transcribed", expressed),
                       ("rna_junction_reads", rna)):
        if mask is not None:
            steps.append((name, qc_passed[mask], "marginal"))
    if "presentable" in qc_passed.columns:
        pres = flag(qc_passed, "presentable")
        robust = qc_passed["presentation_robustness"] == "robust"
        steps.append(("presentable_by_own_MHC", qc_passed[pres], "marginal"))
        steps.append(("presentable_robustly", qc_passed[pres & robust], "marginal"))
    if hc is not None and rna_strong is not None:
        steps.append(("hc_AND_rna", qc_passed[hc & rna_strong], "conjunction"))
        if expressed is not None:
            whole = hc & expressed & rna_strong
            steps.append(("hc_AND_transcribed_AND_rna",
                          qc_passed[whole], "conjunction"))
            if "presentable" in qc_passed.columns:
                steps.append(("AND_presentable",
                              qc_passed[whole & flag(qc_passed, "presentable")],
                              "conjunction"))
                steps.append(("AND_presentable_robustly",
                              qc_passed[whole & flag(qc_passed, "presentable")
                                        & (qc_passed["presentation_robustness"]
                                           == "robust")], "conjunction"))

    rows = []
    for name, frame, kind in steps:
        counts = frame["peptide_origin"].value_counts()
        rows.append({"step": name, "kind": kind, "n": len(frame),
                     "n_somatic": int(counts.get("somatic", 0)),
                     "n_parental": int(counts.get("parental", 0)),
                     "n_both": int(counts.get("both", 0))})
    funnel = pd.DataFrame(rows)
    funnel.to_csv(args.out_dir / "funnel_neoantigens.tsv", sep="\t", index=False)

    with (args.out_dir / "neoantigen_summary.json").open("w") as fh:
        json.dump({"run_dir": str(args.run_dir),
                   "parent_run_dirs": [str(p) for p in args.parent_run_dir],
                   "self_test": bool(proteome),
                   "funnel": rows}, fh, indent=2)

    print()
    print(f"  {'step':<28} {'kind':<12} {'total':>8} {'somatic':>8} "
          f"{'parental':>9} {'both':>6}")
    for row in rows:
        print(f"  {row['step']:<28} {row['kind']:<12} {row['n']:>8,} "
              f"{row['n_somatic']:>8,} {row['n_parental']:>9,} {row['n_both']:>6,}")
    print("  cumulative rows nest; marginal rows do NOT — they are independent "
          "properties of the same QC'd set. Quote a conjunction row.")
    print(f"\n  wrote {out}")
    print(f"        {args.out_dir / 'candidate_peptides.txt'} "
          f"({len(peptides):,} peptides, ready for netMHCpan)")
    print(f"        {args.out_dir / 'funnel_neoantigens.tsv'}")


if __name__ == "__main__":
    main()
