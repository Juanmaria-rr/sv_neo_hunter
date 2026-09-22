#!/usr/bin/env python3
"""
analyse_locus_vs_peptide.py — how often is a shared locus a shared peptide?

THE QUESTION
------------
Junctions recur between patients and cell lines far more often than the peptides
they produce do. That is not a contradiction and it is not a bug, but it is the
single easiest way to overstate a recurrence result: quoting "we found matching
SVs" when what matched was a locus, not a consequence.

This measures the gap, and separates its causes, so the three cross levels can be
justified with numbers rather than asserted:

    same locus      the breakpoint falls near a patient breakpoint
    same mechanism  and in the same gene, with the same SV type
    same consequence  and the peptide sequence is identical

WHAT IT SEPARATES
-----------------
Among junctions landing on or beside a patient breakpoint, three different things
prevent a peptide match, and they are not interchangeable:

    (a) the junction produces no peptide at all — intergenic, intronic without
        protein consequence, or identical to the wild type. It could not match.
    (b) it produces peptides and none matches — the interesting case: same place,
        different sequence.
    (c) it produces peptides and at least one matches.

For (b) the report also tests the mechanism we expect to dominate: inserted bases
at the junction. Two breaks at the same coordinate with different insertions
produce different fusion sequences, and an insertion that is not a multiple of
three shifts the reading frame, so everything downstream translates differently.

    python tools/analyse_locus_vs_peptide.py
    python tools/analyse_locus_vs_peptide.py --branch PON10
"""
from __future__ import annotations

import argparse
import os
import pathlib

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent

#: Evidence levels that mean "seen at this locus" without meaning "same peptide".
LOCUS_LEVELS = ["same_gene_and_svtype", "same_gene", "breakpoint_within_1kb",
                "breakpoint_within_10kb", "breakpoint_within_50kb",
                "breakpoint_within_100kb"]

#: Proximity at which a junction is treated as landing on a patient breakpoint
#: for the cause breakdown. Reported, not a criterion used by the pipeline.
NEAR_BP = 1000


def _num(frame: pd.DataFrame, column: str) -> pd.Series:
    """A numeric column, with absence as NaN rather than a raising cast."""
    if column not in frame.columns:
        return pd.Series(float("nan"), index=frame.index)
    return pd.to_numeric(frame[column], errors="coerce")


def load(branch: str) -> pd.DataFrame:
    """The combined SV table for one branch, or the full one as a fallback."""
    for name in (f"master_sv_{branch}.tsv", "master_sv.tsv"):
        path = REPO / "results" / name
        if path.exists():
            frame = pd.read_csv(path, sep="\t", low_memory=False)
            if "branch" in frame.columns:
                subset = frame[frame["branch"] == branch]
                if len(subset):
                    return subset
            return frame
    raise SystemExit("no master_sv table found; run tools/build_master_table.py")


def analyse(table: pd.DataFrame) -> dict:
    """Every figure the report quotes, computed once, from the table."""
    generated = _num(table, "n_peptides_generated").fillna(0)
    matched = _num(table, "n_peptides_matched").fillna(0)
    distance = _num(table, "patient_bp_dist_bp")
    inserted = _num(table, "insert_len").fillna(0)

    evidence = table.get("patient_evidence", pd.Series(dtype=str))
    locus_only = int(evidence.isin(LOCUS_LEVELS).sum())
    identical = int((evidence == "identical_peptide").sum())

    near = distance <= NEAR_BP
    produces = generated > 0
    result = {
        "junctions": len(table),
        "levels": evidence.value_counts().to_dict(),
        "locus_only": locus_only,
        "identical": identical,
        "ratio": locus_only / identical if identical else float("nan"),
        "near": int(near.sum()),
        "near_no_peptide": int((near & ~produces).sum()),
        "near_peptide_no_match": int((near & produces & (matched == 0)).sum()),
        "near_peptide_match": int((near & produces & (matched > 0)).sum()),
    }

    # The same-coordinate subset: the strongest form of "same place".
    exact = distance == 0
    result["exact"] = int(exact.sum())
    result["exact_with_insert"] = int((exact & (inserted > 0)).sum())
    result["exact_insert_median"] = float(inserted[exact & (inserted > 0)].median()) \
        if (exact & (inserted > 0)).any() else float("nan")
    result["exact_same_svtype"] = int(
        (exact & (table.get("svtype") == table.get("patient_svtype"))).sum())
    result["exact_produces"] = int((exact & produces).sum())
    result["exact_matches"] = int((exact & produces & (matched > 0)).sum())

    # Does the inserted length separate the matching from the non-matching?
    at_locus = exact & produces
    result["insert_median_unmatched"] = float(
        inserted[at_locus & (matched == 0)].median()) if (at_locus & (matched == 0)).any() else float("nan")
    result["insert_median_matched"] = float(
        inserted[at_locus & (matched > 0)].median()) if (at_locus & (matched > 0)).any() else float("nan")

    # Even a matching junction matches only a fraction of its own peptides.
    hit = matched > 0
    fraction = (matched[hit] / generated[hit]).replace([float("inf")], float("nan"))
    result["match_fraction_median"] = float(fraction.median()) if len(fraction) else float("nan")
    result["match_fraction_min"] = float(fraction.min()) if len(fraction) else float("nan")

    examples = table[at_locus & (matched == 0)].copy()
    examples["_insert"] = inserted[examples.index]
    result["examples"] = examples.sort_values("_insert", ascending=False)[
        [c for c in ("cell_line", "gene1", "svtype", "patient_svtype",
                     "event_size", "insert_len", "n_peptides_generated")
         if c in examples.columns]].head(8)
    return result


def steps(table: pd.DataFrame) -> list[dict]:
    """The derivation, one row per step: population, operation, result.

    Written as data rather than prose so the report cannot describe a step the
    code does not take. Each entry states the SET it operates on and the COLUMNS
    it reads, which is what makes a number checkable against the master table.
    """
    generated = _num(table, "n_peptides_generated").fillna(0)
    matched = _num(table, "n_peptides_matched").fillna(0)
    distance = _num(table, "patient_bp_dist_bp")
    inserted = _num(table, "insert_len").fillna(0)
    evidence = table.get("patient_evidence", pd.Series(dtype=str))

    near, produces, hit = distance <= NEAR_BP, generated > 0, matched > 0
    exact = distance == 0
    rows = [
        dict(step="0", population="every admitted junction, all lines",
             operation="rows of master_sv.tsv", columns="—", n=len(table)),
        dict(step="1", population="the same junctions",
             operation="distance to the nearest patient breakpoint, "
                       "min over the two breakends",
             columns="patient_bp_dist_bp", n=int(distance.notna().sum())),
        dict(step="2", population="the same junctions",
             operation="assign the highest evidence level reached, by precedence",
             columns="patient_evidence",
             n=int(evidence.notna().sum())),
        dict(step="3a", population="all junctions",
             operation="reach a locus or mechanism level, but not a peptide",
             columns="patient_evidence",
             n=int(evidence.isin(LOCUS_LEVELS).sum())),
        dict(step="3b", population="all junctions",
             operation="reach an identical peptide",
             columns="patient_evidence",
             n=int((evidence == "identical_peptide").sum())),
        dict(step="4", population="all junctions",
             operation=f"restrict to within {NEAR_BP:,} bp of a patient breakpoint",
             columns="patient_bp_dist_bp", n=int(near.sum())),
        dict(step="5a", population="step 4",
             operation="produce no candidate peptide — cannot match",
             columns="n_peptides_generated", n=int((near & ~produces).sum())),
        dict(step="5b", population="step 4",
             operation="produce peptides, none identical to the catalogue",
             columns="n_peptides_generated, n_peptides_matched",
             n=int((near & produces & ~hit).sum())),
        dict(step="5c", population="step 4",
             operation="produce peptides and at least one matches",
             columns="n_peptides_matched", n=int((near & produces & hit).sum())),
        dict(step="6", population="all junctions",
             operation="restrict to the SAME coordinate as a patient breakpoint",
             columns="patient_bp_dist_bp", n=int(exact.sum())),
        dict(step="7a", population="step 6",
             operation="also the same SV type as the patient event",
             columns="svtype, patient_svtype",
             n=int((exact & (table.get("svtype") == table.get("patient_svtype"))).sum())),
        dict(step="7b", population="step 6", operation="carry inserted bases",
             columns="insert_len", n=int((exact & (inserted > 0)).sum())),
        dict(step="8a", population="step 6", operation="produce peptides",
             columns="n_peptides_generated", n=int((exact & produces).sum())),
        dict(step="8b", population="step 8a", operation="share a peptide",
             columns="n_peptides_matched", n=int((exact & produces & hit).sum())),
    ]
    return rows


def steps_markdown(table: pd.DataFrame) -> str:
    """The derivation as a table."""
    lines = ["| Step | Population | Operation | Columns read | n |",
             "|---|---|---|---|---|"]
    for row in steps(table):
        lines.append(f"| {row['step']} | {row['population']} | {row['operation']} "
                     f"| `{row['columns']}` | **{row['n']:,}** |")
    return "\n".join(lines)


def per_line_markdown(table: pd.DataFrame) -> str:
    """The evidence ladder per cell line, since the lines differ by 80-fold."""
    if "cell_line" not in table.columns:
        return ""
    order = ["identical_peptide", "same_gene_and_svtype", "same_gene",
             "breakpoint_within_1kb", "breakpoint_within_10kb",
             "breakpoint_within_100kb", "not_seen_in_patients"]
    cross = pd.crosstab(table["patient_evidence"], table["cell_line"])
    lines_ = list(cross.columns)
    out = ["| Level | " + " | ".join(lines_) + " |",
           "|---" * (len(lines_) + 1) + "|"]
    for level in order:
        if level not in cross.index:
            continue
        out.append(f"| `{level}` | "
                   + " | ".join(f"{cross.loc[level, c]:,}" for c in lines_) + " |")
    out.append("| **total** | "
               + " | ".join(f"**{int(cross[c].sum()):,}**" for c in lines_) + " |")
    return "\n".join(out)


def distance_markdown(table: pd.DataFrame) -> str:
    """How the junctions distribute over distance, before any level is assigned."""
    distance = _num(table, "patient_bp_dist_bp")
    bands = [("0 (same coordinate)", distance == 0),
             ("1 - 1,000", (distance >= 1) & (distance <= 1_000)),
             ("1,001 - 10,000", (distance > 1_000) & (distance <= 10_000)),
             ("10,001 - 100,000", (distance > 10_000) & (distance <= 100_000)),
             ("> 100,000", distance > 100_000),
             ("no patient breakpoint on that chromosome", distance.isna())]
    total = max(len(table), 1)
    lines = ["| Distance to nearest patient breakpoint | Junctions | % |",
             "|---|---|---|"]
    for label, mask in bands:
        lines.append(f"| {label} | {int(mask.sum()):,} | {100*mask.sum()/total:.1f} |")
    return "\n".join(lines)


def markdown(result: dict) -> str:
    """The tables the executive summary embeds."""
    order = ["identical_peptide", "same_gene_and_svtype", "same_gene",
             "breakpoint_within_1kb", "breakpoint_within_10kb",
             "breakpoint_within_50kb", "breakpoint_within_100kb",
             "not_seen_in_patients"]
    total = max(result["junctions"], 1)
    lines = ["| Level reached | Junctions | % |", "|---|---|---|"]
    for level in order:
        n = result["levels"].get(level)
        if n is None:
            continue
        # Emphasis must sit OUTSIDE the code span, or the asterisks render
        # literally inside it.
        cell = f"**`{level}`**" if level == "identical_peptide" else f"`{level}`"
        lines.append(f"| {cell} | {n:,} | {100 * n / total:.1f} |")
    lines.append(f"| *total* | {result['junctions']:,} | |")

    lines += ["", "| Junctions within "
              f"{NEAR_BP:,} bp of a patient breakpoint | n |", "|---|---|",
              f"| produce no peptide at all — cannot match by construction "
              f"| {result['near_no_peptide']:,} |",
              f"| **produce peptides, none matches** | "
              f"**{result['near_peptide_no_match']:,}** |",
              f"| produce peptides and at least one matches | "
              f"{result['near_peptide_match']:,} |",
              f"| *total* | {result['near']:,} |"]

    lines += ["", "| At the *same coordinate* as a patient breakpoint | n |",
              "|---|---|",
              f"| junctions | {result['exact']:,} |",
              f"| of which same SV type as the patient event | "
              f"{result['exact_same_svtype']:,} |",
              f"| of which carry inserted bases at the junction | "
              f"{result['exact_with_insert']:,} "
              f"(median {result['exact_insert_median']:.0f} bp) |",
              f"| of which produce peptides | {result['exact_produces']:,} |",
              f"| **of which share a peptide** | **{result['exact_matches']:,}** |"]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--branch", default="noPON")
    parser.add_argument("--out", default="results/locus_vs_peptide.tsv")
    args = parser.parse_args()

    table = load(args.branch)
    result = analyse(table)

    print(f"branch {args.branch}: {result['junctions']:,} admitted junctions\n")
    print(markdown(result))
    print(f"\nlocus-level recurrence : {result['locus_only']:,}")
    print(f"peptide-level recurrence: {result['identical']:,}")
    print(f"ratio                   : {result['ratio']:.0f}:1")
    print(f"\ninserted bases at a shared coordinate, median:")
    print(f"  junctions whose peptides do NOT match: "
          f"{result['insert_median_unmatched']:.0f} bp")
    print(f"  junctions that share a peptide       : "
          f"{result['insert_median_matched']:.0f} bp")
    print(f"\nfraction of a matching junction's own peptides that match: "
          f"median {result['match_fraction_median']:.3f}, "
          f"min {result['match_fraction_min']:.3f}")
    print("\nsame coordinate, peptides produced, none matching:")
    print(result["examples"].to_string(index=False))

    out = REPO / args.out
    rows = [{"metric": k, "value": v} for k, v in result.items()
            if not isinstance(v, (dict, pd.DataFrame))]
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"\nwrote {os.path.relpath(out, REPO)}")


if __name__ == "__main__":
    main()
