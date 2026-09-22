#!/usr/bin/env python3
"""
compare_runs.py — what changed between two output directories, and why it matters.

WHY THIS EXISTS
---------------
When a defect in the generator is corrected, replacing the old results answers
the wrong question. The interesting quantity is not "what do we get now" but
"how much was the defect distorting". The second is evidence — for the analysis,
and for the upstream report — and it is lost the moment the old directory is
overwritten.

So both directories are kept and this compares them, run by run, on the
quantities the defect was expected to move:

  cascade counts       admitted junctions through to surviving events
  frame_effect         the share flagged Start-loss, upstream's own
                       low-reliability warning
  junction spanning    peptides that actually cross the breakpoint
  strand composition   the skew that led to the defect being found

A figure that does NOT move is as informative as one that does: the fix was
verified to leave plus-strand transcripts untouched, so anything driven by them
should be identical.

    python tools/compare_runs.py --before results --after results_patch002
"""
from __future__ import annotations

import argparse
import json
import pathlib

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent

#: Counts worth comparing, in cascade order, with a short label.
COUNTS = [
    ("admitted", "admitted breakends"),
    ("junctions", "junctions"),
    ("candidate_peptides", "candidate peptides"),
    ("matches_unique_peptides", "distinct catalogue matches"),
    ("matches_credible", "credible matches"),
    ("events", "genomic events"),
    ("events_private", "not a common variant"),
    ("events_hc", "high-confidence"),
    ("events_rna_supported", "RNA-supported"),
    ("events_private_hc_and_rna", "all three criteria"),
]


def runs(directory: pathlib.Path) -> list[str]:
    """Run directories, i.e. those carrying a summary.json."""
    if not directory.exists():
        return []
    return sorted(p.name for p in directory.iterdir()
                  if (p / "summary.json").exists())


def counts_of(directory: pathlib.Path, run: str) -> dict:
    summary = json.loads((directory / run / "summary.json").read_text())
    return summary.get("counts", {})


def peptide_stats(directory: pathlib.Path, run: str) -> dict:
    """frame_effect shares and junction spanning, from the generator's output."""
    candidates = next((directory / run).glob("*.all_neopeptides.txt"), None)
    if candidates is None:
        return {}
    table = pd.read_csv(candidates, sep="\t", low_memory=False)
    out = {"peptide_rows": len(table)}

    if "frame_effect" in table.columns:
        share = table["frame_effect"].value_counts(normalize=True)
        out["start_loss"] = float(share.get("Start-loss", 0.0))
        # Split by strand: the defect was strand-specific, so a fix should move
        # the minus share and leave the plus share alone.
        if "strand1" in table.columns:
            for strand, label in (("-", "start_loss_minus"),
                                  ("+", "start_loss_plus")):
                subset = table[table["strand1"] == strand]
                out[label] = float(
                    (subset["frame_effect"] == "Start-loss").mean()) \
                    if len(subset) else float("nan")

    if "spans_junction" in table.columns:
        spans = table["spans_junction"].astype(str).str.lower()
        out["spans_junction"] = float((spans == "true").mean())

    if "strand1" in table.columns:
        strands = table["strand1"].dropna()
        out["fraction_minus"] = float((strands == "-").mean()) if len(strands) else float("nan")

    return out


def fmt(value, kind: str = "int") -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "—"
    if value == "NA":
        return "NA"
    if kind == "pct":
        return f"{100 * float(value):.1f}%"
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return str(value)


def delta(before, after, kind: str = "int") -> str:
    """The change, or an em dash when either side is absent or non-numeric."""
    for value in (before, after):
        if value is None or value == "NA" or \
                (isinstance(value, float) and pd.isna(value)):
            return "—"
    try:
        change = float(after) - float(before)
    except (TypeError, ValueError):
        return "—"
    if change == 0:
        return "unchanged"
    return f"{100 * change:+.1f} pp" if kind == "pct" else f"{change:+,.0f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--before", default="results")
    parser.add_argument("--after", default="results_patch002")
    parser.add_argument("--out", default=None, help="write the report to a file")
    args = parser.parse_args()

    before_dir, after_dir = REPO / args.before, REPO / args.after
    shared = [r for r in runs(after_dir) if r in runs(before_dir)]
    if not shared:
        raise SystemExit(
            f"no run appears in both {args.before} and {args.after}.\n"
            f"  {args.before}: {', '.join(runs(before_dir)) or 'nothing'}\n"
            f"  {args.after}: {', '.join(runs(after_dir)) or 'nothing'}")

    lines = [f"# {args.before} vs {args.after}", "",
             f"{len(shared)} run(s) present in both.", ""]

    for run in shared:
        b, a = counts_of(before_dir, run), counts_of(after_dir, run)
        bp, ap = peptide_stats(before_dir, run), peptide_stats(after_dir, run)
        lines += [f"## {run}", "",
                  "| quantity | before | after | change |", "|---|---|---|---|"]
        for key, label in COUNTS:
            if key not in b and key not in a:
                continue
            lines.append(f"| {label} | {fmt(b.get(key))} | {fmt(a.get(key))} "
                         f"| {delta(b.get(key), a.get(key))} |")
        for key, label in (("start_loss", "Start-loss, all"),
                           ("start_loss_minus", "Start-loss, minus strand"),
                           ("start_loss_plus", "Start-loss, plus strand"),
                           ("spans_junction", "peptides spanning the junction"),
                           ("fraction_minus", "candidates on the minus strand")):
            if key in bp or key in ap:
                lines.append(f"| {label} | {fmt(bp.get(key), 'pct')} "
                             f"| {fmt(ap.get(key), 'pct')} "
                             f"| {delta(bp.get(key), ap.get(key), 'pct')} |")
        lines.append("")

    report = "\n".join(lines)
    print(report)
    if args.out:
        (REPO / args.out).write_text(report + "\n")
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
