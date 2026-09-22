#!/usr/bin/env python
"""
compare_netmhcpan_versions.py — how much does netMHCpan 4.1 vs 4.2 change?

WHY
---
An earlier analysis may have been run with a different netMHCpan version. The
binder rule depends on two %Rank cuts, and a rank is read off a version-specific
model, so "the same filter" is not free across versions.

Measuring that does not need re-predicting. The two runs overlap on real
(peptide, allele) pairs -- matched peptides are catalogue peptides, and each
patient's own alleles sit inside the predicted panel -- so this is a join, not a
run.

Which version produced an existing output is not a guess: netMHCpan embeds its
full command line in the header of every plain-text result file.

WHAT IT REPORTS
---------------
1. Pair level — the 2x2 of binder / non-binder under each version, over every
   (peptide, allele) pair present in both, plus per-cut agreement so a
   disagreement can be attributed to affinity, %Rank_BA or %Rank_EL.
2. Distributional — median and 90th-percentile absolute differences in each
   quantity, to distinguish "shifted slightly" from "reordered".
3. Peptide level — the number that matters for the deliverable: for peptides
   whose whole comparable allele set is shared, does "binds at least one of
   them" flip between versions?

WHAT IT CANNOT SAY
------------------
The overlap is whatever the two runs happen to share; it is not a random sample
of the peptide set. A peptide predicted here but absent from the older outputs
contributes nothing, so the concordance describes the shared region only.

Usage
-----
    python tools/compare_netmhcpan_versions.py --limit-files 100   # pilot
    python tools/compare_netmhcpan_versions.py                     # all of them
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import pathlib
import statistics
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from run_netmhcpan import BINDER_THRESHOLDS      # noqa: E402

REPO = pathlib.Path(__file__).resolve().parent.parent
DEFAULT_ROOT = REPO.parent

#: netMHCpan's plain-text table, 0-indexed. VERIFIED to be identical in 4.1 and
#: 4.2 against real output of each -- the header is
#: `Pos MHC Peptide Core Of Gp Gl Ip Il Icore Identity Score_EL %Rank_EL
#:  Score_BA %Rank_BA Aff(nM) [BindLevel]`
#: and the trailing BindLevel columns are absent on non-binders, so every field
#: used here is counted from the left.
COL_ALLELE, COL_PEPTIDE = 1, 2
COL_RANK_EL, COL_RANK_BA, COL_AFF = 12, 14, 15
MIN_FIELDS = 16


def normalise_allele(raw: str) -> str:
    return raw.strip().upper().replace("*", "").replace(" ", "")


def is_binder(aff: float, rank_ba: float, rank_el: float) -> bool:
    return (aff <= BINDER_THRESHOLDS["affinity_nM"]
            and rank_ba <= BINDER_THRESHOLDS["rank_BA"]
            and rank_el <= BINDER_THRESHOLDS["rank_EL"])


def scan_old(paths: list[pathlib.Path], wanted: set[str]) -> dict:
    """(peptide, allele) -> (rank_EL, rank_BA, aff) from netMHCpan 4.1 output.

    The same pair recurs across patients that share an allele. The prediction is
    deterministic, so a repeat must agree; disagreements are counted rather than
    silently overwritten, because they would mean the parse is wrong.
    """
    old: dict[tuple[str, str], tuple[float, float, float]] = {}
    clashes = 0
    for n, path in enumerate(paths, 1):
        try:
            with path.open(errors="replace") as fh:
                for line in fh:
                    if not line.startswith(" "):
                        continue
                    parts = line.split()
                    if len(parts) < MIN_FIELDS:
                        continue
                    peptide = parts[COL_PEPTIDE]
                    if peptide not in wanted:
                        continue
                    try:
                        values = (float(parts[COL_RANK_EL]),
                                  float(parts[COL_RANK_BA]),
                                  float(parts[COL_AFF]))
                    except ValueError:
                        continue
                    key = (peptide, normalise_allele(parts[COL_ALLELE]))
                    if key in old:
                        if old[key] != values:
                            clashes += 1
                    else:
                        old[key] = values
        except OSError:
            continue
        if n % 250 == 0:
            print(f"    {n:,}/{len(paths):,} files, {len(old):,} pairs")
    if clashes:
        print(f"  WARNING: {clashes:,} repeated pairs disagreed with themselves "
              "-- the 4.1 parse is suspect, do not trust the numbers below")
    return old


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", default="parental_noPON")
    parser.add_argument("--root", type=pathlib.Path, default=DEFAULT_ROOT)
    parser.add_argument("--old-dir", type=pathlib.Path,
                        help="directory tree of the older run's .net.out.txt files")
    parser.add_argument("--limit-files", type=int, default=0,
                        help="pilot on this many 4.1 output files (0 = all)")
    parser.add_argument("--out-dir", type=pathlib.Path)
    args = parser.parse_args()

    cross = args.root / "results_cross_patch002" / args.sample
    if args.old_dir is None:
        raise SystemExit("--old-dir is required: the tree of .net.out.txt files "
                         "produced by the older netMHCpan run")
    old_dir = args.old_dir
    out_dir = args.out_dir or cross
    for path in (cross, old_dir):
        if not path.exists():
            raise SystemExit(f"missing input: {path}")

    with (cross / "netmhcpan_peptide_summary.tsv").open(newline="") as fh:
        wanted = {r["peptide"] for r in csv.DictReader(fh, delimiter="\t")}
    print(f"matched peptides: {len(wanted):,}")

    paths = sorted(old_dir.glob("*/*.net.out.txt"))
    if args.limit_files:
        paths = paths[:args.limit_files]
    print(f"scanning {len(paths):,} netMHCpan-4.1 output files")
    old = scan_old(paths, wanted)
    print(f"  4.1 pairs covering a matched peptide: {len(old):,}")
    if not old:
        raise SystemExit("no overlapping pairs -- nothing to compare")

    # -- stream our 4.2 table, keeping only pairs 4.1 also has ---------------
    rows = []
    with (cross / "netmhcpan.tsv").open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row["peptide"], normalise_allele(row["allele"]))
            ref = old.get(key)
            if ref is None:
                continue
            try:
                new = (float(row["rank_EL"]), float(row["rank_BA"]),
                       float(row["affinity_nM"]))
            except (ValueError, TypeError):
                continue
            rows.append((key, ref, new))
    print(f"  pairs present in BOTH versions: {len(rows):,}")
    if not rows:
        raise SystemExit("no overlapping pairs -- nothing to compare")

    # -- pair level ---------------------------------------------------------
    table = collections.Counter()
    per_cut = collections.Counter()
    deltas = {"rank_EL": [], "rank_BA": [], "affinity_nM": []}
    for key, (el1, ba1, aff1), (el2, ba2, aff2) in rows:
        b1, b2 = is_binder(aff1, ba1, el1), is_binder(aff2, ba2, el2)
        table[(b1, b2)] += 1
        per_cut["aff"] += (aff1 <= 500) == (aff2 <= 500)
        per_cut["rank_BA"] += (ba1 <= 2) == (ba2 <= 2)
        per_cut["rank_EL"] += (el1 <= 2) == (el2 <= 2)
        deltas["rank_EL"].append(abs(el2 - el1))
        deltas["rank_BA"].append(abs(ba2 - ba1))
        deltas["affinity_nM"].append(abs(aff2 - aff1))

    n = len(rows)
    both = table[(True, True)]
    only41, only42 = table[(True, False)], table[(False, True)]
    neither = table[(False, False)]
    agree = both + neither

    def pct(x):
        return 100.0 * x / n

    def quantiles(values):
        values = sorted(values)
        return (statistics.median(values), values[int(0.9 * (len(values) - 1))])

    # -- is the disagreement substantive, or is it threshold noise? ---------
    # Every cut is a hard boundary, so a pair sitting at 505 nM under one version
    # and 495 under the other flips the verdict while saying nothing about the
    # models disagreeing. `margin` is the smallest relative distance to any of
    # the three thresholds, taken over BOTH versions: small margin = the call was
    # always going to be fragile.
    def margin(el, ba, aff):
        return min(abs(aff - BINDER_THRESHOLDS["affinity_nM"]) / BINDER_THRESHOLDS["affinity_nM"],
                   abs(ba - BINDER_THRESHOLDS["rank_BA"]) / BINDER_THRESHOLDS["rank_BA"],
                   abs(el - BINDER_THRESHOLDS["rank_EL"]) / BINDER_THRESHOLDS["rank_EL"])

    margins = []
    for key, (el1, ba1, aff1), (el2, ba2, aff2) in rows:
        if is_binder(aff1, ba1, el1) != is_binder(aff2, ba2, el2):
            margins.append(min(margin(el1, ba1, aff1), margin(el2, ba2, aff2)))
    marginal = {f"within_{int(t * 100)}pct_of_a_threshold":
                sum(m <= t for m in margins) for t in (0.10, 0.25, 0.50)}

    # -- peptide level ------------------------------------------------------
    # Only peptides whose entire 4.1 allele set is also in 4.2, so "binds at
    # least one" is asked of the same alleles on both sides.
    old_alleles = collections.defaultdict(set)
    for peptide, allele in old:
        old_alleles[peptide].add(allele)
    shared = collections.defaultdict(set)
    binder41 = collections.defaultdict(bool)
    binder42 = collections.defaultdict(bool)
    for (peptide, allele), (el1, ba1, aff1), (el2, ba2, aff2) in rows:
        shared[peptide].add(allele)
        binder41[peptide] |= is_binder(aff1, ba1, el1)
        binder42[peptide] |= is_binder(aff2, ba2, el2)
    complete = [p for p in shared if shared[p] == old_alleles[p]]
    flips = collections.Counter((binder41[p], binder42[p]) for p in complete)

    # THE DENOMINATOR DECIDES THE HEADLINE. Agreement over all pairs is
    # dominated by peptide-allele combinations that are non-binders under both
    # versions -- the overwhelming majority -- and reads ~99.8% however much the
    # binder calls actually move. The honest denominator is the union of calls
    # that are a binder under AT LEAST ONE version. Both are reported, and the
    # second is the one to quote.
    union_pairs = both + only41 + only42
    union_peps = (flips[(True, True)] + flips[(True, False)]
                  + flips[(False, True)])
    summary = {
        "pairs_compared": n,
        "peptides_touched": len(shared),
        "peptides_without_any_4.1_prediction": len(wanted) - len(shared),
        "peptides_fully_shared": len(complete),
        "agreement_all_pairs_pct": round(pct(agree), 3),
        "agreement_among_binder_calls_pct": (
            round(100.0 * both / union_pairs, 2) if union_pairs else None),
        "binder_calls_union": union_pairs,
        # Three defensible denominators. The union is the strictest; the two
        # directional rates say how much of each version's binder set the other
        # one keeps, and they are what a reader usually means by "agreement".
        "retained_by_4.2_of_4.1_binders_pct": (
            round(100.0 * both / (both + only41), 2) if both + only41 else None),
        "retained_by_4.1_of_4.2_binders_pct": (
            round(100.0 * both / (both + only42), 2) if both + only42 else None),
        "discordant_pairs": len(margins),
        "discordant_marginality": marginal,
        "pair_table": {"binder_in_both": both,
                       "binder_only_4.1": only41,
                       "binder_only_4.2": only42,
                       "binder_in_neither": neither},
        "per_cut_agreement_pct": {k: round(100.0 * v / n, 3)
                                  for k, v in per_cut.items()},
        "abs_delta_median_p90": {k: [round(m, 4), round(p90, 4)]
                                 for k, (m, p90) in
                                 ((k, quantiles(v)) for k, v in deltas.items())},
        "peptide_flips": {"binder_in_both": flips[(True, True)],
                          "binder_only_4.1": flips[(True, False)],
                          "binder_only_4.2": flips[(False, True)],
                          "binder_in_neither": flips[(False, False)]},
        "peptide_agreement_among_binders_pct": (
            round(100.0 * flips[(True, True)] / union_peps, 2)
            if union_peps else None),
        "peptide_binders_4.1": flips[(True, True)] + flips[(True, False)],
        "peptide_binders_4.2": flips[(True, True)] + flips[(False, True)],
        "binder_thresholds": BINDER_THRESHOLDS,
        "files_scanned": len(paths),
        "pilot": bool(args.limit_files),
    }
    # Per-peptide coverage, so the hole in the 4.1 side can be characterised:
    # the old outputs predate patch 002, so peptides the fix created have no 4.1
    # prediction at all, and that absence is not random.
    with (out_dir / "netmhcpan_version_peptides.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["peptide", "n_alleles_compared", "binder_4.1", "binder_4.2",
                    "covered_by_4.1"])
        for peptide in sorted(wanted):
            n_alleles = len(shared.get(peptide, ()))
            w.writerow([peptide, n_alleles, binder41.get(peptide, ""),
                        binder42.get(peptide, ""), bool(n_alleles)])

    with (out_dir / "netmhcpan_version_concordance.json").open("w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"\n  pairs compared              {n:,}")
    print(f"  agreement, ALL pairs        {pct(agree):.2f}%   <- inflated by "
          f"{neither:,} agreed non-binders, do not quote")
    if union_pairs:
        print(f"  agreement, binder calls     "
              f"{100.0 * both / union_pairs:.1f}%   <- {both:,}/{union_pairs:,}"
              " binder in >=1 version")
    print(f"    binder in both            {both:,}")
    print(f"    binder only in 4.1        {only41:,}  ({pct(only41):.3f}%)")
    print(f"    binder only in 4.2        {only42:,}  ({pct(only42):.3f}%)")
    print(f"    binder in neither         {neither:,}")
    if both + only41:
        print(f"  of 4.1 binders, kept by 4.2 {100.0 * both / (both + only41):.1f}%"
              f"   ({both:,}/{both + only41:,})")
    if both + only42:
        print(f"  of 4.2 binders, kept by 4.1 {100.0 * both / (both + only42):.1f}%"
              f"   ({both:,}/{both + only42:,})")
    if margins:
        print(f"  discordant pairs            {len(margins):,}, of which "
              + ", ".join(f"{v:,} ({100.0 * v / len(margins):.0f}%) {k.replace('_', ' ')}"
                          for k, v in marginal.items()))
    print("  agreement per cut          " + ", ".join(
        f"{k} {100.0 * v / n:.2f}%" for k, v in per_cut.items()))
    print("  |delta| median / p90       " + ", ".join(
        f"{k} {m:.3g}/{p90:.3g}" for k, (m, p90) in
        ((k, quantiles(v)) for k, v in deltas.items())))
    print(f"\n  peptides fully shared       {len(complete):,} of {len(shared):,}")
    print(f"    binder under both         {flips[(True, True)]:,}")
    print(f"    binder only under 4.1     {flips[(True, False)]:,}")
    print(f"    binder only under 4.2     {flips[(False, True)]:,}")
    if union_peps:
        print(f"    agreement among binders   "
              f"{100.0 * flips[(True, True)] / union_peps:.1f}%"
              f"   ({flips[(True, True)]:,}/{union_peps:,})")
    print(f"    binder peptides 4.1={flips[(True, True)] + flips[(True, False)]:,}"
          f"  4.2={flips[(True, True)] + flips[(False, True)]:,}")
    print(f"\n  wrote {out_dir}/netmhcpan_version_concordance.json")


if __name__ == "__main__":
    main()
