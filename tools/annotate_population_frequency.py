#!/usr/bin/env python
"""
annotate_population_frequency.py — is this ordinary population variation?

WHY TWO COLUMNS AND NOT ONE
---------------------------
The table previously carried a single `is_private` verdict combining a panel of
normals with gnomAD. Two independent sources were folded into one boolean, and
because an absent gnomAD record counted as a pass, gnomAD's silence overrode
what the panel had actually measured. The result: junctions seen in 60% of the
panel of normals were labelled private to the sample.

The sources disagree too often to be summarised together. Across this catalogue
they agree on fewer than half the junctions: some are common in the panel and
absent from gnomAD, others common in gnomAD and rare in the panel. So each gets
its own column, and the reader sees the disagreement instead of inheriting a
resolution of it.

  is_highfreq_gnomad   gnomad_af_popmax >= POPULATION_COMMON_AF
  is_highfreq_panel    pon_fraction     >= POPULATION_COMMON_AF

THE TWO ABSENCES ARE NOT THE SAME, THOUGH BOTH READ AS FALSE
------------------------------------------------------------
The panel is a closed set of genomes screened by this pipeline, so a junction
absent from it was looked for and not found: a measured zero (`PON_ABSENT_MEANS`).
gnomAD is an external database, where absence can also mean the event is not
represented comparably or that the reciprocal-overlap match failed. Both are
written False, because with two separate columns a False here no longer hides
panel evidence — the panel column carries it. **`gnomad_af_popmax` being empty
is what distinguishes "gnomAD says rare" from "gnomAD has no record".**

Nothing is filtered. These are columns.

Usage
-----
    python tools/annotate_population_frequency.py \\
        --table <dir>/candidate_universe.tsv --out <dir>/candidate_universe.tsv
"""
from __future__ import annotations

import argparse
import pathlib
import sys
from datetime import date

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
from svneo import criteria as C                            # noqa: E402

#: Retired from this table. It combined the two sources into one verdict in
#: which gnomAD's silence outvoted the panel's measurement. Dropped rather than
#: renamed: the replacement columns have the OPPOSITE polarity (True now means
#: common, where `is_private` True meant keep), so anything still filtering on
#: the old name must fail loudly instead of silently selecting the complement.
RETIRED = "is_private"


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--threshold", type=float, default=C.POPULATION_COMMON_AF)
    parser.add_argument("--keep-retired", action="store_true",
                        help="leave `is_private` in place instead of dropping it")
    args = parser.parse_args()

    table = pd.read_csv(args.table, sep="\t", dtype=str, low_memory=False)
    print(f"  table: {len(table):,} rows, {len(table.columns)} columns")
    print(f"  threshold: {args.threshold:.2%} "
          f"(~{args.threshold * 11912:.0f} genomes in a panel of ~11,912)")

    popmax = pd.to_numeric(table.get("gnomad_af_popmax"), errors="coerce")
    fraction = pd.to_numeric(table.get("pon_fraction"), errors="coerce")

    # fillna(0) on each, for different reasons stated in the docstring above.
    table["is_highfreq_gnomad"] = (popmax.fillna(0) >= args.threshold)
    table["is_highfreq_panel"] = (fraction.fillna(0) >= args.threshold)

    junctions = table.drop_duplicates(["chrom1", "pos1", "pos2"])
    g = junctions["is_highfreq_gnomad"]
    p = junctions["is_highfreq_panel"]
    print(f"\n  per distinct junction ({len(junctions):,}):")
    print(f"    common in gnomAD              {g.sum():>5}")
    print(f"      of which no gnomAD record   "
          f"{int((~g & pd.to_numeric(junctions['gnomad_af_popmax'], errors='coerce').isna()).sum()):>5} "
          f"are False for lack of a record, not for being rare")
    print(f"    common in the panel           {p.sum():>5}")
    print(f"    common in BOTH                {(g & p).sum():>5}")
    print(f"    common in NEITHER             {(~g & ~p).sum():>5}")
    print(f"    disagree                      {(g ^ p).sum():>5}")

    if RETIRED in table.columns and not args.keep_retired:
        was = table[RETIRED].astype(str).str.lower().eq("true").sum()
        table = table.drop(columns=[RETIRED])
        print(f"\n  dropped `{RETIRED}` ({was:,} rows had it True)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False, lineterminator="\n")
    print(f"  wrote {args.out} ({len(table.columns)} columns)")

    pd.DataFrame([
        {"key": "generated", "value": date.today().isoformat()},
        {"key": "tool", "value": pathlib.Path(__file__).name},
        {"key": "threshold", "value": args.threshold},
        {"key": "threshold_source", "value":
            "POPULATION_COMMON_AF in criteria.py — 1%, the classical "
            "polymorphism/rare-variant line. ACMG/AMP stand-alone benign sits at "
            "5% but answers a different question."},
        {"key": "retired_column", "value": RETIRED},
        {"key": "note", "value":
            "Absent gnomAD record and absent panel count both read as False, but "
            "they are different states: the panel is a closed screened set so "
            "absence is a measured zero, while gnomAD is external and absence can "
            "mean the event is not comparably represented. gnomad_af_popmax being "
            "empty distinguishes them. Panel coverage is near-complete for DEL and "
            "DUP and partial for TRA/INV; no t2tINV junction has ever matched the "
            "panel, so for that type absence is uninformative."},
    ]).to_csv(args.out.parent / "POPULATION_FREQUENCY_PROVENANCE.tsv",
              sep="\t", index=False)


if __name__ == "__main__":
    main()
