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

WHY THERE IS NO PANEL EQUIVALENT
--------------------------------
There was one, computed from `pon_fraction`, and it has been withdrawn along
with that column. `pon_fraction` was `PON_COUNT` divided by the largest
`PON_COUNT` in the same VCF, because the true size of the panel is not
published: neither the hmftools resource documentation nor the VCF header
states how many samples it was built from, and the field is defined only as
"PON count if in PON".

Two things follow. The denominator is per-sample — it depends on how many
junctions that sample happened to call — so the fraction was not comparable
between lines: the same count gave 8.4% in one and 27.6% in another. And the
largest count observed, 11,912, exceeds every published size of the cohort the
panel is built from (2,399 to ~8,000 patients), so `PON_COUNT` is very unlikely
to be a count of individuals at all.

A number that is neither a frequency nor comparable between samples cannot
support a threshold, so `PON_COUNT` is carried as the raw count it is and the
population-frequency verdict rests on gnomAD, which publishes allele
frequencies with known denominators per ancestry group.

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

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from provenance import code_version                        # noqa: E402
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
from svneo import criteria as C                            # noqa: E402

#: Retired from this table. It combined the two sources into one verdict in
#: which gnomAD's silence outvoted the panel's measurement. Dropped rather than
#: renamed: the replacement columns have the OPPOSITE polarity (True now means
#: common, where `is_private` True meant keep), so anything still filtering on
#: the old name must fail loudly instead of silently selecting the complement.
RETIRED = ["is_private", "pon_fraction", "is_highfreq_panel"]


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
    print(f"  threshold: {args.threshold:.2%} of gnomAD popmax")

    popmax = pd.to_numeric(table.get("gnomad_af_popmax"), errors="coerce")
    table["is_highfreq_gnomad"] = (popmax.fillna(0) >= args.threshold)

    junctions = table.drop_duplicates(["chrom1", "pos1", "pos2"])
    g = junctions["is_highfreq_gnomad"]
    no_record = pd.to_numeric(junctions["gnomad_af_popmax"], errors="coerce").isna()
    print(f"\n  per distinct junction ({len(junctions):,}):")
    print(f"    common in gnomAD              {g.sum():>5}")
    print(f"    not common                    {(~g & ~no_record).sum():>5}")
    print(f"    no gnomAD record at all       {no_record.sum():>5}  "
          f"(False for lack of a record, not for being rare)")

    if not args.keep_retired:
        for column in RETIRED:
            if column in table.columns:
                table = table.drop(columns=[column])
                print(f"  dropped `{column}`")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False, lineterminator="\n")
    print(f"  wrote {args.out} ({len(table.columns)} columns)")

    pd.DataFrame([
        *code_version(),
        {"key": "generated", "value": date.today().isoformat()},
        {"key": "tool", "value": pathlib.Path(__file__).name},
        {"key": "threshold", "value": args.threshold},
        {"key": "threshold_source", "value":
            "POPULATION_COMMON_AF in criteria.py — 1%, the classical "
            "polymorphism/rare-variant line. ACMG/AMP stand-alone benign sits at "
            "5% but answers a different question."},
        {"key": "retired_columns", "value": ", ".join(RETIRED)},
        {"key": "note", "value":
            "False covers both 'gnomAD says rare' and 'gnomAD has no record'; "
            "gnomad_af_popmax being empty distinguishes them. pon_fraction and its "
            "verdict were withdrawn: its denominator is the largest PON_COUNT in "
            "the same VCF, which is per-sample and not the panel size, and the "
            "largest count observed exceeds every published size of the cohort the "
            "panel is built from. PON_COUNT is carried as a raw count."},
    ]).to_csv(args.out.parent / "POPULATION_FREQUENCY_PROVENANCE.tsv",
              sep="\t", index=False)


if __name__ == "__main__":
    main()
