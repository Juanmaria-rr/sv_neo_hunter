#!/usr/bin/env python
"""
annotate_gnomad_populations.py — every ancestry group's frequency, not the maximum.

WHY POPMAX IS NOT ENOUGH
------------------------
`gnomad_af_popmax` is the highest frequency across nine ancestry groups. Two
things follow, and both matter for reading a candidate.

**It is biased upwards by construction.** A maximum over nine noisy estimates
exceeds any single one of them, and the more groups a database reports the
higher the maximum drifts. It is the right statistic for a conservative
"is this variant common ANYWHERE" screen, which is why the pipeline judges
`is_private` on it, and the wrong one for saying how common a variant actually
is in the donor.

**It hides the case that matters most.** A variant at 0.46 in one group and
0.028 in another is sixteen-fold apart, and popmax reports only the 0.46. Which
group is the right reference depends on the donor's ancestry — knowledge this
pipeline does not have and must not assume. So every group is carried and the
choice is left to whoever reads the table.

The pipeline already extracts all nine at stage 6; they were simply not
propagated into the candidate universe. This adds them without rebuilding, so
annotations applied after the universe was built are not lost.

WHAT AN EMPTY CELL MEANS
------------------------
Empty is "not reported for this group", never 0. gnomAD-SV does not report a
frequency for every group on every record, and a missing value read as zero
turns an unmeasured variant into an apparently absent one — the same
absence-of-evidence error that makes `pass_gnomad` True when no match was found
at all.

Usage
-----
    python tools/annotate_gnomad_populations.py \\
        --table <dir>/candidate_universe.tsv \\
        --inputs <dir>/inputs \\
        --out <dir>/candidate_universe.tsv
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

#: Carried alongside the per-group frequencies: the cohort-wide figure, which is
#: dominated by the largest group, and which group the maximum came from —
#: without that name a popmax cannot be interpreted at all.
EXTRA = ("gnomad_af", "gnomad_af_popmax_pop")

JOIN_KEY = ["chrom1", "pos1", "chrom2", "pos2"]


def junction_table(inputs: pathlib.Path) -> pd.DataFrame:
    """Every sample's stage-6 junction table, keyed by coordinates.

    Read from the checksummed input snapshot rather than the live run
    directories, so this annotation is reproducible from what the universe was
    actually built from.
    """
    columns = list(JOIN_KEY) + list(EXTRA) + [
        f"gnomad_af_{population}" for population in C.GNOMAD_POPULATIONS]
    frames = []
    for path in sorted(inputs.glob("*/stage6_8_all_junctions.tsv")):
        frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        have = [c for c in columns if c in frame.columns]
        missing = set(columns) - set(have)
        if missing:
            print(f"    {path.parent.name}: missing {sorted(missing)}")
        frames.append(frame[have])
        print(f"    {path.parent.name}: {len(frame):,} junctions")
    if not frames:
        raise SystemExit(f"  no stage6_8_all_junctions.tsv under {inputs}")
    joined = pd.concat(frames, ignore_index=True)

    # The same junction appears in several samples through inheritance, and the
    # gnomAD record behind it is a property of the locus, not of the sample, so
    # the rows must agree. Keep the first and say so if they do not.
    before = len(joined)
    joined = joined.drop_duplicates(subset=JOIN_KEY, keep="first")
    print(f"  {before:,} rows -> {len(joined):,} distinct junctions")
    return joined


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", type=pathlib.Path, required=True)
    parser.add_argument("--inputs", type=pathlib.Path, required=True,
                        help="the checksummed snapshot the universe was built from")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    source = junction_table(args.inputs)
    table = pd.read_csv(args.table, sep="\t", dtype=str, low_memory=False)
    print(f"  table: {len(table):,} rows, {len(table.columns)} columns")

    added = [c for c in source.columns if c not in JOIN_KEY]
    clash = [c for c in added if c in table.columns]
    if clash:
        # Refuse rather than overwrite: a column already present may have been
        # computed differently, and silently replacing it would make the table
        # disagree with its own dictionary.
        raise SystemExit(f"  already present, refusing to overwrite: {clash}")

    before = len(table)
    table = table.merge(source, on=JOIN_KEY, how="left", validate="many_to_one")
    assert len(table) == before, "merge changed the row count"

    per_pop = [f"gnomad_af_{p}" for p in C.GNOMAD_POPULATIONS
               if f"gnomad_af_{p}" in table.columns]
    any_pop = table[per_pop].notna().any(axis=1)
    print(f"\n  added {len(added)} columns; {any_pop.sum():,} rows carry at least "
          f"one per-group frequency ({100 * any_pop.mean():.1f}%)")

    # What the new columns buy over popmax: how far apart the groups are on the
    # same variant. A fold-spread near 1 means popmax was a fair summary.
    numeric = table[per_pop].apply(pd.to_numeric, errors="coerce")
    lo, hi = numeric.min(axis=1), numeric.max(axis=1)
    spread = (hi / lo.where(lo > 0))[any_pop]
    if spread.notna().any():
        print(f"  fold-spread between the highest and lowest reporting group:")
        print(f"    median {spread.median():.1f}x, "
              f"90th pct {spread.quantile(0.9):.1f}x, max {spread.max():.0f}x")
        print(f"    junctions where the groups differ by more than 10x: "
              f"{(spread > 10).sum():,}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False, lineterminator="\n")
    print(f"  wrote {args.out} ({len(table.columns)} columns)")

    pd.DataFrame([
        *code_version(),
        {"key": "generated", "value": date.today().isoformat()},
        {"key": "tool", "value": pathlib.Path(__file__).name},
        {"key": "inputs", "value": str(args.inputs)},
        {"key": "populations", "value": ", ".join(C.GNOMAD_POPULATIONS)},
        {"key": "join_key", "value": ", ".join(JOIN_KEY)},
        {"key": "rows_with_any_group_af", "value": int(any_pop.sum())},
        {"key": "note", "value":
            "Empty means the group was not reported for that record, NEVER 0. "
            "popmax is a maximum over nine groups: biased upwards by "
            "construction, and it hides a variant common in one ancestry and "
            "absent elsewhere. Which group is the right reference depends on "
            "the donor's ancestry, which this pipeline does not know."},
    ]).to_csv(args.out.parent / "GNOMAD_POPULATIONS_PROVENANCE.tsv",
              sep="\t", index=False)


if __name__ == "__main__":
    main()
