#!/usr/bin/env python
"""
annotate_dna_vaf.py — what fraction of DNA fragments carry the alteration.

WHY A FRACTION AND NOT THE VARIANT COUNT
----------------------------------------
`vf_bp1` says how many DNA fragments support the junction. On its own that
number is bounded by sequencing depth, which varies between samples and across
the genome, so five fragments at a shallow locus and five at a deep one are
different observations. Dividing by the fragments that could have supported it
removes that dependence:

    dna_vaf_bp1 = vf_bp1 / (vf_bp1 + ref_bp1)

Both counts come from the caller, at the same breakend, in the same unit. That
matters more than it sounds: forming a fraction from two quantities bounded by
different things produces a number that tracks the difference between the
bounds rather than anything biological.

WHAT THE VALUE MEANS, AND THE BIAS IT CARRIES
---------------------------------------------
Roughly 0.5 is one altered copy against one intact, 1.0 is every copy altered,
and a few per cent is the shape of a subclonal population or of misalignment.

**It is biased low, and knowing by how much matters more than the value.** A
fragment can only support the variant by spanning the novel junction with enough
anchored sequence on both sides, while a fragment supports the reference merely
by covering the breakpoint position. The variant allele is therefore harder to
observe than the reference one, and the fraction understates it. The size of the
bias depends on read length, insert size and the aligner, none of which this
column knows.

For a deletion there is a second estimator without that asymmetry: read depth
inside the interval against its flanks, which `tools/validate_junction.py`
measures. It needs no read to span anything, so where the two disagree the depth
ratio is the less biased of the two. They are reported separately rather than
combined.

Usage
-----
    python tools/annotate_dna_vaf.py --table <dir>/candidate_universe.tsv \\
        --inputs <dir>/inputs --out <dir>/candidate_universe.tsv
"""
from __future__ import annotations

import argparse
import pathlib
import sys
from datetime import date

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from provenance import code_version                        # noqa: E402

JOIN_KEY = ["chrom1", "pos1", "chrom2", "pos2"]
#: Carried so the fraction can be recomputed, and so a reader can see the
#: denominator rather than trusting the ratio.
COUNTS = ["ref_bp1", "ref_bp2"]


def junction_table(inputs: pathlib.Path) -> pd.DataFrame:
    frames = []
    for path in sorted(inputs.glob("*/stage6_8_all_junctions.tsv")):
        frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        missing = [c for c in JOIN_KEY + COUNTS if c not in frame.columns]
        if missing:
            raise SystemExit(f"  {path} lacks {missing}")
        frames.append(frame[JOIN_KEY + COUNTS])
        print(f"    {path.parent.name}: {len(frame):,} junctions")
    if not frames:
        raise SystemExit(f"  nothing under {inputs}")
    joined = pd.concat(frames, ignore_index=True)
    return joined.drop_duplicates(subset=JOIN_KEY, keep="first")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", type=pathlib.Path, required=True)
    parser.add_argument("--inputs", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    source = junction_table(args.inputs)
    table = pd.read_csv(args.table, sep="\t", dtype=str, low_memory=False)
    print(f"  table: {len(table):,} rows, {len(table.columns)} columns")

    clash = [c for c in COUNTS if c in table.columns]
    if clash:
        raise SystemExit(f"  already present, refusing to overwrite: {clash}")

    before = len(table)
    table = table.merge(source, on=JOIN_KEY, how="left", validate="many_to_one")
    assert len(table) == before, "merge changed the row count"

    for side in ("1", "2"):
        variant = pd.to_numeric(table[f"vf_bp{side}"], errors="coerce")
        reference = pd.to_numeric(table[f"ref_bp{side}"], errors="coerce")
        total = variant + reference
        # A zero denominator is not a VAF of zero: it means nothing was observed
        # at that breakend, and writing 0 would read as "no support for the
        # variant" when the truth is "no fragments either way".
        table[f"dna_vaf_bp{side}"] = (variant / total.where(total > 0)).round(4)
        table[f"dna_depth_bp{side}"] = total

    vaf1 = pd.to_numeric(table["dna_vaf_bp1"], errors="coerce")
    print(f"\n  dna_vaf_bp1 computable on {vaf1.notna().sum():,} rows")
    for label, value in (("median", vaf1.median()), ("mean", vaf1.mean())):
        print(f"    {label:>7}: {value:.3f}")
    for lo, hi, what in ((0, .1, "< 0.10   subclonal or misalignment-shaped"),
                         (.1, .35, "0.10-0.35"),
                         (.35, .65, "0.35-0.65  one-copy-shaped"),
                         (.65, 1.01, "> 0.65   every-copy-shaped")):
        n = vaf1.between(lo, hi, inclusive="left").sum()
        print(f"    {what:<42} {n:>7,}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False, lineterminator="\n")
    print(f"\n  wrote {args.out} ({len(table.columns)} columns)")

    pd.DataFrame([
        *code_version(),
        {"key": "generated", "value": date.today().isoformat()},
        {"key": "tool", "value": pathlib.Path(__file__).name},
        {"key": "inputs", "value": str(args.inputs)},
        {"key": "formula", "value": "vf_bpN / (vf_bpN + ref_bpN), both from the caller"},
        {"key": "rows_computable", "value": int(vaf1.notna().sum())},
        {"key": "median_dna_vaf_bp1", "value": round(float(vaf1.median()), 4)},
        {"key": "note", "value":
            "Biased LOW: a fragment supports the variant only by spanning the "
            "novel junction with anchored sequence on both sides, while it "
            "supports the reference merely by covering the breakpoint. For a "
            "deletion, depth inside the interval against its flanks "
            "(tools/validate_junction.py) has no such asymmetry and is the less "
            "biased estimator. Empty means no fragments at that breakend, not 0."},
    ]).to_csv(args.out.parent / "DNA_VAF_PROVENANCE.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
