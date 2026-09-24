#!/usr/bin/env python
"""
annotate_junction_usage.py — the junction count, with the denominator it needs.

WHY A RAW COUNT IS NOT ENOUGH, AND WHY IT IS STILL THE HEADLINE
--------------------------------------------------------------
`junction_reads` is the evidence that the altered transcript exists, and it
stays the quotable number: it is a direct observation, it needs no model, and
for antigen presentation the absolute amount of altered transcript is closer to
what matters than any ratio.

What it cannot do is be compared. Seven thousand junction reads is overwhelming
support in a gene at 5 TPM and unremarkable in one at 900 TPM, so a threshold on
the raw count is in practice a filter on expression: it drops real events in
quiet genes and admits alignment noise in loud ones. Ranking candidates by it
reproduces the ranking of gene expression, which the TPM column already gives.

So this adds the denominator rather than replacing the count.

THE DENOMINATOR, WHICH IS THE ENTIRE DIFFICULTY
-----------------------------------------------
Three populations could be used and two of them are wrong:

  - **every read at the locus** — wrong. A read lying entirely inside an exon
    never had the opportunity to use either junction. Counting reads that could
    not vote makes the share a function of exon length, not of biology.
  - **every gapped read at the locus** — wrong. A neighbouring intron spliced
    out of every transcript regardless does not compete with the candidate; it
    only dilutes. On a worked example this was the difference between 10% and 28%.
  - **gapped reads whose gap overlaps the SV interval** — right. Those are the
    reads that could have gone either way.

Hence three columns rather than one, because the distinction between them is
exactly what gets misread:

    junction_reads           cross the interval BY the candidate junction
    junction_alt_reads       cross it WITHOUT the candidate junction
    junction_interval_total  the sum — and the denominator, because a share
                             requires the numerator to be inside it
    junction_usage           junction_reads / junction_interval_total

`junction_alt_reads` deliberately EXCLUDES the junction's own reads, so it is
the competing outcome on its own. It is reported for reading, never for
dividing: `junction_reads / junction_alt_reads` is odds, not a percentage, and
is not bounded by 100%.

WHERE IT DOES NOT APPLY
-----------------------
Only a size gap has a canonical alternative competing over the same interval. A
chimeric junction has none, and an insertion is not a choice between two splice
outcomes. Those rows are left EMPTY, never 0 — a 0 would read as "the junction
is never used", the opposite of "the question does not apply here".

Nothing is filtered. These are columns.

Usage
-----
    python tools/annotate_junction_usage.py \\
        --table <dir>/candidate_universe.tsv \\
        --config config/<run>.yaml \\
        --out   <dir>/candidate_universe.tsv
"""
from __future__ import annotations

import argparse
import pathlib
import sys
from datetime import date

import pandas as pd
import yaml

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from validate_junction import junction_usage, MIN_MAPQ     # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
from svneo import criteria as C                            # noqa: E402

#: The only geometry with a competing canonical junction over the same interval.
APPLICABLE_TEST = "sizegap"


def bam_by_sample(config_path: pathlib.Path) -> dict[str, str]:
    """sample name -> RNA BAM, straight from the run config.

    Read from the config rather than passed on the command line so the BAM a
    row is measured against is always the one its calls came from. A row
    measured against another line's BAM would be silently meaningless.
    """
    config = yaml.safe_load(config_path.read_text())
    mapping = {}
    for sample in config.get("samples", []) or []:
        if sample.get("rna_bam"):
            mapping[sample["name"]] = sample["rna_bam"]
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", type=pathlib.Path, required=True)
    parser.add_argument("--config", type=pathlib.Path, required=True,
                        help="the run config, for the per-sample RNA BAM paths")
    parser.add_argument("--sample-column", default="acquired_in",
                        help="which column names the line a row was called in")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    bams = bam_by_sample(args.config)
    print(f"  RNA BAMs for {len(bams)} samples: {', '.join(sorted(bams))}")
    missing = [name for name, path in bams.items() if not pathlib.Path(path).exists()]
    if missing:
        # Refuse rather than silently annotate a subset: a column that is empty
        # because a BAM was absent is indistinguishable from one empty because
        # the test does not apply, and that is the misreading this file exists
        # to prevent.
        raise SystemExit(f"  RNA BAM not found for: {', '.join(missing)} — "
                         f"refusing to write a column that would be silently partial")

    table = pd.read_csv(args.table, sep="\t", dtype=str, low_memory=False)
    print(f"  table: {len(table):,} rows")

    applicable = table["test"] == APPLICABLE_TEST
    print(f"  {applicable.sum():,} rows with test == {APPLICABLE_TEST!r}; "
          f"{(~applicable).sum():,} left empty (not applicable)")

    # One measurement per distinct junction per sample, not per peptide: a locus
    # appears in thousands of rows and the answer cannot differ between them.
    keys = list(zip(table[args.sample_column], table["chrom1"],
                    table["pos1"], table["pos2"], table["test"]))
    distinct = {k for k, ok in zip(keys, applicable) if ok}
    print(f"  {len(distinct):,} distinct junctions to measure")

    cache: dict[tuple, dict] = {}
    for n, key in enumerate(sorted(distinct), 1):
        sample, chrom, pos1, pos2, _ = key
        bam = bams.get(sample)
        if not bam:
            continue
        try:
            cache[key] = junction_usage(bam, str(chrom), int(float(pos1)),
                                        int(float(pos2)))
        except (ValueError, OSError) as error:
            print(f"    {chrom}:{pos1}-{pos2} in {sample}: {error}")
        if n % 50 == 0:
            print(f"    {n:,}/{len(distinct):,}")

    def column(field, fmt=str):
        return [fmt(cache[k][field]) if k in cache else "" for k in keys]

    table["junction_alt_reads"] = column("alt_reads")
    table["junction_interval_total"] = column("total")
    table["junction_usage"] = [
        f"{cache[k]['usage']:.4f}" if k in cache and cache[k]["total"] else ""
        for k in keys]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False, lineterminator="\n")

    measured = table["junction_usage"] != ""
    print(f"\n  wrote {args.out}")
    print(f"  junction_usage present on {measured.sum():,} rows")
    if measured.any():
        values = pd.to_numeric(table.loc[measured, "junction_usage"])
        for label, value in (("min", values.min()), ("median", values.median()),
                             ("max", values.max())):
            print(f"    {label:>7}: {value:.1%}")
        print(f"    below 5% (noise-shaped): {(values < 0.05).sum():,}")
        print(f"    0.35-0.65 (heterozygous-shaped): "
              f"{values.between(0.35, 0.65).sum():,}")

    pd.DataFrame([
        {"key": "generated", "value": date.today().isoformat()},
        {"key": "tool", "value": pathlib.Path(__file__).name},
        {"key": "config", "value": str(args.config)},
        {"key": "sample_column", "value": args.sample_column},
        {"key": "applicable_test", "value": APPLICABLE_TEST},
        {"key": "min_mapq", "value": MIN_MAPQ},
        {"key": "ngap_size_tolerance_bp", "value": C.NGAP_SIZE_TOLERANCE},
        {"key": "ngap_position_tolerance_bp", "value": C.NGAP_POSITION_TOLERANCE},
        {"key": "junctions_measured", "value": len(cache)},
        {"key": "denominator", "value":
            "gapped fragments whose gap overlaps the SV interval "
            "(gap_start < pos2 and gap_end > pos1), the candidate junction "
            "INCLUDED. Ungapped reads are excluded: they could not use either "
            "junction. Gaps elsewhere in the gene are excluded: they do not compete."},
        {"key": "note", "value":
            "junction_alt_reads EXCLUDES junction_reads and is NOT a denominator; "
            "dividing by it gives odds, not a share. Empty means not applicable "
            "(test != sizegap), never zero."},
    ]).to_csv(args.out.parent / "JUNCTION_USAGE_PROVENANCE.tsv",
              sep="\t", index=False)


if __name__ == "__main__":
    main()
