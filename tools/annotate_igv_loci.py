#!/usr/bin/env python
"""
annotate_igv_loci.py — paste-ready IGV coordinates for every candidate.

WHY THIS IS NOT JUST `chrom:pos`
--------------------------------
Three things stand between a table row and actually seeing the reads, and each
one wastes a minute per locus until it is written down.

**The chromosome naming differs between the tables and the BAMs.** These tables
carry `1`; the alignments carry `chr1`. IGV searches the loaded genome, so the
prefix has to match the BAM, not the table. This emits the BAM's form.

**A junction has two ends and one of them is often not on screen.** For an
intra-chromosomal event a single window holds both, padded so the flanks that
make a coverage drop visible are included. For an inter-chromosomal one no single
window can: IGV takes several loci separated by spaces and opens a split view,
which is what `igv_locus` contains in that case.

**Padding has to scale with the lesion.** A fixed window either buries a 100 bp
deletion in empty flank or shows a slice of a 50 kb one. The window here is the
event plus a quarter of its size on each side, with a floor so small events still
show enough flank to judge a depth change, and a ceiling so large ones stay
navigable.

WHAT TO LOOK AT ONCE IT IS OPEN
-------------------------------
For a deletion the decisive observation is in the DNA track, not the RNA one:
depth inside the interval against its flanks, about half for a heterozygous
lesion and unchanged if there is no deletion. The RNA track shows the gapped
reads, but a gap alone cannot distinguish a deleted allele from an alternative
splice junction — both look the same against a reference that still carries the
sequence.

Usage
-----
    python tools/annotate_igv_loci.py --table <dir>/candidate_universe.tsv \\
        --out <dir>/candidate_universe.tsv --chr-prefix chr
"""
from __future__ import annotations

import argparse
import pathlib

import pandas as pd

#: Smallest and largest half-window, in bp. The floor keeps enough flank visible
#: to judge a depth change on a small lesion; the ceiling stops a large one from
#: opening at a zoom where individual reads are not drawn.
MIN_PAD, MAX_PAD = 200, 5000
#: Half-window when looking at one breakend on its own.
BREAKEND_PAD = 150


def fmt(chrom: str, start: int, end: int, prefix: str) -> str:
    name = str(chrom).strip()
    if name.lower().startswith("chr"):
        name = name[3:]
    return f"{prefix}{name}:{max(1, int(start)):,}-{int(end):,}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    parser.add_argument("--chr-prefix", default="chr",
                        help="what the BAM calls its contigs: `chr` for hg38 "
                             "alignments that use chr1, empty for ones that use 1. "
                             "IGV searches the loaded genome, so this must match "
                             "the BAM rather than the table.")
    args = parser.parse_args()

    table = pd.read_csv(args.table, sep="\t", dtype=str, low_memory=False)
    print(f"  {len(table):,} rows")

    loci, bp1s, bp2s, spans = [], [], [], []
    for row in table.itertuples(index=False):
        try:
            c1, p1 = str(row.chrom1), int(float(row.pos1))
            c2, p2 = str(row.chrom2), int(float(row.pos2))
        except (TypeError, ValueError, AttributeError):
            loci.append(""); bp1s.append(""); bp2s.append(""); spans.append("")
            continue

        one = fmt(c1, p1 - BREAKEND_PAD, p1 + BREAKEND_PAD, args.chr_prefix)
        two = fmt(c2, p2 - BREAKEND_PAD, p2 + BREAKEND_PAD, args.chr_prefix)
        bp1s.append(one)
        bp2s.append(two)

        same = str(c1).replace("chr", "") == str(c2).replace("chr", "")
        if same:
            lo, hi = min(p1, p2), max(p1, p2)
            pad = min(MAX_PAD, max(MIN_PAD, (hi - lo) // 4))
            loci.append(fmt(c1, lo - pad, hi + pad, args.chr_prefix))
            spans.append(str(hi - lo))
        else:
            # IGV opens a split view when given several loci separated by spaces;
            # a single window cannot hold an inter-chromosomal junction.
            loci.append(f"{one} {two}")
            spans.append("")

    table["igv_locus"] = loci
    table["igv_bp1"] = bp1s
    table["igv_bp2"] = bp2s

    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False)

    n_split = sum(1 for v in loci if " " in v)
    print(f"  igv_locus: {len(loci) - n_split - loci.count(''):,} single-window, "
          f"{n_split:,} split-view (inter-chromosomal)")
    print(f"  wrote {args.out}")
    for value in list(dict.fromkeys(v for v in loci if v))[:3]:
        print(f"    e.g. {value}")


if __name__ == "__main__":
    main()
