#!/usr/bin/env python3
"""
diagnose_minus_strand_cds.py — minus-strand transcripts get the wrong 5' CDS.

WHAT IT CHECKS
--------------
`neosv.fusion_utils.truncate_cds(transcript, '5', pos)` returns what a
transcript contributes as the HEAD of a fusion protein: its coding sequence from
the start codon up to the breakpoint. This sweeps a breakpoint through every
region of a transcript — 5'UTR, each coding exon, each intron, 3'UTR — and
compares the length returned against the length that should be returned.

It verifies the LENGTH. An earlier version checked that the sequence began with
`ATG` and reported the exonic path as correct; that check is worthless, because
the sequence is always taken from position 0 of the coding sequence and so
begins with ATG at any non-zero length. Starting correctly is not the same as
measuring correctly.

WHAT IT FINDS
-------------
Plus-strand transcripts are correct in every region. Minus-strand transcripts are
wrong in almost all of them, in two distinct ways:

  intronic breakpoint   the head comes back EMPTY — the fusion loses its start
                        codon entirely and is flagged `Start-loss`
  exonic breakpoint     the length is counted from the wrong end, so the head is
                        too short near the start of the transcript and TOO LONG
                        near its end

The second is the dangerous one: sequence the gene does not contribute is added
to the fusion, and the sliding window turns it into peptides. This can fabricate
candidates, not merely lose them.

WHY IT HAPPENS
--------------
`transcript.coding_sequence_position_ranges` is ordered by COORDINATE. Reading
order is the same thing on the plus strand and the reverse of it on the minus
strand. `get_cds_range` documents its output as "from 5' to 3'" and
`get_noncds_range` has a minus-strand branch that compensates for a reversal
that never happened, producing intron intervals with `start > end` — which can
never match — and a first interval that spans most of the gene.

    python tools/diagnose_minus_strand_cds.py            # the sweep
    python tools/diagnose_minus_strand_cds.py --explain   # a toy walkthrough
"""
from __future__ import annotations

import argparse
import os
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent

#: Genes with both strands represented, chosen for being well annotated and
#: multi-exon rather than for any biological relevance here.
GENES = ["COL1A1", "TP53", "BRCA1", "NF1", "KRAS", "BRAF",
         "EGFR", "MYC", "PTEN", "PIK3CA"]


#: A three-exon gene, small enough to follow by hand. Coordinates are arbitrary;
#: what matters is that the exons are separated by introns and that the
#: transcript extends a little beyond the outermost coding exons.
TOY_CDS = [(100, 200), (300, 400), (500, 600)]
TOY_SPAN = (50, 650)
TOY_BREAK = 450                     # in the intron between the 2nd and 3rd exon


def toy_gaps(cds, span, strand):
    """The non-coding intervals, by the arithmetic of `get_noncds_range`.

    Reproduced here rather than imported so the walkthrough shows the formula
    beside its result. It must stay identical to
    `vendor/neosv/transcript_utils.py`; the assertions in `explain()` check it
    against the real function on real transcripts.
    """
    start, end = span
    out = []
    if strand == "+":
        out.append(((start, cds[0][0] - 1),
                    "(transcript.start, cds[0].start-1)"))
        for i in range(1, len(cds)):
            out.append(((cds[i - 1][1] + 1, cds[i][0] - 1),
                        f"(cds[{i-1}].end+1, cds[{i}].start-1)"))
        out.append(((cds[-1][1] + 1, end), "(cds[-1].end+1, transcript.end)"))
    else:
        out.append(((cds[0][1] + 1, end), "(cds[0].end+1, transcript.end)"))
        for i in range(1, len(cds)):
            out.append(((cds[i][1] + 1, cds[i - 1][0] - 1),
                        f"(cds[{i}].end+1, cds[{i-1}].start-1)"))
        out.append(((start, cds[-1][0] - 1),
                    "(transcript.start, cds[-1].start-1)"))
    return out


def explain() -> None:
    """Walk one breakpoint through both strands, showing every interval.

    The exon list arrives sorted by coordinate on BOTH strands. On the plus
    strand that is also reading order; on the minus strand it is the reverse.
    The minus branch compensates for a reversal that never happened.
    """
    cds, span = TOY_CDS, TOY_SPAN
    sizes = [e - s + 1 for s, e in cds]
    print(f"A three-exon gene. Coding exons at {cds} — {sizes[0]} nt each, "
          f"{sum(sizes)} nt of CDS in total.")
    print(f"The transcript spans {span[0]}-{span[1]}. The breakpoint is at "
          f"{TOY_BREAK}, in the intron between the 2nd and 3rd exon.\n")
    print("The exon list arrives sorted by COORDINATE on both strands. What "
          "differs is the\nREADING order, and therefore which exons lie "
          "'before' the break.\n")
    print("A break does not delete anything: it SELECTS what this transcript "
          "contributes as\nthe 5' head of the fusion. The rest is replaced by "
          "the partner's 3' tail.\n")

    for strand in ("+", "-"):
        # Reading order: coordinate order on the plus strand, reversed on the
        # minus. The exons preceding the break are those read before it.
        order = list(range(len(cds))) if strand == "+" else \
            list(reversed(range(len(cds))))
        names = "ABC"
        reading = " -> ".join(names[i] for i in order)
        before = [i for i in order
                  if (cds[i][1] < TOY_BREAK if strand == "+"
                      else cds[i][0] > TOY_BREAK)]
        want_nt = sum(cds[i][1] - cds[i][0] + 1 for i in before)

        print("=" * 72)
        print(f"STRAND {strand}   read {reading}")
        print(f"  exons before the break, in reading order: "
              f"{', '.join(names[i] for i in before) or 'none'}")
        print(f"  so the 5' head SHOULD be {len(before)} exon(s) = {want_nt} nt")
        print("=" * 72)

        gaps = toy_gaps(cds, span, strand)
        print("  the non-coding intervals it builds to locate the break:")
        for index, (interval, formula) in enumerate(gaps):
            note = "  <-- INVERTED: start > end, nothing can fall inside" \
                if interval[0] > interval[1] else ""
            print(f"    gap {index}  {formula:38s} = {interval}{note}")

        match = next((i for i, (iv, _) in enumerate(gaps)
                      if iv[0] <= TOY_BREAK <= iv[1]), None)
        kept = cds[:match] if match is not None else []
        got_nt = sum(e - s + 1 for s, e in kept)
        verdict = "correct" if got_nt == want_nt else \
            f"WRONG — {want_nt} nt expected, the fusion loses its head"
        print(f"\n  first interval containing {TOY_BREAK}: gap {match}")
        print(f"  gap index {match} means '{match} coding exon(s) precede it'")
        print(f"  the 5' head becomes cds[:{match}] = {len(kept)} exon(s) = "
              f"{got_nt} nt   [{verdict}]\n")

    print("On the minus strand the answer should be C alone. Two faults combine "
          "to give\nnothing instead. Gap 0 is measured from the end of the "
          "first exon in the LIST (A)\nrather than the first in reading order "
          "(C), so it collapses leftwards across most\nof the gene; and every "
          "intron comes out with start > end, so none can match.\nThe break "
          "lands in gap 0, which reads as 'before the first exon'.\n")


def expected_head(cds_ranges, strand: str, position: int) -> int:
    """Nucleotides the 5' head SHOULD contain: every coding exon read before the
    breakpoint, plus the part of the exon the breakpoint falls in.

    "Before" means before in READING order — lower coordinates on the plus
    strand, higher on the minus. This is the definition the code is measured
    against, and it is deliberately written from the biology rather than from
    the implementation.
    """
    total = 0
    for start, end in cds_ranges:
        if strand == "+":
            if end < position:
                total += end - start + 1
            elif start <= position <= end:
                total += position - start + 1
        else:
            if start > position:
                total += end - start + 1
            elif start <= position <= end:
                total += end - position + 1
    return total


def sweep(transcript, strand: str) -> list[dict]:
    """One breakpoint per region, in reading order, with expected vs returned."""
    from neosv.fusion_utils import truncate_cds
    ranges = sorted(transcript.coding_sequence_position_ranges)
    reading = ranges if strand == "+" else list(reversed(ranges))

    probes = []
    for index, (start, end) in enumerate(reading):
        probes.append((f"exon {index + 1}", (start + end) // 2))
        if index + 1 < len(reading):
            nxt = reading[index + 1]
            probes.append((f"intron {index + 1}",
                           (end + nxt[0]) // 2 if strand == "+"
                           else (nxt[1] + start) // 2))

    rows = []
    for label, position in probes:
        collection = truncate_cds(transcript, "5", position)
        returned = collection.cut_length if collection else 0
        wanted = expected_head(ranges, strand, position)
        rows.append({"region": label, "position": position,
                     "expected": wanted, "returned": returned,
                     "ok": returned == wanted})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release", type=int, default=115)
    parser.add_argument("--explain", action="store_true",
                        help="walk one breakpoint through both strands on a "
                             "three-exon toy gene, showing every interval")
    parser.add_argument("--detail", action="store_true",
                        help="print every region, not just the summary")
    args = parser.parse_args()

    if args.explain:
        explain()

    sys.path.insert(0, str(REPO / "vendor"))
    from pyensembl import EnsemblRelease
    genome = EnsemblRelease(args.release)

    print("5' head length returned by truncate_cds, swept across every region\n")
    print(f"  {'gene':9s} {'strand':7s} {'regions':>8s} {'correct':>8s} "
          f"{'empty':>6s} {'short':>6s} {'long':>6s}")
    print("  " + "-" * 54)

    tested, failing = [], []
    for name in GENES:
        try:
            gene = genome.genes_by_name(name)[0]
        except Exception:                                   # noqa: BLE001
            continue
        transcript = next((t for t in gene.transcripts
                           if t.is_protein_coding and t.complete), None)
        if transcript is None or \
                len(transcript.coding_sequence_position_ranges) < 3:
            continue
        rows = sweep(transcript, gene.strand)
        tested.append(name)
        good = sum(r["ok"] for r in rows)
        empty = sum(1 for r in rows if not r["ok"] and r["returned"] == 0)
        short = sum(1 for r in rows if not r["ok"] < 0 and
                    0 < r["returned"] < r["expected"])
        long_ = sum(1 for r in rows if r["returned"] > r["expected"])
        if good < len(rows):
            failing.append(name)
        print(f"  {name:9s} {gene.strand:7s} {len(rows):8d} {good:8d} "
              f"{empty:6d} {short:6d} {long_:6d}")
        if args.detail:
            for r in rows:
                mark = "ok" if r["ok"] else f"{r['returned'] - r['expected']:+d}"
                print(f"      {r['region']:12s} {r['position']:12,} "
                      f"expected {r['expected']:6,}  returned {r['returned']:6,}"
                      f"   {mark}")

    print()
    if not tested:
        raise SystemExit(
            "no transcripts were tested — nothing was verified.\n"
            "Set PYENSEMBL_CACHE_DIR to the annotation cache, e.g.\n"
            "  export PYENSEMBL_CACHE_DIR=/path/to/pyensembl_cache")

    if failing:
        print(f"{len(failing)} of {len(tested)} transcripts return a wrong 5' "
              f"head in at least one region: {', '.join(failing)}")
        print("\nRe-run with --detail to see every region, or --explain for a "
              "worked\nthree-exon example of why.")
        sys.exit(1)
    print(f"All {len(tested)} transcripts return the correct 5' head in every "
          f"region.")


if __name__ == "__main__":
    main()
