#!/usr/bin/env python3
"""
inspect_insertion.py — look at the reads behind an insertion call, one by one.

WHY
---
An insertion-driven event is supported by reads that carry the inserted sequence.
A count alone cannot distinguish real support from four familiar artefacts, so
this prints the evidence rather than summarising it:

    PCR/optical duplicates   n copies of one original molecule read as n reads
    low-complexity match     a GC-rich or repetitive insert occurs by chance
    mismapping               reads from a paralogous locus carrying the sequence
    one-sided pile-up        all support at one offset = a single misaligned stack

For each supporting read it reports duplicate flag, mapping quality, start
position, CIGAR, and where inside the read the insertion falls. Reads that are
genuinely independent differ in start position and in the offset of the insert;
duplicates and stacks do not.

A background rate is measured over random loci: how often does this sequence turn
up in reads that have nothing to do with the event? For a low-complexity insert
that number is not zero, and the comparison is the whole point.

USAGE
-----
    python tools/inspect_insertion.py --bam sample.rna.bam \\
        --chrom 12 --pos 132783513 --sequence CGCC...GGG [--window 200] \\
        [--background-loci 200]

    # the same insertion in the DNA BAM: a germline insertion must be there too
    python tools/inspect_insertion.py --bam sample.dna.bam --chrom 12 ...
"""
from __future__ import annotations

import argparse
import collections
import random
import re


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgtN", "TGCAtgcaN"))[::-1]


def complexity(seq: str) -> dict:
    """Composition summary — the first thing to check on any short sequence."""
    import math
    counts = collections.Counter(seq.upper())
    n = len(seq)
    entropy = -sum((c / n) * math.log2(c / n) for c in counts.values())
    longest = best = 1
    for i in range(1, n):
        best = best + 1 if seq[i] == seq[i - 1] else 1
        longest = max(longest, best)
    return {"length": n,
            "gc": round((counts["G"] + counts["C"]) / n, 3),
            "shannon_bits": round(entropy, 2),
            "max_homopolymer": longest,
            "composition": dict(counts)}


def supporting_reads(bam, chrom: str, pos: int, sequence: str, window: int):
    """Every read in the window whose sequence contains the insert (either strand)."""
    contig = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
    reverse = revcomp(sequence)
    found, total = [], 0
    for read in bam.fetch(contig, max(0, pos - window), pos + window):
        if read.is_unmapped or not read.query_sequence:
            continue
        total += 1
        seq = read.query_sequence
        offset = seq.find(sequence)
        strand = "+"
        if offset == -1:
            offset = seq.find(reverse)
            strand = "-" if offset != -1 else "+"
        if offset == -1:
            continue
        found.append({
            "name": read.query_name,
            "start": read.reference_start,
            "mapq": read.mapping_quality,
            "duplicate": read.is_duplicate,
            "secondary": read.is_secondary or read.is_supplementary,
            "cigar": read.cigarstring,
            "insert_offset_in_read": offset,
            "strand_of_match": strand,
            "read_len": len(seq),
        })
    return found, total


def background_rate(bam, sequence: str, loci: int, window: int, seed: int = 0):
    """How often the sequence appears in reads at unrelated loci.

    Without this a raw count is uninterpretable: a low-complexity insert will
    appear in reads anywhere, and the question is whether the event locus is
    enriched over that baseline.
    """
    rng = random.Random(seed)
    reverse = revcomp(sequence)
    main = [(ref, length) for ref, length in zip(bam.references, bam.lengths)
            if re.fullmatch(r"chr\d+", ref) and length > 2_000_000]
    hits = reads = sampled = 0
    for _ in range(loci):
        ref, length = rng.choice(main)
        start = rng.randrange(1_000_000, length - 1_000_000)
        try:
            for read in bam.fetch(ref, start, start + 2 * window):
                if read.is_unmapped or not read.query_sequence:
                    continue
                reads += 1
                if sequence in read.query_sequence or reverse in read.query_sequence:
                    hits += 1
            sampled += 1
        except (ValueError, KeyError):
            continue
    return {"loci_sampled": sampled, "reads_examined": reads,
            "reads_containing_sequence": hits,
            "rate_per_1000_reads": round(1000 * hits / reads, 4) if reads else None}


def main() -> None:
    import pysam
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--pos", type=int, required=True)
    ap.add_argument("--sequence", required=True)
    ap.add_argument("--window", type=int, default=200)
    ap.add_argument("--background-loci", type=int, default=200)
    args = ap.parse_args()

    print(f"insertion: {args.sequence}")
    for key, value in complexity(args.sequence).items():
        print(f"  {key:18} {value}")

    bam = pysam.AlignmentFile(args.bam, "rb")
    reads, total = supporting_reads(bam, args.chrom, args.pos, args.sequence,
                                    args.window)

    print(f"\nlocus {args.chrom}:{args.pos} +/- {args.window} bp — "
          f"{total} reads, {len(reads)} carry the insertion\n")
    if reads:
        print(f"  {'read':32} {'start':>10} {'mapq':>5} {'dup':>4} {'off':>4} "
              f"{'str':>4}  cigar")
        for r in sorted(reads, key=lambda x: x["start"]):
            print(f"  {r['name'][:32]:32} {r['start']:>10} {r['mapq']:>5} "
                  f"{str(r['duplicate'])[0]:>4} {r['insert_offset_in_read']:>4} "
                  f"{r['strand_of_match']:>4}  {r['cigar']}")

        starts = {r["start"] for r in reads}
        offsets = {r["insert_offset_in_read"] for r in reads}
        duplicates = sum(1 for r in reads if r["duplicate"])
        low_mapq = sum(1 for r in reads if r["mapq"] < 20)
        print(f"\n  distinct start positions : {len(starts)} of {len(reads)}")
        print(f"  distinct offsets in read : {len(offsets)} of {len(reads)}")
        print(f"  flagged duplicates       : {duplicates}")
        print(f"  MAPQ < 20                : {low_mapq}")
        print("\n  Independent support looks like many distinct starts AND offsets.")
        print("  One start or one offset repeated = a duplicate stack or a single")
        print("  misaligned pile, not n independent observations.")

    print(f"\nbackground over {args.background_loci} random loci:")
    for key, value in background_rate(bam, args.sequence, args.background_loci,
                                      args.window).items():
        print(f"  {key:28} {value}")
    bam.close()


if __name__ == "__main__":
    main()
