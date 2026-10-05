#!/usr/bin/env python
"""
verify_junction_peptide.py — prove a candidate's junction peptides against the
RNA reads themselves, not against the caller's coordinates.

WHY
---
Every peptide this pipeline emits is computed from a VCF breakpoint plus a
reference transcript. Both can be wrong: a breakpoint can be shifted by
microhomology, and the aligner can place the RNA gap somewhere else again. A
peptide that only exists in the arithmetic is indistinguishable, in every
downstream table, from one the cell actually makes.

This closes that gap. For one event it:

  1. rebuilds the mutant CDS from the transcript's own coding sequence by
     removing the deleted genomic positions, and translates it;
  2. lists the peptides spanning the new junction;
  3. writes the junction NUCLEOTIDE sequence predicted by the breakpoints and
     counts how many RNA reads contain it verbatim (either orientation);
  4. does the same for the junction implied by the N gap the aligner actually
     placed, when that differs.

A peptide supported by reads carrying its own junction sequence is as close to
observed as sequence data gets. A peptide whose junction k-mer appears in zero
reads, while an alternative does, is an artefact of the coordinates.

WHY THE ALIGNER'S GAP IS NOT THE GROUND TRUTH
---------------------------------------------
Spliced aligners score canonical GT-AG motifs, so an N gap is routinely placed a
few bases from the true junction when a nearby motif fits. Comparing the VCF
interval with the BAM gap therefore shows an offset that is not an error in
either. The read SEQUENCE decides between them, which is what step 3-4 use.

Usage
-----
    python tools/verify_junction_peptide.py --rna-bam <bam> --transcript ENST... \\
        --chrom chr1 --pos1 <left breakend> --pos2 <right breakend> \\
        [--release 115] [--flank 20] [--peptide PEPTIDE]
"""
from __future__ import annotations

import argparse
import collections
import sys

from Bio.Seq import Seq

LENGTHS = (8, 9, 10, 11)


def cds_offsets(transcript) -> dict[int, int]:
    """Genomic position -> offset in the transcript's coding sequence."""
    blocks = sorted(transcript.coding_sequence_position_ranges)
    order = blocks if transcript.strand == "+" else blocks[::-1]
    offsets, k = {}, 0
    for start, end in order:
        span = range(start, end + 1) if transcript.strand == "+" else range(end, start - 1, -1)
        for position in span:
            offsets[position] = k
            k += 1
    return offsets


def translate(nt: str) -> str:
    return str(Seq(nt[: len(nt) - len(nt) % 3]).translate(to_stop=True))


def mutant_protein(transcript, offsets, deleted_from: int, deleted_to: int) -> str:
    cds = transcript.coding_sequence
    kept = [p for p in sorted(offsets, key=offsets.get)
            if not deleted_from <= p <= deleted_to]
    return translate("".join(cds[offsets[p]] for p in kept))


def junction_nt(transcript, offsets, last_kept: int, first_kept: int, flank: int) -> str:
    """The flank+flank nucleotides the mutant transcript carries across the seam."""
    cds = transcript.coding_sequence
    left = [p for p in sorted(offsets, key=offsets.get) if p <= last_kept][-flank:] \
        if transcript.strand == "+" else \
        [p for p in sorted(offsets, key=offsets.get) if p >= last_kept][-flank:]
    right = [p for p in sorted(offsets, key=offsets.get) if p >= first_kept][:flank] \
        if transcript.strand == "+" else \
        [p for p in sorted(offsets, key=offsets.get) if p <= first_kept][:flank]
    return "".join(cds[offsets[p]] for p in left + right)


def count_reads(bam_path: str, chrom: str, start: int, end: int, probe: str) -> tuple[int, int]:
    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb")
    contig = chrom if chrom in bam.references else chrom.replace("chr", "")
    other = str(Seq(probe).reverse_complement())
    hits = total = 0
    for read in bam.fetch(contig, max(0, start), end):
        if read.query_sequence is None:
            continue
        total += 1
        if probe in read.query_sequence or other in read.query_sequence:
            hits += 1
    return hits, total


def observed_gaps(bam_path: str, chrom: str, start: int, end: int, min_mapq: int = 20):
    """N-gap populations in the window: (first skipped, last skipped) -> fragments."""
    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb")
    contig = chrom if chrom in bam.references else chrom.replace("chr", "")
    gaps = collections.Counter()
    for read in bam.fetch(contig, max(0, start), end):
        if read.is_unmapped or read.mapping_quality < min_mapq or read.cigartuples is None:
            continue
        position = read.reference_start
        for op, length in read.cigartuples:
            if op in (0, 2, 7, 8):
                position += length
            elif op == 3:
                gaps[(position + 1, position + length)] += 1
                position += length
    return gaps


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--rna-bam", required=True)
    p.add_argument("--transcript", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--pos1", type=int, required=True, help="left breakend, retained")
    p.add_argument("--pos2", type=int, required=True, help="right breakend, retained")
    p.add_argument("--release", type=int, default=115)
    p.add_argument("--flank", type=int, default=20)
    p.add_argument("--peptide", default=None, help="a specific peptide to look for")
    args = p.parse_args()

    import pyensembl
    transcript = pyensembl.EnsemblRelease(args.release).transcript_by_id(args.transcript)
    offsets = cds_offsets(transcript)
    wild = transcript.protein_sequence
    low, high = sorted((args.pos1, args.pos2))
    coding_removed = sum(1 for x in offsets if low < x < high)

    print(f"transcript {transcript.transcript_name} ({transcript.transcript_id}), "
          f"strand {transcript.strand}, protein {len(wild)} aa")
    print(f"deleted genomic interval {low + 1}-{high - 1} ({high - low - 1} bp); "
          f"coding bases removed {coding_removed} "
          f"({'in frame' if coding_removed % 3 == 0 else 'FRAMESHIFT'})")

    mutant = mutant_protein(transcript, offsets, low + 1, high - 1)
    first = next((i for i in range(min(len(wild), len(mutant))) if wild[i] != mutant[i]),
                 min(len(wild), len(mutant)))
    print(f"mutant protein {len(mutant)} aa; first difference at residue {first}")
    print(f"  wild type : …{wild[max(0, first - 10):first + 12]}…")
    print(f"  mutant    : …{mutant[max(0, first - 10):first + 12]}…")

    novel = {mutant[i:i + n] for n in LENGTHS for i in range(len(mutant) - n + 1)} - \
            {wild[i:i + n] for n in LENGTHS for i in range(len(wild) - n + 1)}
    spanning = sorted(x for x in novel
                      if (j := mutant.find(x)) <= first < j + len(x))
    print(f"novel windows {len(novel)}, of which spanning the junction {len(spanning)}")
    for peptide in spanning:
        print(f"    {peptide}")

    probe = junction_nt(transcript, offsets, low, high, args.flank)
    hits, total = count_reads(args.rna_bam, args.chrom, low - 400, high + 400, probe)
    print(f"\njunction sequence predicted by the breakpoints ({args.flank}+{args.flank} nt):")
    print(f"  {probe}")
    print(f"  reads containing it verbatim: {hits} of {total} in the window")

    gaps = observed_gaps(args.rna_bam, args.chrom, low - 400, high + 400)
    print("\nN-gap populations the aligner placed here:")
    for (start, end), n in gaps.most_common(5):
        tag = "  <- matches the breakpoints" if (start, end) == (low + 1, high - 1) else ""
        print(f"  {start}-{end}  len {end - start + 1}  fragments {n}{tag}")
    for (start, end), n in gaps.most_common(5):
        if (start, end) == (low + 1, high - 1) or end - start + 1 != high - low - 1:
            continue
        alt = junction_nt(transcript, offsets, start - 1, end + 1, args.flank)
        alt_hits, _ = count_reads(args.rna_bam, args.chrom, low - 400, high + 400, alt)
        print(f"\nalternative junction implied by the gap at {start}-{end}:")
        print(f"  {alt}")
        print(f"  reads containing it verbatim: {alt_hits}")
        print("  -> the breakpoints are supported; the gap is placed elsewhere by the "
              "aligner" if alt_hits < hits else
              "  -> the GAP is supported and the breakpoints are not: peptides are suspect")

    if args.peptide:
        print(f"\n{args.peptide}: in mutant protein = {args.peptide in mutant}, "
              f"in wild type = {args.peptide in wild}")


if __name__ == "__main__":
    main()
