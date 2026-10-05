#!/usr/bin/env python
"""
audit_peptide_generation.py — rebuild each same-gene deletion/duplication's
mutant protein with an independent model and compare its peptides with stage 2.

WHY
---
Stage 2 (vendored NeoSV) has already produced three silent errors that every
downstream test passed: junction-inserted bases placed in the coding sequence at
intronic breakends, non-AUG start codons translated literally, and transcript
choice by genomic span. Each was found only by comparing the output with an
expectation. This makes that comparison systematic: a second, deliberately
simple implementation of the same biology, run on the same events and the same
transcript, so any disagreement is a bug in one of the two.

THE INDEPENDENT MODEL (spliced coding sequence, genomic coordinates)
-------------------------------------------------------------------
Breakends are read from the admitted VCF. For a DEL joining pos1 to pos2 the
deleted bases are pos1+1 .. pos2-1; for a tandem DUP the duplicated block is
pos1 .. pos2. Each event is classified by where its breakends fall in the
transcript the generator used, and only the cases with an unambiguous spliced
product are rebuilt:

  same intron            mature mRNA unchanged -> no peptides expected
  both in coding exons   CDS bases outside the deletion, joined, with the
                         junction-inserted bases (they are exonic here)
  intron / intron        whole exons removed (DEL) or repeated (DUP) between the
  (different introns)    two introns -> exon skipping / exon duplication
  anything else          one breakend exonic and one intronic, UTR, or not a
                         same-gene DEL/DUP: splice outcome not determined by the
                         sequence alone -> reported as `not_modelled`

The rebuilt CDS (plus 3'UTR, for stop-loss read-through) is translated with a
methionine at a retained native start, and every 8-11-mer absent from the
wild-type protein is a peptide. The comparison is per event: identical, or the
peptides only one implementation produced.

Usage
-----
    python tools/audit_peptide_generation.py --run-dir <cross>/<sample>_<branch> \\
        [--run-dir ...] --release 115 --out audit.tsv
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys

import pandas as pd
from Bio.Seq import Seq

LENGTHS = (8, 9, 10, 11)


def windows(protein: str) -> set[str]:
    return {protein[i:i + n] for n in LENGTHS for i in range(len(protein) - n + 1)}


def read_vcf_events(vcf: pathlib.Path) -> dict[str, dict]:
    """sv_id -> {chrom, pos1, pos2, insertion (forward strand)} for the record
    whose mate lies at the higher coordinate on the same chromosome."""
    out = {}
    for line in open(vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        chrom, pos, sv_id, ref, alt = f[0], int(f[1]), f[2], f[3], f[4]
        m = re.search(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]", alt)
        if not m or m.group(1) != chrom or int(m.group(2)) <= pos:
            continue
        bases = re.sub(r"[\[\]][^\[\]]*[\[\]]", "", alt)
        insertion = bases[len(ref):] if alt.startswith(ref) else bases[:-len(ref)]
        out[sv_id] = {"chrom": chrom, "pos1": pos, "pos2": int(m.group(2)),
                      "insertion": insertion}
    return out


def cds_blocks(t) -> list[tuple[int, int]]:
    """Coding blocks in transcript (5'->3') order, genomic, inclusive."""
    blocks = sorted(t.coding_sequence_position_ranges)
    return blocks[::-1] if t.strand == "-" else blocks


def genomic_seq(t, blocks) -> str:
    """Spliced sequence of the given genomic blocks, in transcript orientation,
    taken from the transcript's own cDNA so no genome FASTA is needed."""
    # Map each genomic CDS base to its offset in coding_sequence.
    cds = t.coding_sequence
    offsets, k = {}, 0
    for s, e in cds_blocks(t):
        rng = range(s, e + 1) if t.strand == "+" else range(e, s - 1, -1)
        for g in rng:
            offsets[g] = k
            k += 1
    out = []
    for s, e in blocks:
        rng = range(s, e + 1) if t.strand == "+" else range(e, s - 1, -1)
        out.append("".join(cds[offsets[g]] for g in rng))
    return "".join(out)


def region(t, pos: int) -> tuple[str, int | None]:
    blocks = sorted(t.coding_sequence_position_ranges)
    for s, e in blocks:
        if s <= pos <= e:
            return "CDS", None
    exons = sorted((x.start, x.end) for x in t.exons)
    for s, e in exons:
        if s <= pos <= e:
            return "UTR", None
    for i in range(len(exons) - 1):
        if exons[i][1] < pos < exons[i + 1][0]:
            return "intron", i
    return "outside", None


def translate(nt: str, native_start: bool) -> str:
    nt = nt[: len(nt) - len(nt) % 3]
    aa = str(Seq(nt).translate(to_stop=True))
    if aa and native_start and not aa.startswith("M"):
        aa = "M" + aa[1:]
    return aa


def rebuild(t, ev: dict, svtype: str) -> tuple[str, str | None]:
    """(status, mutant protein or None)."""
    p1, p2 = ev["pos1"], ev["pos2"]
    r1, i1 = region(t, p1)
    r2, i2 = region(t, p2)
    utr3 = t.three_prime_utr_sequence or ""
    if svtype == "DEL":
        if r1 == "intron" and r2 == "intron" and i1 == i2:
            return "same_intron", t.protein_sequence
        if r1 == "CDS" and r2 == "CDS":
            keep = []
            for s, e in sorted(t.coding_sequence_position_ranges):
                for a, b in ((s, min(e, p1)), (max(s, p2), e)):
                    if a <= b:
                        keep.append((a, b))
            keep = sorted(set(keep))
            keep = keep[::-1] if t.strand == "-" else keep
            # split at the junction to insert the exonic inserted bases
            ins = ev["insertion"] if t.strand == "+" else str(Seq(ev["insertion"]).reverse_complement())
            left = [b for b in keep if (b[1] <= p1 if t.strand == "+" else b[0] >= p2)]
            right = [b for b in keep if b not in left]
            nt = genomic_seq(t, left) + ins + genomic_seq(t, right) + utr3
            return "exonic_deletion", translate(nt, native_start=True)
        if r1 == "intron" and r2 == "intron":
            keep = [(s, e) for s, e in sorted(t.coding_sequence_position_ranges)
                    if e < p1 or s > p2]
            keep = keep[::-1] if t.strand == "-" else keep
            return "exon_skip", translate(genomic_seq(t, keep) + utr3, native_start=True)
    if svtype == "DUP" and r1 == "intron" and r2 == "intron" and i1 != i2:
        blocks = cds_blocks(t)
        inside = [b for b in blocks if b[0] >= p1 and b[1] <= p2]
        if inside:
            idx = blocks.index(inside[-1]) + 1
            dup = blocks[:idx] + inside + blocks[idx:]
            return "exon_duplication", translate(genomic_seq(t, dup) + utr3, native_start=True)
    if svtype == "DUP" and r1 == "intron" and r2 == "intron" and i1 == i2:
        return "same_intron", t.protein_sequence
    return f"not_modelled ({svtype} {r1}/{r2})", None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--run-dir", action="append", required=True, type=pathlib.Path,
                    help="a stage-2 output directory with *.all_neopeptides.txt and "
                         "*.admitted.vcf (or pass --vcf-dir)")
    ap.add_argument("--vcf-dir", action="append", default=None, type=pathlib.Path,
                    help="where each run's admitted VCF lives, if not in --run-dir "
                         "(same order as --run-dir)")
    ap.add_argument("--release", type=int, default=115)
    ap.add_argument("--out", required=True, type=pathlib.Path)
    args = ap.parse_args()

    import pyensembl
    genome = pyensembl.EnsemblRelease(args.release)
    rows = []
    for k, run in enumerate(args.run_dir):
        peps = pd.read_csv(next(run.glob("*.all_neopeptides.txt")), sep="\t", low_memory=False)
        vdir = args.vcf_dir[k] if args.vcf_dir else run
        vcf = read_vcf_events(next(vdir.glob("*.admitted.vcf")))
        same = peps[(peps["gene1"] == peps["gene2"]) & peps["svtype"].isin(["DEL", "DUP"])]
        for sv_id, g in same.groupby("sv_id"):
            ev = vcf.get(str(sv_id))
            t = genome.transcript_by_id(g["transcript_id1"].iloc[0])
            generated = set(g["neopeptide"])
            if ev is None:
                status, protein = "no VCF record", None
            else:
                status, protein = rebuild(t, ev, g["svtype"].iloc[0])
            if protein is None:
                expected, verdict = None, "not compared"
            else:
                expected = windows(protein) - windows(t.protein_sequence)
                verdict = "agree" if expected == generated else "DISAGREE"
            rows.append({
                "sample": run.name, "sv_id": sv_id, "gene": g["gene1"].iloc[0],
                "transcript": t.transcript_id, "svtype": g["svtype"].iloc[0],
                "model": status, "verdict": verdict,
                "n_generated": len(generated),
                "n_expected": None if expected is None else len(expected),
                "only_generated": "" if expected is None else ";".join(sorted(generated - expected)[:5]),
                "only_expected": "" if expected is None else ";".join(sorted(expected - generated)[:5]),
            })
    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    print(out.groupby(["model", "verdict"]).size().to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
