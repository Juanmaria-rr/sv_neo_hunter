#!/usr/bin/env python
"""
run_mhcflurry.py — MHCflurry class I predictions for a peptide list, one row per
(peptide, allele), as a second predictor beside netMHCpan.

WHY A SECOND PREDICTOR, AND WHY ONLY AS AN ANNOTATION
-----------------------------------------------------
Which peptides a predictor calls binders moves between predictors and versions;
the magnitudes are stable, the identities are not. A call that two predictors of
different families agree on is more robust than one only netMHCpan makes. That
is the question this answers. It does NOT redefine `presentable`: two predictors
with no rule for combining them would give two answers, and any rule chosen
after seeing the results could be tuned to them. netMHCpan remains the criterion;
MHCflurry is recorded beside it (tools/annotate_mhcflurry.py).

WHAT IS PREDICTED
-----------------
The presentation model (`Class1PresentationPredictor`), which combines binding
affinity with antigen processing — the counterpart of netMHCpan's eluted-ligand
score. Each allele is predicted as its own one-allele "sample", so every
(peptide, allele) pair gets a row, not just the best allele.

    allele, peptide, affinity_nM, processing_score, presentation_score,
    presentation_percentile

Usage
-----
    python tools/run_mhcflurry.py --peptides peptides.txt \\
        --alleles HLA-A*24:02,HLA-B*35:01,... --out mhcflurry.tsv
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import pandas as pd


#: Binder rule for MHCflurry, the counterpart of run_netmhcpan.BINDER_THRESHOLDS:
#: every cut must hold at once. Affinity <= 500 nM is the classical binder line
#: shared with netMHCpan; presentation percentile <= 2 is MHCflurry's analogue of
#: %Rank_EL <= 2 (affinity and processing combined, ranked against random
#: peptides for that allele). MHCflurry's presentation model has no separate
#: affinity percentile, so the rule has two cuts where netMHCpan's has three.
BINDER_THRESHOLDS = {"affinity_nM": 500.0, "presentation_percentile": 2.0}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--peptides", required=True,
                   help="plain peptide list, or a TSV with --peptide-column")
    p.add_argument("--peptide-column", default=None)
    p.add_argument("--alleles", required=True,
                   help="comma-separated (HLA-A*01:01 or HLA-A01:01), or a file of names")
    p.add_argument("--out", required=True)
    p.add_argument("--chunk-size", type=int, default=20000)
    args = p.parse_args()

    from mhcflurry import Class1PresentationPredictor, __version__

    src = pathlib.Path(args.peptides)
    if args.peptide_column:
        peptides = pd.read_csv(src, sep="\t", usecols=[args.peptide_column])[args.peptide_column]
    else:
        peptides = pd.Series(src.read_text().split())
    peptides = sorted(set(peptides.dropna().astype(str)))
    raw = (pathlib.Path(args.alleles).read_text().split() if pathlib.Path(args.alleles).exists()
           else args.alleles.split(","))
    alleles = []
    for a in raw:
        a = a.strip()
        if a and "*" not in a:                       # HLA-A01:01 -> HLA-A*01:01
            a = a[:5] + "*" + a[5:]
        if a:
            alleles.append(a)

    predictor = Class1PresentationPredictor.load()
    usable = [x for x in peptides if 8 <= len(x) <= 15]
    print(f"  MHCflurry {__version__}: {len(usable):,} peptides x {len(alleles)} alleles",
          file=sys.stderr)
    frames = []
    for i in range(0, len(usable), args.chunk_size):
        chunk = usable[i:i + args.chunk_size]
        out = predictor.predict(chunk, {a: [a] for a in alleles}, verbose=0)
        frames.append(out)
        print(f"    {min(i + args.chunk_size, len(usable)):,}/{len(usable):,}", file=sys.stderr)
    t = pd.concat(frames, ignore_index=True)
    t = t.rename(columns={"sample_name": "allele", "affinity": "affinity_nM"})
    t = t[["allele", "peptide", "affinity_nM", "processing_score",
           "presentation_score", "presentation_percentile"]]
    t.to_csv(args.out, sep="\t", index=False)
    print(f"  wrote {len(t):,} rows -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
