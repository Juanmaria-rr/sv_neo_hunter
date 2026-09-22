#!/usr/bin/env python3
"""
extract_cohort_alleles.py — the HLA class I alleles actually observed in the cohort.

WHY THIS PANEL RATHER THAN A STANDARD ONE
-----------------------------------------
Binding prediction needs an allele set, and the obvious move is a published
reference panel (IEDB's 27, say). That would be defensible in isolation and
wrong here, because it answers a different question from the one the existing
analysis answered. The cohort's patients are typed, so this cohort's allele universe
is a fact to be read off the data, not a choice to be made.

THE SOURCE MATTERS, AND THE OBVIOUS ONE IS BIASED
-------------------------------------------------
There are two places to read alleles from, and they do not agree:

  consensus_HLA.csv        4,439 typed samples -> 221 distinct alleles
  all_neoantigens.tsv      3,026 samples       -> 187 distinct alleles

The second is tempting because it sits beside the peptides, and it is wrong for
this purpose: it lists alleles that appeared *as binders* after filtering, so an
allele carried by patients that never presented an SV neopeptide is absent. 34
alleles are missing for exactly that reason, and the omission is not random —
they are the restrictive ones, so coverage computed over the 187 is optimistic.

The tell is in the per-sample counts: 134 samples show three alleles or fewer,
which no amount of homozygosity can produce, since a fully homozygous individual
still carries three distinct class I alleles. (The 852 samples showing five are
ordinary single-locus homozygosity and mean nothing.)

Reading the genotypes directly avoids all of this. As a consistency check, every
one of the 187 binder alleles does appear in the genotype file — the biased set
is a strict subset of the unbiased one, which is what it should be.

WHAT IT DOES NOT DO
-------------------
It does not weight alleles by frequency, and the resulting list is therefore not
a population. Frequencies live in `top_50_hla_class_I.xlsx` (AFND-derived) and
are applied downstream, where population coverage is computed as
1 - prod_i (1 - f_i)^2 under Hardy-Weinberg. Keeping the two apart means the
prediction does not have to be redone if the population model changes.

Alleles absent from the netMHCpan install are reported and excluded — rare
four-digit types appear in a cohort of thousands that no predictor covers, and
silently dropping them would overstate coverage.

    python tools/extract_cohort_alleles.py \
        --genotypes ../cohort/consensus_HLA.csv \
        --out reference/cohort_observed_alleles.txt
"""
from __future__ import annotations

import argparse
import pathlib

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parent.parent


def supported_alleles(netmhcpan_dir: pathlib.Path) -> set[str]:
    """Allele names this netMHCpan install knows, from data/allelenames."""
    names = netmhcpan_dir / "data" / "allelenames"
    if not names.exists():
        return set()
    return {line.split()[0] for line in names.read_text().splitlines()
            if line.strip()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genotypes", required=True,
                        help="wide CSV of per-sample class I genotypes: "
                             "sampleId, HLA-A_1, HLA-A_2, HLA-B_1, ...")
    parser.add_argument("--out", default="reference/cohort_observed_alleles.txt")
    parser.add_argument("--netmhcpan-dir", default=None,
                        help="to check support; defaults to the local install")
    args = parser.parse_args()

    genotypes = pd.read_csv(args.genotypes)
    allele_columns = [c for c in genotypes.columns if c != "sampleId"]

    # The genotype file writes `A*02:01`; netMHCpan wants `HLA-A02:01`. It
    # rejects the asterisk by returning zero rows without erroring, which is the
    # worst possible failure mode, so both fixes are applied here and the
    # written file is directly usable.
    raw = set()
    for column in allele_columns:
        raw |= set(genotypes[column].dropna().astype(str).str.strip())

    observed = sorted({
        ("HLA-" + a if not a.startswith("HLA-") else a).upper().replace("*", "")
        for a in raw if a and a.lower() != "nan"})

    install = pathlib.Path(args.netmhcpan_dir) if args.netmhcpan_dir else \
        next(iter(sorted(REPO.glob("netMHCpan-*/"))), None)
    known = supported_alleles(install) if install else set()

    if known:
        usable = [a for a in observed if a in known]
        missing = [a for a in observed if a not in known]
    else:
        usable, missing = observed, []
        print("  WARNING: no netMHCpan install found; support not checked")

    samples = genotypes["sampleId"].nunique()
    per_sample = genotypes[allele_columns].nunique(axis=1)

    header = [
        "# HLA class I alleles carried by the reference cohort.",
        "#",
        f"# Extracted by tools/extract_cohort_alleles.py from {pathlib.Path(args.genotypes).name}",
        f"# {samples:,} typed samples; median {int(per_sample.median())} distinct",
        "# alleles per sample (six unless a locus is homozygous).",
        f"# {len(observed)} distinct alleles, {len(usable)} supported by netMHCpan.",
        "#",
        "# Read from the GENOTYPES, not from the binder-filtered neoantigen",
        "# table: the latter omits alleles that never presented anything, which",
        "# biases the panel towards permissive alleles. See the module docstring.",
        "#",
        "# This is the cohort's allele universe, not a population: it carries no",
        "# frequency weighting. Population coverage is computed downstream from",
        "# the AFND top-50 table under Hardy-Weinberg, 1 - prod(1 - f)^2.",
        "#",
        "# Format: netMHCpan names, asterisk stripped (HLA-A*02:01 silently",
        "# returns nothing; HLA-A02:01 works).",
    ]
    if missing:
        header += ["#",
                   f"# EXCLUDED, not in this netMHCpan install ({len(missing)}):",
                   *[f"#   {a}" for a in missing]]

    out = REPO / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(header) + "\n\n" + "\n".join(usable) + "\n")

    print(f"  samples:            {samples:,}")
    print(f"  observed alleles:   {len(observed)}")
    print(f"  supported:          {len(usable)}")
    if missing:
        print(f"  excluded ({len(missing)}):      {', '.join(missing)}")
    print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()
