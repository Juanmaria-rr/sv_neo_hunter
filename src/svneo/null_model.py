"""
null_model.py — how many matches would this many peptides produce by chance?

WHY THIS EXISTS
---------------
Raw match counts are not comparable between samples whose candidate universes
differ by orders of magnitude. An observed run:

    sample A   63,036 candidate peptides -> 288 matches -> 39 events
    sample B        64 candidate peptides ->   0 matches ->  0 events
    sample C       536 candidate peptides ->  21 matches ->  1 event

"39 vs 0 vs 1" invites the reading that A is biologically special, when the
search spaces differ 1,000-fold. Two corrections are needed and both are reported:

    1. a RATE   — matches per 1,000 candidate peptides
    2. a NULL   — the matches expected by chance for a set of that size and
                  composition, giving an empirical p-value

Gene concordance (97.6% of matches landing in the same gene in one run) is strong
evidence against chance convergence, but it is an argument, not a test. This is
the test.

THE PERMUTATION
---------------
Each candidate peptide's residues are shuffled, preserving length and amino-acid
composition exactly. That holds constant the two properties that drive chance
matching — how many peptides there are and what they are made of — while
destroying the biological signal (the actual sequence). The fraction of
permutations reaching the observed match count is the p-value.

WHAT IT DOES NOT CONTROL FOR
----------------------------
Shuffling breaks the reference set's own composition bias only on the sample
side. If the reference catalogue is itself enriched for a residue composition
(as SV-derived frameshift peptides often are), a shuffled sample peptide is
still drawn from a different distribution than a real one. The null is therefore
CONSERVATIVE for the "is this above chance" question and should not be read as a
calibrated false-discovery rate.
"""
from __future__ import annotations

import random

import pandas as pd


def _shuffled(peptide: str, rng: random.Random) -> str:
    residues = list(peptide)
    rng.shuffle(residues)
    return "".join(residues)


def permutation_test(candidate_peptides, reference_peptides,
                     n_permutations: int = 1000, seed: int = 0) -> dict:
    """Empirical null for the number of identical-peptide matches.

    `candidate_peptides`  the sample's unique candidate peptides
    `reference_peptides`  the reference catalogue's unique peptides

    Returns observed count, the null distribution's mean/sd/max, an empirical
    p-value, and the enrichment ratio. A p-value of 0 is reported as
    "< 1/n_permutations" by the caller rather than as zero.
    """
    candidates = [str(p).strip().upper() for p in candidate_peptides if p]
    reference = {str(p).strip().upper() for p in reference_peptides if p}
    if not candidates or not reference:
        return {"observed": 0, "null_mean": 0.0, "p_value": None,
                "note": "empty candidate or reference set"}

    observed = sum(1 for p in set(candidates) if p in reference)

    rng = random.Random(seed)
    null_counts = []
    for _ in range(n_permutations):
        shuffled = {_shuffled(p, rng) for p in candidates}
        null_counts.append(sum(1 for p in shuffled if p in reference))

    mean = sum(null_counts) / len(null_counts)
    variance = sum((c - mean) ** 2 for c in null_counts) / len(null_counts)
    at_least = sum(1 for c in null_counts if c >= observed)

    return {
        "observed": observed,
        "n_candidates": len(set(candidates)),
        "n_reference": len(reference),
        "matches_per_1000_candidates": round(1000.0 * observed / len(set(candidates)), 3),
        "null_mean": round(mean, 3),
        "null_sd": round(variance ** 0.5, 3),
        "null_max": max(null_counts),
        "n_permutations": n_permutations,
        "p_value": at_least / n_permutations,
        "p_value_is_bounded": at_least == 0,   # true -> report as < 1/n
        "enrichment_over_null": round(observed / mean, 2) if mean > 0 else None,
    }


def summarise(results: dict) -> str:
    """One-line human summary, with the p-value reported honestly at the bound."""
    if results.get("p_value") is None:
        return f"null model not evaluated: {results.get('note', 'unknown reason')}"
    n = results["n_permutations"]
    p = (f"< {1 / n:g}" if results["p_value_is_bounded"]
         else f"= {results['p_value']:.3g}")
    return (f"{results['observed']} matches from {results['n_candidates']:,} "
            f"candidates ({results['matches_per_1000_candidates']}/1,000); "
            f"null {results['null_mean']} ± {results['null_sd']}, p {p}"
            + (f", {results['enrichment_over_null']}x over null"
               if results.get("enrichment_over_null") else ""))


def rate_table(per_sample: dict) -> pd.DataFrame:
    """Side-by-side rates for every sample — the table that makes raw counts
    comparable. Always report this next to the raw counts, never instead."""
    rows = []
    for sample, result in per_sample.items():
        rows.append({"sample": sample,
                     "candidates": result.get("n_candidates"),
                     "matches": result.get("observed"),
                     "per_1000": result.get("matches_per_1000_candidates"),
                     "null_mean": result.get("null_mean"),
                     "enrichment": result.get("enrichment_over_null"),
                     "p_value": result.get("p_value")})
    return pd.DataFrame(rows)
