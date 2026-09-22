"""
synthesis.py — stage 9: the funnel, lineage attribution, and cross-sample tables.

THE FUNNEL IS THE PRIMARY DELIVERABLE, NOT THE FINAL NUMBER
-----------------------------------------------------------
A four-figure candidate count becoming a single-digit believable set is the
normal outcome — the field's own benchmark (TESLA, 28 teams, 6 patients) found
only ~6% of top-ranked predictions validate functionally. What matters is WHERE
each order of magnitude is lost, because that is what tells you whether the
analysis is filtering artefacts or filtering signal. So every sample emits the
same rows, in the same order, whatever the numbers turn out to be.

ATTRIBUTION ONLY MEANS SOMETHING IN A LINEAGE
---------------------------------------------
The earliest sample carrying a breakpoint is the one that explains it. In a
derived lineage this is the whole point: an event present in the parent cannot be
attributed to whatever was done to make the child. For unrelated samples
(`parent: null` everywhere) attribution is skipped rather than faked.
"""
from __future__ import annotations

import pandas as pd

from . import criteria

#: The funnel rows, fixed so samples are comparable line by line. Every sample
#: emits all of them; a stage that could not run emits `NA`, never 0 — a zero
#: from a stage that never ran is indistinguishable from a real negative.
FUNNEL_ROWS = [
    ("vcf_records", "SV records in the raw VCF"),
    ("admitted", "after FILTER / pairing / population filters"),
    ("junctions", "paired breakends collapsed to junctions"),
    ("candidate_peptides", "unique candidate neopeptides"),
    ("matches_identical", "identical to a reference peptide (rows: peptide x SV)"),
    ("matches_unique_peptides", "distinct peptides among those matches"),
    ("matches_gene_concordant", "and broken in the same gene"),
    ("matches_credible", "and high-complexity, non-self"),
    ("events", "credible matches collapsed to genomic events"),
    ("events_private", "PON and population-frequency clean"),
    ("events_hc", "high-confidence SV calls"),
    ("events_rna_supported", "junction-crossing reads in RNA"),
    # NOT the strongest set: privacy is absent from this one. Labelling it as
    # such reported an event present in 90% of a gnomAD population as a
    # surviving candidate, so the two are now separate rows and the last row is
    # the one to quote.
    ("events_hc_and_rna", "high-confidence AND transcribed (privacy NOT applied)"),
    ("events_private_hc_and_rna", "private AND high-confidence AND transcribed"),
]


def build_funnel(counts: dict) -> pd.DataFrame:
    """One sample's funnel, in the fixed row order."""
    return pd.DataFrame([{"step": key, "description": description,
                          "n": counts.get(key, "NA")}
                         for key, description in FUNNEL_ROWS])


def attribute(events: pd.DataFrame, sample_name: str, config,
              breakpoint_index: dict) -> pd.DataFrame:
    """Assign each event to the earliest sample in the lineage carrying it.

    `breakpoint_index` maps sample name -> {(chrom, pos), ...} for every admitted
    breakend of that sample. A match within ATTRIBUTION_MATCH_TOLERANCE bp counts
    as the same breakpoint.
    """
    if events.empty:
        return pd.DataFrame()

    ancestors = list(reversed(config.ancestors(sample_name)))   # root first
    tolerance = criteria.ATTRIBUTION_MATCH_TOLERANCE

    def present_in(sample: str, chrom, pos) -> bool:
        index = breakpoint_index.get(sample)
        if not index:
            return False
        chrom = str(chrom).replace("chr", "")
        return any((chrom, int(pos) + d) in index
                   for d in range(-tolerance, tolerance + 1))

    rows = []
    for _, event in events.iterrows():
        earliest, checked = sample_name, []
        for ancestor in ancestors:
            hit = present_in(ancestor, event["chrom1"], event["pos1"]) and \
                  present_in(ancestor, event["chrom2"], event["pos2"])
            checked.append(f"{ancestor}:{'yes' if hit else 'no'}")
            if hit:
                earliest = ancestor
                break
        rows.append({
            "sample_sv_id": event.get("sample_sv_id"), "gene": event.get("gene"),
            "found_in": sample_name, "attributed_to": earliest,
            "lineage_checked": ";".join(checked) if checked else "root sample",
            "interpretation": ("pre-existing in an ancestor — not attributable to "
                               "this sample's derivation"
                               if earliest != sample_name
                               else "acquired in this sample"),
        })
    return pd.DataFrame(rows)


def breakpoint_index(junctions: pd.DataFrame) -> set:
    """All breakend coordinates of a sample, for attribution lookups."""
    index = set()
    for _, row in junctions.iterrows():
        index.add((str(row["chrom1"]).replace("chr", ""), int(row["pos1"])))
        index.add((str(row["chrom2"]).replace("chr", ""), int(row["pos2"])))
    return index


def compare_samples(per_sample: dict) -> pd.DataFrame:
    """Cross-sample comparison table.

    READ THIS ON EVENTS AND RATES, NEVER ON PEPTIDE COUNTS. A high per-candidate
    match rate can be a single locus seen through a sliding window: in one
    observed run 21 "matches" were one 53 bp deletion. Both the raw count and the
    per-1,000-candidate rate are emitted for exactly this reason.
    """
    rows = []
    for sample, result in per_sample.items():
        counts = result.get("counts", {})
        row = {"sample": sample}
        row.update({key: counts.get(key, "NA") for key, _ in FUNNEL_ROWS})
        candidates = counts.get("candidate_peptides")
        matches = counts.get("matches_identical")
        if isinstance(candidates, int) and isinstance(matches, int) and candidates:
            row["matches_per_1000_candidates"] = round(1000.0 * matches / candidates, 3)
        null = result.get("null") or {}
        row["null_mean"] = null.get("null_mean", "NA")
        row["null_p"] = null.get("p_value", "NA")
        rows.append(row)
    return pd.DataFrame(rows)


def interpret(counts: dict) -> list[str]:
    """Plain-language cautions attached to every per-sample report.

    These are emitted with the numbers rather than left to the reader, because
    each corresponds to a documented mistake that produced a wrong conclusion.
    """
    notes = []
    if counts.get("events") == 0 and counts.get("matches_identical", 0) > 0:
        notes.append("Matches collapsed to zero events after credibility QC — "
                     "report the event count, not the match count.")
    if isinstance(counts.get("candidate_peptides"), int) \
            and counts["candidate_peptides"] < 1000:
        notes.append("Small candidate universe: raw counts are not comparable "
                     "with samples having orders of magnitude more candidates. "
                     "Use the per-1,000 rate and the null model.")
    if counts.get("events_rna_supported") == "NA":
        notes.append("No RNA for this sample: RNA rows are NA, not negative. "
                     "Absence of the modality is not absence of the junction.")
    if counts.get("events_private", 0) == 0 and counts.get("events", 0) > 0:
        notes.append("Every event is recurrent in a panel or population — real "
                     "calls, but common polymorphisms rather than private events.")
    return notes
