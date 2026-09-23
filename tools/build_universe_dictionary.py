#!/usr/bin/env python
"""
build_universe_dictionary.py — the candidate universe explained, thresholds included.

WHY NOT JUST THE MEANINGS
-------------------------
`column_meanings.py` says what each column is, and names the constants behind it
— `HC_MIN_SEGMAPQ`, `PON_MAX` — without their values. Anyone reading a table
then has to go and look them up, and anyone quoting a threshold from memory is
one edit away from being wrong.

So this emits the values too, and **imports them from the modules that define
them** rather than restating them. If a threshold moves in `criteria.py`, this
dictionary moves with it; it cannot describe a cut the pipeline is not applying.

It also answers the question the meanings do not: **does this column ever remove
a row?** Most do not. `sv_hc`, `pon_count`, `matched_reference` and
`presentable` are measurements carried alongside the data, and reading them as
filters is the single most common way to misread these tables.

Usage
-----
    python tools/build_universe_dictionary.py \\
        --universe <dir>/candidate_universe.tsv --out <dir>/..._dictionary.tsv
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
import column_meanings                                    # noqa: E402
from svneo import criteria                                # noqa: E402
from run_netmhcpan import BINDER_THRESHOLDS               # noqa: E402

C = criteria
BT = BINDER_THRESHOLDS

#: column -> (stage, kind, rule/threshold). The rule strings interpolate the
#: real constants; none of these numbers is typed out by hand.
#:
#: `kind` is the part people get wrong:
#:   identity     what this row is
#:   annotation   a label carried from upstream
#:   measurement  a number measured on the data
#:   verdict      a boolean the pipeline computed — say whether it filters
SPEC: dict[str, tuple] = {
    # column: (stage, level, kind, derived_from, what it is + any threshold)
    #
    # `level` is the field this documentation was missing. A reader cannot judge
    # a column without knowing whether it is a DNA call, an RNA measurement, an
    # in-silico translation or a prediction — and several columns that look like
    # observations are none of those.

    # -- lineage -----------------------------------------------------------
    "acquired_in": ("lineage", "bookkeeping", "identity", "the run each candidate table came from",
        "The line whose OWN variant calls produce this peptide — where it entered the "
        "lineage. Semicolon-separated if two lines generate the same sequence "
        "independently. **Filter on this to ask what a knockout added, and to stay "
        "inside a line-restricted subset**; filtering on `present_in` instead pulls "
        "back every line that merely inherits it. No threshold."),
    "present_in": ("lineage", "bookkeeping", "identity", "acquired_in + the configured parent/child lineage",
        "Every line carrying it, inheritance included. No threshold."),
    "n_lines_present": ("lineage", "bookkeeping", "measurement", "count of present_in",
        "How many lines of the lineage carry this peptide. 1 means it is confined to "
        "one branch; the maximum equals the number of lines in the run, reached by "
        "anything the root acquired. No threshold."),

    # -- stage 1, the SV call itself (DNA) ---------------------------------
    "sv_id": ("stage 1", "DNA — caller", "identity", "the SV VCF", "Caller's record ID for the junction."),
    "chrom1": ("stage 1", "DNA — caller", "identity", "the SV VCF", "Chromosome of breakend 1."),
    "pos1": ("stage 1", "DNA — caller", "identity", "the SV VCF", "Position of breakend 1, 1-based."),
    "chrom2": ("stage 1", "DNA — caller", "identity", "the SV VCF",
        "Chromosome of the 3' breakend, without a `chr` prefix. Differing from "
        "`chrom1` makes the junction inter-chromosomal, which forces `test` to "
        "`chimeric` and leaves `event_size` undefined."),
    "pos2": ("stage 1", "DNA — caller", "identity", "the SV VCF", "Position of breakend 2."),
    "gene1": ("stage 1", "annotation", "annotation", "Ensembl, by overlap with breakend 1",
        "Gene at the 5' breakend. Reference annotation, not a measurement."),
    "gene2": ("stage 1", "annotation", "annotation", "Ensembl, by overlap with breakend 2",
        "Gene overlapping the 3' breakend, from the reference annotation rather than "
        "from any measurement. It contributes the downstream half of the fusion "
        "protein, so a frameshift makes every residue it contributes novel."),
    "strand1": ("stage 1", "annotation", "annotation", "Ensembl transcript strand",
        "Strand of the 5' transcript. Emitted so results can be stratified by it."),
    "strand2": ("stage 1", "annotation", "annotation", "Ensembl transcript strand",
        "Strand of the 3' transcript, `+` or `-`. Emitted so results can be "
        "stratified by strand pair; note that per-peptide proportions by strand are "
        "dominated by a few prolific loci, so stratify per event."),
    "svtype": ("stage 1", "DNA — caller", "annotation", "the SV VCF",
        "DEL / DUP / INV / BND / TRA, as the caller typed it."),

    # -- stage 2, the predicted fusion protein (IN SILICO) -----------------
    "neopeptide": ("stage 2", "protein — in silico", "identity",
        "translation of the reconstructed fusion transcript",
        "The candidate peptide. Predicted from coordinates and reference sequence; "
        "no read supports it directly."),
    "pep_length": ("stage 2", "protein — in silico", "measurement", "len(neopeptide)",
        f"Length. Windows are generated at {C.PEPTIDE_LENGTHS} (`PEPTIDE_LENGTHS`), "
        "the canonical MHC class I range."),
    "junction_aa": ("stage 2", "protein — in silico", "measurement",
        "(5' CDS nucleotides + inserted bases) // 3",
        "Residue index at which the two partners meet in the fusion protein. "
        "Empty when the 5' coding sequence could not be resolved."),
    "spans_junction": ("stage 2", "protein — in silico", "annotation",
        "position of the peptide within the fusion protein vs `junction_aa`",
        "**Does this peptide straddle the fusion point in the PROTEIN?** True when "
        "`start < junction_aa < start + len(peptide)`. A sliding window also yields "
        "peptides lying wholly on one side: those are ordinary peptides of an intact "
        "reading frame that happen to sit in a fused transcript. "
        "**This is geometry, not evidence** — it says nothing about whether the "
        "junction is transcribed (that is `rna_tier`) or expressed. "
        "`REQUIRE_SPANS_JUNCTION = " + str(C.REQUIRE_SPANS_JUNCTION) + "`: it is NOT "
        "used as a filter, because after a frameshift every residue downstream of the "
        "junction is novel too even though it does not straddle it. Enabling it "
        "discarded 63 of 64 catalogue matches. Empty = junction offset unknown, never "
        "reported as False."),
    "frame_effect": ("stage 2", "protein — in silico", "annotation",
        "reading frame of the fused CDS vs the 5' partner's",
        "What the fusion does to translation: `In-frame` (frame preserved), "
        "`Stop-gain` (a premature stop appears), `Stop-loss` (the normal stop is "
        "removed and translation reads through), `Start-loss` (the 5' segment "
        "contributes no start). Derived, not observed. NOTE: per-peptide proportions "
        "of this column are dominated by a few prolific loci, because one event "
        "yields many overlapping windows — count per event."),

    # -- stage 3, sequence QC ----------------------------------------------
    "is_self": ("stage 3 QC", "protein — vs reference proteome", "verdict — FILTERS",
        "exact substring search against Ensembl `pep.all`",
        "**Does this exact sequence already exist in a normal human protein?** If so "
        "it is not neo-anything and the immune system is tolerant to it. Proteins are "
        "NUL-separated so no peptide can match across a protein boundary, which would "
        "invent a junction sequence. Exact match only: one residue from self is not "
        "flagged. Empty = no proteome configured, never 'passed'."),
    "low_complexity": ("stage 3 QC", "protein — composition", "verdict — FILTERS",
        "composition of the peptide string",
        f"**Is the sequence compositionally trivial?** Such a string matches a "
        f"catalogue by composition rather than by descent. Flagged if ANY holds: "
        f"Shannon entropy < {C.LC_MIN_SHANNON_ENTROPY}; one residue >= "
        f"{C.LC_MAX_SINGLE_AA_FRACTION:.0%} of the peptide; a homopolymer run >= "
        f"{C.LC_MAX_HOMOPOLYMER_RUN}; or <= {C.LC_MIN_DISTINCT_RESIDUES} distinct "
        f"residues (`LC_*` in criteria.py)."),

    # -- stage 6, is the DNA call believable? ------------------------------
    "sv_hc": ("stage 6", "DNA — read evidence", "verdict — does NOT filter",
        "the caller's per-breakend fields",
        f"**Is the SV call itself believable?** ALL of: segmapq >= {C.HC_MIN_SEGMAPQ}, "
        f"variant fragments >= {C.HC_MIN_VF}, caller qual >= {C.HC_MIN_QUAL}; plus, for "
        f"intra-chromosomal events, size >= {C.HC_MIN_SV_SIZE} bp (`HC_MIN_*`). "
        "**Reported, never enforced** — a call can be perfect and still be a common "
        "polymorphism, which is what the privacy columns are for."),
    "vf_bp1": ("stage 6", "DNA — read evidence", "measurement", "caller VF field",
        f"DNA fragments supporting breakend 1. Enters `sv_hc` at >= {C.HC_MIN_VF}. "
        "Compare against `junction_reads`: RNA far above DNA is the splicing-artefact "
        "signature."),
    "vf_bp2": ("stage 6", "DNA — read evidence", "measurement", "caller VF field",
        f"DNA fragments supporting the 3' breakend. Enters `sv_hc` at >= {C.HC_MIN_VF}. "
        "Like `vf_bp1`, compare it against `junction_reads`: RNA support far exceeding "
        "DNA support is the splicing-artefact signature."),
    "qual_bp1": ("stage 6", "DNA — read evidence", "measurement", "caller QUAL",
        f"Caller confidence, breakend 1. Enters `sv_hc` at >= {C.HC_MIN_QUAL}."),
    "qual_bp2": ("stage 6", "DNA — read evidence", "measurement", "caller QUAL",
        f"The SV caller's own confidence score for the 3' breakend, on its Phred-like "
        f"scale. Enters `sv_hc` at >= {C.HC_MIN_QUAL}."),
    "segmapq_bp1": ("stage 6", "DNA — read evidence", "measurement", "caller SEGMAPQ",
        f"Mapping quality of the highest-contributing segment, breakend 1. Low values "
        f"mean the region is repetitive. Enters `sv_hc` at >= {C.HC_MIN_SEGMAPQ}."),
    "segmapq_bp2": ("stage 6", "DNA — read evidence", "measurement", "caller SEGMAPQ",
        f"Mapping quality of the highest-contributing segment at the 3' breakend, "
        f"0-60. Low values mean that end sits in repetitive sequence and the locus "
        f"assignment is doubtful. Enters `sv_hc` at >= {C.HC_MIN_SEGMAPQ}."),

    # -- stage 6, is it the sample's own? ----------------------------------
    "pon_count": ("stage 6", "DNA — cohort of normals", "measurement", "panel-of-normals annotation",
        f"**In how many unrelated normal genomes has this junction been seen?** An "
        f"ABSOLUTE count, not a frequency — read it against `panel_size_estimate`, "
        f"since the same count means opposite things against panels of different "
        f"sizes. Threshold `PON_MAX = {C.PON_MAX}`, applied on the PON10 branch "
        f"only. Empty counts as {C.PON_ABSENT_MEANS}, because treating absent as high "
        "would discard the cleanest breakends."),
    "pon_fraction": ("stage 6", "DNA — cohort of normals", "measurement", "pon_count / panel size",
        "`pon_count` expressed as a fraction of the panel — the interpretable form. "
        "A given raw count means opposite things against panels of different sizes, "
        "so quote this alongside it."),
    "pass_pon": ("stage 6", "DNA — cohort of normals", "verdict — does NOT filter here",
        "pon_count vs PON_MAX", f"`pon_count < {C.PON_MAX}`. On the noPON branch this is "
        "recorded and deliberately excluded from `is_private`."),
    "gnomad_af_used": ("stage 6", "DNA — population database", "measurement",
        "gnomAD-SV v4.1, reciprocal-overlap match",
        f"**How common is this variant in the general population?** Taken from "
        f"`{C.GNOMAD_AF_FIELD}`, matched to a gnomAD-SV record by "
        f"{C.GNOMAD_RECIPROCAL_OVERLAP:.0%} reciprocal overlap. **NaN means not "
        "evaluated**, which is not the same as rare — inter-chromosomal junctions are "
        "left unmatched because gnomAD-SV does not represent them comparably."),
    "gnomad_af_popmax": ("stage 6", "DNA — population database", "measurement", "gnomAD-SV per-ancestry allele frequencies",
        "Highest frequency across ancestry groups. Preferred over a global average, "
        "which is dominated by the largest group and can hide a variant that is rare "
        "worldwide but common where the donor came from."),
    "pass_gnomad": ("stage 6", "DNA — population database", "verdict — does NOT filter",
        "gnomad_af_used vs GNOMAD_MAX_AF",
        f"AF < {C.GNOMAD_MAX_AF} (`GNOMAD_MAX_AF`). **True when AF is NaN** — absence of "
        "evidence reads as a pass, so check `gnomad_af_used` before believing it."),
    "is_private": ("stage 6", "DNA — derived verdict", "verdict — does NOT filter",
        "pass_pon AND pass_gnomad",
        "**Is this the sample's own variant rather than common germline variation?** "
        "On the noPON branch only the gnomAD half counts, so it is half an answer; "
        "`privacy_note` says which halves ran."),
    "privacy_note": ("stage 6", "bookkeeping", "annotation", "which privacy tests ran",
        "Plain text saying which halves of the privacy test actually ran, so "
        "`is_private` can be read for what it is. On a panel-reporting branch it "
        "begins `gnomAD only — panel count reported, not applied`, and gains "
        "`; this event IS in the panel` when the panel has seen the junction in "
        "unrelated normals — i.e. **the panel contradicts `is_private` and was "
        "ignored by design**. On a panel-applying branch it instead reads "
        "`gnomAD not evaluated — PON only`, or is empty when both halves ran."),

    # -- stage 7, expression (RNA) -----------------------------------------
    "gene1_TPM": ("stage 7", "RNA — quantification", "measurement", "Isofox transcript quantification",
        f"Expression of the 5' gene in this sample's RNA. `expressed` needs >= "
        f"{C.TPM_EXPRESSED} TPM."),
    "gene2_TPM": ("stage 7", "RNA — quantification", "measurement", "Isofox transcript quantification",
        f"Expression of the 3' gene in this sample's RNA, in transcripts per million. "
        f"`expressed` requires >= {C.TPM_EXPRESSED} TPM here AND at the 5' gene."),
    "min_side_TPM": ("stage 7", "RNA — quantification", "measurement", f"{C.COVERAGE_SUMMARY}(gene1_TPM, gene2_TPM)",
        "The quieter of the two sides. `min`, never `max` (`COVERAGE_SUMMARY`): taking "
        "the maximum reports a silent locus as covered because its partner is "
        "transcribed."),
    "expressed": ("stage 7", "RNA — quantification", "verdict — does NOT filter", "both TPMs vs TPM_EXPRESSED",
        f"TPM >= {C.TPM_EXPRESSED} at BOTH breakends "
        f"(`EXPRESSION_USES_BOTH_BREAKENDS = {C.EXPRESSION_USES_BOTH_BREAKENDS}`). An "
        "event is in transcribed territory only if both ends are."),

    # -- stage 8, direct read evidence for the junction (RNA) --------------
    "junction_reads": ("stage 8", "RNA — direct read counting", "measurement", "the RNA BAM, counted by this pipeline",
        f"**Reads that actually CROSS the junction** — not coverage, which is context: "
        f"a breakpoint inside a highly expressed gene has thousands of reads whether or "
        f"not the junction exists. Counted by a CIGAR N gap matching the event size "
        f"within {C.NGAP_SIZE_TOLERANCE} bp and starting within "
        f"{C.NGAP_POSITION_TOLERANCE} bp of the breakpoint, or a supplementary "
        f"alignment landing within {C.SA_PARTNER_TOLERANCE} bp of the partner. "
        f"MAPQ >= {C.MIN_READ_MAPQ}; a fragment is counted once however many alignment "
        "records it produces."),
    "rna_tier": ("stage 8", "RNA — direct read counting", "verdict — does NOT filter", "junction_reads",
        f"STRONG >= {C.STRONG_MIN_JUNCTION_READS}, SUGGESTIVE >= "
        f"{C.SUGGESTIVE_MIN_JUNCTION_READS}, WEAK >= 1, NONE = 0, UNTESTABLE when no "
        "test applies. Soft clips and coverage never tier an event "
        "(`SOFTCLIPS_TIER_EVENTS`, `COVERAGE_TIERS_EVENTS` both False): clipped reads "
        "whose supplementary alignments land megabases away end near the breakpoint "
        "and demonstrably do not cross it."),
    "test": ("stage 8", "RNA — direct read counting", "annotation",
        "the event's geometry, decided BEFORE any read is examined",
        "**Which of the three counting mechanisms is applicable to this junction**, "
        "chosen from its geometry alone. One of `insertion` (the junction inserts at "
        f"least {C.MIN_INSERT_LEN} bp), `chimeric` (inter-chromosomal, or BND/TRA), "
        f"`sizegap` (intra-chromosomal with a resolvable size >= "
        f"{C.MIN_TESTABLE_GAP_SIZE} bp), or `none`. **It gates which "
        "`junction_by_*` column can be non-zero**: `sizegap` can only produce ngap "
        "counts, `chimeric` only sa, `insertion` only insert. Order matters — an "
        "insertion is tested as an insertion even when its breakend span is tiny, "
        "because the span is not the lesion. `none` means no mechanism applies and "
        "`rna_tier` is then `UNTESTABLE`, not `NONE`: untested is not tested-negative."),
    "test_reason": ("stage 8", "RNA — direct read counting", "annotation",
        "the geometry that selected `test`, with its measured values",
        "**Why that test and not another, in words, with the numbers that decided "
        "it.** Four shapes: `insertion-driven (insert_len=N)`; "
        "`intra-chromosomal N bp: CIGAR N-gap test`; `inter-chromosomal: needs "
        f"SA-to-partner / fusion caller`; and, when `test` is `none`, either "
        f"`N bp < {C.MIN_TESTABLE_GAP_SIZE} bp: below the aligner's minimum intron, "
        "no N-gap can exist` or `intra-chromosomal with no resolvable size`. Read it "
        "when a junction has no reads: it says whether that is a real negative or a "
        "test that could never have fired."),

    # -- Isofox's own calls, independent of the above ----------------------
    "nearest_alt_sj_bp": ("stage 7", "RNA — Isofox splice calls", "measurement", "Isofox alternative-splice-junction calls",
        "**Distance to the nearest alternative splice junction Isofox calls.** No "
        "threshold — a diagnostic. A breakpoint sitting beside a splice site in a "
        "highly expressed gene inherits that site's spliced reads as apparent junction "
        "evidence. Small here with large `alt_sj_frags` and small `vf_*` is the "
        "artefact signature."),
    "alt_sj_frags": ("stage 7", "RNA — Isofox splice calls", "measurement", "Isofox",
        "Fragments supporting that alternative junction. Compare against `vf_bp1`."),
    "alt_sj_type": ("stage 7", "RNA — Isofox splice calls", "annotation", "Isofox",
        "e.g. NOVEL_INTRON. Independent of this pipeline's own read counting, which is "
        "why it can contradict `junction_reads`."),
    "alt_sj_within_window": ("stage 7", "RNA — Isofox splice calls", "verdict — does NOT filter", "Isofox",
        "Whether an alternative splice junction falls inside the search window."),
    "retained_intron": ("stage 7", "RNA — Isofox splice calls", "measurement", "Isofox", "Retained-intron calls near the junction."),
    "isofox_fusion": ("stage 7", "RNA — Isofox splice calls", "annotation", "Isofox", "Isofox's own fusion call, if any."),
    "isofox_fusion_support": ("stage 7", "RNA — Isofox splice calls", "measurement", "Isofox fusion calling",
        "Fragments Isofox counts behind its own fusion call at this locus. Where both "
        "this and `junction_reads` are populated they are two tools counting the same "
        "junction independently, which is the one place a direct number-against-number "
        "cross-check is possible; usually empty, since Isofox only calls fusions it "
        "recognises."),

    # -- junction evidence by mechanism (RNA, our own counting) ------------
    "junction_by_ngap": ("stage 8", "RNA — direct read counting", "measurement",
        "pysam over the RNA BAM: CIGAR N gaps",
        f"Fragments supporting the junction through a CIGAR `N` gap matching the "
        f"event size within {C.NGAP_SIZE_TOLERANCE} bp and starting within "
        f"{C.NGAP_POSITION_TOLERANCE} bp of the breakpoint. `N` (skip), never `D` "
        "(deletion): `D` found zero reads across a 221 bp deletion where `N` found "
        "31. **Support concentrated here with none by `sa` is the shape splicing "
        "makes** — cross-check `nearest_alt_sj_bp`."),
    "junction_by_sa": ("stage 8", "RNA — direct read counting", "measurement",
        "pysam: SA tag vs the partner breakend",
        f"Fragments supporting it through a supplementary alignment landing within "
        f"{C.SA_PARTNER_TOLERANCE} bp of the partner. The mechanism for chimeric "
        "junctions, where no single read carries a gap."),
    "junction_by_insert": ("stage 8", "RNA — direct read counting", "measurement",
        "pysam: inserted sequence", f"Fragments supporting it through an inserted "
        f"sequence of at least {C.MIN_INSERT_LEN} bp."),
    "coverage_bp1": ("stage 8", "RNA — direct read counting", "measurement",
        f"alignments within {C.COVERAGE_WINDOW} bp of breakend 1",
        "**Context, not evidence.** A breakpoint inside a highly expressed gene has "
        "thousands of reads whether or not the junction exists, which is why "
        "coverage never tiers an event."),
    "coverage_bp2": ("stage 8", "RNA — direct read counting", "measurement",
        f"alignments within {C.COVERAGE_WINDOW} bp of breakend 2",
        "RNA alignments around the 3' breakend. Context for how much signal was "
        "available there, never evidence that the junction exists."),
    "min_coverage": ("stage 8", "RNA — direct read counting", "measurement", "min of the two",
        "The quieter side. `min`, never `max`."),
    "softclip_bp1": ("stage 8", "RNA — direct read counting", "measurement",
        f"soft clips >= {C.MIN_SOFTCLIP_LEN} bp at breakend 1",
        f"**Diagnostic only — soft clips never tier an event** "
        f"(`SOFTCLIPS_TIER_EVENTS = {C.SOFTCLIPS_TIER_EVENTS}`). A clip shows a read "
        "ENDS at the breakpoint; it says nothing about where the rest went. Clipped "
        "reads whose supplementary alignments landed megabases away once produced a "
        "candidate that did not survive."),
    "softclip_bp2": ("stage 8", "RNA — direct read counting", "measurement",
        f"soft clips >= {C.MIN_SOFTCLIP_LEN} bp at breakend 2",
        "Reads whose soft clip sits at the 3' breakend. Diagnostic only: a clip shows "
        "a read ENDS there, not where the rest of it went."),
    "alignments_bp1": ("stage 8", "RNA — direct read counting", "measurement", "pysam",
        "Alignment records examined at breakend 1, before the mapping-quality cut."),
    "alignments_bp2": ("stage 8", "RNA — direct read counting", "measurement", "pysam",
        "Alignment records examined at the 3' breakend before the mapping-quality "
        "cut. The denominator for `low_mapq_bp2`."),
    "low_mapq_bp1": ("stage 8", "RNA — direct read counting", "measurement", "pysam",
        f"Alignments discarded at breakend 1 for MAPQ < {C.MIN_READ_MAPQ}. A high "
        "share means the region is repetitive and the locus assignment is doubtful."),
    "low_mapq_bp2": ("stage 8", "RNA — direct read counting", "measurement", "pysam",
        f"Alignments discarded at the 3' breakend for MAPQ < {C.MIN_READ_MAPQ}. A "
        "high share of `alignments_bp2` means that end is repetitive."),

    # -- lesion size --------------------------------------------------------
    "span": ("stage 1", "DNA — caller", "measurement", "|pos2 - pos1|",
        "Raw coordinate difference. **Not the lesion size** when an insertion is "
        "involved: a 34 bp insertion was once reported as a 2 bp deletion because "
        "the span was taken for the lesion. Compare against `event_size`."),
    "event_size": ("stage 1", "DNA — derived", "measurement", "type-aware size function",
        f"Deleted, duplicated or inserted length as the SV type implies. Empty for "
        f"BND/TRA/SGL, where size is not meaningful. Events below "
        f"{C.MIN_TESTABLE_GAP_SIZE} bp are UNTESTABLE for RNA junction evidence."),
    "insert_len": ("stage 1", "DNA — caller", "measurement", "the VCF's inserted sequence",
        "Inserted bases at the junction. Non-zero here is why `span` and "
        "`event_size` disagree."),

    # -- common fragile sites, two catalogues ------------------------------
    "cfs_narrow_bp1": ("annotation", "genomic interval", "annotation",
        "conservative CFS catalogue, GRCh38",
        "Fragile site containing breakend 1 under the CONSERVATIVE catalogue "
        "(~1.2% of the genome). Empty = not inside one. No threshold: a position "
        "either falls in an interval or it does not."),
    "cfs_narrow_bp2": ("annotation", "genomic interval", "annotation", "conservative CFS catalogue",
        "Fragile site containing the 3' breakend under the conservative catalogue "
        "(~1.2% of the genome). Empty = not inside one. Annotated separately from "
        "`cfs_narrow_bp1` because the two ends can fall on different sides of a "
        "boundary."),
    "cfs_narrow_tier_bp1": ("annotation", "genomic interval", "annotation", "the catalogue's own tier column",
        "`core` (0.5% of the genome) or `extended` (0.7%), at breakend 1."),
    "cfs_narrow_tier_bp2": ("annotation", "genomic interval", "annotation", "the catalogue's tier column",
        "`core` or `extended` region of the conservative catalogue at the 3' "
        "breakend; empty when that end is outside one."),
    "cfs_narrow_status": ("annotation", "genomic interval", "annotation — does NOT filter", "cfs_narrow_bp1 + cfs_narrow_bp2",
        "`none` / `bp1` / `bp2` / `both`. Kept per end rather than collapsed to a "
        "boolean: one end deep inside a fragile site and the other far outside is a "
        "different claim from both ends inside."),
    "cfs_broad_bp1": ("annotation", "genomic interval", "annotation",
        "permissive CFS catalogue, GRCh38",
        "Fragile site containing breakend 1 under the PERMISSIVE catalogue. **That "
        "catalogue covers ~65% of the genome**, so a hit is close to what chance "
        "alone produces and is not evidence on its own."),
    "cfs_broad_bp2": ("annotation", "genomic interval", "annotation", "permissive CFS catalogue",
        "Fragile site containing the 3' breakend under the permissive catalogue, "
        "whose ~65% genome footprint makes a hit close to the chance expectation."),
    "cfs_broad_status": ("annotation", "genomic interval", "annotation — does NOT filter", "cfs_broad_bp1 + cfs_broad_bp2",
        "`none` / `bp1` / `bp2` / `both` under the permissive catalogue."),
    "cfs_agreement": ("annotation", "genomic interval", "annotation — does NOT filter", "the two status columns",
        "`both` / `narrow_only` / `broad_only` / `neither`. **`broad_only` is the "
        "weak class**: a hit the conservative catalogue does not support, and one a "
        "permissive catalogue's genome footprint makes near-certain by chance. Always "
        "test a hit rate against the admitted junctions that produced NO candidate — "
        "a uniform-genome null is wrong here, because these breakends sit in genes by "
        "construction. Expect the permissive catalogue to discriminate poorly."),

    # -- MHC layer, prediction ---------------------------------------------
    "presentable": ("MHC layer", "prediction — netMHCpan", "verdict — does NOT filter",
        "netMHCpan 4.2e against the line's own class I genotype",
        f"**Could this cell present the peptide?** Binds >= 1 allele of the line's OWN "
        f"genotype. Binder = IC50 <= {BT['affinity_nM']:g} nM AND %Rank_BA <= "
        f"{BT['rank_BA']:g} AND %Rank_EL <= {BT['rank_EL']:g}, all three at once — "
        "stricter than netMHCpan's own weak/strong convention. **Predicted binding, "
        "never observed presentation**: says nothing about proteasomal processing, TAP "
        "transport, surface abundance or T-cell recognition."),
    "n_alleles_binding": ("MHC layer", "prediction — netMHCpan", "measurement", "netMHCpan",
        "How many of the cell line's own class I alleles bind this peptide, out of the "
        "six in its genotype. 0 means `presentable` is False. A peptide binding one "
        "allele is presentable by that cell just as surely as one binding several; the "
        "count is breadth, not confidence — for confidence read "
        "`presentation_margin`."),
    "binding_alleles": ("MHC layer", "prediction — netMHCpan", "annotation", "netMHCpan", "Which ones, `;`-separated."),
    "presentation_best_allele": ("MHC layer", "prediction — netMHCpan", "annotation", "netMHCpan",
        "Which of the line's own alleles gives the largest `presentation_margin`, "
        "written as `HLA-A02:01` with no asterisk. This is the allele the robustness "
        "bin refers to; a peptide may bind others less comfortably."),
    "presentation_margin": ("MHC layer", "derived from prediction", "measurement", "distance to the three binder cuts",
        "**How much room is left before the call flips.** Distance to each cut as a "
        "fraction of that cut, the SMALLEST of the three (the constraint holding the "
        "call up), then the LARGEST across binding alleles. 0 = exactly on a boundary. "
        "Scales are mixed on purpose, so it ranks fragility rather than estimating a "
        "flip probability."),
    "presentation_limiting_cut": ("MHC layer", "derived from prediction", "annotation", "argmin of the three distances",
        "Which of the three cuts is closest to failing, and therefore the one holding "
        "the call up: `affinity_nM`, `rank_BA` or `rank_EL`. Empty for non-binders. "
        "Among marginal calls it is spread across all three rather than concentrated "
        "in one, so the fragility comes from the conjunction, not a single badly "
        "placed threshold."),
    "presentation_robustness": ("MHC layer", "derived from prediction", "verdict — does NOT filter", "presentation_margin, binned",
        "One of four labels, from `presentation_margin`: `flippable` (within 10% of a "
        "cut), `marginal` (within 25%), `solid` (within 50%), `robust` (further than "
        "50% inside every cut). Empty for non-binders. **Not a stricter biological "
        "claim** — it is the identical binder criterion minus the calls that a re-run, "
        "a predictor version change or a threshold nudge would flip."),

    # -- cohort-side verdicts, given their own entries so they do not fall
    # -- through to the generic "populated only where matched" text ---------
    "autologous_verdict": ("stages 3-5 (cohort)", "prediction — netMHCpan", "annotation — does NOT filter",
        "the cohort-side binder test over a catalogue patient's genotype",
        "One of four values, and the distinction between them is the point: `binder` "
        "(a patient carrying this peptide has an allele it binds), `non_binder` "
        "(tested and it does not), `unevaluable_no_typed_patient` (every patient "
        "carrying it lacks a linkable genotype), or `unevaluable_allele_outside_panel` "
        "(the carrying patient has an allele that was never predicted). "
        "**`unevaluable` is not `non_binder`** — the binder count is a floor, and the "
        "size of the gap is exactly these rows."),
    "autologous_limiting_cut": ("stages 3-5 (cohort)", "derived from prediction", "annotation — does NOT filter",
        "argmin of the three distances, cohort side",
        "Which of the three cuts holds the cohort-side call up: `affinity_nM`, "
        "`rank_BA` or `rank_EL`. Empty for non-binders and unevaluables."),
    "autologous_robustness": ("stages 3-5 (cohort)", "derived from prediction", "annotation — does NOT filter",
        "the cohort-side margin, binned",
        "`flippable`, `marginal`, `solid` or `robust`, binning the cohort-side margin "
        "at 10%, 25% and 50% of a cut. Empty unless the verdict is `binder`."),
    "panel_limiting_cut": ("stages 3-5 (cohort)", "derived from prediction", "annotation — does NOT filter",
        "argmin of the three distances, panel side",
        "Which cut holds up the PANEL-side call — the union over every cohort allele, "
        "an upper bound rather than a result: `affinity_nM`, `rank_BA` or `rank_EL`."),
    "gene_concordant": ("stages 3-5 (cohort)", "catalogue comparison", "annotation — does NOT filter",
        "this sample's gene vs the catalogue's for the matched peptide",
        "Whether the peptide is broken in the same gene on both sides of the match. "
        "An identical sequence arising from a different gene is convergence rather "
        "than a shared event, so this separates the two. Recorded, not enforced, in "
        "this table."),

    # -- cohort comparison --------------------------------------------------
    "matched_reference": ("stages 3-5", "catalogue comparison", "verdict — does NOT filter",
        "exact sequence identity against the reference catalogue",
        "**Does this exact peptide also occur in the reference cohort?** A sequence "
        "identity test, HLA-independent by construction. **In this table it is an "
        "annotation: the universe is not selected on it.** When False every cohort-side "
        "column is empty because the question was not asked, not because the answer was "
        "no."),
}

#: Suffixed patient-side columns share their base meaning.
PATIENT_SUFFIX = ("_patient",)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--universe", type=pathlib.Path, required=True)
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    table = pd.read_csv(args.universe, sep="\t", dtype=str, low_memory=False)
    rows = []
    for column in table.columns:
        stage, level, kind, source, rule = SPEC.get(column, ("", "", "", "", ""))
        meaning = column_meanings.resolve(column) or ""
        if not stage:
            base = column
            for suffix in PATIENT_SUFFIX:
                if column.endswith(suffix):
                    base = column[:-len(suffix)]
            if base != column and base in SPEC:
                stage, level, kind, source, rule = SPEC[base]
                rule = (rule + " HERE it is the COHORT-side verdict, computed over a "
                               "catalogue patient's genotype — not this line's. See "
                               "`presentable` for the line's own.")
            elif column.startswith(("n_patients", "panel_", "presenting_alleles",
                                    "n_alleles_panel", "binds_panel", "ref_",
                                    "autologous_", "binds_autologous",
                                    "gene_concordant", "svtype_concordant",
                                    "n_rows", "n_patients")):
                stage, level = "stages 3-5 (cohort)", "catalogue comparison"
                kind, source = "annotation — does NOT filter", "the reference catalogue"
                rule = ("Cohort-side column, populated only where "
                        "`matched_reference` is True. Empty means the question "
                        "was not asked, not that the answer was no.")
            else:
                stage, level, kind = "", "", "annotation"
                source, rule = "", "No threshold."
        filled = int(table[column].notna().sum())
        rows.append({
            "column": column,
            "stage": stage,
            "level": level,
            "kind": kind,
            "derived_from": source,
            "what_it_is_and_threshold": rule,
            "meaning": meaning,
            "n_populated": filled,
            "pct_populated": round(100.0 * filled / len(table), 1),
            "n_distinct": int(table[column].nunique(dropna=True)),
            "example": next((str(v) for v in table[column].dropna().head(1)), ""),
        })

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)

    missing = out[out["meaning"] == ""]["column"].tolist()
    print(f"  {len(out)} columns -> {args.out}")
    print(f"  without a meaning: {missing or 'none'}")
    print(f"  with an explicit threshold: "
          f"{int((~out.what_it_is_and_threshold.str.startswith('No threshold')).sum())}")
    print("\n  by level:")
    for level, n in out["level"].value_counts().items():
        print(f"    {level or '(unset)':<34} {n:>3}")
    print("  by kind:")
    for kind, n in out["kind"].value_counts().items():
        print(f"    {kind:<34} {n:>3}")


if __name__ == "__main__":
    main()
