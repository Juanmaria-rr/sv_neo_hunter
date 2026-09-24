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

    # -- IGV navigation ------------------------------------------------------
    "igv_locus": ("navigation", "bookkeeping", "annotation", "chrom/pos, written in the BAM's contig naming",
        "**Paste this into IGV's search box.** Covers the whole junction, padded by "
        "a quarter of the event size each side (floor 200 bp, ceiling 5 kb) so the "
        "flanks that make a depth change visible are on screen. For an "
        "inter-chromosomal junction it holds two loci separated by a space, which "
        "opens IGV's split view. Written as `chr1:…` because the alignments name "
        "their contigs that way while this table uses `1` — IGV searches the genome "
        "the BAM was aligned to."),
    "igv_bp1": ("navigation", "bookkeeping", "annotation", "chrom1/pos1 with a fixed window",
        "The 5' breakend alone, 150 bp each side. Use it when the whole-junction "
        "window is too zoomed out for IGV to draw individual reads."),
    "igv_bp2": ("navigation", "bookkeeping", "annotation", "chrom2/pos2 with a fixed window",
        "The other end of the junction on its own, 150 bp either side. The pair "
        "`igv_bp1`/`igv_bp2` exists because a junction has two ends and they can be "
        "megabases apart or on different chromosomes: when that is the case no single "
        "window shows both, and each has to be looked at in turn."),

    # -- lineage -----------------------------------------------------------
    "acquired_in": ("lineage", "bookkeeping", "identity", "the run each candidate table came from",
        "The line whose OWN variant calls produce this peptide — where it entered the "
        "lineage. Semicolon-separated if two lines generate the same sequence "
        "independently. **Filter on this to ask what a knockout added, and to stay "
        "inside a line-restricted subset**; filtering on `present_in` instead pulls "
        "back every line that merely inherits it. No threshold."),
    "present_in": ("lineage", "bookkeeping", "identity", "acquired_in + the configured parent/child lineage",
        "Every line whose cells carry this peptide, INCLUDING the ones that merely "
        "inherited it from a parent rather than generating it themselves. So a peptide "
        "that arose in the parental line appears here for the parent and for every "
        "clone derived from it. This is the wider of the two lineage columns: use it "
        "to ask 'which lines have it', and `acquired_in` to ask 'which line made it'. "
        "Counting candidates with this column double-counts inheritance. No threshold."),
    "n_lines_present": ("lineage", "bookkeeping", "measurement", "count of present_in",
        "How many of the related cell lines in this run carry the peptide. 1 means it "
        "is confined to a single branch of the family tree — the interesting case, "
        "since it was acquired somewhere specific; the maximum equals the number of "
        "lines and is reached by anything the common ancestor already had, which is "
        "the least interesting case. **It counts LINES, not independent events**: a "
        "peptide present in four lines because one ancestor acquired it happened "
        "once, not four times. No threshold."),

    # -- stage 1, the SV call itself (DNA) ---------------------------------
    "sv_id": ("stage 1", "DNA — caller", "identity", "the SV VCF",
        "The identifier the variant caller gave this junction in its own output. It "
        "is a label for tracing a row back to the caller's file, carries no meaning "
        "of its own, and is NOT stable between runs or callers — do not join tables "
        "on it across pipeline versions."),
    "chrom1": ("stage 1", "DNA — caller", "identity", "the SV VCF",
        "Chromosome holding the FIRST of the junction's two ends. Written without "
        "the `chr` prefix (`1`, not `chr1`), while the alignment files use the "
        "prefix — which is why `igv_locus` exists rather than pasting this straight "
        "into a genome browser."),
    "pos1": ("stage 1", "DNA — caller", "identity", "the SV VCF",
        "Base position of the first end along `chrom1`, counting from 1 as the VCF "
        "does (not from 0, as most programming does — off by one against a BED file "
        "or a Python slice)."),
    "chrom2": ("stage 1", "DNA — caller", "identity", "the SV VCF",
        "Chromosome of the 3' breakend, without a `chr` prefix. Differing from "
        "`chrom1` makes the junction inter-chromosomal, which forces `test` to "
        "`chimeric` and leaves `event_size` undefined."),
    "pos2": ("stage 1", "DNA — caller", "identity", "the SV VCF",
        "Base position of the second end, 1-based like `pos1`. For a deletion on one "
        "chromosome the stretch BETWEEN `pos1` and `pos2` is what is missing, and "
        "its length is `event_size`; when the two ends are on different chromosomes "
        "there is no interval between them and that idea does not apply."),
    "gene1": ("stage 1", "annotation", "annotation", "Ensembl, by overlap with breakend 1",
        "The gene the first breakend lands in, looked up in the reference gene "
        "annotation. It contributes the FIRST part of the fusion protein, so its "
        "reading frame is the one the second partner is read against. Empty means "
        "the breakend fell outside any annotated gene, not that the lookup failed. "
        "**This is annotation, not measurement**: it says where the breakpoint is, "
        "never that the gene is expressed or that the fusion is made."),
    "gene2": ("stage 1", "annotation", "annotation", "Ensembl, by overlap with breakend 2",
        "Gene overlapping the 3' breakend, from the reference annotation rather than "
        "from any measurement. It contributes the downstream half of the fusion "
        "protein, so a frameshift makes every residue it contributes novel."),
    "strand1": ("stage 1", "annotation", "annotation", "Ensembl transcript strand",
        "Which of the DNA's two strands the first partner's gene is read from, `+` "
        "or `-`. It matters because a fusion only makes a coherent protein when the "
        "two partners' orientations are compatible; the pair `strand1`/`strand2` is "
        "what tells them apart. Carried so results can be split by strand pair, not "
        "used to filter."),
    "strand2": ("stage 1", "annotation", "annotation", "Ensembl transcript strand",
        "Strand of the 3' transcript, `+` or `-`. Emitted so results can be "
        "stratified by strand pair; note that per-peptide proportions by strand are "
        "dominated by a few prolific loci, so stratify per event."),
    "svtype": ("stage 1", "DNA — caller", "annotation", "the SV VCF",
        "What kind of rearrangement the caller decided this is: `DEL` a deleted "
        "stretch, `DUP` a duplicated one, `INV` a piece flipped end to end, `BND`/"
        "`TRA` two distant or different chromosomes joined. **It is the caller's "
        "classification, not an independent measurement**, and the pipeline does not "
        "re-derive it — the geometry actually used for read counting is `test`, "
        "which is computed from the coordinates and can disagree with this label."),

    # -- stage 2, the predicted fusion protein (IN SILICO) -----------------
    "neopeptide": ("stage 2", "protein — in silico", "identity",
        "translation of the reconstructed fusion transcript",
        "The candidate peptide sequence itself, in one-letter amino-acid code — the "
        "thing the whole table is about. It is obtained by reconstructing what the "
        "fused transcript would translate to and sliding a window along it, so "
        "**every sequence here is a PREDICTION**: it was computed from breakpoint "
        "coordinates plus the reference genome, and no read, peptide measurement or "
        "experiment shows it exists. One junction yields many overlapping windows, "
        "so counting rows counts windows, not events."),
    "pep_length": ("stage 2", "protein — in silico", "measurement", "len(neopeptide)",
        f"How many amino acids long the peptide is. Only {C.PEPTIDE_LENGTHS} are "
        f"generated (`PEPTIDE_LENGTHS`) because the MHC class I groove is closed at "
        f"both ends and physically holds roughly that many residues; longer or "
        f"shorter sequences are not presented by this pathway and would be wasted "
        f"predictions. Since every length is tried at every position, one junction "
        f"produces several rows differing only by a residue at one end — they are "
        f"not independent candidates."),
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
        "A low count deserves a look, but RNA-to-DNA ratios do not settle anything "
        "here — see `vf_bp2`."),
    "vf_bp2": ("stage 6", "DNA — read evidence", "measurement", "caller VF field",
        f"DNA fragments supporting the 3' breakend. Enters `sv_hc` at >= {C.HC_MIN_VF}. "
        "Low support here is worth noticing, but do NOT read it as a ratio against "
        "`junction_reads`: RNA depth at an expressed gene exceeds DNA depth by orders "
        "of magnitude whether or not the lesion is real. Check DNA depth across the "
        "interval instead."),
    "qual_bp1": ("stage 6", "DNA — read evidence", "measurement", "caller QUAL",
        f"How confident the variant caller is in the first breakend, on a Phred-like "
        f"scale where higher is better. It reflects how well the read evidence fits a "
        f"real junction — **not** how biologically interesting the event is, and not "
        f"whether the junction is ordinary human variation (that is `pon_count`). "
        f"Enters `sv_hc` at >= {C.HC_MIN_QUAL}."),
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
        "pon_count vs PON_MAX", f"True when this junction was seen in fewer than {C.PON_MAX} of the "
        f"normal genomes in the panel (`pon_count < {C.PON_MAX}`) — that is, when it "
        f"is rare enough not to look like ordinary human variation. **On the noPON "
        f"branch it is recorded but deliberately left out of `is_private`**, so a row "
        f"can be False here and still be kept: this table reports the panel rather "
        f"than filtering on it, which is what lets you decide the threshold "
        f"afterwards instead of inheriting one."),
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
        f"True when the junction is rarer in the general population than "
        f"{C.GNOMAD_MAX_AF} (`GNOMAD_MAX_AF`), using the frequency in a large public "
        f"database of human variation. A common variant is something the person was "
        f"born with, not something the tumour or the cell line acquired. **Careful: "
        f"it is also True when no frequency was found at all** — absence of evidence "
        f"is scored as a pass, so a True here can mean 'rare' or 'never looked up'. "
        f"Check `gnomad_af_used` before believing it."),
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
        f"How strongly the first partner gene is transcribed in this sample, in "
        f"transcripts per million — normalised for gene length and library size, so "
        f"it is comparable across genes and samples as a raw read count is not. It "
        f"describes the GENE, not the junction: a gene can be highly expressed while "
        f"the altered transcript is never made. `expressed` needs >= "
        f"{C.TPM_EXPRESSED} TPM here and at the partner."),
    "gene2_TPM": ("stage 7", "RNA — quantification", "measurement", "Isofox transcript quantification",
        f"Expression of the 3' gene in this sample's RNA, in transcripts per million. "
        f"`expressed` requires >= {C.TPM_EXPRESSED} TPM here AND at the 5' gene."),
    "min_side_TPM": ("stage 7", "RNA — quantification", "measurement", f"{C.COVERAGE_SUMMARY}(gene1_TPM, gene2_TPM)",
        "The LOWER of `gene1_TPM` and `gene2_TPM` — expression of whichever partner "
        "gene is the quieter of the two. It is the minimum and never the maximum "
        "(`COVERAGE_SUMMARY`) because a fusion transcript needs both halves: taking "
        "the maximum would report a locus as transcribed on the strength of its "
        "partner alone, while the side that actually limits the product is silent."),
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
    "junction_alt_reads": ("stage 8", "RNA — direct read counting", "measurement",
        "the RNA BAM: gapped reads over the SV interval, minus `junction_reads`",
        "**Reads that cross the same interval WITHOUT using the candidate junction** "
        "— the alternative, which for a deletion is normally the canonical intron. "
        "This column DELIBERATELY EXCLUDES `junction_reads`: it is the count of the "
        "competing outcome on its own. **It is not the denominator of any percentage "
        "here** — dividing by it gives odds, not a share, and is not bounded by 100%. "
        "It is reported so the two outcomes can be read side by side. A read is "
        "included when it carries a CIGAR `N` gap overlapping the SV interval "
        "(`gap_start < pos2 and gap_end > pos1`) whose size or position does NOT "
        "match the event. Reads with no gap are excluded throughout: a read lying "
        "entirely inside an exon never had the opportunity to use either junction, "
        "and counting non-voters in a vote share measures exon length, not biology."),
    "junction_interval_total": ("stage 8", "RNA — direct read counting", "measurement",
        "`junction_reads` + `junction_alt_reads`",
        "**Every read that crosses the SV interval with a gap, the candidate junction "
        "INCLUDED.** The word `total` is load-bearing: unlike `junction_alt_reads` "
        "this one CONTAINS `junction_reads`, which is precisely what makes it a valid "
        "denominator — a percentage requires the numerator to be inside it. It is the "
        "population of reads that could have gone either way, so a share of it reads "
        "as 'what fraction of transcripts through this interval used the lesion'. "
        "Not to be confused with total coverage at the locus, nor with all gapped "
        "reads in the gene: gaps elsewhere (a neighbouring intron spliced in every "
        "transcript regardless) are excluded because they do not compete."),
    "junction_usage": ("stage 8", "RNA — direct read counting", "measurement",
        "`junction_reads` / `junction_interval_total`",
        "**The fraction of transcripts crossing this interval that use the candidate "
        "junction**, 0-1. Its purpose is to normalise for expression, which the raw "
        "count cannot: 7,000 junction reads mean one thing in a gene at 900 TPM and "
        "another in a gene at 5 TPM, and any threshold on the raw count is in effect "
        "a filter on expression. Roughly 0.5 is what a heterozygous deletion "
        "transcribed as readily as the intact allele would give; well below that is "
        "allelic imbalance or a less stable transcript; a few per cent is the shape "
        "of alignment noise in a highly expressed gene. **Limitation:** only defined "
        "where `test` == `sizegap`. A chimeric junction has no canonical alternative "
        "competing for the same interval and an insertion is not a choice between two "
        "splice outcomes, so the column is EMPTY for those — empty meaning "
        "not applicable, never 0, which would read as 'the junction is never used'. "
        "**The denominator is the confidence, so read it too.** A usage of 100% over "
        "a `junction_interval_total` of 15 says only that nothing competes at that "
        "locus — there is no alternative to measure against, so the figure is "
        "trivially high and means nothing. Observed here: the junctions above 65% "
        "have a median denominator of 15 and almost all have `junction_alt_reads` "
        "of 0-2, while those below 5% have a median denominator of 209. A high "
        "share is only interpretable when there was something to lose it to. "
        "In this respect the metric fails the same way the raw count does, in the "
        "opposite direction: neither column is self-sufficient. "
        "It does NOT filter and it is not a measure of protein abundance: use it "
        "alongside `junction_reads`, which stays the quotable evidence that the "
        "altered transcript exists at all."),
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
        "threshold — a diagnostic, and one that is easy to over-read. A small "
        "distance is NOT evidence against the junction: a transcript from a "
        "deleted allele, aligned to a reference that still carries the deleted "
        "bases, looks exactly like a novel intron, so a splice-aware tool calls "
        "one whether or not the lesion is genomic. It is the expected observation "
        "either way. **Only DNA depth across the interval separates the two** — "
        "see `tools/validate_junction.py`."),
    "alt_sj_frags": ("stage 7", "RNA — Isofox splice calls", "measurement", "Isofox",
        "Fragments Isofox counts behind that alternative junction. Do NOT compare "
        "it against `vf_bp1` as a ratio: RNA depth at an expressed gene exceeds "
        "DNA depth by orders of magnitude whether or not the lesion is real, so "
        "the ratio measures expression rather than artefact."),
    "alt_sj_type": ("stage 7", "RNA — Isofox splice calls", "annotation", "Isofox",
        "How the splice-aware tool Isofox classified the gap it sees here, e.g. "
        "`NOVEL_INTRON` for a splice junction not in the reference annotation. "
        "**A `NOVEL_INTRON` call is NOT evidence against the SV**: a transcript from "
        "a deleted allele, aligned to a reference that still contains the deleted "
        "bases, looks exactly like a novel intron, so the call is what you would see "
        "either way and cannot tell the two apart. Only DNA depth across the interval "
        "can. This column comes from a different tool than `junction_reads` and is "
        "allowed to disagree with it."),
    "alt_sj_within_window": ("stage 7", "RNA — Isofox splice calls", "verdict — does NOT filter", "Isofox",
        "True when Isofox reports some alternative splice junction close enough to "
        "this breakpoint to be confusable with it. It is a flag to LOOK, not a "
        "verdict: the SV and a nearby splice site can coexist, and the fact that "
        "both are in the window is exactly the situation where read counts alone "
        "cannot attribute the gap. Does not filter."),
    "retained_intron": ("stage 7", "RNA — Isofox splice calls", "measurement", "Isofox",
        "How many retained-intron events Isofox reports at this locus — places where "
        "an intron was NOT spliced out and stayed in the transcript. It matters "
        "because a retained intron changes what the transcript actually contains, so "
        "the protein this pipeline predicted from the spliced sequence may not be the "
        "one made. A count of the locus, not of this peptide; empty means Isofox was "
        "not run or reported nothing here, not zero retention."),
    "isofox_fusion": ("stage 7", "RNA — Isofox splice calls", "annotation", "Isofox",
        "The fusion Isofox called at this locus from the RNA, if it called one — a "
        "SECOND, independent opinion, since Isofox works from transcripts while "
        "stage 1 works from DNA. Agreement is reassuring; **disagreement is not "
        "decisive either way**, because the two tools are answering slightly "
        "different questions and Isofox only sees what is transcribed. Empty means "
        "no call here, which for an unexpressed gene is expected."),
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
        f"Fragments that support the junction by being split across BOTH of its "
        f"ends: the aligner matched part of the read at one breakend and recorded "
        f"the rest as a second (supplementary) alignment, which landed within "
        f"{C.SA_PARTNER_TOLERANCE} bp of the partner. This is the counting mechanism "
        f"for junctions joining places too distant to describe as a gap — two "
        f"different chromosomes, say — where no single read can carry one. Which "
        f"mechanism is allowed to fire is fixed by `test` before any read is "
        f"examined, so a zero here is usually 'not applicable', not 'searched and "
        f"found nothing'."),
    "junction_by_insert": ("stage 8", "RNA — direct read counting", "measurement",
        "pysam: inserted sequence",
        f"Fragments that support the junction by carrying the INSERTED sequence "
        f"itself — bases present in the read that are absent from the reference, at "
        f"least {C.MIN_INSERT_LEN} bp of them. This is the counting mechanism for "
        f"insertions, one of three; which one is allowed to fire is decided by `test` "
        f"before any read is looked at, so a zero here usually means this mechanism "
        f"was not applicable, not that the search failed."),
    "coverage_bp1": ("stage 8", "RNA — direct read counting", "measurement",
        f"alignments within {C.COVERAGE_WINDOW} bp of breakend 1",
        "**Context, not evidence.** A breakpoint inside a highly expressed gene has "
        "thousands of reads whether or not the junction exists, which is why "
        "coverage never tiers an event."),
    "coverage_bp2": ("stage 8", "RNA — direct read counting", "measurement",
        f"alignments within {C.COVERAGE_WINDOW} bp of breakend 2",
        "How many RNA reads sit near the second breakend. Like `coverage_bp1` this "
        "is CONTEXT — how much signal was available to detect anything at all — and "
        "never evidence: a breakpoint inside a busy gene has deep coverage whether "
        "or not the junction is real. Only reads CROSSING the junction are evidence, "
        "and coverage deliberately never tiers an event."),
    "min_coverage": ("stage 8", "RNA — direct read counting", "measurement", "min of the two",
        "The LOWER of `coverage_bp1` and `coverage_bp2` — RNA read depth at "
        "whichever end of the junction has less of it. The minimum and never the "
        "maximum, because detecting a junction needs signal at BOTH ends: the "
        "maximum would report a locus as well covered when one side is silent and "
        "nothing could have been seen there. Still context, not evidence."),
    "softclip_bp1": ("stage 8", "RNA — direct read counting", "measurement",
        f"soft clips >= {C.MIN_SOFTCLIP_LEN} bp at breakend 1",
        f"**Diagnostic only — soft clips never tier an event** "
        f"(`SOFTCLIPS_TIER_EVENTS = {C.SOFTCLIPS_TIER_EVENTS}`). A clip shows a read "
        "ENDS at the breakpoint; it says nothing about where the rest went. Clipped "
        "reads whose supplementary alignments landed megabases away once produced a "
        "candidate that did not survive."),
    "softclip_bp2": ("stage 8", "RNA — direct read counting", "measurement",
        f"soft clips >= {C.MIN_SOFTCLIP_LEN} bp at breakend 2",
        "Reads whose alignment stops abruptly at the second breakend, leaving a tail "
        "the aligner set aside. Clips PILE UP at a real breakpoint, so this is a "
        "useful hint — but only a hint: a clip shows a read ENDS there and says "
        "nothing about where the rest of it belongs. Clipped reads whose remainder "
        "turned out to map megabases away once produced a candidate that did not "
        "survive, which is why soft clips never tier an event here."),
    "alignments_bp1": ("stage 8", "RNA — direct read counting", "measurement", "pysam",
        "How many alignment records were looked at near the first breakend before "
        "any quality filtering. It is the DENOMINATOR for `low_mapq_bp1`: on its own "
        "the number of discarded reads means nothing, since 100 discarded out of 120 "
        "is a broken locus and 100 out of 50,000 is routine. Records, not fragments, "
        "so one molecule may appear more than once here."),
    "alignments_bp2": ("stage 8", "RNA — direct read counting", "measurement", "pysam",
        "How many alignment records were looked at near the second breakend before "
        "quality filtering — the denominator for `low_mapq_bp2`, for the same reason "
        "`alignments_bp1` is for the first end: a count of discarded reads is "
        "uninterpretable without knowing how many there were to start with."),
    "low_mapq_bp1": ("stage 8", "RNA — direct read counting", "measurement", "pysam",
        f"How many alignments near the first breakend were thrown out because the "
        f"aligner was not confident where they came from (MAPQ < {C.MIN_READ_MAPQ}). "
        f"**Read it as a SHARE of `alignments_bp1`, never on its own.** A high share "
        f"means the sequence here occurs in several places in the genome, so reads "
        f"cannot be assigned to this locus and everything measured at this end is "
        f"doubtful — including a junction count of zero, which may mean 'could not "
        f"look' rather than 'looked and found nothing'."),
    "low_mapq_bp2": ("stage 8", "RNA — direct read counting", "measurement", "pysam",
        f"How many alignments near the second breakend were discarded for low "
        f"mapping confidence (MAPQ < {C.MIN_READ_MAPQ}). As with `low_mapq_bp1`, "
        f"divide by `alignments_bp2` before reading it: a large raw count at a deeply "
        f"covered locus can be a small and unremarkable fraction."),

    # -- lesion size --------------------------------------------------------
    "span": ("stage 1", "DNA — caller", "measurement", "|pos2 - pos1|",
        "How far apart the junction's two ends are along the chromosome, in bases — "
        "simply `pos2` minus `pos1`. **It is NOT the size of the lesion**, and using "
        "it as one is a mistake this pipeline has actually made: where bases are "
        "INSERTED at the junction the two ends can be almost adjacent while the "
        "change is large, and a 34 bp insertion was once reported as a 2 bp deletion "
        "for exactly that reason. `event_size` is the size-aware figure; this column "
        "is kept beside it so the two can be compared. Undefined when the ends are on "
        "different chromosomes, where the distance has no meaning."),
    "event_size": ("stage 1", "DNA — derived", "measurement", "type-aware size function",
        f"Deleted, duplicated or inserted length as the SV type implies. Empty for "
        f"BND/TRA/SGL, where size is not meaningful. Events below "
        f"{C.MIN_TESTABLE_GAP_SIZE} bp are UNTESTABLE for RNA junction evidence."),
    "insert_len": ("stage 1", "DNA — caller", "measurement", "the VCF's inserted sequence",
        "How many bases are INSERTED at the seam — sequence present in the sample "
        "but absent from the reference, added rather than removed. 0 for a plain "
        "deletion. It is the reason `span` and `event_size` can disagree: a junction "
        "whose two ends are 2 bp apart but which carries 34 inserted bases is a 34 bp "
        "change, not a 2 bp one. It also decides how the RNA evidence can be counted "
        "at all, since an insertion is found by matching the inserted sequence inside "
        "reads rather than by looking for a gap."),

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
        "How central this fragile site is under the conservative catalogue, at the "
        "first breakend: `core` (those regions cover 0.5% of the genome) or "
        "`extended` (0.7%). The percentages are the point — they are the rate at "
        "which a breakpoint dropped at random would land here anyway, so they are "
        "the baseline any claim about fragile sites has to beat. Empty means this "
        "end is outside the catalogue."),
    "cfs_narrow_tier_bp2": ("annotation", "genomic interval", "annotation", "the catalogue's tier column",
        "The same tiering as `cfs_narrow_tier_bp1` but for the second breakend: "
        "`core` or `extended` under the conservative catalogue, empty when that end "
        "falls outside it. The two ends are tiered separately because they can land "
        "on different sides of a boundary, and one end inside a fragile site with the "
        "other far outside is a different claim from both ends inside."),
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
        "Which fragile site the second breakend falls in under the PERMISSIVE "
        "catalogue. Common fragile sites are not points but regions of megabases, and "
        "published catalogues disagree about them by orders of magnitude — this one "
        "covers about 65% of the genome. **At that size a hit is very nearly what "
        "chance alone produces**, so a name here is close to meaningless on its own "
        "and only `cfs_agreement` makes it interpretable."),
    "cfs_broad_status": ("annotation", "genomic interval", "annotation — does NOT filter", "cfs_broad_bp1 + cfs_broad_bp2",
        "Which ENDS of the junction fall inside a fragile site under the permissive "
        "catalogue: `none`, `bp1` (only the first), `bp2` (only the second), or "
        "`both`. Kept as four values rather than one true/false because collapsing "
        "them would throw away the distinction that matters — both ends inside is a "
        "much stronger claim than one. Given that catalogue's ~65% genome footprint, "
        "expect `both` by chance most of the time; read it against `cfs_agreement`."),
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
    "binding_alleles": ("MHC layer", "prediction — netMHCpan", "annotation", "netMHCpan",
        "Which of this cell line's own class I alleles are predicted to bind the "
        "peptide, `;`-separated, written as `HLA-` then the gene, group and protein "
        "with no asterisk. Empty where `presentable` is "
        "False — no allele bound — so the two columns must agree. Everyone carries "
        "six class I alleles and each binds a different set of peptides, so a peptide "
        "presentable in one individual may be invisible in another; this column names "
        "the alleles responsible in THIS line. **Predicted binding, not observed "
        "presentation.**"),
    "presentation_best_allele": ("MHC layer", "prediction — netMHCpan", "annotation", "netMHCpan",
        "The single allele of this cell line's own six that is predicted to bind the "
        "peptide most comfortably — the one furthest from failing the binder cuts — "
        "written as `HLA-` plus gene, group and protein, without the asterisk some "
        "notations use. It names "
        "the allele that the robustness bin is about, so the two are read together; "
        "the peptide may bind others as well, less comfortably. **This is the best "
        "case, not a typical one**: summarising a row by its strongest allele "
        "flatters it, which is why `binding_alleles` lists them all."),
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

#: Terms the columns use that a reader outside this field has no reason to know.
#: Defined ONCE here and written beside the dictionary, rather than re-explained
#: in thirty entries or — worse — assumed. Ordered so an earlier term never
#: depends on a later one.
GLOSSARY: list[tuple[str, str]] = [
    ("structural variant (SV)",
     "A rearrangement that joins two pieces of a genome that are not normally "
     "adjacent: a deleted stretch, a duplicated one, an inversion, or two "
     "different chromosomes stuck together. Bigger than a single-letter mutation."),
    ("breakend",
     "One of the two ENDS an SV joins. Every junction has exactly two, so every "
     "positional column in this table comes in a pair: `chrom1`/`pos1` is one end "
     "and `chrom2`/`pos2` the other. They can be far apart, or on different "
     "chromosomes. `bp1` and `bp2` in a column name mean breakend 1 and 2."),
    ("junction",
     "The seam itself — the point where the two breakends are joined. The new "
     "sequence spanning that seam is what can make a protein the body has never "
     "seen, which is the whole reason this table exists."),
    ("5' and 3'",
     "Which way round the two partners sit along the gene. 5' ('five prime') is "
     "the start-of-the-gene side, 3' the end. In a fusion the 5' partner "
     "contributes the first half of the protein and the 3' partner the second."),
    ("read, fragment",
     "A read is one stretch of sequence the machine produced. Sequencing produces "
     "them in PAIRS from the same physical molecule, and the pair is a FRAGMENT. "
     "Counts here are per fragment, so one molecule is one observation however "
     "many records it left behind — otherwise a molecule would be counted twice."),
    ("coverage, depth",
     "How many reads sit over a given base. It is CONTEXT, not evidence: a "
     "breakpoint inside a busy gene has deep coverage whether or not the SV is "
     "real. What counts as evidence is reads crossing the junction."),
    ("MAPQ (mapping quality)",
     "How confident the aligner is that a read came from where it placed it, "
     "0-60. Low MAPQ means the sequence occurs in several places in the genome, "
     "so the position is a guess. Repetitive regions are where callers go wrong."),
    ("CIGAR, and the `N` gap",
     "A compact description of how a read lines up against the reference: so many "
     "matching bases, so many inserted, so many skipped. A skip is written `N`, "
     "and a read carrying one crossed something that is in the reference but NOT "
     "in this sample's transcript — either a spliced-out intron or a deletion."),
    ("supplementary alignment (SA)",
     "When one read matches two separate places, the aligner records the second "
     "as a supplementary alignment. That is the signature of a read crossing a "
     "junction too abrupt to describe as a gap — two chromosomes joined, say."),
    ("soft clip",
     "The aligner gave up on the tail of a read and set it aside. Clipped reads "
     "PILE UP near a breakpoint, which is suggestive, but a clipped read has not "
     "been shown to cross anything — which is why they never tier an event here."),
    ("splicing, intron, exon",
     "A gene is transcribed whole and then edited: introns are cut out and exons "
     "joined. The difficulty running through this table is that a deleted stretch "
     "of DNA and a spliced-out intron look IDENTICAL in RNA — both are an `N` gap "
     "against a reference that still contains the sequence."),
    ("reading frame, frameshift",
     "Protein is read three letters at a time. An SV that removes a number of "
     "bases not divisible by three shifts the frame, and every residue after the "
     "junction becomes novel — not just the ones straddling the seam."),
    ("TPM",
     "Transcripts Per Million: how strongly a gene is expressed, normalised for "
     "gene length and library size. Comparable across genes and samples, which "
     "raw read counts are not."),
    ("peptide, neopeptide",
     "A short fragment of protein, here 8-11 residues. A NEOpeptide is one the "
     "SV creates that no normal human protein contains — the candidate target."),
    ("HLA / MHC class I, allele",
     "The molecule that carries peptide fragments to the cell surface for T cells "
     "to inspect. Everyone has six class I versions (ALLELES), and each binds a "
     "different set of peptides, so whether a peptide can be shown at all depends "
     "on the individual's genotype."),
    ("IC50 (nM), %Rank",
     "Two outputs of a binding predictor. IC50 is predicted binding strength — "
     "LOWER is tighter. %Rank is where the peptide sits against a background of "
     "random peptides for that same allele — also lower is better — and it is the "
     "fairer of the two across alleles, because some alleles bind everything "
     "tightly and would otherwise dominate."),
    ("panel of normals (PON)",
     "A collection of genomes from people without the disease. A junction seen "
     "there is ordinary human variation rather than something this sample "
     "acquired, however clean the call looks."),
    ("this is a PREDICTION, not an observation",
     "Most protein-level and binding columns here are computed from coordinates "
     "and reference sequence. No read shows the peptide, and no experiment shows "
     "it reaching the cell surface. The RNA columns are the only direct "
     "measurements of whether the altered transcript is really made."),
]


#: Suffixed patient-side columns share their base meaning.
PATIENT_SUFFIX = ("_patient",)


def adds_nothing(rule: str, meaning: str) -> bool:
    """Is `meaning` already contained in `rule`?

    Not symmetric on purpose. The failure this catches is not two cells saying
    the same thing — it is a SHORTER, weaker cell sitting beside a fuller one,
    which leaves a reader deciding which to trust and usually reading the worse.
    So the test is whether nearly every word of `meaning` already appears in
    `rule`, not whether the two resemble each other.
    """
    A = set(meaning.lower().split())
    B = set(rule.lower().split())
    if not A:
        return True
    if len(A & B) / len(A) > 0.6:
        return True
    # A note less than half the length of the definition beside it is a legacy
    # summary of that definition, phrased differently. Keeping it re-creates the
    # problem in another form: the reader has to decide which cell is current,
    # and the shorter one is the one they read.
    return len(meaning.split()) * 2 < len(rule.split())


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
        # `meaning` comes from the shared module and `rule` from this file; on
        # many columns they say the same thing twice, which makes the table look
        # padded and leaves a reader guessing which cell is authoritative. Drop
        # it when it adds nothing.
        if meaning and adds_nothing(rule, meaning):
            meaning = ""

        filled = int(table[column].notna().sum())
        rows.append({
            "column": column,
            "stage": stage,
            "level": level,
            "kind": kind,
            "derived_from": source,
            "what_it_is": rule,
            "also_noted": meaning,
            "n_populated": filled,
            "pct_populated": round(100.0 * filled / len(table), 1),
            "n_distinct": int(table[column].nunique(dropna=True)),
            "example": next((str(v) for v in table[column].dropna().head(1)), ""),
        })

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)

    # The glossary is what makes the entries readable without re-explaining a
    # breakend in thirty of them. Written beside the dictionary rather than
    # inside it so the dictionary stays one row per column.
    glossary_path = args.out.parent / (args.out.stem + "_GLOSSARY.tsv")
    pd.DataFrame([{"term": term, "plain_english": text}
                  for term, text in GLOSSARY]).to_csv(
        glossary_path, sep="\t", index=False)

    missing = out[out["what_it_is"] == ""]["column"].tolist()
    print(f"  {len(out)} columns -> {args.out}")
    print(f"  {len(GLOSSARY)} glossary terms -> {glossary_path}")
    print(f"  second note kept only where it adds something: "
          f"{int((out['also_noted'] != '').sum())}/{len(out)}")
    print(f"  without a definition: {missing or 'none'}")
    print(f"  with an explicit threshold: "
          f"{int((~out.what_it_is.str.startswith('No threshold')).sum())}")
    print("\n  by level:")
    for level, n in out["level"].value_counts().items():
        print(f"    {level or '(unset)':<34} {n:>3}")
    print("  by kind:")
    for kind, n in out["kind"].value_counts().items():
        print(f"    {kind:<34} {n:>3}")


if __name__ == "__main__":
    main()
