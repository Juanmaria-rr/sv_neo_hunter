"""
column_meanings.py — what every column of the master tables means, in our use.

WHY THIS IS A FILE AND NOT A COMMENT
------------------------------------
The master tables are the artefact shared when someone wants to check a result
themselves, and a column whose meaning has to be guessed is worse than an absent
one: it will be guessed wrong, silently. `pon_count` looks like a frequency,
`coverage_bp1` looks like evidence, `sf_bp1` looks like it might be the split
reads that support the junction in RNA. None of those readings is right.

Each entry says what the value IS, and where it matters, what it is NOT. Columns
recording a threshold's verdict point at the constant in `criteria.py`, so the
number and its justification are one hop away.

`tools/build_master_table.py` emits these into `<table>_column_dictionary.tsv`
beside each table; `check()` below fails if a column has no entry, so a new
column cannot ship undocumented.
"""
from __future__ import annotations

MEANINGS: dict[str, str] = {

    # -- identity --------------------------------------------------------
    "sample": "Run directory this row came from, `<cell_line>_<branch>`. Unique "
              "per row set; matches a folder under `results/` holding the "
              "criteria manifest that produced it.",
    "cell_line": "The biological sample alone, without the branch suffix. Use "
                 "this to group across branches, or to match a row to a report, "
                 "which names the line only.",
    "branch": "Which thresholds produced the row. `PON10` enforces the panel of "
              "normals; `noPON` measures and reports it without filtering on it, "
              "at both admission and the population-frequency verdict.",
    "junction_key": "The two breakend record IDs joined by `|`. Identifies the "
                    "junction within this VCF; not stable across call sets.",
    "sample_sv_id": "ID of the first breakend of the pair, from the VCF `ID` "
                    "column. The key everything downstream joins on, and the key "
                    "events are deduplicated by.",
    "prefix": "Output-file prefix used by the peptide generator for this sample.",
    "sv_id": "SV identifier recovered from the VCF by coordinate, linking a "
             "peptide back to the call that produced it. NeoSV carries no "
             "identifier of its own.",

    # -- breakend coordinates and lesion ---------------------------------
    "chrom1": "Chromosome of breakend 1, without the `chr` prefix.",
    "pos1": "Position of breakend 1 (1-based, as in the VCF).",
    "chrom2": "Chromosome of breakend 2 (the mate).",
    "pos2": "Position of breakend 2 (1-based). With `pos1`, defines the span; "
            "the lesion size is `event_size`, not the difference.",
    "svtype": "SV class as called: DEL, DUP, INS, INV, or BND for anything "
              "inter-chromosomal.",
    "span": "Raw distance between the two breakends, `pos2 - pos1`. NOT the "
            "lesion size — see `event_size`. Empty for inter-chromosomal events.",
    "event_size": "The lesion, computed type-aware: DEL `span-1` bases deleted, "
                  "DUP the `span` duplicated, INS the inserted length, BND "
                  "undefined (empty). Using `span` instead made 248 insertions "
                  "look like 1 bp events.",
    "insert_len": "Bases inserted at the junction. For an insertion-driven event "
                  "this, not the span, is the lesion, and it decides which RNA "
                  "test applies (`MIN_INSERT_LEN`).",
    "insert_seq": "The inserted sequence itself. Used to test RNA reads for the "
                  "inserted bases when the geometry makes a gap test impossible.",
    "alt_bp1": "The raw VCF ALT of breakend 1, e.g. `A[chr1:123456[`. The "
               "bracket orientation encodes which side of each breakend is "
               "joined, which is what makes the junction orientable.",
    "filter": "The caller's own FILTER value. `PASS`; `PON` = the caller found "
              "the breakpoint in its panel of normals; `INFERRED` = deduced from "
              "a copy-number transition with no read support, so no junction "
              "sequence exists to translate.",

    # -- SV call support -------------------------------------------------
    "vf_bp1": "Variant fragments supporting breakend 1 in the DNA, from the "
              "caller's `VF`. Compared against `HC_MIN_VF`.",
    "vf_bp2": "Variant fragments supporting breakend 2 in the DNA (`VF`).",
    "sf_bp1": "Split fragments at breakend 1 (`SF`) — fragments whose alignment "
              "is broken across the junction, in the DNA. A component of `VF`, "
              "not RNA evidence.",
    "sf_bp2": "Split fragments at breakend 2 (`SF`), in the DNA.",
    "df_bp1": "Discordant fragments at breakend 1 (`DF`) — pairs whose mates map "
              "at an unexpected distance or orientation, supporting the junction "
              "without spanning it.",
    "df_bp2": "Discordant fragments at breakend 2 (`DF`).",
    "ref_bp1": "Reference-supporting fragments at breakend 1 (`REF`) — the "
               "denominator against which `vf_bp1` is a variant allele fraction.",
    "ref_bp2": "Reference-supporting fragments at breakend 2 (`REF`).",
    "qual_bp1": "The caller's quality score for breakend 1. Compared against "
                "`HC_MIN_QUAL`.",
    "qual_bp2": "The caller's quality score for breakend 2.",
    "segmapq_bp1": "Mapping quality of the assembled segment at breakend 1 "
                   "(`SEGMAPQ`). Low values mean repetitive sequence, which is "
                   "where callers make their mistakes. Against `HC_MIN_SEGMAPQ`.",
    "segmapq_bp2": "Mapping quality of the assembled segment at breakend 2.",
    "homseq_bp1": "Microhomology sequence at breakend 1 (`HOMSEQ`). Homology "
                  "makes the exact breakpoint ambiguous by its own length, so a "
                  "long value means the coordinate is approximate.",
    "homseq_bp2": "Microhomology sequence at breakend 2 (`HOMSEQ`).",
    "imprecise_bp1": "The caller flagged breakend 1 as `IMPRECISE`: the position "
                     "is an estimate, so a base-resolution junction peptide built "
                     "from it is unreliable.",
    "imprecise_bp2": "The caller flagged breakend 2 as `IMPRECISE`.",
    "sv_hc": "Our own high-confidence verdict, combining `HC_MIN_SEGMAPQ`, "
             "`HC_MIN_VF`, `HC_MIN_QUAL` and `HC_MIN_SV_SIZE`. **Reported, not "
             "enforced** — no event is dropped for failing it. Copy-number "
             "fields are excluded (`TRUST_COPY_NUMBER_FIELDS`).",

    # -- annotation ------------------------------------------------------
    "gene1": "Gene broken at breakend 1, from the peptide generator's transcript "
             "annotation (Ensembl 115). Empty when the breakend is intergenic.",
    "gene2": "Gene broken at breakend 2.",
    "transcript_id1": "Ensembl transcript used to build the 5'/3' side at "
                      "breakend 1. Peptide sequence depends on this choice.",
    "transcript_id2": "Ensembl transcript used at breakend 2.",
    "strand1": "Strand of `transcript_id1`. Emitted so results can be stratified "
               "by transcript strand.",
    "strand2": "Strand of `transcript_id2`.",
    "frameshift": "RETIRED COLUMN NAME — present only in tables built before "
                  "2026-09-07, and always EMPTY in them. This pipeline read "
                  "NeoSV's frame verdict under the wrong attribute name, so "
                  "`getattr` returned the default for every row and the column "
                  "shipped blank without erroring. Superseded by `frame_effect`; "
                  "if a table you hold has this column, its frame information "
                  "was never populated. Rebuild with "
                  "`tools/build_master_table.py`.",
    "frame_effect": "NeoSV's verdict on the reading frame: `In-frame`, "
                    "`Stop-gain`, `Stop-loss`, or `Start-loss`. The last is a "
                    "RELIABILITY WARNING, not a biological class — upstream "
                    "loses the start codon when the 5' CDS is empty or under "
                    "three residues, falls back to the next ATG, and states "
                    "such predictions are of low reliability. Decisive for any "
                    "peptide that does NOT span the junction, since only an "
                    "altered downstream frame makes such a peptide neo.",
    "junction_nt": "Nucleotide offset of the junction within the fusion "
                   "sequence: the length of everything 5' of the breakpoint, "
                   "inserted bases included.",
    "junction_aa": "Residue index of the junction in the fusion protein, "
                   "`junction_nt // 3`.",
    "spans_junction": "Whether this peptide actually crosses the breakpoint. "
                      "**Annotation only** — `REQUIRE_SPANS_JUNCTION` is False. "
                      "Applying it as a filter discarded 63 of 64 catalogue "
                      "matches, because a sliding window legitimately yields "
                      "peptides whose novelty comes from the reading frame "
                      "rather than from straddling the breakpoint. Empty when "
                      "the junction offset is unknown, so 'not determined' is "
                      "never reported as 'does not span'.",
    "neopeptide": "The candidate peptide sequence, after subtracting anything "
                  "already present in the wild-type protein of either transcript.",
    "pep_length": "Length of `neopeptide`, 8-11 (`PEPTIDE_LENGTHS`).",

    # -- how far this junction got ---------------------------------------
    "n_peptides_generated": "Candidate peptides this junction produced. 0 means "
                            "the junction hit no coding transcript, or produced "
                            "only wild-type sequence — the denominator that makes "
                            "every later count interpretable.",
    "n_peptides_matched": "How many of them are identical to a catalogue peptide. "
                          "Exact string match, no alignment, no similarity score.",
    "n_peptides_credible": "Matched peptides that are also gene-concordant, "
                           "high-complexity and absent from the normal proteome.",
    "catalogue_genes": "Gene(s) the catalogue attributes the matched peptide(s) "
                       "to, `;`-separated. Compared against BOTH breakend genes: "
                       "testing only the first lost 8 of 64 matches, all "
                       "translocations.",
    "gene_concordant": "Whether a catalogue gene matches either breakend's gene. "
                       "Coincidental peptide convergence would not preferentially "
                       "land in the same gene, so this separates a shared locus "
                       "from a shared accident.",
    "low_complexity": "Peptide flagged as compositionally trivial by entropy, "
                      "homopolymer run, distinct residues or single-residue "
                      "fraction (`LC_*` constants). Such a string matches by "
                      "composition rather than by descent.",
    "is_self": "Peptide found in the normal proteome, so not neo-anything. Empty "
               "when no proteome was configured — recorded as not evaluated, "
               "never as passed.",
    "credible": "Matched AND gene-concordant AND high-complexity AND non-self.",

    # -- was it seen in patients? ----------------------------------------
    "patient_evidence": "Strongest level of patient recurrence this junction "
                        "reaches: `identical_peptide` (shared consequence) > "
                        "`same_gene_and_svtype` (shared mechanism) > `same_gene` "
                        "> `breakpoint_within_1kb|10kb|50kb|100kb` (shared locus "
                        "only) > `not_seen_in_patients`. The levels are different "
                        "kinds of evidence, not degrees of one.",
    "patient_bp_dist_bp": "Distance in bp from this junction to the nearest "
                          "patient breakpoint — the smaller of the two breakend "
                          "distances. Unthresholded, so any window can be applied.",
    "patient_bp_dist_bp1": "Distance from breakend 1 to the nearest patient "
                           "breakpoint.",
    "patient_bp_dist_bp2": "Distance from breakend 2 to the nearest patient "
                           "breakpoint.",
    "patient_bp_both_within_10kb": "Both breakends are within 10 kb of a patient "
                                   "breakpoint, not just one. One end landing "
                                   "near a patient breakpoint is far weaker "
                                   "evidence than both ends doing so.",
    "patient_svtype": "SV type of the nearest patient event, so a shared locus "
                      "can be told apart from a shared mechanism.",
    "patient_gene": "Gene annotated to the nearest patient breakpoint.",
    "patient_same_gene": "Patient gene matching this junction's gene, if any.",
    "patient_same_gene_and_svtype": "The patient event is in the same gene AND of "
                                    "the same SV type — the level-2 cross.",

    # -- panel of normals ------------------------------------------------
    "pon_count": "How many samples of the caller's panel of normals carry this "
                 "breakpoint. **An absolute COUNT, not a frequency** — it means "
                 "nothing without `panel_size_estimate`. Empty means absent from "
                 "the panel, which passes (`PON_ABSENT_MEANS`).",
    "pon_fraction": "`pon_count` as a fraction of the panel — the interpretable "
                    "form. The same `PON_MAX = 10` means 'seen in well under "
                    "1% of normals' against a panel of many thousands, and 'seen "
                    "in a fifth of them' against a panel of fifty. Always read "
                    "the count against `panel_size_estimate`.",
    "panel_size_estimate": "Largest `PON_COUNT` observed anywhere in this VCF, "
                           "used as the panel-size denominator. An estimate: the "
                           "true panel size is not recorded in the VCF.",
    "pass_pon": "Whether `pon_count < PON_MAX`. Computed and reported in every "
                "branch; whether it votes on `is_private` is `pon_in_privacy`.",
    "pon_in_privacy": "Whether the panel count participated in the `is_private` "
                      "verdict for this run. False means the panel was measured "
                      "and reported but did not remove anything, so a row can be "
                      "`is_private` while carrying a high `pon_count`.",

    # -- population frequency --------------------------------------------
    "gnomad_af": "Global allele frequency of the matching gnomAD-SV v4.1 record, "
                 "by reciprocal overlap at `GNOMAD_RECIPROCAL_OVERLAP`. Empty "
                 "means no match, or that no gnomAD resource was configured — "
                 "see `privacy_note`, since the two are not the same.",
    "gnomad_af_popmax": "The highest of the nine ancestry-group frequencies. The "
                        "conservative choice for a FILTER, since nothing passes "
                        "for being rare in one group — and the WRONG number to "
                        "describe a given sample, because that highest group may "
                        "be unrelated to the donor.",
    "gnomad_af_popmax_pop": "Which ancestry group `gnomad_af_popmax` came from. "
                            "Without it the maximum cannot be judged against a "
                            "donor's ancestry.",
    "gnomad_af_afr": "Frequency in African/African-American.",
    "gnomad_af_ami": "Frequency in the Amish group. A small, bottlenecked "
                     "population: frequencies here are noisy and often drive "
                     "`gnomad_af_popmax` without being informative elsewhere.",
    "gnomad_af_amr": "Frequency in Admixed American. A plausible reference for "
                     "the the cell line, whose donor ancestry is not documented.",
    "gnomad_af_asj": "Frequency in Ashkenazi Jewish.",
    "gnomad_af_eas": "Frequency in the East Asian group.",
    "gnomad_af_fin": "Frequency in the Finnish group. Bottlenecked, so it can "
                     "differ substantially from Non-Finnish European.",
    "gnomad_af_mid": "Frequency in Middle Eastern.",
    "gnomad_af_nfe": "Frequency in Non-Finnish European. A plausible reference "
                     "for the the cell line, whose donor ancestry is not documented.",
    "gnomad_af_sas": "Frequency in South Asian.",
    "gnomad_af_used": "The value the automatic verdict actually compared against "
                      "`GNOMAD_MAX_AF`. All nine groups are in the table so any "
                      "row can be re-judged against a different reference without "
                      "re-running anything.",
    "gnomad_af_field": "Which column `gnomad_af_used` was taken from, set by "
                       "`GNOMAD_AF_FIELD`.",
    "pass_gnomad": "Whether `gnomad_af_used < GNOMAD_MAX_AF` (0.001). A missing "
                   "frequency passes, so read it with `privacy_note`.",
    "is_private": "Clears BOTH recurrence filters. **Not 'unique to this cell "
                  "line'** — a match with the patient catalogue is the point of "
                  "the analysis. It asks whether the variant is common enough in "
                  "the population that patient and cell line would share it "
                  "anyway, which would make the match uninformative. Reports "
                  "label this column 'Not a common variant'.",
    "privacy_note": "Which basis produced `is_private`: gnomAD not evaluated "
                    "(panel only), panel reported but not applied, or both "
                    "applied. Distinguishes three states a bare boolean hides.",

    # -- expression (stage 7) --------------------------------------------
    "gene1_TPM": "Isofox gene-level TPM for `gene1`. Context: says the gene is "
                 "transcribed, not that the junction is.",
    "gene2_TPM": "Isofox gene-level TPM for `gene2`.",
    "isoform1_TPM": "Isofox TPM for `transcript_id1` specifically. A gene can be "
                    "well expressed through an isoform that does not contain the "
                    "breakpoint.",
    "isoform2_TPM": "Isofox TPM for `transcript_id2`.",
    "min_side_TPM": "The lower of the two sides' TPM. `min`, never `max`: an "
                    "event is in transcribed territory only if BOTH ends are. "
                    "Taking the maximum once reported a silent locus as being "
                    "covered by thousands of reads — the count belonged to its "
                    "transcribed partner, three orders of magnitude away.",
    "expressed": "`min_side_TPM >= TPM_EXPRESSED` (1.0). Context for the junction "
                 "test, not evidence for it.",

    # -- independent RNA evidence from the quantifier ---------------------
    "nearest_alt_sj_bp": "Distance in bp to the nearest novel splice junction "
                         "Isofox called. A breakpoint coinciding with one is "
                         "transcribed through a rearranged structure that a CIGAR "
                         "gap at a fixed offset can miss.",
    "alt_sj_frags": "Fragments supporting that novel splice junction, as counted "
                    "by Isofox — independent of our own read counting.",
    "alt_sj_type": "Isofox's classification of it, e.g. `NOVEL_INTRON`, "
                   "`NOVEL_EXON`.",
    "alt_sj_within_window": "The novel splice junction is within "
                            "`RNA_JUNCTION_MAX_DIST` (100 bp). The raw distance "
                            "travels alongside so a different cutoff can be used.",
    "isofox_fusion": "Isofox's own fusion call ID at this locus. A fusion "
                     "CALLER's verdict, which is stronger evidence for an "
                     "inter-chromosomal junction than hand-counted supplementary "
                     "alignments.",
    "isofox_fusion_support": "Fragments supporting that fusion call.",
    "retained_intron": "Retained introns Isofox reports at this locus, which "
                       "change what the transcript actually contains.",

    # -- junction evidence (stage 8) --------------------------------------
    "test": "Which RNA test the event's geometry allows: `sizegap` (CIGAR N gap), "
            "`insertion` (inserted bases in the read), `chimeric` (reads mapping "
            "to both loci), or `none`. Chosen from geometry, not from what would "
            "give an answer.",
    "test_reason": "Why that test, in words, including the size or insert length "
                   "that decided it.",
    "coverage_bp1": "Fragments overlapping breakend 1 in the RNA, MAPQ >= "
                    "`MIN_READ_MAPQ`. **Context, not evidence** — a breakpoint in "
                    "a highly expressed gene has thousands of reads whether or "
                    "not the junction exists.",
    "coverage_bp2": "Fragments overlapping breakend 2 in the RNA.",
    "min_coverage": "The lower of the two. `COVERAGE_SUMMARY = min`, never max.",
    "alignments_bp1": "Raw alignment records seen at breakend 1, before "
                      "collapsing to fragments. One fragment can emit several "
                      "records; support that read as 7 alignments was 3 "
                      "fragments (`COUNT_UNIQUE_FRAGMENTS`).",
    "alignments_bp2": "Raw alignment records seen at breakend 2.",
    "low_mapq_bp1": "Records at breakend 1 discarded for MAPQ < `MIN_READ_MAPQ` "
                    "(20). STAR emits only 0, 1, 3 and 255, so 3 marks a "
                    "multi-mapping read; a stack of them once looked like "
                    "support.",
    "low_mapq_bp2": "Records at breakend 2 discarded for low MAPQ.",
    "softclip_bp1": "Soft-clipped reads at breakend 1. **Does not tier an event** "
                    "(`SOFTCLIPS_TIER_EVENTS` is False): a clip shows a read ENDS "
                    "at the breakpoint, not that it crosses. Five such 'split "
                    "reads' once had supplementary alignments 3-156 Mb from the "
                    "partner breakend.",
    "softclip_bp2": "Soft-clipped reads at breakend 2.",
    "junction_reads": "**The evidence.** Unique fragments actually crossing the "
                      "junction, by whichever test applies. This is the only "
                      "quantity that tiers an event.",
    "junction_by_ngap": "Of those, how many were found as a CIGAR `N` gap "
                        "matching the deletion, within "
                        "`NGAP_POSITION_TOLERANCE`/`NGAP_SIZE_TOLERANCE`.",
    "junction_by_sa": "How many were found via an SA tag landing on the partner "
                      "breakend, within `SA_PARTNER_TOLERANCE` (1 kb). The right "
                      "chromosome alone is not enough.",
    "junction_by_insert": "How many were found by matching the inserted sequence "
                          "inside the read.",
    "rna_tier": "STRONG >= `STRONG_MIN_JUNCTION_READS` (5), SUGGESTIVE >= 3, WEAK "
                ">= 1, NONE = tested and zero, UNTESTABLE = geometry allows no "
                "valid test, NA = no RNA data. UNTESTABLE and NA are NOT "
                "negatives: recording a 2 bp deletion as NONE fabricates a "
                "result, since aligners do not emit `N` skips below their "
                "minimum intron.",
    # master_peptides.tsv only — one row per matched peptide
    # ================================================================

    "peptide": "The matched peptide sequence — identical in this sample and in "
               "the catalogue. Exact string match: a near-match is a different "
               "peptide and would be presented differently, or not at all.",

    # -- carried unchanged from the reference catalogue -------------------
    # These are the catalogue's own analysis columns. They are joined through
    # untouched so a match can be read against what the catalogue already knew
    # about that peptide; the pipeline neither recomputes nor validates them,
    # and their precise definitions belong to the catalogue, not to this repo.
    "ref_peptide": "The catalogue's peptide sequence — identical to `peptide` by "
                   "construction, kept so the join is auditable.",
    "ref_gene": "Gene the catalogue attributes the peptide to. What "
                "`gene_concordant` compares against.",
    "ref_svtype": "SV type the catalogue recorded for this peptide (DEL/DUP/TRA). "
                  "Used by the level-2 cross.",
    "or_original": "Catalogue column, carried unchanged: the odds ratio from the "
                   "cohort's immune-selection analysis as originally computed. "
                   "Not recomputed here.",
    "or_clean": "Catalogue column, carried unchanged: the same odds ratio after "
                "the catalogue's own cleaning step. Not recomputed here.",
    "significant_p05": "Catalogue column: whether that odds ratio was significant "
                       "at p < 0.05 in the cohort analysis.",
    "any_top50": "Catalogue column: whether the peptide reached the cohort "
                 "analysis's top-50 ranking.",
    "confirmed_clean": "Catalogue column: patients in which the peptide was "
                       "confirmed, after cleaning.",
    "confirmed_ge1": "Catalogue column: confirmed in at least 1 patient.",
    "confirmed_ge5": "Catalogue column: confirmed in at least 5 patients.",
    "present_not_conf_clean": "Catalogue column: patients where the peptide was "
                              "present but not confirmed.",
    "compat_no_pep_clean": "Catalogue column: patients HLA-compatible but "
                           "carrying no peptide — a denominator of the cohort's "
                           "2x2 depletion table.",
    "incompat_no_pep_clean": "Catalogue column: patients HLA-incompatible and "
                             "carrying no peptide — the other denominator.",
    "unique_breaks": "Catalogue column: distinct breakpoints producing this "
                     "peptide across the cohort. A peptide reachable by many "
                     "breakpoints is easier to hit by chance.",
    "is_cfs": "Catalogue column: the locus lies in a common fragile site, where "
              "breakpoints recur for reasons unrelated to selection.",
    "hla_pres_cov_European": "Catalogue column: fraction of the European HLA "
                             "population predicted to present this peptide.",
    "hla_pres_cov_NorthAm": "Catalogue column: the same for North American HLA "
                            "frequencies.",
    "hla_pres_cov_EastAsian": "Catalogue column: the same for East Asian.",
    "hla_pres_cov_Australian": "Catalogue column: the same for Australian.",
    "hla_pres_cov_mean4": "Catalogue column: mean of the four coverages above. "
                          "Population-level presentability of the peptide — NOT "
                          "a statement about this cell line, whose HLA type has "
                          "not been determined.",

    # -- the SV annotation joined onto each peptide -----------------------
    # Suffixed `_anno` because the peptide table already carries the generator's
    # own copy of the same fields; both are kept so a mismatch between the two is
    # visible rather than silently resolved.
    "prefix_anno": "Sample prefix from the annotation table (joined on `sv_id`).",
    "chrom1_anno": "Breakend 1 chromosome, from the annotation table.",
    "pos1_anno": "Breakend 1 position, from the annotation table.",
    "gene1_anno": "Breakend 1 gene, from the annotation table.",
    "transcript_id1_anno": "Breakend 1 transcript, from the annotation table.",
    "strand1_anno": "Breakend 1 transcript strand, from the annotation table.",
    "chrom2_anno": "Breakend 2 chromosome, from the annotation table.",
    "pos2_anno": "Breakend 2 position, from the annotation table.",
    "gene2_anno": "Breakend 2 gene, from the annotation table.",
    "transcript_id2_anno": "Breakend 2 transcript, from the annotation table.",
    "strand2_anno": "Breakend 2 transcript strand, from the annotation table.",
    "svtype_anno": "SV type, from the annotation table.",

    # -- the cross, per peptide -------------------------------------------
    "concordant_gene1": "The catalogue's gene matches breakend 1's gene.",
    "concordant_gene2": "The catalogue's gene matches breakend 2's gene. Tested "
                        "separately because a junction joins two genes and the "
                        "catalogue may annotate either side — testing only "
                        "breakend 1 lost 8 of 64 matches, all translocations.",
    "svtype_concordant": "This junction's SV type matches the catalogue's for "
                         "that peptide. Level-2 evidence: a shared mechanism.",

    # -- event-level context on a peptide row ------------------------------
    "n_peptides_in_event": "How many matched peptides come from the SAME "
                           "junction as this one. A sliding window over one "
                           "fusion yields up to ~40 overlapping sequences, so "
                           "this is the factor by which counting peptide rows "
                           "would overstate the number of findings — 21 in one "
                           "observed case.",
    "n_peptides": "The same count, carried from the event table.",
    "peptides": "All matched peptides of this event, `;`-separated. Makes the "
                "collapse from peptides to one event inspectable.",
    "ref_genes": "Catalogue gene(s) for this event's matched peptides, "
                 "`;`-separated.",

    # -- catalogue columns, carried through the cross ----------------------
    # From neoantigen_ranking_patch002.tsv, the aggregated patient catalogue.
    # Two of them mislead if read at face value; see the entries.
    "gene2_x": "The 3' partner gene from THIS sample's annotation. A pandas "
               "merge suffix, and an exact duplicate of `gene2` — verified "
               "identical in every row of a full run. Kept only because removing "
               "it would change existing column positions; use `gene2`.",
    "gene2_y": "The 3' partner gene recorded in the CATALOGUE for this peptide, "
               "as opposed to `gene2_x`/`gene2` which come from this sample. "
               "Distinct from `ref_gene`, which is `;`-separated and can list "
               "several; the two agree on the 52% of rows where the catalogue "
               "entry names a single gene.",
    "svpattern": "A NeoSV column preserved for schema compatibility and **always "
                 "empty here** — `build_patient_catalogue.py` does not carry it "
                 "through the regeneration path and emits it blank so the column "
                 "exists. Empty means not populated, not missing data.",
    "n_rows": "Rows in the catalogue contributing this peptide, i.e. how many "
              "(patient, SV) combinations produced it. Always >= `n_patients`.",
    "n_patients": "Distinct catalogue patients carrying this peptide. The "
                  "recurrence count on the patient side.",
    "n_samples": "**A copy of `n_patients`, not an independent count.** "
                 "`build_patient_catalogue.py` sets `n_samples = n_patients` "
                 "outright, because the regeneration path works one VCF per "
                 "patient and has no second sample grain to count. Reading it as "
                 "corroborating evidence double-counts the same number.",

    # -- MHC binding -------------------------------------------------------
    # Written by tools/propagate_binding_to_events.py. Present only for samples
    # where netMHCpan was actually run; blank elsewhere means NOT MEASURED.
    # Binder throughout = IC50 <= 500 nM AND %Rank_BA <= 2 AND %Rank_EL <= 2,
    # all three at once, the earlier analysis's own conjunction.
    "n_alleles_panel": "How many of the 218 cohort alleles this peptide binds. "
                       "A breadth statistic over the panel, NOT population "
                       "coverage: it is unweighted by allele frequency.",
    "binds_panel": "Binds at least one of the 218 alleles observed anywhere in "
                   "the reference cohort. **An upper bound, never a result.** It asks "
                   "whether *someone* could present the peptide, which is not "
                   "the question the earlier analysis answered; quoting it "
                   "against that analysis inflates the count roughly four-fold.",
    "panel_margin": "How far inside the binder definition the best panel allele "
                    "sits, as a fraction of the tightest cut. See "
                    "`autologous_margin` for the full definition.",
    "panel_limiting_cut": "Which of the three cuts is closest to failing for "
                          "`panel_margin` — the one actually holding the call up.",
    "panel_best_allele": "The panel allele giving `panel_margin`.",
    "n_patients_typed": "Catalogue patients carrying this peptide who have a "
                        "linkable class I genotype. Zero means no autologous "
                        "verdict is possible, NOT that the peptide fails.",
    "n_patients_presenting": "Of those, how many have an allele the peptide "
                             "binds — the autologous evidence count.",
    "n_patients_unevaluable": "Carrying patients whose genotype includes an "
                              "allele outside the predicted panel, so a "
                              "non-binding result there is unproven.",
    "binds_autologous": "Binds an allele of at least one patient who actually "
                        "carries this peptide. **This is the comparable "
                        "quantity** — it reproduces the rule used by the "
                        "earlier analysis (`annotate_hla.py` upstream), which "
                        "filtered each peptide against its own patient's ~6 "
                        "alleles rather than against a pooled panel.",
    "autologous_verdict": "`binder`, `non_binder`, or one of two `unevaluable_*` "
                          "states. **`unevaluable` is not `non_binder`**: 14% of "
                          "matched peptides come only from patients with no "
                          "linkable HLA, so the binder count is a floor.",
    "presenting_alleles": "The alleles that actually presented it in a carrying "
                          "patient, `;`-separated.",
    "autologous_margin": "How far inside the binder definition this call sits: "
                         "for each of the three cuts, the distance to it as a "
                         "fraction of the cut, then the SMALLEST of the three, "
                         "then the LARGEST across presenting alleles (a peptide "
                         "is as robust as its most comfortable allele). 0 means "
                         "exactly on a boundary. Scales are mixed on purpose "
                         "(nM against percentile), so it is a heuristic, not a "
                         "calibrated probability. Its use: 95% of the calls that "
                         "disagreed between netMHCpan 4.1 and 4.2 sat within 50% "
                         "of a threshold, so this column separates calls that "
                         "survive a re-run from calls that do not.",
    "autologous_limiting_cut": "Which cut is closest to failing — `affinity_nM`, "
                               "`rank_BA` or `rank_EL`. Among marginal calls it "
                               "splits roughly 43 / 12 / 44, so no single "
                               "threshold explains the fragility.",
    "autologous_best_allele": "The presenting allele giving `autologous_margin`.",
    "autologous_rank_EL": "%Rank_EL at `autologous_best_allele`.",
    "autologous_rank_BA": "%Rank_BA at `autologous_best_allele`.",
    "autologous_affinity_nM": "Predicted IC50 at `autologous_best_allele`.",
    "autologous_rank_Neo": "%Rank from netMHCpan 4.2's CEDAR neoepitope-finetuned "
                           "head at `autologous_best_allele`. Emitted only by a "
                           "`neo_*` criterion, since the default three-way rule "
                           "does not read that column.",
    "autologous_robustness": "`autologous_margin` binned: `flippable` (<=10% of "
                             "a cut), `marginal` (<=25%), `solid` (<=50%), "
                             "`robust` (>50%). Not a stricter biological claim — "
                             "the same criterion, minus the calls a perturbation "
                             "would flip.",
    "n_peptides_binding_panel": "Peptides of this event binding >= 1 panel "
                                "allele.",
    "n_peptides_binding_autologous": "Peptides of this event binding an allele "
                                     "of a patient who carries them.",
    "n_peptides_binding_autologous_robust": "Of those, the ones with "
                                            "`autologous_robustness` = `robust`.",
    "n_peptides_unevaluable": "Peptides of this event with no autologous "
                              "verdict available.",
    "binds_autologous_robust": "At least one robust autologous binder in this "
                               "event.",

    # -- junction evidence, broken down by mechanism -----------------------
    # Counted by this pipeline directly from the RNA BAM with pysam, NOT taken
    # from Isofox. Isofox's own calls (`alt_sj_*`, `isofox_fusion*`) are separate
    # columns and are independent of these — which is why they can disagree, and
    # why that disagreement is informative rather than a bug.
    "junction_by_ngap": "Fragments supporting the junction via a CIGAR `N` gap "
                        "matching the event size and position. `N` (skip), not "
                        "`D` (deletion): searching for `D` found zero reads "
                        "across a 221 bp deletion where `N` found 31. **Reads "
                        "counted only here, with none by `sa`, is the shape "
                        "splicing makes** — check `nearest_alt_sj_bp`.",
    "junction_by_sa": "Fragments supporting it via a supplementary alignment "
                      "landing near the partner breakend. The mechanism for "
                      "chimeric junctions, where no single read carries a gap.",
    "junction_by_insert": "Fragments supporting it via an inserted sequence.",
    "coverage_bp1": "Alignments in the window around breakend 1. **Context, not "
                    "evidence**: a breakpoint inside a highly expressed gene has "
                    "thousands of reads whether or not the junction exists.",
    "coverage_bp2": "Alignments in the window around breakend 2. Context, not evidence, for the same reason as `coverage_bp1`.",
    "min_coverage": "The quieter of the two coverages. `min`, never `max`.",
    "softclip_bp1": "Reads whose soft clip sits at breakend 1. **Diagnostic "
                    "only — soft clips never tier an event.** A clip shows a "
                    "read ENDS at the breakpoint; it says nothing about where "
                    "the rest of it went.",
    "softclip_bp2": "Reads whose soft clip sits at breakend 2. Diagnostic only, like `softclip_bp1`: a clip shows where a read ends, not where it goes.",
    "alignments_bp1": "Alignment records examined at breakend 1, before the "
                      "mapping-quality cut.",
    "alignments_bp2": "Alignment records examined at breakend 2, before the mapping-quality cut.",
    "low_mapq_bp1": "Alignments discarded at breakend 1 for ambiguous placement. "
                    "A high share means the region is repetitive.",
    "low_mapq_bp2": "Alignments discarded at breakend 2 for ambiguous placement. A high share means that end sits in repetitive sequence.",

    # -- lesion size -------------------------------------------------------
    "span": "Raw coordinate difference between the breakends, |pos2 - pos1|. "
            "**Not the lesion size** when an insertion is involved: a 34 bp "
            "insertion once read as a 2 bp deletion because the span was taken "
            "for the lesion. Compare with `event_size`.",
    "event_size": "Type-aware lesion size: the deleted, duplicated or inserted "
                  "length as the SV type implies, rather than the coordinate "
                  "span. Undefined (empty) for BND/TRA/SGL, where size is not a "
                  "meaningful quantity.",
    "insert_len": "Inserted bases at the junction. Non-zero here is why `span` "
                  "and `event_size` can disagree.",

    # -- common fragile sites ----------------------------------------------
    # Annotated against TWO catalogues on purpose: published definitions differ
    # by two orders of magnitude in how much genome they call fragile, so a hit
    # under one and not the other is the useful signal.
    "cfs_narrow_bp1": "Fragile site containing breakend 1 under the CONSERVATIVE "
                      "catalogue (~1% of the genome). Empty = not inside one.",
    "cfs_narrow_bp2": "Fragile site containing breakend 2 under the conservative catalogue. Empty = not inside one. The two ends are annotated separately because they can fall on different sides of a boundary.",
    "cfs_narrow_tier_bp1": "`core` or `extended` region of that catalogue, at "
                           "breakend 1.",
    "cfs_narrow_tier_bp2": "`core` or `extended` region of the conservative catalogue, at breakend 2.",
    "cfs_narrow_status": "`none` / `bp1` / `bp2` / `both` — which ends fall "
                         "inside. Kept per end rather than collapsed to a "
                         "boolean: one end deep inside a fragile site and the "
                         "other far outside is a different claim from both ends "
                         "inside.",
    "cfs_broad_bp1": "Fragile site containing breakend 1 under the PERMISSIVE "
                     "catalogue. **That catalogue covers ~65% of the genome, so "
                     "a hit here is close to the null expectation and is not "
                     "evidence on its own.**",
    "cfs_broad_bp2": "Fragile site containing breakend 2 under the permissive catalogue, whose large genome footprint makes a hit close to the chance expectation.",
    "cfs_broad_status": "`none` / `bp1` / `bp2` / `both` under the permissive "
                        "catalogue.",
    "cfs_agreement": "`both`, `narrow_only`, `broad_only` or `neither`. "
                     "**`broad_only` is the weak class** — a hit the "
                     "conservative catalogue does not support, and one the "
                     "permissive catalogue's genome footprint makes near-certain "
                     "by chance. Test any hit rate against the junctions that "
                     "produced no candidate, never against a uniform genome.",

    # -- lineage and candidate universe ------------------------------------
    # Written by tools/build_candidate_universe.py. Kept here rather than in that
    # script so a column cannot mean one thing in the universe and another in a
    # table derived from it.
    "acquired_in": "The line whose OWN variant calls produce this peptide — where "
                   "it entered the lineage. Semicolon-separated if two lines "
                   "generate the same sequence independently. Use this to ask "
                   "'what did this knockout add'.",
    "present_in": "Every line carrying the peptide: the acquiring line plus all "
                  "its descendants, since a clone inherits its parent's genome. "
                  "**Use this to ask 'what can line X present'** — filtering on "
                  "`acquired_in` silently drops everything the line inherited.",
    "n_lines_present": "How many lines carry it; 1 means private to one branch "
                       "of the lineage.",
    "matched_reference": "Whether the peptide also appears in the reference "
                         "cohort. **In the candidate universe this is an "
                         "annotation, never a selection.** When False every "
                         "cohort-side column is empty because the question was "
                         "not asked, not because the answer was no.",

    # -- presentability by the sample's OWN genotype ------------------------
    # Distinct from the cohort-side `autologous_*` columns, which ask the same
    # thing of a catalogue patient's genotype. Sharing a name would invite the
    # conflation the panel/autologous distinction exists to prevent.
    "presentable": "Binds at least one allele of the cell line's OWN class I "
                   "genotype. For a cell line there is no panel/autologous "
                   "ambiguity: the alleles ARE its own, so this means the cell "
                   "could present the peptide. Predicted binding, never observed "
                   "presentation.",
    "n_alleles_binding": "How many of the line's own alleles it binds.",
    "binding_alleles": "Which of the line's own alleles bind it, `;`-separated.",
    "presentation_best_allele": "The line's allele giving the largest margin.",
    "presentation_margin": "How far inside the binder definition the call sits: "
                           "the distance to each of the three cuts as a fraction "
                           "of that cut, the smallest of the three, then the "
                           "largest across binding alleles. 0 is exactly on a "
                           "boundary. Scales are mixed on purpose, so it ranks "
                           "fragility rather than estimating a flip probability.",
    "presentation_limiting_cut": "Which of the three cuts is closest to failing.",
    "presentation_robustness": "`presentation_margin` binned: `flippable` (<=10% "
                               "of a cut), `marginal` (<=25%), `solid` (<=50%), "
                               "`robust` (>50%). Not a stricter biological claim "
                               "— the same criterion minus the calls a re-run or "
                               "a predictor change would flip.",

    # -- cohort-side verdicts, renamed apart -------------------------------
    "binds_autologous_patient": "The COHORT-side verdict: binds an allele of a "
                                "catalogue patient who carries this peptide. "
                                "Distinct from `presentable`, which is this cell "
                                "line's own genotype. Empty unless "
                                "`matched_reference`.",
    "autologous_verdict_patient": "Cohort-side verdict: `binder`, `non_binder` or "
                                  "one of the `unevaluable_*` states. "
                                  "Unevaluable is NOT non-binder.",
    "autologous_margin_patient": "Cohort-side margin, defined as "
                                 "`presentation_margin` but over a patient's "
                                 "alleles.",
    "autologous_limiting_cut_patient": "Cohort-side limiting cut.",
    "autologous_best_allele_patient": "Cohort-side presenting allele.",
    "autologous_robustness_patient": "Cohort-side robustness bin.",
}


#: Alternative binder criteria reuse every MHC column name with a suffix. Their
#: meanings are the base meaning plus which rule produced them, so they are
#: derived rather than restated — a new criterion needs one line here, not
#: eighteen entries.
CRITERION_SUFFIXES = {
    "_neo_strong": "under the **neo_strong** rule (`%Rank_Neo <= 0.5`, netMHCpan "
                   "4.2's CEDAR neoepitope head) instead of the three-way "
                   "EL/BA rule. No cut is established for that head; 0.5 is "
                   "netMHCpan's generic strong-binder convention applied to it, "
                   "chosen because it passes a comparable number of "
                   "peptide-allele rows to the three-way rule. **Not comparable "
                   "with the earlier analysis.**",
    "_neo_weak": "under the **neo_weak** rule (`%Rank_Neo <= 2`), which is "
                 "several times more permissive than the three-way rule. **Not "
                 "comparable with the earlier analysis.**",
}


def check(columns) -> list[str]:
    """Columns with no entry. Used to fail a build rather than ship a guess."""
    return [c for c in columns if resolve(c) is None]


def resolve(column: str) -> str | None:
    """The meaning of a column, expanding criterion suffixes."""
    if column in MEANINGS:
        return MEANINGS[column]
    for suffix, note in CRITERION_SUFFIXES.items():
        if column.endswith(suffix):
            base = MEANINGS.get(column[:-len(suffix)])
            if base:
                return f"{base} Computed {note}"
    return None
