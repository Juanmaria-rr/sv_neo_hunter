"""
criteria.py — every threshold and decision rule, in one importable place.

WHY ONE MODULE
--------------
The scientific claim of a recurrence analysis is that *identical scrutiny* was
applied to every sample. That is only credible if thresholds cannot drift between
samples, so every cutoff lives here, is imported by the stages, and is echoed
into each run's JSON manifest. Changing a threshold means changing it once.

Every constant carries (a) its rationale and (b) the observation that set it.
Where a value differs from the workflow it came from, the change is marked
`CHANGED:` with the reason. See ../../docs/PROVENANCE.md for the full history.

Thresholds may be overridden per run from the config file (`criteria:` block);
`apply_overrides()` is the only sanctioned way to do it, and it records what was
overridden so the manifest never silently disagrees with the code.
"""
from __future__ import annotations

import math
from collections import Counter

# ============================================================================
# STAGE 1 — VCF admission
# ============================================================================

#: Somatic VCFs: keep only these FILTER values.
#: `PON`      = found in the panel of normals (recurrent artefact / common germline)
#: `INFERRED` = breakend deduced by the caller from a copy-number transition, with
#:              no read evidence and often no mate; it cannot yield a trustworthy
#:              junction peptide because peptide generation needs two precise
#:              coordinates.
SOMATIC_KEEP_FILTERS = frozenset({"PASS"})

#: Somatic FILTER values admitted by a PANEL-REPORTING branch, where the panel
#: count is annotated but never removes a candidate.
#:
#: WHY THIS SET IS NEEDED AT ALL. The panel filter reaches a candidate by two
#: different routes, and relaxing one does nothing to the other:
#:
#:   germline VCFs  the panel is an INFO field, `PON_COUNT`, and it is this
#:                  pipeline that thresholds it (`PON_MAX`). Setting the
#:                  threshold to None is enough to stop filtering.
#:   somatic VCFs   the caller has ALREADY applied its panel filter and recorded
#:                  the verdict as `FILTER=PON`. `PON_MAX` never sees those
#:                  records, because admission drops them one step earlier for
#:                  not being PASS.
#:
#: A branch that only sets `pon_max=None` is therefore a no-op on somatic call
#: sets — it produces output identical to the filtered branch while being
#: labelled as unfiltered, which is worse than not running it. Measured on the
#: the sample set: 82, 58 and 96 records carry `FILTER=PON`, none of which
#: `pon_max` could ever have reached.
#:
#: `INFERRED` stays out. That is not a panel judgement: the caller deduced the
#: breakend from a copy-number transition with no read support, so there is no
#: junction sequence to translate. Admitting it would add records that can only
#: leave again at the peptide-generation step, while making the branch
#: incomparable with the filtered one in a second, unrelated respect.
SOMATIC_KEEP_FILTERS_REPORT_PON = frozenset({"PASS", "PON"})

#: Germline/parental PON threshold. A record is kept when PON_COUNT < this.
#:
#: THIS IS AN ABSOLUTE COUNT, NOT A FREQUENCY. Interpret it against the panel
#: size: with a ~12,000-sample panel `< 10` means "seen in < 0.08% of normals"
#: (strict); the same number against a 50-sample panel would mean "< 20%"
#: (permissive). Always report `pon_fraction()` alongside the raw threshold.
PON_MAX = 10

#: A record with no PON_COUNT is not in the panel, so it is kept. Treating
#: "absent" as a high count would discard the cleanest breakends.
PON_ABSENT_MEANS = 0

#: CHANGED (new): population-frequency filter, applied *in addition* to the PON.
#: A panel of normals has false negatives for inherited variation — observed:
#: CMSS1 PON_COUNT=1 but gnomAD-SV AF 0.287; AGMO PON_COUNT=3, AF 0.263. Both
#: pass PON<10 while being documented common polymorphisms.
GNOMAD_MAX_AF = 0.001

#: Reciprocal-overlap fraction required to call a cell-line SV the same event as
#: a gnomAD-SV record. 0.5 is the community default for SV matching.
GNOMAD_RECIPROCAL_OVERLAP = 0.5

#: Evaluate stages 6-8 (privacy, expression, junction reads) on EVERY admitted
#: junction, not only on those that produced a credible peptide match.
#:
#: The verdict-bearing stages only need the credible ones, but a table meant for
#: judging the call set is useless when 99.3% of its rows are empty: measuring
#: only what already passed cannot show what the rest looks like. The cost is one
#: BAM lookup per breakend, which is minutes, and nothing downstream changes —
#: `credible_events.tsv` is still built from the credible subset.
EVALUATE_ALL_JUNCTIONS = True

#: Ancestry groups reported by gnomAD-SV v4.1. All of them are extracted and
#: carried into the tables: which one is the right reference depends on the
#: sample's own ancestry, which is knowledge the pipeline does not have and must
#: not assume. A global average is dominated by the largest group (nfe) and can
#: hide a variant that is rare worldwide but common where the donor came from —
#: observed here: LPP at 0.028 in afr and 0.456 in eas, sixteen-fold apart.
GNOMAD_POPULATIONS = ("afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "sas")

#: Which frequency the automatic `is_private` verdict uses.
#:   "popmax" the highest frequency across GNOMAD_POPULATIONS — conservative,
#:            and the convention in clinical genetics: a variant common anywhere
#:            is not private
#:   "global" the cohort-wide AF
#:   "afr" | "nfe" | ...  one named group, when the donor's ancestry is known
#: Whatever is chosen, every population's frequency is still emitted, so the
#: verdict can be recomputed for a different reference without re-running.
GNOMAD_AF_FIELD = "popmax"


# ============================================================================
# STAGE 2 — peptide generation
# ============================================================================

#: Sliding-window peptide lengths spanning the junction. 8-11 covers the length
#: distribution of MHC-I ligands.
PEPTIDE_LENGTHS = (8, 9, 10, 11)

#: MHC binding prediction is NOT part of the recurrence test. Whether the sample
#: line's HLA would present a peptide has no bearing on whether the peptide
#: exists in it; an IC50/rank filter here would discard true sequence matches.
#: Presentability is a separate, downstream question (see `hla` in the config).
RUN_MHC_PREDICTION = False

#: Peptide column in NeoSV-Trace's all_neopeptides output. Getting this wrong
#: silently collapses the unique count (4 instead of 20,710 in one run).
PEPTIDE_COLUMN = "neopeptide"

#: Whether to keep ONLY peptides that span the junction. Default False.
#:
#: An SV-derived peptide need NOT cross the breakpoint. When the variant causes a
#: frameshift, every residue translated downstream exists only because of the SV,
#: and is as neoantigenic as one straddling the junction. The generator already
#: enforces novelty by subtracting the wild-type peptides of both transcripts
#: (`set(mut) - set(wt)`), so whatever survives is absent from the normal
#: proteome of those transcripts by construction — spanning is not what makes a
#: peptide neo.
#:
#: MEASURED: enabling this filter on a parental germline call set discarded
#: 63 of 64 catalogue matches (the surviving set collapsed from 13 events to 1),
#: while the null model put those same 64 matches at 47.7x chance (p < 0.001).
#: An earlier sample suggested only 6.3% of candidates did not span, which
#: generalised badly: that sample had few frameshift fusions.
#:
#: The flag is still emitted per peptide as an annotation, and is worth reporting
#: — a junction-spanning peptide is mechanistically tighter evidence — but it
#: must not gate the recurrence test.
REQUIRE_SPANS_JUNCTION = False
SPANS_JUNCTION_COLUMN = "spans_junction"


# ============================================================================
# STAGE 4 — sequence QC
# ============================================================================

#: A peptide is low-complexity if ANY of these fires. Short repetitive peptides
#: ("FFFFFFFFF" is present in real curated lists) match between unrelated genomes
#: by chance and are classic prediction artefacts.
LC_MIN_SHANNON_ENTROPY = 2.0
LC_MAX_SINGLE_AA_FRACTION = 0.5
LC_MAX_HOMOPOLYMER_RUN = 4
LC_MIN_DISTINCT_RESIDUES = 3


def shannon_entropy(peptide: str) -> float:
    """Shannon entropy (bits) over a peptide's residue composition."""
    if not peptide:
        return 0.0
    n = len(peptide)
    return -sum((c / n) * math.log2(c / n) for c in Counter(peptide).values())


def longest_homopolymer(peptide: str) -> int:
    """Length of the longest run of identical residues."""
    if not peptide:
        return 0
    best = run = 1
    for i in range(1, len(peptide)):
        run = run + 1 if peptide[i] == peptide[i - 1] else 1
        best = max(best, run)
    return best


def is_low_complexity(peptide: str) -> bool:
    """Any criterion firing is sufficient."""
    if not peptide:
        return True
    if len(set(peptide)) <= LC_MIN_DISTINCT_RESIDUES:
        return True
    if max(Counter(peptide).values()) / len(peptide) >= LC_MAX_SINGLE_AA_FRACTION:
        return True
    if longest_homopolymer(peptide) >= LC_MAX_HOMOPOLYMER_RUN:
        return True
    return shannon_entropy(peptide) < LC_MIN_SHANNON_ENTROPY


def is_credible(low_complexity: bool, self_peptide: bool,
                gene_concordant: bool) -> bool:
    """credible = high-complexity AND non-self AND gene-concordant.

    Gene concordance is a credibility criterion, not just an annotation: a match
    in a *different* gene is more likely short-sequence convergence than shared
    biology. NOTE the cost, which must be stated whenever "credible" is quoted:
    this excludes by construction a genuine convergent peptide arising from a
    different gene. It trades sensitivity for specificity.
    """
    return (not low_complexity) and (not self_peptide) and gene_concordant


# ============================================================================
# STAGE 5 — event deduplication
# ============================================================================

#: The deduplication key: the sample's own SV identifier. One genomic event =
#: one SV call, no matter how many sliding-window peptides it produced (up to
#: ~40 from one junction). Counting peptide rows overstates findings — in one
#: observed case by 21-fold (21 peptide rows, 1 deletion).
#:
#: NEVER key on a per-peptide value such as an odds ratio: window peptides from
#: one SV carry different values and survive as separate rows (178 vs 39).
EVENT_DEDUP_COLUMN = "sample_sv_id"

#: CHANGED: report the junction-level count alongside the caller's grouped id.
#: A caller's grouped-id string merges clustered junctions (3 junctions at one
#: locus -> 1) and can double-count a breakend shared by two groups. Observed:
#: 46 honest VCF junctions reported as 39 events. Both are emitted; the junction
#: count is the honest genomic unit and the grouped count is kept for continuity.
REPORT_BOTH_EVENT_GRAINS = True


# ============================================================================
# STAGE 6 — SV call confidence
# ============================================================================

#: High-confidence SV: BOTH breakends must satisfy all of these.
HC_MIN_SEGMAPQ = 30   # mapping quality of the highest-contributing segment
HC_MIN_VF = 5         # variant fragments supporting the breakend
HC_MIN_QUAL = 20      # caller confidence
HC_MIN_SV_SIZE = 100  # bp, intra-chromosomal only: sub-100 bp intra-chromosomal
                      # calls are frequently alignment artefacts around indels
                      # and homopolymers

#: Copy-number fields are EXCLUDED from every credibility judgement by default.
#: They depend on the caller's purity/ploidy fit, which is unreliable for clonal
#: cell lines (an observed fit of purity 0.37 / ploidy 3.3 on a clonal pair that
#: should sit near 1.0). Raw evidence fields (VF, SF, REF, SEGMAPQ) are unaffected.
TRUST_COPY_NUMBER_FIELDS = False


def breakend_is_hc(segmapq, vf, qual, sv_size=None) -> bool:
    """HC test for ONE breakend. `sv_size=None` for inter-chromosomal, which
    skips the size floor (a translocation has no meaningful size)."""
    if segmapq is None or vf is None or qual is None:
        return False
    ok = segmapq >= HC_MIN_SEGMAPQ and vf >= HC_MIN_VF and qual >= HC_MIN_QUAL
    if sv_size is not None:
        ok = ok and sv_size >= HC_MIN_SV_SIZE
    return ok


def event_is_hc(breakends: list[dict]) -> bool:
    """An event is HC only if EVERY breakend is HC.

    HC is REPORTED AS A SUBSET, NOT APPLIED AS A FILTER: a non-HC event with
    orthogonal RNA evidence can be more believable than an HC event with none.
    """
    if not breakends:
        return False
    return all(breakend_is_hc(b.get("SEGMAPQ"), b.get("VF"), b.get("QUAL"),
                              b.get("size")) for b in breakends)


def pon_fraction(pon_count: int, panel_size: int) -> float:
    """Absolute PON_COUNT -> fraction of the panel. `panel_size` may be estimated
    as max(PON_COUNT) over the VCF, which is a lower bound on the true size."""
    if panel_size <= 0:
        raise ValueError("panel_size must be positive")
    return pon_count / panel_size


def _absent(value) -> bool:
    """True when a value means "not recorded".

    Both `None` and NaN mean absent, and they are not interchangeable in
    comparisons: NaN < 10 is False, so a record with no PON_COUNT — i.e. one that
    is NOT in the panel of normals, the cleanest kind — would be judged to fail
    the panel filter. Values read from a DataFrame arrive as NaN, never None.
    """
    if value is None:
        return True
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return number != number      # NaN != NaN


def passes_pon(pon_count, pon_max: int = None) -> bool:
    """Stage-1/6 PON rule. An absent PON_COUNT means the breakend is not in the
    panel, which passes."""
    limit = PON_MAX if pon_max is None else pon_max
    value = PON_ABSENT_MEANS if _absent(pon_count) else float(pon_count)
    return value < limit


def is_private(pon_count, gnomad_af, pon_max: int = None,
               gnomad_max: float = None) -> bool:
    """CHANGED (new): an event is private only if it clears BOTH recurrence
    filters. Either one alone has documented false negatives."""
    # An unmatched event is absent from gnomAD, which passes. Note the asymmetry
    # this creates and report it: "not matched" and "not evaluated" look the same
    # here, which is why annotate_privacy() carries a `privacy_note` column.
    af = 0.0 if _absent(gnomad_af) else float(gnomad_af)
    return passes_pon(pon_count, pon_max) and \
        af < (GNOMAD_MAX_AF if gnomad_max is None else gnomad_max)


# ============================================================================
# STAGE 7 — expression
# ============================================================================

#: Expression is reported as a GRADIENT with its background rate, never as a bare
#: cutoff: in one reference transcriptome 55.3% of ALL genes clear TPM > 0, so
#: "expressed" at that threshold carries almost no information.
TPM_GRADIENT = (0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0)

#: Threshold when a single boolean is unavoidable. Applied to the DISRUPTED
#: ISOFORM, not the gene: gene TPM sums all isoforms and over-credits events (in
#: 13 of 13 expressed events the disrupted isoform was < 50% of gene TPM).
TPM_EXPRESSED = 1.0

#: CHANGED: check BOTH breakends' transcripts. Using only the first side ignores
#: the partner gene of an inter-chromosomal event.
EXPRESSION_USES_BOTH_BREAKENDS = True

#: Max distance (bp) from a breakpoint to an annotated splice junction for that
#: junction to count as support.
RNA_JUNCTION_MAX_DIST = 100


# ============================================================================
# STAGE 8 — RNA junction evidence  (see rna.py for the full rationale)
# ============================================================================

#: CHANGED: 3 -> 5, to match pVACfuse's default read-support cutoff so counts are
#: comparable with the field.
STRONG_MIN_JUNCTION_READS = 5

#: The previous headline threshold, kept as a visible intermediate tier so the
#: change in counts is auditable rather than silent.
SUGGESTIVE_MIN_JUNCTION_READS = 3

#: Smallest intra-chromosomal event for which a CIGAR `N`-gap test is valid.
#: Aligners do not emit `N` operations shorter than their minimum intron (~20 bp
#: for STAR), so below this the test is inapplicable, not negative.
MIN_TESTABLE_GAP_SIZE = 20

NGAP_SIZE_TOLERANCE = 10      # bp, matching a read's gap to the event size
NGAP_POSITION_TOLERANCE = 15  # bp, how close the gap starts to the breakpoint
SA_PARTNER_TOLERANCE = 1000   # bp, how close an SA must land to the partner
COVERAGE_WINDOW = 200         # bp, window for local coverage
MIN_SOFTCLIP_LEN = 10         # bp, informative soft-clip (diagnostic only)
MIN_INSERT_LEN = 10           # bp, insertion long enough to search for

#: Minimum mapping quality for a read to count as junction support.
#: A read placed ambiguously is not evidence for a specific locus.
#:
#: WHY 20, AND WHY IT IS BARELY A CHOICE
#: -------------------------------------
#: Field convention: GRIDSS2 treats alignments below MAPQ 20 as unmapped, having
#: raised the bar from GRIDSS 1's < 10. That lineage matters here — GRIDSS/ESVEE
#: is what produced the SV calls this pipeline consumes, so 20 keeps read-level
#: validation on the same footing as the caller.
#:
#: In RNA the threshold is not tunable at all. STAR does not emit a continuous
#: scale: MAPQ = 255 for unique alignments and int(-10*log10(1 - 1/Nmap))
#: otherwise, which yields only 3 (2 loci), 1 (3-4 loci) and 0 (>=5 loci).
#: Measured over 73,232 RNA alignments: exactly four distinct values (0, 1, 3,
#: 255) and ZERO reads between 10 and 20. Any cutoff in 4..254 is the same
#: decision — keep unique alignments, drop multi-mappers.
#:
#: In DNA (BWA-MEM2) the scale is continuous and the threshold does bite, but
#: lightly: of 74,825 alignments, 393 (0.5%) sit between 10 and 19, i.e. reads
#: with a >10% chance of being misplaced. 93.5% sit at MAPQ 60.
#:
#: Lowering to 10 would therefore recover nothing in RNA and only doubtful reads
#: in DNA. Raising it above 60 would silently discard every BWA alignment.
MIN_READ_MAPQ = 20

#: Count distinct FRAGMENTS, not alignments. A single fragment can contribute
#: several alignment records (supplementary/chimeric), and counting records
#: inflates support: 7 alignments, 3 fragments in the case above.
COUNT_UNIQUE_FRAGMENTS = True

#: THE TWO RULES THAT DEFINE THIS STAGE.
#: Soft-clips never tier an event: a clip shows a read ENDS at a breakpoint, not
#: that it crosses. Coverage never tiers an event: a breakpoint inside an
#: expressed gene has thousands of reads whether or not the junction exists.
SOFTCLIPS_TIER_EVENTS = False
COVERAGE_TIERS_EVENTS = False

#: Coverage is summarised across breakends with `min`, never `max`. An event is
#: in transcribed territory only if BOTH ends are. Using `max` reported a silent
#: locus (3 reads) as having 8,913 because its partner was highly expressed.
COVERAGE_SUMMARY = "min"


# ============================================================================
# STAGE 9 — null model and synthesis
# ============================================================================

#: Permutation iterations for the chance-match null. Raw match counts are not
#: comparable between samples whose candidate-peptide universes differ by orders
#: of magnitude (63,036 vs 64 in one observed pair).
NULL_PERMUTATIONS = 1000

#: Permutation strategy: shuffle each candidate peptide's residues, preserving
#: length and amino-acid composition. This holds composition bias constant, so
#: the null answers "how many matches would this many peptides of this
#: composition produce by chance".
NULL_STRATEGY = "shuffle_residues"

#: Breakpoint-proximity windows (kb) for cross level 3. Reported as a GRADIENT
#: because the count grows smoothly with window size and any single cutoff is
#: arbitrary. Both "at least one breakend within" and "both within" are emitted:
#: they answer different questions and the first is far more permissive.
PROXIMITY_WINDOWS_KB = (1, 10, 50, 100)

#: Tolerance (bp) when asking whether the same breakpoint exists in another
#: sample's VCF, for lineage attribution.
ATTRIBUTION_MATCH_TOLERANCE = 20


# ============================================================================
# Overrides and manifest
# ============================================================================

_OVERRIDDEN: dict = {}


def apply_overrides(overrides: dict) -> dict:
    """Apply `criteria:` overrides from the config, recording what changed.

    The only sanctioned way to change a threshold at run time. Every override is
    stored and echoed into the manifest, so a result can never be traced to a
    threshold that is not written down next to it.
    """
    applied = {}
    for key, value in (overrides or {}).items():
        name = key.upper()
        if not hasattr(globals(), "__contains__") or name not in globals():
            raise KeyError(f"unknown criterion '{key}' (expected one of the "
                           f"constants defined in criteria.py)")
        applied[name] = {"from": globals()[name], "to": value}
        globals()[name] = value
    _OVERRIDDEN.update(applied)
    return applied


def manifest() -> dict:
    """Every threshold as a flat dict, written into each run's manifest."""
    keys = [k for k, v in globals().items()
            if k.isupper() and not k.startswith("_")
            and isinstance(v, (int, float, str, bool, tuple, frozenset))]
    out = {k: (sorted(globals()[k]) if isinstance(globals()[k], frozenset)
               else list(globals()[k]) if isinstance(globals()[k], tuple)
               else globals()[k]) for k in sorted(keys)}
    if _OVERRIDDEN:
        out["_overridden"] = _OVERRIDDEN
    return out
