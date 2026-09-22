# Derivation of the criteria

Every non-obvious threshold in `src/svneo/criteria.py` was set by a specific
observation, not chosen by convention. This document records which, so that a
reader can distinguish a considered decision from an arbitrary one, and so that a
threshold can be re-derived rather than merely copied when the inputs change.

The observations come from applying two independently developed implementations
to the same dataset — a set of near-diploid cell lines tested against a
patient-derived SV-neoantigen catalogue — and reconciling where they disagreed.

## The reconciliation that produced the method

Three candidates were scored as strong RNA-supported recurrences by one
implementation. Re-measured against the same alignment file with the rules now in
this repository:

| Candidate | Geometry | Coverage, breakend 1 | Coverage, breakend 2 | Reads crossing the junction | Verdict |
|---|---|---|---|---|---|
| A | inter-chromosomal | 3 | 8,913 | 0 | not supported |
| B | 2 bp deletion, 34 bp insertion | 468 | 467 | 0 | untestable by a gap test |
| C | 220 bp deletion | 2,372 | 2,367 | 31 | supported, but `PON_COUNT` 3,513 |

None survives. Candidate A's supporting reads were soft-clipped reads whose
supplementary alignments landed 3–156 Mb from the partner breakend — they end
near the breakpoint and demonstrably do not cross it. Candidate B is a 2 bp event
that no CIGAR gap test can evaluate, and whose real lesion is a 34 bp insertion
misrepresented by breakend geometry. Candidate C is genuine: 31 crossing reads,
and a homozygous-deletion depth profile in the matched DNA — and it is present in
3,513 unrelated normals, so it is a common polymorphism rather than a private
event.

A general lesson followed, and is enforced structurally: **a conclusion must not
outlive the code that produced it.** The discrepancy arose because an
implementation was corrected while conclusions computed by its earlier version
were carried forward. Every output directory therefore carries the criteria
manifest and generator version that produced it.

## Threshold origins

| Constant | Observation that set it |
|---|---|
| `SOFTCLIPS_TIER_EVENTS = False` | "split reads" supporting candidate A were clipped reads whose supplementary alignments pointed megabases away from the partner |
| `COVERAGE_TIERS_EVENTS = False` | promoting on coverage alone marks any breakpoint inside an expressed gene as validated |
| `COVERAGE_SUMMARY = "min"` | taking the maximum reported candidate A's silent locus (3 reads) as having 8,913, the partner's count |
| `MIN_TESTABLE_GAP_SIZE = 20` | 2–3 bp events were recorded as negative by a test that cannot apply; aligners emit no `N` operation below their minimum intron |
| insertion test precedes the gap test | candidate B: breakend span 3, deleted length 2, inserted 34 bases — the span is not the lesion |
| `SA_PARTNER_TOLERANCE = 1000` | matching only the partner chromosome accepted alignments 3–156 Mb from the breakend |
| CIGAR `N`, not `D` | searching for `D` found zero reads across a 221 bp deletion; `N` found 31 (e.g. `112M218N39M`) |
| `STRONG_MIN_JUNCTION_READS = 5` | aligned with the default read-support cutoff of an established fusion-neoantigen tool; the previous value of 3 is retained as a lower tier |
| `GNOMAD_MAX_AF = 0.001` | candidates with `PON_COUNT` 1 and 3 carried population allele frequencies of 0.287 and 0.263 |
| `EVENT_DEDUP_COLUMN` | deduplicating on a per-peptide odds ratio left sliding-window peptides separate, inflating an event count more than fourfold |
| `REPORT_BOTH_EVENT_GRAINS` | a caller's grouped identifier merged clustered junctions at one locus and double-counted a shared breakend: 46 junctions reported as 39 events |
| `TRUST_COPY_NUMBER_FIELDS = False` | a clonal sample pair that should fit near purity 1.0 was fitted at purity 0.37 / ploidy 3.3, making copy-number-derived fields unreliable |
| `EXPRESSION_USES_BOTH_BREAKENDS` | using only the first breakend's transcript ignores the partner gene of an inter-chromosomal event |
| `TPM_EXPRESSED` applied to the isoform | in every expressed event examined, the disrupted isoform carried under 50% of the gene's TPM, and often under 1% |
| `REQUIRE_SPANS_JUNCTION = False` | tried as a filter and withdrawn: spanning the breakpoint is not what makes a peptide SV-derived, because after a frameshift every downstream residue exists only because of the variant. Enabling it discarded 63 of 64 catalogue matches (13 events collapsed to 1) while the null model put those same matches at 49x chance. Retained as an annotation |
| `MIN_READ_MAPQ = 20` | matches GRIDSS2's unmapped threshold (raised from GRIDSS 1's 10), the lineage that produced the SV calls. In RNA it is not tunable: STAR emits only 0, 1, 3 and 255, with zero reads between 10 and 20 over 73,232 alignments, so any cutoff in 4..254 is the same decision. In DNA it affects 0.5% of alignments |
| `COUNT_UNIQUE_FRAGMENTS` | one fragment can emit several alignment records; an insertion whose support read as 7 alignments was 3 fragments, one contributing 4 records at MAPQ 3 |
| `NULL_PERMUTATIONS` | candidate universes differing by three orders of magnitude between samples make raw match counts incomparable |
| germline vs somatic admission | every derived sample's "germline" call set was the same parental genome, 86.2% of breakends at identical positions; analysing it repeatedly would inflate any shared-across-samples claim |

## Regression tests

`tests/test_svneo.py` pins each of the above. They are not coverage tests: each
guards a rule whose violation yields a plausible but wrong result. A failure
there means a conclusion is about to become false.

## Retained design principles

Carried in unchanged because they were right: thresholds centralised in one
importable module with per-run manifests; parents analysed before descendants,
with attribution to the earliest sample carrying a breakpoint; a shared parental
call set analysed once; events rather than peptide rows as the reported unit;
`PON_COUNT` treated as an absolute count and reported as a fraction of panel
size; MHC prediction excluded from a sequence-identity test; and a missing
modality recorded as `NA` rather than as a negative result.
