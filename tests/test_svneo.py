"""
test_svneo.py — the rules that must not regress.

These are not coverage tests. Each one pins a decision that was got wrong at
least once in the two source workflows and produced a published wrong number.
If one of these fails, a conclusion somewhere is about to become false.

    pytest tests/ -v        (or: python tests/test_svneo.py)
"""
import os
import pathlib
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from svneo import confidence, criteria, generators, rna, null_model, vcf  # noqa: E402


# ---------------------------------------------------------------------------
# RULE 1 — only junction-crossing reads tier an event
# ---------------------------------------------------------------------------

def test_coverage_alone_never_reaches_a_supported_tier():
    """A breakpoint in a highly expressed gene has thousands of reads whether or
    not the junction exists. The workflow this replaces promoted such events to
    MODERATE on coverage alone, which is how a silent locus whose PARTNER carried
    thousands of reads became 'RNA-validated'."""
    assert rna.rna_tier(junction_reads=0, testable=True) == "NONE"


def test_tier_thresholds():
    assert rna.rna_tier(5, True) == "STRONG"       # aligned with pVACfuse
    assert rna.rna_tier(3, True) == "SUGGESTIVE"   # the old STRONG, kept visible
    assert rna.rna_tier(1, True) == "WEAK"
    assert rna.rna_tier(0, True) == "NONE"


def test_softclips_and_coverage_are_declared_non_tiering():
    """The manifest must state these, because a reader cannot infer them."""
    assert criteria.SOFTCLIPS_TIER_EVENTS is False
    assert criteria.COVERAGE_TIERS_EVENTS is False
    assert criteria.COVERAGE_SUMMARY == "min"


# ---------------------------------------------------------------------------
# RULE 3 — no valid test means UNTESTABLE, never NONE
# ---------------------------------------------------------------------------

def test_tiny_event_is_untestable_not_negative():
    """A 3 bp deletion cannot produce a CIGAR N gap: aligners do not emit N
    operations below their minimum intron. Recording it as NONE fabricates a
    negative result."""
    test, reason = rna.select_test("DEL", "3", "3", event_size=3, insert_len=0)
    assert test == "none"
    assert "minimum intron" in reason
    assert rna.rna_tier(0, testable=False) == "UNTESTABLE"


def test_insertion_driven_event_gets_the_insertion_test():
    """Observed case: span 3, deleted length 2, but 34 inserted bases. A gap test
    can never confirm it; it must be tested as an insertion."""
    test, reason = rna.select_test("DEL", "12", "12", event_size=2, insert_len=34)
    assert test == "insertion"
    assert "34" in reason


def test_interchromosomal_uses_the_chimeric_test():
    test, _ = rna.select_test("BND", "1", "8", event_size=None, insert_len=0)
    assert test == "chimeric"


def test_large_intrachromosomal_uses_the_gap_test():
    test, _ = rna.select_test("DEL", "15", "15", event_size=220, insert_len=0)
    assert test == "sizegap"


def test_support_is_counted_per_fragment_with_a_mapq_floor():
    """One fragment can emit several alignment records. An insertion whose RNA
    support read as 7 alignments was 3 fragments, one of them contributing 4
    records at MAPQ 3 (STAR's multi-mapping value)."""
    assert criteria.COUNT_UNIQUE_FRAGMENTS is True
    assert criteria.MIN_READ_MAPQ >= 10


def test_missing_event_size_does_not_crash():
    """NaN is truthy in Python: `if size:` lets a missing size through and
    int(float(nan)) then raises. Inter-chromosomal events have no size by
    definition, so every size arriving from a DataFrame can be NaN."""
    nan = float("nan")
    test, _ = rna.select_test("DEL", "3", "3", event_size=nan, insert_len=nan)
    assert test == "none"
    test, _ = rna.select_test("BND", "1", "8", event_size=nan, insert_len=nan)
    assert test == "chimeric"


def test_missing_rna_is_na_not_negative():
    assert rna.rna_tier(0, testable=True, has_rna=False) == "NA"


# ---------------------------------------------------------------------------
# SA tags: the right chromosome is not enough
# ---------------------------------------------------------------------------

def test_sa_tag_must_land_on_the_partner_not_just_its_chromosome():
    """The reads that produced a false STRONG carried SA tags to the correct
    chromosome but 3-156 Mb from the partner breakend."""
    assert rna.sa_hits_partner("chr8,127700000,+,50M100S,60,0;", "8", 127700100)
    assert not rna.sa_hits_partner("chr8,3000000,+,50M100S,60,0;", "8", 127700100)
    assert not rna.sa_hits_partner("chr12,125065102,+,50M,60,0;", "15", 68420999)


# ---------------------------------------------------------------------------
# Type-aware event size
# ---------------------------------------------------------------------------

def test_admitted_vcf_selects_by_id_not_coordinate():
    """Two breakends can share a position and belong to different junctions.
    Selecting by (chrom, pos) writes out a REJECTED record because its neighbour
    passed, letting a common germline junction reach the peptide generator."""
    import io, tempfile, pandas as pd
    from svneo import vcf as vcf_mod
    body = ("##fileformat=VCFv4.2\n"
            "chr17\t1000000\t26224\tN\t]chr17:1000037]CAG\t50\tPASS\tSVTYPE=BND\n"
            "chr17\t1000000\t26225\tN\t]chr17:999894]CAG\t50\tPASS\tSVTYPE=BND\n")
    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as fh:
        fh.write(body)
        source = fh.name
    admitted = pd.DataFrame([{"chrom": "17", "pos": 1000000, "id": "26225"}])
    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as fh:
        out = fh.name
    written = vcf_mod.write_admitted_vcf(source, admitted, out)
    assert written == 1, f"wrote {written} records, expected only the admitted one"
    assert "26224" not in open(out).read(), "a rejected breakend was written out"


def test_event_size_is_type_aware():
    """The breakend span is not the lesion. Using the span made 248 insertions
    look like 1 bp events."""
    assert vcf.event_size("DEL", 1000, 1221) == 220     # span - 1 bases deleted
    assert vcf.event_size("DUP", 1000, 1221) == 221     # the duplicated tract
    assert vcf.event_size("INS", 1000, 1001, insert_len=34) == 34
    assert vcf.event_size("BND", 1000, 999999) is None  # undefined


# ---------------------------------------------------------------------------
# Privacy: both filters, and "not evaluated" is not "clean"
# ---------------------------------------------------------------------------

def test_privacy_requires_both_filters():
    """A panel of normals has false negatives for inherited variation: an
    observed event with PON_COUNT=1 had a population AF of 0.287."""
    assert criteria.is_private(pon_count=1, gnomad_af=0.0001)
    assert not criteria.is_private(pon_count=1, gnomad_af=0.287)
    assert not criteria.is_private(pon_count=3513, gnomad_af=0.0)


def test_absent_pon_count_means_not_in_panel():
    """Absent means NOT in the panel, which passes — and absent arrives as NaN
    from a DataFrame, never as None. NaN < 10 is False, so without an explicit
    check the cleanest events (no PON_COUNT at all) were judged NOT private."""
    assert criteria.passes_pon(None)
    assert criteria.passes_pon(float("nan"))
    assert criteria.is_private(float("nan"), float("nan"))
    assert not criteria.is_private(3513, float("nan"))


# ---------------------------------------------------------------------------
# Credibility and the null model
# ---------------------------------------------------------------------------

def test_low_complexity_catches_the_classic_artefact():
    assert criteria.is_low_complexity("FFFFFFFFF")
    assert not criteria.is_low_complexity("SLYNTVATL")


def test_credible_requires_all_three():
    assert criteria.is_credible(False, False, True)
    assert not criteria.is_credible(True, False, True)    # low complexity
    assert not criteria.is_credible(False, True, True)    # self peptide
    assert not criteria.is_credible(False, False, False)  # different gene


def test_null_model_reports_a_bounded_p_value():
    """A real signal must not be reported as p = 0; it is p < 1/n."""
    reference = ["SLYNTVATL", "GILGFVFTL", "NLVPMVATV"]
    candidates = list(reference) + [f"AAAAAAAA{i}" for i in range(50)]
    result = null_model.permutation_test(candidates, reference, n_permutations=100)
    assert result["observed"] == 3
    assert result["p_value_is_bounded"] or result["p_value"] < 0.05
    assert "matches_per_1000_candidates" in result


def test_null_model_makes_sample_sizes_comparable():
    """The point of the rate: 3 matches from 50 candidates is not the same
    finding as 3 from 50,000."""
    reference = ["SLYNTVATL"]
    small = null_model.permutation_test(["SLYNTVATL"] + ["ACDEFGHIK"] * 1,
                                        reference, n_permutations=50)
    assert small["matches_per_1000_candidates"] > 100


def test_spans_junction_annotates_but_does_not_filter():
    """Spanning the junction is NOT what makes a peptide SV-derived: after a
    frameshift every downstream residue exists only because of the variant.
    Filtering on it discarded 63 of 64 real catalogue matches in a run where the
    null model put those same matches at 47.7x chance."""
    assert criteria.REQUIRE_SPANS_JUNCTION is False
    assert criteria.SPANS_JUNCTION_COLUMN == "spans_junction"


# ---------------------------------------------------------------------------
# Stage 2 is pluggable, and its contract is enforced
# ---------------------------------------------------------------------------

def test_generator_registry_and_unknown_name_fails_loudly():
    """A typo in `peptide_generator:` must not silently fall back to a default."""
    assert "neosv_trace" in generators.REGISTRY
    assert "precomputed" in generators.REGISTRY
    try:
        generators.get_generator("neosvtrace")     # missing underscore
        raise AssertionError("unknown backend was accepted")
    except generators.GeneratorError as error:
        assert "Available" in str(error)


def test_generator_contract_rejects_a_missing_required_column():
    """Reading the wrong peptide column once collapsed a count from 20,710 to 4.
    A backend that omits a required column must fail, not degrade."""
    import pandas as pd
    good_anno = pd.DataFrame([{"sv_id": "1", "chrom1": "1", "pos1": 100,
                               "gene1": "A", "chrom2": "1", "pos2": 200}])
    bad = pd.DataFrame([{"sv_id": "1", "gene1": "A"}])          # no `neopeptide`
    try:
        generators.PeptideGenerator.validate_output(bad, good_anno, "test")
        raise AssertionError("missing column was accepted")
    except generators.GeneratorError as error:
        assert "neopeptide" in str(error)


def test_vendored_generator_is_present_and_licensed():
    """The repo must run without an external checkout, and must ship the licence
    that permits redistributing the vendored tool."""
    from svneo.generators.neosv import VENDOR_DIR
    assert os.path.isdir(os.path.join(VENDOR_DIR, "neosv")), \
        f"vendored NeoSV missing from {VENDOR_DIR}"
    assert os.path.exists(os.path.join(VENDOR_DIR, "LICENSE.NeoSV")), \
        "MIT licence must ship with the vendored copy"
    assert os.path.exists(os.path.join(VENDOR_DIR, "VENDOR.md"))


def test_frame_patch_is_applied():
    """Patch 001 removes a stop-codon offset that is wrong for pyensembl > 2.3.13.
    Unpatched, the tool emits 170 peptides where the patched one emits 64, with
    none in common. Reverting it silently invalidates every identity match."""
    from svneo.generators.neosv import VENDOR_DIR
    source = pathlib.Path(VENDOR_DIR) / "neosv" / "fusion_class.py"
    text = source.read_text()
    assert "start = end - self.cut_length + 1 - 3" not in text, \
        "the -3 stop-codon offset is back: see vendor/patches/001-*.md"
    assert "start = end - self.cut_length + 1" in text
    assert (pathlib.Path(VENDOR_DIR) / "patches" /
            "001-pyensembl-stop-codon-frame.patch").exists()


# ---------------------------------------------------------------------------
# The panel filter acts at two points, and a branch must relax both
#
# These guard a specific silent failure: a branch declared with `pon_max: null`
# and nothing else produced output IDENTICAL to the filtered branch on somatic
# call sets, while being labelled unfiltered. The caller had already written its
# panel verdict into FILTER, one step before `pon_max` is consulted.
# ---------------------------------------------------------------------------

def _breakends(*filters) -> "pd.DataFrame":
    """A minimal admitted-breakend frame: mated, no panel count of our own."""
    import pandas as pd
    return pd.DataFrame({
        "filter": list(filters),
        "paired": [True] * len(filters),
        "pon_count": [None] * len(filters),
    })


def test_pon_max_alone_cannot_reach_a_caller_filtered_somatic_record():
    """The no-op that motivated `admit_panel_filtered`."""
    breakends = _breakends("PASS", "PON", "PON", "INFERRED")
    strict, _ = vcf.admit(breakends, "somatic", pon_max=None)
    assert len(strict) == 1, \
        "pon_max=None must NOT admit FILTER=PON on a somatic set; if it does, " \
        "the two branches are no longer distinguishable by admission alone"


def test_admit_panel_filtered_admits_pon_and_reports_how_many():
    breakends = _breakends("PASS", "PON", "PON", "INFERRED")
    relaxed, funnel = vcf.admit(breakends, "somatic", pon_max=None,
                               admit_panel_filtered=True)
    assert len(relaxed) == 3, f"expected PASS + 2 PON, got {len(relaxed)}"
    assert funnel["caller_pon_admitted"] == 2, \
        "the funnel must state how many panel-flagged records were admitted, " \
        "or the relaxation is invisible in the output"
    # INFERRED is not a panel judgement and must stay out in both branches.
    assert "INFERRED" not in set(relaxed["filter"])


def test_pon_still_votes_on_privacy_unless_told_otherwise():
    """Relaxing admission alone leaves `is_private` filtering by panel count."""
    import pandas as pd
    events = pd.DataFrame({"sample_sv_id": ["a"], "pon_count": [3513.0]})

    default = confidence.annotate_privacy(events, panel_size=12000, gnomad_af={})
    assert not bool(default.is_private.iloc[0]), \
        "a high PON_COUNT must fail privacy while the panel votes"

    reported = confidence.annotate_privacy(events, panel_size=12000, gnomad_af={},
                                           pon_in_privacy=False)
    assert bool(reported.is_private.iloc[0]), \
        "with the panel reporting only, privacy rests on gnomAD alone"
    # Reported, not discarded: the deciding value stays in the table either way.
    assert reported.pon_count.iloc[0] == 3513.0
    assert not bool(reported.pass_pon.iloc[0]), \
        "pass_pon must still record that the event IS in the panel"
    assert "not applied" in str(reported.privacy_note.iloc[0]), \
        "the note must say the panel was reported rather than applied, or a " \
        "reader cannot tell which basis produced is_private"


def test_minus_strand_cds_is_read_in_transcript_order():
    """Patch 002. Coding exon ranges arrive ordered by COORDINATE, which is
    reading order on the plus strand and its reverse on the minus strand. Every
    consumer assumes reading order — `get_cds_range` says so in its own
    docstring — so on a minus-strand transcript the introns came out with
    start > end (unmatchable) and exons were counted from the wrong end.

    Effect before the patch: 0 of 188 minus-strand regions returned the correct
    5' head, against 87 of 87 on the plus strand. An intronic breakpoint gave an
    empty head; an exonic one gave the wrong length, and near the end of a
    transcript gave TOO MUCH — adding sequence the gene does not contribute,
    which the sliding window turns into peptides.

    This asserts the source, not the behaviour, so it runs without an annotation
    cache. `tools/diagnose_minus_strand_cds.py` checks the behaviour.
    """
    from svneo.generators.neosv import VENDOR_DIR
    source = pathlib.Path(VENDOR_DIR) / "neosv" / "transcript_utils.py"
    text = source.read_text()

    assert "cds_ranges = sorted(transcript.coding_sequence_position_ranges)" in text, \
        "get_cds_range no longer sorts: see vendor/patches/002-*.md"
    assert "cds_ranges = cds_ranges[::-1]" in text, \
        "the minus-strand reversal is gone: see vendor/patches/002-*.md"
    # get_noncds_range must consume the ordered accessor, not the raw one, or
    # it rebuilds the inverted intervals the patch exists to remove. Count in
    # CODE only: the patch comments name the accessor they replaced, and
    # counting those would make this assertion depend on prose.
    code = "\n".join(line.split("#", 1)[0] for line in text.splitlines())
    assert code.count("transcript.coding_sequence_position_ranges") == 1, \
        "something reads the raw accessor again; route it through get_cds_range"
    assert (pathlib.Path(VENDOR_DIR) / "patches" /
            "002-minus-strand-cds-order.patch").exists()


# ---------------------------------------------------------------------------
# A freshly copied template must fail with something the user can act on
#
# Copying config/template.yaml and running --dry-run is the first thing anyone
# does with this repository. It used to raise FileNotFoundError on the first
# placeholder it met, behind a stack trace — which reads as a broken tool rather
# than an unconfigured one, and reveals one bad path per invocation.
# ---------------------------------------------------------------------------

def test_unedited_template_reports_every_placeholder_at_once():
    import tempfile
    from svneo import config as config_mod

    template = pathlib.Path(__file__).parent.parent / "config" / "template.yaml"
    assert template.exists(), "config/template.yaml is referenced by the README"

    with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as handle:
        handle.write(template.read_text())
        copy = handle.name
    try:
        config_mod.load(copy)
    except config_mod.ConfigError as error:
        message = str(error)
    except Exception as error:                                  # noqa: BLE001
        raise AssertionError(
            f"expected ConfigError, got {type(error).__name__}: {error}")
    else:
        raise AssertionError("an unedited template must not load")
    finally:
        os.unlink(copy)

    assert "reference.peptides" in message
    # Every bad path, not just the first: otherwise fixing the config is a
    # one-error-per-run guessing game.
    assert message.count("/path/to/") > 1, \
        f"only one placeholder reported:\n{message}"
    assert "template.yaml" in message, \
        "the message should say what to do, not only what is wrong"


# ---------------------------------------------------------------------------
# Every column of a shared table must carry an explanation
#
# The master tables are what someone is handed when they want to check a result
# themselves. A column whose meaning has to be guessed is worse than an absent
# one: `pon_count` looks like a frequency, `coverage_bp1` looks like evidence,
# `softclip_bp1` looks like junction support. None of those readings is right.
# ---------------------------------------------------------------------------

def test_every_documented_column_has_a_written_meaning():
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "tools"))
    import column_meanings

    empty = [name for name, text in column_meanings.MEANINGS.items()
             if not text or len(text) < 25]
    assert not empty, f"placeholder meanings: {empty}"

    # Any table already produced must be fully covered. Nothing to check on a
    # fresh clone, which has no results — the assertion is on what exists.
    import glob
    for path in glob.glob(os.path.join(
            os.path.dirname(__file__), "..", "results", "master_*.tsv")):
        if "column_dictionary" in path:
            continue
        with open(path) as handle:
            columns = handle.readline().rstrip("\n").split("\t")
        missing = column_meanings.check(columns)
        assert not missing, \
            f"{os.path.basename(path)} has undocumented columns: {missing}"


# ---------------------------------------------------------------------------
# The manifest must be complete
# ---------------------------------------------------------------------------

def test_manifest_carries_every_threshold_used_in_a_verdict():
    manifest = criteria.manifest()
    for key in ("STRONG_MIN_JUNCTION_READS", "MIN_TESTABLE_GAP_SIZE", "PON_MAX",
                "GNOMAD_MAX_AF", "HC_MIN_SEGMAPQ", "COVERAGE_SUMMARY",
                "NULL_PERMUTATIONS"):
        assert key in manifest, f"{key} missing from the manifest"


if __name__ == "__main__":
    failures = 0
    for name, function in sorted(globals().items()):
        if name.startswith("test_") and callable(function):
            try:
                function()
                print(f"  PASS  {name}")
            except AssertionError as error:
                failures += 1
                print(f"  FAIL  {name}: {error}")
    print(f"\n{failures} failure(s)")
    sys.exit(1 if failures else 0)
