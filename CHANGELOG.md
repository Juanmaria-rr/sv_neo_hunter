# Changelog

Notable changes to this project. Format follows [Keep a Changelog](https://keepachangelog.com/);
versioning is [semantic](https://semver.org/).

## [0.1.0] — unreleased

First working version: the analysis as configurable software rather than a
per-dataset script.

### Added

- **Panel-reporting branches.** A branch may now carry the panel of normals as
  an annotation without filtering on it, via `admit_panel_filtered:` (admit
  records the caller rejected as `FILTER=PON`) and `pon_in_privacy:` (keep the
  panel count out of the `is_private` verdict). Both are needed because the
  panel reaches a candidate twice, and on a somatic call set the caller has
  already applied it as a `FILTER` value — a branch setting only `pon_max: null`
  is a no-op there, producing output identical to the filtered branch under an
  unfiltered label.
- **`events_private_hc_and_rna`** in `summary.json`: the count of events
  satisfying every criterion. The pre-existing `events_hc_and_rna` does not
  include privacy, and both are now reported separately — an event can be
  confidently called and transcribed while being a common polymorphism.
- **HTML validation reports** (`tools/build_html_report.py`), one per sample and
  branch: each filter rendered with what it measures, how, why, its threshold and
  the threshold's justification, beside the count it removed and the names of
  what it removed. Grain changes along the cascade are marked explicitly.

### Fixed

- The HTML report described `events_hc_and_rna` as "satisfies every criterion",
  which overstated a common germline event as a surviving candidate. Both reports
  now name the criteria each figure covers.

- **Config-driven runs.** A run is defined by a reference peptide catalogue and a
  sample set; no cohort, tissue or cell line is named in the code. Samples form a
  lineage via `parent:`, which drives execution order and event attribution.
- **Nine stages** — admission and breakend pairing (`vcf`), peptide generation
  (`generators`), the three-level cross with sequence QC and event
  deduplication (`cross`), SV confidence and privacy (`confidence`), expression
  and junction evidence (`rna`), funnel and attribution (`synthesis`), plus a
  permutation null (`null_model`).
- **Pluggable peptide generation.** Stage 2 sits behind a documented contract
  (`generators/base.py`) with three backends: `neosv` (default), `neosv_trace`
  (external checkout), `precomputed` (reuse a stored output). No other module
  imports a generation tool.
- **Vendored NeoSV** under its MIT licence, so the repository runs end to end
  without an external checkout. See `vendor/VENDOR.md`.
- **`tools/compare_generators.py`** — runs two generator packages over the same
  input and reports whether their peptide sets are interchangeable. Swapping a
  generator without this check can silently invalidate every identity match.
- **Regression tests** pinning the rules whose violation produces plausible but
  wrong results, rather than pursuing coverage.
- **`tools/build_master_table.py`** — one auditable row per candidate peptide,
  carrying every catalogue column, every generator column, and one explicit
  pass/fail column per criterion, plus `first_failed_gate`.
- **`tools/inspect_insertion.py`** — prints the individual reads behind an
  insertion call with mapping quality, duplicate flag and in-read offset, against
  a background rate over random loci. A count cannot distinguish real support
  from a duplicate stack; this can.
- **`notebooks/pipeline_walkthrough_pyspark.ipynb`** — every stage re-derived
  independently in PySpark and asserted against the packaged pipeline, 16 checks.
  Generated from `notebooks/build_walkthrough.py`.

### Fixed

- **Stage 7 joined the annotation positionally** (first column, `prefix`) instead
  of on `sv_id`, so no expression was ever attached and every event returned with
  no gene, transcript or TPM.
- **NaN handling in two filters.** NaN is truthy in Python and every comparison
  against it is False, so an event with no size crashed the junction test, and an
  event with no `PON_COUNT` — absent from the panel, the cleanest possible
  record — was judged **not** private, inverting the filter for the best
  candidates.
- **The peptide generator received the raw VCF** rather than the stage-1 admitted
  set, building peptides from panel and copy-number-inferred breakends.
- **Reading-frame correction in the vendored generator**
  (`vendor/patches/001-pyensembl-stop-codon-frame`). NeoSV subtracts 3 from the
  3' fusion start to compensate for pyensembl behaviour that changed at v2.3.13
  (openvax/pyensembl#176). With later releases the subtraction shifts the fusion
  one amino acid out of frame: measured on 51 SVs with pyensembl 2.10.1, the
  unpatched tool emits 170 peptides against 64 patched, **with none in common**,
  because frame-shifted sequences survive the wild-type subtraction as spurious
  neopeptides.

### Method changes relative to the workflows this replaces

- RNA evidence is tiered on **junction-crossing reads only**. Coverage is
  context and soft-clips are diagnostic; neither can promote a tier.
- Coverage is summarised across breakends with **`min`, never `max`** — taking
  the maximum attributes a partner locus's expression to a silent one.
- An event with no applicable test is **`UNTESTABLE`**, not `NONE`. The test is
  selected from event geometry: chimeric, gap, or inserted-sequence.
- Privacy requires **two** filters, panel of normals **and** population
  frequency; either alone has documented false negatives.
- `spans_junction` is emitted per peptide as an **annotation, not a filter**
  (`REQUIRE_SPANS_JUNCTION = False`). Spanning the breakpoint is not what makes a
  peptide SV-derived: after a frameshift every downstream residue exists only
  because of the variant. Enabling it as a filter discarded 63 of 64 catalogue
  matches in a run whose null model put those matches at 49x chance.
- Gene concordance is evaluated against **both** breakends. A junction joins two
  genes and the catalogue peptide may be annotated to either; testing only the
  first side lost 8 of 64 matches, all translocations.
- Junction support is counted per **fragment**, not per alignment record
  (`COUNT_UNIQUE_FRAGMENTS`), with a mapping-quality floor
  (`MIN_READ_MAPQ = 20`, matching GRIDSS2). One insertion whose support read as
  7 alignments was 3 fragments, one contributing 4 records at MAPQ 3.
- The strong-evidence threshold is **5 junction reads**, aligned with pVACfuse;
  the previous 3 is retained as a visible lower tier.
- Match counts are reported at both grains — rows (peptide x SV) and distinct
  peptides — because one peptide can be produced by several SVs and mixing the
  two makes runs incomparable.
- Every output directory carries the criteria manifest and generator version
  that produced it, so a result cannot outlive the code behind it.

### Known limitations

Read evidence supports "the junction is transcribed", never "the peptide is
presented" — that needs immunopeptidome mass spectrometry. Gene concordance is
part of the credibility definition, so convergent peptides arising from a
different gene are excluded by construction. No liftover is performed: reference
and samples must share a genome build and annotation release.
