# Running, reading and reviewing a run

How to go from a config to tables you can audit, and how to check that a run
did what it claims. See [CONFIGURATION.md](CONFIGURATION.md) for the config
itself.

## What a run looks like

Console output is one block per sample and branch, with the attrition printed as
it happens:

```
=== parental_PON10  [germline] ===
  stage1: 17,953 records -> 2,286 admitted -> 1,143 junctions
  stage2 [neosv]: 24,548 peptide rows -> 20,710 unique
  stages3-5: 64 matches -> 57 credible -> 13 events
  stage6: 4 HC, 9 private (PON+population)
  stage8: tiers {'NONE': 12, 'WEAK': 1}
  null: 64 matches from 20,710 candidates (3.09/1,000); null 1.3 ± 0.8,
        p < 0.001, 49x over null
```

`funnel.tsv` is the same cascade as a table, with a fixed row order so two
samples can be read side by side:

| step | description | n |
|---|---|---|
| `vcf_records` | SV records in the raw VCF | 17,953 |
| `admitted` | after FILTER / pairing / population filters | 2,286 |
| `junctions` | paired breakends collapsed to junctions | 1,143 |
| `candidate_peptides` | unique candidate neopeptides | 20,710 |
| `matches_identical` | identical to a reference peptide (rows: peptide × SV) | 64 |
| `matches_unique_peptides` | distinct peptides among those matches | 64 |
| `matches_gene_concordant` | and broken in the same gene | 64 |
| `matches_credible` | and high-complexity, non-self | 57 |
| `events` | credible matches collapsed to genomic events | **13** |
| `events_private` | PON and population-frequency clean | 9 |
| `events_hc` | high-confidence SV calls | 4 |
| `events_rna_supported` | junction-crossing reads in RNA | 1 |
| `events_hc_and_rna` | high-confidence AND transcribed (privacy NOT applied) | 0 |
| `events_private_hc_and_rna` | private AND high-confidence AND transcribed | **0** |

Read the last three rows carefully: they are **overlapping sets, not a chain**.
All three start from the 13 events, and one event can fail more than one
criterion. `events_hc_and_rna` deliberately excludes privacy, so it is not the
headline figure — an event can be confidently called and transcribed while being
a germline polymorphism carried by most of the population. Quote
`events_private_hc_and_rna`.

Note also that the counting unit changes along the cascade: breakends →
junctions → peptide rows → genomic events. A peptide count read as a finding
count overstates the evidence by up to ~40×, which is why `events` is the
reported unit.

The numbers above come from access-controlled inputs and are shown to make the
output shape concrete; they are not reproducible from a clone. See
[Reproducing results](#reproducing-results).


## Reproducing results

**The inputs of the published analysis are not public.** The reference catalogue
is access-controlled, and the sample VCFs and BAMs belong to the group that
generated them. Cloning this repository therefore does not reproduce our figures;
it reproduces the *method*. Concretely:

| What | Availability |
|---|---|
| This code, thresholds, and the vendored generator | in the repository |
| gnomAD-SV v4.1 sites VCF (1.7 GB) | public download |
| Ensembl reference via pyensembl | public download |
| Reference peptide catalogue | access-controlled (cohort data agreement) |
| Sample SV calls, RNA and DNA alignments | belong to the originating group |

To confirm that an installation computes what ours does, run the integration
test. It executes the whole chain on a small synthetic dataset shipped in
`tests/data/` and asserts exact values:

```bash
python tests/test_integration.py
```

A pass means the environment reproduces the reference behaviour. A failure means
it does not — find out why before running real data. The fixture deliberately
contains the cases that have broken this pipeline: an insertion-driven event
whose span misrepresents its size, an inter-chromosomal junction with no defined
size, a panel-of-normals record that must be rejected, and a rejected breakend
sharing a coordinate with an accepted one.

To run the analysis on your own data, copy `config/template.yaml`, point it at a
peptide catalogue and a set of samples, and follow [Quick start](#quick-start).

## Reviewing a run

| Tool | Question it answers |
|---|---|
| `tools/build_master_table.py` | two tables of **values, not verdicts**: `master_sv.tsv` (one row per admitted junction) and `master_peptides.tsv` (one row per matched peptide), each with a column dictionary. `--runs`/`--suffix` restrict the combined table to the runs a given report covers |
| `tools/build_report.py` | a technical results report per sample, in markdown, generated from its outputs |
| `tools/build_html_report.py` | the same run as a **filter-by-filter validation page**: each filter with what it measures, how, why, its threshold and the threshold's justification, beside the count it removed and the names of what it removed. `--all` for every run |
| `tools/build_summary_html.py` | renders a `docs/` markdown document to a self-contained, shareable HTML page. The markdown stays canonical; the page is generated, never edited |
| `tools/build_focused_summary.py` | derives a single-branch version of a document — drops branch-comparison blocks, removes the other branch's table columns, and can drop whole sections by heading. Avoids a second hand-written document that would diverge from the first |
| `tools/sync_event_tables.py` | fills the events tables in `docs/` from the run outputs, so a results table in prose cannot drift from the run it describes. `--check` gates a commit |
| `tools/check_docs.py` | verifies the documentation still describes the code; exits 1 if stale |
| `tools/column_meanings.py` | one written explanation per master-table column, emitted into `<table>_column_dictionary.tsv`. A column with no entry is marked `UNDOCUMENTED` in the artefact and fails the test suite, so a new column cannot ship unexplained |
| `tools/analyse_locus_vs_peptide.py` | measures the gap between locus-level and peptide-level recurrence, and separates its causes. Answers "how can junctions match while peptides do not?" with numbers |
| `tools/check_strand_bias.py` | strand composition against a length-weighted background (see open questions) |
| `tools/inspect_insertion.py` | the individual reads behind an insertion call, with mapping quality, duplicate flag and in-read offset, against a background rate — is this support real or a duplicate stack? |
| `tools/compare_generators.py` | do two peptide generators produce interchangeable output? Run before swapping one. |
| `tools/strand_stratified_comparison.py` | tests *why* patch 002 removed catalogue matches instead of adding them, by stratifying the losses on the strand pair of each fusion. The plus/plus class is the internal control: if a shared methodological artefact explains the losses, that class must be untouched |
| `tools/compare_runs.py` | what changed between two output directories, on the quantities a generator defect was expected to move. Written so a correction can be *measured* rather than silently applied — replacing the old results answers the wrong question |
| `tools/diagnose_minus_strand_cds.py` | sweeps a breakpoint through every region of a transcript and checks the 5' coding segment returned. Exits 1 if patch 002 is missing or reverted |
| `notebooks/pipeline_walkthrough_pyspark.ipynb` | every stage re-derived independently in PySpark and asserted against the pipeline. Change `SAMPLE` in the first cell to walk through any line. |

The notebook is generated from `notebooks/build_walkthrough.py`; edit that and
re-run it rather than editing the `.ipynb`.

