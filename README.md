# sv_neo_hunter

Recurrence testing for structural-variant-derived neoantigens.

Given a catalogue of neopeptides and a set of samples with SV calls, `svneo`
answers one question: **do those peptides reappear in those samples, and is the
recurrence real or an artefact?** An apparent match can arise from a genuine
shared genomic event, or from mapping noise, low-complexity sequence,
coincidental peptide convergence, or a common germline polymorphism. The pipeline
is built to tell them apart and to make the attrition at every step explicit.

> The repository is `sv_neo_hunter`; the Python package inside it is `svneo`, so
> every command below is `python -m svneo.run …`, run from the repository root.

**Start here:** [How it works](#how-it-works) ·
[Installation](#installation) · [Quick start](#quick-start) ·
[Design rationale](#design-rationale) · [Limitations](#limitations)

**Reference:** [Configuration](docs/CONFIGURATION.md) ·
[Peptide generation](docs/PEPTIDE_GENERATION.md) ·
[Running and reviewing](docs/USAGE.md) ·
[Threshold provenance](docs/PROVENANCE.md) ·
[MHC binding method](docs/METHOD_mhc_propagation.md)

## How it works

The analysis is defined by two inputs, and nothing in `src/svneo` names a
specific cohort, tissue or cell line:

```
REFERENCE   a catalogue of peptides to look for   (any cohort, any caller)
SAMPLES     anything with an SV VCF               (+ RNA and DNA if available)
```

| Stage | Module | Question |
|---|---|---|
| 1 | `vcf.py` | Which SV records are admissible, and what is the real lesion size? |
| **2** | **`generators/`** | What neopeptides could this sample's SVs produce? *(pluggable)* |
| 3–5 | `cross.py` | Which reference peptides match, how tightly, and how many *events* is that? |
| 6 | `confidence.py` | Is the SV call believable, and is it private to the sample? |
| 7–8 | `rna.py` | Is the locus transcribed, and do reads **cross** the junction? |
| 9 | `synthesis.py` | What survives, and which sample explains it? |
| — | `null_model.py` | How many matches would chance alone produce? |

Samples form a lineage via `parent:`, which controls execution order (parents
first) and attribution (the earliest sample carrying a breakpoint explains it).
For unrelated samples, leave `parent: null` and the lineage collapses to a flat
list.

A sample missing a modality is recorded as `NA` for the affected stages and
continues. Missing data is never scored as a negative result.

## Repository layout

```
src/svneo/            the pipeline; nothing here names a cohort or sample
    criteria.py       EVERY threshold, one constant per decision, each with why
    config.py         the config schema (Sample, Branch, Reference, RunConfig)
    vcf.py            stage 1   admission, breakend pairing, event size
    generators/       stage 2   peptide generation (the one pluggable stage)
    cross.py          stages 3-5  catalogue cross, sequence QC, event collapse
    confidence.py     stage 6   SV confidence, PON and gnomAD privacy
    rna.py            stages 7-8  expression, junction-crossing read evidence
    synthesis.py      stage 9   funnel, cross-sample comparison, attribution
    null_model.py     how many matches chance alone would produce
    run.py            the orchestrator; entry point

vendor/neosv/         NeoSV (MIT), vendored; two documented patches
config/               template.yaml + a complete worked example
tools/                reporting and verification scripts
tests/                unit tests + an integration test on a synthetic fixture
docs/                 PROVENANCE.md, METHOD_mhc_propagation.md
notebooks/            an independent PySpark re-derivation of every stage
results/              run outputs (git-ignored; regenerated, not tracked)
results_<variant>/    a variant of the same analysis (a patch applied, a
                      different catalogue) kept in its own directory, so a
                      change can be measured rather than silently applied.
                      Each carries its own PROVENANCE.md
```

`README.md` (this file) is how to run it.
[`docs/PROVENANCE.md`](docs/PROVENANCE.md) derives every non-obvious threshold in
`src/svneo/criteria.py` from the observation that set it, so a considered
decision can be told from an arbitrary one.

## Installation

```bash
git clone https://github.com/Juanmaria-rr/sv_neo_hunter.git
cd sv_neo_hunter
conda create -n svneo python=3.10 -y && conda activate svneo
pip install -r requirements.txt
```

The peptide generator ships with the repository, so no separate install is
needed. One reference download is required:

```bash
export PYENSEMBL_CACHE_DIR=$PWD/pyensembl_cache
pyensembl install --release 115 --species homo_sapiens     # ~3.2 GB, once
```

> **The Ensembl release used to generate peptides must match the one the
> reference catalogue was built with.** Peptide sequences depend on the
> annotation; crossing catalogues built on different releases compares strings
> that were never comparable.

<details>
<summary><b>Version constraints that are not optional</b></summary>

| Package | Constraint | Why |
|---|---|---|
| **pyensembl** | **> 2.3.13** | The vendored generator carries a patch that is correct only for releases *after* the coding-sequence accessor changed. On ≤ 2.3.13 the patch must be reverted, or every fusion is one amino acid out of frame. See `vendor/patches/001-*`. |
| Python | ≥ 3.10 | union types in annotations |
| PySpark | 3.5.x on **Java 11** | review notebook only; PySpark 4.x requires Java 17 |

Reference environment used in development: Python 3.10.20, pandas 2.3.3,
pysam 0.24.0, **pyensembl 2.10.1**, PySpark 3.5.3, OpenJDK 11.0.26. Every run
records its own versions in `summary.json`.

</details>

## Quick start

```bash
# 1. validate the config without computing anything
python -m svneo.run --config config/my_run.yaml --validate-only

# 2. run it
python -m svneo.run --config config/my_run.yaml --out-dir results

# 3. turn the outputs into auditable tables
python tools/build_master_table.py --results-dir results
```

`config/lineage.example.yaml` is a complete worked example with every field
filled in; `config/template.yaml` is the same with the values blank. See
[Configuration](docs/CONFIGURATION.md) for what each field does, and
[Running and reviewing](docs/USAGE.md) for what a run prints and how to check it.

## Outputs

Per sample and branch, under `<out_dir>/<sample>[_<branch>]/`:

| File | Content |
|---|---|
| `stage1_junctions.tsv` | admitted junctions with type-aware `event_size` and `insert_len` |
| `<sample>.admitted.vcf` | the admitted call set handed to stage 2 |
| `stage3_matches.tsv` | peptide matches with concordance and QC flags |
| `stage3_gene_svtype.tsv`, `stage3_proximity.tsv` | cross levels 2 and 3 |
| `credible_events.tsv` | **the reported unit** — one row per genomic event |
| `stage7_expression.tsv`, `stage7_background.json` | expression with its background gradient |
| `stage8_rna_evidence.tsv` | per-breakend counts, test applied, evidence tier |
| `stage6_8_all_junctions.tsv` | the same evidence for **every** admitted junction, not only those that matched (`criteria.EVALUATE_ALL_JUNCTIONS`) |
| `stage3_junction_recurrence.tsv` | per junction: distance to the nearest patient breakpoint, and what was found there |
| `master_sv.tsv`, `master_peptides.tsv` | the auditable tables — values, not verdicts |
| `funnel.tsv` | the 14-row evidence funnel, fixed row order so samples are comparable line by line |
| `summary.json` | counts, null model, generator provenance, full threshold manifest |

Cross-sample, under `<out_dir>/`: `comparison.tsv`, `null_model_rates.tsv`,
`attribution.tsv`, the combined `master_sv.tsv` / `master_peptides.tsv` with their
column dictionaries, and `reports/<sample>.md`.

**Was a junction seen in patients?** `master_sv.tsv` carries `patient_evidence`,
naming the strongest level of recurrence each junction reaches — `identical_peptide`,
`same_gene_and_svtype`, `same_gene`, `breakpoint_within_1kb|10kb|100kb`, or
`not_seen_in_patients` — alongside the unthresholded distance in
`patient_bp_dist_bp`. The three levels are not the same kind of evidence: a shared
peptide is a shared consequence, a shared gene and SV type a shared mechanism, and
a nearby breakpoint only a shared locus.

**Population frequency** is carried for all nine gnomAD ancestry groups
(`gnomad_af_afr` … `gnomad_af_sas`) plus `gnomad_af_popmax`. Which one is the
right reference depends on the donor's ancestry, which the pipeline does not know;
`criteria.GNOMAD_AF_FIELD` selects what the automatic verdict uses, and every
group's frequency is in the table so it can be re-judged without re-running.

Read `comparison.tsv` on events and per-1,000-candidate rates, not on peptide
counts: a single locus seen through a sliding window can produce dozens of
matches.

## Design rationale

Five rules distinguish this implementation from a naive one. Each addresses a
failure mode that produces confident, wrong results, illustrated below with
observed cases.

**1. Only junction-crossing reads count as RNA evidence.** Coverage is context, and
a soft-clip shows that a read *ends* at a breakpoint, not that it crosses. An
inter-chromosomal candidate scored "5 split reads" on clipped reads whose
supplementary alignments landed 3–156 Mb from the partner breakend; no read
connected the two loci at all.

**2. Coverage is summarised across breakends with `min`, never `max`.** An event is in
transcribed territory only if both ends are. Taking the maximum reported a locus
carrying a handful of reads as having thousands — the count belonged to its highly
expressed partner, three orders of magnitude away.

**3. An event with no valid test is `UNTESTABLE`, not `NONE`.** `NONE` asserts
"tested and negative". A 2 bp deletion cannot produce a CIGAR `N` gap, because
aligners do not emit skips shorter than their minimum intron; scoring it `NONE`
fabricates a result. The test is selected from event geometry: chimeric evidence
for inter-chromosomal junctions, gap matching for intra-chromosomal events above
the floor, and inserted-sequence matching for insertion-driven events, which the
breakend span misrepresents.

**4. Privacy requires two independent filters.** A panel of normals has false
negatives for inherited variation: candidates with `PON_COUNT` of 1 and 3 carried
population allele frequencies of 0.287 and 0.263. An event must clear both the
panel and a population-frequency threshold. A well-called, well-expressed,
read-supported deletion present in thousands of unrelated normals is real — and
disqualified, because the property that matters is not whether the call is
correct but whether it is the sample's own.

**5. Counting is done on events, never peptide rows.** One SV yields up to ~40
overlapping sliding-window peptides. In one case 21 apparent matches were a
single 53 bp deletion, a 21-fold overstatement. A related trap: deduplicating on
a per-peptide value such as an odds ratio leaves window peptides separate.

Two further points shape interpretation. Raw match counts are not comparable
between samples whose candidate universes differ by orders of magnitude, so every
sample reports a per-1,000-candidate rate and a permutation null alongside the
raw count. And a funnel that reduces a four-figure candidate count to a
single-digit set is the expected outcome, not a failure: a 28-team benchmark
found ~6% of top-ranked neoantigen predictions validate functionally. What
matters is where each order of magnitude is lost.

Full derivation, including which observation set each threshold, is in
[`docs/PROVENANCE.md`](docs/PROVENANCE.md).

## Limitations

- **Presentation is out of scope.** Read evidence supports "the junction is
  transcribed". Claiming presentation requires immunopeptidome mass spectrometry.
- **Gene concordance is part of the credibility definition**, so a genuine
  convergent peptide arising from a different gene is excluded by construction.
  This trades sensitivity for specificity and should be stated whenever
  "credible" counts are quoted.
- **A negative RNA result is not absence of the lesion.** The gene may not be
  expressed in that sample, and nonsense-mediated decay of an aberrant transcript
  is a real possibility.
- **Copy-number fields are excluded from credibility judgements** by default,
  because they depend on a caller's purity/ploidy fit that is unreliable for
  clonal samples.
- **No liftover is performed.** The reference catalogue and every sample VCF must
  be on the same genome build and annotation release.
- **A defect in the vendored generator, corrected by patch 002 but present in every pre-patch figure.** On **minus-strand
  transcripts the 5' coding segment is wrong in every region**, while every
  plus-strand region is correct. An intronic breakpoint
  returns an empty head, so the fusion becomes the 3' partner alone, is flagged
  `Start-loss`, places the junction at residue ~0 and can never yield a
  junction-spanning peptide. An exonic breakpoint returns the wrong length, and
  near the end of a transcript returns **too much** — adding sequence the gene
  does not contribute, which the sliding window turns into peptides. The visible
  symptoms are a large excess of `Start-loss` calls on the minus strand and a
  deficit of junction-spanning peptides among matches. Root cause, evidence and
  the fix are in
  [`vendor/patches/002-minus-strand-cds-order.md`](vendor/patches/002-minus-strand-cds-order.md);
  reproduce with `python tools/diagnose_minus_strand_cds.py`, which exits 1 if
  the patch is missing or reverted. Keep patched and unpatched runs in separate
  output directories so the difference can be measured.
- **MHC binding is predicted with netMHCpan-4.2e, the earlier analysis used
  4.1.** Among peptide-allele calls that bind under at least one version, only
  69.4% agree (measured, not assumed — `tools/compare_netmhcpan_versions.py`).
  Counts are robust to this, peptide identities are not. See
  [`docs/METHOD_mhc_propagation.md`](docs/METHOD_mhc_propagation.md).

## MHC binding

Binding prediction is **not** part of the recurrence test, which is a sequence
identity test and HLA-independent by construction. It is a separate, optional
layer, needed only to compare these counts with an analysis that had already
filtered for presentation.

| script | what it does |
|---|---|
| `tools/run_netmhcpan.py` | netMHCpan over a peptide list; one row per (peptide, allele), no collapsing to a verdict |
| `tools/propagate_binding_to_events.py` | carries binding up to peptides and events, in two branches |

**The two branches are not interchangeable.** `binds_panel` asks whether a
peptide could be presented by *someone* in the cohort — an upper bound, never a
result. `binds_autologous` asks whether it binds an allele of a patient who
actually carries it, which is the only quantity comparable with a per-sample
presentation filter. Quoting the first against such a count inflates it.

Method, denominators and limitations:
[`docs/METHOD_mhc_propagation.md`](docs/METHOD_mhc_propagation.md).

## Testing

```bash
python tests/test_svneo.py        # unit tests: the rules that must not regress
python tests/test_integration.py  # the whole chain on the synthetic fixture
```

The suite pins the decisions above rather than pursuing coverage: each test
guards a rule whose violation produces a plausible but wrong result.

## Licence

**MIT** — see [`LICENSE`](LICENSE).

The vendored dependency, NeoSV, is also MIT (Copyright (c) 2022 Yang Shi); its
licence ships as [`vendor/LICENSE.NeoSV`](vendor/LICENSE.NeoSV) and applies to
`vendor/neosv/`. Modifications to it are limited to two documented patches under
[`vendor/patches/`](vendor/patches/).

The reading-frame bug that patch fixes was identified in
[NeoSV-Trace](https://github.com/winterga/NeoSV-Trace) by Greyson Wintergerst;
no code from that fork is redistributed here.

## Citation

If this pipeline contributes to published work, cite this repository, the
peptide generator it vendors — **NeoSV** (Shi, Jing & Xi, *Genome Biol* 24:169,
2023), see [`vendor/VENDOR.md`](vendor/VENDOR.md) — and the reference peptide
catalogue used. No code from NeoSV-Trace is redistributed here; that fork is
credited for identifying the reading-frame bug, not vendored.
