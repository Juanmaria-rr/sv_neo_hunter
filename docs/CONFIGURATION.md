# Configuration

A run is defined entirely by a config file. `config/lineage.example.yaml` is a
complete worked example; `config/template.yaml` is the same with every field
blank. Validate one without computing anything:

```bash
python -m svneo.run --config config/my_run.yaml --validate-only
```


`config/template.yaml` documents every field. The essentials:

```yaml
run_name: my_run
peptide_generator: neosv              # stage-2 backend
out_dir: results

genome:
  build: GRCh38                       # recorded, not enforced: no liftover
  ensembl_release: 115                # MUST match the catalogue's release

criteria:                             # optional per-run threshold overrides
  PON_MAX: 10                         # any constant from src/svneo/criteria.py

reference:
  name: my_catalogue                  # appears in the run manifest
  peptides: /path/to/reference_peptides.tsv
  peptide_column: neoantigen          # required: the peptide sequence
  gene_column: gene1                  # null disables gene concordance
  svtype_column: svtype               # null disables cross level 2
  proteome: /path/to/proteome.fa.gz   # null records the self test as NA
  pool: /path/to/catalogue_with_coords.tsv   # optional; enables cross level 3
  pool_peptide_column: neoantigen

samples:
  - name: sample_A
    vcf: /path/to/sample_A.sv.vcf.gz  # required; everything else is optional
    vcf_kind: somatic                 # germline | somatic
    parent: null                      # name of another sample, or null
    rna_bam: /path/to/sample_A.bam    # stages 7-8; NA without it
    isofox_dir: /path/to/isofox/A     # transcript-level expression context
    isofox_prefix: sample_A.isf       # filename prefix inside isofox_dir
    dna_bam: /path/to/sample_A.dna.bam  # carried for follow-up; no stage reads it yet
    notes: |                          # free text, echoed into summary.json
      Anything a reader needs to interpret this sample.
```

A sample missing a modality gets `NA` for the affected stages and continues, so
the minimum viable sample is a name and an SV VCF.

<details>
<summary><b>External resources</b></summary>

```yaml
resources:
  pyensembl_cache: /path/to/pyensembl_cache        # required by stage 2
  gnomad_sv: /path/to/gnomad.v4.1.sv.sites.vcf.gz  # strongly recommended
  gnomad_helper: /path/to/annotator.py             # optional
  neosv_path: /path/to/an/external/neosv           # optional; vendored by default
  peptides_from:                                   # optional; skip stage 2
    sample_A: /path/to/previous/run/sample_A
```

| Key | Effect if absent |
|---|---|
| `pyensembl_cache` | stage 2 falls back to pyensembl's default location; set it explicitly so runs are reproducible across machines |
| `gnomad_sv` | **privacy degrades silently to the panel of normals alone.** Every `gnomad_af_*` column is written as `NA` and `privacy_note` records `"gnomAD not evaluated — PON only"`. A panel has false negatives for inherited variation, so this is a real loss of specificity, not a cosmetic one — see design rationale 4 |
| `gnomad_helper` | the reciprocal-overlap join runs directly; the helper only exists to reuse a site's existing annotator and its local-copy discovery |
| `neosv_path` | the vendored copy under `vendor/neosv/` is used, which is the intended default |
| `peptides_from` | stage 2 runs normally; point it at a previous run's output prefix to reproduce a run without regenerating peptides |

Download gnomAD-SV v4.1 sites (~1.7 GB) from
[gnomad.broadinstitute.org](https://gnomad.broadinstitute.org/downloads#v4-structural-variants).
It must be on the same genome build as the samples; no liftover is performed.

The three cross levels answer progressively weaker questions, and the config
controls which are available:

| Level | Needs | Question |
|---|---|---|
| 1 identical peptide | `peptide_column` | a shared **consequence** |
| 2 gene + SV type | `gene_column`, `svtype_column` | a shared **mechanism** |
| 3 breakpoint proximity | `pool` (a catalogue **with coordinates**) | a shared **locus** only |

`pool` is separate from `peptides` because the two differ in what they carry: the
curated peptide set may have no coordinates, while the larger pool has
breakpoints but is not peptide-filtered. Without `pool`, level 3 and the
`patient_evidence` proximity tiers are simply absent, not zero.

`vcf_kind` is a property of the sample, not of the filename: `germline` means
"the background genome shared with descendants" and receives the
population-frequency filters, `somatic` means "acquired relative to the parent"
and is assumed already panel-filtered by the caller.

</details>

<details>
<summary><b>Sensitivity branches</b></summary>

A branch re-runs the same samples under a different threshold, so that any
difference in the result is attributable to that threshold and nothing else.
Branches differ in stage-1 admission and — for the panel of normals only — in
whether the panel count votes on the privacy verdict at stage 6.

```yaml
branches:
  - name: PON10
    pon_max: 10                 # keep records with PON_COUNT < 10
  - name: noPON
    pon_max: null               # this pipeline's own PON threshold: off
    admit_panel_filtered: true  # also admit somatic FILTER=PON records
    pon_in_privacy: false       # PON reported, is_private = gnomAD only
```

A branch may also set `gnomad_max_af:` to override the population-frequency
threshold for that branch alone; `null` disables the population filter. Any
threshold in `criteria.py` can likewise be overridden per run from the config's
`criteria:` block, and every override is echoed into the run manifest.

**The panel of normals reaches a candidate twice, and a branch meaning "report
the panel, do not filter on it" has to relax both:**

| Where | Germline call set | Somatic call set |
|---|---|---|
| Admission | `PON_COUNT` is an INFO field this pipeline thresholds → `pon_max: null` is enough | the caller already decided and wrote `FILTER=PON`; those records are dropped for not being `PASS` **before** `pon_max` is consulted → needs `admit_panel_filtered: true` |
| Privacy (stage 6) | `is_private` combines the panel count with population frequency → `pon_in_privacy: false` removes the panel's vote | same |

Setting only `pon_max: null` is therefore a **no-op on somatic call sets**: it
produces output identical to the filtered branch under an unfiltered label. The
panel count is measured and carried into every table in both branches; what
changes is whether it removes anything.

Thresholds live in `src/svneo/criteria.py`, one constant per decision, each with
its rationale. They can be overridden per run from the config's `criteria:`
block; every override is echoed into the run manifest, so a result can never be
traced to a threshold that is not recorded beside it.

`config/lineage.example.yaml` is a complete worked example with every field
filled in and paths as placeholders. Local run configurations (`config/*.yaml`)
are git-ignored, since they carry absolute paths and the filenames of
access-controlled inputs.

</details>
