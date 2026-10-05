# HLA typing of the samples

`presentable` and every column derived from it are computed against a class I
genotype. This step measures that genotype from each sample's own reads, so the
genotype used downstream has a recorded origin, and so a derived sample that has
lost an HLA haplotype is not credited with alleles it no longer carries.

## Tool and version

**LILAC 1.6** (hmftools), the version run by nf-core/oncoanalyser 2.0.0. LILAC is
also how the reference cohort's patients are typed, so using it here keeps the
sample side and the cohort side of every presentability comparison on the same
method. `tools/run_lilac.py` reproduces oncoanalyser's `LILAC_CALLING`
subworkflow step for step.

## Method

### Why the reads are realigned first

On an ALT-aware GRCh38 (such as `GRCh38_masked_exclusions_alts_hlas`, HMF 25.1),
reads from the HLA genes are distributed between chr6 and the ~500 HLA and alt
contigs. LILAC works on chr6 coordinates. So, per DNA BAM:

| step | command | why |
|---|---|---|
| slice | `samtools view --regions-file grch38_alt.plus_homologous.bed` | chr6 HLA region, the HLA/alt contigs and homologous loci |
| fastq | `samtools sort -n` → `samtools fastq` | recover read pairs |
| realign | `bwa-mem2 mem -Y` against chr6 only → sort → index | put every HLA read back on chr6 |

The chr6 sequence must be the one the BAMs were aligned against. For HMF 25.1 it
can be taken from the published FASTA with an HTTP byte-range request using its
`.fai` (offset, line length), without downloading the whole genome.

### How each sample is typed

| sample | LILAC inputs | what it reports |
|---|---|---|
| lineage root (`parent: null`) | `-reference_bam` own realigned, `-rna_bam` own | the germline genotype, with per-allele DNA and RNA fragment support |
| derived sample | `-reference_bam` **root's** realigned, `-tumor_bam` own realigned, `-rna_bam` own, `-purple_dir` own | the same genotype plus, per allele, support in the derived sample, its copy number and somatic variants — i.e. loss of an allele relative to the root |

Derived samples are typed against the lineage **root**, not their immediate
parent, because the root is the germline their somatic (PURPLE) calls were made
against.

### Outputs

`--out-dir` is required and has no default; the tool refuses to write over an
existing `hla_genotypes.tsv`.

| file | content |
|---|---|
| `lilac/<sample>/` | LILAC's own output, untouched (`*.lilac.tsv` per allele, `*.lilac.qc.tsv`) |
| `hla_genotypes.tsv` | all samples collected: one row per (sample, allele), LILAC's columns plus `mode`, `reference_sample`, `qc_status` |
| `realigned/<sample>.hla.realigned.bam` | step-3 output, reusable with `--reuse-realigned` |
| `PROVENANCE.tsv` | LILAC version, resources, realignment reference, per-sample inputs, code commit |

## Reading the result

- **DNA decides the genotype; RNA is support.** An allele with DNA fragments and
  no RNA fragments is typed but may not be expressed. Expression matters for
  presentation; the genotype does not change with it.
- **A derived sample's allele with tumour copy number near 0** has been lost.
  Its peptides cannot be presented by that allele in that sample, whatever the
  parental genotype says.
- **`qc_status` other than `PASS`** means LILAC itself flags the call (low
  coverage, unmatched fragments, discordant pairs). Read its QC file before using
  the genotype.

## Limitations

- RNA BAMs are used as aligned (as oncoanalyser does). RNA reads that the aligner
  placed on HLA/alt contigs are not seen, so RNA support is an underestimate.
- Typing resolution is LILAC's: two-field (protein-level) class I alleles. Class
  II is not typed.
- Allele-specific copy number needs `purple_dir`; without it a derived sample
  is typed but loss is not assessed.

## Running

```bash
conda env create -f envs/hla_typing.yml          # bwa-mem2, samtools
python tools/run_lilac.py --config config/<run>.yaml --out-dir <new dir> --realign-only
python tools/run_lilac.py --config config/<run>.yaml --out-dir <same dir> --reuse-realigned
```

Configuration: the `resources.lilac` block and per-sample `purple_dir` in
[`config/template.yaml`](../config/template.yaml).
