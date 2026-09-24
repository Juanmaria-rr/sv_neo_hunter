# External resources

Every input this pipeline does not produce itself, and the version of it that
the current results were built from. **None of these files is in this
repository**: they are either licensed, access-controlled, or too large. This
file records what to obtain and how to tell whether you have the same thing.

Where a resource cannot be pinned by checksum, that is stated rather than left
to look settled.

---

## Software

| resource | version | notes |
|---|---|---|
| Python | 3.10 | `environment.yml` pins the full environment |
| pysam | 0.24.0 | read counting from BAMs |
| pandas | 2.3.3 | |
| pyensembl | 2.10.1 | transcript annotation lookups |
| netMHCpan | **4.2e** | licensed (DTU Health Tech), academic use, redistribution not permitted. Obtain per machine; `tools/run_netmhcpan.py` locates it and explains the setup if absent. Earlier results in this project used **4.1** — the two are not interchangeable, see `docs/`. |

Recreate the environment with:

    conda env create -f environment.yml

## Variant calls (upstream of this pipeline)

| resource | version |
|---|---|
| SV caller | **Esvee 1.0.2** (from the VCF header `##EsveeVersion`) |
| Purity/ploidy and SV post-processing | **PURPLE 4.1** (`##purpleVersion`) |
| Reference genome | GRCh38 |

## Annotation

| resource | version | how it is pinned |
|---|---|---|
| Ensembl | **release 115**, GRCh38 | `pyensembl install --release 115 --species homo_sapiens`. The proteome used for the self test is the same release's `pep.all`. |
| gnomAD-SV | **v4.1** | matched by 0.5 reciprocal overlap; the nine ancestry groups in `criteria.GNOMAD_POPULATIONS` |
| Common fragile sites | two catalogues, conservative and permissive | passed as BED paths; their genome footprints are measured at run time and written into `CFS_ANNOTATION_PROVENANCE.tsv`, because a hit rate cannot be read without them |

## Panel of normals — **size unknown, and it matters**

The PON arrives inside the SV VCF as the `PON_COUNT` INFO field, defined there
only as *"PON count if in PON"*.

**How many samples the panel was built from is not documented anywhere we could
find**: not in the VCF header, not in the hmftools resource documentation, not
in the Esvee documentation. Two consequences, both of which shaped the columns
this pipeline emits:

- A derived `pon_fraction` existed and was **withdrawn**. Its denominator was
  the largest `PON_COUNT` in the same VCF — a lower bound on the panel size, and
  one that varies per sample because it depends on how many junctions that
  sample called. The same count gave 8.4% in one line and 27.6% in another.
- `PON_COUNT` may not count individuals at all. The largest value observed here
  is **11,912**, above every published size of the cohort the panel is built
  from (2,399 patients in 2019; 4,375 in 2023; over 8,000 in the database). It
  is more likely counting breakpoint observations than carriers.

`pon_count` is therefore carried as a raw count and is sound **as an ordinal** —
thousands against a minimum of 1 does separate recurrent from private — and not
as a frequency. Population frequency rests on gnomAD, which publishes allele
frequencies with known denominators per ancestry group.

**If this is ever resolved** (an issue on hmftools would settle it), a real
`pon_fraction` becomes computable and comparable, and the population verdict can
rest on both sources instead of one.

## Sample data

Access-controlled and never in this repository: the SV VCFs, the RNA and DNA
BAMs, the per-patient HLA genotypes, and the reference cohort catalogue. Paths
live in `config/<run>.yaml`, which is gitignored; `config/template.yaml` shows
the shape.

Each run copies the stage-6/8 tables it depends on into `inputs/` with sha256
checksums, so a later run can prove it used the same inputs:

    python tools/snapshot_candidate_inputs.py --check --out-dir <run>/inputs
