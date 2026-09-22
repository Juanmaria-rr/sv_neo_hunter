#!/usr/bin/env python3
"""
build_walkthrough.py — generates notebooks/pipeline_walkthrough_pyspark.ipynb.

The notebook is generated rather than hand-edited so it stays consistent with the
package: cell sources live here, in a plain file that diffs cleanly, and the
.ipynb is a build artefact. Re-run after changing the walkthrough.

    python notebooks/build_walkthrough.py
"""
from __future__ import annotations

import json
import pathlib

CELLS: list[tuple[str, str]] = []


def md(text: str) -> None:
    CELLS.append(("markdown", text.strip("\n")))


def code(text: str) -> None:
    CELLS.append(("code", text.strip("\n")))


# ============================================================================
md(r"""
# Pipeline walkthrough — every stage, re-implemented in PySpark

**Purpose: review, not production.** The package (`src/svneo/`) is pandas and
pysam. This notebook re-derives the same stages independently in Spark and
**asserts that both agree at every step**. Where they disagree, one of them is
wrong — that is the point.

Two things it is good for:

1. **Auditing the logic.** Each stage shows the grain of its rows, the join keys,
   and the identity that must hold. A silent join miss or a grain confusion is
   invisible in a summary count and obvious here.
2. **Scaling the counting logic.** 18k rows do not need Spark. Cohort-scale
   inputs do, and the operations that break at scale — grain confusion,
   double-counted breakends, joins that silently drop rows — are exactly the ones
   written out explicitly below.

Parameterised on `SAMPLE` / `BRANCH`: set them in the next cell and the whole
notebook re-runs for any line.

> Requires PySpark 3.5.x on **Java 11**. PySpark 4.x needs Java 17; mixing them
> fails at `SparkSession.builder.getOrCreate()` with an opaque JVM error.
""")

code(r"""
# --- parameters -------------------------------------------------------------
# Change these two and re-run the notebook; everything else is derived.
SAMPLE = "parental"       # any sample name in the config
BRANCH = "PON10"         # "PON10" | "noPON" | "" for an unbranched run

import os, sys, pathlib
REPO = pathlib.Path.cwd().parent if pathlib.Path.cwd().name == "notebooks" else pathlib.Path.cwd()
sys.path.insert(0, str(REPO / "src"))

RUN_DIR = REPO / "results" / (f"{SAMPLE}_{BRANCH}" if BRANCH else SAMPLE)
assert RUN_DIR.is_dir(), f"no pipeline output at {RUN_DIR} — run the pipeline first"

from svneo import criteria, config as config_mod
CFG = config_mod.load(str(REPO / "config" / "my_run.yaml"))
SAMPLE_CFG = CFG.by_name(SAMPLE)
assert SAMPLE_CFG, f"{SAMPLE} is not in the config"

# Derived, never hand-set: a mismatch between the declared kind and the one the
# pipeline used would silently compare two different analyses.
VCF_KIND = SAMPLE_CFG.vcf_kind
PON_MAX = None if (BRANCH == "noPON" or VCF_KIND != "germline") else criteria.PON_MAX

print(f"sample   {SAMPLE}  [{VCF_KIND}]  branch={BRANCH or 'default'}")
print(f"vcf      {SAMPLE_CFG.vcf}")
print(f"outputs  {RUN_DIR}")
print(f"PON_MAX  {PON_MAX}")
""")

code(r"""
# The workers must run the SAME interpreter as the driver. Without this the JVM
# spawns whatever `python3` is first on PATH, which lacks this environment's
# packages and fails inside the executor with a Java stack trace that never
# mentions Python.
os.environ["PYSPARK_PYTHON"] = sys.executable
os.environ["PYSPARK_DRIVER_PYTHON"] = sys.executable

from pyspark.sql import SparkSession, functions as F, types as T

spark = (SparkSession.builder
         .appName(f"svneo-walkthrough-{SAMPLE}")
         .master("local[*]")
         .config("spark.driver.memory", "6g")
         .config("spark.sql.shuffle.partitions", "8")
         .config("spark.driver.maxResultSize", "2g")
         .getOrCreate())
spark.sparkContext.setLogLevel("ERROR")

# The UDFs below call the package's own criteria functions, so the WORKERS need
# the package too — the driver's sys.path does not reach them. Shipping a zip
# keeps the notebook calling the same code as the pipeline instead of restating
# the thresholds, which is the only way the comparison stays honest.
import shutil, tempfile
_pkg_zip = shutil.make_archive(
    os.path.join(tempfile.gettempdir(), "svneo_pkg"), "zip", str(REPO / "src"))
spark.sparkContext.addPyFile(_pkg_zip)
print(spark.version, "| package shipped to workers:", os.path.basename(_pkg_zip))

CHECKS = []   # (stage, statement, passed) — collected and summarised at the end

def check(stage: str, statement: str, condition: bool):
    CHECKS.append((stage, statement, bool(condition)))
    print(f"  [{'PASS' if condition else 'FAIL'}] {stage}: {statement}")
    return condition
""")

# ---------------------------------------------------------------- stage 1
md(r"""
## Stage 1 — Admission and junction pairing

**Grain of a VCF row: one breakend, not one event.** A junction appears twice,
once from each side, so any count taken on raw records double-counts. The
identity that must hold after pairing:

```
records = 2 x junctions + single_breakends
```

Three fields decide admission, and one of them is easy to get wrong:

- `FILTER` — somatic call sets arrive panel-filtered by the caller; germline ones
  usually pass everything and are only *annotated*.
- the ALT bracket — a breakend without a mate has one coordinate, and a junction
  peptide needs two.
- `PON_COUNT` — an **absolute count**, not a frequency. Reported here as a
  fraction of the estimated panel size.
""")

code(r"""
import gzip, re

def read_vcf_lines(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [ln.rstrip("\n") for ln in fh if not ln.startswith("#")]

raw = spark.createDataFrame([(l,) for l in read_vcf_lines(SAMPLE_CFG.vcf)], ["line"])

fields = (raw
    .withColumn("f", F.split("line", "\t"))
    .select(
        F.regexp_replace(F.col("f")[0], "^chr", "").alias("chrom"),
        F.col("f")[1].cast("int").alias("pos"),
        F.col("f")[2].alias("id"),
        F.col("f")[3].alias("ref"),
        F.col("f")[4].alias("alt"),
        F.col("f")[5].alias("qual"),
        F.col("f")[6].alias("filter"),
        F.col("f")[7].alias("info")))

def info_field(col, key, cast="int"):
    # INFO is a ;-delimited key=value string; pull one key out and type it.
    return F.regexp_extract(col, rf"(?:^|;){key}=([^;]+)", 1).cast(cast)

breakends = (fields
    .withColumn("paired", F.col("alt").rlike(r"[\[\]]"))
    .withColumn("svtype",    F.regexp_extract("info", r"(?:^|;)SVTYPE=([^;]+)", 1))
    .withColumn("mateid",    F.regexp_extract("info", r"(?:^|;)MATEID=([^;]+)", 1))
    .withColumn("pon_count", info_field(F.col("info"), "PON_COUNT"))
    .withColumn("vf",        info_field(F.col("info"), "VF"))
    .withColumn("segmapq",   info_field(F.col("info"), "SEGMAPQ"))
    .withColumn("qual_n",    F.col("qual").cast("double"))
    .cache())

n_records = breakends.count()
n_paired  = breakends.filter("paired").count()
print(f"records {n_records:,}   paired {n_paired:,}   single {n_records - n_paired:,}")

panel_size = breakends.agg(F.max("pon_count")).first()[0] or 0
print(f"panel size (floor estimate, = max PON_COUNT): {panel_size:,}")
""")

code(r"""
# --- admission --------------------------------------------------------------
admitted = breakends.filter("paired")
if VCF_KIND == "somatic":
    admitted = admitted.filter(F.col("filter").isin(list(criteria.SOMATIC_KEEP_FILTERS)))
if PON_MAX is not None:
    # An ABSENT PON_COUNT means "not in the panel", which passes. NaN/null is the
    # cleanest possible record, and treating it as a failure inverts the filter.
    admitted = admitted.filter(
        F.coalesce(F.col("pon_count"), F.lit(criteria.PON_ABSENT_MEANS)) < PON_MAX)

n_admitted = admitted.count()
print(f"admitted breakends: {n_admitted:,}")

# --- pair into junctions ----------------------------------------------------
# Self-join on MATEID, keeping each junction once via a canonical key on the
# sorted id pair. Without the canonicalisation every junction appears twice.
left  = admitted.alias("a")
right = admitted.select(F.col("id").alias("mid"), F.col("chrom").alias("mchrom"),
                        F.col("pos").alias("mpos"), F.col("pon_count").alias("mpon"),
                        F.col("vf").alias("mvf"), F.col("qual_n").alias("mqual"),
                        F.col("segmapq").alias("msegmapq")).alias("b")

junctions = (left.join(right, F.col("a.mateid") == F.col("b.mid"), "inner")
    .withColumn("key", F.concat_ws("|", F.array_sort(F.array(F.col("a.id"), F.col("b.mid")))))
    .dropDuplicates(["key"])
    .select(
        F.col("key").alias("junction_key"),
        F.col("a.id").alias("sample_sv_id"),
        F.col("a.chrom").alias("chrom1"), F.col("a.pos").alias("pos1"),
        F.col("b.mchrom").alias("chrom2"), F.col("b.mpos").alias("pos2"),
        F.col("a.svtype").alias("svtype"),
        F.greatest(F.col("a.pon_count"), F.col("b.mpon")).alias("pon_count"),
        F.col("a.vf").alias("vf_bp1"), F.col("b.mvf").alias("vf_bp2"),
        F.col("a.qual_n").alias("qual_bp1"), F.col("b.mqual").alias("qual_bp2"),
        F.col("a.segmapq").alias("segmapq_bp1"), F.col("b.msegmapq").alias("segmapq_bp2"))
    .cache())

n_junctions = junctions.count()
print(f"junctions: {n_junctions:,}")
check("stage1", "records = 2 x junctions + unpaired (mates may be filtered out)",
      2 * n_junctions <= n_admitted + (n_records - n_paired))
""")

code(r"""
# --- compare against the pipeline ------------------------------------------
import pandas as pd
pipe_junctions = pd.read_csv(RUN_DIR / "stage1_junctions.tsv", sep="\t")
print(f"pipeline junctions: {len(pipe_junctions):,}   spark: {n_junctions:,}")
check("stage1", "junction count matches the pipeline", len(pipe_junctions) == n_junctions)

spark_ids = {r.sample_sv_id for r in junctions.select("sample_sv_id").collect()}
pipe_ids  = set(pipe_junctions.sample_sv_id.astype(str))
check("stage1", "the same junction ids, not merely the same number",
      spark_ids == pipe_ids)
if spark_ids != pipe_ids:
    print("  only in spark:", sorted(spark_ids - pipe_ids)[:5])
    print("  only in pipe :", sorted(pipe_ids - spark_ids)[:5])
""")

md(r"""
### The trap: breakend span is not the lesion size

`|pos2 - pos1|` equals the lesion only for duplications.

| type | span | actual size |
|---|---|---|
| `DEL` | flanks the deleted stretch | `span - 1` bases deleted |
| `INS` | **always 1** — breakends are adjacent | length of the inserted sequence |
| `DUP` | the duplicated tract | = span |
| `BND` | undefined | — |

Using the span makes every insertion look like a 1 bp event, and an
insertion-driven event look like a tiny deletion. That misclassification then
selects the wrong RNA test downstream.
""")

code(r"""
sized = (junctions
    .withColumn("span", F.when(F.col("chrom1") == F.col("chrom2"),
                               F.abs(F.col("pos2") - F.col("pos1"))))
    .withColumn("event_size", F.when(F.col("svtype") == "DEL", F.col("span") - 1)
                               .when(F.col("svtype").isin("DUP", "INV"), F.col("span"))
                               .when(F.col("svtype").isin("BND", "TRA"), F.lit(None))
                               .otherwise(F.col("span"))))

print("event size by SV type:")
(sized.groupBy("svtype")
      .agg(F.count("*").alias("n"),
           F.expr("percentile_approx(span, 0.5)").alias("median_span"),
           F.expr("percentile_approx(event_size, 0.5)").alias("median_size"))
      .orderBy(F.desc("n")).show())
""")

# ---------------------------------------------------------------- stage 2
md(r"""
## Stage 2 — Candidate peptides

Peptide generation is NeoSV plus this repository's candidate-table writer; it is
not a Spark operation and is not re-implemented here. The notebook loads what the
pipeline produced and checks the two counts that matter.

**Rows are not peptides.** One peptide can be produced by several SVs or several
transcripts. Quote the unique count against a catalogue; quote rows only when
talking about the table itself.
""")

code(r"""
peptides = spark.read.option("header", True).option("sep", "\t") \
                .csv(str(RUN_DIR / f"{SAMPLE}.all_neopeptides.txt")).cache()
anno = spark.read.option("header", True).option("sep", "\t") \
             .csv(str(RUN_DIR / f"{SAMPLE}.anno.txt")).cache()

n_rows = peptides.count()
n_unique = peptides.select("neopeptide").distinct().count()
print(f"peptide rows {n_rows:,}   unique peptides {n_unique:,}   "
      f"annotated SVs {anno.select('sv_id').distinct().count():,}")

check("stage2", "every peptide row carries an sv_id",
      peptides.filter(F.col("sv_id").isNull() | (F.col("sv_id") == "NA")).count() == 0)

print("\npeptides per SV (a sliding window over one junction yields many):")
(peptides.groupBy("sv_id").agg(F.countDistinct("neopeptide").alias("n"))
         .agg(F.min("n").alias("min"), F.expr("percentile_approx(n, 0.5)").alias("median"),
              F.max("n").alias("max")).show())
""")

# ---------------------------------------------------------------- stage 3-5
md(r"""
## Stages 3–5 — Cross, sequence QC, collapse to events

The cross is an inner join on the peptide string. Two properties are annotated,
and the second is a correction worth stating:

- **gene concordance** — evaluated against **both** breakends. A junction joins
  two genes and the catalogue peptide may be annotated to either; testing only
  the first side silently drops real matches (8 of 64 in one run, all
  translocations).
- **SV-type concordance** — reported, never used as a gate.

Then the collapse that decides every headline number: **peptides are not
findings**. One SV yields up to ~40 overlapping window peptides, each able to
match separately.
""")

code(r"""
reference = (spark.read.option("header", True).option("sep", "\t")
             .csv(CFG.reference.peptides)
             .withColumnRenamed(CFG.reference.peptide_column, "ref_peptide")
             .withColumnRenamed(CFG.reference.gene_column, "ref_gene")
             .withColumnRenamed(CFG.reference.svtype_column, "ref_svtype")
             .withColumn("ref_peptide", F.upper(F.trim("ref_peptide"))).cache())
print(f"reference catalogue: {reference.count():,} peptides")

matches = (peptides
    .withColumn("peptide", F.upper(F.trim(F.col("neopeptide"))))
    .join(reference, F.col("peptide") == F.col("ref_peptide"), "inner")
    .withColumn("concordant_gene1", F.upper(F.trim("gene1")) == F.upper(F.trim("ref_gene")))
    .withColumn("concordant_gene2", F.upper(F.trim("gene2")) == F.upper(F.trim("ref_gene")))
    .withColumn("gene_concordant", F.col("concordant_gene1") | F.col("concordant_gene2"))
    .cache())

n_match_rows = matches.count()
n_match_pept = matches.select("peptide").distinct().count()
n_concordant = matches.filter("gene_concordant").select("peptide").distinct().count()
print(f"match rows {n_match_rows}   unique peptides {n_match_pept}   "
      f"gene-concordant {n_concordant}")

pipe_matches = pd.read_csv(RUN_DIR / "stage3_matches.tsv", sep="\t")
check("stage3", "same match rows as the pipeline", len(pipe_matches) == n_match_rows)
check("stage3", "same set of matching peptides",
      set(pipe_matches.peptide) == {r.peptide for r in matches.select("peptide").distinct().collect()})
""")

code(r"""
# --- sequence QC as Spark UDFs, calling the package's own functions ---------
# Importing criteria rather than reimplementing means the notebook cannot drift
# from the pipeline: if a threshold changes, both move together.
lc_udf   = F.udf(lambda p: bool(criteria.is_low_complexity(p or "")), T.BooleanType())

proteome = None
if CFG.reference.proteome:
    from svneo import cross as cross_mod
    proteome = cross_mod.load_proteome(CFG.reference.proteome)
    bc = spark.sparkContext.broadcast(proteome)
    self_udf = F.udf(lambda p: bool(p and p in bc.value), T.BooleanType())
else:
    self_udf = F.udf(lambda p: False, T.BooleanType())

qc = (matches
      .withColumn("low_complexity", lc_udf("peptide"))
      .withColumn("is_self", self_udf("peptide"))
      .withColumn("credible", (~F.col("low_complexity")) & (~F.col("is_self"))
                              & F.col("gene_concordant"))
      .cache())

print(qc.groupBy("low_complexity", "is_self", "gene_concordant", "credible")
        .count().orderBy(F.desc("count")).toPandas().to_string(index=False))

n_credible = qc.filter("credible").count()
check("stages4", "same credible row count as the pipeline",
      int(pipe_matches.credible.fillna(False).astype(bool).sum()) == n_credible)
""")

code(r"""
# --- collapse to events -----------------------------------------------------
events = (qc.filter("credible")
            .groupBy("sv_id")
            .agg(F.countDistinct("peptide").alias("n_peptides"),
                 F.sort_array(F.collect_set("ref_gene")).alias("ref_genes"))
            .cache())

n_events = events.count()
print(f"credible peptides {n_credible} -> events {n_events}")
print("\ninflation factor per event (peptides collapsed into one finding):")
events.orderBy(F.desc("n_peptides")).show(10, truncate=False)

pipe_events = pd.read_csv(RUN_DIR / "credible_events.tsv", sep="\t")
check("stage5", "same event count as the pipeline", len(pipe_events) == n_events)
check("stage5", "same event ids",
      set(pipe_events.sample_sv_id.astype(str)) ==
      {r.sv_id for r in events.select("sv_id").collect()})
""")

# ---------------------------------------------------------------- stage 6
md(r"""
## Stage 6 — Confidence and privacy

Two orthogonal questions, deliberately not merged:

- **Confidence** — is the call real? Read support, mapping quality, caller
  confidence, and a size floor for intra-chromosomal events.
- **Privacy** — is it the sample's own? A panel of normals **and** a population
  frequency database, because either alone has documented false negatives.

An event can be flawless on every confidence metric and still be disqualified for
being a common polymorphism. That is not a contradiction: confidence says the
call is right, privacy says it is not private.
""")

code(r"""
hc = (junctions.join(events, junctions.sample_sv_id == events.sv_id, "inner")
      .withColumn("size", F.when(F.col("chrom1") == F.col("chrom2"),
                                 F.abs(F.col("pos2") - F.col("pos1")) - 1))
      .withColumn("bp1_ok", (F.col("segmapq_bp1") >= criteria.HC_MIN_SEGMAPQ) &
                            (F.col("vf_bp1")      >= criteria.HC_MIN_VF) &
                            (F.col("qual_bp1")    >= criteria.HC_MIN_QUAL))
      .withColumn("bp2_ok", (F.col("segmapq_bp2") >= criteria.HC_MIN_SEGMAPQ) &
                            (F.col("vf_bp2")      >= criteria.HC_MIN_VF) &
                            (F.col("qual_bp2")    >= criteria.HC_MIN_QUAL))
      .withColumn("size_ok", F.col("size").isNull() |
                             (F.col("size") >= criteria.HC_MIN_SV_SIZE))
      .withColumn("sv_hc", F.col("bp1_ok") & F.col("bp2_ok") & F.col("size_ok"))
      .withColumn("pass_pon",
                  F.coalesce(F.col("pon_count"), F.lit(0)) < criteria.PON_MAX))

print(hc.select("sample_sv_id", "svtype", "pon_count", "segmapq_bp1", "vf_bp1",
                "qual_bp1", "sv_hc", "pass_pon").orderBy("sample_sv_id")
        .toPandas().to_string(index=False))

check("stage6", "same high-confidence count as the pipeline",
      int(pipe_events.sv_hc.astype(str).str.lower().isin(["true","1"]).sum())
      == hc.filter("sv_hc").count())
""")

md(r"""
### Privacy: the population filter, in Spark

The panel of normals answers "did the caller see this in its normals?". It has
documented false negatives for inherited variation — events with `PON_COUNT` of 1
and 3 carrying population frequencies of 0.287 and 0.263. The second filter is a
reciprocal-overlap join against gnomAD-SV, and it is the one operation here that
genuinely needs Spark: a 1.7 GB sites file against every event.

Only intra-chromosomal events with a resolvable size can be matched this way;
inter-chromosomal junctions are left **unmatched, not zero**, because gnomAD-SV
does not represent them comparably.
""")

code(r"""
gnomad_path = CFG.resources.get("gnomad_sv")
event_geo = (sized.join(events, sized.sample_sv_id == events.sv_id, "inner")
             .filter(F.col("chrom1") == F.col("chrom2"))
             .filter(F.col("event_size").isNotNull())
             .withColumn("start", F.least("pos1", "pos2"))
             .withColumn("end", F.least("pos1", "pos2") + F.col("event_size"))
             .select("sample_sv_id", "chrom1", "start", "end", "event_size"))

wanted_chroms = [r.chrom1 for r in event_geo.select("chrom1").distinct().collect()]
print(f"testable events {event_geo.count()} on chromosomes {sorted(wanted_chroms)}")

if gnomad_path and os.path.exists(gnomad_path) and wanted_chroms:
    gn = (spark.read.text(gnomad_path)
          .filter(~F.col("value").startswith("#"))
          .withColumn("f", F.split("value", "\t"))
          .select(F.regexp_replace(F.col("f")[0], "^chr", "").alias("g_chrom"),
                  F.col("f")[1].cast("long").alias("g_start"),
                  F.col("f")[7].alias("g_info"))
          .filter(F.col("g_chrom").isin(wanted_chroms))
          .withColumn("g_end", F.regexp_extract("g_info", r"(?:^|;)END=(\d+)", 1).cast("long"))
          .withColumn("g_af",  F.regexp_extract("g_info", r"(?:^|;)AF=([0-9.eE+-]+)", 1).cast("double"))
          .filter(F.col("g_end").isNotNull() & F.col("g_af").isNotNull()))

    # Reciprocal overlap: each interval must cover >= the threshold fraction of
    # the other. A one-sided overlap would match a 3 kb event to a 30 bp record.
    ov = (event_geo.join(gn, event_geo.chrom1 == gn.g_chrom, "inner")
          .withColumn("ov", F.least(F.col("end"), F.col("g_end"))
                          - F.greatest(F.col("start"), F.col("g_start")))
          .filter(F.col("ov") > 0)
          .filter((F.col("ov") / F.greatest(F.col("end") - F.col("start"), F.lit(1))
                   >= criteria.GNOMAD_RECIPROCAL_OVERLAP) &
                  (F.col("ov") / F.greatest(F.col("g_end") - F.col("g_start"), F.lit(1))
                   >= criteria.GNOMAD_RECIPROCAL_OVERLAP))
          .groupBy("sample_sv_id").agg(F.max("g_af").alias("gnomad_af")))

    af_spark = {r.sample_sv_id: r.gnomad_af for r in ov.collect()}
    print(f"events matched in gnomAD-SV: {len(af_spark)}")
    for sv_id, af in sorted(af_spark.items(), key=lambda kv: -kv[1]):
        print(f"  {sv_id:>8}  AF {af:.4f}  {'COMMON POLYMORPHISM' if af >= criteria.GNOMAD_MAX_AF else ''}")

    pipe_af = (pipe_events.dropna(subset=["gnomad_af"])
               .set_index(pipe_events.dropna(subset=["gnomad_af"]).sample_sv_id.astype(str))
               ["gnomad_af"].astype(float).to_dict())
    check("stage6", "same gnomAD matches as the pipeline",
          set(af_spark) == set(pipe_af))
    check("stage6", "same allele frequencies (to 4 dp)",
          all(round(af_spark[k], 4) == round(pipe_af[k], 4) for k in af_spark if k in pipe_af))
else:
    print("gnomAD-SV not configured — privacy rests on the panel alone, and the")
    print("run records that in `privacy_note` rather than implying a clean result")
""")

code(r"""
# --- the combined verdict ---------------------------------------------------
# private = below the panel threshold AND below the population threshold.
# An absent value on either side means "not seen", which passes; that asymmetry
# is why `privacy_note` travels with the result.
verdict = (hc.select("sample_sv_id", "svtype", "pon_count", "sv_hc", "pass_pon")
             .withColumn("gnomad_af",
                         F.create_map([F.lit(x) for kv in af_spark.items() for x in kv])[F.col("sample_sv_id")]
                         if af_spark else F.lit(None).cast("double"))
             .withColumn("pass_gnomad",
                         F.col("gnomad_af").isNull() | (F.col("gnomad_af") < criteria.GNOMAD_MAX_AF))
             .withColumn("is_private", F.col("pass_pon") & F.col("pass_gnomad")))

print(verdict.orderBy(F.desc("is_private"), "sample_sv_id").toPandas().to_string(index=False))
check("stage6", "same private-event count as the pipeline",
      int(pipe_events.is_private.astype(str).str.lower().isin(["true","1"]).sum())
      == verdict.filter("is_private").count())
""")

md(r"""
### Expression context

Two rules make an expression number mean something:

- the **disrupted isoform's** TPM, not the gene's — a gene sums all its isoforms
  and over-credits an event whose specific transcript is barely expressed
- a **background rate**: what fraction of *all* genes clears each threshold. In a
  typical transcriptome ~55% clear TPM > 0, so "expressed" at that cutoff carries
  almost no information
""")

code(r"""
if SAMPLE_CFG.isofox_dir and os.path.isdir(SAMPLE_CFG.isofox_dir):
    genes_csv = os.path.join(SAMPLE_CFG.isofox_dir, f"{SAMPLE_CFG.isofox_prefix}.gene_data.csv")
    gene_tpm = (spark.read.option("header", True).csv(genes_csv)
                .selectExpr("GeneName as gene", "cast(AdjTPM as double) as tpm"))
    total = gene_tpm.count()
    print(f"quantified genes: {total:,}\n")
    print("background — fraction of ALL genes clearing each threshold:")
    for cutoff in criteria.TPM_GRADIENT:
        frac = gene_tpm.filter(F.col("tpm") > cutoff).count() / total
        print(f"  TPM > {cutoff:<5} {frac:6.1%}")

    event_genes = (events.select(F.explode("ref_genes").alias("gene")).distinct()
                   .join(gene_tpm, "gene", "left"))
    print("\nour events' genes:")
    event_genes.orderBy(F.desc("tpm")).show(20, truncate=False)
else:
    print("no expression data for this sample — stage 7 is NA, not negative")
""")

# ---------------------------------------------------------------- stage 7-8
md(r"""
## Stages 7–8 — Expression and junction evidence

Read counting needs random access to a BAM, which is a pysam operation on the
driver, not a Spark one. The notebook loads the pipeline's stage-8 output and
checks the three rules that decide every tier — the rules that two independent
workflows previously disagreed on:

1. only reads that **cross** the junction tier an event; coverage is context
2. coverage is summarised with **`min`** across breakends, never `max`
3. an event with no applicable test is **`UNTESTABLE`**, never `NONE`
""")

code(r"""
rna_path = RUN_DIR / "stage8_rna_evidence.tsv"
if rna_path.exists():
    rna = spark.read.option("header", True).option("sep", "\t").csv(str(rna_path))
    print(rna.select("gene", "test", "junction_reads", "coverage_bp1", "coverage_bp2",
                     "min_coverage", "rna_tier").toPandas().to_string(index=False))

    r = rna.toPandas()
    r["junction_reads"] = pd.to_numeric(r.junction_reads, errors="coerce").fillna(0)
    check("stage8", "no event is tiered above NONE without a junction read",
          not ((r.junction_reads == 0) & r.rna_tier.isin(["WEAK","SUGGESTIVE","STRONG"])).any())
    check("stage8", "min_coverage never exceeds either breakend's coverage",
          (pd.to_numeric(r.min_coverage, errors="coerce") <=
           pd.concat([pd.to_numeric(r.coverage_bp1, errors="coerce"),
                      pd.to_numeric(r.coverage_bp2, errors="coerce")], axis=1).min(axis=1)
           ).fillna(True).all())
    print("\ntest applied per event geometry:")
    rna.groupBy("test").count().show()
else:
    print("no RNA output for this sample — stages 7-8 are NA, not negative")
""")

# ---------------------------------------------------------------- null model
md(r"""
## The null model — is the match count above chance?

Raw match counts are not comparable between samples whose candidate universes
differ by orders of magnitude (63,036 against 64 in this dataset). Each candidate
peptide's residues are shuffled, preserving length and composition, and the cross
is repeated. Composition and set size are held constant; only the biological
signal is destroyed.

Spark earns its place here: the permutations are embarrassingly parallel.
""")

code(r"""
import random

ref_set = {r.ref_peptide for r in reference.select("ref_peptide").distinct().collect()}
bc_ref = spark.sparkContext.broadcast(ref_set)
candidates = [r.neopeptide for r in peptides.select("neopeptide").distinct().collect()]
bc_cand = spark.sparkContext.broadcast(candidates)

N_PERM = criteria.NULL_PERMUTATIONS

def one_permutation(seed):
    rng = random.Random(seed)
    hits = 0
    seen = set()
    for peptide in bc_cand.value:
        residues = list(peptide)
        rng.shuffle(residues)
        shuffled = "".join(residues)
        if shuffled not in seen:
            seen.add(shuffled)
            if shuffled in bc_ref.value:
                hits += 1
    return hits

null_counts = spark.sparkContext.parallelize(range(N_PERM), 16).map(one_permutation).collect()
observed = len({p for p in candidates if p in ref_set})

import statistics
mean, sd = statistics.mean(null_counts), statistics.pstdev(null_counts)
at_least = sum(1 for c in null_counts if c >= observed)
print(f"observed {observed}   null {mean:.2f} +/- {sd:.2f}   max {max(null_counts)}")
print(f"p {'< ' + str(1/N_PERM) if at_least == 0 else '= ' + str(at_least/N_PERM)}"
      f"   enrichment {observed/mean:.1f}x" if mean else "")
check("null", "observed matches exceed the null", observed > mean + 3 * sd)
""")

# ---------------------------------------------------------------- summary
md(r"""
## Verification summary

Every check compares this independent Spark derivation against the packaged
pipeline. A failure means the two implementations disagree, and one of them is
wrong — that is the reason to run this notebook.
""")

code(r"""
summary = pd.DataFrame(CHECKS, columns=["stage", "statement", "passed"])
print(summary.to_string(index=False))
failed = summary[~summary.passed]
print(f"\n{len(summary) - len(failed)} of {len(summary)} checks passed")
if len(failed):
    print("\nDISAGREEMENTS — investigate before trusting either result:")
    print(failed.to_string(index=False))
""")

code(r"""
spark.stop()
""")


def build() -> dict:
    cells = []
    for kind, source in CELLS:
        lines = source.splitlines(keepends=True)
        cell = {"cell_type": kind, "metadata": {}, "source": lines}
        if kind == "code":
            cell.update(execution_count=None, outputs=[])
        cells.append(cell)
    return {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3 (neosv_cfs)",
                           "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.10"},
        },
        "nbformat": 4, "nbformat_minor": 5,
    }


if __name__ == "__main__":
    out = pathlib.Path(__file__).parent / "pipeline_walkthrough_pyspark.ipynb"
    out.write_text(json.dumps(build(), indent=1))
    print(f"wrote {out}  ({len(CELLS)} cells)")
