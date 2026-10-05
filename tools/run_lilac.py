#!/usr/bin/env python
"""
run_lilac.py — class I HLA typing of every sample in a run config, with LILAC.

WHY THIS EXISTS
---------------
`presentable` is decided entirely by the genotype it is computed against. A
genotype with no recorded origin makes every presentability call unauditable,
and a genotype that is right for the parent can be wrong for a derived line
that has lost an HLA haplotype. So the genotype is measured per sample, from
that sample's own reads, and the measurement is kept beside the result.

WHY LILAC, AND WHY DNA FIRST
----------------------------
LILAC (hmftools) is what types the reference cohort's patients, so typing the
samples with it keeps the two sides comparable. It types from DNA, where depth
is even, and uses RNA only as supporting evidence per allele: an allele with DNA
support and no RNA support is typed but probably not expressed. With a tumour
(here: derived-line) BAM and a PURPLE directory it also reports allele-specific
copy number, i.e. loss of an HLA haplotype relative to the reference sample.

WHAT IT DOES, PER SAMPLE — replicating nf-core/oncoanalyser 2.0.0
--------------------------------------------------------------------
On an ALT-aware GRCh38 (reads at HLA genes split across chr6 and the HLA/alt
contigs), LILAC needs the reads gathered back onto chr6 first:

  1. slice   samtools view --regions-file <hla slice bed>  (chr6 + alt/HLA contigs)
  2. fastq   samtools sort -n | samtools fastq
  3. realign bwa-mem2 mem -Y against chr6 alone, sort, index
  4. LILAC   root sample:     -reference_bam <own realigned> [-rna_bam]
             derived sample:  -reference_bam <root's realigned>
                              -tumor_bam <own realigned> [-rna_bam] [-purple_dir]

The derived samples are typed against the lineage ROOT, not their immediate
parent, because that is the germline the somatic (PURPLE) calls were made
against. RNA BAMs are used as they are; oncoanalyser does not realign them
either, so RNA reads placed on HLA contigs are not counted (a known loss of RNA
support, not of DNA typing).

Configuration — the run config's `resources.lilac` block:

    resources:
      lilac:
        jar:          /path/to/lilac_v1.6.jar
        resource_dir: /path/to/hmf_pipeline_resources/misc/lilac
        slice_bed:    /path/to/grch38_alt.plus_homologous.bed
        chr6_fasta:   /path/to/chr6.fa      # bwa-mem2 indexed; same build as the BAMs
        java_xmx:     8G                    # optional

and per sample `dna_bam` (required to type), `rna_bam`, `purple_dir` (optional).

Outputs, under --out-dir (required, no default; refuses to overwrite):

    realigned/<sample>.hla.realigned.bam     step 3, reusable with --reuse-realigned
    lilac/<sample>/                          LILAC's own output, untouched
    hla_genotypes.tsv                        one row per (sample, allele), collected
    PROVENANCE.tsv                           tools, versions, inputs, code commit

Usage
-----
    python tools/run_lilac.py --config config/<run>.yaml --out-dir <new dir>
    python tools/run_lilac.py ... --samples <name> --dry-run
"""
from __future__ import annotations

import argparse
import os
import pathlib
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from provenance import code_version                          # noqa: E402
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
from svneo.config import load as load_config                 # noqa: E402

LILAC_GENES = ("A", "B", "C")


def tool(name: str, configured: str | None = None) -> str:
    """Resolve an executable: configured path, else the active environment's
    bin directory, else PATH. Failing here is better than failing mid-pipeline."""
    for candidate in (configured, str(pathlib.Path(sys.executable).parent / name),
                      shutil.which(name)):
        if candidate and os.access(candidate, os.X_OK):
            return candidate
    sys.exit(f"  required tool not found: {name} (set resources.lilac.{name} or "
             f"activate the environment in envs/hla_typing.yml)")


def run(cmd: str, dry: bool, log: pathlib.Path | None = None) -> None:
    print(f"   $ {cmd}", file=sys.stderr)
    if dry:
        return
    with open(log, "a") if log else open(os.devnull, "w") as handle:
        res = subprocess.run(["bash", "-o", "pipefail", "-c", cmd],
                             stdout=handle if log else None, stderr=subprocess.STDOUT if log else None)
    if res.returncode:
        sys.exit(f"  FAILED (exit {res.returncode}): {cmd}" + (f"\n  see {log}" if log else ""))


def realign(sample, bam: str, out: pathlib.Path, cfg: dict, tools: dict,
            threads: int, dry: bool, reuse: bool) -> pathlib.Path:
    final = out / "realigned" / f"{sample}.hla.realigned.bam"
    if reuse and final.exists() and pathlib.Path(str(final) + ".bai").exists():
        print(f"  {sample}: reusing {final}", file=sys.stderr)
        return final
    tmp = out / "realigned" / f"{sample}.tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    log = out / "realigned" / f"{sample}.log"
    st, bw = tools["samtools"], tools["bwa-mem2"]
    q = shlex.quote
    run(f"{st} view -@{threads} -Obam --regions-file {q(cfg['slice_bed'])} {q(bam)} "
        f"| {st} sort -n -@{threads} -T {q(str(tmp / 'n'))} -o {q(str(tmp / 'byname.bam'))}",
        dry, log)
    run(f"{st} fastq -@{threads} {q(str(tmp / 'byname.bam'))} "
        f"-1 {q(str(tmp / 'R1.fq.gz'))} -2 {q(str(tmp / 'R2.fq.gz'))} "
        f"-0 {q(str(tmp / 'other.fq.gz'))} -s {q(str(tmp / 'single.fq.gz'))}", dry, log)
    rg = f"@RG\\tID:{sample}\\tSM:{sample}"
    run(f"{bw} mem -Y -t {threads} -R {q(rg)} {q(cfg['chr6_fasta'])} "
        f"{q(str(tmp / 'R1.fq.gz'))} {q(str(tmp / 'R2.fq.gz'))} "
        f"| {st} sort -@{threads} -T {q(str(tmp / 's'))} -o {q(str(final))} - "
        f"&& {st} index {q(str(final))}", dry, log)
    if not dry:
        shutil.rmtree(tmp)
    return final


def collect(out: pathlib.Path, typed: list[dict]) -> pd.DataFrame:
    """LILAC's per-allele table and QC verdict, one row per (sample, allele)."""
    rows = []
    for t in typed:
        d = out / "lilac" / t["sample"]
        alleles = d / f"{t['sample']}.lilac.tsv"
        qc = d / f"{t['sample']}.lilac.qc.tsv"
        if not alleles.exists():
            continue
        a = pd.read_csv(alleles, sep="\t")
        q = pd.read_csv(qc, sep="\t").iloc[0].to_dict() if qc.exists() else {}
        a.insert(0, "sample", t["sample"])
        a.insert(1, "mode", t["mode"])
        a.insert(2, "reference_sample", t["reference"])
        a["qc_status"] = q.get("Status", "")
        rows.append(a)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--config", required=True)
    p.add_argument("--out-dir", required=True, type=pathlib.Path)
    p.add_argument("--samples", nargs="*", default=None, help="subset (default: all with a DNA BAM)")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--reuse-realigned", action="store_true")
    p.add_argument("--realign-only", action="store_true",
                   help="steps 1-3 only (LILAC resources not needed yet)")
    p.add_argument("--force", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    cfg = load_config(args.config)
    lilac = cfg.resources.get("lilac") or sys.exit("  config has no resources.lilac block")
    needed = ("slice_bed", "chr6_fasta") if args.realign_only else \
        ("jar", "resource_dir", "slice_bed", "chr6_fasta")
    for key in needed:
        if not lilac.get(key) or not os.path.exists(lilac[key]):
            sys.exit(f"  resources.lilac.{key} missing or not found: {lilac.get(key)}")
    tools = {n: tool(n, lilac.get(n)) for n in ("samtools", "bwa-mem2", "java")}

    out = args.out_dir
    if (out / "hla_genotypes.tsv").exists() and not (args.force or args.reuse_realigned):
        sys.exit(f"  REFUSING: {out} already holds hla_genotypes.tsv. Use a new --out-dir "
                 f"(or --reuse-realigned to re-run LILAC on existing realignments).")
    out.mkdir(parents=True, exist_ok=True)

    by = {s.name: s for s in cfg.samples}
    chosen = [s for s in cfg.ordered_samples()
              if s.dna_bam and (args.samples is None or s.name in args.samples)]
    roots = {s.name: ([s.name] + cfg.ancestors(s.name))[-1] for s in cfg.samples}

    realigned: dict[str, pathlib.Path] = {}

    def get_realigned(name: str) -> pathlib.Path:
        if name not in realigned:
            realigned[name] = realign(name, by[name].dna_bam, out, lilac, tools,
                                      args.threads, args.dry_run, args.reuse_realigned)
        return realigned[name]

    if args.realign_only:
        for s in chosen:
            get_realigned(s.name)
        return

    typed = []
    for s in chosen:
        root = roots[s.name]
        q = shlex.quote
        own = get_realigned(s.name)
        if root == s.name or not by[root].dna_bam:
            mode, ref_bam, extra = "germline", own, ""
        else:
            mode, ref_bam = "tumor_vs_root", get_realigned(root)
            extra = f"-tumor_bam {q(str(own))} "
            if s.purple_dir:
                extra += f"-purple_dir {q(s.purple_dir)} "
        if s.rna_bam:
            extra += f"-rna_bam {q(s.rna_bam)} "
        sdir = out / "lilac" / s.name
        sdir.mkdir(parents=True, exist_ok=True)
        print(f"\n── {s.name}  ({mode}{', reference ' + root if mode != 'germline' else ''})",
              file=sys.stderr)
        # Locale pinned: under a comma-decimal locale (e.g. es_ES) LILAC writes
        # "0,00", which every downstream parser reads as text, not a number.
        run(f"{tools['java']} -Xmx{lilac.get('java_xmx', '8G')} "
            f"-Duser.language=en -Duser.country=US -jar {q(lilac['jar'])} "
            f"-sample {q(s.name)} -reference_bam {q(str(ref_bam))} {extra}"
            f"-ref_genome {q(lilac['chr6_fasta'])} -ref_genome_version V38 "
            f"-resource_dir {q(lilac['resource_dir'])} -threads {args.threads} "
            f"-output_dir {q(str(sdir))}", args.dry_run, sdir / "lilac.log")
        typed.append({"sample": s.name, "mode": mode,
                      "reference": s.name if mode == "germline" else root})

    if args.dry_run:
        return
    table = collect(out, typed)
    table.to_csv(out / "hla_genotypes.tsv", sep="\t", index=False)

    version = subprocess.run([tools["java"], "-jar", lilac["jar"], "-version"],
                             capture_output=True, text=True).stdout.strip()
    tv = {n: subprocess.run([tools[n], "--version"], capture_output=True, text=True)
          .stdout.splitlines()[0:1] for n in ("samtools",)}
    prov = [{"key": "generated_utc", "value": datetime.now(timezone.utc).isoformat(timespec="seconds")},
            {"key": "tool", "value": "run_lilac.py"},
            {"key": "lilac", "value": version or lilac["jar"]},
            {"key": "lilac_resources", "value": lilac["resource_dir"]},
            {"key": "slice_bed", "value": lilac["slice_bed"]},
            {"key": "realign_reference", "value": lilac["chr6_fasta"]},
            {"key": "samtools", "value": " ".join(tv["samtools"])},
            {"key": "bwa-mem2", "value": tools["bwa-mem2"]},
            {"key": "method", "value": "nf-core/oncoanalyser 2.0.0 LILAC_CALLING (alt genome): "
                                       "slice -> name-sort -> fastq -> bwa-mem2 -Y on chr6 -> LILAC"}]
    for t in typed:
        s = by[t["sample"]]
        prov.append({"key": f"sample[{s.name}]",
                     "value": f"mode={t['mode']}; reference={t['reference']}; dna={s.dna_bam}; "
                              f"rna={s.rna_bam or '-'}; purple={s.purple_dir or '-'}"})
    pd.DataFrame(prov + code_version()).to_csv(out / "PROVENANCE.tsv", sep="\t", index=False)
    print(f"\nwrote {out / 'hla_genotypes.tsv'} ({len(table)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
