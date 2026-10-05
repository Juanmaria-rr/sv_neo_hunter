#!/usr/bin/env python
"""
run_optitype.py — an independent class I HLA call for every sample, with OptiType.

WHY A SECOND TYPER
------------------
LILAC (tools/run_lilac.py) is the primary call: it is what types the reference
cohort and it reports allele loss. But one tool cannot reveal its own errors.
OptiType reaches its call by a different algorithm (integer linear programming
over reads mapped to an HLA exon reference) and is among the best-benchmarked
class I typers, so agreement between the two is evidence and a disagreement
marks the allele that needs an external typing.

It is run on the SAME reads LILAC sees, so a disagreement is about the
algorithm, not the input:
  DNA  the read pairs from run_lilac.py's realigned/<sample>.hla.realigned.bam
  RNA  the sample's RNA BAM sliced to the same HLA regions (resources.lilac.slice_bed)
DNA decides the genotype; the RNA call says which alleles are expressed well
enough to be typed from transcripts, and is weaker by nature where expression is low.

When the LILAC table exists in --hla-dir, a concordance table is written:
per (sample, gene), each method's two alleles, and whether they agree.

Configuration — the run config's `resources.optitype` block:

    resources:
      optitype:
        bin: /path/to/envs/svneo_optitype/bin   # optitype, razers3, glpsol
        threads: 8                               # optional

Outputs, under --hla-dir (the run_lilac.py output directory):

    optitype/<sample>/{dna,rna}/<sample>_result.tsv   OptiType's own output
    optitype_genotypes.tsv                             all calls, one row per (sample, source)
    hla_concordance.tsv                                LILAC vs OptiType DNA vs OptiType RNA
    OPTITYPE_PROVENANCE.tsv

Usage
-----
    python tools/run_optitype.py --config config/<run>.yaml --hla-dir <run_lilac out dir>
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

GENES = ("A", "B", "C")


def sh(cmd: str, log: pathlib.Path, env: dict) -> None:
    print(f"   $ {cmd}", file=sys.stderr)
    with open(log, "a") as handle:
        res = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], stdout=handle,
                             stderr=subprocess.STDOUT, env=env)
    if res.returncode:
        sys.exit(f"  FAILED (exit {res.returncode}): {cmd}\n  see {log}")


def to_fastq(bam: str, work: pathlib.Path, samtools: str, threads: int,
             log: pathlib.Path, env: dict, regions: str | None = None) -> tuple[str, str]:
    q = shlex.quote
    work.mkdir(parents=True, exist_ok=True)
    view = (f"{samtools} view -@{threads} -Obam --regions-file {q(regions)} {q(bam)} | "
            if regions else f"cat {q(bam)} | ")
    sh(f"{view}{samtools} sort -n -@{threads} -T {q(str(work / 'n'))} -o {q(str(work / 'n.bam'))} -",
       log, env)
    r1, r2 = work / "R1.fq", work / "R2.fq"
    sh(f"{samtools} fastq -@{threads} -1 {q(str(r1))} -2 {q(str(r2))} -0 /dev/null "
       f"-s /dev/null {q(str(work / 'n.bam'))}", log, env)
    return str(r1), str(r2)


def optitype(sample: str, source: str, r1: str, r2: str, out: pathlib.Path,
             bin_dir: str, threads: int, log: pathlib.Path, env: dict) -> pathlib.Path:
    q = shlex.quote
    flag = "--rna" if source == "rna" else "--dna"
    sh(f"{q(os.path.join(bin_dir, 'optitype'))} run -i {q(r1)} -i {q(r2)} {flag} "
       f"-o {q(str(out))} -p {q(sample)} --threads {threads}", log, env)
    return out / f"{sample}_result.tsv"


def two_field(a) -> str:
    """'HLA-A*24:02:01' / 'A*24:02' / 'A*24:02:01:02G' -> 'A*24:02'.

    An expression suffix on the second field is KEPT ('A*24:09N' stays null):
    stripping it would equate a non-expressed allele with an expressed one."""
    s = str(a).replace("HLA-", "").strip()
    if "*" not in s:
        return s
    gene, fields = s.split("*", 1)
    parts = fields.split(":")
    two = ":".join(parts[:2])
    if len(parts) > 2 and parts[-1][-1:].isalpha() and parts[-1][-1] in "NLSCAQ":
        two += parts[-1][-1]
    return f"{gene}*{two}"


def lilac_alleles(hla_dir: pathlib.Path) -> pd.DataFrame:
    path = hla_dir / "hla_genotypes.tsv"
    if not path.exists():
        return pd.DataFrame()
    t = pd.read_csv(path, sep="\t")
    col = next(c for c in t.columns if c.lower() == "allele")
    t["allele2"] = t[col].map(two_field)
    t["gene"] = t["allele2"].str[0]
    return t


def concordance(geno: pd.DataFrame, lilac: pd.DataFrame,
                roots: dict[str, str] | None = None) -> pd.DataFrame:
    rows = []
    for sample in sorted(set(geno["sample"]) | set(lilac.get("sample", []))):
        for gene in GENES:
            rec = {"sample": sample, "gene": f"HLA-{gene}"}
            calls = {}
            if len(lilac):
                l = lilac[(lilac["sample"] == sample) & (lilac["gene"] == gene)]
                calls["lilac"] = sorted(l["allele2"])
            for source in ("dna", "rna"):
                g = geno[(geno["sample"] == sample) & (geno["source"] == source)]
                if len(g):
                    calls[f"optitype_{source}"] = sorted(two_field(g.iloc[0][f"{gene}{i}"])
                                                         for i in (1, 2))
            for k, v in calls.items():
                rec[k] = " / ".join(v)
            primary = calls.get("lilac") or calls.get("optitype_dna")
            rec["dna_concordant"] = (calls.get("lilac") == calls.get("optitype_dna")
                                     if "lilac" in calls and "optitype_dna" in calls else None)
            rec["rna_concordant"] = (primary == calls.get("optitype_rna")
                                     if primary and "optitype_rna" in calls else None)
            rec["genotype"] = " / ".join(primary) if primary else ""
            rec["status"] = ("confirmed" if rec["dna_concordant"] else
                             "DISCORDANT — external typing needed" if rec["dna_concordant"] is False
                             else "single method")
            rows.append(rec)
    out = pd.DataFrame(rows)
    return lineage_check(out, roots or {})


def lineage_check(conc: pd.DataFrame, roots: dict[str, str]) -> pd.DataFrame:
    """A derived sample cannot change a germline allele except by somatic mutation,
    so its genotype should equal the root's. One method disagreeing, while the other
    methods on that sample AND the confirmed root genotype agree, is that method's
    error — not an uncertain genotype. Anything else stays flagged."""
    if not len(conc):
        return conc
    root_geno = {(r.sample, r.gene): r.genotype for r in conc.itertuples()
                 if r.status == "confirmed"}
    conc["root_sample"] = conc["sample"].map(lambda s: roots.get(s, s))
    conc["matches_root"] = [root_geno.get((r.root_sample, r.gene)) == r.genotype
                            if (r.root_sample, r.gene) in root_geno else None
                            for r in conc.itertuples()]
    fix = ((conc["dna_concordant"] == False) & (conc["rna_concordant"] == True)  # noqa: E712
           & (conc["matches_root"] == True) & (conc["root_sample"] != conc["sample"]))  # noqa: E712
    conc.loc[fix, "status"] = ("resolved — one method discordant; the sample's other "
                               "methods and the root's confirmed genotype agree")
    return conc


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--config", required=True)
    p.add_argument("--hla-dir", required=True, type=pathlib.Path)
    p.add_argument("--samples", nargs="*", default=None)
    p.add_argument("--no-rna", action="store_true")
    p.add_argument("--reuse", action="store_true",
                   help="keep existing <sample>_result.tsv files; type only what is missing")
    p.add_argument("--concordance-only", action="store_true",
                   help="rebuild hla_concordance.tsv from existing outputs, e.g. after LILAC finishes")
    args = p.parse_args()

    cfg = load_config(args.config)
    roots = {x.name: ([x.name] + cfg.ancestors(x.name))[-1] for x in cfg.samples}
    if args.concordance_only:
        geno = pd.read_csv(args.hla_dir / "optitype_genotypes.tsv", sep="\t")
        conc = concordance(geno, lilac_alleles(args.hla_dir), roots)
        conc.to_csv(args.hla_dir / "hla_concordance.tsv", sep="\t", index=False)
        print(conc.to_string(index=False), file=sys.stderr)
        return

    ot = cfg.resources.get("optitype") or sys.exit("  config has no resources.optitype block")
    lil = cfg.resources.get("lilac") or {}
    bin_dir = ot.get("bin") or sys.exit("  resources.optitype.bin is required")
    threads = int(ot.get("threads", 8))
    samtools = lil.get("samtools") or shutil.which("samtools") or sys.exit("  samtools not found")
    env = dict(os.environ, PATH=f"{bin_dir}:{os.environ.get('PATH', '')}")

    rows, prov = [], []
    for s in cfg.ordered_samples():
        if args.samples and s.name not in args.samples:
            continue
        base = args.hla_dir / "optitype" / s.name
        sources = []
        realigned = args.hla_dir / "realigned" / f"{s.name}.hla.realigned.bam"
        if realigned.exists():
            sources.append(("dna", str(realigned), None))
        else:
            print(f"  {s.name}: no realigned DNA BAM — run run_lilac.py first", file=sys.stderr)
        if s.rna_bam and not args.no_rna:
            sources.append(("rna", s.rna_bam, lil.get("slice_bed")))
        for source, bam, regions in sources:
            out = base / source
            log = out / "optitype.log"
            result = out / f"{s.name}_result.tsv"
            if args.reuse and result.exists():
                print(f"  {s.name} {source}: reusing {result}", file=sys.stderr)
                r = pd.read_csv(result, sep="\t", index_col=0).iloc[0].to_dict()
                rows.append({"sample": s.name, "source": source, **r})
                prov.append({"key": f"input[{s.name},{source}]", "value": bam})
                continue
            if (out / "work").exists():
                shutil.rmtree(out / "work")      # left by an interrupted run
            out.mkdir(parents=True, exist_ok=True)
            print(f"\n── {s.name} {source}", file=sys.stderr)
            r1, r2 = to_fastq(bam, out / "work", samtools, threads, log, env, regions)
            result = optitype(s.name, source, r1, r2, out, bin_dir, threads, log, env)
            shutil.rmtree(out / "work")
            r = pd.read_csv(result, sep="\t", index_col=0).iloc[0].to_dict()
            rows.append({"sample": s.name, "source": source, **r})
            prov.append({"key": f"input[{s.name},{source}]", "value": bam})

    geno = pd.DataFrame(rows)
    geno.to_csv(args.hla_dir / "optitype_genotypes.tsv", sep="\t", index=False)
    conc = concordance(geno, lilac_alleles(args.hla_dir), roots)
    conc.to_csv(args.hla_dir / "hla_concordance.tsv", sep="\t", index=False)

    version = subprocess.run([os.path.join(bin_dir, "optitype"), "--version"],
                             capture_output=True, text=True, env=env).stdout.strip()
    pd.DataFrame([{"key": "generated_utc",
                   "value": datetime.now(timezone.utc).isoformat(timespec="seconds")},
                  {"key": "tool", "value": "run_optitype.py"},
                  {"key": "optitype", "value": version},
                  {"key": "dna_input", "value": "run_lilac.py realigned BAM (same reads as LILAC)"},
                  {"key": "rna_input", "value": "RNA BAM sliced with resources.lilac.slice_bed"}]
                 + prov + code_version()).to_csv(
        args.hla_dir / "OPTITYPE_PROVENANCE.tsv", sep="\t", index=False)
    print(f"\nwrote {args.hla_dir / 'optitype_genotypes.tsv'} and hla_concordance.tsv",
          file=sys.stderr)
    print(conc.to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
