#!/usr/bin/env python
"""
build_hla_report.py — the HLA typing of a run, explained: an HTML report and a
Markdown file, generated from the typing outputs.

Reads what run_lilac.py and run_optitype.py wrote into --hla-dir, plus each
derived sample's PURPLE fit for context, and writes

    HLA_TYPING_REPORT.html   tables, a figure, a method block on every table
    HLA_TYPING_REPORT.md     the same content as Markdown (no figure; its numbers
                             are in the tables)

Every number comes from those files; every verdict comes from a named rule with
its threshold stated below. Nothing names a sample: names come from the inputs.

HOW EACH QUESTION IS ANSWERED
-----------------------------
genotype      LILAC's call, checked against OptiType on DNA and on RNA
              (hla_concordance.tsv; statuses defined in docs/HLA_TYPING.md).
allele lost?  From DNA FRAGMENTS, not copy number: in each derived sample, the
              share of the gene's DNA fragments that support each allele. An
              allele lost from the sample has a share near 0; below
              LOSS_FRAGMENT_SHARE it is called lost. LILAC's allele copy number
              is shown but not used for the verdict, because it inherits PURPLE's
              purity/ploidy fit, which is not trustworthy for a clonal sample
              (criteria.TRUST_COPY_NUMBER_FIELDS) — the report prints that fit.
expressed?    RNA fragments supporting the allele (LILAC). At least
              EXPRESSED_MIN_RNA_FRAGMENTS, and a share of the gene's RNA
              fragments of at least LOW_EXPRESSION_SHARE, else flagged.
matches what  Optionally (--predictions), the alleles a netMHCpan prediction table
presentable   was computed against are compared with the typed genotype.
used?

Usage
-----
    python tools/build_hla_report.py --config config/<run>.yaml --hla-dir <dir> \\
        [--predictions <netmhcpan table with an 'allele' column>]
"""
from __future__ import annotations

import argparse
import glob
import html
import pathlib
import sys
from datetime import date

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from provenance import code_version                          # noqa: E402
from run_optitype import two_field                           # noqa: E402
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
from svneo import criteria as C                              # noqa: E402
from svneo.config import load as load_config                 # noqa: E402

LOSS_FRAGMENT_SHARE = 0.10        # allele's share of the gene's DNA fragments
EXPRESSED_MIN_RNA_FRAGMENTS = 10
LOW_EXPRESSION_SHARE = 0.15       # allele's share of the gene's RNA fragments
IMBALANCE_DISAGREEMENT = 0.08     # |expected - observed| minor-allele DNA share
GENES = ("A", "B", "C")
COLOURS = [("#2a78d6", "#3987e5"), ("#eb6834", "#d95926"), ("#1baf7a", "#199e70"),
           ("#eda100", "#c98500"), ("#e87ba4", "#d55181"), ("#008300", "#008300")]


# --------------------------------------------------------------------------- data
def load(hla_dir: pathlib.Path, cfg):
    lilac = pd.read_csv(hla_dir / "hla_genotypes.tsv", sep="\t")
    lilac["allele2"] = lilac["Allele"].map(two_field)
    lilac["gene"] = lilac["allele2"].str[0]
    opti = pd.read_csv(hla_dir / "optitype_genotypes.tsv", sep="\t")
    conc = pd.read_csv(hla_dir / "hla_concordance.tsv", sep="\t")
    qc = {}
    for s in lilac["sample"].unique():
        f = hla_dir / "lilac" / s / f"{s}.lilac.qc.tsv"
        if f.exists():
            qc[s] = pd.read_csv(f, sep="\t").iloc[0]
    prov = {}
    for name in ("PROVENANCE.tsv", "OPTITYPE_PROVENANCE.tsv"):
        f = hla_dir / name
        if f.exists():
            prov[name] = pd.read_csv(f, sep="\t")
    purple = {}
    for s in cfg.samples:
        if not s.purple_dir:
            continue
        pur = glob.glob(f"{s.purple_dir}/{s.name}.purple.purity.tsv")
        gene = glob.glob(f"{s.purple_dir}/{s.name}.purple.cnv.gene.tsv")
        rec = {}
        if pur:
            p = pd.read_csv(pur[0], sep="\t").iloc[0]
            rec.update(purity=p.get("purity"), ploidy=p.get("ploidy"))
        if gene:
            g = pd.read_csv(gene[0], sep="\t")
            g = g[g["gene"].isin([f"HLA-{x}" for x in GENES])]
            rec["hla_cn"] = dict(zip(g["gene"], g["minCopyNumber"].round(2)))
        purple[s.name] = rec
    return lilac, opti, conc, qc, prov, purple


def allele_evidence(lilac: pd.DataFrame) -> pd.DataFrame:
    t = lilac.copy()
    derived = t["mode"] != "germline"
    # DNA support in the sample itself: tumour fragments for a derived sample,
    # reference fragments for the root (whose own reads are the reference).
    t["dna_frags"] = t["TumorTotal"].where(derived, t["RefTotal"])
    for col, share in (("dna_frags", "dna_share"), ("RnaTotal", "rna_share")):
        tot = t.groupby(["sample", "gene"])[col].transform("sum")
        t[share] = (t[col] / tot).where(tot > 0)
    t["lost"] = t["dna_share"] < LOSS_FRAGMENT_SHARE
    t["expressed"] = t["RnaTotal"] >= EXPRESSED_MIN_RNA_FRAGMENTS
    t["low_expression"] = t["expressed"] & (t["rna_share"] < LOW_EXPRESSION_SHARE)
    t["somatic"] = t[[c for c in t.columns if c.startswith("Somatic")]].sum(axis=1)
    return t


def final_genotype(conc: pd.DataFrame) -> dict[str, dict[str, str]]:
    out = {}
    for r in conc.itertuples():
        out.setdefault(r.sample, {})[r.gene] = (r.genotype, r.status)
    return out


# --------------------------------------------------------------------------- shared text
def verdicts(ev, conc, order, roots, predicted):
    lines = []
    n = len(conc)
    ok = conc["status"].eq("confirmed").sum()
    res = conc["status"].str.startswith("resolved").sum()
    bad = conc["status"].str.startswith("DISCORDANT").sum()
    lines.append(f"**Genotype.** {ok} of {n} sample–gene calls confirmed by LILAC and OptiType "
                 f"on DNA; {res} resolved by lineage consistency; {bad} discordant.")
    root_geno = conc[conc["sample"].isin(set(roots.values()))]
    same = all(conc[conc["sample"] == s].set_index("gene")["genotype"].to_dict()
               == root_geno[root_geno["sample"] == roots[s]].set_index("gene")["genotype"].to_dict()
               for s in order)
    lines.append("**Lineage.** Every sample carries its root's genotype." if same else
                 "**Lineage.** At least one sample's genotype differs from its root's — see the table.")
    lost = ev[ev["lost"]]
    lines.append("**Loss.** No allele is lost in any sample (lowest share of a gene's DNA "
                 f"fragments: {ev['dna_share'].min():.2f}; loss threshold {LOSS_FRAGMENT_SHARE})."
                 if lost.empty else
                 "**Loss.** Lost: " + ", ".join(f"{r.allele2} in {r.sample}" for r in lost.itertuples()) + ".")
    # Does LILAC's allele copy number (from the PURPLE fit) imply an imbalance that
    # the reads do not show? Expected DNA share of the lower-CN allele = cn_lo/(cn_lo+cn_hi).
    der = ev[ev["mode"] != "germline"]
    if len(der):
        g = der.groupby(["sample", "gene"])
        exp = g["TumorCopyNumber"].min() / g["TumorCopyNumber"].sum()
        obs = g["dna_share"].min()
        gap = (obs - exp).abs()
        worst = gap.idxmax()
        if gap.max() > IMBALANCE_DISAGREEMENT:
            lines.append(f"**Copy number vs reads.** LILAC's allele copy numbers imply a minor-allele "
                         f"DNA share as low as {exp.min():.2f}, but the observed shares are "
                         f"{obs.min():.2f}–{obs.max():.2f} (largest gap: {worst[0]} HLA-{worst[1]}, "
                         f"expected {exp[worst]:.2f}, observed {obs[worst]:.2f}). The reads do not show "
                         "the allelic imbalance the copy numbers suggest: those copy numbers reflect the "
                         "PURPLE fit (section 3), not the sample.")
        else:
            lines.append("**Copy number vs reads.** LILAC's allele copy numbers and the observed DNA "
                         "fragment shares agree.")
    notexp = ev[~ev["expressed"]]
    low = ev[ev["low_expression"]]
    lines.append("**Expression.** Every allele has RNA support in every sample "
                 f"(minimum {int(ev['RnaTotal'].min()):,} fragments)."
                 if notexp.empty else
                 "**Expression.** No RNA support: " + ", ".join(f"{r.allele2} in {r.sample}" for r in notexp.itertuples()) + ".")
    if not low.empty:
        lines.append("**Low relative expression** (share of the gene's RNA fragments below "
                     f"{LOW_EXPRESSION_SHARE}): " + ", ".join(f"{r.allele2} in {r.sample}" for r in low.itertuples()) + ".")
    if predicted is not None:
        typed = {f"HLA-{a}" for a in ev["allele2"]}
        lines.append("**Predictions.** The binding predictions were computed against exactly "
                     "the typed alleles, so no prediction changes." if predicted == typed else
                     f"**Predictions.** The binding predictions used {sorted(predicted)}, the typing "
                     f"gives {sorted(typed)}: predictions must be re-run for "
                     f"{sorted(typed - predicted)} and dropped for {sorted(predicted - typed)}.")
    return lines


METHOD_TEXT = {
    "genotype": ("Per sample and gene, the two alleles at two-field (protein) resolution from three calls: LILAC on DNA (RNA as support), OptiType on DNA, OptiType on RNA.",
                 "LILAC 1.6 on HLA reads realigned to chr6 (nf-core/oncoanalyser 2.0.0 procedure); the root typed on its own reads, derived samples as tumour against the root. OptiType 1.5.0 on the same realigned DNA reads, and on the RNA BAM sliced to the same regions. Status rules are in docs/HLA_TYPING.md.",
                 "In LILAC's tumour mode the genotype is solved on the ROOT's reads; the derived sample contributes allele support, copy number and somatic variants. The per-sample independent DNA typing is OptiType's. Two-field resolution: synonymous and non-coding differences are not resolved."),
    "evidence": ("Per allele and sample: DNA fragments supporting the allele in that sample, their share of the gene's fragments, RNA fragments and share, LILAC allele copy number, somatic variants in the allele.",
                 f"Fragments from LILAC (`RefTotal` for the root, `TumorTotal` for derived samples; `RnaTotal`). Lost: DNA share < {LOSS_FRAGMENT_SHARE}. Expressed: ≥ {EXPRESSED_MIN_RNA_FRAGMENTS} RNA fragments; low relative expression: RNA share < {LOW_EXPRESSION_SHARE}.",
                 "Shares are within a gene: a share near 0.5 is what two equally represented alleles give, and departures reflect both biology and the alleles' differing numbers of distinguishing positions (unique vs shared fragments). Copy number is printed but not used: it inherits the PURPLE fit (see next table). RNA reads the aligner placed on HLA contigs are not counted, so RNA support is an underestimate."),
    "purple": ("PURPLE's fitted purity and ploidy for each derived sample, and its total copy number at the three class I genes.",
               "Read from <code>&lt;sample&gt;.purple.purity.tsv</code> and <code>.purple.cnv.gene.tsv</code> in the sample's configured PURPLE directory.",
               "A clonal cell line typed against its own parent is ~100% 'tumour'; a fitted purity far below 1 means the fit is not describing the sample, and the allele copy numbers LILAC derives from it should not be read as absolute counts. This pipeline does not trust copy-number fields for clonal samples (TRUST_COPY_NUMBER_FIELDS = False)."),
    "qc": ("LILAC's own quality record per sample.",
           "Copied from <code>&lt;sample&gt;.lilac.qc.tsv</code>. Score margin: how far the chosen solution beats the next; next solution: the closest alternative genotype.",
           "WARN_UNMATCHED_AMINO_ACID means some fragments carry amino-acid sequences matching no candidate allele; it is common and, with a large score margin and an alternative differing in one rare allele, does not by itself question the call."),
    "optitype": ("OptiType's six alleles per sample and source, with the number of reads it used.",
                 "OptiType 1.5.0 (RazerS3 mapping to its HLA exon reference, integer linear programming).",
                 "OptiType on DNA works from hundreds of reads here, so a call between near-identical alleles can flip; on RNA it uses tens of thousands, but only for expressed alleles."),
}

LIMITS = [
    "Class I only (A, B, C). Class II is not typed.",
    "Two-field resolution; null alleles are kept distinct (`N` suffix), but a non-coding variant that silences an allele would not be seen in the genotype — RNA support is the check for that.",
    "LILAC's allele copy number depends on PURPLE's purity/ploidy fit; for clonal samples that fit is not reliable, so loss is judged from fragment shares instead.",
    "RNA BAMs are used as aligned; RNA reads on HLA/alt contigs are not counted.",
    "No external (clinical) typing was available to compare against; agreement is between two algorithms on the same reads, not with a gold standard.",
]


# --------------------------------------------------------------------------- html
def e(x) -> str:
    return html.escape("" if x is None else str(x))


def tbl(headers, rows) -> str:
    th = "".join(f"<th>{h}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f"<div class=tw><table><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table></div>"


def method_html(key) -> str:
    a, b, c = METHOD_TEXT[key]
    return (f"<div class=method><div><h4>What is counted</h4><p>{a}</p></div>"
            f"<div><h4>How it is computed</h4><p>{b}</p></div>"
            f"<div><h4>What it does not capture</h4><p>{c}</p></div></div>")


def md_inline(s: str) -> str:
    import re
    s = e(s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", s)
    return re.sub(r"`([^`]+)`", r"<code>\1</code>", s)


def figure(ev: pd.DataFrame, order: list[str], colour: dict) -> str:
    """Share of each gene's DNA fragments per allele, per sample (grouped bars)."""
    alleles = sorted(ev["allele2"].unique())
    W, H, left, top, bottom = 760, 300, 48, 20, 60
    gw = (W - left - 10) / len(alleles)
    bw = min(16, (gw - 14) / len(order))
    y = lambda v: top + (1 - v) * (H - top - bottom)              # noqa: E731
    parts = [f"<svg viewBox='0 0 {W} {H}' role=img aria-label='DNA fragment share per allele' style='width:100%;height:auto'>"]
    for v in (0, 0.25, 0.5, 0.75, 1.0):
        parts.append(f"<line x1='{left}' x2='{W-10}' y1='{y(v):.1f}' y2='{y(v):.1f}' style='stroke:var(--border);stroke-width:1'></line>"
                     f"<text x={left-6} y={y(v)+4:.1f} text-anchor=end class=ax>{v:.2f}</text>")
    parts.append(f"<line x1='{left}' x2='{W-10}' y1='{y(LOSS_FRAGMENT_SHARE):.1f}' y2='{y(LOSS_FRAGMENT_SHARE):.1f}' "
                 f"style='stroke:var(--stop-bd);stroke-width:1.5;stroke-dasharray:4 3'></line>"
                 f"<text x={left+6} y={y(LOSS_FRAGMENT_SHARE)-4:.1f} class=ax>loss threshold ({LOSS_FRAGMENT_SHARE})</text>")
    for i, a in enumerate(alleles):
        x0 = left + i * gw + (gw - bw * len(order) - 2 * (len(order) - 1)) / 2
        for j, s in enumerate(order):
            r = ev[(ev["allele2"] == a) & (ev["sample"] == s)]
            if r.empty or pd.isna(r["dna_share"].iloc[0]):
                continue
            v = float(r["dna_share"].iloc[0])
            x = x0 + j * (bw + 2)
            parts.append(f"<rect x='{x:.1f}' y='{y(v):.1f}' width='{bw:.1f}' height='{y(0)-y(v):.1f}' rx='3' "
                         f"style='fill:var({colour[s]})'><title>{e(s)} {e(a)}: {v:.2f} of the gene's DNA fragments "
                         f"({int(r['dna_frags'].iloc[0])} fragments)</title></rect>")
        parts.append(f"<text x={left + i*gw + gw/2:.1f} y={H-bottom+18} text-anchor=middle class=lb>{e(a)}</text>")
    parts.append("</svg>")
    legend = " ".join(f"<span class=chip><i style='background:var({colour[s]})'></i>{e(s)}</span>" for s in order)
    return f"<div class=fig>{''.join(parts)}<div class=legend>{legend}</div></div>"


CSS = """
:root{color-scheme:light;--s0:#f4f4f1;--s1:#fcfcfb;--s2:#eceae4;--border:#dcd9d0;--t1:#0b0b0b;--t2:#52514e;--t3:#76746d;
--stop-bd:#e34948;--ok:#1baf7a;--warn-bg:#fdf3e7;--warn-bd:#eda100;__L__}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;--s0:#121211;--s1:#1a1a19;--s2:#232321;
--border:#35342f;--t1:#fff;--t2:#c3c2b7;--t3:#94928a;--stop-bd:#e66767;--ok:#199e70;--warn-bg:#2b2415;--warn-bd:#c98500;__D__}}
:root[data-theme="dark"]{color-scheme:dark;--s0:#121211;--s1:#1a1a19;--s2:#232321;--border:#35342f;--t1:#fff;--t2:#c3c2b7;
--t3:#94928a;--stop-bd:#e66767;--ok:#199e70;--warn-bg:#2b2415;--warn-bd:#c98500;__D__}
*{box-sizing:border-box}body{margin:0;background:var(--s0);color:var(--t1);font:15px/1.6 ui-sans-serif,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
main{max-width:1040px;margin:0 auto;padding:32px 16px 80px}h1{font-size:27px;line-height:1.2;margin:0 0 6px}
h2{font-size:20px;margin:44px 0 8px;padding-top:14px;border-top:1px solid var(--border)}h3{font-size:16px;margin:24px 0 6px}
h4{font-size:12px;text-transform:uppercase;letter-spacing:.04em;color:var(--t3);margin:0 0 4px}p,li{color:var(--t2)}b{color:var(--t1)}
code{font:12.5px ui-monospace,Menlo,monospace;background:var(--s2);padding:1px 4px;border-radius:4px;color:var(--t1)}
.meta{color:var(--t3);font-size:13px}.tw{overflow-x:auto;margin:10px 0}
table{border-collapse:collapse;width:100%;font-size:13.5px;background:var(--s1)}th,td{text-align:left;padding:6px 9px;border-bottom:1px solid var(--border);vertical-align:top}
th{font-size:12px;color:var(--t3);background:var(--s2)}.method{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:14px;font-size:13px;background:var(--s2);border-radius:8px;padding:12px 14px;margin:12px 0}
.method p{margin:0}.summary{background:var(--s1);border:1px solid var(--border);border-radius:10px;padding:8px 18px}
.ok{color:var(--ok);font-weight:600}.bad{color:var(--stop-bd);font-weight:600}.warn{background:var(--warn-bg);border-left:3px solid var(--warn-bd);padding:8px 12px;border-radius:6px}
.fig{background:var(--s1);border:1px solid var(--border);border-radius:10px;padding:12px}.ax{font-size:11px;fill:var(--t3)}.lb{font-size:12px;fill:var(--t2)}
.legend{display:flex;flex-wrap:wrap;gap:14px;font-size:13px;color:var(--t2);margin-top:6px}.chip{display:inline-flex;align-items:center;gap:6px}
.chip i{width:10px;height:10px;border-radius:3px;display:inline-block}
"""


def status_cell(status: str) -> str:
    if status == "confirmed":
        return "<span class=ok>confirmed</span>"
    if status.startswith("resolved"):
        return f"<span class=ok>resolved</span><br><small>{e(status.split('—', 1)[-1].strip())}</small>"
    return f"<span class=bad>{e(status)}</span>"


# --------------------------------------------------------------------------- build
def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--config", required=True)
    p.add_argument("--hla-dir", required=True, type=pathlib.Path)
    p.add_argument("--predictions", default=None,
                   help="optional netMHCpan table (column 'allele') to check against the typing")
    p.add_argument("--out-html", type=pathlib.Path, default=None)
    p.add_argument("--out-md", type=pathlib.Path, default=None)
    args = p.parse_args()
    out_html = args.out_html or args.hla_dir / "HLA_TYPING_REPORT.html"
    out_md = args.out_md or args.hla_dir / "HLA_TYPING_REPORT.md"

    cfg = load_config(args.config)
    lilac, opti, conc, qc, prov, purple = load(args.hla_dir, cfg)
    order = [s.name for s in cfg.ordered_samples() if s.name in set(lilac["sample"])]
    roots = {s.name: ([s.name] + cfg.ancestors(s.name))[-1] for s in cfg.samples}
    parents = {s.name: s.parent for s in cfg.samples}
    by = {s.name: s for s in cfg.samples}
    ev = allele_evidence(lilac)
    geno = final_genotype(conc)
    predicted = None
    if args.predictions:
        alle = pd.read_csv(args.predictions, sep="\t", usecols=["allele"])["allele"].unique()
        predicted = {f"HLA-{two_field(a)}" for a in alle}
    summary = verdicts(ev, conc, order, roots, predicted)
    stamp = code_version()
    commit = next((r["value"] for r in stamp if r["key"] == "code_commit"), "?")
    clean = next((r["value"] for r in stamp if r["key"] == "code_tree_clean"), "?")
    lv = prov.get("PROVENANCE.tsv", pd.DataFrame(columns=["key", "value"])).set_index("key")["value"]
    ov = prov.get("OPTITYPE_PROVENANCE.tsv", pd.DataFrame(columns=["key", "value"])).set_index("key")["value"]

    # ---------- table rows shared by both outputs
    samples_rows = [[s, parents.get(s) or "— (root)", "germline" if roots[s] == s else f"tumour vs {roots[s]}",
                     "yes" if by[s].dna_bam else "no", "yes" if by[s].rna_bam else "no",
                     "yes" if by[s].purple_dir else "no"] for s in order]
    geno_rows = []
    for s in order:
        for g in GENES:
            r = conc[(conc["sample"] == s) & (conc["gene"] == f"HLA-{g}")]
            if r.empty:
                continue
            r = r.iloc[0]
            geno_rows.append([s, f"HLA-{g}", r.get("lilac", ""), r.get("optitype_dna", ""),
                              r.get("optitype_rna", ""), r["genotype"], r["status"]])
    ev_rows = []
    for s in order:
        for r in ev[ev["sample"] == s].sort_values("allele2").itertuples():
            ev_rows.append([s, r.allele2, int(r.dna_frags), f"{r.dna_share:.2f}", int(r.RnaTotal),
                            f"{r.rna_share:.2f}" if r.rna_share == r.rna_share else "—",
                            f"{r.TumorCopyNumber:.2f}" if r.mode != "germline" else "—",
                            int(r.somatic),
                            "LOST" if r.lost else ("expressed" if r.expressed else "no RNA")
                            + (" (low share)" if r.low_expression else "")])
    purple_rows = [[s, f"{purple[s].get('purity', float('nan')):.2f}", f"{purple[s].get('ploidy', float('nan')):.2f}",
                    ", ".join(f"{k} {v}" for k, v in purple[s].get("hla_cn", {}).items())]
                   for s in order if s in purple]
    qc_rows = [[s, q["Status"], q["ScoreMargin"], q["NextSolutionAlleles"], int(q["TotalFragments"]),
                int(q["FittedFragments"]), int(q["UnmatchedFragments"]),
                f"{int(q['A_LowCoverageBases'])}/{int(q['B_LowCoverageBases'])}/{int(q['C_LowCoverageBases'])}",
                int(q["SomaticVariantsMatched"])] for s, q in qc.items()]
    ot_rows = [[r.sample, r.source.upper(), " ".join(two_field(getattr(r, f"{g}{i}")) for g in GENES for i in (1, 2)),
                f"{int(r.Reads):,}"] for r in opti.itertuples()]
    final = conc[conc["sample"] == order[0]].set_index("gene")["genotype"].to_dict() if order else {}

    tools = [("LILAC", lv.get("lilac", "?")), ("LILAC resources", "hmf_pipeline_resources.38 v2.0.0 (misc/lilac)"),
             ("Realignment", "samtools 1.21 + bwa-mem2 2.2.1 -Y onto chr6 of GRCh38_masked_exclusions_alts_hlas (HMF 25.1)"),
             ("OptiType", ov.get("optitype", "?")), ("Procedure", lv.get("method", "?")),
             ("Code", f"sv_neo_hunter {commit[:7]} (clean tree: {clean})")]

    # ---------- HTML
    colour = {s: f"--c{i}" for i, s in enumerate(order)}
    light = "".join(f"--c{i}:{COLOURS[i % 6][0]};" for i in range(len(order)))
    dark = "".join(f"--c{i}:{COLOURS[i % 6][1]};" for i in range(len(order)))
    H = [f"<h1>HLA class I typing</h1><p class=meta>{e(' · '.join(order))} · generated {date.today()} · "
         f"code {e(commit[:7])} · contains results, not for publication</p>",
         "<div class=summary><h3>Summary</h3>",
         tbl(["gene", "genotype"], [[g, f"<b>{e(v)}</b>"] for g, v in final.items()]),
         "<ul>" + "".join(f"<li>{md_inline(x)}</li>" for x in summary) + "</ul></div>",
         "<h2>What was analysed</h2><p>Each sample's own sequencing: whole-genome DNA (typing), RNA "
         "(expression of each allele, and an independent call) and, for derived samples, the PURPLE "
         "copy-number fit against the lineage root.</p>",
         tbl(["sample", "parent", "LILAC mode", "DNA", "RNA", "PURPLE"], [[e(c) for c in r] for r in samples_rows]),
         "<h2>Method</h2><ol>"
         "<li><b>Gather HLA reads onto chr6.</b> The BAMs were aligned to an alt-aware GRCh38 in which reads "
         "from the HLA genes are split between chr6 and ~500 HLA/alt contigs. Reads in the HLA regions and on "
         "those contigs are sliced out, converted back to read pairs and realigned to chr6 alone.</li>"
         "<li><b>Type with LILAC</b> (primary). The root sample on its own reads; each derived sample as "
         "'tumour' against the root, which adds per-allele support in the derived sample, allele copy number "
         "from PURPLE, and somatic variants inside HLA. RNA is passed as support.</li>"
         "<li><b>Check with OptiType</b> on the same DNA reads and, separately, on the sample's RNA.</li>"
         "<li><b>Compare.</b> Two-field alleles per gene; agreement of LILAC and OptiType-DNA confirms; a lone "
         "discordant method is resolved only if the sample's other methods and the root's confirmed genotype "
         "agree.</li>"
         "<li><b>Assess loss and expression</b> from fragment counts per allele.</li></ol>",
         tbl(["component", "version / setting"], [[e(a), e(b)] for a, b in tools]),
         "<h2>Results</h2><h3>1 · Genotype per sample and method</h3>",
         tbl(["sample", "gene", "LILAC", "OptiType DNA", "OptiType RNA", "genotype", "status"],
             [[e(r[0]), e(r[1]), e(r[2]), e(r[3]), e(r[4]), f"<b>{e(r[5])}</b>", status_cell(r[6])] for r in geno_rows]),
         method_html("genotype"),
         "<h3>2 · Is any allele lost? Is each expressed?</h3>",
         figure(ev, order, colour),
         "<p class=meta>Figure: for each allele, its share of the gene's DNA fragments in each sample. Two "
         "equally represented alleles give ~0.5; a lost allele falls to ~0 (dashed line: loss threshold). "
         "Hover a bar for the fragment count.</p>",
         tbl(["sample", "allele", "DNA fragments", "DNA share", "RNA fragments", "RNA share", "allele CN (LILAC)",
              "somatic variants", "verdict"], [[e(c) for c in r] for r in ev_rows]),
         method_html("evidence"),
         "<h3>3 · Copy-number context: the PURPLE fit</h3>",
         tbl(["sample", "fitted purity", "fitted ploidy", "total copy number at HLA genes"],
             [[e(c) for c in r] for r in purple_rows]) if purple_rows else "<p>No PURPLE directory configured.</p>",
         method_html("purple"),
         "<h3>4 · LILAC quality</h3>",
         tbl(["sample", "status", "score margin", "next solution", "fragments", "fitted", "unmatched",
              "low-coverage bases A/B/C", "somatic variants"], [[e(c) for c in r] for r in qc_rows]),
         method_html("qc"),
         "<h3>5 · OptiType calls</h3>",
         tbl(["sample", "source", "alleles", "reads used"], [[e(c) for c in r] for r in ot_rows]),
         method_html("optitype"),
         "<h2>Limitations</h2><ul>" + "".join(f"<li>{md_inline(x)}</li>" for x in LIMITS) + "</ul>",
         "<h2>Files</h2>" + tbl(["file", "content"], [
             ["<code>hla_concordance.tsv</code>", "the genotype table (section 1)"],
             ["<code>hla_genotypes.tsv</code>", "LILAC per allele and sample (section 2)"],
             ["<code>optitype_genotypes.tsv</code>", "OptiType calls (section 5)"],
             ["<code>lilac/&lt;sample&gt;/</code>", "LILAC's own output and QC (section 4)"],
             ["<code>PROVENANCE.tsv</code>, <code>OPTITYPE_PROVENANCE.tsv</code>", "inputs, versions, code commit"]]),
         "<h2>Provenance</h2>" + tbl(["key", "value"], [[e(r["key"]), e(r["value"])] for r in stamp])]
    css = CSS.replace("__L__", light).replace("__D__", dark)
    out_html.write_text("<!doctype html><html lang=en><head><meta charset=utf-8>"
                        "<meta name=viewport content='width=device-width,initial-scale=1'>"
                        f"<title>HLA Typing Report</title><style>{css}</style></head><body><main>"
                        + "".join(H) + "</main></body></html>")

    # ---------- Markdown
    def mdt(headers, rows) -> str:
        esc = lambda x: str(x).replace("|", "\\|").replace("*", "\\*")   # noqa: E731
        return ("| " + " | ".join(headers) + " |\n|" + "---|" * len(headers) + "\n"
                + "".join("| " + " | ".join(esc(c) for c in r) + " |\n" for r in rows))

    def mdm(key) -> str:
        import re
        a, b, c = (re.sub(r"<[^>]+>", "", x).replace("&lt;", "<").replace("&gt;", ">")
                   for x in METHOD_TEXT[key])
        return f"\n> **Counted.** {a}\n>\n> **Computed.** {b}\n>\n> **Not captured.** {c}\n"

    M = [f"# HLA class I typing\n\n{' · '.join(order)} · generated {date.today()} · code `{commit[:7]}` · "
         "contains results, not for publication\n\n## Summary\n\n",
         mdt(["gene", "genotype"], [[g, v] for g, v in final.items()]), "\n",
         "".join(f"- {x}\n" for x in summary),
         "\n## What was analysed\n\nEach sample's own sequencing: whole-genome DNA (typing), RNA (expression "
         "of each allele, and an independent call) and, for derived samples, the PURPLE copy-number fit "
         "against the lineage root.\n\n", mdt(["sample", "parent", "LILAC mode", "DNA", "RNA", "PURPLE"], samples_rows),
         "\n## Method\n\n"
         "1. **Gather HLA reads onto chr6.** The BAMs were aligned to an alt-aware GRCh38 in which HLA reads are "
         "split between chr6 and ~500 HLA/alt contigs; reads in the HLA regions and on those contigs are sliced "
         "out, converted back to pairs and realigned to chr6 alone.\n"
         "2. **Type with LILAC** (primary): the root on its own reads; each derived sample as 'tumour' against "
         "the root (allele support in the derived sample, allele copy number from PURPLE, somatic HLA "
         "variants). RNA as support.\n"
         "3. **Check with OptiType** on the same DNA reads and, separately, on the sample's RNA.\n"
         "4. **Compare** two-field alleles per gene: LILAC = OptiType-DNA confirms; a lone discordant method is "
         "resolved only if the sample's other methods and the root's confirmed genotype agree.\n"
         "5. **Assess loss and expression** from fragment counts per allele.\n\n",
         mdt(["component", "version / setting"], tools),
         "\n## Results\n\n### 1 · Genotype per sample and method\n\n",
         mdt(["sample", "gene", "LILAC", "OptiType DNA", "OptiType RNA", "genotype", "status"], geno_rows), mdm("genotype"),
         "\n### 2 · Is any allele lost? Is each expressed?\n\n",
         mdt(["sample", "allele", "DNA fragments", "DNA share", "RNA fragments", "RNA share", "allele CN (LILAC)",
              "somatic variants", "verdict"], ev_rows), mdm("evidence"),
         "\n### 3 · Copy-number context: the PURPLE fit\n\n",
         mdt(["sample", "fitted purity", "fitted ploidy", "total copy number at HLA genes"], purple_rows), mdm("purple"),
         "\n### 4 · LILAC quality\n\n",
         mdt(["sample", "status", "score margin", "next solution", "fragments", "fitted", "unmatched",
              "low-coverage bases A/B/C", "somatic variants"], qc_rows), mdm("qc"),
         "\n### 5 · OptiType calls\n\n", mdt(["sample", "source", "alleles", "reads used"], ot_rows), mdm("optitype"),
         "\n## Limitations\n\n", "".join(f"- {x}\n" for x in LIMITS),
         "\n## Files\n\n", mdt(["file", "content"], [
             ["hla_concordance.tsv", "the genotype table (section 1)"],
             ["hla_genotypes.tsv", "LILAC per allele and sample (section 2)"],
             ["optitype_genotypes.tsv", "OptiType calls (section 5)"],
             ["lilac/<sample>/", "LILAC's own output and QC (section 4)"],
             ["PROVENANCE.tsv, OPTITYPE_PROVENANCE.tsv", "inputs, versions, code commit"]]),
         f"\n## Provenance\n\nGenerated by `tools/build_hla_report.py`, sv_neo_hunter `{commit[:7]}` "
         f"(clean tree: {clean}). Regenerate rather than edit.\n"]
    out_md.write_text("".join(M))
    print(f"wrote {out_html}\nwrote {out_md}", file=sys.stderr)


if __name__ == "__main__":
    main()
