#!/usr/bin/env python
"""
build_lineage_report.py — one report for every line of a lineage, generated.

WHY GENERATED
-------------------
The first results report covered two lines and was written by hand. Extending it
to more lines would mean re-deriving each figure by hand and keeping prose and
tables in step with a universe that keeps being re-annotated. A hand-edited page
has already been half-destroyed by one bad string replacement and had no source
to rebuild it from. So this builds the report from the run's own outputs: every
number in it is computed here from the files it names, and refining the report
means editing this file and re-running it.

WHAT IT READS
-------------
  --config         the run config: the lines, their lineage and their DNA BAMs
  --universe-dir   candidate_universe.tsv, its column dictionary and glossary
  --cross-dir      per-sample summary.json and funnel.tsv (stage-1 admission)
  --judgements     OPTIONAL. Per-event verdicts made by a person (e.g. "not a
                   lesion", with the reason). No column of the universe
                   encodes them, so they are kept in a separate file that
                   names who decided what and why, and the report always shows
                   the pipeline's count beside the count after them.
                   Columns: acquired_in, sv_id, judgement, reason.

WHAT IT MEASURES ITSELF
-----------------------
DNA depth across the interval of every surviving and near-miss event, the
discriminator between a genomic deletion and a splicing event (see
validate_junction.py). It is written to <out>.dna_depth.tsv so the figure in
the report can be checked against the table. --no-dna-depth skips it.

Nothing here names a line, gene or cohort: every name comes from the inputs.

Usage
-----
    python tools/build_lineage_report.py --config config/<run>.yaml \\
        --universe-dir <results>/universe --cross-dir <cross> \\
        --judgements <private>/event_judgements.tsv --out <dir>/REPORT.html
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import pathlib
import re
import sys
from datetime import date

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from provenance import code_version                          # noqa: E402
from validate_junction import copy_number, verdict, FLANK_BP, EDGE_TRIM  # noqa: E402
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
from svneo import criteria as C                              # noqa: E402
from svneo.config import load as load_config                 # noqa: E402

# Categorical slots in fixed order (validated palette), light / dark. A line
# keeps its slot by lineage position, never by rank in any table.
LINE_COLOURS = [("#2a78d6", "#3987e5"), ("#eb6834", "#d95926"),
                ("#1baf7a", "#199e70"), ("#eda100", "#c98500"),
                ("#e87ba4", "#d55181"), ("#008300", "#008300"),
                ("#4a3aa7", "#9085e9"), ("#e34948", "#e66767")]

RNA_CROSSING = ["STRONG", "SUGGESTIVE"]      # the credible_and_presentable recipe
CRITERIA = {                                  # name -> (label, column test)
    "hc":   ("high-confidence SV call", lambda t: flag(t, "sv_hc")),
    "expr": ("both partner genes transcribed", lambda t: flag(t, "expressed")),
    "rna":  ("RNA reads cross the junction",
             lambda t: t["rna_tier"].isin(RNA_CROSSING)),
    "pres": ("presentable by the line's own HLA (netMHCpan)", lambda t: flag(t, "presentable")),
}
DEPTH_MAX_INTERVAL = 2_000_000   # beyond this the depth scan is slow and the
                                 # "interval" is a chromosome arm, not a lesion
DEPTH_TYPES = {"DEL": "≈0.5 heterozygous, ≈0 homozygous, ≈1 no deletion",
               "DUP": "≈1.5 one extra copy, ≈1 no duplication"}


# --------------------------------------------------------------------------- data
def flag(t: pd.DataFrame, col: str) -> pd.Series:
    if col not in t:
        return pd.Series(False, index=t.index)
    return t[col].astype(str).str.lower().isin(["true", "1", "yes"])


def acquired_by(t: pd.DataFrame, line: str) -> pd.Series:
    return t["acquired_in"].astype(str).str.split(";").apply(lambda x: line in x)


def carried_by(t: pd.DataFrame, line: str) -> pd.Series:
    return t["present_in"].astype(str).str.split(";").apply(lambda x: line in x)


def sha16(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def qc_mask(t):
    return ~flag(t, "is_self") & ~flag(t, "low_complexity")


def conj_without_pres(t: pd.DataFrame) -> pd.Series:
    """QC-clean and every criterion except presentability."""
    m = qc_mask(t)
    for key, (_, test) in CRITERIA.items():
        if key != "pres":
            m &= test(t)
    return m


def funnel(t: pd.DataFrame, lines: list[str]) -> pd.DataFrame:
    qc = qc_mask(t)
    conj = qc.copy()
    for _, test in CRITERIA.values():
        conj &= test(t)
    rows = []
    for line in lines:
        a = acquired_by(t, line)
        row = {"line": line, "candidates": int(a.sum()),
               "non_self": int((a & ~flag(t, "is_self")).sum()),
               "qc": int((a & qc).sum())}
        for key, (_, test) in CRITERIA.items():
            row[key] = int((a & qc & test(t)).sum())
        row["conjunction"] = int((a & conj).sum())
        if "presentable_mhcflurry" in t:
            mf = flag(t, "presentable_mhcflurry")
            base = conj_without_pres(t)
            row["pres_mf"] = int((a & qc & mf).sum())
            row["pres_any"] = int((a & qc & (mf | flag(t, "presentable"))).sum())
            row["pres_both"] = int((a & qc & mf & flag(t, "presentable")).sum())
            row["conj_mf"] = int((a & base & mf).sum())
            row["conj_both"] = int((a & conj & mf).sum())
            row["conj_any"] = int((a & (conj | (base & mf))).sum())
        row["matched_reference"] = int((a & flag(t, "matched_reference")).sum())
        row["carried"] = int(carried_by(t, line).sum())
        rows.append(row)
    return pd.DataFrame(rows).set_index("line")


def stage1(cross: pathlib.Path, lines: list[str], branch: str) -> pd.DataFrame:
    rows = []
    for line in lines:
        d = cross / f"{line}_{branch}"
        s = json.loads((d / "summary.json").read_text()).get("stage1_funnel", {})
        f = pd.read_csv(d / "funnel.tsv", sep="\t").set_index("step")["n"]
        rows.append({"line": line, "records": s.get("records"),
                     "after_filter": s.get("after_filter_pass"),
                     "pon_admitted": s.get("caller_pon_admitted"),
                     "paired": s.get("after_paired_only"),
                     "junctions": f.get("junctions"),
                     "peptides": f.get("candidate_peptides")})
    return pd.DataFrame(rows).set_index("line")


def events(t: pd.DataFrame, mask: pd.Series) -> pd.DataFrame:
    """Peptide rows -> one row per (acquiring line, junction)."""
    sub = t[mask].copy()
    sub["line"] = sub["acquired_in"].astype(str)
    agg = {
        "genes": ("gene1", lambda s: "/".join(sorted(set(s.dropna().astype(str))))),
        "gene2": ("gene2", lambda s: "/".join(sorted(set(s.dropna().astype(str))))),
        "chrom1": ("chrom1", "first"), "pos1": ("pos1", "first"),
        "chrom2": ("chrom2", "first"), "pos2": ("pos2", "first"),
        "svtype": ("svtype", "first"), "event_size": ("event_size", "first"),
        "peptides": ("neopeptide", "nunique"),
        "junction_reads": ("junction_reads", "max"),
        "junction_usage": ("junction_usage", "max"),
        "junction_interval_total": ("junction_interval_total", "max"),
        "test": ("test", "first"),
        "vf_bp1": ("vf_bp1", "max"), "ref_bp1": ("ref_bp1", "max"),
        "dna_vaf_bp1": ("dna_vaf_bp1", "max"), "dna_vaf_bp2": ("dna_vaf_bp2", "max"),
        "pon_count": ("pon_count", "max"),
        "gnomad_af_popmax": ("gnomad_af_popmax", "max"),
        "nearest_alt_sj_bp": ("nearest_alt_sj_bp", "min"),
        "alt_sj_type": ("alt_sj_type", "first"),
        "min_side_TPM": ("min_side_TPM", "max"),
        "cfs_narrow_status": ("cfs_narrow_status", "first"),
        "cfs_narrow_tier_bp1": ("cfs_narrow_tier_bp1", "first"),
        "cfs_narrow_tier_bp2": ("cfs_narrow_tier_bp2", "first"),
        "cfs_broad_status": ("cfs_broad_status", "first"),
        "cfs_agreement": ("cfs_agreement", "first"),
        "matched_reference": ("matched_reference",
                              lambda s: s.astype(str).str.lower().eq("true").any()),
    }
    if "missing" in sub:
        agg["missing"] = ("missing", "first")
    if "predictors" in sub:
        agg["predictors"] = ("predictors", lambda s: (
            f"netMHCpan {int(s.isin(['both', 'netMHCpan only']).sum())} · "
            f"MHCflurry {int(s.isin(['both', 'MHCflurry only']).sum())} · "
            f"both {int((s == 'both').sum())}") if (s != "").any() else "")
    return sub.groupby(["line", "sv_id"]).agg(**agg).reset_index()


def near_miss_mask(t: pd.DataFrame) -> pd.Series:
    qc = qc_mask(t)
    tests = {k: test(t) for k, (_, test) in CRITERIA.items()}
    n_fail = sum((~v).astype(int) for v in tests.values())
    t["missing"] = ""
    for k, v in tests.items():
        t.loc[~v, "missing"] = k
    return qc & (n_fail == 1)


def dna_depth(ev: pd.DataFrame, bams: dict[str, str]) -> pd.DataFrame:
    out = []
    for r in ev.itertuples():
        rec = {"line": r.line, "sv_id": r.sv_id, "depth_inside": None,
               "depth_flank": None, "depth_ratio": None, "depth_verdict": "",
               "depth_note": ""}
        bam = bams.get(r.line)
        size = abs(int(r.pos2) - int(r.pos1)) if str(r.chrom1) == str(r.chrom2) else None
        if r.svtype not in DEPTH_TYPES:
            rec["depth_note"] = f"not measured: depth does not test a {r.svtype}"
        elif size is None or size <= 2 * EDGE_TRIM:
            rec["depth_note"] = "not measured: interval too small to trim"
        elif size > DEPTH_MAX_INTERVAL:
            rec["depth_note"] = f"not measured: interval {size:,} bp"
        elif not bam or not pathlib.Path(bam).exists():
            rec["depth_note"] = "not measured: no DNA BAM for this line"
        else:
            lo, hi = sorted((int(r.pos1), int(r.pos2)))
            cn = copy_number(bam, str(r.chrom1), lo, hi)
            rec.update(depth_inside=round(cn["inside"], 1),
                       depth_flank=round(cn["flank"], 1),
                       depth_ratio=round(cn["ratio"], 2))
            if r.svtype == "DEL":
                rec["depth_verdict"] = verdict(cn["ratio"])[0]
            else:
                ratio = cn["ratio"]
                rec["depth_verdict"] = ("undetermined" if ratio != ratio else
                                        "genomic gain" if ratio > 1.3 else
                                        "NOT a genomic gain" if ratio < 1.15 else
                                        "ambiguous")
        out.append(rec)
    return pd.DataFrame(out)


# --------------------------------------------------------------------------- html
def esc(x) -> str:
    return html.escape("" if x is None else str(x))


def md_inline(text: str) -> str:
    """The dictionary's light markdown: `code` and **bold**."""
    s = esc(text)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", s)
    return s


def num(x, digits=0) -> str:
    if x is None or (isinstance(x, float) and x != x) or x == "":
        return "<span class=muted>—</span>"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return esc(x)
    return f"{v:,.{digits}f}"


def cfs_cell(r, null: dict) -> str:
    """Per-breakend fragile-site status. Only the conservative catalogue can flag:
    the permissive one covers most of the genome, so a hit there is near chance."""
    narrow = str(getattr(r, "cfs_narrow_status", "") or "")
    broad = str(getattr(r, "cfs_broad_status", "") or "")
    if narrow and narrow not in ("none", "nan"):
        tiers = sorted({str(x) for x in (getattr(r, "cfs_narrow_tier_bp1", None),
                                         getattr(r, "cfs_narrow_tier_bp2", None))
                        if isinstance(x, str) and x})
        return (f"<span class='pon pon-lo'>CFS · {esc(narrow)}</span>"
                f"<br><small>conservative{' · ' + esc('/'.join(tiers)) if tiers else ''}"
                f"; permissive {esc(broad)}</small>")
    if broad and broad not in ("none", "nan"):
        chance = null.get("broad_expected_hit_rate_pct")
        return (f"<span class=muted>permissive only · {esc(broad)}</span>"
                f"<br><small>{f'chance rate {chance}%' if chance else 'near chance'}</small>")
    return "<span class=muted>none</span>"


def pon_cell(value) -> str:
    """A junction seen in the panel of normals is flagged, never left as a bare number:
    on a branch that does not apply the panel, this is the only place it shows."""
    v = pd.to_numeric(value, errors="coerce")
    if v != v or v <= 0:
        return "<span class=muted>no PON record</span>"
    cls = "pon-hi" if v >= C.PON_MAX else "pon-lo"
    verdict = f"≥ PON_MAX {C.PON_MAX}" if v >= C.PON_MAX else f"&lt; PON_MAX {C.PON_MAX}"
    return f"<span class='pon {cls}'>IN PON · {v:,.0f}</span><br><small>{verdict}</small>"


def plural(n: int, word: str) -> str:
    return f"{n:,} {word}{'' if n == 1 else 's'}"


def method(counted: str, computed: str, limits: str) -> str:
    return (f"<div class=method><div><h4>What is counted</h4><p>{counted}</p></div>"
            f"<div><h4>How it is computed</h4><p>{computed}</p></div>"
            f"<div><h4>What it does not capture</h4><p>{limits}</p></div></div>")


def table(headers: list[str], rows: list[list[str]], cls: str = "") -> str:
    th = "".join(f"<th>{h}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f"<div class=tw><table class='{cls}'><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table></div>"


def chip(line: str, colour_var: str) -> str:
    return f"<span class=chip><i style='background:var({colour_var})'></i>{esc(line)}</span>"


def bar(value: int, total: int, colour_var: str, title: str) -> str:
    pct = 100 * value / total if total else 0
    width = max(pct, 0.6) if value else 0
    return (f"<div class=bar title='{esc(title)}'><span style='width:{width:.2f}%;"
            f"background:var({colour_var})'></span></div><small>{pct:.1f}%</small>")


CSS = """
:root{color-scheme:light;--surface-0:#f4f4f1;--surface-1:#fcfcfb;--surface-2:#eceae4;
--border:#dcd9d0;--text-primary:#0b0b0b;--text-secondary:#52514e;--text-muted:#76746d;
--warn-bg:#fdf3e7;--warn-bd:#eda100;--stop-bg:#fdecec;--stop-bd:#e34948;--ok-bd:#1baf7a;__LIGHT__}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;
--surface-0:#121211;--surface-1:#1a1a19;--surface-2:#232321;--border:#35342f;
--text-primary:#fff;--text-secondary:#c3c2b7;--text-muted:#94928a;--warn-bg:#2b2415;
--warn-bd:#c98500;--stop-bg:#2e1b1b;--stop-bd:#e66767;--ok-bd:#199e70;__DARK__}}
:root[data-theme="dark"]{color-scheme:dark;--surface-0:#121211;--surface-1:#1a1a19;
--surface-2:#232321;--border:#35342f;--text-primary:#fff;--text-secondary:#c3c2b7;
--text-muted:#94928a;--warn-bg:#2b2415;--warn-bd:#c98500;--stop-bg:#2e1b1b;
--stop-bd:#e66767;--ok-bd:#199e70;__DARK__}
*{box-sizing:border-box}
body{margin:0;background:var(--surface-0);color:var(--text-primary);
font:15px/1.6 ui-sans-serif,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
main{max-width:1080px;margin:0 auto;padding:32px 16px 80px}
h1{font-size:28px;line-height:1.2;margin:0 0 8px}h2{font-size:21px;margin:48px 0 8px;
padding-top:16px;border-top:1px solid var(--border)}h3{font-size:16px;margin:28px 0 6px}
h4{font-size:12px;text-transform:uppercase;letter-spacing:.04em;color:var(--text-muted);margin:0 0 4px}
p,li{color:var(--text-secondary)}b,strong{color:var(--text-primary)}
code{font:12.5px ui-monospace,SFMono-Regular,Menlo,monospace;background:var(--surface-2);
padding:1px 4px;border-radius:4px;color:var(--text-primary)}
.muted,small{color:var(--text-muted)}
.meta{color:var(--text-muted);font-size:13px}
nav{display:flex;flex-wrap:wrap;gap:6px 14px;font-size:13px;margin:16px 0}
nav a{color:var(--text-secondary)}
.card{background:var(--surface-1);border:1px solid var(--border);border-radius:10px;padding:16px 18px;margin:14px 0}
.note{border-left:3px solid var(--warn-bd);background:var(--warn-bg);padding:10px 14px;border-radius:6px;margin:14px 0}
.stop{border-left:3px solid var(--stop-bd);background:var(--stop-bg);padding:10px 14px;border-radius:6px;margin:14px 0}
.ok{border-left:3px solid var(--ok-bd);padding:10px 14px;margin:14px 0;background:var(--surface-1);border-radius:6px}
.tw{overflow-x:auto;margin:10px 0}
table{border-collapse:collapse;width:100%;font-size:13.5px;background:var(--surface-1)}
th,td{text-align:left;padding:6px 9px;border-bottom:1px solid var(--border);vertical-align:top}
th{font-size:12px;color:var(--text-muted);font-weight:600;background:var(--surface-2)}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
.chip{display:inline-flex;align-items:center;gap:6px;white-space:nowrap}
.chip i{width:10px;height:10px;border-radius:3px;display:inline-block}
.bar{display:inline-block;width:110px;height:8px;background:var(--surface-2);border-radius:4px;
vertical-align:middle;margin-right:6px;overflow:hidden}
.bar span{display:block;height:100%;border-radius:0 4px 4px 0}
.method{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:14px;
font-size:13px;background:var(--surface-2);border-radius:8px;padding:12px 14px;margin:12px 0}
.method p{margin:0}
.pon{display:inline-block;font-size:12px;font-weight:700;padding:1px 7px;border-radius:5px;white-space:nowrap;color:var(--text-primary)}
.pon-hi{background:var(--stop-bg);border:1px solid var(--stop-bd)}
.pon-lo{background:var(--warn-bg);border:1px solid var(--warn-bd)}
.tag{font-size:11px;border:1px solid var(--border);border-radius:10px;padding:0 7px;
color:var(--text-muted);margin-left:6px;vertical-align:middle;font-weight:500}
details{margin:10px 0}summary{cursor:pointer;color:var(--text-secondary);font-weight:600}
input[type=search]{width:100%;max-width:420px;padding:7px 10px;border:1px solid var(--border);
border-radius:6px;background:var(--surface-1);color:var(--text-primary);font:inherit}
.dict td:first-child{white-space:nowrap}
.dict td:nth-child(3){min-width:320px}
"""


def build(args) -> None:
    cfg = load_config(args.config)
    lines = [s.name for s in cfg.ordered_samples()]
    parents = {s.name: s.parent for s in cfg.samples}
    bams = {s.name: s.dna_bam for s in cfg.samples}
    colour = {l: f"--line{i}" for i, l in enumerate(lines)}
    light = "".join(f"--line{i}:{LINE_COLOURS[i % 8][0]};" for i in range(len(lines)))
    dark = "".join(f"--line{i}:{LINE_COLOURS[i % 8][1]};" for i in range(len(lines)))

    udir = pathlib.Path(args.universe_dir)
    upath = udir / "candidate_universe.tsv"
    t = pd.read_csv(upath, sep="\t", low_memory=False)
    dictionary = pd.read_csv(udir / "candidate_universe_column_dictionary.tsv", sep="\t")
    glossary_path = udir / "candidate_universe_column_dictionary_GLOSSARY.tsv"
    glossary = pd.read_csv(glossary_path, sep="\t") if glossary_path.exists() else pd.DataFrame()

    fun = funnel(t, lines)
    s1 = stage1(pathlib.Path(args.cross_dir), lines, args.branch)

    qc = qc_mask(t)
    conj = qc.copy()
    for _, test in CRITERIA.values():
        conj &= test(t)
    mf_conj = (conj_without_pres(t) & flag(t, "presentable_mhcflurry")
               if "presentable_mhcflurry" in t else pd.Series(False, index=t.index))
    t["predictors"] = [("both" if a and b else "netMHCpan only" if a else
                        "MHCflurry only" if b else "") for a, b in zip(conj, mf_conj)]
    surv = events(t, conj | mf_conj)
    nm_mask = near_miss_mask(t)
    near = events(t, nm_mask)
    # A surviving event usually also has peptides that fail one criterion (most
    # often presentability). Those are not near misses: the event already survives.
    surviving = set(zip(surv["line"], surv["sv_id"].astype(str)))
    near = near[[(l, str(i)) not in surviving
                 for l, i in zip(near["line"], near["sv_id"])]]

    judged = pd.DataFrame(columns=["acquired_in", "sv_id", "judgement", "reason"])
    if args.judgements:
        judged = pd.read_csv(args.judgements, sep="\t", comment="#")
    jkey = {(str(r.acquired_in), str(r.sv_id)): (r.judgement, r.reason)
            for r in judged.itertuples()}

    cfs_prov = udir / "CFS_ANNOTATION_PROVENANCE.tsv"
    cfs_null = (dict(pd.read_csv(cfs_prov, sep="\t").values) if cfs_prov.exists() else {})

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    depth_path = out.with_suffix(".dna_depth.tsv")
    measure = pd.concat([surv, near[near["missing"] != "pres"]]
                        ).drop_duplicates(["line", "sv_id"])
    if args.no_dna_depth:
        depth = pd.DataFrame(columns=["line", "sv_id", "depth_inside", "depth_flank",
                                      "depth_ratio", "depth_verdict", "depth_note"])
    else:
        print(f"measuring DNA depth across {len(measure)} intervals …", file=sys.stderr)
        depth = dna_depth(measure, bams)
        depth.to_csv(depth_path, sep="\t", index=False)
    dkey = {(str(r.line), str(r.sv_id)): r for r in depth.itertuples()}

    S = []   # sections
    N7 = 6 if args.omit_reference_cross else 7     # sections after the optional cross
    roots = [l for l in lines if not parents.get(l)]

    # ---------------------------------------------------------------- header
    lineage_rows = [[chip(l, colour[l]), esc(parents.get(l) or "— (root)"),
                     f"<td class=n>{fun.loc[l,'candidates']:,}", f"{fun.loc[l,'carried']:,}"]
                    for l in lines]
    S.append(f"""
<h1>{esc(' · '.join(lines))}</h1>
<p class=meta>Candidate SV neoantigens per line, from SV call to verified candidate ·
generated {date.today().isoformat()} · branch <code>{esc(args.branch)}</code> ·
universe <code>{sha16(upath)}</code> · contains results, not for publication</p>
<nav><a href=#lineage>Lineage</a><a href=#admission>1 Admission</a><a href=#funnel>2 Funnel</a>
<a href=#verdict>3 Verdict per line</a><a href=#survivors>4 Survivors</a>
<a href=#nearmiss>5 Near misses</a>{'' if args.omit_reference_cross else '<a href=#reference>6 Reference cohort</a>'}
<a href=#notfiltered>{N7} Not filtered</a><a href=#limits>{N7+1} Limitations</a>
<a href=#columns>{N7+2} Column guide</a><a href=#sources>Sources</a></nav>

<h2 id=lineage>Lineage, and the two counting bases</h2>
<p>A derived line inherits its parent's genome, so every count exists on two bases,
and they are not interchangeable:</p>
<ul><li><code>acquired_in</code> — peptides produced by the line's <b>own</b> variant calls,
i.e. what this line added relative to its parent. <b>Every funnel in this report is on
this basis.</b></li>
<li><code>present_in</code> — everything the line carries, inherited included. This is the
basis of the per-line view files. For a root line the two coincide.</li></ul>
{table(["line", "parent", "<span class=n>acquired</span>", "carried (present_in)"],
       [[r[0], r[1], f"{fun.loc[l,'candidates']:,}", f"{fun.loc[l,'carried']:,}"]
        for r, l in zip(lineage_rows, lines)])}
{method("Candidate peptides (rows of the universe, one per distinct 8–11-mer), per line, on the two bases.",
        "<code>acquired_in</code> and <code>present_in</code> are split on <code>;</code> and tested for membership, so a peptide generated independently by two lines counts for both.",
        "A peptide count is not an event count: one junction yields dozens of overlapping windows. The somatic lines are called against their parent, so their acquired sets are what the caller saw as new; an event the caller missed in the parent can surface as 'acquired' in a child.")}
""")

    # ---------------------------------------------------------------- stage 1
    rows = []
    for key, label in [("records", "records in the SV VCF"),
                       ("after_filter", "after FILTER"),
                       ("pon_admitted", "…of which FILTER=PON, admitted deliberately"),
                       ("paired", "after requiring a mate"),
                       ("junctions", "junctions (paired breakends)"),
                       ("peptides", "candidate peptides generated")]:
        rows.append([label] + [num(s1.loc[l, key]) for l in lines])
    S.append(f"""
<h2 id=admission>1 · Admission and peptide generation <span class=tag>filters</span></h2>
<p>Each line's SV VCF is read, records are admitted on FILTER and on having a mate,
and paired breakends are collapsed into junctions. A root line is analysed from its
germline VCF; a derived line from its somatic VCF alone, because its germline genome
is its parent's. On the <code>{esc(args.branch)}</code> branch the panel of normals is
measured and not applied (<code>PON_MAX = {C.PON_MAX}</code> is not enforced).
Every junction then yields a fusion protein and every 8–11-residue window across it
(<code>PEPTIDE_LENGTHS = {tuple(C.PEPTIDE_LENGTHS) if hasattr(C, 'PEPTIDE_LENGTHS') else '(8, 9, 10, 11)'}</code>).</p>
{table(["step"] + [chip(l, colour[l]) for l in lines], rows)}
{method("VCF records (one per breakend), then junctions (two paired breakends), then peptides. The unit changes twice down the table.",
        "Counts read from <code>&lt;cross&gt;/&lt;line&gt;_" + esc(args.branch) + "/summary.json</code> (<code>stage1_funnel</code>) and <code>funnel.tsv</code>. Peptides are translated from the reference sequence in the 5′ partner's frame; nothing is re-derived from reads here.",
        "Admission says a record is well formed, not that the rearrangement is real, acquired, or expressed. Peptides here are generated per sample before lineage attribution, so they can differ from the acquired counts in the funnel below.")}
""")

    # ---------------------------------------------------------------- funnel
    frows = []
    spec = [("candidates", "candidates acquired by this line", ""),
            ("non_self", "not present in the normal proteome", "cumulative"),
            ("qc", "… and not low-complexity", "cumulative")]
    spec += [(k, lab, "marginal") for k, (lab, _) in CRITERIA.items()]
    two = "pres_mf" in fun.columns
    if two:
        spec += [("pres_mf", "presentable by the line's own HLA (MHCflurry)", "marginal"),
                 ("pres_any", "presentable by either predictor", "marginal"),
                 ("pres_both", "presentable by both predictors", "marginal")]
    spec += [("conjunction", "all at once — presentable by netMHCpan" if two
              else "all of the above at once", "conjunction")]
    if two:
        spec += [("conj_mf", "all at once — presentable by MHCflurry", "conjunction"),
                 ("conj_any", "all at once — presentable by either predictor", "conjunction"),
                 ("conj_both", "all at once — presentable by both predictors", "conjunction")]
    for key, label, kind in spec:
        cells = [f"{label} {f'<span class=tag>{kind}</span>' if kind else ''}"]
        for l in lines:
            v, tot = int(fun.loc[l, key]), int(fun.loc[l, "candidates"])
            cells.append(f"<b>{v:,}</b><br>{bar(v, tot, colour[l], f'{l}: {v:,} of {tot:,}')}")
        frows.append(cells)
    S.append(f"""
<h2 id=funnel>2 · Funnel, acquired basis</h2>
<p>Only the <b>conjunction</b> row is quotable as "survivors". Sequence QC nests (each
row a subset of the one above). Confidence, transcription, RNA evidence and
presentability are <b>independent properties of the QC-passed set</b>: each is measured
against that set, not against the row above it, so the marginal rows do not form a
pipeline and one can exceed another.</p>
{table(["step"] + [chip(l, colour[l]) for l in lines], frows)}
{method("Peptides acquired by each line. Bars are the share of that line's own acquired candidates, so they compare proportions between lines and never absolute numbers.",
        f"<code>is_self</code>: exact substring of the Ensembl proteome. <code>low_complexity</code>: entropy, dominant residue, homopolymer or distinct-residue rule. <code>sv_hc</code>: segment MAPQ ≥ {C.HC_MIN_SEGMAPQ}, ≥ {C.HC_MIN_VF} supporting fragments, QUAL ≥ {C.HC_MIN_QUAL}. <code>expressed</code>: both partners above the TPM floor. RNA: <code>rna_tier</code> ∈ {{{', '.join(RNA_CROSSING)}}}, i.e. reads that <i>cross</i> the junction by the mechanism the junction's geometry allows. <code>presentable</code>: IC50 ≤ 500 nM AND %Rank_BA ≤ 2 AND %Rank_EL ≤ 2 for ≥ 1 of the line's own class I alleles (netMHCpan 4.2e). <code>presentable_mhcflurry</code>: affinity ≤ 500 nM AND presentation percentile ≤ 2 for ≥ 1 own allele (MHCflurry 2.2.1 presentation model). Four presentability criteria, each with its conjunction view: netMHCpan (<code>credible_and_presentable</code>), MHCflurry (<code>_mhcflurry</code>), either predictor (<code>_either</code>) and both (<code>_both</code>). None is ranked above the others: requiring both is stricter, not better — the two predictors disagree on many peptides and neither is a measurement, so a peptide presented by only one is a candidate, not a reject. The per-line verdict and section 4 use the either-predictor conjunction and state which predictor presents each peptide.",
        "Overlapping windows are not independent: a share per peptide should be re-checked per event (section 4 does). The marginal rows say nothing about what survives alongside them. Predicted binding is not presentation.")}
""")

    # ---------------------------------------------------------------- verdict per line
    vrows = []
    for l in lines:
        q = int(fun.loc[l, "qc"])
        n = int(fun.loc[l, "conj_any"]) if "conj_any" in fun.columns else int(fun.loc[l, "conjunction"])
        split = (f" (netMHCpan {int(fun.loc[l, 'conjunction'])}, MHCflurry "
                 f"{int(fun.loc[l, 'conj_mf'])}, both {int(fun.loc[l, 'conj_both'])})"
                 if "conj_any" in fun.columns else "")
        s_ev = surv[surv.line == l]
        kept = [r for r in s_ev.itertuples() if (l, str(r.sv_id)) not in jkey]
        if n:
            def in_dna(r):
                d = dkey.get((l, str(r.sv_id)))
                return d is not None and str(d.depth_verdict).startswith("genomic")
            def common(r):
                pon = pd.to_numeric(r.pon_count, errors="coerce")
                af = pd.to_numeric(r.gnomad_af_popmax, errors="coerce")
                return (pon == pon and pon >= C.PON_MAX) or (af == af and af >= C.POPULATION_COMMON_AF)
            confirmed = [r for r in kept if in_dna(r)]
            private = [r for r in confirmed if not common(r)]
            text = (f"<b>{n}</b> peptides{split} from <b>{plural(len(s_ev), 'event')}</b> → "
                    f"<b>{len(kept)}</b> after judgements → <b>{len(confirmed)}</b> with the "
                    f"lesion confirmed by DNA depth → <b>{len(private)}</b> of those not common "
                    f"(PON &lt; {C.PON_MAX} and gnomAD popmax &lt; {C.POPULATION_COMMON_AF})")
        else:
            zero = [CRITERIA[k][0] for k in CRITERIA if fun.loc[l, k] == 0]
            smallest = min(CRITERIA, key=lambda k: fun.loc[l, k])
            why = (f"no QC-passed peptide has {', '.join(zero)}" if zero else
                   f"each criterion is met by some peptides (smallest: {CRITERIA[smallest][0]}, "
                   f"{fun.loc[l, smallest]:,} of {q:,}), but never all four by the same one")
            text = f"<b>no candidate satisfies every criterion</b> — {why}"
        vrows.append([chip(l, colour[l]), f"{q:,}", text])
    S.append(f"""
<h2 id=verdict>3 · Verdict per line</h2>
{table(["line", "QC-passed (acquired)", "result"], vrows)}
{method("Per line: the conjunction count in peptides and in events (junctions), then three successive event counts: after the per-event judgements of section 4, with DNA depth across the interval matching a genomic lesion, and of those, the ones that are not common variation. Each arrow is a subset of the one before.",
        "An event is one <code>(acquired_in, sv_id)</code> pair. Where nothing survives, the marginal counts of section 2 say whether one criterion is empty on its own or whether the four simply never coincide; section 5 lists the events that miss by exactly one.",
        "Zero survivors is a statement about this call set and these thresholds, not proof that the line carries no SV neoantigen: an expressed junction below the RNA tier, an unannotated allele or a missed call would all be invisible here.")}
""")

    # ---------------------------------------------------------------- survivors
    def ev_rows(ev: pd.DataFrame, with_missing=False) -> list[list[str]]:
        rows = []
        for r in ev.sort_values(["line", "genes", "sv_id"]).itertuples():
            d = dkey.get((r.line, str(r.sv_id)))
            j = jkey.get((r.line, str(r.sv_id)))
            depth_cell = ("<span class=muted>" + esc(d.depth_note) + "</span>"
                          if d is not None and d.depth_note else
                          f"{num(d.depth_inside,1)} / {num(d.depth_flank,1)} = <b>{num(d.depth_ratio,2)}</b><br><small>{esc(d.depth_verdict)}</small>"
                          if d is not None else "<span class=muted>not measured</span>")
            row = [chip(r.line, colour.get(r.line, "--line0")),
                   f"<b>{esc(r.genes)}</b>" + (f"→{esc(r.gene2)}" if r.gene2 != r.genes else "")
                   + (f"<br><small>presentable peptides: {esc(r.predictors)}</small>"
                      if getattr(r, "predictors", "") else ""),
                   f"{esc(r.sv_id)}<br><small>{esc(r.svtype)} {num(r.event_size)} bp</small>",
                   num(r.peptides)]
            if with_missing:
                row.append(esc(CRITERIA[r.missing][0]) if r.missing in CRITERIA else "")
            row += [num(r.junction_reads) + f"<br><small>{esc(r.test)}</small>",
                    num(r.junction_usage, 2) + f"<br><small>of {num(r.junction_interval_total)}</small>",
                    f"{num(r.vf_bp1)} / {num(r.ref_bp1)}",
                    f"{num(r.dna_vaf_bp1, 2)} · {num(r.dna_vaf_bp2, 2)}",
                    depth_cell,
                    pon_cell(r.pon_count) + f"<br><small>gnomAD {num(r.gnomad_af_popmax, 3)}</small>",
                    num(r.nearest_alt_sj_bp),
                    cfs_cell(r, cfs_null)]
            if not with_missing:
                row.append(f"<b>{esc(j[0])}</b>" if j else "<span class=muted>none</span>")
            rows.append(row)
        return rows

    head = ["line", "gene", "sv_id", "pep"]
    tail = ["RNA reads", "usage", "DNA var / ref", "DNA VAF bp1 · bp2",
            "DNA depth inside / flanks", "panel of normals", "alt SJ (bp)",
            "fragile site"]
    S.append(f"""
<h2 id=survivors>4 · Verifying the survivors <span class=tag>subtracts</span></h2>
<p>This is where the report stops adding and starts subtracting. Each surviving event is
checked against the DNA: <b>is the lesion in the genome?</b> Read counts in RNA cannot
answer that — a deleted allele and an alternative splice junction both produce the same
gapped read — but DNA depth inside the interval against its flanks can, because
expression does not enter it.</p>
<p class=note><b>Panel of normals.</b> On this branch the panel is measured and not applied, so an
event seen in normal genomes can survive the funnel. Every event with <code>pon_count</code> &gt; 0 is
marked <span class='pon pon-hi'>IN PON</span> (red: at or above <code>PON_MAX</code>, the cut the PON branch
would apply; amber: below it). Such an event is inherited human variation, not something the line acquired.</p>
{table(head + tail + ["judgement"], ev_rows(surv)) if len(surv) else "<p class=stop>No line has a surviving event.</p>"}
{method("One row per surviving event: its peptides, RNA support, DNA allele fraction at the breakpoint, and DNA depth across the interval. The judgement column is a person's verdict, not a pipeline output.",
        f"Depth: mean per-base DNA coverage (no quality threshold) inside the interval trimmed {EDGE_TRIM} bp at each edge, over the mean of two {FLANK_BP} bp flanks starting 100 bp outside each breakend, on the acquiring line's own DNA BAM. Expectations: DEL {DEPTH_TYPES['DEL']}; DUP {DEPTH_TYPES['DUP']}. Fragile site: breakend status (<code>none/bp1/bp2/both</code>) under the conservative catalogue ({cfs_null.get('narrow_regions', '?')} regions, {cfs_null.get('narrow_genome_coverage_pct', '?')}% of the genome, {cfs_null.get('narrow_expected_hit_rate_pct', '?')}% chance hit rate for two random breakends) and the permissive one ({cfs_null.get('broad_regions', '?')} regions, {cfs_null.get('broad_genome_coverage_pct', '?')}%, {cfs_null.get('broad_expected_hit_rate_pct', '?')}%); only a conservative hit is flagged, a permissive-only hit is at chance level. Written to <code>{esc(depth_path.name)}</code>. Usage: <code>junction_reads ÷ junction_interval_total</code>, in fragments. 'alt SJ' is the distance to the nearest splice junction Isofox calls independently. Judgements come from <code>{esc(pathlib.Path(args.judgements).name) if args.judgements else 'no judgements file'}</code>.",
        "Depth needs DNA coverage and a lesion large enough to trim; otherwise it is undetermined, never negative. It establishes the lesion, not that the transcript is translated or the peptide presented. A ratio near 1 on a cell line can also mean the deletion is subclonal. The judgements are applied by hand and nothing in the universe is filtered on them, so the first count in section 3 is the pipeline's and the second is a reading of it.")}
""")
    if len(judged):
        grouped = judged.groupby(["acquired_in", "judgement", "reason"], sort=False)["sv_id"] \
            .apply(lambda s: ", ".join(map(str, s))).reset_index()
        S.append("<h3>Judgements applied</h3>" + table(
            ["line", "events (sv_id)", "judgement", "reason"],
            [[esc(r.acquired_in), esc(r.sv_id), f"<b>{esc(r.judgement)}</b>", md_inline(r.reason)]
             for r in grouped.itertuples()]))

    # ---------------------------------------------------------------- near misses
    nm_rows = near[near["missing"] != "pres"]
    pres_only = near[near["missing"] == "pres"]
    pres_summary = ", ".join(f"{l} {int((pres_only.line == l).sum())}" for l in lines)
    S.append(f"""
<h2 id=nearmiss>5 · Near misses: events failing exactly one criterion</h2>
<p>Where a line has no survivor, the useful question is how close it came. These are
QC-passed peptides that meet three of the four criteria of the conjunction, collapsed to
events, excluding events that already survive through another of their peptides. Events that miss <b>only on presentability</b> are counted but not listed
(per line: {esc(pres_summary)}): they are real transcribed junctions whose peptides
this genotype cannot present, and DNA depth is not measured for them.</p>
{table(head + ["fails on"] + tail, ev_rows(nm_rows, with_missing=True)) if len(nm_rows) else "<p>None.</p>"}
{method("Events with ≥ 1 QC-passed peptide that fails exactly one of the four criteria. The criterion is named per event.",
        "Same columns and depth measurement as section 4. A peptide failing two or more criteria is not listed.",
        "A near miss is not a candidate. Relaxing the one failing criterion would admit it, and that is a decision about the threshold, to be taken knowingly — not a finding.")}
""")

    # ---------------------------------------------------------------- reference
    rrows = [[chip(l, colour[l]), f"{fun.loc[l,'candidates']:,}",
              f"{fun.loc[l,'matched_reference']:,}",
              f"{100*fun.loc[l,'matched_reference']/fun.loc[l,'candidates']:.1f}%" if fun.loc[l, 'candidates'] else "—",
              ", ".join(sorted({f"{g}" for g in surv[(surv.line == l) & surv.matched_reference].genes})) or "—"]
             for l in lines]
    if not args.omit_reference_cross:
        S.append(f"""
<h2 id=reference>6 · Cross against the reference patient cohort <span class=tag>no filter</span></h2>
{table(["line", "acquired", "also in the reference cohort", "share", "survivor events matched"], rrows)}
{method("Acquired peptides whose exact sequence also occurs in the reference patient catalogue (<code>matched_reference</code>).",
        "Exact string identity, HLA-independent. The matched event's gene and SV type are carried in <code>gene_concordant</code> / <code>svtype_concordant</code> to separate a shared event from convergent sequence.",
        "Absence from one finite catalogue is not rarity in patients, and a match is not evidence any patient presents it — that needs their genotype (the <code>*_patient</code> columns). A high share in a root line mostly reflects common germline variation shared with patients.")}
""")

    # ---------------------------------------------------------------- not filtered + limitations
    pon_line = "; ".join(
        f"{l}: {int((flag(t, 'pass_pon') == False)[acquired_by(t, l)].sum()):,} of {fun.loc[l,'candidates']:,} peptides on a junction failing the panel"
        for l in lines)
    S.append(f"""
<h2 id=notfiltered>{N7} · What none of this filters on</h2>
<ul>
<li><b>The panel of normals.</b> Measured (<code>pon_count</code>, <code>pass_pon</code>) and not applied
on this branch — {esc(pon_line)}. Filtering on it is a decision to take with <code>pon_count</code> in view.</li>
<li><b>Population frequency.</b> <code>gnomad_af_popmax</code> and <code>is_highfreq_gnomad</code> are recorded;
empty <code>gnomad_af_popmax</code> means no gnomAD record, not rarity.</li>
{'' if args.omit_reference_cross else '<li><b>Membership of the patient cohort</b> (<code>matched_reference</code>).</li>'}
<li><b>Fragile sites</b> (<code>cfs_*</code>) — annotated per breakend under two catalogues, neither selects.</li>
<li><b>Junction usage</b> — a share whose denominator must be read with it.</li>
</ul>
<h2 id=limits>{N7+1} · Limitations</h2>
<ul>
<li><code>presentable</code> is predicted binding, not observed presentation: processing, transport,
surface abundance and T-cell recognition are not addressed. Only immunopeptidomics would.</li>
<li>The self test is exact identity against the reference proteome; a peptide one residue from self,
or self in this line through its own germline variation, is not caught.</li>
<li>A negative RNA result is not absence of the lesion: the gene may be silent or the transcript
degraded by nonsense-mediated decay.</li>
<li>The derived lines are called somatically against their parent. An event the caller missed in the
parent can appear as acquired in a child; DNA depth in the parent is the check.</li>
<li>Binding predictions are netMHCpan 4.2e; individual binder identities are sensitive to the predictor version,
aggregate magnitudes much less so.</li>
</ul>
""")

    # ---------------------------------------------------------------- column guide
    gl = ""
    if len(glossary):
        gl = "<details><summary>Glossary — terms used throughout the column guide</summary>" + table(
            ["term", "meaning"], [[f"<b>{esc(r.term)}</b>", md_inline(r.plain_english)]
                                  for r in glossary.itertuples()]) + "</details>"
    groups = []
    for stage, g in dictionary.groupby("stage", sort=False):
        rows = []
        for r in g.itertuples():
            also = f"<br><small>{md_inline(r.also_noted)}</small>" if isinstance(r.also_noted, str) and r.also_noted else ""
            rows.append([f"<code>{esc(r.column)}</code><br><small>{esc(r.kind)} · {esc(r.level)}</small>",
                         f"{num(r.pct_populated, 1)}%",
                         md_inline(r.what_it_is) + also,
                         f"<small>{md_inline(r.derived_from)}</small>"])
        groups.append(f"<h3>{esc(stage)} <span class=tag>{len(g)} columns</span></h3>"
                      + table(["column", "populated", "what it is", "derived from"], rows, "dict"))
    S.append(f"""
<h2 id=columns>{N7+2} · Column guide — <code>candidate_universe.tsv</code></h2>
<p>All {len(dictionary)} columns, grouped by the pipeline stage that produces them, rendered from
<code>candidate_universe_column_dictionary.tsv</code> (which imports every threshold from the
module that applies it, so the values cannot drift from the code). <i>Populated</i> is the share
of all {len(t):,} rows with a value; an empty cell means the question was not asked
(e.g. no cohort match, no RNA BAM), never "no".</p>
{gl}
<p><input type=search id=colfilter placeholder="Filter columns, e.g. junction, gnomad, presentable"></p>
<div id=coltables>{''.join(groups)}</div>
""")

    # ---------------------------------------------------------------- sources
    prov = code_version()
    S.append(f"""
<h2 id=sources>Sources and provenance</h2>
{table(["file", "what it holds"], [
    [f"<code>{esc(upath)}</code>", f"{len(t):,} peptide rows × {t.shape[1]} columns; sha256 {sha16(upath)}…"],
    [f"<code>{esc(udir / 'views' / 'MANIFEST.tsv')}</code>", "the named per-line views on the carried basis"],
    [f"<code>{esc(pathlib.Path(args.cross_dir))}/&lt;line&gt;_{esc(args.branch)}/</code>", "stage-1 admission (summary.json, funnel.tsv)"],
    [f"<code>{esc(depth_path)}</code>", "DNA depth measured for this report"],
    [f"<code>{esc(args.judgements or '—')}</code>", "per-event judgements (private)"],
    [f"<code>{esc(args.config)}</code>", "run config: lines, lineage, BAMs"],
])}
{table(["key", "value"], [[esc(r['key']), esc(r['value'])] for r in prov])}
""")

    script = """<script>
const f=document.getElementById('colfilter');
f&&f.addEventListener('input',()=>{const q=f.value.toLowerCase();
document.querySelectorAll('#coltables tbody tr').forEach(tr=>{
tr.style.display=tr.textContent.toLowerCase().includes(q)?'':'none'});});
</script>"""
    css = CSS.replace("__LIGHT__", light).replace("__DARK__", dark)
    page = (f"<!doctype html><html lang=en><head><meta charset=utf-8>"
            f"<meta name=viewport content='width=device-width,initial-scale=1'>"
            f"<title>Lineage neoantigen report</title><style>{css}</style></head>"
            f"<body><main>{''.join(S)}</main>{script}</body></html>")
    out.write_text(page)
    print(f"wrote {out}  ({len(page)/1024:.0f} KB)", file=sys.stderr)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--config", required=True)
    p.add_argument("--universe-dir", required=True)
    p.add_argument("--cross-dir", required=True)
    p.add_argument("--branch", default="noPON")
    p.add_argument("--judgements", default=None)
    p.add_argument("--out", required=True)
    p.add_argument("--no-dna-depth", action="store_true")
    p.add_argument("--omit-reference-cross", action="store_true",
                   help="leave out the patient-cohort cross (section 6), e.g. while the "
                        "reference catalogue was built with a generator the samples' "
                        "peptides no longer share; the universe keeps the columns")
    build(p.parse_args())


if __name__ == "__main__":
    main()
