#!/usr/bin/env python3
"""
build_report.py — a technical results report for one sample, from its outputs.

Generated rather than written so that every figure in the prose comes from the
run it describes, and so that reporting a second sample is not a second writing
job. Refining the report means editing this file and re-running it.

Nothing here recomputes anything: it reads what the pipeline wrote. If a number
in a report is wrong, the pipeline produced it, and the report is faithful.

    python tools/build_report.py --run-dir results/parental_PON10
    python tools/build_report.py --all          # every run in results/
"""
from __future__ import annotations

import argparse
import json
import os
import pathlib

import pandas as pd

RNA_SUPPORTED = ("STRONG", "SUGGESTIVE", "WEAK")


def read(path, **kwargs):
    return pd.read_csv(path, sep="\t", low_memory=False, **kwargs) \
        if os.path.exists(path) else pd.DataFrame()


def fmt(value, nd=0):
    """Numbers with thousands separators; NA stays visibly NA."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "NA"
    if isinstance(value, str):
        return value
    return f"{value:,.{nd}f}"


def pct(part, whole):
    return "—" if not whole else f"{100.0 * part / whole:.1f}%"


def section(title, level=2):
    return f"\n{'#' * level} {title}\n"


def table(frame: pd.DataFrame, columns=None, rename=None) -> str:
    """Markdown table from a DataFrame, empty-safe."""
    if frame is None or frame.empty:
        return "_No rows._\n"
    view = frame[[c for c in (columns or frame.columns) if c in frame.columns]].copy()
    if rename:
        view = view.rename(columns=rename)
    header = "| " + " | ".join(str(c) for c in view.columns) + " |"
    rule = "|" + "|".join("---" for _ in view.columns) + "|"
    rows = []
    for _, row in view.iterrows():
        cells = []
        for value in row:
            if isinstance(value, float):
                if pd.isna(value):
                    cells.append("—")
                elif float(value).is_integer():
                    cells.append(f"{value:,.0f}")
                elif abs(value) >= 1000:
                    cells.append(f"{value:,.0f}")      # 2.06e+04 reads as noise
                else:
                    cells.append(f"{value:,.4g}")
            else:
                cells.append("—" if pd.isna(value) else str(value))
        rows.append("| " + " | ".join(cells) + " |")
    return "\n".join([header, rule, *rows]) + "\n"


def build(run_dir: str) -> str:
    run = pathlib.Path(run_dir)
    summary = json.loads((run / "summary.json").read_text()) \
        if (run / "summary.json").exists() else {}
    counts = summary.get("counts", {})
    sample = summary.get("sample", run.name)
    branch = summary.get("branch", "")
    kind = summary.get("vcf_kind", "")

    junctions = read(run / "stage1_junctions.tsv")
    peptides = read(run / f"{sample}.all_neopeptides.txt")
    matches = read(run / "stage3_matches.tsv")
    events = read(run / "credible_events.tsv")
    rna = read(run / "stage8_rna_evidence.tsv")
    expression = read(run / "stage7_expression.tsv")
    sv_master = read(run / "master_sv.tsv")
    background = json.loads((run / "stage7_background.json").read_text()) \
        if (run / "stage7_background.json").exists() else {}
    funnel = read(run / "funnel.tsv")

    out = []
    W = out.append

    # ---------------------------------------------------------------- header
    W(f"# {sample} — SV-neoantigen recurrence report\n")
    W(f"**Call set:** {kind} · **Branch:** {branch or 'default'} · "
      f"**Reference catalogue:** {fmt(summary.get('null_model', {}).get('n_reference'))} "
      f"curated patient peptides\n")

    modalities = summary.get("modalities", {})
    available = [k for k, v in modalities.items() if v]
    missing = [k for k, v in modalities.items() if not v]
    W(f"**Modalities available:** {', '.join(available) or 'none'}"
      + (f" · **absent:** {', '.join(missing)}" if missing else "") + "\n")

    generator = summary.get("peptide_generator", {})
    W(f"**Peptide generator:** `{generator.get('generator')}` "
      f"({generator.get('version')})\n")

    # ---------------------------------------------------------------- headline
    events_n = counts.get("events", 0)
    private = counts.get("events_private", "NA")
    hc = counts.get("events_hc", "NA")
    rna_n = counts.get("events_rna_supported", "NA")
    hc_and_rna = counts.get("events_hc_and_rna", "NA")
    all_three = counts.get("events_private_hc_and_rna", "NA")

    W(section("Result"))
    W(f"Of {fmt(counts.get('vcf_records'))} structural-variant records, "
      f"{fmt(counts.get('matches_unique_peptides', counts.get('matches_identical')))} "
      f"candidate peptides are identical to a catalogue entry, collapsing to "
      f"**{fmt(events_n)} credible genomic events**. Of those, {fmt(private)} are "
      f"private to this sample, {fmt(hc)} are high-confidence SV calls, "
      f"{fmt(rna_n)} carry junction-crossing reads in RNA, "
      f"{fmt(hc_and_rna)} are high-confidence *and* transcribed, and "
      f"**{fmt(all_three)} satisfy all three**.\n\n"
      f"The last two figures are different questions and must not be conflated: "
      f"a confidently called, transcribed junction that is common in the "
      f"population is a polymorphism, not a recurrent tumour neoantigen.\n")

    null = summary.get("null_model", {})
    if null.get("p_value") is not None:
        p = (f"< {1 / null['n_permutations']:g}" if null.get("p_value_is_bounded")
             else f"= {null['p_value']:.3g}")
        W(f"The match count is **{null.get('enrichment_over_null', '—')}× the "
          f"permutation null** ({null.get('null_mean')} ± {null.get('null_sd')} "
          f"expected, p {p}), so the recurrence itself is not chance. "
          f"Match rate: **{null.get('matches_per_1000_candidates')} per 1,000 "
          f"candidate peptides** — the figure to compare across samples, since "
          f"candidate universes differ by orders of magnitude.\n")

    # ---------------------------------------------------------------- funnel
    W(section("Evidence funnel"))
    W("Where each order of magnitude is lost. `NA` means a stage did not run, "
      "never that it failed.\n")
    W(table(funnel, ["step", "description", "n"],
            {"step": "Step", "description": "What it asks", "n": "n"}))

    # ---------------------------------------------------------------- stage 1
    W(section("Stage 1 — Admission"))
    s1 = summary.get("stage1_funnel", {})
    W(table(pd.DataFrame([{"filter": k, "records remaining": v} for k, v in s1.items()]),
            rename={"filter": "Rule applied", "records remaining": "Records remaining"}))

    if not junctions.empty:
        # BND has no defined size, so its median is over an all-NaN column;
        # numpy warns about the empty slice. The NA is the correct answer.
        import warnings
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        by_type = (junctions.groupby("svtype")
                   .agg(junctions=("sample_sv_id", "count"),
                        median_size=("event_size", "median"),
                        with_insert=("insert_len", lambda s:
                                     int((pd.to_numeric(s, errors="coerce") >= 10).sum())))
                   .reset_index().sort_values("junctions", ascending=False))
        W(f"\n{fmt(len(junctions))} admitted junctions, by type. `with_insert` counts "
          f"those carrying ≥10 inserted bases — these need the insertion test, not a "
          f"gap test, and the breakend span misrepresents their size.\n")
        W(table(by_type, rename={"svtype": "SV type", "junctions": "n",
                                 "median_size": "median size (bp)",
                                 "with_insert": "insertion-driven"}))

        pon = pd.to_numeric(junctions.get("pon_count"), errors="coerce")
        W(f"\nPanel of normals: {fmt(int(pon.notna().sum()))} of {fmt(len(junctions))} "
          f"junctions carry a `PON_COUNT` ({pct(int(pon.notna().sum()), len(junctions))}); "
          f"median {fmt(pon.median())}, maximum {fmt(pon.max())}.\n")

    # ---------------------------------------------------------------- stage 2
    W(section("Stage 2 — Candidate peptide generation"))
    if not peptides.empty:
        per_sv = peptides.groupby("sv_id")["neopeptide"].nunique()
        W(f"{fmt(len(peptides))} peptide rows → "
          f"{fmt(peptides.neopeptide.nunique())} distinct peptides, from "
          f"{fmt(peptides.sv_id.nunique())} junctions.\n")
        if not sv_master.empty:
            productive = int((pd.to_numeric(sv_master.n_peptides_generated,
                                            errors="coerce") > 0).sum())
            W(f"\n**Only {fmt(productive)} of {fmt(len(sv_master))} admitted junctions "
              f"({pct(productive, len(sv_master))}) produce any peptide at all** — the "
              f"rest do not disrupt a coding transcript. This is the denominator that "
              f"makes a match count interpretable.\n")
        W(f"\nPeptides per productive junction: median {fmt(per_sv.median())}, "
          f"maximum {fmt(per_sv.max())}. A sliding window over one junction yields many "
          f"overlapping peptides, which is why findings are counted as events.\n")
        if "spans_junction" in peptides.columns:
            spans = peptides.spans_junction.astype(str).str.lower().eq("true").sum()
            W(f"\n{fmt(int(spans))} of {fmt(len(peptides))} peptide rows span the "
              f"breakpoint ({pct(int(spans), len(peptides))}). Reported, not filtered: "
              f"after a frameshift, downstream residues are SV-derived without crossing "
              f"the junction.\n")

    # ------------------------------------------------- patient recurrence
    if not sv_master.empty and "patient_evidence" in sv_master.columns:
        W(section("Was this call set seen in patients?"))
        W("Every admitted junction, by the strongest level of patient recurrence "
          "it reaches. The three levels are **not the same kind of evidence**: an "
          "identical peptide is a shared consequence, a shared gene and SV type is "
          "a shared mechanism, and a nearby breakpoint is only a shared locus — at "
          "a fragile site, possibly no more than a shared propensity to break.\n")
        order = ["identical_peptide", "same_gene_and_svtype", "same_gene",
                 "breakpoint_within_1kb", "breakpoint_within_10kb",
                 "breakpoint_within_100kb", "not_seen_in_patients"]
        counts_by = sv_master.patient_evidence.value_counts()
        rows = [{"level": level, "junctions": int(counts_by.get(level, 0)),
                 "share": pct(int(counts_by.get(level, 0)), len(sv_master))}
                for level in order if counts_by.get(level, 0)]
        W(table(pd.DataFrame(rows), rename={"level": "Strongest level reached",
                                            "junctions": "Junctions", "share": "Share"}))
        distances = pd.to_numeric(sv_master.get("patient_bp_dist_bp"), errors="coerce")
        if distances.notna().any():
            W(f"\nDistance to the nearest patient breakpoint, over "
              f"{fmt(int(distances.notna().sum()))} junctions with one on the same "
              f"chromosome: median {fmt(distances.median())} bp, "
              f"{fmt(int((distances <= 1000).sum()))} within 1 kb.\n")
            W("\nThe distance is reported unthresholded in `master_sv.tsv` "
              "(`patient_bp_dist_bp`): 300 bp and 90 kb both fall inside "
              "'within 100 kb' and do not mean the same thing.\n")

    # ---------------------------------------------------------------- stage 3-5
    W(section("Stages 3–5 — Cross, sequence QC, collapse to events"))
    if matches.empty:
        W("No candidate peptide is identical to a catalogue entry. The cross ends here; "
          "downstream stages have nothing to evaluate.\n")
    else:
        concordant = matches.gene_concordant.astype(str).str.lower().eq("true").sum()
        lc = matches.low_complexity.astype(str).str.lower().eq("true").sum()
        selfp = matches.is_self.astype(str).str.lower().eq("true").sum()
        credible = matches.credible.astype(str).str.lower().eq("true").sum()
        W(f"{fmt(len(matches))} match rows ({fmt(matches.peptide.nunique())} distinct "
          f"peptides). One peptide can be produced by several junctions, which is why "
          f"both grains are reported.\n")
        W(f"\n- **gene-concordant:** {fmt(int(concordant))} "
          f"({pct(int(concordant), len(matches))}) — evaluated against *both* breakends\n"
          f"- **low-complexity:** {fmt(int(lc))} removed\n"
          f"- **self-proteome:** {fmt(int(selfp))} removed\n"
          f"- **credible:** {fmt(int(credible))} rows → "
          f"**{fmt(counts.get('events'))} genomic events**\n")
        if not events.empty and "n_peptides" in events.columns:
            top = events.nlargest(5, "n_peptides")[["sample_sv_id", "gene", "n_peptides"]]
            W(f"\nPeptides collapsed per event (the inflation a peptide count would "
              f"report as independent findings):\n")
            W(table(top, rename={"sample_sv_id": "SV id", "gene": "catalogue gene",
                                 "n_peptides": "peptides"}))

    # ---------------------------------------------------------------- stage 6-8
    if not events.empty:
        W(section("Stages 6–8 — Confidence, privacy, and RNA evidence"))
        merged = events.copy()
        for extra in (expression, rna):
            if not extra.empty:
                cols = [c for c in extra.columns if c not in ("gene",)]
                merged = merged.merge(extra[cols], on="sample_sv_id", how="left",
                                      suffixes=("", "_dup"))
        merged = merged.drop(columns=[c for c in merged.columns if c.endswith("_dup")])

        W("Every credible event, with the evidence behind each verdict. "
          "`is_private` requires clearing **both** the panel of normals and the "
          "population frequency; an event can be a flawless call and still be a "
          "common polymorphism.\n")
        W(table(merged.sort_values("rna_tier" if "rna_tier" in merged.columns
                                   else "sample_sv_id"),
                ["gene", "svtype", "event_size", "insert_len", "pon_count", "gnomad_af",
                 "sv_hc", "is_private", "min_side_TPM", "test", "junction_reads",
                 "rna_tier"],
                {"gene": "Gene", "svtype": "Type", "event_size": "Size (bp)",
                 "insert_len": "Insert", "pon_count": "PON", "gnomad_af": "gnomAD AF",
                 "sv_hc": "HC", "is_private": "Private", "min_side_TPM": "TPM (min side)",
                 "test": "RNA test", "junction_reads": "Junction reads",
                 "rna_tier": "Tier"}))

        if not rna.empty and "test" in rna.columns:
            by_test = rna.groupby("test").size().reset_index(name="events")
            W(f"\nRNA test selected by event geometry — an event with no applicable "
              f"test is `UNTESTABLE`, never `NONE`:\n")
            W(table(by_test, rename={"test": "Test", "events": "Events"}))

        if background:
            W(f"\n**Expression background.** Of "
              f"{fmt(background.get('n_genes_quantified'))} quantified genes, "
              f"{background.get('TPM>0.0')}% clear TPM > 0 and "
              f"{background.get('TPM>1.0')}% clear TPM > 1. A bare 'expressed' claim at "
              f"TPM > 0 therefore carries almost no information; the disrupted "
              f"isoform's TPM is the honest figure.\n")

    # ---------------------------------------------------------------- notes
    notes = summary.get("notes", [])
    if notes:
        W(section("Cautions emitted with these numbers"))
        for note in notes:
            W(f"- {note}\n")

    # ---------------------------------------------------------------- provenance
    W(section("Provenance"))
    criteria = summary.get("criteria", {})
    W(f"Produced by `svneo` from `{run.name}`. The complete threshold manifest "
      f"({len([k for k in criteria if not k.startswith('_')])} criteria) is in "
      f"`summary.json` alongside these outputs, so no figure here can outlive the "
      f"criteria that produced it.\n")
    key = ["PON_MAX", "GNOMAD_MAX_AF", "STRONG_MIN_JUNCTION_READS", "MIN_READ_MAPQ",
           "MIN_TESTABLE_GAP_SIZE", "HC_MIN_SEGMAPQ", "HC_MIN_VF", "HC_MIN_QUAL",
           "HC_MIN_SV_SIZE", "NULL_PERMUTATIONS"]
    W("\nThresholds that decided the verdicts above:\n")
    W(table(pd.DataFrame([{"criterion": k, "value": criteria.get(k)}
                          for k in key if k in criteria]),
            rename={"criterion": "Criterion", "value": "Value"}))

    W("\n**Reading these numbers.** Counts are events, not peptide rows, except where "
      "stated. `NA` distinguishes 'not evaluated' from 'evaluated and failed'. Read "
      "evidence supports the claim that a junction is *transcribed*, never that a "
      "peptide is *presented* — that requires immunopeptidome mass spectrometry. "
      "Results here are not stratified by transcript strand.\n")
    return "".join(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", help="a single results/<sample> directory")
    ap.add_argument("--all", action="store_true", help="every run under --results-dir")
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--out-dir", default=None, help="default: <results-dir>/reports")
    args = ap.parse_args()

    runs = []
    if args.run_dir:
        runs = [args.run_dir]
    elif args.all:
        runs = [os.path.join(args.results_dir, d)
                for d in sorted(os.listdir(args.results_dir))
                if os.path.isdir(os.path.join(args.results_dir, d))]
    else:
        ap.error("give --run-dir or --all")

    out_dir = args.out_dir or os.path.join(args.results_dir, "reports")
    os.makedirs(out_dir, exist_ok=True)
    for run in runs:
        if not os.path.exists(os.path.join(run, "summary.json")):
            print(f"  skipped {os.path.basename(run)} (no summary.json)")
            continue
        text = build(run)
        path = os.path.join(out_dir, f"{os.path.basename(run)}.md")
        with open(path, "w") as handle:
            handle.write(text)
        print(f"  wrote {path}  ({len(text.splitlines())} lines)")


if __name__ == "__main__":
    main()
