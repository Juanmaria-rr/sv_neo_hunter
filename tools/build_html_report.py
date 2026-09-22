#!/usr/bin/env python3
"""
build_html_report.py — a filter-by-filter validation report, as a web page.

WHAT THIS IS FOR
----------------
The markdown report says how many candidates survive each stage. This one says
WHY each stage exists, HOW it is measured, and WHY the threshold sits where it
does — beside the count of what it removed. That is what makes a result
checkable by someone who did not write the pipeline: a number without its
criterion cannot be argued with.

Every filter is rendered as: what it measures, how it is measured, the threshold
and its justification, and how many candidates entered, passed and were removed.

Where nothing reaches the RNA stage, the report says so explicitly and explains
which earlier filter emptied the set — a blank RNA section is otherwise
indistinguishable from an RNA stage that ran and found nothing.

    python tools/build_html_report.py --run-dir results/parental_PON10
    python tools/build_html_report.py --all
"""
from __future__ import annotations

import argparse
import html
import json
import os
import pathlib

import pandas as pd

# ============================================================================
# The rationale for each filter. Kept next to the code that counts it, so a
# threshold cannot drift away from its justification.
# ============================================================================

FILTERS = {
    "filter_pass": dict(
        name="Caller FILTER",
        measures="Whether the SV caller flagged the record as passing its own quality checks.",
        how="Read straight from the VCF <code>FILTER</code> column. Somatic call sets arrive "
            "already panel-filtered by the caller; germline sets are usually all-PASS by "
            "construction and are only <em>annotated</em> with a panel count, which is why "
            "the panel filter below has to be applied separately.",
        why="Records the caller itself rejected carry no usable evidence. Two classes matter: "
            "<code>PON</code> means the breakpoint was seen in a panel of normal samples, and "
            "<code>INFERRED</code> means the caller deduced the breakend from a copy-number "
            "transition with no read support at all — there is no junction sequence to build "
            "a peptide from.",
        threshold="Keep <code>FILTER == PASS</code>.",
        threshold_why="Not a tunable number: it is the caller's own verdict. The alternative "
                      "is to second-guess the caller with worse information than it had."),
    "paired": dict(
        name="Paired breakends",
        measures="Whether the breakend has a mate, i.e. whether both ends of the junction are known.",
        how="The ALT field of a mated breakend contains a bracket with the partner coordinate "
            "(<code>N[chr5:123456[</code>). Single breakends have none.",
        why="A junction peptide is built by joining two sequences. One coordinate is not enough: "
            "there is nothing to join it to. Single breakends can be real rearrangements, but "
            "they cannot yield a junction peptide, which is what this analysis is about.",
        threshold="Keep records whose ALT contains <code>[</code> or <code>]</code>.",
        threshold_why="Structural, not statistical. Dropping these is a limitation to state, "
                      "not a quality judgement — a single breakend may be a perfectly real SV."),
    "panel": dict(
        name="Panel of normals",
        measures="How many unrelated normal samples show the same breakpoint.",
        how="The caller annotates <code>PON_COUNT</code> against its own panel. An absent value "
            "means the breakpoint is not in the panel at all, which passes.",
        why="A breakpoint present in thousands of unrelated normals is either a recurrent "
            "artefact of the caller and reference, or a common germline polymorphism. Either "
            "way it is not private to this sample — and privacy is the property the whole "
            "analysis depends on. A recurrent site is the single strongest reason to disbelieve "
            "a candidate.",
        threshold="<code>PON_COUNT &lt; 10</code>.",
        threshold_why="This is an ABSOLUTE COUNT, not a frequency, so it only means something "
                      "against the panel size. The largest count observed in this panel is "
                      "~12,000, so &lt;10 means <em>seen in under 0.08% of normals</em> — strict. "
                      "The same number against a 50-sample panel would mean &lt;20%, which would "
                      "be permissive. The pipeline reports the count as a fraction of panel size "
                      "for exactly this reason, and runs an unfiltered branch alongside so the "
                      "filter's effect is visible rather than assumed."),
    "generated": dict(
        name="Produces a candidate peptide",
        measures="Whether the junction, once annotated against transcripts, yields any peptide at all.",
        how="The junction sequence is reconstructed, checked for frame, translated, and cut into "
            "8–11mers. Peptides already present in the wild-type protein of either transcript "
            "are subtracted.",
        why="Most junctions do not disrupt a coding transcript: they fall in intergenic space, "
            "in introns without changing the protein, or produce a sequence identical to the "
            "wild type. Those cannot generate a neoantigen by definition.",
        threshold="At least one peptide survives the wild-type subtraction.",
        threshold_why="No threshold to tune. This is the denominator that makes every later "
                      "count interpretable: a match rate quoted against admitted junctions "
                      "rather than against productive ones understates itself by an order of "
                      "magnitude."),
    "identical": dict(
        name="Identical to a catalogue peptide",
        measures="Whether a candidate peptide string appears verbatim in the patient catalogue.",
        how="Exact string intersection, case- and whitespace-normalised. No alignment, no "
            "similarity score.",
        why="This is the recurrence question itself. Exactness is deliberate: a near-match is a "
            "different peptide, and would be presented differently or not at all. Anything "
            "looser would answer a different question.",
        threshold="Exact match on the peptide sequence.",
        threshold_why="MHC binding prediction is deliberately NOT applied here. Whether the "
                      "sample's HLA would present the peptide has no bearing on whether the "
                      "peptide exists in it, and a binding filter would discard true sequence "
                      "matches. Presentability is a separate question, asked later and only of "
                      "survivors."),
    "concordant": dict(
        name="Gene concordance",
        measures="Whether the gene broken in this sample is the gene the catalogue attributes "
                 "the peptide to.",
        how="Compared against <strong>both</strong> breakends: a junction joins two genes and the "
            "catalogue may annotate the peptide to either side.",
        why="Short peptides can coincide between unrelated genomes by chance. Coincidental "
            "convergence would not preferentially land in the same gene, so concordance "
            "separates a shared locus from a shared accident.",
        threshold="The catalogue's gene matches either breakend's gene.",
        threshold_why="Testing only the first breakend loses real matches — 8 of 64 in one run, "
                      "all translocations whose catalogue gene sat on the second side. The cost "
                      "of the criterion is that a genuine convergent peptide arising in a "
                      "different gene is excluded by construction: sensitivity traded for "
                      "specificity, and worth stating whenever the count is quoted."),
    "complexity": dict(
        name="Sequence complexity",
        measures="Whether the peptide is repetitive enough to match by chance.",
        how="Four independent tests, any of which is disqualifying: Shannon entropy below 2 bits, "
            "one residue occupying half the peptide or more, a homopolymer run of 4 or more, or "
            "3 or fewer distinct residues.",
        why="Repetitive peptides such as <code>FFFFFFFFF</code> occur in unrelated proteomes by "
            "chance and are a classic prediction artefact. They match everything and mean nothing.",
        threshold="Entropy ≥ 2 bits AND max single residue &lt; 50% AND longest run &lt; 4 AND "
                  "&gt; 3 distinct residues.",
        threshold_why="Four criteria rather than one because low complexity takes several forms: "
                      "a homopolymer, a two-residue alternation and a low-entropy mixture all "
                      "match promiscuously but fail different single tests."),
    "self": dict(
        name="Not in the normal proteome",
        measures="Whether the peptide already exists in the healthy human proteome.",
        how="Exact substring search against the full Ensembl proteome, with sequences joined by "
            "a NUL separator so no peptide can match across a protein boundary.",
        why="A peptide present verbatim in normal protein is not neo. It is subject to central "
            "tolerance and would not be recognised as foreign, whatever the SV did.",
        threshold="The peptide is not a substring of any reference protein.",
        threshold_why="The separator matters: concatenating proteins without one invents junction "
                      "sequences that do not exist, and would wrongly mark real neopeptides as self."),
    "private": dict(
        name="Private to the sample",
        measures="Whether the event is specific to this sample or present in the general population.",
        how="Two independent filters, both of which must pass: the caller's panel of normals, "
            "and population allele frequency from gnomAD-SV by reciprocal overlap.",
        why="The panel and the population database disagree in both directions. Events with a "
            "panel count of 1 and 3 have carried population frequencies of 0.287 and 0.263 — "
            "common polymorphisms the panel missed. A well-called, well-expressed, "
            "read-supported deletion present in most of the population is real and still "
            "disqualified: what matters is not whether the call is correct but whether it is "
            "the sample's own.",
        threshold="Panel count &lt; 10 AND population frequency &lt; 0.001.",
        threshold_why="Frequency is reported for all nine gnomAD ancestry groups plus the "
                      "maximum across them, because which one is the right reference depends on "
                      "the donor's ancestry — knowledge the pipeline does not have. The automatic "
                      "verdict uses the maximum, the conservative choice: common anywhere is not "
                      "private. An unmatched event is treated as absent, and the run records that "
                      "it was unmatched rather than implying a clean result."),
    "confidence": dict(
        name="SV call confidence",
        measures="Whether the read evidence supports the junction being real.",
        how="Both breakends must independently satisfy mapping quality, fragment support and "
            "caller quality, plus a size floor for intra-chromosomal events.",
        why="Peptide credibility says nothing about whether the underlying SV exists. These are "
            "orthogonal questions and are kept apart: an event can be a flawless call and a "
            "common polymorphism, or a private event on a shaky call.",
        threshold="SEGMAPQ ≥ 30, VF ≥ 5, QUAL ≥ 20 at both breakends; intra-chromosomal size ≥ 100 bp.",
        threshold_why="Mapping quality 30 means under 1 in 1,000 chance of misplacement — below "
                      "that, repetitive sequence is the likelier explanation than a rearrangement. "
                      "The 100 bp floor exists because sub-100 bp intra-chromosomal calls are "
                      "frequently alignment artefacts around indels and homopolymers. This is "
                      "REPORTED, NOT APPLIED as a filter: a non-confident event with orthogonal "
                      "RNA evidence is more believable than a confident one with none."),
    "rna": dict(
        name="Transcribed across the junction",
        measures="Whether RNA reads actually cross the breakpoint.",
        how="The test is chosen by the event's geometry: supplementary alignments landing on the "
            "partner breakend for inter-chromosomal junctions, a CIGAR gap of matching size for "
            "deletions, or the inserted sequence itself for insertion-driven events. Support is "
            "counted per fragment, not per alignment record, and reads placed ambiguously do not "
            "count.",
        why="A genomic event that is never transcribed cannot produce a peptide. This is the "
            "only stage that observes the consequence rather than the cause.",
        threshold="≥ 5 junction-crossing fragments at mapping quality ≥ 20.",
        threshold_why="Five matches the default read-support cutoff of established fusion-"
                      "neoantigen tools, so counts are comparable with the field; 3–4 is reported "
                      "as a lower tier rather than discarded. Coverage never promotes a tier: a "
                      "breakpoint inside an expressed gene has thousands of reads whether or not "
                      "the junction exists. A soft-clip shows a read <em>ends</em> at the "
                      "breakpoint, not that it crosses, so soft-clips are diagnostic only. "
                      "Counting alignment records instead of fragments inflated one event's "
                      "support from 3 to 7."),
}


def esc(text) -> str:
    return html.escape(str(text))


def bar(entered: int, passed: int) -> str:
    """A proportional bar: what passed, what was removed."""
    if not entered:
        return ""
    kept = 100.0 * passed / entered
    return (f'<div class="bar" role="img" aria-label="{passed} of {entered} passed">'
            f'<span class="kept" style="width:{kept:.2f}%"></span></div>')


def grain_change(from_unit: str, to_unit: str, explanation: str) -> str:
    """A visible marker where the counting unit changes.

    Without it the cascade reads as if numbers jump arbitrarily — 2,286 breakends
    followed by 1,143 junctions looks like an error. Every grain change is a
    place where a reader can be misled, so each one is announced.
    """
    return f"""
<div class="grain">
  <span class="from">{esc(from_unit)}</span>
  <span class="arrow">&rarr;</span>
  <span class="to">{esc(to_unit)}</span>
  <p>{explanation}</p>
</div>"""


def losses(caption: str, rows: list[dict], columns: list[str]) -> str:
    """Name what a filter removed, with the value that removed it.

    A count alone is not auditable. "Panel of normals removed 12,054" cannot be
    checked, argued with, or reconciled against a previous analysis; naming the
    event and the value that removed it can. This block exists because a reader who
    remembers a gene from an earlier run has no other way to find out where it
    went, and its silent absence reads as a bug rather than as a filter working.
    """
    if not rows:
        return ""
    head = "".join(f"<th>{esc(c)}</th>" for c in columns)
    body = "".join(
        "<tr>" + "".join(f"<td>{esc(str(r.get(c, '')))}</td>" for c in columns) + "</tr>"
        for r in rows)
    return f"""
  <details class="losses">
    <summary>{caption}</summary>
    <div class="scroll"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>
  </details>"""


def unfiltered_losses(run: pathlib.Path) -> list[dict]:
    """Catalogue-matching events the panel filter removed, recovered from the
    parallel unfiltered run.

    The panel filter acts BEFORE peptide generation, so a filtered run holds no
    record of what it discarded — it never annotated those junctions and never
    learned whether they matched the catalogue. The unfiltered branch exists
    precisely so the filter's cost is measurable instead of assumed, and this is
    where that measurement is read back.
    """
    sibling = run.parent / f"{run.name.split('_')[0]}_noPON"
    target = sibling / "master_sv.tsv"
    if sibling == run or not target.exists():
        return []
    frame = pd.read_csv(target, sep="\t", low_memory=False)
    if "n_peptides_matched" not in frame.columns:
        return []
    frame = frame[pd.to_numeric(frame.n_peptides_matched, errors="coerce") > 0]

    kept = set()
    mine = run / "master_sv.tsv"
    if mine.exists():
        kept = set(pd.read_csv(mine, sep="\t", low_memory=False)
                   .get("sv_id", pd.Series(dtype=str)).astype(str))

    out = []
    for _, row in frame.iterrows():
        if str(row.get("sv_id")) in kept:
            continue
        gene = row.get("gene1") or ""
        if row.get("gene2") and row.get("gene2") != gene:
            gene = f"{gene}–{row.get('gene2')}"
        reads = row.get("junction_reads")
        out.append({
            "gene": _text(gene), "svtype": _text(row.get("svtype")),
            "PON_COUNT": _num(row.get("pon_count")) or "—",
            "gnomAD popmax": _num(row.get("gnomad_af_popmax"), 3) or "—",
            "population": _text(row.get("gnomad_af_popmax_pop")),
            "RNA reads": _num(reads) if not pd.isna(reads) else "—",
            "RNA tier": _text(row.get("rna_tier")),
        })
    def rna_reads(row) -> float:
        try:
            return float(row["RNA reads"].replace(",", ""))
        except (ValueError, AttributeError):
            return -1.0
    out.sort(key=lambda r: -rna_reads(r))
    return out


def _num(value, decimals: int = 0) -> str:
    """Format a number for display; empty string when absent."""
    if value is None or (isinstance(value, float) and pd.isna(value)) or value == "":
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    return f"{number:,.{decimals}f}" if decimals else f"{number:,.0f}"


def _truth(frame, column: str):
    """A boolean Series from a column that may be text, missing, or absent."""
    if not len(frame) or column not in frame.columns:
        return pd.Series([False] * len(frame), index=frame.index)
    return frame[column].astype(str).str.lower().isin(("true", "1"))


def _pon_votes(events) -> bool | None:
    """Did the panel count participate in `is_private` for this run?

    Written by `annotate_privacy()`. None for a run that predates the column, so
    a report never asserts which basis was used when it cannot know.
    """
    if not len(events) or "pon_in_privacy" not in events.columns:
        return None
    return bool(_truth(events, "pon_in_privacy").iloc[0])


def _text(value) -> str:
    """A cell's text, with absence rendered as an em dash.

    NEVER write `value or ""` for these: pandas yields float('nan') for a missing
    cell and NaN is TRUTHY, so the fallback never fires and the literal string
    "nan" reaches the page. That mistake has shipped in this codebase repeatedly.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "—"
    text = str(value).strip()
    return "—" if text.lower() in ("nan", "none", "") else text


def _zero(value) -> bool:
    """True only for a real, measured zero — never for a missing value.

    An absent coverage figure means the test did not run; reporting that as
    "not transcribed" is the confusion this whole report exists to prevent.
    """
    if value is None or value == "" or (isinstance(value, float) and pd.isna(value)):
        return False
    try:
        return float(value) == 0
    except (TypeError, ValueError):
        return False


def named_losses(frame, column: str, keep_true: bool, caption: str,
                 extra: list[str] | None = None) -> str:
    """Name the rows a boolean criterion removed, with their deciding values."""
    if not len(frame) or column not in frame.columns:
        return ""
    truth = frame[column].astype(str).str.lower().isin(("true", "1"))
    dropped = frame[~truth] if keep_true else frame[truth]
    if not len(dropped):
        return ""
    columns = ["gene"] + (extra or [])
    rows = []
    for _, row in dropped.iterrows():
        gene = row.get("gene1")
        entry = {"gene": _text(gene if _text(gene) != "—" else row.get("gene"))}
        for name in (extra or []):
            value = row.get(name)
            if "af" in name:
                entry[name] = _num(value, 3) or "—"
            elif name in ("pon_count", "event_size"):
                entry[name] = _num(value) or "—"
            else:
                entry[name] = _text(value)
        rows.append(entry)
    return losses(caption, rows[:40], columns)


def panel_reported(admitted_pon, pon_voted, kind: str) -> str:
    """The panel filter, deliberately not applied.

    This block occupies the place a filter card would have. Without it the
    cascade simply lacks a panel step, which reads as "the panel removed nothing"
    rather than "the panel was reported and not enforced" — the opposite
    conclusion from the same numbers.
    """
    routes = []
    if kind == "germline":
        routes.append("this pipeline's own <code>PON_COUNT</code> threshold "
                      "(<code>PON_MAX</code>) is switched off")
    if admitted_pon is not None:
        routes.append(f"<b>{admitted_pon:,}</b> record(s) the caller had already "
                      f"rejected as <code>FILTER=PON</code> were admitted and "
                      f"annotated rather than dropped")
    if pon_voted is False:
        routes.append("the panel count does not vote in <code>is_private</code>, "
                      "which rests on population frequency alone")
    return f"""
<section class="filter reported">
  <header>
    <h3>Panel of normals — reported, not applied</h3>
    <div class="counts"><span class="n"><b>0</b><small>removed</small></span></div>
  </header>
  <dl>
    <dt>What this branch does</dt>
    <dd>The panel count is measured and carried into every table, but it removes
    nothing: {"; ".join(routes)}.</dd>
    <dt>Why a branch like this exists</dt>
    <dd>The panel threshold is the single most consequential filter in the
    cascade, and its cost cannot be read off a filtered run — the discarded
    records are never annotated, so nothing records what they would have
    matched. Running the same sample with the panel reported instead of enforced
    makes that cost measurable rather than assumed.</dd>
    <dt>How to read the tables that follow</dt>
    <dd>With <code>pon_count</code>, <code>pon_fraction</code> and the nine
    <code>gnomad_af_*</code> columns in view. An event surviving here is
    <em>not</em> a candidate neoantigen by itself: a junction can be real,
    expressed and confidently called while being a germline polymorphism carried
    by most of the population. Measured on this lineage, the events with the
    strongest RNA support are exactly the common ones.</dd>
    <dt>Why not just lower the threshold</dt>
    <dd>Because the panel reaches a candidate twice — once at admission and once
    inside <code>is_private</code> — and on a somatic call set it arrives as a
    <code>FILTER</code> value the caller already applied. Relaxing only the
    numeric threshold changes nothing there: it would produce output identical to
    the filtered branch under an unfiltered label.</dd>
  </dl>
</section>"""


def filter_card(key: str, entered, passed, note: str = "", removed_detail: str = "") -> str:
    spec = FILTERS[key]
    removed = None if entered is None or passed is None else entered - passed
    counts = (f'<div class="counts">'
              f'<span class="n"><b>{entered:,}</b><small>entered</small></span>'
              f'<span class="n pass"><b>{passed:,}</b><small>passed</small></span>'
              f'<span class="n drop"><b>{removed:,}</b><small>removed</small></span>'
              f'</div>{bar(entered, passed)}'
              if removed is not None else
              '<div class="counts"><span class="n"><b>NA</b><small>not evaluated</small></span></div>')
    return f"""
<section class="filter">
  <header>
    <h3>{esc(spec['name'])}</h3>
    {counts}
  </header>
  <dl>
    <dt>What it measures</dt><dd>{spec['measures']}</dd>
    <dt>How it is measured</dt><dd>{spec['how']}</dd>
    <dt>Why this filter exists</dt><dd>{spec['why']}</dd>
    <dt>Threshold</dt><dd class="thr">{spec['threshold']}</dd>
    <dt>Why this threshold</dt><dd>{spec['threshold_why']}</dd>
  </dl>
  {f'<p class="note">{note}</p>' if note else ''}
  {removed_detail}
</section>"""


def build(run_dir: str) -> str:
    run = pathlib.Path(run_dir)
    summary = json.loads((run / "summary.json").read_text())
    counts = summary.get("counts", {})
    s1 = summary.get("stage1_funnel", {})
    sample, branch = summary.get("sample"), summary.get("branch", "")
    kind = summary.get("vcf_kind", "")

    sv = pd.read_csv(run / "master_sv.tsv", sep="\t", low_memory=False) \
        if (run / "master_sv.tsv").exists() else pd.DataFrame()
    matches = pd.read_csv(run / "stage3_matches.tsv", sep="\t", low_memory=False) \
        if (run / "stage3_matches.tsv").exists() else pd.DataFrame()
    events = pd.read_csv(run / "credible_events.tsv", sep="\t") \
        if (run / "credible_events.tsv").exists() else pd.DataFrame()
    rna = pd.read_csv(run / "stage8_rna_evidence.tsv", sep="\t") \
        if (run / "stage8_rna_evidence.tsv").exists() else pd.DataFrame()

    def flag(frame, col):
        return int(frame[col].astype(str).str.lower().isin(("true", "1")).sum()) \
            if (len(frame) and col in frame.columns) else 0

    cards = []
    records = s1.get("records", 0)
    cards.append(filter_card("filter_pass", records, s1.get("after_filter_pass", records),
                             "This is a germline call set: every record passes FILTER by "
                             "construction, so the panel filter below does the work."
                             if kind == "germline" else ""))

    # A branch that reports the panel instead of filtering on it has to say so
    # HERE, where the filter would otherwise have been, or its absence is
    # indistinguishable from a filter that removed nothing.
    admitted_pon = s1.get("caller_pon_admitted")
    pon_voted_here = _pon_votes(events) if len(events) else _pon_votes(sv)
    if admitted_pon is not None or pon_voted_here is False:
        cards.append(panel_reported(admitted_pon, pon_voted_here, kind))
    cards.append(filter_card("paired", s1.get("after_filter_pass", records),
                             s1.get("after_paired_only", 0)))
    pon_key = next((k for k in s1 if k.startswith("after_pon")), None)
    if pon_key:
        lost = unfiltered_losses(run)
        detail = losses(
            f"What this filter cost: {len(lost)} junction(s) that DID match the "
            f"catalogue in the unfiltered branch", lost,
            ["gene", "svtype", "PON_COUNT", "gnomAD popmax", "population",
             "RNA reads", "RNA tier"]) if lost else ""
        cards.append(filter_card(
            "panel", s1.get("after_paired_only", 0), s1[pon_key],
            "This filter acts before annotation, so this run holds no record of "
            "what it discarded. The names below are recovered from the parallel "
            "unfiltered run — without it, the cost of this threshold would be "
            "invisible rather than merely large."
            if lost else "", detail))
    junctions = counts.get("junctions", 0)
    cards.append(grain_change(
        f"{s1.get(pon_key, 0):,} breakends" if pon_key else "breakends",
        f"{junctions:,} junctions",
        "Each junction has two breakends, so the count roughly halves. Records "
        "whose mate did not survive the filters above are dropped here rather "
        "than paired against a rejected partner."))
    productive = int((pd.to_numeric(sv.n_peptides_generated, errors="coerce") > 0).sum()) \
        if "n_peptides_generated" in sv.columns else 0
    cards.append(filter_card("generated", junctions, productive))

    matched_sv = int((pd.to_numeric(sv.n_peptides_matched, errors="coerce") > 0).sum()) \
        if "n_peptides_matched" in sv.columns else 0
    cards.append(filter_card("identical", productive, matched_sv,
                             f"At peptide grain this is "
                             f"<b>{counts.get('matches_unique_peptides', 0)} distinct peptides</b> "
                             f"across {matched_sv} junctions — one junction yields many "
                             f"overlapping window peptides, which is why findings are counted "
                             f"as events."))

    n_matches = len(matches)
    if n_matches:
        cards.append(grain_change(
            f"{matched_sv:,} junctions", f"{n_matches:,} peptide rows",
            "The next three filters judge peptide sequences, not junctions, so "
            "the unit changes. One junction contributes several peptides — the "
            "sliding window produces many overlapping ones — and each is tested "
            "on its own."))
        cards.append(filter_card(
            "concordant", n_matches, flag(matches, "gene_concordant"), "",
            named_losses(matches, "gene_concordant", True,
                         "Peptides removed as gene-discordant",
                         ["peptide", "catalogue_gene"])))
        after_conc = flag(matches, "gene_concordant")
        lc = flag(matches, "low_complexity")
        concordant = matches[matches.gene_concordant.astype(str).str.lower()
                             .isin(("true", "1"))] if "gene_concordant" in matches else matches
        cards.append(filter_card(
            "complexity", after_conc, after_conc - lc, "",
            named_losses(concordant, "low_complexity", False,
                         "Peptides removed as low-complexity", ["peptide"])))
        cards.append(filter_card("self", after_conc - lc,
                                 after_conc - lc - flag(matches, "is_self")))

    n_events = counts.get("events", 0)
    if n_events:
        cards.append(grain_change(
            f"{flag(matches, 'credible'):,} credible peptide rows",
            f"{n_events:,} genomic events",
            "Back to genomic grain, and the most consequential collapse in the "
            "pipeline: many peptides from one junction are ONE finding, not many. "
            "Reporting the peptide count as if it were a finding count has "
            "overstated evidence twenty-fold elsewhere."))
        cards.append(filter_card(
            "private", n_events, flag(events, "is_private"), "",
            named_losses(events, "is_private", True,
                         "Events removed as non-private",
                         ["pon_count", "gnomad_af_popmax", "gnomad_af_popmax_pop"])))
        cards.append(filter_card(
            "confidence", n_events, flag(events, "sv_hc"),
            "Reported, not applied — no event is dropped here.",
            named_losses(events, "sv_hc", True,
                         "Events below the confidence bar (retained, but flagged)",
                         ["svtype", "event_size"])))
        supported = int(rna.rna_tier.isin(["STRONG", "SUGGESTIVE", "WEAK"]).sum()) \
            if "rna_tier" in rna.columns else 0
        tests = rna.test.value_counts().to_dict() if "test" in rna.columns else {}
        note = ("Tests applied by geometry: "
                + ", ".join(f"<b>{v}</b> {k}" for k, v in tests.items()) + ".") if tests else ""
        silent = []
        if "rna_tier" in rna.columns:
            for _, row in rna[~rna.rna_tier.isin(["STRONG", "SUGGESTIVE", "WEAK"])].iterrows():
                silent.append({
                    "gene": row.get("gene1", "") or row.get("gene", ""),
                    "test": row.get("test", ""),
                    "junction reads": _num(row.get("junction_reads")),
                    "min coverage": _num(row.get("min_coverage")),
                    "why": "test not applicable to this geometry"
                           if str(row.get("test", "")) == "none"
                           else ("locus not transcribed"
                                 if _zero(row.get("min_coverage"))
                                 else "expressed locus, but nothing crosses the junction"),
                })
        cards.append(filter_card(
            "rna", n_events, supported, note,
            losses("Events with no RNA support, and the reason for each", silent,
                   ["gene", "test", "junction reads", "min coverage", "why"])))

    # ---------------------------------------------------------------- verdict
    # `events_hc_and_rna` does NOT include privacy — the name says so, but an
    # earlier version of this report described it as "satisfies every criterion",
    # which overstated a common germline event as a survivor.
    # The all-criteria count is its own figure; fall back to recomputing it when
    # reading a run that predates it, rather than silently substituting.
    hc_and_rna = counts.get("events_hc_and_rna", 0)
    both = counts.get("events_private_hc_and_rna")
    if both is None:
        both = int((_truth(events, "is_private") & _truth(events, "sv_hc")
                    & events.sample_sv_id.isin(
                        rna[rna.rna_tier.isin(["STRONG", "SUGGESTIVE"])].sample_sv_id)).sum()) \
            if len(events) and len(rna) and "rna_tier" in rna.columns else 0
    if not n_matches:
        emptied = ("no candidate peptide matched the catalogue" if productive
                   else "no admitted junction produced a candidate peptide at all")
        verdict = f"""
<section class="verdict empty">
  <h2>Nothing reaches the RNA stage — and why that is not a null result</h2>
  <p>The chain stops before RNA because <b>{esc(emptied)}</b>. The RNA stage was
  therefore never run: there was nothing to test. That is a different statement
  from "the junctions are not transcribed", and the tables record it as
  <code>NA</code>, not as zero.</p>
  <p>Of {junctions:,} admitted junctions, {productive:,} produced candidate peptides
  ({100.0 * productive / junctions if junctions else 0:.1f}%), and none of those
  peptides appears in the reference catalogue. The sample's SVs are real and may
  well be transcribed; what is absent is any overlap with the patient catalogue,
  which is the question this analysis asks.</p>
  <p class="caution">A small candidate universe makes raw counts incomparable
  with samples that have orders of magnitude more. Compare the per-1,000 rate and
  the permutation null instead.</p>
</section>"""
    elif not both:
        reasons = []
        if n_events:
            not_private = n_events - flag(events, "is_private")
            not_hc = n_events - flag(events, "sv_hc")
            supported = int(rna.rna_tier.isin(["STRONG", "SUGGESTIVE", "WEAK"]).sum()) \
                if "rna_tier" in rna.columns else 0
            if not_private:
                reasons.append(f"<b>{not_private}</b> are common in the population or the "
                               f"panel — real events, but not this sample's own")
            if not_hc:
                reasons.append(f"<b>{not_hc}</b> do not meet the SV-confidence bar")
            reasons.append(f"<b>{n_events - supported}</b> have no read crossing the junction "
                           f"in RNA")
        verdict = f"""
<section class="verdict">
  <h2>No candidate satisfies every criterion</h2>
  <p>{n_events} credible genomic events enter the final stages, and none is
  simultaneously private, confidently called and transcribed:</p>
  <ul>{''.join(f'<li>{r}</li>' for r in reasons)}</ul>
  <p>These are overlapping sets, not a sequence: an event can fail more than one.
  The per-event table below shows which.</p>
  <p class="caution">This is not a failed analysis. Reducing a four-figure
  candidate count to a single-digit believable set is the expected outcome — a
  28-team benchmark found roughly 6% of top-ranked neoantigen predictions
  validate functionally. What matters is that each removal is attributable to a
  stated criterion.</p>
</section>"""
    else:
        pon_voted = _pon_votes(events)
        privacy_basis = ("population frequency alone — the panel count is "
                         "reported in the table but does not vote in this branch"
                         if pon_voted is False else
                         "both the panel count and population frequency")
        verdict = f"""
<section class="verdict survives">
  <h2>{both} candidate(s) satisfy every criterion</h2>
  <p>Private, confidently called, <em>and</em> transcribed across the junction.
  Privacy here rests on {privacy_basis}. See the per-event table below.</p>
  <p class="note">A further <b>{hc_and_rna}</b> event(s) are confidently called
  and transcribed but <b>not</b> private. Those two counts must not be conflated:
  a real, expressed junction that is common in the population is a polymorphism,
  not a recurrent tumour neoantigen.</p>
</section>"""

    # ---------------------------------------------------------------- events
    table_html = ""
    if n_events and len(events):
        merged = events.copy()
        if len(rna):
            merged = merged.merge(rna[[c for c in rna.columns if c != "gene"]],
                                  on="sample_sv_id", how="left")
        cols = [("gene", "Gene"), ("svtype", "Type"), ("event_size", "Size bp"),
                ("insert_len", "Insert"), ("pon_count", "Panel"),
                ("gnomad_af_popmax", "gnomAD max"), ("gnomad_af_popmax_pop", "in"),
                ("sv_hc", "Confident"), ("is_private", "Private"),
                ("test", "RNA test"), ("junction_reads", "Crossing reads"),
                ("rna_tier", "Tier")]
        head = "".join(f"<th>{esc(label)}</th>" for _, label in cols if _ in merged.columns)
        rows = []
        for _, row in merged.iterrows():
            cells = []
            for key, _label in cols:
                if key not in merged.columns:
                    continue
                value = row[key]
                text = "—" if pd.isna(value) else (
                    f"{value:,.0f}" if isinstance(value, float) and float(value).is_integer()
                    else f"{value:.3g}" if isinstance(value, float) else esc(value))
                css = ""
                if key in ("sv_hc", "is_private"):
                    css = " class=\"yes\"" if str(value).lower() in ("true", "1") else " class=\"no\""
                if key == "rna_tier" and str(value) in ("STRONG", "SUGGESTIVE", "WEAK"):
                    css = " class=\"yes\""
                cells.append(f"<td{css}>{text}</td>")
            rows.append("<tr>" + "".join(cells) + "</tr>")
        table_html = f"""
<section>
  <h2>The surviving events, one by one</h2>
  <p>Every credible event with the evidence behind each verdict. Values, not
  judgements: the thresholds above say how to read them.</p>
  <div class="scroll"><table><thead><tr>{head}</tr></thead>
  <tbody>{''.join(rows)}</tbody></table></div>
</section>"""

    generator = summary.get("peptide_generator", {})
    criteria = summary.get("criteria", {})

    return f"""<title>{esc(sample)} — filter-by-filter validation</title>
<style>
:root {{
  --ground:#F6F8FA; --surface:#FFFFFF; --surface-alt:#EDF1F5; --ink:#141A21;
  --ink-soft:#55606E; --ink-faint:#8A94A2; --rule:#D8DFE7; --steel:#215D8C;
  --steel-soft:#DDE9F3; --amber:#B26B12; --amber-soft:#F7EEDF; --garnet:#8F3A3A;
  --garnet-soft:#F5E4E4; --moss:#2C6A4E;
  --serif:"Iowan Old Style","Palatino Linotype",Palatino,"Book Antiqua",Georgia,serif;
  --mono:ui-monospace,"SF Mono",SFMono-Regular,Menlo,Consolas,monospace;
}}
@media (prefers-color-scheme:dark) {{ :root:not([data-theme="light"]) {{
  --ground:#0E1319; --surface:#161D25; --surface-alt:#1D2731; --ink:#E6EBF1;
  --ink-soft:#A3AEBC; --ink-faint:#6D7887; --rule:#2A3540; --steel:#6FA8D4;
  --steel-soft:#17303F; --amber:#D9A05B; --amber-soft:#33260F; --garnet:#D08585;
  --garnet-soft:#3A2020; --moss:#74B695;
}} }}
:root[data-theme="dark"] {{
  --ground:#0E1319; --surface:#161D25; --surface-alt:#1D2731; --ink:#E6EBF1;
  --ink-soft:#A3AEBC; --ink-faint:#6D7887; --rule:#2A3540; --steel:#6FA8D4;
  --steel-soft:#17303F; --amber:#D9A05B; --amber-soft:#33260F; --garnet:#D08585;
  --garnet-soft:#3A2020; --moss:#74B695;
}}
*{{box-sizing:border-box}}
body{{background:var(--ground);color:var(--ink);font-family:var(--serif);
 font-size:17px;line-height:1.6;margin:0;padding:0 1.25rem 5rem}}
.wrap{{max-width:52rem;margin:0 auto}}
header.page{{padding:3.5rem 0 2rem;border-bottom:2px solid var(--ink)}}
.eyebrow{{font-family:var(--mono);font-size:.72rem;letter-spacing:.16em;
 text-transform:uppercase;color:var(--steel);margin:0 0 .75rem}}
h1{{font-size:clamp(1.9rem,4.5vw,2.6rem);line-height:1.15;margin:0 0 .75rem;
 font-weight:600;letter-spacing:-.015em;text-wrap:balance}}
h2{{font-size:1.4rem;margin:3rem 0 .75rem;font-weight:600;text-wrap:balance}}
h3{{font-size:1.1rem;margin:0;font-weight:600}}
p{{margin:0 0 .9rem}} .lede{{font-size:1.1rem;color:var(--ink-soft);margin:0}}
code{{font-family:var(--mono);font-size:.85em;background:var(--surface-alt);
 padding:.1em .35em;border-radius:2px}}
.filter{{background:var(--surface);border:1px solid var(--rule);border-radius:3px;
 margin:1rem 0;overflow:hidden}}
.filter>header{{display:flex;flex-wrap:wrap;gap:1rem;align-items:center;
 justify-content:space-between;padding:.9rem 1.2rem;background:var(--surface-alt);
 border-bottom:1px solid var(--rule)}}
.counts{{display:flex;gap:1.4rem;font-family:var(--mono);
 font-variant-numeric:tabular-nums}}
.n{{display:flex;flex-direction:column;align-items:flex-end;line-height:1.15}}
.n b{{font-size:1.15rem}} .n small{{font-size:.68rem;color:var(--ink-faint);
 text-transform:uppercase;letter-spacing:.08em}}
.n.pass b{{color:var(--moss)}} .n.drop b{{color:var(--garnet)}}
.bar{{width:100%;height:4px;background:var(--garnet-soft);flex-basis:100%}}
.bar .kept{{display:block;height:100%;background:var(--moss)}}
dl{{margin:0;padding:1.1rem 1.2rem;display:grid;grid-template-columns:auto 1fr;
 gap:.5rem 1.1rem}}
dt{{font-family:var(--mono);font-size:.68rem;letter-spacing:.09em;
 text-transform:uppercase;color:var(--ink-faint);padding-top:.28rem;white-space:nowrap}}
dd{{margin:0;color:var(--ink-soft);font-size:.96rem}}
dd.thr{{font-family:var(--mono);font-size:.85rem;color:var(--amber)}}
.note{{margin:0;padding:.75rem 1.2rem;background:var(--steel-soft);
 font-size:.9rem;border-top:1px solid var(--rule)}}
.filter.reported{{border-left:3px solid var(--moss)}}
.filter.reported h3{{color:var(--moss)}}
.losses{{margin-top:1rem;border-top:1px solid var(--rule);padding-top:.8rem}}
.losses summary{{cursor:pointer;font-size:.85rem;font-weight:600;color:var(--steel);
 letter-spacing:.01em}}
.losses summary::marker{{color:var(--ink-faint)}}
.losses .scroll{{overflow-x:auto;margin-top:.7rem}}
.losses table{{border-collapse:collapse;width:100%;font-size:.82rem;
 font-family:var(--mono);font-variant-numeric:tabular-nums}}
.losses th{{text-align:left;font-weight:600;color:var(--ink-soft);
 border-bottom:1px solid var(--rule);padding:.35rem .7rem .35rem 0;white-space:nowrap}}
.losses td{{padding:.3rem .7rem .3rem 0;border-bottom:1px solid var(--rule);
 white-space:nowrap}}
.grain{{display:flex;flex-wrap:wrap;align-items:baseline;gap:.6rem;
 margin:1.6rem 0 .6rem;padding:.7rem 1rem;border-left:3px solid var(--amber);
 background:var(--amber-soft);border-radius:0 3px 3px 0}}
.grain .from,.grain .to{{font-family:var(--mono);font-size:.82rem;font-weight:600}}
.grain .arrow{{color:var(--amber);font-size:1.1rem}}
.grain p{{flex-basis:100%;margin:.35rem 0 0;font-size:.88rem;color:var(--ink-soft)}}
.verdict{{background:var(--surface);border:2px solid var(--amber);border-radius:3px;
 padding:1.4rem 1.5rem;margin:2.5rem 0}}
.verdict.empty{{border-color:var(--ink-faint)}}
.verdict.survives{{border-color:var(--moss)}}
.verdict h2{{margin-top:0}}
.verdict ul{{margin:.5rem 0 1rem;padding-left:1.2rem;color:var(--ink-soft)}}
.caution{{font-size:.92rem;color:var(--ink-soft);border-left:3px solid var(--steel);
 padding-left:.9rem;margin-bottom:0}}
.scroll{{overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-size:.86rem;
 font-variant-numeric:tabular-nums}}
th,td{{text-align:left;padding:.45rem .6rem;border-bottom:1px solid var(--rule);
 white-space:nowrap}}
th{{font-family:var(--mono);font-size:.66rem;letter-spacing:.08em;
 text-transform:uppercase;color:var(--ink-soft);border-bottom:1px solid var(--ink)}}
td.yes{{color:var(--moss);font-weight:600}} td.no{{color:var(--ink-faint)}}
footer{{margin-top:3.5rem;padding-top:1.25rem;border-top:1px solid var(--rule);
 font-size:.85rem;color:var(--ink-faint)}}
</style>
<div class="wrap">
<header class="page">
  <p class="eyebrow">svneo · filter-by-filter validation</p>
  <h1>{esc(sample)}</h1>
  <p class="lede">{esc(kind)} call set{f', branch {esc(branch)}' if branch else ''} ·
  {records:,} SV records in · {n_events} credible events out ·
  <b>{both if isinstance(both, int) else 'NA'}</b> satisfying every criterion.</p>
</header>

<h2>Every filter, and why</h2>
<p>Each stage below states what it measures, how, and why its threshold sits
where it does — beside the count of what it removed. A number without its
criterion cannot be checked.</p>
{''.join(cards)}
{verdict}
{table_html}
<footer>
  Generated from <code>{esc(run.name)}</code> by <code>tools/build_html_report.py</code>.
  Peptides from <code>{esc(generator.get('generator'))}</code>
  ({esc(generator.get('version'))}). The full threshold manifest
  ({len([k for k in criteria if not k.startswith('_')])} criteria) is in
  <code>summary.json</code> beside these outputs.
  Read evidence supports the claim that a junction is <em>transcribed</em>, never
  that a peptide is <em>presented</em>. An unresolved strand skew
  means no result here is strand-controlled.
</footer>
</div>"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--results-dir", default="results")
    args = ap.parse_args()

    runs = [args.run_dir] if args.run_dir else (
        [os.path.join(args.results_dir, d) for d in sorted(os.listdir(args.results_dir))
         if os.path.exists(os.path.join(args.results_dir, d, "summary.json"))]
        if args.all else [])
    if not runs:
        ap.error("give --run-dir or --all")

    out_dir = os.path.join(args.results_dir, "reports")
    os.makedirs(out_dir, exist_ok=True)
    for run in runs:
        path = os.path.join(out_dir, f"{os.path.basename(run)}.html")
        with open(path, "w") as handle:
            handle.write(build(run))
        print(f"  wrote {path}")


if __name__ == "__main__":
    main()
