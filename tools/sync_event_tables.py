#!/usr/bin/env python3
"""
sync_event_tables.py — keep the events tables in `docs/` generated, not typed.

A results table transcribed by hand into prose is the most reliable way for a
document to stop describing its own run: nothing errors, and the reader has no
way to tell. This fills marked regions of a markdown file from the run outputs,
so the tables in the documentation are machine-maintained while the surrounding
argument stays hand-written.

Mark a region in the document:

    <!-- events:parental_noPON -->
    ...anything here is replaced...
    <!-- /events:parental_noPON -->

then run this. The markdown stays readable on its own — the table is really in
the file, not a placeholder resolved at render time — and `--check` verifies the
file is current without writing, so it can gate a commit.

    python tools/sync_event_tables.py
    python tools/sync_event_tables.py --check
"""
from __future__ import annotations

import argparse
import os
import pathlib
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

REPO = pathlib.Path(__file__).resolve().parent.parent
#: Documents to keep in sync. Passed on the command line, because which
#: write-ups exist is a property of a given analysis, not of this tool.
DOCS: list[str] = []

#: Column key in the results, and its heading in the document.
EVENT_COLUMNS = [
    ("gene", "Gene"), ("svtype", "Type"), ("event_size", "Size"),
    ("n_peptides", "Peptides"), ("pon_count", "PON"),
    ("gnomad_af_nfe", "gnomAD nfe"), ("gnomad_af_amr", "gnomAD amr"),
    ("gnomad_af_popmax", "gnomAD max"), ("is_private", "Not a common variant"),
    ("sv_hc", "HC"), ("test", "RNA test"), ("min cov", "min cov"),
    ("junction_reads", "Crossing"), ("rna_tier", "Tier"),
]
#: `min cov` above reads from `min_coverage`; kept separate so the heading can
#: differ from the column name without a second mapping.
SOURCE = {"min cov": "min_coverage"}

COUNTED = {"event_size", "n_peptides", "pon_count", "min_coverage",
           "junction_reads"}

#: Per-group frequencies shown beside the maximum. The maximum alone is the
#: right default for a FILTER (it is the conservative bound) and the wrong thing
#: to read as a description of a given sample: it reports whichever of the nine
#: ancestry groups happens to be highest, which may be unrelated to the donor.
#: All nine are in `master_sv.tsv` either way.


def cell(key: str, value, population=None) -> str:
    """One cell, with absence as an em dash.

    `value or "—"` does NOT work here: pandas yields float('nan') for a missing
    cell, and NaN is truthy, so the fallback never fires and the literal string
    "nan" reaches the document.
    """
    source = SOURCE.get(key, key)
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "—"
    text = str(value).strip()
    if text.lower() in ("nan", "none", ""):
        return "—"
    if source in COUNTED:
        return f"{float(value):,.0f}"
    if source.startswith("gnomad_af"):
        rendered = f"{float(value):.3f}"
        # The maximum carries the group it came from; a bare maximum invites the
        # reader to apply a frequency from an ancestry unrelated to the sample.
        # The named per-group columns need no such label.
        if source == "gnomad_af_popmax" and population:
            return f"{rendered} ({population})"
        return rendered
    if source in ("is_private", "sv_hc"):
        yes = text.lower() in ("true", "1")
        if source == "is_private":
            return "yes" if yes else "**no**"
        return "yes" if yes else "no"
    return text


def events_table(run: str) -> str:
    """The credible events of one run as a markdown table, from its own outputs."""
    directory = REPO / "results" / run
    events_path = directory / "credible_events.tsv"
    if not events_path.exists():
        return f"*No `credible_events.tsv` for `{run}`.*"

    frame = pd.read_csv(events_path, sep="\t", low_memory=False)
    rna_path = directory / "stage8_rna_evidence.tsv"
    if rna_path.exists():
        rna = pd.read_csv(rna_path, sep="\t", low_memory=False)
        wanted = [c for c in ("sample_sv_id", "rna_tier", "junction_reads",
                              "test", "min_coverage") if c in rna.columns]
        frame = frame.merge(rna[wanted], on="sample_sv_id", how="left")

    # Strongest evidence first. File order is genomic and carries no meaning to
    # a reader asking "what survived".
    rank = {"STRONG": 0, "SUGGESTIVE": 1, "WEAK": 2}
    tiers = frame["rna_tier"] if "rna_tier" in frame.columns \
        else pd.Series(index=frame.index, dtype=object)
    frame["_tier"] = tiers.map(lambda t: rank.get(str(t), 3))
    private = frame["is_private"] if "is_private" in frame.columns \
        else pd.Series(False, index=frame.index)
    frame["_notpriv"] = ~private.astype(str).str.lower().isin(("true", "1"))
    frame = frame.sort_values(["_tier", "_notpriv", "gene"], kind="stable")

    lines = ["| " + " | ".join(label for _, label in EVENT_COLUMNS) + " |",
             "|" + "---|" * len(EVENT_COLUMNS)]
    for _, row in frame.iterrows():
        population = row.get("gnomad_af_popmax_pop")
        if isinstance(population, float) and pd.isna(population):
            population = None
        lines.append("| " + " | ".join(
            cell(key, row.get(SOURCE.get(key, key)), population)
            for key, _ in EVENT_COLUMNS) + " |")
    return "\n".join(lines)


def analysis_block(name: str) -> str:
    """A generated analysis table, by name. Extend here, not in the document."""
    import analyse_locus_vs_peptide as analysis
    builders = {
        "locus_vs_peptide":
            lambda t: analysis.markdown(analysis.analyse(t)),
        "locus_vs_peptide_steps": analysis.steps_markdown,
        "locus_vs_peptide_per_line": analysis.per_line_markdown,
        "locus_vs_peptide_distance": analysis.distance_markdown,
    }
    if name not in builders:
        return f"*Unknown analysis block `{name}`.*"
    return builders[name](analysis.load("noPON"))


def sync(text: str) -> tuple[str, int]:
    """Replace every marked region with a freshly generated table."""
    # The body is matched as its own group, and may be empty: a freshly added
    # marker pair has nothing between it, and must still be filled on first run.
    pattern = re.compile(
        r"(<!--\s*events:(\S+?)\s*-->)(.*?)(<!--\s*/events:\2\s*-->)", re.S)
    filled = 0

    def replace(match: re.Match) -> str:
        nonlocal filled
        filled += 1
        return (f"{match.group(1)}\n{events_table(match.group(2))}\n"
                f"{match.group(4)}")

    text, filled = pattern.sub(replace, text), filled

    analysis = re.compile(
        r"(<!--\s*analysis:(\S+?)\s*-->)(.*?)(<!--\s*/analysis:\2\s*-->)", re.S)

    def replace_analysis(match: re.Match) -> str:
        nonlocal filled
        filled += 1
        return (f"{match.group(1)}\n{analysis_block(match.group(2))}\n"
                f"{match.group(4)}")

    return analysis.sub(replace_analysis, text), filled


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true",
                        help="exit 1 if a document is out of date; write nothing")
    parser.add_argument("--doc", action="append", default=[], metavar="PATH",
                        help="a document whose event tables to sync; repeatable")
    args = parser.parse_args()
    docs = args.doc or DOCS
    if not docs:
        raise SystemExit("no documents given — pass --doc PATH at least once")

    stale = []
    for name in docs:
        path = REPO / name
        if not path.exists():
            continue
        before = path.read_text()
        after, filled = sync(before)
        if not filled:
            continue
        if after == before:
            print(f"  {name}: {filled} table(s) already current")
            continue
        if args.check:
            stale.append(name)
            print(f"  {name}: {filled} table(s) OUT OF DATE")
        else:
            path.write_text(after)
            print(f"  {name}: {filled} table(s) regenerated")

    if stale:
        print("\nrun `python tools/sync_event_tables.py` to refresh")
        sys.exit(1)


if __name__ == "__main__":
    main()
