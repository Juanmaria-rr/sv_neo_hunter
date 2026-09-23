#!/usr/bin/env python
"""
build_universe_views.py — named subsets of the candidate universe, derived not made.

WHY NAMED VIEWS RATHER THAN HAND-CUT FILES
------------------------------------------
Three generations of master table already exist in this project, and the one
that looked most recent turned out to be stale, 255x smaller than the current
one, with nothing in the file saying so. Every extra hand-made subset multiplies
that surface.

So a subset here is never cut by hand. It is a named recipe applied to a frozen
universe file, and every output carries a `PROVENANCE.tsv` recording which
universe it came from, that file's checksum, the exact filter, and the resulting
counts. Re-running reproduces it; a stale copy is detectable by comparing one
checksum.

Add a view by adding one entry to `VIEWS`, not by writing another script.

Usage
-----
    python tools/build_universe_views.py --universe <dir>/candidate_universe.tsv \\
        --out-dir <dir>/views --view all
"""
from __future__ import annotations

import argparse
import hashlib
import pathlib
from datetime import date

import pandas as pd


def flag(frame: pd.DataFrame, column: str) -> pd.Series:
    """A TSV boolean, read as one. Missing column -> all False, never a crash."""
    if column not in frame.columns:
        return pd.Series(False, index=frame.index)
    return frame[column].astype(str).str.strip().str.lower().eq("true")


def carries(frame: pd.DataFrame, line: str) -> pd.Series:
    """Rows whose `present_in` includes this line — inherited peptides included."""
    return frame["present_in"].astype(str).str.split(";").apply(lambda v: line in v)


#: name -> (what it selects, one-line description)
#: Each predicate takes the universe and returns a boolean mask.
VIEWS: dict[str, tuple] = {
    "sequence_qc_passed": (
        lambda t: ~flag(t, "is_self") & ~flag(t, "low_complexity"),
        "Candidates that are neither proteome sequences nor compositionally "
        "trivial. The honest starting point for anything downstream."),
    "presentable": (
        lambda t: flag(t, "presentable"),
        "Binds at least one allele of the lineage's own class I genotype. "
        "Predicted binding, not observed presentation."),
    "presentable_robust": (
        lambda t: flag(t, "presentable") & t.get(
            "presentation_robustness", pd.Series(dtype=str)).eq("robust"),
        "Presentable with a comfortable margin on every cut — the subset that "
        "survives a re-run, a predictor version change or a threshold nudge."),
    "credible_and_presentable": (
        lambda t: (~flag(t, "is_self") & ~flag(t, "low_complexity")
                   & flag(t, "sv_hc") & flag(t, "expressed")
                   & t.get("rna_tier", pd.Series(dtype=str)).isin(
                       ["STRONG", "SUGGESTIVE"])
                   & flag(t, "presentable")),
        "Everything at once: QC-clean, on a high-confidence SV, transcribed, "
        "with junction-crossing RNA reads, and presentable by the line's own "
        "MHC. Check each survivor against `nearest_alt_sj_bp` before believing "
        "it — a breakpoint beside a splice site inherits that site's reads."),
    "in_reference_cohort": (
        lambda t: flag(t, "matched_reference"),
        "Candidates that ALSO appear in the reference cohort. The recurrence "
        "question, asked of the same rows."),
    "not_in_reference": (
        lambda t: ~flag(t, "matched_reference"),
        "Candidates absent from the reference cohort. NOT the default view: the "
        "universe is everything, and cohort membership is a column."),
}


def checksum(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()[:16]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--universe", type=pathlib.Path, required=True)
    parser.add_argument("--view", action="append", default=[],
                        help=f"one of: {', '.join(VIEWS)}, or `all`; repeatable")
    parser.add_argument("--lines", action="append", default=[], metavar="LINE",
                        help="restrict the universe to peptides ACQUIRED by these "
                             "lines before computing any view; repeatable. A "
                             "peptide acquired by a line outside the set is "
                             "dropped entirely — including one a selected line "
                             "would inherit, since it cannot be inherited from a "
                             "line that is not there. The exclusion is recorded "
                             "in every provenance file, because a subset that "
                             "does not say which lines it dropped is "
                             "indistinguishable from the full universe.")
    parser.add_argument("--subset-name", default=None,
                        help="directory name for a --lines restriction; defaults "
                             "to the lines joined by '+'")
    parser.add_argument("--per-line", action="store_true",
                        help="also write one file per line in `present_in`, "
                             "carrying everything that line has including what "
                             "it inherited")
    parser.add_argument("--out-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()

    names = list(VIEWS) if "all" in args.view or not args.view else args.view
    unknown = [n for n in names if n not in VIEWS]
    if unknown:
        raise SystemExit(f"unknown view(s): {', '.join(unknown)}")

    universe = pd.read_csv(args.universe, sep="\t", dtype=str, low_memory=False)
    source_sum = checksum(args.universe)
    full_rows = len(universe)
    print(f"  universe: {full_rows:,} peptides  sha256:{source_sum}")

    restriction = ""
    if args.lines:
        known = {line for value in universe["acquired_in"].astype(str)
                 for line in value.split(";") if line}
        unknown = [l for l in args.lines if l not in known]
        if unknown:
            raise SystemExit(f"unknown line(s): {', '.join(unknown)}. "
                             f"Present: {', '.join(sorted(known))}")
        selected = set(args.lines)
        keep = universe["acquired_in"].astype(str).apply(
            lambda v: set(v.split(";")) <= selected)
        dropped = sorted(known - selected)
        universe = universe[keep]
        restriction = (f"acquired_in within {{{', '.join(sorted(selected))}}}; "
                       f"excluded {{{', '.join(dropped) or 'none'}}}. "
                       "NOTE: `present_in` still names the excluded lines where "
                       "they inherit a kept peptide — that is a true fact about "
                       "the peptide and is not rewritten. Filter on "
                       "`acquired_in` to stay inside this subset.")
        args.out_dir = args.out_dir / (args.subset_name
                                       or "+".join(sorted(selected)))
        print(f"  restricted to {', '.join(sorted(selected))}: "
              f"{len(universe):,} peptides "
              f"(dropped {full_rows - len(universe):,} acquired by "
              f"{', '.join(dropped)})")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest = []

    def emit(subdir: pathlib.Path, name: str, subset: pd.DataFrame,
             description: str, filter_text: str) -> None:
        subdir.mkdir(parents=True, exist_ok=True)
        path = subdir / f"{name}.tsv"
        subset.to_csv(path, sep="\t", index=False)
        pd.DataFrame([
            {"key": "view", "value": name},
            {"key": "description", "value": description},
            {"key": "filter", "value": filter_text},
            {"key": "generated", "value": date.today().isoformat()},
            {"key": "tool", "value": pathlib.Path(__file__).name},
            {"key": "source", "value": str(args.universe)},
            {"key": "source_sha256_16", "value": source_sum},
            {"key": "line_restriction", "value": restriction or "none — all lines"},
            {"key": "rows", "value": len(subset)},
            {"key": "rows_after_restriction", "value": len(universe)},
            {"key": "rows_in_full_universe", "value": full_rows},
        ]).to_csv(subdir / f"{name}.PROVENANCE.tsv", sep="\t", index=False)
        manifest.append({"file": str(path.relative_to(args.out_dir)),
                         "view": name, "rows": len(subset),
                         "pct_of_restricted": round(100.0 * len(subset) / len(universe), 2),
                         "line_restriction": restriction or "none",
                         "source_sha256_16": source_sum,
                         "generated": date.today().isoformat()})
        print(f"    {name:<26} {len(subset):>7,}")

    # The restricted universe itself, not only its views: a subset that ships
    # only derived cuts forces anyone asking a new question back to the full
    # table, where the excluded lines are present again.
    if args.lines:
        emit(args.out_dir, "candidate_universe", universe,
             "The candidate universe restricted to the selected lines. This is "
             "the subset's master table; every view beside it derives from it.",
             restriction)

    print("\n  views over the restricted universe:"
          if args.lines else "\n  views over the whole universe:")
    for name in names:
        predicate, description = VIEWS[name]
        emit(args.out_dir, name, universe[predicate(universe)], description,
             f"VIEWS['{name}'] in {pathlib.Path(__file__).name}")

    if args.per_line:
        lines = sorted({line for value in universe["present_in"].astype(str)
                        for line in value.split(";") if line})
        if args.lines:
            # A line outside the restriction still appears in `present_in`,
            # because it inherits what the selected lines acquired — and after
            # the restriction it carries exactly what its parent does, so its
            # directory would be a duplicate under a name the subset excludes.
            lines = [line for line in lines if line in set(args.lines)]
        for line in lines:
            print(f"\n  {line} (everything it carries, inherited included):")
            held = universe[carries(universe, line)]
            subdir = args.out_dir / line
            emit(subdir, "all_candidates", held,
                 f"Every candidate {line} carries, acquired or inherited.",
                 f"present_in contains '{line}'")
            for name in names:
                predicate, description = VIEWS[name]
                emit(subdir, name, held[predicate(held)],
                     f"{description} Restricted to what {line} carries.",
                     f"present_in contains '{line}' AND VIEWS['{name}']")

    pd.DataFrame(manifest).to_csv(args.out_dir / "MANIFEST.tsv",
                                  sep="\t", index=False)
    print(f"\n  wrote {len(manifest)} files + MANIFEST.tsv under {args.out_dir}")
    print("  every file carries a .PROVENANCE.tsv naming its source checksum — "
          "compare it against the universe to tell current from stale.")


if __name__ == "__main__":
    main()
