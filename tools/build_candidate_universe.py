#!/usr/bin/env python
"""
build_candidate_universe.py — every candidate neopeptide a lineage can produce.

WHAT THIS IS, AND WHAT IT DELIBERATELY IS NOT
---------------------------------------------
The candidate universe of a whole lineage: one row per distinct peptide, across
every line, **regardless of whether it also appears in a reference cohort**.
Whether it does is carried as a column (`matched_reference`) and never as a
filter, because the question here is what these cells can make, not what they
share with anyone.

That is the opposite selection from the recurrence test, which keeps only the
peptides that DO match. Both views come from the same rows; only the column you
read changes. Nothing in this file is filtered on cohort membership.

ONE ROW PER PEPTIDE, NOT ONE PER (PEPTIDE, LINE)
------------------------------------------------
A derived clone inherits its parent's genome, so a peptide acquired by the
parental line is present in every descendant. Repeating it per line would
multiply the parental set by the number of clones and invite double-counting in
any downstream statistic. Instead:

  `acquired_in`  the line whose OWN variant calls produce it — where it entered
                 the lineage. Semicolon-separated if more than one line
                 generates the same sequence independently.
  `present_in`   every line that carries it: the acquiring line plus all its
                 descendants. This is the column to filter on for "what can line
                 X present".

So "candidates in clone_A" is `present_in` containing clone_A, which
includes what it inherited; "what a given knockout added" is `acquired_in`.

Usage
-----
    python tools/build_candidate_universe.py \\
        --line parental:<dir>:                       \\
        --line clone_A:<dir>:parental              \\
        --line clone_B:<dir>:clone_A      \\
        --master <results>/master_peptides.tsv      \\
        --sample parental_noPON --sample clone_A_noPON \\
        --out-dir <results>/universe
"""
from __future__ import annotations

import argparse
import hashlib
import pathlib
import sys
from datetime import date

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import column_meanings                                   # noqa: E402
from build_combined_master import CROSS_COLUMNS, RENAME  # noqa: E402

#: Columns that describe the peptide and its junction. Taken from the line that
#: acquired the peptide, since that is where its junction actually is.
PEPTIDE_COLUMNS = [
    "neopeptide", "pep_length", "sv_id", "chrom1", "pos1", "gene1", "strand1",
    "chrom2", "pos2", "gene2", "strand2", "svtype", "frame_effect",
    "spans_junction", "junction_aa", "low_complexity", "is_self",
    "sv_hc", "vf_bp1", "vf_bp2", "qual_bp1", "qual_bp2", "segmapq_bp1",
    "segmapq_bp2", "pon_count", "pon_fraction", "pass_pon", "gnomad_af_used",
    "gnomad_af_popmax", "pass_gnomad", "is_private", "privacy_note",
    "gene1_TPM", "gene2_TPM", "min_side_TPM", "expressed", "rna_tier",
    "junction_reads", "test", "test_reason",
    # Junction evidence broken down by mechanism, and the lesion size. Which
    # mechanism produced a count says which geometry the junction has, not
    # whether it is real: an N-gap population is what BOTH a deleted allele and
    # an alternative splice junction produce. `event_size` is type-aware where
    # `span` is the raw coordinate difference, so the two disagree whenever an
    # insertion is involved.
    "junction_by_ngap", "junction_by_sa", "junction_by_insert",
    "coverage_bp1", "coverage_bp2", "min_coverage",
    "softclip_bp1", "softclip_bp2",
    "alignments_bp1", "alignments_bp2", "low_mapq_bp1", "low_mapq_bp2",
    "span", "event_size", "insert_len",
    "nearest_alt_sj_bp", "alt_sj_frags", "alt_sj_type", "alt_sj_within_window",
    "retained_intron", "isofox_fusion", "isofox_fusion_support",
    "presentable", "n_alleles_binding", "binding_alleles",
    "presentation_best_allele", "presentation_margin",
    "presentation_limiting_cut", "presentation_robustness",
]


def descendants(line: str, parents: dict[str, str]) -> list[str]:
    """Every line downstream of `line`, itself included."""
    out = [line]
    changed = True
    while changed:
        changed = False
        for child, parent in parents.items():
            if parent in out and child not in out:
                out.append(child)
                changed = True
    return out


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
    parser.add_argument("--line", action="append", required=True,
                        metavar="NAME:DIR:PARENT",
                        help="a line, its acquired-catalogue directory, and its "
                             "parent line (empty for the root); repeatable")
    parser.add_argument("--master", type=pathlib.Path,
                        help="master_peptides.tsv, to record which candidates "
                             "also matched a reference cohort. Optional: "
                             "without it `matched_reference` is absent rather "
                             "than False, because unasked is not no.")
    parser.add_argument("--sample", action="append", default=[],
                        help="runs of --master to consider; repeatable")
    parser.add_argument("--out-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()

    lines, parents = {}, {}
    for spec in args.line:
        name, directory, parent = (spec.split(":") + ["", ""])[:3]
        lines[name] = pathlib.Path(directory)
        if parent:
            parents[name] = parent

    frames = []
    for name, directory in lines.items():
        path = directory / "candidate_neoantigens.tsv"
        if not path.exists():
            raise SystemExit(f"missing {path}")
        frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        frame = frame.drop_duplicates("neopeptide")
        frame["acquired_in"] = name
        frames.append(frame)
        print(f"  {name:<22} {len(frame):>7,} acquired")

    stacked = pd.concat(frames, ignore_index=True)

    # A sequence can be generated independently by two lines. That is one
    # peptide with two origins, not two peptides: collapse, keeping every label.
    acquired = (stacked.groupby("neopeptide")["acquired_in"]
                       .agg(lambda s: ";".join(sorted(set(s)))))
    table = stacked.drop_duplicates("neopeptide").drop(columns=["acquired_in"])
    table = table.merge(acquired.rename("acquired_in"), on="neopeptide")

    reach = {name: descendants(name, parents) for name in lines}
    table["present_in"] = table["acquired_in"].map(
        lambda labels: ";".join(sorted(
            {line for label in labels.split(";") for line in reach[label]})))
    table["n_lines_present"] = table["present_in"].str.count(";") + 1
    print(f"\n  distinct peptides in the universe: {len(table):,}")

    # -- did it also match a cohort? a column, never a filter ------------------
    if args.master:
        master = pd.read_csv(args.master, sep="\t", dtype=str, low_memory=False)
        if args.sample:
            master = master[master["sample"].isin(args.sample)]
        keep = ["neopeptide"] + [c for c in CROSS_COLUMNS if c in master.columns]
        cross = master[keep].drop_duplicates("neopeptide").rename(columns=RENAME)
        cross["matched_reference"] = True
        table = table.merge(cross, on="neopeptide", how="left")
        table["matched_reference"] = table["matched_reference"].notna()
        print(f"  of those, also in the reference cohort: "
              f"{int(table['matched_reference'].sum()):,} "
              f"({100 * table['matched_reference'].mean():.1f}%) "
              "— reported, NOT filtered on")

    ordered = (["acquired_in", "present_in", "n_lines_present"]
               + [c for c in PEPTIDE_COLUMNS if c in table.columns]
               + [c for c in table.columns
                  if c.startswith(("binds_", "autologous_", "panel_", "n_patients",
                                   "n_alleles_panel", "presenting_alleles",
                                   "ref_", "matched_reference", "gene_concordant",
                                   "svtype_concordant", "n_rows"))])
    seen, columns = set(), []
    for c in ordered:
        if c in table.columns and c not in seen:
            seen.add(c)
            columns.append(c)
    table = table[columns]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out_dir / "candidate_universe.tsv"
    table.to_csv(out, sep="\t", index=False)

    # -- dictionary -----------------------------------------------------------
    local = {
        "acquired_in": "The line whose OWN variant calls produce this peptide — "
                       "where it entered the lineage. Semicolon-separated if two "
                       "lines generate the same sequence independently. This is "
                       "the column for 'what did this knockout add'.",
        "present_in": "Every line carrying the peptide: the acquiring line plus "
                      "all its descendants, since a clone inherits its parent's "
                      "genome. **This is the column for 'what can line X "
                      "present'** — filtering on `acquired_in` would silently "
                      "drop everything the line inherited.",
        "n_lines_present": "How many lines carry it; 1 means private to one "
                           "branch of the lineage.",
        "matched_reference": "Whether this peptide also appears in the reference "
                             "cohort. **Reported, never filtered on: this table "
                             "is the candidate universe, so cohort membership is "
                             "a property of a peptide here, not an entry "
                             "condition.** When False, every cross and "
                             "patient-side column is empty because the question "
                             "was not asked, not because the answer was no.",
        "presentable": "Binds at least one allele of the lineage's own class I "
                       "genotype — the cells could present it. Predicted "
                       "binding, never observed presentation.",
        "n_alleles_binding": "How many of the line's own alleles it binds.",
        "binding_alleles": "Which ones, `;`-separated.",
        "presentation_best_allele": "The allele giving the largest margin.",
        "presentation_margin": "Distance to each of the three binder cuts as a "
                               "fraction of that cut, the smallest of the three, "
                               "then the largest across binding alleles. 0 is "
                               "exactly on a boundary. A fragility heuristic, "
                               "not a calibrated probability.",
        "presentation_limiting_cut": "Which cut is closest to failing.",
        "presentation_robustness": "`flippable` (<=10% of a cut), `marginal` "
                                   "(<=25%), `solid` (<=50%), `robust` (>50%).",
        "binds_autologous_patient": "The cohort-side verdict: binds an allele of "
                                    "a catalogue patient carrying this peptide. "
                                    "Distinct from `presentable`, which is this "
                                    "lineage's own genotype. Empty unless "
                                    "`matched_reference`.",
        "autologous_verdict_patient": "Cohort-side verdict: `binder`, "
                                      "`non_binder` or `unevaluable_*`. "
                                      "Unevaluable is NOT non-binder.",
        "autologous_margin_patient": "Cohort-side margin.",
        "autologous_limiting_cut_patient": "Cohort-side limiting cut.",
        "autologous_best_allele_patient": "Cohort-side presenting allele.",
        "autologous_robustness_patient": "Cohort-side robustness bin.",
    }
    rows = []
    for column in table.columns:
        meaning = local.get(column) or column_meanings.resolve(column)
        filled = int(table[column].notna().sum())
        rows.append({"column": column,
                     "meaning": meaning or "UNDOCUMENTED — add an entry to "
                                           "tools/column_meanings.py",
                     "n_populated": filled,
                     "pct_populated": round(100.0 * filled / len(table), 1),
                     "example": next((str(v) for v in table[column].dropna().head(1)), "")})
    dictionary = pd.DataFrame(rows)
    dictionary.to_csv(args.out_dir / "candidate_universe_column_dictionary.tsv",
                      sep="\t", index=False)

    # -- provenance -----------------------------------------------------------
    # Stamped so a derived subset can always be told apart from a stale one:
    # three generations of master table already exist on disk and only the
    # checksum distinguishes them.
    provenance = [
        {"key": "generated", "value": date.today().isoformat()},
        {"key": "tool", "value": pathlib.Path(__file__).name},
        {"key": "rows", "value": len(table)},
        {"key": "columns", "value": len(table.columns)},
        {"key": "cohort_filter_applied", "value": "NO — matched_reference is a column"},
        {"key": "output_sha256_16", "value": checksum(out)},
    ]
    for name, directory in lines.items():
        source = directory / "candidate_neoantigens.tsv"
        provenance.append({"key": f"source[{name}]",
                           "value": f"{source} sha256:{checksum(source)}"})
    if args.master:
        provenance.append({"key": "master",
                           "value": f"{args.master} sha256:{checksum(args.master)}"})
    pd.DataFrame(provenance).to_csv(args.out_dir / "PROVENANCE.tsv",
                                    sep="\t", index=False)

    missing = int(dictionary["meaning"].str.startswith("UNDOCUMENTED").sum())
    print(f"\n  wrote {out}")
    print(f"        {len(table):,} rows x {len(table.columns)} columns, "
          f"{missing} undocumented")
    print(f"        PROVENANCE.tsv, column dictionary")
    print("\n  acquired_in:")
    for label, n in table["acquired_in"].value_counts().items():
        print(f"    {label:<40} {n:>7,}")
    print("  present_in (what each line carries):")
    for name in lines:
        n = int(table["present_in"].str.split(";").apply(lambda v: name in v).sum())
        print(f"    {name:<40} {n:>7,}")


if __name__ == "__main__":
    main()
