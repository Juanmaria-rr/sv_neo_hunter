#!/usr/bin/env python
"""
build_combined_master.py — one table per candidate peptide, across cell lines.

WHY A THIRD TABLE
-----------------
Two analyses have been run over the same samples and they live in different
files at different grains:

  master_peptides.tsv         one row per peptide that MATCHED a reference
                              catalogue, carrying the cross and the MHC layer
  candidate_neoantigens.tsv   one row per CANDIDATE peptide, carrying sequence
                              QC and presentability by the line's own genotype

The candidate universe is the superset: every matched peptide is a candidate,
but most candidates never matched. Reproducing either analysis from a single
file therefore means starting from the candidate grain and joining the cross
onto it, which is what this script does.

A peptide that never matched gets EMPTY cross columns, not False. Those two are
different facts and the distinction is the whole point of keeping the columns
rather than a verdict: `matched_reference` says whether the cross had anything
to say, and only then do its columns mean anything.

WHICH CELL LINE A PEPTIDE COMES FROM
------------------------------------
A derived clone is analysed from its somatic VCF alone, because its germline
genome IS the parent's. So a peptide seen in the clone's catalogue was either
acquired by the clone or inherited from the parent, and `cell_line` records
which run produced it. `peptide_origin` adds a third value, `both`, for a
peptide that arises on both sides — inherited AND regenerated, which is a fact
about the peptide rather than a duplicate row to drop.

Usage
-----
    python tools/build_combined_master.py \\
        --catalogue <results>/<clone>/neoantigens/candidate_neoantigens.tsv \\
        --master    <results>/master_peptides.tsv \\
        --sample parental_noPON --sample clone_A_noPON \\
        --out-dir   <results>/combined_master
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import column_meanings                                   # noqa: E402

#: Cross and MHC columns worth carrying onto a candidate. Everything else in the
#: matched table is either a duplicate of a candidate column or an artefact of
#: the merge that produced it.
CROSS_COLUMNS = [
    "ref_peptide", "ref_gene", "ref_svtype", "n_patients", "n_rows",
    "gene_concordant", "svtype_concordant",
    "n_alleles_panel", "binds_panel", "panel_margin", "panel_limiting_cut",
    "panel_best_allele", "n_patients_typed", "n_patients_presenting",
    "n_patients_unevaluable", "binds_autologous", "autologous_verdict",
    "presenting_alleles", "autologous_margin", "autologous_limiting_cut",
    "autologous_best_allele", "autologous_robustness",
    "binds_autologous_neo_strong", "autologous_verdict_neo_strong",
    "autologous_margin_neo_strong", "autologous_robustness_neo_strong",
]

#: Renamed on the way in, so a cross column cannot be mistaken for the
#: cell-line-autologous one that shares its concept.
RENAME = {
    "binds_autologous": "binds_autologous_patient",
    "autologous_verdict": "autologous_verdict_patient",
    "autologous_margin": "autologous_margin_patient",
    "autologous_limiting_cut": "autologous_limiting_cut_patient",
    "autologous_best_allele": "autologous_best_allele_patient",
    "autologous_robustness": "autologous_robustness_patient",
}

#: Column order in the output: identity first, then what the peptide is, then
#: how good the call is, then what it means immunologically.
ORDER = [
    "cell_line", "branch", "peptide_origin", "neopeptide", "pep_length",
    "sv_id", "chrom1", "pos1", "gene1", "strand1", "chrom2", "pos2", "gene2",
    "strand2", "svtype", "frame_effect", "spans_junction", "junction_aa",
    "low_complexity", "is_self",
    "sv_hc", "vf_bp1", "vf_bp2", "qual_bp1", "qual_bp2", "segmapq_bp1",
    "segmapq_bp2", "pon_count", "pon_fraction", "pass_pon", "gnomad_af_used",
    "gnomad_af_popmax", "pass_gnomad", "is_private", "privacy_note",
    "gene1_TPM", "gene2_TPM", "min_side_TPM", "expressed", "rna_tier",
    "junction_reads", "test", "test_reason",
    "nearest_alt_sj_bp", "alt_sj_frags", "alt_sj_type", "alt_sj_within_window",
    "retained_intron", "isofox_fusion", "isofox_fusion_support",
    "presentable", "n_alleles_binding", "binding_alleles",
    "presentation_best_allele", "presentation_margin",
    "presentation_limiting_cut", "presentation_robustness",
    "matched_reference",
]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--catalogue", type=pathlib.Path, required=True,
                        help="candidate_neoantigens.tsv covering the samples")
    parser.add_argument("--master", type=pathlib.Path, required=True,
                        help="master_peptides.tsv — the cross and MHC layer")
    parser.add_argument("--sample", action="append", required=True,
                        metavar="NAME",
                        help="a run to include, e.g. parental_noPON; repeatable")
    parser.add_argument("--origin-map", action="append", default=[],
                        metavar="ORIGIN=CELL_LINE",
                        help="name the cell line behind each origin, e.g. "
                             "parental=parental somatic=clone_A. Without it "
                             "`cell_line` is the origin label itself.")
    parser.add_argument("--out-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()

    mapping = dict(pair.split("=", 1) for pair in args.origin_map)

    table = pd.read_csv(args.catalogue, sep="\t", dtype=str, low_memory=False)
    table = table.drop_duplicates("neopeptide")
    print(f"  candidates: {len(table):,}")

    table["cell_line"] = table["peptide_origin"].map(
        lambda o: mapping.get(o, o))
    # `both` has no single cell line; say so rather than picking one.
    table.loc[table["peptide_origin"] == "both", "cell_line"] = "both"

    master = pd.read_csv(args.master, sep="\t", dtype=str, low_memory=False)
    master = master[master["sample"].isin(args.sample)]
    if master.empty:
        raise SystemExit(f"no rows in {args.master.name} for "
                         f"{', '.join(args.sample)}")
    print(f"  matched peptides in {', '.join(args.sample)}: "
          f"{master['neopeptide'].nunique():,}")

    branch = master["branch"].dropna().unique()
    table["branch"] = branch[0] if len(branch) == 1 else ";".join(sorted(branch))

    keep = ["neopeptide"] + [c for c in CROSS_COLUMNS if c in master.columns]
    cross = master[keep].drop_duplicates("neopeptide").rename(columns=RENAME)
    cross["matched_reference"] = True
    table = table.merge(cross, on="neopeptide", how="left",
                        suffixes=("", "_cross"))
    table["matched_reference"] = table["matched_reference"].notna()
    matched = int(table["matched_reference"].astype(bool).sum())
    print(f"  of the candidates, {matched:,} matched the reference catalogue")

    ordered = [c for c in ORDER if c in table.columns]
    rest = [c for c in table.columns
            if c not in ordered and not c.endswith("_cross")
            and c not in ("origin", "origin_run", "sample_sv_id", "prefix",
                          "transcript_id1", "transcript_id2", "junction_nt")]
    table = table[ordered + rest]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out_dir / "master_candidates.tsv"
    table.to_csv(out, sep="\t", index=False)

    # -- the dictionary -------------------------------------------------------
    # Generated from the same source the other master tables use, so a column
    # cannot mean one thing here and another there. Columns this script creates
    # are described here because they exist nowhere else.
    local = {
        "cell_line": "Which line the peptide's junction comes from. A "
                     "derived clone is analysed from its somatic VCF alone "
                     "(its germline genome IS the parent's), so a peptide is "
                     "either acquired by the clone or inherited from the "
                     "parent. `both` means it arises on both sides — inherited "
                     "AND regenerated, which is a fact about the peptide, not a "
                     "duplicate row.",
        "branch": "The sensitivity branch these rows come from. `noPON` means "
                  "the panel-of-normals filter was NOT applied: `pon_count` is "
                  "reported as a column so it can be filtered on deliberately, "
                  "and `is_private` is therefore a gnomAD-only judgement.",
        "peptide_origin": "`somatic`, `parental` or `both` — the raw origin "
                          "label that `cell_line` names.",
        "matched_reference": "Whether this peptide matched the reference "
                             "catalogue. **When False, every cross and "
                             "patient-side MHC column is empty because the "
                             "question was never asked — not because the answer "
                             "was no.**",
        "presentable": "Binds at least one allele of the cell line's OWN class I "
                       "genotype. For a cell line this is the autologous "
                       "question and there is no panel branch to confuse it "
                       "with: the alleles ARE the line's, so this means the cell "
                       "could present the peptide. Predicted binding, never "
                       "observed presentation.",
        "n_alleles_binding": "How many of the line's own alleles it binds.",
        "binding_alleles": "Which ones, `;`-separated.",
        "presentation_best_allele": "The line's allele giving the largest margin.",
        "presentation_margin": "How far inside the binder definition the call "
                               "sits: the distance to each of the three cuts as "
                               "a fraction of that cut, the smallest of the "
                               "three, then the largest across binding alleles. "
                               "0 is exactly on a boundary. A heuristic for "
                               "fragility, not a calibrated probability.",
        "presentation_limiting_cut": "Which cut is closest to failing.",
        "presentation_robustness": "`presentation_margin` binned: `flippable` "
                                   "(<=10% of a cut), `marginal` (<=25%), "
                                   "`solid` (<=50%), `robust` (>50%). The robust "
                                   "subset is the same criterion minus the calls "
                                   "a re-run or a version change would flip.",
        "binds_autologous_patient": "The PATIENT-side autologous verdict from "
                                    "the cross — binds an allele of a catalogue "
                                    "patient who carries this peptide. Distinct "
                                    "from `presentable`, which is about this "
                                    "cell line. Empty unless "
                                    "`matched_reference`.",
        "autologous_verdict_patient": "Patient-side verdict: `binder`, "
                                      "`non_binder`, or `unevaluable_*`. "
                                      "Unevaluable is NOT non-binder.",
        "autologous_margin_patient": "Patient-side margin.",
        "autologous_limiting_cut_patient": "Patient-side limiting cut.",
        "autologous_best_allele_patient": "Patient-side presenting allele.",
        "autologous_robustness_patient": "Patient-side robustness bin.",
    }

    rows = []
    for column in table.columns:
        meaning = local.get(column) or column_meanings.resolve(column)
        filled = table[column].notna().sum()
        rows.append({
            "column": column,
            "meaning": meaning or "UNDOCUMENTED — add an entry to "
                                  "tools/column_meanings.py",
            "n_populated": filled,
            "pct_populated": round(100.0 * filled / len(table), 1),
            "example": next((str(v) for v in table[column].dropna().head(3)), ""),
        })
    dictionary = pd.DataFrame(rows)
    dictionary.to_csv(args.out_dir / "master_candidates_column_dictionary.tsv",
                      sep="\t", index=False)

    missing = int((dictionary["meaning"].str.startswith("UNDOCUMENTED")).sum())
    print(f"\n  wrote {out}")
    print(f"        {len(table):,} rows x {len(table.columns)} columns")
    print(f"        {args.out_dir / 'master_candidates_column_dictionary.tsv'}"
          f" ({missing} undocumented)")
    counts = table["cell_line"].value_counts()
    for line, n in counts.items():
        m = int(table.loc[table["cell_line"] == line, "matched_reference"]
                .astype(bool).sum())
        print(f"        {line:<14} {n:>7,} candidates, {m:>6,} matched")


if __name__ == "__main__":
    main()
