#!/usr/bin/env python3
"""
build_master_table.py — the two auditable tables, per sample and combined.

DATA, NOT VERDICTS
------------------
These tables carry the VALUES each criterion was judged on, not the judgement.
`pass_sv_confidence = False` hides whether the mapping quality, the fragment
support, the caller quality or the size floor was the reason; `segmapq_bp1 = 12,
vf_bp1 = 3, qual_bp1 = 51, event_size = 34` says exactly which. The verdicts live
in `credible_events.tsv` and `funnel.tsv`, computed from these same values, so
nothing is lost by keeping them apart — and anyone can re-judge a call, or apply
a different threshold, without re-running anything.

TWO GRAINS, BECAUSE THEY ANSWER DIFFERENT QUESTIONS
---------------------------------------------------
`master_peptides.tsv`  one row per candidate peptide that matched the reference.
                       Use it to ask why a *candidate* did or did not survive.
`master_sv.tsv`        one row per ADMITTED junction, whether or not it produced
                       a peptide or a match. Use it to ask what the sample's SV
                       call set looks like and which parts of it went anywhere.

The SV table deliberately includes junctions that produced nothing. A table of
only the successful ones cannot answer "how much of the call set is inert", and
that denominator is what makes a match count interpretable.

WHY BOTH, RATHER THAN ONE COLLAPSED TABLE
-----------------------------------------
A peptide-grain table cannot be collapsed to SVs without losing the per-peptide
criteria (complexity, self-proteome, catalogue annotation); an SV-grain table
cannot be expanded without inventing rows. Each is emitted at the grain its
criteria are evaluated at, and `sample_sv_id` joins them.

`NA` means not measured, never "failed": a junction that never reached the RNA
stage has no `rna_tier`, which is not the same as having no reads.

USAGE
-----
    python tools/build_master_table.py --results-dir results
"""
from __future__ import annotations

import argparse
import os

import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import column_meanings  # noqa: E402

#: Peptide-grain gates, in the order the pipeline applies them.
PEPTIDE_GATES = [
    ("pass_identical", "peptide string identical to a catalogue entry"),
    ("pass_gene_concordant", "the gene broken here is the catalogue's gene"),
    ("pass_complexity", "not a low-complexity sequence"),
    ("pass_not_self", "not present verbatim in the normal proteome"),
    ("pass_private_panel", "breakpoint below the panel-of-normals threshold"),
    ("pass_private_population", "below the population allele-frequency threshold"),
    ("pass_sv_confidence", "both breakends meet mapping quality, support and size"),
    ("pass_rna_junction", "reads cross the junction in RNA"),
]

#: SV-grain gates. An admitted junction can fail earlier than any peptide gate —
#: by producing no peptide at all, which no peptide-grain table can show.
SV_GATES = [
    ("pass_admitted", "survived FILTER, pairing and the panel filter"),
    ("pass_generated_peptides", "the annotator produced at least one peptide"),
    ("pass_matched", "at least one peptide matched the reference catalogue"),
    ("pass_credible", "at least one match survived sequence QC"),
    ("pass_private_panel", "below the panel-of-normals threshold"),
    ("pass_private_population", "below the population allele-frequency threshold"),
    ("pass_sv_confidence", "both breakends meet mapping quality, support and size"),
    ("pass_rna_junction", "reads cross the junction in RNA"),
]

RNA_SUPPORTED = ("STRONG", "SUGGESTIVE", "WEAK")


def _read(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False) \
        if os.path.exists(path) else pd.DataFrame()


def _flag(frame: pd.DataFrame, column: str, default=False) -> pd.Series:
    """Read a stringly-typed boolean column; absent means `default`."""
    if column not in frame.columns:
        return pd.Series([default] * len(frame), index=frame.index)
    return frame[column].astype(str).str.lower().isin(("true", "1", "yes"))


def _numeric(frame: pd.DataFrame, columns) -> pd.DataFrame:
    for column in columns:
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


NUMERIC = ["or_original", "or_clean", "hla_pres_cov_mean4", "pon_count", "gnomad_af",
           "event_size", "insert_len", "span", "vf_bp1", "vf_bp2", "qual_bp1",
           "qual_bp2", "segmapq_bp1", "segmapq_bp2", "junction_reads", "coverage_bp1",
           "coverage_bp2", "min_coverage", "min_side_TPM", "pep_length",
           "junction_aa", "n_peptides", "pos1", "pos2", "pon_fraction"]


def _first_failure(row, gates) -> str:
    for column, _ in gates:
        value = row.get(column)
        if pd.isna(value):
            return f"{column} (not evaluated)"
        if not bool(value):
            return column
    return "-"


def build_peptide_table(run_dir: str, label: str) -> pd.DataFrame:
    """One row per matched candidate peptide, with every criterion."""
    matches = _read(os.path.join(run_dir, "stage3_matches.tsv"))
    if matches.empty:
        return pd.DataFrame()

    table = matches.copy()
    identify(table, label)
    table["sample_sv_id"] = table["sv_id"].astype(str)

    for extra in (_read(os.path.join(run_dir, "credible_events.tsv")),
                  _read(os.path.join(run_dir, "stage7_expression.tsv")),
                  _read(os.path.join(run_dir, "stage8_rna_evidence.tsv"))):
        if extra.empty:
            continue
        extra = extra.drop(columns=[c for c in ("gene",) if c in extra.columns])
        extra["sample_sv_id"] = extra["sample_sv_id"].astype(str)
        table = table.merge(extra, on="sample_sv_id", how="left", suffixes=("", "_dup"))
    table = table.drop(columns=[c for c in table.columns if c.endswith("_dup")])

    # MHC binding, when it was run. It joins on the PEPTIDE, not on
    # `sample_sv_id` like everything else above: binding is a property of the
    # sequence and one peptide can arise from several junctions. Optional —
    # binding is not part of the recurrence test, and most samples have no
    # prediction at all, so its absence leaves the columns out rather than
    # filling them with a value that would read as measured.
    # Alternative binder criteria land in `peptide_binding_<criterion>.tsv` and
    # are merged with a suffixed name, so a candidate carries every verdict it
    # has rather than only the default one. `three_way` keeps the bare names
    # because it is the comparable rule and everything downstream expects them.
    for suffix, filename in [("", "peptide_binding.tsv")] + [
            (f"_{c}", f"peptide_binding_{c}.tsv")
            for c in ("neo_strong", "neo_weak")]:
        binding = _read(os.path.join(run_dir, filename))
        if binding.empty or "peptide" not in binding.columns:
            continue
        binding = binding.drop_duplicates("peptide")
        if suffix:
            binding = binding.rename(columns={c: f"{c}{suffix}"
                                              for c in binding.columns
                                              if c != "peptide"})
        overlap = [c for c in binding.columns
                   if c != "peptide" and c in table.columns]
        table = table.merge(binding.drop(columns=overlap), on="peptide", how="left")

    # How many peptides the same junction contributed — context for reading one
    # row, since a sliding window over one locus produces many.
    table["n_peptides_in_event"] = table.groupby("sample_sv_id")["peptide"] \
                                        .transform("nunique")
    # Derived booleans that ARE data rather than verdicts (they record what the
    # sequence is, not whether it passed): low_complexity, is_self,
    # gene_concordant, spans_junction all stay as the generator/QC emitted them.
    drop = [c for c in table.columns if c.startswith("pass_") or c == "credible"]
    return _numeric(table.drop(columns=drop), NUMERIC)


def build_sv_table(run_dir: str, label: str) -> pd.DataFrame:
    """One row per admitted junction, including those that produced nothing."""
    junctions = _read(os.path.join(run_dir, "stage1_junctions.tsv"))
    if junctions.empty:
        return pd.DataFrame()

    table = junctions.copy()
    identify(table, label)
    table["sample_sv_id"] = table["sample_sv_id"].astype(str)

    # Gene and transcript annotation, from the generator's own output.
    prefix = label.split("_")[0]
    anno = _read(os.path.join(run_dir, f"{prefix}.anno.txt"))
    if not anno.empty and "sv_id" in anno.columns:
        keep = [c for c in ("sv_id", "gene1", "transcript_id1", "strand1", "gene2",
                            "transcript_id2", "strand2") if c in anno.columns]
        anno = anno[keep].drop_duplicates("sv_id").rename(columns={"sv_id": "sample_sv_id"})
        table = table.merge(anno, on="sample_sv_id", how="left", suffixes=("", "_anno"))

    # How far each junction got: peptides generated, matched, credible.
    peptides = _read(os.path.join(run_dir, f"{prefix}.all_neopeptides.txt"))
    if not peptides.empty and "sv_id" in peptides.columns:
        generated = (peptides.groupby("sv_id")["neopeptide"].nunique()
                     .rename("n_peptides_generated"))
        table = table.merge(generated, left_on="sample_sv_id", right_index=True, how="left")
    table["n_peptides_generated"] = table.get(
        "n_peptides_generated", pd.Series(dtype=float)).fillna(0).astype(int)

    matches = _read(os.path.join(run_dir, "stage3_matches.tsv"))
    if not matches.empty:
        matches["sv_id"] = matches["sv_id"].astype(str)
        matched = matches.groupby("sv_id")["peptide"].nunique().rename("n_peptides_matched")
        credible = (matches[_flag(matches, "credible")].groupby("sv_id")["peptide"]
                    .nunique().rename("n_peptides_credible"))
        genes = matches.groupby("sv_id")["ref_gene"].agg(
            lambda s: ";".join(sorted(set(s.dropna().astype(str))))).rename("catalogue_genes")
        for extra in (matched, credible, genes):
            table = table.merge(extra, left_on="sample_sv_id", right_index=True, how="left")
    for column in ("n_peptides_matched", "n_peptides_credible"):
        table[column] = table.get(column, pd.Series(dtype=float)).fillna(0).astype(int)

    # --- was this junction seen in patients? ------------------------------
    # Three independent levels, from strongest to weakest. They are reported
    # separately AND summarised, because they mean different things: an identical
    # peptide is a shared consequence, a shared gene and SV type is a shared
    # mechanism, and a nearby breakpoint is a shared locus — which at a fragile
    # site may be no more than a shared propensity to break.
    recurrence = _read(os.path.join(run_dir, "stage3_junction_recurrence.tsv"))
    if not recurrence.empty:
        recurrence["sample_sv_id"] = recurrence["sample_sv_id"].astype(str)
        # Everything is read as text; distances must be numeric before they can
        # be compared against a threshold.
        recurrence = _numeric(recurrence, ["patient_bp_dist_bp", "patient_bp_dist_bp1",
                                           "patient_bp_dist_bp2"])
        table = table.merge(recurrence, on="sample_sv_id", how="left")

    level2 = _read(os.path.join(run_dir, "stage3_gene_svtype.tsv"))
    if not level2.empty and "sv_id" in level2.columns:
        hit = level2.copy()
        hit["sample_sv_id"] = hit["sv_id"].astype(str)
        flag = (hit.groupby("sample_sv_id")
                .apply(lambda g: bool(_flag(g, "svtype_match").any()), include_groups=False)
                .rename("patient_same_gene_and_svtype").reset_index())
        genes = (hit.groupby("sample_sv_id")["ref_gene"]
                 .agg(lambda s: ";".join(sorted(set(s.dropna().astype(str)))[:5]))
                 .rename("patient_same_gene").reset_index()) \
            if "ref_gene" in hit.columns else None
        table = table.merge(flag, on="sample_sv_id", how="left")
        if genes is not None:
            table = table.merge(genes, on="sample_sv_id", how="left")
    table["patient_same_gene_and_svtype"] = table.get(
        "patient_same_gene_and_svtype", pd.Series(dtype=object)).fillna(False)

    def evidence(row) -> str:
        """The strongest level of patient recurrence this junction reaches."""
        if row.get("n_peptides_matched", 0) > 0:
            return "identical_peptide"
        if bool(row.get("patient_same_gene_and_svtype", False)):
            return "same_gene_and_svtype"
        if pd.notna(row.get("patient_same_gene")):
            return "same_gene"
        distance = pd.to_numeric(row.get("patient_bp_dist_bp"), errors="coerce")
        if pd.notna(distance):
            for limit, label in ((1_000, "breakpoint_within_1kb"),
                                 (10_000, "breakpoint_within_10kb"),
                                 (100_000, "breakpoint_within_100kb")):
                if distance <= limit:
                    return label
        return "not_seen_in_patients"
    table["patient_evidence"] = table.apply(evidence, axis=1)

    # Evidence for stages 6-8. Prefer the full-call-set pass when it exists:
    # the per-stage files cover only events that reached them, which leaves the
    # population and RNA columns empty for the ~99% of junctions that produced no
    # match — and a table meant for judging the call set cannot be mostly blank.
    full = _read(os.path.join(run_dir, "stage6_8_all_junctions.tsv"))
    sources = [full] if not full.empty else [
        _read(os.path.join(run_dir, f)) for f in
        ("credible_events.tsv", "stage7_expression.tsv", "stage8_rna_evidence.tsv")]

    # Per-event MHC counts, when binding was run. Only the binding columns are
    # taken: the rest of that file is credible_events.tsv, already merged above,
    # and re-merging it would collide column for column.
    events_binding = _read(os.path.join(run_dir, "events_binding.tsv"))
    if not events_binding.empty and "sample_sv_id" in events_binding.columns:
        keep = ["sample_sv_id"] + [c for c in events_binding.columns
                                   if c.startswith(("n_peptides_binding",
                                                    "n_peptides_unevaluable",
                                                    "binds_"))]
        sources = sources + [events_binding.loc[:, keep]]

    for extra in sources:
        if extra.empty or "sample_sv_id" not in extra.columns:
            continue
        extra = extra.drop(columns=[c for c in ("gene", "n_peptides", "peptides",
                                                "ref_genes", "junction_key")
                                    if c in extra.columns])
        extra = extra.loc[:, ~extra.columns.duplicated()].drop_duplicates("sample_sv_id")
        extra["sample_sv_id"] = extra["sample_sv_id"].astype(str)
        overlap = [c for c in extra.columns
                   if c != "sample_sv_id" and c in table.columns]
        extra = extra.drop(columns=overlap)
        table = table.merge(extra, on="sample_sv_id", how="left")

    drop = [c for c in table.columns
            if c.startswith("pass_") or c in ("is_private", "sv_hc")]
    return _numeric(table.drop(columns=drop), NUMERIC)


def column_dictionary(table: pd.DataFrame, gates) -> pd.DataFrame:
    """What each column means, where it came from, and how populated it is.

    The meanings live in `tools/column_meanings.py`, one written entry per
    column. A column with no entry is reported as `UNDOCUMENTED` rather than
    left blank, so it is visible in the artefact itself instead of being
    mistaken for a column whose meaning is obvious.
    """
    gate_help = dict(gates)
    origins = [
        (lambda c: c in ("sample", "cell_line", "branch", "sample_sv_id",
                         "junction_key", "prefix", "sv_id"), "identity"),
        (lambda c: c.startswith(("ref_", "or_", "hla_pres", "confirmed_", "is_cfs",
                                 "any_top50", "significant_", "unique_breaks",
                                 "compat_", "incompat_", "present_not_")),
         "reference catalogue"),
        (lambda c: c.startswith("patient_"), "stage 3 — patient cross"),
        (lambda c: c.startswith("gnomad_") or c in ("pon_count", "pon_fraction",
                                                    "panel_size_estimate",
                                                    "pass_pon", "pass_gnomad",
                                                    "is_private", "privacy_note",
                                                    "pon_in_privacy"),
         "stage 6 — population frequency"),
        (lambda c: c in ("sv_hc",) or c.endswith(("_bp1", "_bp2")) and
         any(c.startswith(k) for k in ("vf", "sf", "df", "ref", "qual", "segmapq",
                                       "homseq", "imprecise", "alt")),
         "SV call"),
        (lambda c: c.endswith("TPM") or c in ("expressed", "min_side_TPM",
                                              "retained_intron", "isofox_fusion",
                                              "isofox_fusion_support")
         or c.startswith(("alt_sj", "nearest_alt_sj")), "stage 7 — expression"),
        (lambda c: c in ("test", "test_reason", "junction_reads", "rna_tier")
         or c.startswith(("coverage_", "min_coverage", "softclip_", "junction_by_",
                          "alignments_", "low_mapq_")), "stage 8 — junction reads"),
        # MHC binding is its own layer, run after the pipeline proper, so it gets
        # its own label rather than falling through to "call and annotation".
        # Placed BEFORE the n_peptides rule, which would otherwise swallow
        # n_peptides_binding_*.
        (lambda c: c.startswith(("binds_", "panel_", "autologous_",
                                 "n_peptides_binding", "n_peptides_unevaluable",
                                 "n_patients_typed", "n_patients_presenting",
                                 "n_patients_unevaluable", "n_alleles_panel",
                                 "presenting_alleles")),
         "MHC binding — optional layer, not part of the recurrence test"),
        (lambda c: c.startswith("n_peptides") or c in ("gene_concordant",
                                                       "low_complexity", "is_self",
                                                       "credible",
                                                       "catalogue_genes"),
         "stages 3-5 — cross and QC"),
        (lambda c: True, "stage 1-2 — call and annotation"),
    ]

    rows = []
    for column in table.columns:
        origin = next(label for test, label in origins if test(column))
        meaning = column_meanings.resolve(column)
        if meaning is None:
            meaning = gate_help.get(column, "UNDOCUMENTED — add an entry to "
                                            "tools/column_meanings.py")
        rows.append({"column": column, "origin": origin, "meaning": meaning,
                     "non_null": int(table[column].notna().sum()),
                     "pct_populated": round(
                         100 * table[column].notna().sum() / max(len(table), 1), 1),
                     "example": next((str(v) for v in table[column].dropna()[:1]), "")})
    return pd.DataFrame(rows)


def identify(table, label: str) -> None:
    """Insert the run identifier, split into the two things it conflates.

    `sample` is the run directory, which is how a row is traced back to the
    outputs and the criteria manifest that produced it. But reports name the
    cell line alone, so cross-referencing a table row against a report meant
    reading past a suffix. Both are now columns:

        sample     parental_noPON     the run — unique, matches results/<dir>
        cell_line  parental           the biological sample, as reports name it
        branch     noPON             which thresholds produced this row

    Splitting them also makes the table groupable by line across branches, which
    a single concatenated string does not allow.
    """
    cell_line, _, branch = label.partition("_")
    table.insert(0, "branch", branch or "default")
    table.insert(0, "cell_line", cell_line)
    table.insert(0, "sample", label)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-dir", default="results")
    ap.add_argument("--runs", nargs="*", metavar="RUN",
                    help="restrict the COMBINED table to these run directories "
                         "(e.g. parental_noPON). Per-run tables are still written "
                         "for every run.")
    ap.add_argument("--suffix", default="",
                    help="suffix for the combined filenames, so a restricted "
                         "table cannot overwrite the full one "
                         "(e.g. --suffix _noPON -> master_sv_noPON.tsv)")
    args = ap.parse_args()

    available = sorted(
        e for e in os.listdir(args.results_dir)
        if os.path.exists(os.path.join(args.results_dir, e, "summary.json")))
    if args.runs:
        # A run name that matches nothing must fail loudly. Silently combining
        # fewer runs than asked for produces a table that looks complete and is
        # not — and this table is the one shared as the record of a report.
        unknown = [r for r in args.runs if r not in available]
        if unknown:
            raise SystemExit(
                "these runs do not exist in " + args.results_dir + ":\n  "
                + "\n  ".join(unknown) + "\n\navailable:\n  "
                + "\n  ".join(available))

    peptide_tables, sv_tables = [], []
    for entry in sorted(os.listdir(args.results_dir)):
        run_dir = os.path.join(args.results_dir, entry)
        if not os.path.isdir(run_dir) or \
                not os.path.exists(os.path.join(run_dir, "summary.json")):
            continue      # e.g. reports/ — a directory, but not a run

        peptides = build_peptide_table(run_dir, entry)
        svs = build_sv_table(run_dir, entry)
        # Per-run tables are always written; only the combined one is filtered.
        in_scope = (not args.runs) or entry in args.runs
        if not peptides.empty:
            peptides.to_csv(os.path.join(run_dir, "master_peptides.tsv"),
                            sep="\t", index=False)
            if in_scope:
                peptide_tables.append(peptides)
        if not svs.empty:
            svs.to_csv(os.path.join(run_dir, "master_sv.tsv"), sep="\t", index=False)
            if in_scope:
                sv_tables.append(svs)
        written = [n for n, frame in (("master_sv.tsv", svs),
                                      ("master_peptides.tsv", peptides))
                   if not frame.empty]
        print(f"  {entry:26} SVs {len(svs):>5}  matched peptides {len(peptides):>4}"
              f"   -> {', '.join(written) or 'nothing (no admitted junctions)'}")

    for tables, name, gates in ((peptide_tables, "master_peptides", PEPTIDE_GATES),
                                (sv_tables, "master_sv", SV_GATES)):
        if not tables:
            continue
        combined = pd.concat(tables, ignore_index=True)
        path = os.path.join(args.results_dir, f"{name}{args.suffix}.tsv")
        combined.to_csv(path, sep="\t", index=False)
        column_dictionary(combined, gates).to_csv(
            path.replace(".tsv", "_column_dictionary.tsv"), sep="\t", index=False)
        scope = ", ".join(args.runs) if args.runs else "every run"
        print(f"\ncombined {name}{args.suffix}: {len(combined)} rows x "
              f"{len(combined.columns)} columns  [{scope}]")
        print(f"  {path}")
        print(f"  {path.replace('.tsv', '_column_dictionary.tsv')}")

        print(f"  attrition:")
        for column, description in gates:
            if column not in combined.columns:
                continue
            passed = int(combined[column].astype("boolean").fillna(False).sum())
            evaluated = int(combined[column].notna().sum())
            print(f"    {column:26} {passed:>5} / {evaluated:<5} evaluated  {description}")


if __name__ == "__main__":
    main()
