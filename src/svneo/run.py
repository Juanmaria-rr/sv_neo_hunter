#!/usr/bin/env python3
"""
run.py — the orchestrator: every sample, every branch, same stages, same order.

    python -m svneo.run --config config/my_run.yaml --out-dir results
    python -m svneo.run --config ... --samples parental --branches PON10
    python -m svneo.run --config ... --dry-run          # plan only, no compute

WHAT THIS GUARANTEES
--------------------
1. Parents run before children, so a hit in a parent is attributable to the
   parent (config.ordered_samples()).
2. Every sample gets the same stages with the same thresholds; a sample lacking a
   modality gets `NA` for those stages and continues. Missing data is never
   recorded as a negative result.
3. Every output directory contains the criteria manifest that produced it, so no
   number can outlive the thresholds it came from. This rule exists because a
   previous analysis published conclusions its own corrected code no longer
   reproduced.
4. Failures are per sample, not per run: one sample's crash does not discard the
   others' results.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import traceback

import pandas as pd

from . import config as config_mod
from . import confidence, criteria, cross, generators, null_model, rna, synthesis, vcf


def _write(frame: pd.DataFrame, path: str) -> None:
    if frame is None or (hasattr(frame, "empty") and frame.empty):
        return
    frame.to_csv(path, sep="\t", index=False)


def run_sample(sample, branch, cfg, reference, proteome, pool,
               out_root: str, breakpoint_indices: dict) -> dict:
    """One sample x one branch, all stages. Returns its counts and artefacts."""
    tag = f"{sample.name}" + (f"_{branch.name}" if branch.name != "default" else "")
    out_dir = os.path.join(out_root, tag)
    os.makedirs(out_dir, exist_ok=True)
    counts: dict = {}
    print(f"\n=== {tag}  [{sample.vcf_kind}] ===")

    # -- stage 1: admission ------------------------------------------------
    breakends = vcf.read_breakends(sample.vcf)
    counts["vcf_records"] = len(breakends)
    # The panel threshold only applies where the panel is an INFO field this
    # pipeline thresholds itself. On a somatic set the caller already decided,
    # and `admit_panel_filtered` is the only lever that can reopen it.
    pon_max = branch.pon_max if sample.vcf_kind == "germline" else None
    admitted, funnel1 = vcf.admit(breakends, sample.vcf_kind, pon_max,
                                  admit_panel_filtered=branch.admit_panel_filtered)
    counts["admitted"] = len(admitted)
    junctions = vcf.pair_junctions(admitted)
    counts["junctions"] = len(junctions)
    panel_size = vcf.panel_size_estimate(breakends)
    _write(junctions, os.path.join(out_dir, "stage1_junctions.tsv"))
    print(f"  stage1: {counts['vcf_records']:,} records -> "
          f"{counts['admitted']:,} admitted -> {counts['junctions']:,} junctions")

    breakpoint_indices[sample.name] = synthesis.breakpoint_index(junctions)

    # The generator applies no admission logic of its own, so it must be handed
    # the ADMITTED call set. Feeding it the raw VCF would build peptides from PON
    # and INFERRED breakends and leave artefact-derived peptides in the candidate
    # universe, where they can match a reference peptide by chance.
    admitted_vcf = os.path.join(out_dir, f"{sample.name}.admitted.vcf")
    n_written = vcf.write_admitted_vcf(sample.vcf, admitted, admitted_vcf)
    print(f"  stage1: wrote {n_written:,} admitted records -> "
          f"{os.path.basename(admitted_vcf)}")

    # -- stage 2: candidate peptides (the one pluggable stage) -------------
    # The backend is chosen by config (`peptide_generator:`); nothing else in the
    # pipeline imports a peptide-generation tool. See generators/base.py for the
    # contract every backend must satisfy.
    prefix = os.path.join(out_dir, sample.name)
    stored = cfg.resources.get("peptides_from", {}).get(sample.name)
    backend_name = "precomputed" if stored else cfg.peptide_generator
    try:
        backend = generators.get_generator(
            backend_name, **({"source_prefix": stored} if stored
                             else {"neosv_path": cfg.resources.get("neosv_path")}))
        generated = backend.generate(
            admitted_vcf, out_dir, sample.name,
            release=cfg.ensembl_release,
            cache_dir=cfg.resources.get("pyensembl_cache"))
    except generators.GeneratorError as error:
        print(f"  stage2 SKIPPED: {error}")
        counts["candidate_peptides"] = "NA"
        return {"counts": counts, "out_dir": out_dir, "error": str(error)}

    peptide_table, anno = generated.peptides, generated.annotation
    counts["candidate_peptides"] = int(peptide_table["neopeptide"].nunique())
    print(f"  stage2 [{generated.generator}]: {len(peptide_table):,} peptide rows "
          f"-> {counts['candidate_peptides']:,} unique   {generated.stats}")

    # -- stages 3-5: cross, QC, events -------------------------------------
    matches = cross.level1_identical(peptide_table, anno, reference)
    # Rows are (peptide x SV x transcript): one peptide can be produced by several
    # SVs and appear more than once. Both grains are reported because they answer
    # different questions and mixing them makes runs incomparable — the unique
    # count is the one to quote against a catalogue, and the one the null model
    # tests.
    counts["matches_identical"] = len(matches)
    counts["matches_unique_peptides"] = int(matches["peptide"].nunique()) \
        if len(matches) else 0
    if not matches.empty:
        counts["matches_gene_concordant"] = int(
            matches["gene_concordant"].fillna(False).sum())
        matches = cross.sequence_qc(matches, proteome)
        counts["matches_credible"] = int(matches["credible"].fillna(False).sum())
        _write(matches, os.path.join(out_dir, "stage3_matches.tsv"))

    level2 = cross.level2_gene_svtype(anno, reference)
    _write(level2, os.path.join(out_dir, "stage3_gene_svtype.tsv"))
    if pool is not None and not pool.empty:
        _write(cross.level3_proximity(junctions, pool),
               os.path.join(out_dir, "stage3_proximity.tsv"))
        # The same question keyed on OUR junctions, so every admitted junction
        # can state whether anything like it was seen in patients.
        _write(cross.junction_recurrence(junctions, pool),
               os.path.join(out_dir, "stage3_junction_recurrence.tsv"))

    events = cross.to_events(matches, junctions) if not matches.empty else pd.DataFrame()
    counts["events"] = len(events)
    print(f"  stages3-5: {counts['matches_identical']} matches -> "
          f"{counts.get('matches_credible', 0)} credible -> {counts['events']} events")

    # -- stage 6: confidence and privacy -----------------------------------
    if not events.empty:
        events = confidence.annotate_confidence(events)
        af = confidence.annotate_gnomad(events, cfg.resources.get("gnomad_sv", ""),
                                        cfg.resources.get("gnomad_helper"))
        events = confidence.annotate_privacy(events, panel_size, af,
                                             pon_in_privacy=branch.pon_in_privacy)
        counts["events_hc"] = int(events["sv_hc"].sum())
        counts["events_private"] = int(events["is_private"].sum())
        basis = "PON+population" if branch.pon_in_privacy else "population only"
        print(f"  stage6: {counts['events_hc']} HC, "
              f"{counts['events_private']} private ({basis})")

    # -- stage 7: expression -----------------------------------------------
    if not events.empty and sample.has_expression:
        expr, background = rna.expression(events, anno, sample.isofox_dir,
                                          sample.isofox_prefix)
        # Isofox's own calls — splice junctions, fusions, retained introns — are
        # independent of our read counting and were previously unused.
        context = rna.isofox_context(events, anno, sample.isofox_dir,
                                     sample.isofox_prefix)
        if not context.empty:
            expr = expr.merge(context, on="sample_sv_id", how="left")
        _write(expr, os.path.join(out_dir, "stage7_expression.tsv"))
        counts["events_expressed"] = int(expr["expressed"].sum())
        with open(os.path.join(out_dir, "stage7_background.json"), "w") as fh:
            json.dump(background, fh, indent=2)
    else:
        counts["events_expressed"] = "NA"

    # -- stages 6-8 over EVERY admitted junction ---------------------------
    # The verdict stages above run on credible events only. This pass measures
    # the same evidence for the whole call set, so the per-SV table can be used
    # to judge junctions that never produced a match — otherwise its RNA and
    # population columns are empty for 99% of rows.
    if criteria.EVALUATE_ALL_JUNCTIONS and not junctions.empty:
        every = junctions.copy()
        every["gene"] = None
        print(f"  stages6-8: measuring all {len(every):,} admitted junctions")
        every = confidence.annotate_confidence(every)
        af_all = confidence.annotate_gnomad(every, cfg.resources.get("gnomad_sv", ""),
                                            cfg.resources.get("gnomad_helper"))
        every = confidence.annotate_privacy(every, panel_size, af_all,
                                            pon_in_privacy=branch.pon_in_privacy)
        if sample.has_expression:
            expr_all, _ = rna.expression(every, anno, sample.isofox_dir,
                                         sample.isofox_prefix)
            context_all = rna.isofox_context(every, anno, sample.isofox_dir,
                                             sample.isofox_prefix)
            if not context_all.empty:
                expr_all = expr_all.merge(context_all, on="sample_sv_id", how="left")
            every = every.merge(expr_all.drop(columns=["gene"], errors="ignore"),
                                on="sample_sv_id", how="left", suffixes=("", "_dup"))
        evidence_all = rna.evaluate(every, sample.rna_bam if sample.has_rna else None)
        every = every.merge(evidence_all.drop(columns=["gene"], errors="ignore"),
                            on="sample_sv_id", how="left", suffixes=("", "_dup"))
        every = every.drop(columns=[c for c in every.columns if c.endswith("_dup")])
        _write(every, os.path.join(out_dir, "stage6_8_all_junctions.tsv"))

    # -- stage 8: junction evidence ----------------------------------------
    if not events.empty:
        evidence = rna.evaluate(events, sample.rna_bam if sample.has_rna else None)
        _write(evidence, os.path.join(out_dir, "stage8_rna_evidence.tsv"))
        if sample.has_rna:
            supported = evidence["rna_tier"].isin(["STRONG", "SUGGESTIVE", "WEAK"])
            counts["events_rna_supported"] = int(supported.sum())
            merged = events.merge(evidence[["sample_sv_id", "rna_tier"]],
                                  on="sample_sv_id", how="left")
            counts["events_hc_and_rna"] = int(
                (merged["sv_hc"] & merged["rna_tier"].isin(["STRONG", "SUGGESTIVE"])).sum())
            # Privacy is NOT in the count above, and a report that says "satisfies
            # every criterion" over it overstates the result: an event can be
            # high-confidence and transcribed while still being common in a
            # population. The all-criteria count is therefore computed and
            # reported separately.
            counts["events_private_hc_and_rna"] = int(
                (merged["is_private"] & merged["sv_hc"]
                 & merged["rna_tier"].isin(["STRONG", "SUGGESTIVE"])).sum())
            print(f"  stage8: tiers {evidence['rna_tier'].value_counts().to_dict()}")
        else:
            counts["events_rna_supported"] = "NA"
            counts["events_hc_and_rna"] = "NA"
            counts["events_private_hc_and_rna"] = "NA"
        _write(events, os.path.join(out_dir, "credible_events.tsv"))

    # -- null model ---------------------------------------------------------
    null = null_model.permutation_test(
        peptide_table[criteria.PEPTIDE_COLUMN].dropna().unique(),
        reference["ref_peptide"].dropna().unique(),
        n_permutations=criteria.NULL_PERMUTATIONS)
    print(f"  null: {null_model.summarise(null)}")

    # -- funnel, notes, manifest -------------------------------------------
    _write(synthesis.build_funnel(counts), os.path.join(out_dir, "funnel.tsv"))
    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        json.dump({"sample": sample.name, "branch": branch.name,
                   "vcf_kind": sample.vcf_kind, "counts": counts,
                   "stage1_funnel": funnel1, "null_model": null,
                   "modalities": sample.modalities(),
                   "peptide_generator": generated.manifest(),
                   "notes": synthesis.interpret(counts),
                   "criteria": criteria.manifest()}, fh, indent=2, default=str)

    for note in synthesis.interpret(counts):
        print(f"  NOTE: {note}")
    return {"counts": counts, "null": null, "events": events, "out_dir": out_dir}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-dir", default=None, help="overrides the config value")
    ap.add_argument("--samples", nargs="*", help="subset by name")
    ap.add_argument("--branches", nargs="*", help="subset by name")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the plan and exit without computing")
    args = ap.parse_args()

    try:
        cfg = config_mod.load(args.config)
    except (config_mod.ConfigError, FileNotFoundError) as error:
        # A misconfiguration is the user's to fix, so print what is wrong and
        # nothing else. A stack trace here only obscures the message.
        sys.exit(f"config error in {args.config}:\n{error}")
    print(config_mod.describe(cfg))
    if cfg.criteria_overrides:
        applied = criteria.apply_overrides(cfg.criteria_overrides)
        print(f"criteria overrides: {applied}")

    out_root = args.out_dir or cfg.out_dir
    # A name that matches nothing must fail loudly. Selecting zero samples and
    # exiting 0 looks exactly like a successful run, and sample names use hyphens
    # while the directories they came from use underscores — an easy mismatch.
    def select(available, requested, what):
        if not requested:
            return list(available)
        names = {getattr(x, "name") for x in available}
        unknown = [r for r in requested if r not in names]
        if unknown:
            raise SystemExit(
                f"error: no {what} named {unknown} in the config.\n"
                f"       available: {sorted(names)}\n"
                f"       (sample names use the config's spelling, which may differ "
                f"from directory names)")
        return [x for x in available if x.name in requested]

    samples = select(cfg.ordered_samples(), args.samples, "sample")
    branches = select(cfg.branches, args.branches, "branch")

    if args.dry_run:
        print(f"\nDRY RUN — would process {len(samples)} sample(s) x "
              f"{len(branches)} branch(es) into {out_root}/")
        return

    os.makedirs(out_root, exist_ok=True)
    reference = cross.load_reference(cfg.reference)
    proteome = (cross.load_proteome(cfg.reference.proteome)
                if cfg.reference.proteome else None)
    if proteome is None:
        print("NOTE: no proteome configured — the self-peptide test is recorded "
              "as NA, not as passing.")
    pool = None
    if cfg.reference.pool and os.path.exists(cfg.reference.pool):
        pool = pd.read_csv(cfg.reference.pool, sep="\t", dtype=str)
        pool = pool.rename(columns={cfg.reference.pool_peptide_column: "ref_peptide"})

    results, indices = {}, {}
    for branch in branches:
        for sample in samples:
            key = f"{sample.name}" + (f"_{branch.name}" if branch.name != "default" else "")
            try:
                results[key] = run_sample(sample, branch, cfg, reference, proteome,
                                          pool, out_root, indices)
            except Exception:                      # one sample must not sink the run
                print(f"  FAILED: {key}", file=sys.stderr)
                traceback.print_exc()
                results[key] = {"counts": {}, "error": "see traceback"}

    # -- cross-sample synthesis --------------------------------------------
    _write(synthesis.compare_samples(results),
           os.path.join(out_root, "comparison.tsv"))
    _write(null_model.rate_table({k: v.get("null", {}) for k, v in results.items()}),
           os.path.join(out_root, "null_model_rates.tsv"))

    attributions = []
    for key, result in results.items():
        events = result.get("events")
        if events is None or events.empty:
            continue
        sample_name = key.split("_")[0] if "_" in key else key
        sample = cfg.by_name(sample_name) or cfg.by_name(key)
        if sample:
            attributions.append(synthesis.attribute(events, sample.name, cfg, indices))
    if attributions:
        _write(pd.concat(attributions, ignore_index=True),
               os.path.join(out_root, "attribution.tsv"))

    print(f"\nwrote {out_root}/comparison.tsv, null_model_rates.tsv"
          + (", attribution.tsv" if attributions else ""))
    print("Read comparison.tsv on EVENTS and per-1,000 rates, never on peptide "
          "counts: one locus seen through a sliding window can produce dozens.")


if __name__ == "__main__":
    main()
