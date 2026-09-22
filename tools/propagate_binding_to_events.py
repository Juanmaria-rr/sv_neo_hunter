#!/usr/bin/env python
"""
propagate_binding_to_events.py — carry MHC binding from peptides up to events.

WHY THIS EXISTS
---------------
A cross reports distinct recurrent peptides; a published recurrence count from a
presentation-filtered analysis is not the same kind of object, because one has
already been filtered for MHC presentation and the other has not. Running
netMHCpan on the unfiltered side is what makes them comparable -- but only if the
*same question* is asked of both.

TWO QUESTIONS, ONLY ONE OF WHICH IS COMPARABLE
----------------------------------------------
`tools/run_netmhcpan.py` predicts every matched peptide against the whole cohort
allele panel. That answers "could this peptide be presented by *someone* in the
cohort" -- a union over hundreds of alleles.

An analysis that filters each
peptide against the HLA genotype of the *patient in whom it was observed* --
~6 alleles. That answers "was it presented where we actually saw it", and it is
the question that produced the 130.

So this script computes both and keeps them apart:

  binds_panel        peptide binds >= 1 allele of the cohort panel.
                     UPPER BOUND. Not comparable with the earlier analysis.
  binds_autologous   peptide binds >= 1 allele of at least one HLA-typed patient
                     who actually carries that peptide.
                     THIS is the comparable quantity; quote this one.

The binder definition itself is imported from `run_netmhcpan.py` rather than
restated, so the three-way cut (IC50 <= 500 nM AND %Rank_BA <= 2 AND
%Rank_EL <= 2) cannot drift between the two scripts.

COVERAGE IS NOT COMPLETE, AND SILENCE IS NOT A NEGATIVE
-------------------------------------------------------
Not every catalogue patient carries a linkable HLA genotype. A peptide
contributed only by untyped patients cannot
be evaluated autologously -- it is reported as `unevaluable`, never as a
non-binder, because those are different facts. The same applies to a patient
carrying an allele outside the predicted panel.

Usage
-----
    python tools/propagate_binding_to_events.py \
        --sample parental_noPON
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from run_netmhcpan import BINDER_THRESHOLDS      # noqa: E402  single source of truth

REPO = pathlib.Path(__file__).resolve().parent.parent
#: The cross outputs and the patient catalogue live beside the repository, not
#: inside it: they are derived from patient genotypes and must stay unpublished.
DEFAULT_ROOT = REPO.parent

#: RNA tiers. `events_rna_supported` admits WEAK; the conjunction rows do not.
#: Copied deliberately from src/svneo/run.py -- if that changes, this must too.
RNA_SUPPORTED = {"STRONG", "SUGGESTIVE", "WEAK"}
RNA_CONJUNCTION = {"STRONG", "SUGGESTIVE"}


def normalise_allele(raw: str) -> str:
    """The three allele spellings in play, reduced to one.

    They genuinely differ and nothing warns you: `reference/cohort_observed_alleles.txt`
    holds netMHCpan's *input* form with no asterisk (`HLA-A01:01`, because
    netMHCpan returns zero rows without erroring if you feed it one), while both
    netMHCpan's *output* table and `<cohort>/HLAs/*.hla.txt` carry it
    (`HLA-A*01:01`). Comparing the two forms directly yields an empty
    intersection and therefore a clean, entirely wrong set of zeros.
    """
    return raw.strip().upper().replace(" ", "").replace("*", "")


#: Selectable binder definitions. Each is a set of `column <= threshold` cuts
#: that must ALL hold. `three_way` is the earlier analysis's own rule and the
#: default; nothing else is comparable with its numbers.
#:
#: The `neo_*` entries use netMHCpan 4.2's CEDAR neoepitope-finetuned head, which
#: has NO established cut — these are netMHCpan's generic strong/weak rank
#: conventions applied to it, which is a choice, not a standard. Measured on this
#: conventions applied to it. `neo_strong` is of comparable stringency to the
#: three-way rule; `neo_weak` is several times more permissive.
CRITERIA = {
    # Imported, never restated: the three-way rule has one definition, in
    # run_netmhcpan.py, so it cannot drift between the run and its propagation.
    "three_way": BINDER_THRESHOLDS,
    "neo_strong": {"rank_Neo": 0.5},
    "neo_weak": {"rank_Neo": 2.0},
}

#: How far inside a cut a binder sits, as a fraction of that cut. 0 means exactly
#: on the boundary, 1 means the value is 0 (or, for affinity, essentially so).
#: Scales are mixed on purpose -- nM against percentile -- so this is a usable
#: heuristic, NOT a calibrated probability. What it is good for: separating calls
#: that survive any perturbation from calls that a re-run, a version bump or a
#: threshold nudge would flip. Measured motivation: 95% of the calls that
#: disagreed between netMHCpan 4.1 and 4.2 sat within 50% of a threshold.
def margin_of(values: dict, cuts: dict) -> tuple[float, str]:
    """(smallest relative distance inside any cut, which cut that was)."""
    distances = {column: (threshold - values[column]) / threshold
                 for column, threshold in cuts.items()}
    cut = min(distances, key=distances.get)
    return distances[cut], cut


def read_binding_alleles(path: pathlib.Path,
                         cuts: dict) -> dict[str, dict[str, dict]]:
    """peptide -> {allele: {column: value}} for rows passing every cut.

    Streamed: the prediction table is ~7.3M rows / 461-578 MB, of which well
    under 2% pass, so only the survivors are ever held in memory. The numbers are
    kept, not just the allele names, because the margin cannot be recovered from
    a set membership.
    """
    binding: dict[str, dict[str, dict]] = collections.defaultdict(dict)
    seen = 0
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in cuts if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(
                f"{path.name} has no column(s) {', '.join(missing)} — that "
                "criterion needs a prediction table produced with the matching "
                "flags (`run_netmhcpan.py --neo` for the rank_Neo columns)")
        for row in reader:
            seen += 1
            try:
                values = {column: float(row[column]) for column in cuts}
            except (ValueError, KeyError, TypeError):
                continue
            if all(values[c] <= t for c, t in cuts.items()):
                binding[row["peptide"]][normalise_allele(row["allele"])] = values
    print(f"  netMHCpan: {seen:,} predictions, "
          f"{len(binding):,} peptides bind >= 1 panel allele")
    return dict(binding)


def read_genotype(path: pathlib.Path) -> set[str]:
    return {normalise_allele(line) for line in path.read_text().splitlines()
            if line.strip()}


def autologous_pass(peptides_of_interest: set[str],
                    binding: dict[str, dict[str, dict]],
                    catalogue_dir: pathlib.Path,
                    hla_dir: pathlib.Path,
                    panel: set[str],
                    cuts: dict) -> dict[str, dict]:
    """One pass over the HLA-typed patients, accumulating per peptide.

    Deliberately incremental: a peptide's contributing patients are never all
    held at once. Each patient's genotype is intersected with the already-known
    binding alleles and the result folded in, so memory stays O(peptides).
    """
    stats = {p: {"n_patients_typed": 0, "n_patients_presenting": 0,
                 "n_patients_unevaluable": 0, "presenting_alleles": set(),
                 "best": None}
             for p in peptides_of_interest}

    genotypes = sorted(hla_dir.glob("*.hla.txt"))
    print(f"  patients with HLA: {len(genotypes):,}")
    linked = 0
    for n, hla_path in enumerate(genotypes, 1):
        patient = hla_path.name.split(".")[0]
        peptide_path = catalogue_dir / f"{patient}.peptides.tsv"
        if not peptide_path.exists():
            continue
        linked += 1
        genotype = read_genotype(hla_path)
        # An allele we never predicted cannot be used to call a non-binder.
        blind = bool(genotype - panel)

        with peptide_path.open(newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            next(reader, None)                       # header
            for row in reader:
                if not row:
                    continue
                peptide = row[0]
                entry = stats.get(peptide)
                if entry is None:
                    continue
                entry["n_patients_typed"] += 1
                predicted = binding.get(peptide) or {}
                hit = set(predicted) & genotype
                if hit:
                    entry["n_patients_presenting"] += 1
                    entry["presenting_alleles"] |= hit
                    # A peptide is as robust as its most comfortable presenting
                    # allele, so keep the maximum margin seen, not the first.
                    for allele in hit:
                        values = predicted[allele]
                        margin, cut = margin_of(values, cuts)
                        if entry["best"] is None or margin > entry["best"][0]:
                            entry["best"] = (margin, cut, allele, values)
                elif blind:
                    entry["n_patients_unevaluable"] += 1
        if n % 500 == 0:
            print(f"    {n:,}/{len(genotypes):,} patients")
    print(f"  patients linked to a peptide file: {linked:,}")
    return stats


def classify(entry: dict, binds_panel: bool) -> str:
    """The autologous verdict, with 'unevaluable' kept distinct from 'no'."""
    if entry["n_patients_presenting"]:
        return "binder"
    if entry["n_patients_typed"] == 0:
        return "unevaluable_no_typed_patient"
    if entry["n_patients_unevaluable"] and binds_panel:
        return "unevaluable_allele_outside_panel"
    return "non_binder"


def funnel(events: list[dict], rna: dict[str, str], keep) -> dict[str, int]:
    """The event rows of the standard funnel, restricted by `keep`."""
    kept = [e for e in events if keep(e)]

    def flag(event, column):
        return str(event.get(column, "")).strip().lower() == "true"

    tier = lambda e: rna.get(e["sample_sv_id"], "")          # noqa: E731
    return {
        "events": len(kept),
        "events_private": sum(flag(e, "is_private") for e in kept),
        "events_hc": sum(flag(e, "sv_hc") for e in kept),
        "events_rna_supported": sum(tier(e) in RNA_SUPPORTED for e in kept),
        "events_hc_and_rna": sum(flag(e, "sv_hc") and tier(e) in RNA_CONJUNCTION
                                 for e in kept),
        "events_private_hc_and_rna": sum(
            flag(e, "is_private") and flag(e, "sv_hc")
            and tier(e) in RNA_CONJUNCTION for e in kept),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", default="parental_noPON")
    parser.add_argument("--root", type=pathlib.Path, default=DEFAULT_ROOT,
                        help="directory holding results_cross/, "
                             "the catalogue peptide files and the genotypes")
    parser.add_argument("--cross-dir", type=pathlib.Path,
                        help="overrides <root>/results_cross/<sample>")
    parser.add_argument("--catalogue-dir", type=pathlib.Path,
                        help="one peptide file per catalogue patient")
    parser.add_argument("--hla-dir", type=pathlib.Path,
                        help="one class I genotype file per typed patient")
    parser.add_argument("--panel", type=pathlib.Path,
                        default=REPO / "reference" / "cohort_observed_alleles.txt")
    parser.add_argument("--criterion", choices=sorted(CRITERIA),
                        default="three_way",
                        help="binder definition. `three_way` is the earlier "
                             "analysis's rule and the only one comparable with "
                             "its numbers; the neo_* rules use 4.2's CEDAR head, "
                             "for which no cut is established.")
    parser.add_argument("--predictions", type=pathlib.Path,
                        help="prediction table; defaults to netmhcpan.tsv, or "
                             "netmhcpan_neo.tsv for a neo_* criterion")
    parser.add_argument("--peptide-list", type=pathlib.Path,
                        help="the matched peptide set, one per line or a TSV "
                             "with a `peptide` column. Defaults to this "
                             "sample's netmhcpan_peptide_summary.tsv. Use it "
                             "when several samples share one prediction table — "
                             "a lineage's peptides overlap heavily, so "
                             "predicting each sample separately repeats work.")
    parser.add_argument("--label", default="",
                        help="suffix for output filenames, so a second criterion "
                             "does not overwrite the first")
    parser.add_argument("--out-dir", type=pathlib.Path,
                        help="defaults to the cross directory")
    args = parser.parse_args()

    cross = args.cross_dir or args.root / "results_cross" / args.sample
    catalogue = args.catalogue_dir or args.root / "catalogue_peptides"
    hla_dir = args.hla_dir or args.root / "cohort" / "HLAs"
    out_dir = args.out_dir or cross
    for path in (cross, catalogue, hla_dir):
        if not path.exists():
            raise SystemExit(f"missing input: {path}")

    cuts = CRITERIA[args.criterion]
    default_predictions = ("netmhcpan_neo.tsv" if args.criterion.startswith("neo")
                           else "netmhcpan.tsv")
    predictions = args.predictions or cross / default_predictions
    label = args.label or ("" if args.criterion == "three_way"
                           else f"_{args.criterion}")

    print(f"sample {args.sample}")
    print(f"  criterion: {args.criterion} — "
          + " AND ".join(f"{c} <= {t:g}" for c, t in cuts.items()))
    print(f"  predictions: {predictions.name}")
    binding = read_binding_alleles(predictions, cuts)

    panel = {normalise_allele(line) for line in
             args.panel.read_text().splitlines() if line.strip()
             and not line.lstrip().startswith("#")}
    # A spelling mismatch here would silently make every peptide a non-binder,
    # which is exactly the failure mode netMHCpan itself has with `*`.
    predicted = set().union(*binding.values()) if binding else set()
    if predicted and not (predicted & panel):
        raise SystemExit("no allele spelling in common between the prediction "
                         "table and the panel file -- refusing to report zeros")

    with (cross / "credible_events.tsv").open(newline="") as fh:
        events = list(csv.DictReader(fh, delimiter="\t"))
    with (cross / "stage8_rna_evidence.tsv").open(newline="") as fh:
        rna = {r["sample_sv_id"]: r.get("rna_tier", "")
               for r in csv.DictReader(fh, delimiter="\t")}

    event_peptides = {e["sample_sv_id"]: [p for p in e["peptides"].split(";") if p]
                      for e in events}
    credible_peptides = set().union(*event_peptides.values()) if event_peptides else set()

    # TWO PEPTIDE LEVELS, AND THEY ANSWER DIFFERENT QUESTIONS.
    # `matches_unique_peptides` is the number that stands opposite a published
    # recurrence count, so it is the one to filter for comparability.
    # The peptides carried by credible events (1,809) are a much smaller set,
    # downstream of gene concordance and the complexity/self cuts -- useful for
    # the event funnel, useless for the comparison. Both are reported.
    peptide_list = args.peptide_list or cross / "netmhcpan_peptide_summary.tsv"
    lines = [line.rstrip("\n") for line in
             peptide_list.read_text().splitlines() if line.strip()]
    header = lines[0].split("\t") if lines else []
    if "peptide" in header:                      # a TSV with a named column
        column = header.index("peptide")
        matched_peptides = {line.split("\t")[column] for line in lines[1:]}
    else:                                        # a bare list, no header
        matched_peptides = {line.split("\t")[0].strip() for line in lines}
    # The peptide set is the same whichever criterion is applied — it is what was
    # submitted for prediction, not what passed — so it is read from the original
    # summary regardless of `--criterion`.
    missing = credible_peptides - matched_peptides
    if missing:
        raise SystemExit(f"{len(missing):,} peptides in credible events were "
                         "never submitted to netMHCpan -- the two tables do not "
                         "describe the same run")
    print(f"  credible events: {len(events):,}, peptides in them: "
          f"{len(credible_peptides):,}; matched peptides overall: "
          f"{len(matched_peptides):,}")

    stats = autologous_pass(matched_peptides, binding, catalogue, hla_dir,
                            panel, cuts)
    peptides = matched_peptides

    # -- peptide level ------------------------------------------------------
    peptide_rows = []
    for peptide in sorted(peptides):
        entry = stats[peptide]
        alleles = binding.get(peptide) or {}
        verdict = classify(entry, bool(alleles))
        panel_best = max((margin_of(values, cuts) + (allele,)
                          for allele, values in alleles.items()),
                         default=None)
        best = entry["best"]
        peptide_rows.append({
            "peptide": peptide,
            "n_alleles_panel": len(alleles),
            "binds_panel": bool(alleles),
            "panel_margin": round(panel_best[0], 4) if panel_best else "",
            "panel_limiting_cut": panel_best[1] if panel_best else "",
            "panel_best_allele": panel_best[2] if panel_best else "",
            "n_patients_typed": entry["n_patients_typed"],
            "n_patients_presenting": entry["n_patients_presenting"],
            "n_patients_unevaluable": entry["n_patients_unevaluable"],
            "binds_autologous": verdict == "binder",
            "autologous_verdict": verdict,
            "presenting_alleles": ";".join(sorted(entry["presenting_alleles"])),
            "autologous_margin": round(best[0], 4) if best else "",
            "autologous_limiting_cut": best[1] if best else "",
            "autologous_best_allele": best[2] if best else "",
            **{f"autologous_{column}": (best[3][column] if best else "")
               for column in cuts},
            "autologous_robustness": (
                "" if not best else
                "flippable" if best[0] <= 0.10 else
                "marginal" if best[0] <= 0.25 else
                "solid" if best[0] <= 0.50 else "robust"),
        })
    _write(out_dir / f"peptide_binding{label}.tsv", peptide_rows)

    panel_binders = {r["peptide"] for r in peptide_rows if r["binds_panel"]}
    auto_binders = {r["peptide"] for r in peptide_rows if r["binds_autologous"]}
    unevaluable = {r["peptide"] for r in peptide_rows
                   if r["autologous_verdict"].startswith("unevaluable")}
    # The subset that survives any perturbation: comfortably inside every cut.
    # Not a stricter biological claim -- the same criterion, minus the calls that
    # a re-run or a version bump would flip.
    robust = {r["peptide"] for r in peptide_rows
              if r["autologous_robustness"] == "robust"}

    # -- event level --------------------------------------------------------
    event_rows = []
    for event in events:
        peps = event_peptides[event["sample_sv_id"]]
        row = dict(event)
        row["n_peptides_binding_panel"] = sum(p in panel_binders for p in peps)
        row["n_peptides_binding_autologous"] = sum(p in auto_binders for p in peps)
        row["n_peptides_unevaluable"] = sum(p in unevaluable for p in peps)
        row["n_peptides_binding_autologous_robust"] = sum(p in robust for p in peps)
        row["binds_panel"] = row["n_peptides_binding_panel"] > 0
        row["binds_autologous"] = row["n_peptides_binding_autologous"] > 0
        row["binds_autologous_robust"] = row["n_peptides_binding_autologous_robust"] > 0
        event_rows.append(row)
    _write(out_dir / f"events_binding{label}.tsv", event_rows)

    # -- the funnel, three branches ----------------------------------------
    branches = {
        "all": funnel(event_rows, rna, lambda e: True),
        "mhc_panel": funnel(event_rows, rna, lambda e: e["binds_panel"]),
        "mhc_autologous": funnel(event_rows, rna, lambda e: e["binds_autologous"]),
        "mhc_autologous_robust": funnel(event_rows, rna,
                                        lambda e: e["binds_autologous_robust"]),
    }
    def scope_counts(subset: set[str]) -> dict[str, int]:
        return {"n_all": len(subset),
                "n_mhc_panel": len(subset & panel_binders),
                "n_mhc_autologous": len(subset & auto_binders),
                "n_mhc_autologous_robust": len(subset & robust),
                "n_unevaluable": len(subset & unevaluable)}

    peptide_counts = {"matches_unique_peptides": scope_counts(matched_peptides),
                      "peptides_in_credible_events": scope_counts(credible_peptides)}

    # One schema for every row, so the three branches read down the same columns.
    # The peptide row carries the branch counts directly; `n_unevaluable` is only
    # meaningful for the autologous branch, and is blank elsewhere rather than 0.
    funnel_rows = [{"step": step, "scope": "peptide", **counts}
                   for step, counts in peptide_counts.items()]
    for step in ("events", "events_private", "events_hc", "events_rna_supported",
                 "events_hc_and_rna", "events_private_hc_and_rna"):
        funnel_rows.append({"step": step, "scope": "event",
                            "n_all": branches["all"][step],
                            "n_mhc_panel": branches["mhc_panel"][step],
                            "n_mhc_autologous": branches["mhc_autologous"][step],
                            "n_mhc_autologous_robust":
                                branches["mhc_autologous_robust"][step],
                            "n_unevaluable": ""})
    _write(out_dir / f"funnel_mhc{label}.tsv", funnel_rows)

    with (out_dir / f"mhc_propagation{label}.json").open("w") as fh:
        json.dump({"sample": args.sample,
                   "criterion": args.criterion,
                   "binder_thresholds": cuts,
                   "predictions": str(predictions),
                   "panel_alleles": len(panel),
                   "peptides": peptide_counts,
                   "funnel": branches}, fh, indent=2)

    print()
    for step, counts in peptide_counts.items():
        print(f"  {step:<28} all={counts['n_all']:>6}"
              f"  panel={counts['n_mhc_panel']:>5}"
              f"  autologous={counts['n_mhc_autologous']:>5}"
              f"  robust={counts['n_mhc_autologous_robust']:>5}"
              f"  (unevaluable {counts['n_unevaluable']:,})")
    for step in ("events", "events_private", "events_hc", "events_rna_supported",
                 "events_hc_and_rna", "events_private_hc_and_rna"):
        print(f"  {step:<28} all={branches['all'][step]:>5}"
              f"  panel={branches['mhc_panel'][step]:>5}"
              f"  autologous={branches['mhc_autologous'][step]:>5}"
              f"  robust={branches['mhc_autologous_robust'][step]:>5}")
    print(f"\n  wrote {out_dir}/peptide_binding{label}.tsv, "
          f"events_binding{label}.tsv, funnel_mhc{label}.tsv, "
          f"mhc_propagation{label}.json")


def _write(path: pathlib.Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        # lineterminator is explicit: csv defaults to \r\n, which turns every
        # last column into "value\r" for awk, cut and grep. These files are meant
        # to be read with shell tools.
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
