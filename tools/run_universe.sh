#!/usr/bin/env bash
# run_universe.sh — build the candidate universe end to end, in order.
#
# WHY THIS FILE EXISTS
# --------------------
# The annotators that build the universe are order-dependent and not
# idempotent: several refuse to overwrite a column they would produce, and one
# retires columns an earlier step created. Run them in a different order and
# you get a different table, or an error. That order lived only in whoever ran
# them last. It lives here now.
#
# SAFETY
# ------
# --out-dir is REQUIRED and there is no default, so no invocation can land on a
# working directory by accident. The script refuses to write into a directory
# that already holds a candidate_universe.tsv unless --force is given, which is
# the case that would destroy work in progress. To verify a change, build into
# a NEW directory and diff against the old one; never rebuild in place.
#
# Usage
# -----
#   tools/run_universe.sh --config config/<run>.yaml --out-dir <new dir>
#   tools/run_universe.sh --config ... --out-dir ... --dry-run
set -euo pipefail

CONFIG=""; OUT=""; FORCE=0; DRY=0
CROSS=""   # directory holding the per-sample pipeline runs
MASTER=""  # master_peptides.tsv, for the matched_reference annotation

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)  CONFIG="$2"; shift 2 ;;
    --out-dir) OUT="$2";    shift 2 ;;
    --cross)   CROSS="$2";  shift 2 ;;
    --master)  MASTER="$2"; shift 2 ;;
    --force)   FORCE=1;     shift ;;
    --dry-run) DRY=1;       shift ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

[[ -n "$CONFIG" ]] || { echo "  --config is required" >&2; exit 2; }
[[ -n "$OUT"    ]] || { echo "  --out-dir is required (deliberately no default)" >&2; exit 2; }
[[ -f "$CONFIG" ]] || { echo "  config not found: $CONFIG" >&2; exit 2; }

if [[ -e "$OUT/rpe1_universe/candidate_universe.tsv" && $FORCE -eq 0 ]]; then
  cat >&2 <<MSG
  REFUSING TO WRITE: $OUT already holds a candidate_universe.tsv.

  Rebuilding in place would overwrite work in progress, and if a step failed
  half way you would be left with a table that is neither the old one nor a
  complete new one. Build into a new directory and compare:

      tools/run_universe.sh --config $CONFIG --out-dir ${OUT}_rebuild
      diff <(head -1 $OUT/rpe1_universe/candidate_universe.tsv | tr '\t' '\n') \\
           <(head -1 ${OUT}_rebuild/rpe1_universe/candidate_universe.tsv | tr '\t' '\n')

  --force overrides this. It is not recoverable.
MSG
  exit 1
fi

U="$OUT/rpe1_universe"
run() {
  echo; echo "── $1"; shift
  if [[ $DRY -eq 1 ]]; then printf '   %q' "$@"; echo; else "$@"; fi
}

echo "config:   $CONFIG"
echo "out-dir:  $OUT"
echo "code:     $(git rev-parse --short HEAD 2>/dev/null || echo 'not a git checkout')$(git diff --quiet 2>/dev/null || echo ' (DIRTY)')"
[[ $DRY -eq 1 ]] && echo "MODE:     dry run, nothing will be written"

# The lineage and the per-sample run directories come from the config; they are
# passed explicitly rather than discovered, so a rerun cannot silently pick up
# a different set of inputs than the one recorded in this invocation.
: "${CROSS:?  --cross is required: the directory holding the per-sample pipeline runs}"
: "${MASTER:?  --master is required: master_peptides.tsv}"

# The sample names and their parent/child relationships come from the config,
# which is not in this repository: those names identify the material. Deriving
# them keeps this script free of them, so it can be tracked publicly and still
# describe the real run.
LINES=$(python tools/lineage_args.py --config "$CONFIG" --cross "$CROSS" --what lines)
SAMPLES=$(python tools/lineage_args.py --config "$CONFIG" --cross "$CROSS" --what samples)
FROMS=$(python tools/lineage_args.py --config "$CONFIG" --cross "$CROSS" --what from)

run "1/11  candidate universe, one row per peptide across the lineage" \
    bash -c "python tools/build_candidate_universe.py $LINES \
      --master '$MASTER' $SAMPLES --out-dir '$U'"

run "2/11  checksummed snapshot of the stage-6/8 inputs" \
    bash -c "python tools/snapshot_candidate_inputs.py $FROMS --out-dir '$OUT/inputs'"

run "3/11  common fragile sites, two catalogues, per breakend" \
    python tools/annotate_fragile_sites.py --table "$U/candidate_universe.tsv" \
      --narrow "${CFS_NARROW:?set CFS_NARROW to the conservative catalogue BED}" \
      --broad  "${CFS_BROAD:?set CFS_BROAD to the permissive catalogue BED}" \
      --out "$U/candidate_universe.tsv"

run "4/11  IGV loci" \
    python tools/annotate_igv_loci.py --table "$U/candidate_universe.tsv" \
      --out "$U/candidate_universe.tsv" --chr-prefix chr

run "5/11  junction usage from the RNA BAMs" \
    python tools/annotate_junction_usage.py --table "$U/candidate_universe.tsv" \
      --config "$CONFIG" --out "$U/candidate_universe.tsv"

run "6/11  gnomAD per-ancestry frequencies" \
    python tools/annotate_gnomad_populations.py --table "$U/candidate_universe.tsv" \
      --inputs "$OUT/inputs" --out "$U/candidate_universe.tsv"

run "7/11  DNA variant allele fraction" \
    python tools/annotate_dna_vaf.py --table "$U/candidate_universe.tsv" \
      --inputs "$OUT/inputs" --out "$U/candidate_universe.tsv"

# Must come after 6: it reads gnomad_af_popmax, and it retires columns.
run "8/11  population-frequency verdict, and retire the superseded columns" \
    python tools/annotate_population_frequency.py --table "$U/candidate_universe.tsv" \
      --out "$U/candidate_universe.tsv"

# Must come last of the table steps: it documents whatever columns exist.
run "9/11  column dictionary and glossary" \
    python tools/build_universe_dictionary.py --universe "$U/candidate_universe.tsv" \
      --out "$U/candidate_universe_column_dictionary.tsv"

run "10/11 views, all lines and per line" \
    python tools/build_universe_views.py --universe "$U/candidate_universe.tsv" \
      --out-dir "$U/views" --view all --per-line

# Optional: a subset restricted to some of the lines. --subset asks for it as
# "name=line,line", so the names stay in the caller's hands and out of here.
if [[ -n "${SUBSET:-}" ]]; then
  NAME="${SUBSET%%=*}"; MEMBERS="${SUBSET#*=}"
  ARGS=""; for m in ${MEMBERS//,/ }; do ARGS="$ARGS --lines $m"; done
  run "11/11 subset: $NAME" \
      bash -c "python tools/build_universe_views.py \
        --universe '$U/candidate_universe.tsv' --out-dir '$U/subsets' \
        --view all --per-line $ARGS --subset-name '$NAME'"
else
  echo; echo "── 11/11 no subset requested (set SUBSET='name=lineA,lineB')"
fi

if [[ $DRY -eq 0 ]]; then
  echo; echo "── verifying every derived file matches the master it claims"
  python tools/check_derived_freshness.py --master "$U/candidate_universe.tsv"
fi
echo; echo "done: $U"
