#!/usr/bin/env python
"""
snapshot_candidate_inputs.py — give the candidate analysis its own inputs.

THE PROBLEM THIS SOLVES, AND THE ONE IT COULD CREATE
----------------------------------------------------
Candidate neopeptides are stage-2 output: generated from a sample's own SV
junctions, long before any reference cohort is consulted. They happen to be
written into a run directory whose purpose was the cross, so a candidate
analysis that reads them from there silently depends on that directory. Re-run
the cross, or edit anything under it, and every downstream candidate table
changes with nothing to say so.

Copying the two files the candidate analysis actually needs breaks that
coupling. But a copy creates the mirror problem: it goes stale when the source
is legitimately regenerated, and again nothing says so. A snapshot without drift
detection is not an improvement, it is the same failure pointing the other way.

So every copied file is recorded with its origin path and checksum, and
`--check` re-reads the origins and reports three states per file:

    current   the origin still hashes to what was copied
    DRIFTED   the origin changed — the snapshot is stale, re-snapshot decide
    MISSING   the origin is gone — the snapshot is now the only copy

WHAT IS COPIED, AND WHAT DELIBERATELY IS NOT
--------------------------------------------
Only the two files the candidate universe reads:

    *.all_neopeptides.txt        stage 2 — the candidates themselves
    stage6_8_all_junctions.tsv   stages 6-8 — confidence, privacy, expression
                                 and RNA evidence for EVERY admitted junction,
                                 not just the ones that matched

Everything else in a run directory is either cross output (which the candidate
analysis must not consult) or too large to duplicate for no gain. The snapshot
is a fraction of a percent of the run it comes from.

Usage
-----
    python tools/snapshot_candidate_inputs.py \\
        --from <cross>/parental_noPON --from <cross>/clone_A_noPON \\
        --out-dir <candidates>/inputs

    python tools/snapshot_candidate_inputs.py --out-dir <candidates>/inputs --check
"""
from __future__ import annotations

import argparse
import hashlib
import pathlib
import shutil
import sys
from datetime import datetime

import pandas as pd

#: Filenames the candidate analysis reads. A glob so the stage-2 file can carry
#: whatever prefix its run used.
WANTED = ["*.all_neopeptides.txt", "stage6_8_all_junctions.tsv"]

SOURCE_FILE = "SOURCE.tsv"


def checksum(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def snapshot(run_dir: pathlib.Path, out_dir: pathlib.Path) -> list[dict]:
    target = out_dir / run_dir.name
    target.mkdir(parents=True, exist_ok=True)
    rows = []
    for pattern in WANTED:
        hits = sorted(run_dir.glob(pattern))
        if not hits:
            print(f"    WARNING: {run_dir.name} has no {pattern} — "
                  "downstream columns from it will be absent")
            continue
        for source in hits:
            digest = checksum(source)
            shutil.copy2(source, target / source.name)
            rows.append({"file": source.name,
                         "origin": str(source.resolve()),
                         "origin_sha256": digest,
                         "bytes": source.stat().st_size,
                         "copied": datetime.now().isoformat(timespec="seconds")})
            print(f"    {source.name:<44} {source.stat().st_size / 1e6:>7.1f} MB")
    pd.DataFrame(rows).to_csv(target / SOURCE_FILE, sep="\t", index=False)
    return rows


def check(out_dir: pathlib.Path) -> int:
    """Compare every recorded origin against its current state. Returns drift."""
    drifted = 0
    records = sorted(out_dir.glob(f"*/{SOURCE_FILE}"))
    if not records:
        raise SystemExit(f"no {SOURCE_FILE} under {out_dir} — nothing to check")
    for record in records:
        sample = record.parent.name
        table = pd.read_csv(record, sep="\t")
        for row in table.itertuples(index=False):
            origin = pathlib.Path(row.origin)
            if not origin.exists():
                state = "MISSING"
            elif checksum(origin) == row.origin_sha256:
                state = "current"
            else:
                state = "DRIFTED"
            if state != "current":
                drifted += 1
            print(f"  {state:<8} {sample}/{row.file}")
            if state == "DRIFTED":
                print(f"           origin changed since {row.copied}: {origin}")
            elif state == "MISSING":
                print(f"           origin gone: {origin}")
    return drifted


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--from", dest="runs", action="append", default=[],
                        type=pathlib.Path, metavar="RUN_DIR",
                        help="a pipeline run directory to take inputs from; "
                             "repeatable")
    parser.add_argument("--out-dir", type=pathlib.Path, required=True)
    parser.add_argument("--check", action="store_true",
                        help="do not copy; report whether each recorded origin "
                             "still matches what was copied. Exits 1 on drift, "
                             "so it can gate a rebuild.")
    args = parser.parse_args()

    if args.check:
        drifted = check(args.out_dir)
        print(f"\n  {drifted} file(s) drifted or missing")
        if drifted:
            print("  the snapshot no longer reflects its sources — re-snapshot, "
                  "or decide deliberately to keep the old inputs")
        sys.exit(1 if drifted else 0)

    if not args.runs:
        raise SystemExit("--from is required unless --check")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    total = 0
    for run_dir in args.runs:
        if not run_dir.exists():
            raise SystemExit(f"missing run directory: {run_dir}")
        print(f"  {run_dir.name}")
        total += len(snapshot(run_dir, args.out_dir))

    print(f"\n  snapshotted {total} file(s) into {args.out_dir}")
    print("  each sample carries a SOURCE.tsv with its origin and checksum.")
    print("  run with --check before trusting a rebuild: it tells you whether "
          "the sources moved under you.")


if __name__ == "__main__":
    main()
