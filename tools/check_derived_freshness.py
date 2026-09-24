#!/usr/bin/env python
"""
check_derived_freshness.py — is every derived file current with its source?

THE FAILURE THIS CATCHES
------------------------
Annotations are applied to the master table one at a time, each adding columns.
Views and subsets are built FROM that master and are not rebuilt when it
changes. Nothing errors when they fall behind: the stale files are still there,
still readable, still named as though they were current.

That is worse than a missing file. A subset carrying 108 columns beside a master
carrying 125 looks authoritative, and the columns it lacks are exactly the ones
most recently added — the ones someone is most likely to have come looking for.
Observed here: views and a subset silently missing two whole annotation passes.

The detection mechanism already existed. Every generated file records the
sha256 of the table it was built from, precisely so current can be told from
stale. What was missing was anything that actually ran the comparison, which
left it depending on someone remembering. This is that check.

WHAT IT DOES NOT DO
-------------------
It does not rebuild anything, because it cannot know the arguments each file was
built with — a view may carry a line restriction that is not recoverable from
its output. It reports, and exits non-zero so a pipeline step can refuse to
continue on stale inputs.

A file whose provenance records no source checksum is reported separately rather
than passed over. Not checkable is not the same as current, and folding the two
together would be the same mistake this script exists to prevent.

Usage
-----
    python tools/check_derived_freshness.py --master <dir>/candidate_universe.tsv
    python tools/check_derived_freshness.py --master <dir>/table.tsv --root <dir>
"""
from __future__ import annotations

import argparse
import hashlib
import pathlib
import sys

#: The provenance key every generator writes. Kept here rather than inlined so a
#: change to the convention is a one-line change in one place.
SOURCE_KEY = "source_sha256_16"
DIGEST_CHARS = 16


def digest(path: pathlib.Path) -> str:
    """sha256 of a file, truncated as the provenance files record it."""
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(chunk)
    return hasher.hexdigest()[:DIGEST_CHARS]


def provenance(path: pathlib.Path) -> dict[str, str]:
    rows = {}
    for line in path.read_text().splitlines()[1:]:
        if "\t" in line:
            key, value = line.split("\t", 1)
            rows[key] = value
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--master", type=pathlib.Path, required=True,
                        help="the table whose checksum derived files must match")
    parser.add_argument("--root", type=pathlib.Path, default=None,
                        help="directory to scan, recursively "
                             "(default: the master's own directory)")
    parser.add_argument("--quiet", action="store_true",
                        help="print only what is wrong")
    args = parser.parse_args()

    if not args.master.exists():
        raise SystemExit(f"  master not found: {args.master}")
    root = args.root or args.master.parent

    current = digest(args.master)
    columns = len(args.master.open().readline().rstrip("\n").split("\t"))
    print(f"  master {args.master.name}: {current}, {columns} columns")

    # Two naming conventions are in use: `<output>.PROVENANCE.tsv` beside a file
    # derived from the master, and `<STEP>_PROVENANCE.tsv` written by an
    # annotator that modified the master in place. Both are scanned. The second
    # kind records no source checksum and so lands in `unchecked` — which is the
    # point: passing over it silently is the failure this script exists to stop.
    found = {path for pattern in ("*.PROVENANCE.tsv", "*_PROVENANCE.tsv")
             for path in root.rglob(pattern)}

    fresh, stale, unchecked = [], [], []
    for path in sorted(found):
        claimed = provenance(path).get(SOURCE_KEY)
        relative = path.relative_to(root)
        if claimed is None:
            unchecked.append(relative)
        elif claimed == current:
            fresh.append(relative)
        else:
            stale.append((relative, claimed))

    if not args.quiet:
        print(f"  {len(fresh)} current")
    if unchecked:
        print(f"  {len(unchecked)} record no {SOURCE_KEY} — NOT checkable, which is "
              f"not the same as current:")
        for relative in unchecked:
            print(f"      {relative}")
    if stale:
        print(f"\n  {len(stale)} STALE — built from a different version of the master:")
        for relative, claimed in stale:
            print(f"      {relative}")
            print(f"        claims {claimed}, master is {current}")
        print("\n  Rebuild these before reading them. Their column set is behind the "
              "master's, so a column added since will be missing rather than empty.")
        sys.exit(1)

    print("  nothing stale")


if __name__ == "__main__":
    main()
