#!/usr/bin/env python
"""
query_sv_pon.py — look up any junction in the SV panel of normals the caller used.

WHY
---
`pon_count` reaches our tables only for the junctions the caller happened to
annotate, and only as a bare number: a junction with no `PON_COUNT` is
indistinguishable from one the caller never checked. The panel itself is a
table, so it can be asked directly — for a junction in our samples, for one that
produced no peptide, or for any coordinates at all.

THE RESOURCE
------------
hmftools' SV panel of normals, the same files the SV caller is given by
nf-core/oncoanalyser (`hmf_pipeline_resources.38_v2.0.0--3.tar.gz`):

    dna/sv/sv_pon.38.bedpe.gz    junctions: chrom1 start1 end1 chrom2 start2 end2 . count strand1 strand2
    dna/sv/sgl_pon.38.bed.gz     single breakends: chrom start end . count strand

Both store an INTERVAL per side, not a point, because callers disagree about a
breakpoint by a few bases; a query matches when each side falls inside the
record's interval (widened by --slop).

WHAT THE COUNT IS, AND IS NOT
-----------------------------
The number of panel samples carrying the junction. It is NOT a frequency: the
panel size is not stated in the resource or its documentation. What the file
does give is a lower bound — the largest count it contains — which this tool
prints with `--describe`, so a count can at least be read against it.

Usage
-----
    python tools/query_sv_pon.py --bedpe <sv_pon.38.bedpe.gz> --describe
    python tools/query_sv_pon.py --bedpe ... --junction chr1:3820477-chr1:3820524
    python tools/query_sv_pon.py --bedpe ... --table junctions.tsv --out annotated.tsv
"""
from __future__ import annotations

import argparse
import gzip
import pathlib
import sys


def parse_junction(text: str) -> tuple[str, int, str, int]:
    """'chr1:3820477-chr1:3820524' -> ('chr1', 3820477, 'chr1', 3820524)."""
    try:
        left, right = text.split("-chr") if "-chr" in text else text.split("-")
        if not right.startswith("chr"):
            right = "chr" + right
        c1, p1 = left.rsplit(":", 1)
        c2, p2 = right.rsplit(":", 1)
        return c1, int(p1), c2, int(p2)
    except ValueError:
        raise SystemExit(f"  cannot parse junction '{text}'; "
                         f"expected chrom:pos-chrom:pos")


def norm(chrom: str) -> str:
    return str(chrom).replace("chr", "")


def describe(path: pathlib.Path) -> None:
    biggest = total = 0
    with gzip.open(path, "rt") as handle:
        for line in handle:
            f = line.split("\t")
            if len(f) < 8:
                continue
            total += 1
            count = int(f[7])
            biggest = max(biggest, count)
    print(f"{path.name}: {total:,} records, largest count {biggest:,}")
    print(f"  -> the panel holds AT LEAST {biggest:,} samples; the resource does not "
          f"state its size, so this is a lower bound and a count is only an ordinal.")


def lookup(path: pathlib.Path, queries: list[tuple[str, int, str, int]],
           slop: int) -> list[list[dict]]:
    """One pass over the panel for every query; all matches per query."""
    wanted = {}
    for i, (c1, p1, c2, p2) in enumerate(queries):
        wanted.setdefault(tuple(sorted((norm(c1), norm(c2)))), []).append(i)
    hits: list[list[dict]] = [[] for _ in queries]
    with gzip.open(path, "rt") as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            a, b = norm(f[0]), norm(f[3])
            for i in wanted.get(tuple(sorted((a, b))), ()):
                c1, p1, c2, p2 = queries[i]
                s1, e1, s2, e2 = int(f[1]), int(f[2]), int(f[4]), int(f[5])
                forward = (norm(c1) == a and s1 - slop <= p1 <= e1 + slop
                           and s2 - slop <= p2 <= e2 + slop)
                reverse = (norm(c1) == b and s2 - slop <= p1 <= e2 + slop
                           and s1 - slop <= p2 <= e1 + slop)
                if forward or reverse:
                    hits[i].append({"count": int(f[7]),
                                    "record": f"{f[0]}:{s1}-{e1} {f[3]}:{s2}-{e2}",
                                    "strands": f"{f[8]}{f[9]}" if len(f) > 9 else ""})
    return hits


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--bedpe", required=True, type=pathlib.Path)
    p.add_argument("--junction", action="append", default=[],
                   metavar="chrom:pos-chrom:pos", help="repeatable")
    p.add_argument("--table", type=pathlib.Path,
                   help="TSV with chrom1, pos1, chrom2, pos2 columns")
    p.add_argument("--out", type=pathlib.Path, help="annotated copy of --table")
    p.add_argument("--slop", type=int, default=0,
                   help="widen each side's interval before matching (default 0; "
                        "the panel's intervals already allow for imprecision)")
    p.add_argument("--describe", action="store_true")
    args = p.parse_args()

    if args.describe:
        describe(args.bedpe)
        if not (args.junction or args.table):
            return

    queries, frame = [q and parse_junction(q) for q in args.junction], None
    if args.table:
        import pandas as pd
        frame = pd.read_csv(args.table, sep="\t", low_memory=False)
        missing = [c for c in ("chrom1", "pos1", "chrom2", "pos2") if c not in frame.columns]
        if missing:
            raise SystemExit(f"  --table lacks {missing}")
        queries += [(str(r.chrom1), int(r.pos1), str(r.chrom2), int(r.pos2))
                    for r in frame.itertuples()]
    if not queries:
        raise SystemExit("  nothing to look up: pass --junction or --table")

    print(f"  querying {len(queries):,} junction(s) against {args.bedpe.name}",
          file=sys.stderr)
    hits = lookup(args.bedpe, queries, args.slop)
    counts = [max((h["count"] for h in hs), default=0) for hs in hits]

    for (c1, p1, c2, p2), hs, count in zip(queries, hits, counts):
        if hs:
            best = max(hs, key=lambda h: h["count"])
            print(f"{c1}:{p1}-{c2}:{p2}\tpon_count={count}\t{best['record']}\t"
                  f"{best['strands']}" + (f"\t({len(hs)} records)" if len(hs) > 1 else ""))
        else:
            print(f"{c1}:{p1}-{c2}:{p2}\tnot in the panel")

    if frame is not None and args.out:
        offset = len(args.junction)
        frame["pon_count_panel"] = counts[offset:]
        frame["in_panel"] = [bool(h) for h in hits[offset:]]
        frame.to_csv(args.out, sep="\t", index=False)
        print(f"  wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
