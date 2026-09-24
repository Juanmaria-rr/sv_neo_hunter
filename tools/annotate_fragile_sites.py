#!/usr/bin/env python
"""
annotate_fragile_sites.py — do the breakpoints fall in a common fragile site?

TWO DEFINITIONS, BECAUSE THERE IS NO SINGLE ONE
-----------------------------------------------
Common fragile sites are not points. They are regions of megabases whose limits
depend on the study, the cell type and the assay, and published catalogues
disagree by two orders of magnitude about how much genome they cover. So this
annotates against BOTH a narrow and a broad definition and records what each
says, rather than picking one and hiding the choice.

**The coverage of each definition is the fact that makes its column readable.**
A breakpoint landing in a definition that covers most of the genome is not
evidence of anything: it is the expected outcome. This script therefore measures
each catalogue's genome coverage and writes it into the provenance, so a hit
rate can be compared against the rate chance alone would produce. Reporting
"N% of our candidates are in fragile sites" without that denominator is the
mistake this design exists to prevent.

PER BREAKEND, NOT PER JUNCTION
------------------------------
A junction has two ends and they can fall on different sides of a boundary. Each
end is annotated separately and a status column says which: `none`, `bp1`,
`bp2`, or `both`. Collapsing that to one boolean would throw away the case that
matters most — one end inside a fragile site and the other far outside is a
different claim from both ends inside.

Nothing is filtered. These are columns.

Usage
-----
    python tools/annotate_fragile_sites.py \\
        --table <dir>/candidate_universe.tsv \\
        --narrow <cfs>/cfs_canonical_grch38.bed \\
        --broad  <cfs>/cfs_grch38.bed \\
        --out    <dir>/candidate_universe.tsv
"""
from __future__ import annotations

import argparse
import pathlib
import sys
from datetime import date

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from provenance import code_version                        # noqa: E402

#: Assumed genome size for the coverage figure. Only used to express a
#: catalogue's footprint as a percentage, never in the annotation itself.
GENOME_BP = 3.1e9


def normalise_chrom(value) -> str:
    """`chr1`, `1`, `x` -> `1`, `X`. Catalogues and VCFs disagree on all three."""
    text = str(value).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    return text.upper()


def load_bed(path: pathlib.Path) -> tuple[pd.DataFrame, float]:
    """A BED-ish catalogue plus the fraction of genome it covers.

    Accepts a header or none, and takes the name from column 4 when present.
    Extra columns (`gene`, `band`, `tier`) are kept when the file has a header.
    """
    first = path.read_text().split("\n", 1)[0].split("\t")
    has_header = not first[1].strip().isdigit()
    frame = pd.read_csv(path, sep="\t", header=0 if has_header else None,
                        dtype=str, comment="#")
    if not has_header:
        names = ["chrom", "start", "end", "name", "score", "strand"]
        frame.columns = names[:len(frame.columns)]
    frame = frame.rename(columns={frame.columns[0]: "chrom",
                                  frame.columns[1]: "start",
                                  frame.columns[2]: "end"})
    frame["chrom"] = frame["chrom"].map(normalise_chrom)
    frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
    frame["end"] = pd.to_numeric(frame["end"], errors="coerce")
    frame = frame.dropna(subset=["chrom", "start", "end"])
    if "name" not in frame.columns:
        frame["name"] = frame["chrom"] + ":" + frame["start"].astype(int).astype(str)

    # Merge overlaps before measuring coverage, or a region counted twice
    # inflates the footprint and with it the null rate.
    total = 0
    for _, group in frame.groupby("chrom"):
        spans = sorted(zip(group["start"], group["end"]))
        lo, hi = spans[0]
        for start, end in spans[1:]:
            if start <= hi:
                hi = max(hi, end)
            else:
                total += hi - lo
                lo, hi = start, end
        total += hi - lo
    return frame, total / GENOME_BP


def lookup(frame: pd.DataFrame, chrom: str, pos) -> dict | None:
    """The catalogue row containing this position, or None."""
    try:
        position = float(pos)
    except (TypeError, ValueError):
        return None
    hit = frame[(frame["chrom"] == chrom)
                & (frame["start"] <= position) & (position < frame["end"])]
    if hit.empty:
        return None
    return hit.iloc[0].to_dict()


def annotate(table: pd.DataFrame, catalogues: dict) -> pd.DataFrame:
    # One lookup per distinct breakend rather than per peptide: a locus appears
    # in thousands of rows and the answer cannot differ between them.
    ends = pd.concat([
        table[["chrom1", "pos1"]].rename(columns={"chrom1": "chrom", "pos1": "pos"}),
        table[["chrom2", "pos2"]].rename(columns={"chrom2": "chrom", "pos2": "pos"}),
    ]).drop_duplicates()
    ends["chrom"] = ends["chrom"].map(normalise_chrom)
    print(f"  distinct breakends to look up: {len(ends):,}")

    for label, (frame, coverage) in catalogues.items():
        cache = {}
        for chrom, pos in zip(ends["chrom"], ends["pos"]):
            cache[(chrom, str(pos))] = lookup(frame, chrom, pos)

        for side in ("1", "2"):
            keys = list(zip(table[f"chrom{side}"].map(normalise_chrom),
                            table[f"pos{side}"].astype(str)))
            hits = [cache.get(key) for key in keys]
            table[f"cfs_{label}_bp{side}"] = [h["name"] if h else "" for h in hits]
            if any(h and "tier" in h for h in hits if h):
                table[f"cfs_{label}_tier_bp{side}"] = [
                    h.get("tier", "") if h else "" for h in hits]

        bp1 = table[f"cfs_{label}_bp1"].astype(bool)
        bp2 = table[f"cfs_{label}_bp2"].astype(bool)
        table[f"cfs_{label}_status"] = [
            "both" if a and b else "bp1" if a else "bp2" if b else "none"
            for a, b in zip(bp1, bp2)]
        rate = 100.0 * (table[f"cfs_{label}_status"] != "none").mean()
        print(f"  {label:<8} covers {100 * coverage:5.1f}% of the genome; "
              f"{rate:5.1f}% of rows have at least one end inside "
              f"(chance alone would give about "
              f"{100 * (1 - (1 - coverage) ** 2):.0f}%)")

    # Where the two definitions disagree is the informative column: a hit under
    # the broad catalogue alone is weak by construction.
    if len(catalogues) == 2:
        a, b = list(catalogues)
        table["cfs_agreement"] = [
            "both" if x != "none" and y != "none"
            else f"{a}_only" if x != "none"
            else f"{b}_only" if y != "none" else "neither"
            for x, y in zip(table[f"cfs_{a}_status"], table[f"cfs_{b}_status"])]
    return table


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", type=pathlib.Path, required=True,
                        help="any table with chrom1/pos1/chrom2/pos2")
    parser.add_argument("--narrow", type=pathlib.Path, required=True,
                        help="the conservative catalogue")
    parser.add_argument("--broad", type=pathlib.Path, required=True,
                        help="the permissive catalogue")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()

    catalogues = {}
    for label, path in (("narrow", args.narrow), ("broad", args.broad)):
        frame, coverage = load_bed(path)
        catalogues[label] = (frame, coverage)
        print(f"  {label}: {len(frame)} regions from {path.name}")

    table = pd.read_csv(args.table, sep="\t", dtype=str, low_memory=False)
    print(f"  table: {len(table):,} rows")
    table = annotate(table, catalogues)
    table.to_csv(args.out, sep="\t", index=False)

    provenance = [*code_version(),
                  {"key": "generated", "value": date.today().isoformat()},
                  {"key": "tool", "value": pathlib.Path(__file__).name}]
    for label, (frame, coverage) in catalogues.items():
        source = args.narrow if label == "narrow" else args.broad
        provenance += [
            {"key": f"{label}_catalogue", "value": str(source)},
            {"key": f"{label}_regions", "value": len(frame)},
            {"key": f"{label}_genome_coverage_pct", "value": round(100 * coverage, 2)},
            {"key": f"{label}_expected_hit_rate_pct",
             "value": round(100 * (1 - (1 - coverage) ** 2), 1)},
        ]
    provenance.append({"key": "note", "value":
                       "Expected hit rate is the rate two independent random "
                       "breakends would produce. A observed rate at or below it "
                       "is NOT evidence of fragile-site involvement."})
    pd.DataFrame(provenance).to_csv(
        args.out.parent / "CFS_ANNOTATION_PROVENANCE.tsv", sep="\t", index=False)

    print(f"\n  wrote {args.out}")
    for label in catalogues:
        print(f"\n  cfs_{label}_status:")
        print(table[f"cfs_{label}_status"].value_counts().to_string())
    if "cfs_agreement" in table.columns:
        print("\n  cfs_agreement:")
        print(table["cfs_agreement"].value_counts().to_string())


if __name__ == "__main__":
    main()
