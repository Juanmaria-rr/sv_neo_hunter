#!/usr/bin/env python3
"""
run_netmhcpan.py — MHC class I binding predictions for a peptide list.

WHAT THIS EMITS, AND WHY THAT SHAPE
-----------------------------------
One row per (peptide, allele), carrying netMHCpan's raw numbers: eluted-ligand
score and rank, binding-affinity score and rank, and affinity in nM.

It deliberately does NOT collapse that to a per-peptide verdict. Any such
collapse — "best rank across the panel", or a population coverage weighted by
allele frequencies — encodes a choice about what question is being asked, and
that choice is not settled here. Aggregation is
cheap arithmetic on this table; the prediction is not. Keeping them separate
means a change of convention costs seconds instead of a re-run.

`--summary` writes an additional per-peptide file, but only as a convenience
view over the same rows, and it says in its header which rule produced it.

TWO PROPERTIES OF netMHCpan THAT COST A DEBUGGING SESSION
---------------------------------------------------------
1. `TMPDIR`. macOS hands every process a per-user temporary directory whose
   path is ~49 characters (`/var/folders/3l/h7dyr.../T/`). netMHCpan crashes on
   it with `Trace/BPT trap` — no message, no exit code worth reading. A short
   path works. This module forces one rather than trusting the environment,
   because the failure mode looks like a corrupt install and is not.

2. Allele naming. `HLA-A02:01` is accepted; `HLA-A*02:01` silently yields zero
   rows rather than an error. Asterisks are stripped on the way in.

Mixed peptide lengths in one call are fine, and several alleles per call are
fine; both were verified against this install rather than assumed.

SCALE
-----
Cost is peptides x alleles. The full patient catalogue is ~10.4M peptides, which
is not a sensible target for a laptop; the peptides that survived a cross
(tens of thousands) are. Chunking keeps memory flat but does not make the
former cheap — choose the input set deliberately.

    python tools/run_netmhcpan.py \
        --peptides results_cross_patch002/parental_noPON/stage3_matches.tsv \
        --peptide-column peptide \
        --alleles HLA-A02:01,HLA-B07:02 \
        --out results_cross_patch002/parental_noPON/netmhcpan.tsv
"""
from __future__ import annotations

import argparse
import pathlib
import re
import shutil
import subprocess
import tempfile
import time

REPO = pathlib.Path(__file__).resolve().parent.parent

#: netMHCpan-4.x predicts class I over this length range. Anything else is
#: dropped with a count, never silently.
MIN_LENGTH, MAX_LENGTH = 8, 14

#: The 20 standard residues. Peptides carrying anything else (`X` from an
#: ambiguous codon, `*` from a stop) are not predictable and are dropped.
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

#: The binder definition used by the earlier analysis, recovered from
#: NeoSV-Trace's `mhc_filter` and its argument defaults, and confirmed against
#: its output: in 506,324 rows the maxima are exactly 500.0, 2.0 and 2.0 with no
#: violation, so the defaults were applied unmodified.
#:
#: All three must hold SIMULTANEOUSLY — `mhc_filter` ands them. That is stricter
#: than netMHCpan's own convention, where %Rank_EL <= 2 alone defines a weak
#: binder and <= 0.5 a strong one. No 0.5 threshold appears anywhere in that
#: analysis: there is one binder class, not two.
#:
#: The affinity cut overlaps with the BA rank cut without being implied by it —
#: they are different scales — so per-criterion attrition is worth reporting
#: separately to see which one actually bites.
BINDER_THRESHOLDS = {"affinity_nM": 500.0, "rank_BA": 2.0, "rank_EL": 2.0}

#: A data line from netMHCpan's plain-text output: rank, allele, peptide, ...
#: The header and the dashed rules are skipped by requiring a leading integer.
DATA_LINE = re.compile(r"^\s*\d+\s+\S+\s+[A-Z]+\s")

#: Columns taken from netMHCpan's output, by position, and the name each gets
#: here. Positions are stable across 4.1/4.2 for the `-p -BA` invocation.
FIELDS = [(1, "allele"), (2, "peptide"), (11, "score_EL"), (12, "rank_EL"),
          (13, "score_BA"), (14, "rank_BA"), (15, "affinity_nM")]

#: `-neo` (4.2 only) APPENDS two columns and moves none, so the positions above
#: stay valid; verified against this install by running the bundled test set with
#: and without the flag and confirming Score_EL / %Rank_EL / Score_BA / %Rank_BA
#: / Aff(nM) come back bit-identical. The extra head is the CEDAR neoepitope
#: fine-tuning described in the 4.2 paper (Nilsson et al., Front Immunol 2025).
NEO_FIELDS = [(16, "score_Neo"), (17, "rank_Neo")]


def find_netmhcpan(explicit: str | None) -> pathlib.Path:
    """Locate the launcher, or explain how to install it.

    netMHCpan is licensed and cannot live in this repository (see .gitignore),
    so it is found per-machine rather than vendored.
    """
    candidates = []
    if explicit:
        candidates.append(pathlib.Path(explicit))
    env = shutil.which("netMHCpan")
    if env:
        candidates.append(pathlib.Path(env))
    candidates += sorted(REPO.glob("netMHCpan-*/netMHCpan"))

    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()

    raise SystemExit(
        "netMHCpan not found.\n"
        "  It is licensed third-party software (DTU Health Tech, free for\n"
        "  academic use) and is deliberately not tracked in this repository.\n\n"
        "  Install it:\n"
        "    1. Request it at https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/\n"
        "    2. Unpack it inside this repository as netMHCpan-4.2/\n"
        "    3. Set NMHOME in its `netMHCpan` script to that absolute path\n"
        "    4. chmod +x netMHCpan-4.2/netMHCpan\n\n"
        "  Or pass --netmhcpan /path/to/netMHCpan.")


def read_peptides(path: pathlib.Path, column: str | None) -> list[str]:
    """Distinct peptides from a plain list or a named column of a TSV."""
    lines = path.read_text().splitlines()
    if not lines:
        raise SystemExit(f"{path} is empty")

    if column:
        header = lines[0].split("\t")
        if column not in header:
            raise SystemExit(
                f"column {column!r} not in {path.name}.\n"
                f"  available: {', '.join(header[:20])}")
        index = header.index(column)
        values = [line.split("\t")[index] for line in lines[1:]
                  if len(line.split("\t")) > index]
    else:
        # A peptide list may carry extra whitespace-separated fields — the
        # netMHCpan distribution's own test/test.pep is `peptide score`, and
        # `-p` reads only the first token. Match that, or every line looks
        # "too long" and the run aborts with nothing predictable.
        values = [line.split()[0] for line in lines if line.split()]

    # Order-preserving deduplication: the same peptide reached by different
    # junctions is one prediction, but stable order keeps runs diffable.
    seen, peptides = set(), []
    for value in values:
        peptide = value.strip().upper()
        if peptide and peptide not in seen:
            seen.add(peptide)
            peptides.append(peptide)
    return peptides


def partition(peptides: list[str]) -> tuple[list[str], dict[str, int]]:
    """Split into predictable peptides and a tally of why the rest were dropped.

    Reported rather than silently filtered: a large drop count usually means the
    wrong column was passed, and that should be visible.
    """
    keep, dropped = [], {"too_short": 0, "too_long": 0, "non_standard_aa": 0}
    for peptide in peptides:
        if len(peptide) < MIN_LENGTH:
            dropped["too_short"] += 1
        elif len(peptide) > MAX_LENGTH:
            dropped["too_long"] += 1
        elif not set(peptide) <= STANDARD_AA:
            dropped["non_standard_aa"] += 1
        else:
            keep.append(peptide)
    return keep, dropped


def normalise_alleles(raw: str) -> list[str]:
    """Accept a comma list or a file; strip the asterisk netMHCpan rejects.

    A panel file carries its own provenance in `#` comments, and those must be
    stripped before splitting — netMHCpan does not reject an unknown token like
    `#`, it aborts partway through with "cannot be found in hla_pseudo list".
    """
    path = pathlib.Path(raw)
    if path.is_file():
        text = "\n".join(line.split("#", 1)[0]
                         for line in path.read_text().splitlines())
    else:
        text = raw

    # Order-preserving deduplication: a repeated allele would otherwise be
    # predicted twice and double-count in any per-peptide allele tally.
    seen, alleles = set(), []
    for token in re.split(r"[,\s]+", text):
        allele = token.strip().replace("*", "")
        if allele and allele not in seen:
            seen.add(allele)
            alleles.append(allele)

    if not alleles:
        raise SystemExit(f"no alleles parsed from {raw!r}")
    return alleles


#: netMHCpan caps the `-a` argument, so a cohort panel of ~220 alleles (about
#: 2,400 characters) has to be batched. There is no allele-file option; `-h`
#: confirms `-a` is the only route.
#:
#: The cap it reports is not the cap it has. Exceeding 1024 characters gives a
#: clean "PTYPE_LINE parameter option -a too long. Max size 1024", but well
#: below that it dies with `Trace/BPT trap` instead. Measured on this install
#: with 20 peptides: 88 alleles / 970 chars succeeds, 90 / 992 traps. Whether
#: the true limit counts alleles or characters was not worth pinning down —
#: 800 characters (~66 alleles) sits clear of both, and the cost of a smaller
#: batch is only a few more subprocess calls.
MAX_ALLELE_ARG_CHARS = 800


def allele_batches(alleles: list[str]) -> list[list[str]]:
    """Split the panel into groups whose joined form fits netMHCpan's `-a`."""
    batches, current, length = [], [], 0
    for allele in alleles:
        extra = len(allele) + (1 if current else 0)
        if current and length + extra > MAX_ALLELE_ARG_CHARS:
            batches.append(current)
            current, length = [allele], len(allele)
        else:
            current.append(allele)
            length += extra
    if current:
        batches.append(current)
    return batches


def run_chunk(launcher: pathlib.Path, peptides: list[str], alleles: list[str],
              tmpdir: pathlib.Path, extra_flags: list[str] | None = None) -> list[list[str]]:
    """One netMHCpan call; returns parsed rows."""
    peptide_file = tmpdir / "chunk.pep"
    peptide_file.write_text("\n".join(peptides) + "\n")

    # TMPDIR is forced to a short path here — see the module docstring. Passing
    # it through the environment is the only lever: the launcher honours an
    # already-set TMPDIR and only falls back to /tmp when it is unset.
    environment = {"TMPDIR": str(tmpdir), "PATH": "/usr/bin:/bin:/usr/sbin:/sbin",
                   "HOME": str(pathlib.Path.home())}

    result = subprocess.run(
        [str(launcher), "-p", str(peptide_file), "-a", ",".join(alleles), "-BA",
         *(extra_flags or [])],
        capture_output=True, text=True, env=environment)

    rows = [line.split() for line in result.stdout.splitlines()
            if DATA_LINE.match(line)]

    if not rows:
        raise SystemExit(
            f"netMHCpan returned no predictions for {len(peptides)} peptide(s).\n"
            f"  exit code: {result.returncode}\n"
            f"  stderr: {result.stderr.strip()[:500] or '(empty)'}\n"
            f"  stdout tail: {result.stdout.strip()[-500:] or '(empty)'}\n"
            "  A `Trace/BPT trap` here means the temporary path is too long.")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="netMHCpan class I predictions, one row per peptide-allele.")
    parser.add_argument("--peptides", required=True,
                        help="plain peptide list, or a TSV with --peptide-column")
    parser.add_argument("--peptide-column", default=None)
    parser.add_argument("--alleles", required=True,
                        help="comma-separated, or a file of allele names. No "
                             "default: the panel is an analytical choice.")
    parser.add_argument("--out", required=True)
    parser.add_argument("--summary", default=None,
                        help="also write a per-peptide best-rank view")
    parser.add_argument("--neo", action="store_true",
                        help="netMHCpan 4.2 only: also emit the CEDAR "
                             "neoepitope-finetuned score and %%rank as "
                             "score_Neo / rank_Neo. Adds columns; the EL and BA "
                             "columns are unchanged, so an output produced with "
                             "this flag is a superset of one produced without.")
    parser.add_argument("--chunk-size", type=int, default=5000)
    parser.add_argument("--netmhcpan", default=None)
    args = parser.parse_args()

    launcher = find_netmhcpan(args.netmhcpan)
    alleles = normalise_alleles(args.alleles)
    requested = read_peptides(pathlib.Path(args.peptides), args.peptide_column)
    peptides, dropped = partition(requested)

    print(f"  netMHCpan: {launcher}")
    print(f"  alleles:   {len(alleles)} ({', '.join(alleles[:6])}"
          f"{', ...' if len(alleles) > 6 else ''})")
    print(f"  peptides:  {len(peptides):,} predictable of {len(requested):,}")
    for reason, count in dropped.items():
        if count:
            print(f"             dropped {count:,} ({reason.replace('_', ' ')})")
    if not peptides:
        raise SystemExit("nothing predictable — check --peptide-column")
    print(f"  predictions to make: {len(peptides) * len(alleles):,}")

    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fields = FIELDS + (NEO_FIELDS if args.neo else [])
    extra_flags = ["-neo"] if args.neo else []
    if args.neo:
        print("  -neo: adding score_Neo / rank_Neo (4.2 only)")

    batches = allele_batches(alleles)
    if len(batches) > 1:
        print(f"  allele batches: {len(batches)} "
              f"(netMHCpan caps -a at {MAX_ALLELE_ARG_CHARS} chars)")

    written, started = 0, time.time()
    with tempfile.TemporaryDirectory(prefix="nmp_", dir="/tmp") as raw_tmp:
        tmpdir = pathlib.Path(raw_tmp)
        with out_path.open("w") as out:
            out.write("\t".join(name for _, name in fields) + "\n")
            for start in range(0, len(peptides), args.chunk_size):
                chunk = peptides[start:start + args.chunk_size]
                for batch in batches:
                    for row in run_chunk(launcher, chunk, batch, tmpdir,
                                         extra_flags):
                        out.write("\t".join(row[i] for i, _ in fields) + "\n")
                        written += 1
                done = min(start + args.chunk_size, len(peptides))
                elapsed = time.time() - started
                eta = elapsed * (len(peptides) - done) / done if done else 0
                print(f"    {done:,}/{len(peptides):,} peptides -> "
                      f"{written:,} rows  [{elapsed / 60:.1f} min elapsed, "
                      f"~{eta / 60:.0f} min left]", flush=True)

    print(f"  wrote {written:,} rows -> {args.out}")

    if args.summary:
        write_summary(out_path, pathlib.Path(args.summary), len(alleles))


def write_summary(rows_path: pathlib.Path, summary_path: pathlib.Path,
                  n_alleles: int) -> None:
    """Per-peptide binder counts under the earlier analysis's own definition."""
    import pandas as pd

    table = pd.read_csv(rows_path, sep="\t")
    for column in BINDER_THRESHOLDS:
        table[column] = pd.to_numeric(table[column], errors="coerce")

    passes = {column: table[column] <= threshold
              for column, threshold in BINDER_THRESHOLDS.items()}
    binder = passes["affinity_nM"] & passes["rank_BA"] & passes["rank_EL"]
    table["is_binder"] = binder

    aggregations = dict(
        best_rank_EL=("rank_EL", "min"),
        best_rank_BA=("rank_BA", "min"),
        best_affinity_nM=("affinity_nM", "min"),
        n_alleles_binding=("is_binder", "sum"),
    )
    # The neoepitope head is carried through but NOT folded into `is_binder`:
    # there is no established rank cut for it, and inventing one here would
    # quietly change what "binder" means between runs. It is reported alongside
    # so a criterion can be chosen downstream, deliberately.
    if "rank_Neo" in table.columns:
        table["rank_Neo"] = pd.to_numeric(table["rank_Neo"], errors="coerce")
        aggregations["best_rank_Neo"] = ("rank_Neo", "min")
    summary = table.groupby("peptide").agg(**aggregations).reset_index()
    summary["n_alleles_binding"] = summary["n_alleles_binding"].astype(int)
    summary["fraction_alleles_binding"] = \
        summary["n_alleles_binding"] / n_alleles

    summary.to_csv(summary_path, sep="\t", index=False)

    print(f"  wrote {len(summary):,} peptides -> {summary_path}")
    print(f"    binder = affinity <= {BINDER_THRESHOLDS['affinity_nM']:g} nM"
          f" AND %Rank_BA <= {BINDER_THRESHOLDS['rank_BA']:g}"
          f" AND %Rank_EL <= {BINDER_THRESHOLDS['rank_EL']:g}, all three")

    # Per-criterion attrition. The three cuts overlap without implying one
    # another, so the marginal pass rates say which one is actually doing the
    # filtering — and the conjunction is necessarily below every margin.
    total = len(table)
    print(f"    peptide-allele rows: {total:,}")
    for column, threshold in BINDER_THRESHOLDS.items():
        share = passes[column].mean()
        print(f"      {column:12} <= {threshold:6g}: "
              f"{int(passes[column].sum()):>10,}  ({100 * share:5.2f}%)")
    print(f"      {'all three':12}          : {int(binder.sum()):>10,}  "
          f"({100 * binder.mean():5.2f}%)")

    with_binder = int((summary["n_alleles_binding"] > 0).sum())
    print(f"    peptides with >=1 binding allele: {with_binder:,}"
          f" of {len(summary):,} ({100 * with_binder / len(summary):.1f}%)")
    print("    NOTE: a panel statistic, not population coverage — unweighted "
          "by allele frequency")


if __name__ == "__main__":
    main()
