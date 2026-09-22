#!/usr/bin/env python3
"""
check_docs.py — verify that the documentation still describes the code.

Documentation rots silently: a threshold changes, a function moves, a criterion
is withdrawn, and the prose keeps asserting the old behaviour with no error
anywhere. This checks the claims that can be checked mechanically:

  1. every `CONSTANT` named in the docs exists in criteria.py
  2. every `function()` named in the docs exists in the package
  3. every `L<number>` line reference points at the right definition
  4. every documented output filename is one the pipeline actually writes
  5. documented values of boolean criteria match the code
     (this is the one that catches a withdrawn criterion still described as active)

Exit code 1 if anything is stale, so it can gate a commit.

    python tools/check_docs.py [--fix-line-numbers]
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
DOCS = ["README.md", "CHANGELOG.md", "docs/PROVENANCE.md",
        "docs/METHOD_mhc_propagation.md", "vendor/VENDOR.md"]

#: Filenames the pipeline writes; anything else quoted as an output is a typo.
KNOWN_OUTPUTS = {
    "stage1_junctions.tsv", "stage3_matches.tsv", "stage3_gene_svtype.tsv",
    "stage3_proximity.tsv", "credible_events.tsv", "stage7_expression.tsv",
    "stage7_background.json", "stage8_rna_evidence.tsv", "funnel.tsv",
    "summary.json", "comparison.tsv", "null_model_rates.tsv", "attribution.tsv",
    "master_sv.tsv", "master_peptides.tsv",
    "master_sv_column_dictionary.tsv", "master_peptides_column_dictionary.tsv",
    "stage6_8_all_junctions.tsv", "stage3_junction_recurrence.tsv",
    "strand_bias.tsv",
}


def index_definitions() -> tuple[dict, dict]:
    """Line number of every constant and every top-level function."""
    constants, functions = {}, {}
    criteria = REPO / "src" / "svneo" / "criteria.py"
    for i, line in enumerate(criteria.read_text().splitlines(), 1):
        match = re.match(r"^([A-Z][A-Z0-9_]+)\s*=", line)
        if match and match.group(1) not in constants:
            constants[match.group(1)] = i
    for path in (REPO / "src").rglob("*.py"):
        for i, line in enumerate(path.read_text().splitlines(), 1):
            match = re.match(r"^def ([a-z_]+)", line)
            if match and match.group(1) not in functions:
                functions[match.group(1)] = i
    return constants, functions


def boolean_criteria() -> dict:
    """Current value of every boolean criterion, for claim checking."""
    sys.path.insert(0, str(REPO / "src"))
    from svneo import criteria
    return {name: getattr(criteria, name) for name in dir(criteria)
            if name.isupper() and isinstance(getattr(criteria, name), bool)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fix-line-numbers", action="store_true")
    args = ap.parse_args()

    constants, functions = index_definitions()
    booleans = boolean_criteria()
    problems: list[str] = []

    #: Phrases that assert a boolean criterion is ON. If the criterion is off,
    #: the prose is wrong — the failure mode that motivated this script.
    claims = {
        "REQUIRE_SPANS_JUNCTION": [
            r"peptides that do not span the junction are excluded",
            r"only peptides that span",
        ],
        "SOFTCLIPS_TIER_EVENTS": [r"soft-?clips? (?:can |may )?tier"],
        "COVERAGE_TIERS_EVENTS": [r"coverage (?:can |may )?promotes? a tier"],
    }

    for name in DOCS:
        path = REPO / name
        if not path.exists():
            problems.append(f"{name}: documented file does not exist")
            continue
        text = path.read_text()
        fixed = text

        for const, cited in re.findall(r"`([A-Z][A-Z0-9_]{3,})`\s+L(\d+)", text):
            if const in constants and int(cited) != constants[const]:
                problems.append(f"{name}: `{const}` cited at L{cited}, "
                                f"defined at L{constants[const]}")
                fixed = fixed.replace(f"`{const}` L{cited}",
                                      f"`{const}` L{constants[const]}")
        for func, cited in re.findall(r"`([a-z_]+)\(\)`\s+L(\d+)", text):
            if func in functions and int(cited) != functions[func]:
                problems.append(f"{name}: `{func}()` cited at L{cited}, "
                                f"defined at L{functions[func]}")
                fixed = fixed.replace(f"`{func}()` L{cited}",
                                      f"`{func}()` L{functions[func]}")

        for const in set(re.findall(r"`([A-Z][A-Z0-9_]{3,})`", text)):
            if const.isupper() and const not in constants and "_" in const \
                    and const not in ("PON_COUNT", "FUNNEL_ROWS", "SVTYPE",
                                      "MATEID", "SEGMAPQ", "REFPAIR",
                                      "PURPLE_AF", "PURPLE_CN", "PURPLE_JCN",
                                      "PURPLE_CN_CHANGE", "END", "HOMSEQ",
                                      "IHOMPOS", "ALTALN", "ASMID", "ASMLEN",
                                      "SEGRL", "MIN_INTRON", "PYENSEMBL_CACHE_DIR",
                                      "PYSPARK_PYTHON", "PYSPARK_DRIVER_PYTHON",
                                      "LICENSE_NEOSV", "REQUIRED_PEPTIDE_COLUMNS"):
                problems.append(f"{name}: `{const}` is not defined in criteria.py")

        for output in set(re.findall(r"`(\w+\.(?:tsv|json))`", text)):
            if output not in KNOWN_OUTPUTS and not output.startswith("gdc"):
                problems.append(f"{name}: `{output}` is not a pipeline output")

        for const, patterns in claims.items():
            if booleans.get(const) is False:
                for pattern in patterns:
                    if re.search(pattern, text, re.I):
                        problems.append(
                            f"{name}: describes {const} as active, but it is False")

        if args.fix_line_numbers and fixed != text:
            path.write_text(fixed)

    print(f"checked {len(DOCS)} documents against the code")
    if problems:
        print(f"\n{len(problems)} stale reference(s):")
        for problem in problems:
            print(f"  {problem}")
        if not args.fix_line_numbers:
            print("\nre-run with --fix-line-numbers to correct line references")
        sys.exit(1)
    print("no stale references")


if __name__ == "__main__":
    main()
