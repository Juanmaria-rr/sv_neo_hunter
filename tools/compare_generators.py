#!/usr/bin/env python3
"""
compare_generators.py — are two peptide generators interchangeable?

WHY THIS MATTERS
----------------
A recurrence test compares peptide strings between a reference catalogue and a
sample. That comparison is only valid if both sides were produced by generators
that emit the same peptides. Swapping a generator — for a newer version, a
different licence situation, or a different tool — therefore needs evidence, not
assumption: a single changed residue in the sliding-window logic silently breaks
every identity match.

This script runs two generator packages over the same admitted VCF, up to and
including peptide generation, and compares the resulting peptide sets per SV and
in aggregate.

USAGE
-----
    python tools/compare_generators.py \\
        --vcf work/sample.admitted.vcf \\
        --package-a /path/to/NeoSV-main   --module-a neosv \\
        --package-b ../vendor             --module-b neosv_trace \\
        --release 115 --cache-dir /path/to/pyensembl_cache

Both packages are imported into one process under their own module names, so
they must not share a package name.

INTERPRETING THE RESULT
-----------------------
`identical`      the generators are interchangeable for this input
`a_only/b_only`  peptides one produces and the other does not — any non-empty
                 set means catalogues built with the two are NOT comparable
"""
from __future__ import annotations

import argparse
import importlib
import json
import sys


def load_package(path: str, module: str):
    """Import a generator package from an arbitrary directory."""
    if path not in sys.path:
        sys.path.insert(0, path)
    for name in list(sys.modules):
        if name == module or name.startswith(f"{module}."):
            del sys.modules[name]
    return importlib.import_module(module)


def peptides_from(package, vcf_path: str, release: str, cache_dir: str | None,
                  lengths: str = "8-11") -> dict:
    """Run one generator up to peptide generation. Returns {sv_key: {peptides}}.

    Mirrors the tool's own main() ordering: load -> infer pattern -> dedup ->
    build fusions -> set sequences -> generate peptides. Stops before any MHC
    step, which is not part of peptide generation.
    """
    inp = importlib.import_module(f"{package.__name__}.input")
    sv_utils = importlib.import_module(f"{package.__name__}.sv_utils")
    fusion_utils = importlib.import_module(f"{package.__name__}.fusion_utils")
    seq_utils = importlib.import_module(f"{package.__name__}.sequence_utils")

    ensembl = inp.ensembl_load(str(release), None, None, cache_dir)
    svs = [sv_utils.sv_pattern_infer_vcf(v) for v in inp.vcf_load(vcf_path)]
    svs = sv_utils.remove_duplicate(svs)
    window = inp.get_window_range(lengths)

    out: dict = {}
    for sv in svs:
        fusion = fusion_utils.sv_to_svfusion(sv, ensembl)
        if fusion.is_empty():
            continue
        fusion.nt_sequence = seq_utils.set_nt_seq(fusion)
        fusion.aa_sequence = seq_utils.set_aa_seq(fusion)
        peptides = set(seq_utils.generate_neoepitopes(fusion, window))
        # Key on coordinates, not on any id: only one of the two packages
        # necessarily carries an sv_id, and the comparison must not depend on it.
        key = f"{sv.chrom1}:{sv.pos1}-{sv.chrom2}:{sv.pos2}"
        out.setdefault(key, set()).update(peptides)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, help="admitted VCF (plain, not .gz)")
    ap.add_argument("--package-a", required=True)
    ap.add_argument("--module-a", default="neosv")
    ap.add_argument("--package-b", required=True)
    ap.add_argument("--module-b", default="neosv_trace")
    ap.add_argument("--release", default="115")
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args()

    a = peptides_from(load_package(args.package_a, args.module_a),
                      args.vcf, args.release, args.cache_dir)
    b = peptides_from(load_package(args.package_b, args.module_b),
                      args.vcf, args.release, args.cache_dir)

    all_a = {p for s in a.values() for p in s}
    all_b = {p for s in b.values() for p in s}
    per_sv_differences = {k: {"a_only": sorted(a.get(k, set()) - b.get(k, set())),
                              "b_only": sorted(b.get(k, set()) - a.get(k, set()))}
                          for k in set(a) | set(b)
                          if a.get(k, set()) != b.get(k, set())}

    result = {
        "module_a": args.module_a, "module_b": args.module_b,
        "svs_a": len(a), "svs_b": len(b),
        "peptides_a": len(all_a), "peptides_b": len(all_b),
        "shared": len(all_a & all_b),
        "a_only": sorted(all_a - all_b), "b_only": sorted(all_b - all_a),
        "svs_with_differences": len(per_sv_differences),
        "identical": all_a == all_b and not per_sv_differences,
    }

    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("a_only", "b_only")}, indent=2))
    if result["identical"]:
        print("\nIDENTICAL — the generators are interchangeable for this input, "
              "so catalogues built with either are comparable.")
    else:
        print(f"\nDIFFERENT — {len(result['a_only'])} peptides only in "
              f"{args.module_a}, {len(result['b_only'])} only in {args.module_b}.")
        print("Catalogues built with the two are NOT comparable by peptide identity.")
        for peptide in result["a_only"][:5]:
            print(f"  only in {args.module_a}: {peptide}")
        for peptide in result["b_only"][:5]:
            print(f"  only in {args.module_b}: {peptide}")

    if args.json_out:
        with open(args.json_out, "w") as handle:
            json.dump(result, handle, indent=2)
        print(f"\nwrote {args.json_out}")

    sys.exit(0 if result["identical"] else 1)


if __name__ == "__main__":
    main()
