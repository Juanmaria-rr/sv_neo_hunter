#!/usr/bin/env python
"""
lineage_args.py — turn a run config into the arguments the universe build needs.

The lineage is a property of the samples, and the samples are named in the run
config, which is not in this repository because those names identify the
material. Deriving the arguments here keeps them out of the tracked pipeline
script: `tools/run_universe.sh` never contains a sample name.

Emits shell-quoted arguments on one line, for `eval` or command substitution.

Usage
-----
    python tools/lineage_args.py --config config/<run>.yaml \\
        --cross <dir> --branch noPON --what lines
"""
from __future__ import annotations

import argparse
import pathlib
import shlex

import yaml


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=pathlib.Path, required=True)
    parser.add_argument("--cross", required=True,
                        help="directory holding the per-sample pipeline runs")
    parser.add_argument("--branch", default="noPON",
                        help="branch suffix on each run directory")
    parser.add_argument("--what", choices=("lines", "samples", "from", "names"),
                        default="lines")
    args = parser.parse_args()

    config = yaml.safe_load(args.config.read_text())
    samples = config.get("samples") or []
    if not samples:
        raise SystemExit(f"  no samples in {args.config}")

    out = []
    for sample in samples:
        name = sample["name"]
        run_dir = f"{args.cross}/{name}_{args.branch}"
        if args.what == "lines":
            out += ["--line", f"{name}:{run_dir}:{sample.get('parent') or ''}"]
        elif args.what == "samples":
            out += ["--sample", f"{name}_{args.branch}"]
        elif args.what == "from":
            out += ["--from", run_dir]
        else:
            out.append(name)
    print(" ".join(shlex.quote(x) for x in out))


if __name__ == "__main__":
    main()
