#!/usr/bin/env python
"""
provenance.py — stamp every generated file with the code that produced it.

WHY
---
Each generated table already records its inputs and their checksums, which
answers "what went in". It did not record what the code was, which leaves the
other half of the question open: a table and the repository can disagree, and
nothing says whether the table predates a fix or follows it.

`git rev-parse HEAD` closes that, on one condition — the working tree has to be
clean. A commit id on a table built from edited-but-uncommitted files is worse
than no id at all, because it names a version that never produced it. So the
dirty state is recorded explicitly rather than hidden, and anything built from
a dirty tree says so in its own provenance.
"""
from __future__ import annotations

import pathlib
import subprocess
import sys
from datetime import datetime, timezone

REPO = pathlib.Path(__file__).resolve().parent.parent


def _git(*args: str) -> str | None:
    try:
        out = subprocess.run(("git", "-C", str(REPO)) + args,
                             capture_output=True, text=True, timeout=10)
        return out.stdout.strip() if out.returncode == 0 else None
    except (OSError, subprocess.SubprocessError):
        return None


def code_version() -> list[dict]:
    """Provenance rows naming the code, ready to concatenate with the caller's."""
    commit = _git("rev-parse", "HEAD")
    if commit is None:
        return [{"key": "code_commit", "value": "not a git checkout"}]

    dirty = _git("status", "--porcelain") or ""
    rows = [
        {"key": "code_commit", "value": commit},
        {"key": "code_branch", "value": _git("rev-parse", "--abbrev-ref", "HEAD") or ""},
        {"key": "code_committed", "value": _git("log", "-1", "--format=%cI") or ""},
        {"key": "code_tree_clean", "value": str(not dirty)},
        {"key": "generated_utc",
         "value": datetime.now(timezone.utc).isoformat(timespec="seconds")},
        {"key": "python", "value": sys.version.split()[0]},
    ]
    if dirty:
        # Name the files, not just the fact. Without them "dirty" is a warning
        # nobody can act on six months later.
        changed = sorted(line[3:] for line in dirty.splitlines())
        rows.append({"key": "code_tree_modified",
                     "value": ", ".join(changed[:20])
                              + (f" (+{len(changed) - 20} more)" if len(changed) > 20 else "")})
        rows.append({"key": "WARNING", "value":
                     "Built from a MODIFIED working tree: code_commit does not "
                     "describe the code that ran. Commit before a run whose "
                     "output you intend to keep."})
    return rows


def package_versions(*names: str) -> list[dict]:
    """Versions of the packages a tool actually used."""
    rows = []
    for name in names:
        try:
            module = __import__(name)
            rows.append({"key": f"version_{name}",
                         "value": getattr(module, "__version__", "unknown")})
        except ImportError:
            rows.append({"key": f"version_{name}", "value": "not installed"})
    return rows


if __name__ == "__main__":
    for row in code_version():
        print(f"{row['key']}\t{row['value']}")
