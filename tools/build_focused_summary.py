#!/usr/bin/env python3
"""
build_focused_summary.py — the executive summary for a single branch.

WHY THIS IS GENERATED AND NOT WRITTEN
-------------------------------------
The focused summary drops the panel-filtered branch and the permutation null, so
that a reader who only wants "what did the full call set show" is not navigating
a branch comparison. Everything it says is already in
the source document, which stays canonical.

Writing it by hand would create a second document asserting the same thresholds,
and the two would diverge — which is not hypothetical here: the canonical summary
once asserted that the panel branch did not apply to the derived lines, months
after that stopped being true. So this script derives it instead, and the derived
file carries a header saying so.

HOW THE CANONICAL DOCUMENT MARKS WHAT IS BRANCH-SPECIFIC
--------------------------------------------------------
    <!-- focused:drop -->    ...   <!-- /focused:drop -->
        a block that exists only to compare branches, or to report the null

Tables that report run outputs are generated separately, by
`tools/sync_event_tables.py`, so both documents inherit the same machine-
maintained figures.

Table columns whose header names the dropped branch are removed automatically,
so a two-branch cascade table becomes a one-branch one without being rewritten.

    python tools/build_focused_summary.py
    python tools/build_focused_summary.py --drop-branch PON10 --keep-branch noPON
"""
from __future__ import annotations

import argparse
import pathlib
import re
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent

#: Two forms, applied in this order.
#:
#: BLOCK — markers alone on their lines: the whole span, newlines included, goes.
#: INLINE — markers inside a line: only the span between them goes, and the
#: line's own newline is left alone. Without the distinction, an inline fence at
#: the end of a bullet swallows the newline and welds the next bullet onto it,
#: which is how a list silently becomes a paragraph.
DROP_BLOCK = re.compile(
    r"^[ \t]*<!--\s*focused:drop\s*-->[ \t]*\n"
    r".*?"
    r"^[ \t]*<!--\s*/focused:drop\s*-->[ \t]*\n",
    re.S | re.M)
DROP_INLINE = re.compile(
    r"[ \t]*<!--\s*focused:drop\s*-->.*?<!--\s*/focused:drop\s*-->", re.S)

def drop_branch_columns(text: str, branch: str) -> str:
    """Remove every markdown table column whose header names `branch`."""
    out, table, in_table = [], [], False

    def flush() -> None:
        if not table:
            return
        rows = [[c.strip() for c in line.strip().strip("|").split("|")]
                for line in table]
        width = len(rows[0])
        keep = [i for i, head in enumerate(rows[0])
                if branch.lower() not in head.lower()]
        if len(keep) == width:
            out.extend(table)                       # nothing to drop
        else:
            for row in rows:
                if len(row) != width:               # ragged: leave it alone
                    out.extend(table)
                    break
            else:
                for row in rows:
                    out.append("| " + " | ".join(row[i] for i in keep) + " |\n")
        table.clear()

    for line in text.splitlines(keepends=True):
        stripped = line.strip()
        if stripped.startswith("|") and stripped.endswith("|"):
            table.append(line)
            in_table = True
            continue
        if in_table:
            flush()
            in_table = False
        out.append(line)
    flush()
    return "".join(out)


def strip_branch_suffix(text: str, branch: str) -> str:
    """Drop the kept branch's name from sample identifiers.

    In a one-branch document `parental_noPON` is just parental, and carrying the
    suffix invites the opposite of the intended reading: a suffixed name looks
    like a *variant* of the sample, as if there were another parental elsewhere
    with different data. There is not — this document describes one call set.

    Only identifiers are touched: `<sample>_noPON`, `<sample><br>noPON`,
    `<sample> (noPON)` and the parenthesised heading form. Standalone prose
    mentions survive, because Part A still has to explain what reporting the
    panel without filtering on it means, and that explanation names the setting.
    """
    escaped = re.escape(branch)
    patterns = [
        (rf"(?<=[A-Za-z0-9])_{escaped}\b", ""),          # parental_noPON
        (rf"<br>\s*`?{escaped}`?", ""),                  # WT<br>noPON
        (rf"\s*\(\s*`?{escaped}`?\s*\)", ""),            # ... (`noPON`)
        (rf"(?<=[A-Za-z0-9])\s+{escaped}\b(?=\s*\|)", ""),  # `parental noPON |`
        # A header cell left holding only the branch name, once the other
        # branch's column has gone. `| Step | noPON |` names a comparison that
        # no longer exists in this document; the column is simply the count.
        (rf"\|\s*`?{escaped}`?\s*\|", "| n |"),
    ]
    for pattern, replacement in patterns:
        text = re.sub(pattern, replacement, text)
    return text


def drop_section(text: str, title: str) -> tuple[str, bool]:
    """Remove a whole section, from its heading to the next one of equal or
    higher level.

    Matched on the heading text rather than on a marker, so a caller can prune a
    document without editing it — the source stays the one place the section is
    written. Returns the text and whether anything was found, so a typo in a
    section name fails loudly instead of silently pruning nothing.
    """
    pattern = re.compile(rf"^(#{{1,6}})\s+{re.escape(title)}\s*$", re.M)
    match = pattern.search(text)
    if not match:
        return text, False
    level = len(match.group(1))
    following = re.compile(rf"^#{{1,{level}}}\s+", re.M)
    nxt = following.search(text, match.end())
    end = nxt.start() if nxt else len(text)
    return text[:match.start()] + text[end:], True


#: Marks the provenance note this script inserts, so re-deriving a document that
#: already carries one replaces it instead of stacking a second.
NOTE_MARKER = "**Generated document — do not edit.**"


def strip_existing_note(text: str) -> str:
    """Remove a provenance note left by a previous derivation."""
    return re.sub(
        rf"\n*^> .*?{re.escape(NOTE_MARKER)}.*?(?=\n[^>\n])|"
        rf"\n*(?:^>.*\n)*?^>.*{re.escape(NOTE_MARKER)}.*\n(?:^>.*\n)*",
        "\n\n", text, count=1, flags=re.M | re.S)


def collapse_blank_runs(text: str) -> str:
    """Three or more blank lines become one, after blocks are removed."""
    return re.sub(r"\n{4,}", "\n\n\n", text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True,
                        help="the summary document to condense")
    parser.add_argument("--out", required=True)
    parser.add_argument("--drop-branch", default="PON10")
    parser.add_argument("--keep-branch", default="noPON")
    parser.add_argument("--drop-section", action="append", default=[],
                        metavar="TITLE",
                        help="remove a whole section by its heading text; "
                             "repeatable. Fails if the heading is not found.")
    parser.add_argument("--no-html", action="store_true",
                        help="write the markdown only")
    args = parser.parse_args()

    source = REPO / args.source
    if not source.exists():
        raise SystemExit(f"no such document: {source}")
    text = source.read_text()

    dropped = len(DROP_BLOCK.findall(text)) + len(DROP_INLINE.findall(text))
    text = DROP_BLOCK.sub("", text)
    text = DROP_INLINE.sub("", text)
    text = drop_branch_columns(text, args.drop_branch)
    text = strip_branch_suffix(text, args.keep_branch)

    missing = []
    for title in args.drop_section:
        text, found = drop_section(text, title)
        if not found:
            missing.append(title)
    if missing:
        # A section name that matches nothing must not pass quietly: the caller
        # asked for it to be gone, and a silent no-op ships it anyway.
        raise SystemExit(
            "these sections were not found in " + args.source + ":\n  "
            + "\n  ".join(missing)
            + "\n\nheadings available:\n  "
            + "\n  ".join(re.findall(r"^#{1,3} (.+)$", text, re.M)))

    text = strip_existing_note(text)
    text = collapse_blank_runs(text)

    # Say what this file is, in the file. A derived document that does not
    # announce itself gets edited, and the edit is lost on the next build.
    title_end = text.index("\n", text.index("# "))
    omitted = ["the panel-filtered comparison", "the permutation null"]
    omitted += [f"\u201c{title}\u201d" for title in args.drop_section]
    if len(omitted) > 1:
        omissions = ", ".join(omitted[:-1]) + " and " + omitted[-1]
    else:
        omissions = omitted[0]

    canonical = args.source
    note = (
        f"\n\n> **This describes the full SV call set**, with the panel of "
        f"normals measured and reported at every step but never used to remove a "
        f"candidate. Each sample therefore appears once, under its own name.\n"
        f">\n"
        f"> **{NOTE_MARKER[2:-2]}** Built by `tools/build_focused_summary.py` "
        f"from [`{args.source}`]({pathlib.Path(args.source).name}). "
        f"[`{canonical}`]({pathlib.Path(canonical).name}) is the canonical "
        f"version and additionally carries {omissions}. Edit the source and "
        f"rebuild.\n"
        f">\n"
        f"> No *result* is omitted here: the panel-filtered call set is a strict "
        f"subset of this one, so every event it reports appears below.\n")
    text = text[:title_end + 1] + note + text[title_end + 1:]

    out = REPO / args.out
    out.write_text(text)
    print(f"  wrote {out.relative_to(REPO)}  "
          f"({len(text.splitlines()):,} lines, {dropped} block(s) dropped)")

    if not args.no_html:
        subprocess.run(
            [sys.executable, str(REPO / "tools" / "build_summary_html.py"),
             "--doc", args.out], check=True)


if __name__ == "__main__":
    main()
