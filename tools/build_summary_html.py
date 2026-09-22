#!/usr/bin/env python3
"""
build_summary_html.py — render a repository document as a shareable web page.

ONE SOURCE, TWO OUTPUTS
-----------------------
The markdown in `docs/` is canonical: it is versioned, reviewed in diffs, and
checked against the code by `tools/check_docs.py`. This script renders it to a
self-contained HTML file for people who will read the analysis but not clone the
repository.

The page is generated, never edited. A hand-edited HTML copy is a second source
that drifts from the first, and the reader has no way to tell which one is
current — the failure this repository has already had once, when a conclusions
document outlived the code that produced its numbers.

    python tools/build_summary_html.py
    python tools/build_summary_html.py --doc <your-summary>.md
"""
from __future__ import annotations

import argparse
import html
import pathlib
import re

import markdown

REPO = pathlib.Path(__file__).resolve().parent.parent

# ============================================================================
# Design
#
# A reference document, read in long sittings and returned to: the treatment is
# a printed technical report, not a landing page. Serif for running text because
# the prose carries the argument; a monospaced face for every measured figure so
# a number is visibly a number. The one accent is reserved for thresholds and
# criteria names, which is what a reader scans for when checking a claim.
# ============================================================================

CSS = """
:root{
  --ground:#FBFAF7; --surface:#FFFFFF; --ink:#1A1D21; --ink-soft:#4A5158;
  --ink-faint:#868D95; --rule:#E2DFD8; --rule-soft:#EFEDE8;
  --accent:#8A4B2A; --accent-soft:#F6EDE6;
  --steel:#2E5C7E; --moss:#3A6B4A; --amber:#9A6B18;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --ground:#15171A; --surface:#1C1F23; --ink:#E8E6E1; --ink-soft:#B0B6BD;
    --ink-faint:#767D85; --rule:#2C3137; --rule-soft:#23272C;
    --accent:#D89468; --accent-soft:#2B2018;
    --steel:#7FAFD4; --moss:#7FBE93; --amber:#D6A548;
  }
}
:root[data-theme="dark"]{
  --ground:#15171A; --surface:#1C1F23; --ink:#E8E6E1; --ink-soft:#B0B6BD;
  --ink-faint:#767D85; --rule:#2C3137; --rule-soft:#23272C;
  --accent:#D89468; --accent-soft:#2B2018;
  --steel:#7FAFD4; --moss:#7FBE93; --amber:#D6A548;
}

*{box-sizing:border-box}
body{
  background:var(--ground); color:var(--ink);
  font-family:"Iowan Old Style","Palatino Linotype",Palatino,Georgia,serif;
  font-size:17px; line-height:1.62; margin:0;
  -webkit-font-smoothing:antialiased;
}
.wrap{max-width:60rem; margin:0 auto; padding:4rem 1.5rem 6rem}

h1{
  font-size:2.1rem; line-height:1.18; margin:0 0 .6rem; text-wrap:balance;
  letter-spacing:-.015em; font-weight:600;
}
h2{
  font-size:1.42rem; margin:3.4rem 0 1rem; text-wrap:balance; font-weight:600;
  padding-bottom:.45rem; border-bottom:2px solid var(--rule);
}
h3{
  font-size:1.1rem; margin:2.3rem 0 .7rem; text-wrap:balance; font-weight:600;
  color:var(--accent);
}
h4{font-size:.97rem; margin:1.6rem 0 .5rem; font-weight:600; color:var(--ink-soft)}
p{margin:0 0 1rem; max-width:68ch}
ul,ol{margin:0 0 1rem; padding-left:1.35rem; max-width:68ch}
li{margin-bottom:.35rem}
li>p{margin-bottom:.4rem}
strong{font-weight:600}
em{font-style:italic; color:var(--ink-soft)}
a{color:var(--steel); text-decoration:underline; text-underline-offset:2px}
hr{border:0; border-top:1px solid var(--rule); margin:2.6rem 0}

code{
  font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
  font-size:.855em; background:var(--accent-soft); color:var(--accent);
  padding:.09em .34em; border-radius:3px; word-break:break-word;
}
pre{
  background:var(--surface); border:1px solid var(--rule); border-radius:4px;
  padding:.9rem 1rem; overflow-x:auto; margin:0 0 1.2rem;
}
pre code{background:none; color:var(--ink-soft); padding:0; font-size:.82rem;
  line-height:1.55}

blockquote{
  margin:0 0 1.2rem; padding:.75rem 1.1rem; border-left:3px solid var(--amber);
  background:var(--surface); color:var(--ink-soft); border-radius:0 3px 3px 0;
}
blockquote p:last-child{margin-bottom:0}

.tablewrap{overflow-x:auto; margin:0 0 1.4rem}
table{
  border-collapse:collapse; width:100%; font-size:.86rem;
  font-variant-numeric:tabular-nums;
  font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
}
th{
  text-align:left; font-weight:600; color:var(--ink-soft); white-space:nowrap;
  border-bottom:2px solid var(--rule); padding:.45rem .8rem .45rem 0;
  vertical-align:bottom;
}
td{padding:.38rem .8rem .38rem 0; border-bottom:1px solid var(--rule-soft);
  vertical-align:top}
td:first-child,th:first-child{font-family:inherit}
tr:last-child td{border-bottom:none}
table code{font-size:.9em}

/* Table of contents, built from the h2/h3 structure. */
.toc{
  background:var(--surface); border:1px solid var(--rule); border-radius:5px;
  padding:1.1rem 1.4rem; margin:2.2rem 0 3rem; font-size:.9rem;
}
.toc p{font-weight:600; margin:0 0 .55rem; font-family:inherit;
  text-transform:uppercase; letter-spacing:.06em; font-size:.72rem;
  color:var(--ink-faint)}
.toc ol{list-style:none; padding:0; margin:0; max-width:none}
.toc li{margin-bottom:.28rem}
.toc li.sub{padding-left:1.2rem; font-size:.93em}
.toc li.top{margin-top:.6rem; font-weight:600}
.toc li.top:first-child{margin-top:0}
.toc a{color:var(--ink); text-decoration:none}
.toc a:hover{color:var(--accent); text-decoration:underline}

.masthead{border-bottom:2px solid var(--ink); padding-bottom:1.2rem;
  margin-bottom:.4rem}
.masthead .kicker{
  font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
  font-size:.7rem; text-transform:uppercase; letter-spacing:.14em;
  color:var(--accent); margin:0 0 .7rem;
}
.masthead .meta{
  font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
  font-size:.74rem; color:var(--ink-faint); margin:.7rem 0 0;
}
footer{
  margin-top:4rem; padding-top:1.2rem; border-top:1px solid var(--rule);
  font-size:.78rem; color:var(--ink-faint);
  font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace;
}
:focus-visible{outline:2px solid var(--steel); outline-offset:2px}
@media (max-width:640px){
  body{font-size:16px}
  .wrap{padding:2.2rem 1.1rem 4rem}
  h1{font-size:1.65rem}
}
@media print{
  body{background:#fff; color:#000}
  .toc{break-inside:avoid}
  h2,h3{break-after:avoid}
}
"""


def slug(text: str) -> str:
    """A stable anchor id from a heading."""
    plain = re.sub(r"<[^>]+>", "", text)
    return re.sub(r"[^a-z0-9]+", "-", plain.lower()).strip("-") or "section"


def render(source: str) -> tuple[str, str, str]:
    """Markdown -> (title, table of contents, body). Headings gain anchors."""
    converted = markdown.markdown(
        source, extensions=["tables", "fenced_code", "sane_lists", "attr_list"])

    title = "Document"
    first = re.search(r"<h1>(.*?)</h1>", converted, re.S)
    if first:
        title = re.sub(r"<[^>]+>", "", first.group(1)).strip()

    entries: list[tuple[int, str, str]] = []

    def anchor(match: re.Match) -> str:
        level, text = int(match.group(1)), match.group(2)
        identifier = slug(text)
        # h1-h3 are anchored, so a cross-reference to a top-level section
        # resolves; anything deeper is detail and stays out of the contents.
        entries.append((level, identifier, re.sub(r"<[^>]+>", "", text)))
        return f'<h{level} id="{identifier}">{text}</h{level}>'

    body = re.sub(r"<h([123])>(.*?)</h\1>", anchor, converted, flags=re.S)
    # The document title is the page masthead, not a contents entry.
    entries[:] = entries[1:] if entries and entries[0][0] == 1 else entries

    # Tables scroll inside their own container: a wide cascade table must never
    # make the page itself scroll sideways.
    body = body.replace("<table>", '<div class="tablewrap"><table>')
    body = body.replace("</table>", "</table></div>")

    items = "".join(
        f'<li class="{"sub" if level == 3 else ""}">'
        f'<a href="#{identifier}">{html.escape(text)}</a></li>'
        for level, identifier, text in entries)
    toc = f'<nav class="toc"><p>Contents</p><ol>{items}</ol></nav>' if items else ""
    return title, toc, body


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--doc", required=True,
                        help="markdown source, relative to the repository root")
    parser.add_argument("--out", default=None,
                        help="output path (default: results/reports/<name>.html)")
    args = parser.parse_args()

    source_path = REPO / args.doc
    if not source_path.exists():
        raise SystemExit(f"no such document: {source_path}")
    source = source_path.read_text()

    title, toc, body = render(source)

    # The first paragraph after the title serves as the standfirst; strip it from
    # the body so it is not printed twice.
    #
    # Both patterns are anchored at the start and forbidden from crossing a
    # closing tag. A single `<h1>.*?</h1>` here is NOT safe with re.DOTALL: a
    # document with several level-1 headings (Part A, Part B, …) lets the lazy
    # quantifier run to a *later* `</h1>` when the title is not immediately
    # followed by a paragraph, silently swallowing the body. That deleted 47,550
    # of 49,774 characters and produced a page that still looked plausible.
    lead = ""
    title_match = re.match(r"\s*<h1>(?:(?!</h1>).)*</h1>\s*", body, re.S)
    if title_match:
        body = body[title_match.end():]
    lead_match = re.match(r"<p>(?:(?!</p>).)*</p>\s*", body, re.S)
    if lead_match:
        lead = re.match(r"<p>(.*)</p>\s*", lead_match.group(0), re.S).group(1)
        body = body[lead_match.end():]

    revision = "uncommitted"
    try:
        import subprocess
        revision = subprocess.run(
            ["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, check=True).stdout.strip() or revision
    except Exception:
        pass

    page = f"""<title>{html.escape(title)}</title>
<style>{CSS}</style>
<div class="wrap">
  <header class="masthead">
    <p class="kicker">svneo &middot; technical report</p>
    <h1>{html.escape(title)}</h1>
    {f"<p>{lead}</p>" if lead else ""}
    <p class="meta">Generated from <code>{html.escape(args.doc)}</code> at
    revision <code>{html.escape(revision)}</code>. The markdown in the repository
    is canonical; this page is rendered from it and is never edited directly.</p>
  </header>
  {toc}
  {body}
  <footer>
    svneo &middot; rendered by <code>tools/build_summary_html.py</code> from
    <code>{html.escape(args.doc)}</code> &middot; revision
    <code>{html.escape(revision)}</code>
  </footer>
</div>"""

    out = pathlib.Path(args.out) if args.out else \
        REPO / "results" / "reports" / (source_path.stem + ".html")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page)
    print(f"  wrote {out.relative_to(REPO)}  ({len(page) / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
