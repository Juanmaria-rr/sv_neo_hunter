#!/usr/bin/env python
"""
validate_junction.py — is a called junction a genomic lesion or a splicing event?

THE QUESTION THIS ANSWERS, AND THE ONE IT REPLACES
--------------------------------------------------
An SV whose breakpoints sit inside a transcribed gene produces RNA reads with a
gap at the deletion, and so does an alternative splice junction. Read counts
alone cannot tell them apart, and neither can the ratio of RNA support to DNA
support: at a gene expressed at hundreds of TPM the RNA depth exceeds the DNA
depth by orders of magnitude whether or not the lesion is real. That ratio looks
like a discriminator and is not one — it is a measure of expression.

**The discriminator is DNA copy number across the interval.** A heterozygous
deletion removes one of two copies, so read depth inside the deleted span falls
to about half of the flanks. A splicing event does not touch the DNA at all, so
the ratio stays near one. Nothing about expression enters this comparison.

Two further things are worth seeing at the same time, and this prints both:

  - the SIZES of the N gaps at the locus, which separate the junction's own
    population from the gene's ordinary introns, and
  - where each gap STARTS, since a gap belonging to the junction begins at the
    breakpoint while an intron begins at its own splice site.

A tool that calls splicing will label a transcript from a deleted allele a novel
intron, because from its point of view that is exactly what it looks like. Such a
call is therefore not evidence against the lesion; it is the expected observation
either way.

Usage
-----
    python tools/validate_junction.py --rna-bam <bam> --dna-bam <bam> \\
        --chrom <chrom> --pos1 <breakend1> --pos2 <breakend2> \\
        --out <dir>/validation.html
"""
from __future__ import annotations

import argparse
import collections
import pathlib
import sys

import numpy as np
import pysam

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))
from svneo import criteria as C                           # noqa: E402

#: Reads below this are ambiguously placed and are not evidence for a locus.
MIN_MAPQ = 20
#: Flank used for the copy-number comparison, and the margin trimmed from the
#: interval's own edges so partly-overlapping reads do not blur the boundary.
FLANK_BP, EDGE_TRIM = 500, 20


def contig(bam: pysam.AlignmentFile, chrom: str) -> str:
    return chrom if chrom in bam.references else f"chr{chrom}"


def gap_profile(bam_path: str, chrom: str, pos1: int, pos2: int,
                window: int = 400) -> tuple[collections.Counter, collections.Counter, int]:
    """N-gap sizes and start positions around the junction."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref = contig(bam, chrom)
    sizes, starts, total = collections.Counter(), collections.Counter(), 0
    for read in bam.fetch(ref, max(0, pos1 - window), pos2 + window):
        if read.is_unmapped or not read.cigartuples:
            continue
        if read.mapping_quality < MIN_MAPQ:
            continue
        total += 1
        cursor = read.reference_start
        for op, length in read.cigartuples:
            if op == 3:                      # N — reference skipped
                sizes[length] += 1
                starts[cursor] += 1
            if op in (0, 2, 3, 7, 8):        # reference-consuming
                cursor += length
    return sizes, starts, total


def junction_usage(bam_path: str, chrom: str, pos1: int, pos2: int,
                   window: int = 400) -> dict:
    """What fraction of transcripts crossing this interval use the junction?

    THE DENOMINATOR IS THE WHOLE POINT
    ----------------------------------
    A raw junction count cannot be interpreted without knowing what it is a
    count out of. Seven thousand reads is overwhelming support in a gene at
    5 TPM and unremarkable in one at 900 TPM, so any threshold placed on the
    raw number is in practice a filter on expression.

    Two candidate denominators are wrong and one is right:

      - **every read at the locus** is wrong, because a read lying entirely
        inside an exon never had the opportunity to use either junction.
        Counting non-voters in a vote share measures exon length.
      - **every gapped read at the locus** is wrong, because a neighbouring
        intron spliced out of every transcript regardless does not compete
        with anything; including it only dilutes.
      - **gapped reads whose gap overlaps the SV interval** is right: those
        are the reads that could have gone either way.

    `alt` here EXCLUDES the junction's own reads, so it is the competing
    outcome on its own and is deliberately not the denominator — dividing by
    it would give odds rather than a share. `total` is the sum, and it is the
    denominator precisely because it contains the numerator.

    Fragments, not alignment records: a read pair spanning the junction is one
    observation however many records it produced.
    """
    size = pos2 - pos1
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref = contig(bam, chrom)

    hit, alt = set(), set()
    alt_sizes: collections.Counter = collections.Counter()
    for read in bam.fetch(ref, max(0, pos1 - window), pos2 + window):
        if read.is_unmapped or not read.cigartuples:
            continue
        if read.mapping_quality < MIN_MAPQ:
            continue
        cursor = read.reference_start
        for op, length in read.cigartuples:
            if op == 3:                                   # N — reference skipped
                start, end = cursor, cursor + length
                if start < pos2 and end > pos1:           # overlaps the interval
                    matches = (abs(length - size) <= C.NGAP_SIZE_TOLERANCE
                               and abs(start - pos1) <= C.NGAP_POSITION_TOLERANCE)
                    if matches:
                        hit.add(read.query_name)
                    else:
                        alt.add(read.query_name)
                        alt_sizes[length] += 1
            if op in (0, 2, 3, 7, 8):                     # reference-consuming
                cursor += length

    # A fragment carrying both a matching and a non-matching gap is the
    # junction's: it demonstrably used it.
    alt -= hit
    total = len(hit) + len(alt)
    return {"junction_reads": len(hit), "alt_reads": len(alt),
            "total": total, "usage": len(hit) / total if total else float("nan"),
            "alt_sizes": alt_sizes}


def copy_number(bam_path: str, chrom: str, pos1: int, pos2: int) -> dict:
    """Depth inside the interval against its flanks — the actual discriminator."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref = contig(bam, chrom)

    def depth(start: int, end: int) -> np.ndarray:
        counts = bam.count_coverage(ref, max(0, start), end, quality_threshold=0)
        return np.array(counts).sum(axis=0)

    left = depth(pos1 - FLANK_BP - 100, pos1 - 100)
    inside = depth(pos1 + EDGE_TRIM, pos2 - EDGE_TRIM)
    right = depth(pos2 + 100, pos2 + FLANK_BP + 100)
    flank = (float(left.mean()) + float(right.mean())) / 2
    return {"left": float(left.mean()), "inside": float(inside.mean()),
            "right": float(right.mean()), "flank": flank,
            "ratio": float(inside.mean()) / flank if flank else float("nan")}


def verdict(ratio: float) -> tuple[str, str]:
    if ratio != ratio:
        return "undetermined", "no DNA coverage to compare"
    if ratio < 0.25:
        return "genomic, homozygous", "depth inside is near zero"
    if ratio < 0.70:
        return "genomic, heterozygous", "depth inside is about half the flanks"
    if ratio > 0.85:
        return "NOT a genomic deletion", "depth inside matches the flanks"
    return "ambiguous", "depth is between the two expectations"


def bar(value: float, largest: float, width: int = 46) -> str:
    return "#" * max(1, int(round(width * value / largest))) if value else ""


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rna-bam", required=True)
    parser.add_argument("--dna-bam", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--pos1", type=int, required=True)
    parser.add_argument("--pos2", type=int, required=True)
    parser.add_argument("--label", default="")
    parser.add_argument("--out", type=pathlib.Path,
                        help="write an HTML figure here as well as printing")
    args = parser.parse_args()

    size = args.pos2 - args.pos1
    print(f"  junction {args.chrom}:{args.pos1}-{args.pos2}  ({size} bp)"
          + (f"  {args.label}" if args.label else ""))

    sizes, starts, total = gap_profile(args.rna_bam, args.chrom, args.pos1, args.pos2)
    print(f"\n  RNA: {total:,} reads at MAPQ >= {MIN_MAPQ}, "
          f"{sum(sizes.values()):,} carrying an N gap")
    top = sizes.most_common(6)
    largest = top[0][1] if top else 1
    print(f"\n  {'gap size':>10}  {'reads':>8}  {'':46}")
    for length, count in top:
        mark = "  <- the junction" if abs(length - size) <= 10 else ""
        print(f"  {length:>8} bp  {count:>8,}  {bar(count, largest)}{mark}")

    print(f"\n  {'gap starts at':>14}  {'reads':>8}   offset from breakpoint")
    for position, count in starts.most_common(5):
        delta = position - args.pos1
        mark = "  <- the junction" if abs(delta) <= 15 else ""
        print(f"  {position:>14}  {count:>8,}   {delta:+d} bp{mark}")

    use = junction_usage(args.rna_bam, args.chrom, args.pos1, args.pos2)
    print(f"\n  JUNCTION USAGE — normalised for expression, which the raw count is not")
    print(f"  {'junction_reads':>24}  {use['junction_reads']:>8,}   cross the interval BY the junction")
    print(f"  {'junction_alt_reads':>24}  {use['alt_reads']:>8,}   cross it WITHOUT the junction")
    print(f"  {'junction_interval_total':>24}  {use['total']:>8,}   the sum — the denominator")
    if use["total"]:
        print(f"  {'junction_usage':>24}  {use['usage']:>8.1%}   ~50% if one deleted allele "
              f"is transcribed like the intact one")

    cn = copy_number(args.dna_bam, args.chrom, args.pos1, args.pos2)
    call, why = verdict(cn["ratio"])
    print(f"\n  DNA depth: flanks {cn['flank']:.1f}x, inside {cn['inside']:.1f}x")
    print(f"  ratio inside/flanks = {cn['ratio']:.2f}")
    print(f"  expected 0.50 if heterozygous, ~1.00 if no genomic deletion")
    print(f"\n  VERDICT: {call} — {why}")

    if args.out:
        write_html(args.out, args, size, sizes, starts, total, cn, call, why, use)
        print(f"\n  wrote {args.out}")


def write_html(path, args, size, sizes, starts, total, cn, call, why, use) -> None:
    top = sizes.most_common(6)
    largest = top[0][1] if top else 1
    HIT = ' class="hit"'
    rows = []
    for length, count in top:
        is_junction = abs(length - size) <= 10
        cls = HIT if is_junction else ""
        jbar = " j" if is_junction else ""
        what = "the junction" if is_junction else "an ordinary intron of the gene"
        width = 100 * count / largest
        rows.append(
            f'<tr{cls}><td class="num">{length:,} bp</td>'
            f'<td class="num">{count:,}</td>'
            f'<td class="barcell"><span class="track">'
            f'<span class="bar{jbar}" style="width:{width:.1f}%"></span></span></td>'
            f'<td>{what}</td></tr>')
    start_rows = []
    for position, count in starts.most_common(5):
        delta = position - args.pos1
        near = abs(delta) <= 15
        cls = HIT if near else ""
        where = "at the breakpoint" if near else "elsewhere"
        start_rows.append(
            f'<tr{cls}><td class="num">{position:,}</td>'
            f'<td class="num">{count:,}</td><td class="num">{delta:+,} bp</td>'
            f'<td>{where}</td></tr>')
    ratio_pct = max(0.0, min(1.0, cn["ratio"])) * 100
    if use["total"]:
        pct = 100 * use["usage"]
        alt_rows = "".join(
            f'<tr><td class="num">{length:,} bp</td><td class="num">{count:,}</td></tr>'
            for length, count in use["alt_sizes"].most_common(4))
        usage_block = f"""
<h2>Junction usage — the count with its denominator</h2>
<p>The raw count answers <em>does this transcript exist</em>. It cannot answer
<em>is it the dominant form</em>, because a count scales with expression: seven
thousand reads is overwhelming in a gene at 5&nbsp;TPM and unremarkable at
900&nbsp;TPM. The share below normalises that away.</p>
<div class="stack">
  <span class="seg j" style="width:{pct:.1f}%"></span>
  <span class="seg a" style="width:{100 - pct:.1f}%"></span>
</div>
<p class="legend"><span class="key j"></span>uses the junction —
<strong>{use['junction_reads']:,}</strong> ({pct:.1f}%)
&nbsp;&nbsp;<span class="key a"></span>crosses the interval another way —
<strong>{use['alt_reads']:,}</strong> ({100 - pct:.1f}%)</p>
<p>Roughly <strong>50%</strong> is what one deleted allele transcribed as readily as
the intact one would give. Materially below that is allelic imbalance or a less
stable transcript; a few per cent is the shape of alignment noise in a
highly expressed gene.</p>
<div class="scroll"><table>
<thead><tr><th>competing gap size</th><th class="num">reads</th></tr></thead>
<tbody>{alt_rows}</tbody></table></div>
<div class="card" style="margin-top:14px">
<p class="foot"><strong>Method.</strong> Numerator <code>junction_reads</code>:
fragments with a CIGAR <code>N</code> gap within {C.NGAP_SIZE_TOLERANCE}&nbsp;bp of the
event size starting within {C.NGAP_POSITION_TOLERANCE}&nbsp;bp of the breakpoint.
Denominator <code>junction_interval_total</code>: those <em>plus</em>
<code>junction_alt_reads</code> — fragments whose gap overlaps the interval
(<code>gap_start &lt; pos2 and gap_end &gt; pos1</code>) without matching the event.
MAPQ&nbsp;&ge;&nbsp;{MIN_MAPQ}; each fragment counted once.</p>
<p class="foot"><strong>What is excluded and why.</strong> Reads with no gap: one lying
inside an exon never had the opportunity to use either junction, so counting it would
measure exon length. Gaps elsewhere in the gene: an intron spliced from every transcript
regardless does not compete. Both exclusions raise the figure, and omitting them is the
usual way this number is understated.</p>
<p class="foot"><strong>Limitations.</strong> Defined only where the event is a size gap.
A chimeric junction has no canonical alternative over the same interval and an insertion
is not a choice between two splice outcomes — there the share is undefined, not zero.
It is a share of transcripts, <em>not</em> a measure of peptide abundance, and it does
not by itself establish that the lesion is genomic: that is the DNA panel above.</p>
</div>"""
    else:
        usage_block = ""


    path.write_text(f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">
<title>Junction validation</title><style>
:root{{color-scheme:light;--s0:#f4f4f1;--s1:#fcfcfb;--s2:#eceae4;--bd:#dcd9d0;
--bs:#c4c0b4;--t1:#0b0b0b;--t2:#52514e;--t3:#76746d;--a:#2a78d6;--j:#eb6834;
--ok:#1baf7a;--warn:#eda100}}
@media(prefers-color-scheme:dark){{:root:not([data-theme=light]){{color-scheme:dark;
--s0:#121211;--s1:#1a1a19;--s2:#232321;--bd:#35342f;--bs:#4a4943;--t1:#fff;
--t2:#c3c2b7;--t3:#94928a;--a:#3987e5;--j:#d95926;--ok:#199e70;--warn:#c98500}}}}
:root[data-theme=dark]{{color-scheme:dark;--s0:#121211;--s1:#1a1a19;--s2:#232321;
--bd:#35342f;--bs:#4a4943;--t1:#fff;--t2:#c3c2b7;--t3:#94928a;--a:#3987e5;
--j:#d95926;--ok:#199e70;--warn:#c98500}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--s0);color:var(--t1);
font:15px/1.6 ui-sans-serif,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}}
.w{{max-width:820px;margin:0 auto;padding:40px 20px 70px}}
@media(max-width:560px){{.w{{padding:24px 16px 48px}}}}
h1{{font-size:24px;margin:0 0 4px;letter-spacing:-.01em}}
h2{{font-size:17px;margin:38px 0 10px;padding-bottom:6px;border-bottom:1px solid var(--bd)}}
p{{margin:0 0 12px;color:var(--t2)}}.meta{{font-size:13px;color:var(--t3);margin:0 0 6px}}
code{{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:.875em;background:var(--s2);
border:1px solid var(--bd);border-radius:4px;padding:1px 5px;color:var(--t1)}}
.card{{background:var(--s1);border:1px solid var(--bd);border-radius:10px;padding:16px 18px;margin:0 0 14px}}
.verdict{{border-left:3px solid var(--ok)}}
.scroll{{overflow-x:auto}}table{{border-collapse:collapse;width:100%;font-size:13.5px;min-width:480px}}
th,td{{padding:7px 10px;text-align:left;border-bottom:1px solid var(--bd)}}
th{{font-size:11px;text-transform:uppercase;letter-spacing:.05em;color:var(--t3);
font-weight:600;border-bottom:1px solid var(--bs);white-space:nowrap}}
td.num{{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}}
tr.hit td{{font-weight:600;color:var(--t1)}}
.barcell{{width:40%;min-width:120px}}
.track{{background:var(--s2);border-radius:3px;height:10px;width:100%;display:block}}
.bar{{height:10px;border-radius:3px;background:var(--a);display:block;min-width:2px}}
.bar.j{{background:var(--j)}}
.gauge{{position:relative;height:30px;background:var(--s2);border-radius:6px;margin:14px 0 6px}}
.fill{{position:absolute;left:0;top:0;bottom:0;border-radius:6px;background:var(--j)}}
.tick{{position:absolute;top:-6px;bottom:-6px;width:2px;background:var(--bs)}}
.ticklab{{position:absolute;top:26px;font-size:11px;color:var(--t3);transform:translateX(-50%)}}
.foot{{font-size:12.5px;color:var(--t3);margin:0 0 10px}}
.stack{{display:flex;gap:2px;height:30px;margin:16px 0 10px}}
.seg{{height:30px;border-radius:4px;display:block;min-width:2px}}
.seg.j{{background:var(--j)}}.seg.a{{background:var(--a)}}
.legend{{font-size:13px}}.key{{display:inline-block;width:10px;height:10px;border-radius:3px;margin-right:6px;vertical-align:baseline}}
.key.j{{background:var(--j)}}.key.a{{background:var(--a)}}
</style></head><body><div class="w">

<h1>Junction validation</h1>
<p class="meta"><code>{args.chrom}:{args.pos1:,}-{args.pos2:,}</code> · {size} bp
{" · " + args.label if args.label else ""} · internal, contains results</p>

<div class="card verdict">
<p><strong>{call}</strong> — {why}.</p>
<p>DNA depth inside the interval is <strong>{cn['inside']:.1f}×</strong> against
<strong>{cn['flank']:.1f}×</strong> in the flanks: a ratio of
<strong>{cn['ratio']:.2f}</strong>.</p>
</div>

{usage_block}

<h2>The discriminator: DNA copy number</h2>
<p>A heterozygous deletion removes one of two copies, so depth inside falls to about
half the flanks. Splicing does not touch the DNA, so the ratio stays near one.
<strong>Expression does not enter this comparison</strong>, which is why it works where
read counts do not.</p>
<div class="gauge">
  <span class="fill" style="width:{ratio_pct:.1f}%"></span>
  <span class="tick" style="left:50%"></span><span class="ticklab" style="left:50%">0.50 heterozygous</span>
  <span class="tick" style="left:100%"></span><span class="ticklab" style="left:100%">1.00 no deletion</span>
</div>
<p class="foot" style="margin-top:26px">Left flank {cn['left']:.1f}× ·
inside {cn['inside']:.1f}× · right flank {cn['right']:.1f}×. Flanks are
{FLANK_BP} bp starting 100 bp outside each breakpoint; the interval is trimmed
{EDGE_TRIM} bp at each edge so partly-overlapping reads do not blur the boundary.</p>

<h2>RNA gap sizes at the locus</h2>
<p>{total:,} reads at MAPQ ≥ {MIN_MAPQ}, {sum(sizes.values()):,} carrying an N gap.
The junction's own population is separate from the gene's ordinary introns — different
size, different start.</p>
<div class="scroll"><table>
<thead><tr><th>gap size</th><th class="num">reads</th><th class="barcell">&nbsp;</th><th>what it is</th></tr></thead>
<tbody>{''.join(rows)}</tbody></table></div>

<h2>Where each gap starts</h2>
<div class="scroll"><table>
<thead><tr><th class="num">position</th><th class="num">reads</th><th class="num">offset</th><th>&nbsp;</th></tr></thead>
<tbody>{''.join(start_rows)}</tbody></table></div>

<div class="card">
<p><strong>A splice-aware tool calling this a novel intron is not evidence against the
lesion.</strong> A transcript from a deleted allele, aligned to a reference that still
carries the deleted bases, looks exactly like a novel intron. The call is the expected
observation either way, so it cannot discriminate — only the DNA can.</p>
<p><strong>Nor is the ratio of RNA reads to DNA fragments a discriminator.</strong> At a
highly expressed gene the RNA depth exceeds the DNA depth by orders of magnitude whether
or not the lesion is real; that ratio measures expression, not artefact.</p>
</div>

</div></body></html>
""")


if __name__ == "__main__":
    main()
