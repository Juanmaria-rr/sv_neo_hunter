# Patch 002 — coding exons must reach consumers in reading order

**File:** `neosv/transcript_utils.py`
**Functions:** `get_cds_range`, `get_noncds_range`
**Status:** applied to the vendored copy; not yet submitted upstream
**Reproduce:** `python tools/diagnose_minus_strand_cds.py` — exits 1 unpatched,
0 patched

## What was wrong

`truncate_cds(transcript, '5', pos)` returns what a transcript contributes as the
**head** of a fusion protein: its coding sequence from the start codon to the
breakpoint. On minus-strand transcripts it returned the wrong length in **every**
region of the transcript.

Swept across every coding exon and intron of nine transcripts, comparing the
length returned against the length that should be returned:

| | regions | correct before | correct after |
|---|---|---|---|
| EGFR, PTEN, PIK3CA (plus) | 87 | 87 | 87 |
| ITGA11, TP53, BRCA1, GOLGA3, KRAS, BRAF (minus) | 188 | **0** | **188** |

Two distinct failure modes on the minus strand:

| Breakpoint in | Result before the patch |
|---|---|
| an **intron** | the head came back **empty** — no start codon, so the fusion was flagged `Start-loss`, the junction sat at residue 0, and no peptide could span it |
| a **coding exon** | the length was counted from the wrong end — too short near the start of the transcript, **too long** near its end |

The second mode is the more damaging: sequence the gene does not contribute
entered the fusion, and the sliding window turned it into peptides. On a worked
three-exon example a breakpoint inside the last coding exon returned 278 nt where
26 were correct — roughly 84 residues the gene never contributed. The defect
could **fabricate** candidates, not merely lose them.

## Why it happened

`transcript.coding_sequence_position_ranges` is ordered by **coordinate**.
Reading order is the same thing on the plus strand and the reverse of it on the
minus strand.

`get_cds_range` documents its output as "from 5' to 3', cds 1, cds 2, ...", and
both consumers rely on that:

- `get_noncds_range` builds the introns by subtracting neighbouring exons. Its
  minus-strand branch compensates for a reversal that never happened, so the
  intervals come out with `start > end` — no position can ever fall inside one —
  and the first interval, meant to be the 5'UTR, is measured from the end of the
  first exon *in the list* rather than the first in reading order, collapsing
  across most of the gene.
- `truncate_cds` slices the list to count the exons preceding a breakpoint. With
  coordinate order on a minus-strand transcript it counts from the wrong end.

The invariant was documented and never enforced. A breakpoint in an intron then
matched only the oversized first interval, whose index means "no exons precede
this", so the head was empty.

### Worked example

Three coding exons — A 100–200, B 300–450, C 550–600 — in a transcript spanning
50–650. On the minus strand it is read C, B, A.

| Breakpoint | Correct head | Returned before |
|---|---|---|
| inside A (150) | 253 nt | 51 nt |
| intron A–B (250) | 202 nt | **0 nt** |
| inside B (375) | 127 nt | 177 nt |
| intron B–C (500) | 51 nt | **0 nt** |
| inside C (575) | 26 nt | **278 nt** |

`python tools/diagnose_minus_strand_cds.py --explain` prints this walkthrough
with every interval and the formula that produced it.

## The change

Sort the ranges into reading order in `get_cds_range` — descending coordinate for
minus-strand transcripts — and route `get_noncds_range` through it instead of
reading the raw accessor. This restores the invariant the docstring already
claims, for every consumer at once.

Deliberately *not* done: patching `truncate_cds`. Its minus-strand branches are
already written for reading order and become correct as soon as they receive it.
Correcting them individually would add a second exception on top of the first,
and leave the documented invariant still false.

## Verification

- `tools/diagnose_minus_strand_cds.py` — 275 of 275 regions correct, exit 0
- **Zero regression on the plus strand**: 90 probe positions, none changed. The
  necessary condition — repair the minus strand without disturbing the half that
  already worked.
- 189 of 194 minus-strand positions changed, as expected. The five unchanged are
  positions where the wrong answer happened to coincide with the right one.
- `tests/test_svneo.py::test_minus_strand_cds_is_read_in_transcript_order`

## Consequence for results produced before this patch

Every figure generated before 2026-09-08 was produced with the defect present.
Minus-strand genes were systematically degraded, and the effects are documented
in `docs/OPEN_QUESTIONS.md`: 84.0% of minus-strand fusions flagged `Start-loss`
against 13.8% of plus-strand ones, a median `junction_aa` of 3 among
catalogue-matching peptides, and 1 of 408 matching peptides spanning its junction
where 67 were expected.

Re-running is a decision about published numbers rather than a refactor, so the
old results are retained until that decision is taken. A reference catalogue
built with the same tool family would carry the same artefact; correcting only
one side may *reduce* apparent overlap, which would be informative rather than a
regression.
