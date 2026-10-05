# Patch 005 — choose the transcript a cell most likely makes

**File:** `neosv/transcript_utils.py` (`get_transcript`, new
`PREFERRED_TRANSCRIPTS` / `set_preferred_transcripts`)
**Status:** applied to the vendored copy; not yet submitted upstream

## What was wrong

At each breakend NeoSV used the complete transcript with the largest GENOMIC
SPAN. Span has no biological meaning and favours rare isoforms with distant
exons or N-terminal extensions. The peptides of an SV depend on the isoform: an
exon-skipping isoform can turn an exonic in-frame deletion into an apparent
intronic one, with entirely different junction peptides; and the expression call
of stage 7 uses the chosen isoform's TPM, so an unexpressed isoform can make an
expressed event look silent.

## The fix

Among complete transcripts at the breakend: MANE Select first, then Ensembl
canonical, then the longest as before. The tags are read from the release's GTF
(pyensembl does not load them) by `svneo/generators/_neosv_extensions.py`, cached
beside it, and installed with `set_preferred_transcripts`. Each output row
records the rule used (`transcript_rule1/2`). With no tags available the
behaviour is unchanged and every row says `longest`.

Test: `test_preferred_transcript_order`.
