# Patch 006 — minus-strand 3' side cut inside a coding exon

**File:** `neosv/fusion_utils.py` (`truncate_cds`, minus strand, direction 3)
**Status:** applied to the vendored copy; not yet submitted upstream

## What was wrong

For the 3' side of a fusion whose breakpoint falls inside a coding exon of a
minus-strand transcript, the retained part of that exon is [exon start,
breakpoint] — the transcript reads from high to low coordinates. NeoSV set
`start = pos`, keeping [breakpoint, exon end] instead. The 3' length was
miscounted, so the 3' sequence was taken from the wrong offset and the frame
shifted. Plus-strand events and intronic breakpoints were unaffected.

Found by `tools/audit_peptide_generation.py`, an independent reconstruction of
same-gene deletions and duplications: every plus-strand exonic deletion agreed,
11 of 12 minus-strand ones did not (e.g. a deletion removing an intron exactly at
both exon edges, whose correct product is the wild-type protein, produced dozens
of peptides). After the patch all modelled events agree.

## The fix

`end = pos` for the minus-strand 3' side. Test: `test_minus_strand_3prime_cut`.
