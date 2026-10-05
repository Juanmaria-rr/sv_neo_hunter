# Patch 004 — a non-AUG start codon still starts with methionine

**Files:** `neosv/sequence_utils.py` (`set_aa_seq`), `neosv/fusion_class.py`
(new `SVFusion.starts_at_native_start`, `frame_effect`)
**Status:** applied to the vendored copy; not yet submitted upstream

## What was wrong

Some Ensembl transcripts are annotated with a non-AUG initiation codon (CUG,
GUG, ...), typically N-terminally extended isoforms. Translation still begins
with methionine and Ensembl's protein starts with M, but NeoSV translated the
first codon literally (CUG -> L). The fusion protein then differed from wild type
at residue 1, so the four N-terminal windows (8-11 residues) passed the wild-type
subtraction as spurious neopeptides, for ANY SV hitting such a transcript. The
same literal test (`startswith('ATG')`) labelled those events `Start-loss`.

## The fix

When the fusion keeps the 5' transcript's own start codon, its first residue is
methionine; `Start-loss` now means the 5' side lost its own start codon,
compared with the transcript's annotated CDS rather than the string `ATG`.

Test: `test_non_aug_start_is_methionine`.
