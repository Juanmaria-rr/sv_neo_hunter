# Patch 001 — pyensembl stop-codon offset shifts the 3' fusion frame

**File:** `vendor/neosv/fusion_class.py` · **Patch:** `001-pyensembl-stop-codon-frame.patch`
**Status:** required for correctness with pyensembl > 2.3.13 · **Sent upstream:** not yet

## The bug

NeoSV compensates for pyensembl's older behaviour, in which
`coding_sequence_position_ranges` excluded the stop codon while `coding_sequence`
included it, by subtracting 3 when computing the start of the 3' side of a
fusion:

```python
start = end - self.cut_length + 1 - 3
```

pyensembl changed that accessor at **v2.3.13** ([openvax/pyensembl#176](https://github.com/openvax/pyensembl/issues/176)).
With any later release the subtraction is no longer a correction but an error: it
moves the 3' coding sequence three nucleotides — one amino acid — out of frame.

The patch removes the subtraction:

```python
start = end - self.cut_length + 1
```

## Why it matters, measured

Frame-shifted fusion protein sequences do not match the wild-type protein, so
they survive NeoSV's `set(mut) − set(wt)` filter and are emitted as neopeptides.
The result is not a handful of wrong peptides but a wholly different, larger, and
spurious set. On 51 SVs from one sample, with pyensembl 2.10.1:

| | Unpatched | Patched |
|---|---|---|
| Non-empty fusions | 15 | 15 |
| Fusions with identical AA sequence to the patched run | 1 of 15 | — |
| Unique peptides | **170** | **64** |
| Peptides shared with the patched run | **0** | — |

2.7× more peptides, none of them in common. Any analysis run unpatched on a
modern pyensembl is comparing against sequences that do not exist.

## Verification

The patched library reproduces, byte for byte, the peptide set of NeoSV-Trace
(github.com/winterga/NeoSV-Trace), an independent fork that carries the same
correction:

```
rows 64 vs 64 · peptides 64 vs 64, 64 shared, 0 either-only
(peptide, sv_id) keys identical · spans_junction agrees on 64/64
```

Reproduce with `tools/compare_generators.py`.

## Attribution

The issue was identified and first corrected in NeoSV-Trace by Greyson
Wintergerst; the comment there records the diagnosis and the pyensembl issue
number. This repository applies the same one-line correction to the published
MIT-licensed NeoSV rather than depending on the fork.

## Applying and reverting

```bash
# revert to upstream behaviour (only correct with pyensembl <= 2.3.13)
cd vendor/neosv && patch -R -p0 < ../patches/001-pyensembl-stop-codon-frame.patch
```

If a future NeoSV release fixes this upstream, drop the patch and re-run
`tools/compare_generators.py` to confirm the peptide set is unchanged before
using it.
