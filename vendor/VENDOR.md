# Vendored dependency: NeoSV

`vendor/neosv/` is a copy of the published NeoSV package, included so this
repository runs end to end without an external checkout or network access.

> Shi, Y., Jing, B. & Xi, R. Comprehensive analysis of neoantigens derived from
> structural variation across whole genomes from 2528 tumors.
> *Genome Biology* **24**, 169 (2023).
> https://github.com/ysbioinfo/NeoSV · `pip install neosv`

| | |
|---|---|
| Licence | **MIT**, Copyright (c) 2022 Yang Shi — see [`LICENSE.NeoSV`](LICENSE.NeoSV) |
| Source | `NeoSV-main/neosv`, line endings normalised to LF |
| Local modifications | **two patches** — see [`patches/`](patches/) |

MIT permits modification and redistribution provided the copyright notice is
retained. `LICENSE.NeoSV` is included for that purpose and must accompany any
distribution of this repository.

## The patches

[`001-pyensembl-stop-codon-frame`](patches/001-pyensembl-stop-codon-frame.md) —
three characters, and required for correctness. pyensembl changed its
coding-sequence accessor at v2.3.13; NeoSV's compensating `- 3` then shifts the
3' side of a fusion one amino acid out of frame. Unpatched, with pyensembl
2.10.1, the tool emits 170 peptides where the patched version emits 64, **with
none in common** — the frame-shifted sequences pass the wild-type subtraction as
spurious neopeptides.

The issue was identified and first corrected in
[NeoSV-Trace](https://github.com/winterga/NeoSV-Trace) by Greyson Wintergerst.
This repository applies the same correction to the MIT-licensed original rather
than depending on the fork.

## What this repository adds, and where

Peptide sequences come from NeoSV, untouched. Three things it does not provide
are implemented in `src/svneo/generators/_neosv_extensions.py` — this
repository's own code, not the tool's:

| Need | Why NeoSV cannot supply it |
|---|---|
| **SV identifiers** on every peptide row | NeoSV carries no `sv_id`; without it peptides cannot be collapsed to genomic events, which is the reported unit |
| **Candidate peptides before MHC filtering** | NeoSV writes only post-netMHCpan binders; a sequence-identity test is HLA-independent and must not require a licensed predictor |
| **`spans_junction`** per peptide | the sliding window also emits peptides lying wholly on one side of the breakpoint, which are not SV-derived at all |

`sv_id` is recovered by matching breakend coordinates back to the VCF rather
than threaded through the library, so the diff against upstream stays at the
single frame patch.

Equivalence is verified rather than assumed: the patched library plus these
extensions reproduces NeoSV-Trace's output exactly — 64 of 64 peptides shared,
identical `(peptide, sv_id)` keys, `spans_junction` agreeing on every row.
Reproduce with `tools/compare_generators.py`.

## Patch policy

Keep the diff against upstream minimal: additive needs belong in
`_neosv_extensions.py`, not in `vendor/neosv/`. Before adding a patch, ask
whether the same result can be obtained outside the vendored tree.

If a patch is unavoidable:

1. Apply it to `vendor/neosv/`.
2. Add `patches/NNN-short-name.patch` and a matching `.md` explaining the bug,
   the measured effect, and how to revert.
3. Record it in the table below.
4. Re-run `tools/compare_generators.py` against the previous behaviour and state
   whether peptides changed.

Never patch silently: a peptide-level difference invalidates comparisons against
catalogues built elsewhere, which is precisely what patch 001 demonstrates.

| Patch | File | Effect on peptides | Upstream? |
|---|---|---|---|
| 001-pyensembl-stop-codon-frame | `fusion_class.py` | corrects them; unpatched output is spurious with pyensembl > 2.3.13 | not yet submitted |
| 002-minus-strand-cds-order | `transcript_utils.py` | corrects them; unpatched, minus-strand transcripts return the wrong 5' CDS in every region — 0 of 188 correct | not yet submitted |

## Known upstream behaviour worked around in the backend

Not bugs, but constraints that shape `src/svneo/generators/neosv.py`:

1. **Only `.vcf` / `.bedpe` are accepted**, checked by file extension — a
   `.vcf.gz` is rejected outright, so the backend decompresses to a temporary
   file.
2. **No admission logic of its own.** Given a raw VCF it loads panel-of-normals
   and copy-number-inferred records too, so it is handed the stage-1 admitted
   call set, never the raw file.
3. **Output is post-MHC only**, hence the candidate table written by this
   repository.

## Related: the NeoSV-Trace fork

The fork is not vendored here. It produces peptides identical to the patched
NeoSV, so it adds nothing this pipeline needs, and its additions carry no
declared licence. The `neosv_trace` backend remains available for reproducing
analyses run with it, pointing at an external checkout:

```yaml
peptide_generator: neosv_trace
resources:
  neosv_path: /path/to/NeoSV-Trace     # directory containing neosv_trace/
```
