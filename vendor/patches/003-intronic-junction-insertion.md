# Patch 003 — bases inserted at an intronic junction are not coding sequence

**File:** `neosv/fusion_class.py`
**Members:** new `SVFusion.junction_in_cds`; `SVFusion.nt_sequence_ins`
**Status:** applied to the vendored copy; not yet submitted upstream

## What was wrong

SV callers report the bases inserted at a junction as part of the ALT allele
(e.g. an ALT of `AGA[chrN:<pos>[` carries two bases, `GA`, between the breakends). NeoSV
builds the fusion coding sequence as

    5' coding half  +  inserted bases  +  3' coding half

for EVERY junction, whatever region the breakends fall in. When a breakend is
intronic, `truncate_cds` already does the right thing with the halves: it keeps
whole exons up to and from the intron, as splicing would. But the inserted bases
lie in that intron and are spliced out with it. Placing them between the halves
inserts sequence the mature transcript does not contain:

- insertion length not a multiple of 3 -> a frameshift the cell never makes;
- multiple of 3 -> an in-frame insertion of residues that do not exist.

Either way the "mutant" protein differs from wild type, its windows pass the
wild-type subtraction, and they are emitted as neopeptides. For a deletion or
duplication inside a single intron, the correct result — and the result after
this patch — is the wild-type protein, hence no peptides.

## The fix

The inserted bases are used only when BOTH breakends cut inside a coding exon
(`truncate_cds` marks that exon `intact=False`). Otherwise `nt_sequence_ins`
returns an empty string.

## Effect

On a cell-line lineage, re-running the patched generator on the same admitted VCFs
removed the large majority of peptides from deletions and tandem duplications that
lie inside a single intron; almost all of them owed their peptides to a junction
insertion. Inter-gene fusions with intronic breakpoints changed sequence but were
kept, and coding-exon junctions were unchanged. Run-specific counts are recorded
with the results, not here.

## Why RNA evidence does not catch it

The lesion is real and the intron is transcribed in pre-mRNA, so reads spanning
the deletion exist. RNA confirms the lesion, not the protein.

## Not addressed here

A junction with one breakend in a coding exon and the other intronic is still
modelled by joining to the next exon boundary; what the cell actually transcribes
there depends on splice-site usage that this model does not represent.
