# Method — MHC binding propagated to peptides and events

Accompanies the outputs of `tools/propagate_binding_to_events.py`:
`peptide_binding.tsv`, `events_binding.tsv`, `funnel_mhc.tsv` and
`mhc_propagation.json`. Read this before quoting any number from them.

Binding prediction is **not** part of the recurrence test, which is a sequence
identity test and HLA-independent by construction. It is a separate, optional
layer, needed only when these counts are to be compared with an analysis that
had already filtered for presentation.

## 1. What is counted

Three columns run down every table. They are **not** three estimates of one
quantity; they answer three different questions, and only one of them is
comparable with a presentation-filtered analysis.

| column | question |
|---|---|
| `n_all` | recurrence, no MHC filter at all |
| `n_mhc_panel` | could this be presented by **someone** in the cohort? |
| `n_mhc_autologous` | was it presentable in the **patient where it was observed**? |

- **`n_mhc_panel` is an upper bound, not a result.** It takes the union of every
  class I allele observed anywhere in the cohort. A peptide binding one rare
  allele carried by one patient counts, even if the patient who actually carried
  the peptide could never have presented it. On a real cohort the two branches
  differ severalfold, so quoting this one against an earlier presentation-filtered
  count inflates it.
- **`n_mhc_autologous` is the comparable quantity.** A peptide counts only if at
  least one HLA-typed patient both carries it and has, in their own genotype, an
  allele it binds — roughly six alleles rather than a few hundred. This is what a
  per-sample presentation filter measures.

Two peptide levels are reported and they are easy to confuse:

- **matched unique peptides** — distinct peptides identical to a reference
  peptide. This is the level that stands opposite a published recurrence count.
- **peptides carried by credible events** — a much smaller subset, downstream of
  gene concordance and the complexity / self cuts. Correct input for the event
  funnel, wrong input for a comparison.

## 2. How it is calculated

**Inputs.** The peptide x allele prediction table from `tools/run_netmhcpan.py`;
one peptide file per catalogue patient; one class I genotype per typed patient;
`credible_events.tsv` and `stage8_rna_evidence.tsv` from the same run.

**Binder definition.** Imported from `tools/run_netmhcpan.py`, not restated:
`IC50 <= 500 nM` **and** `%Rank_BA <= 2` **and** `%Rank_EL <= 2`, all three
simultaneously. This is stricter than netMHCpan's own convention, where
`%Rank_EL <= 2` alone is a weak binder and `<= 0.5` a strong one. It recognises
one binder class, not two. `--criterion` selects a different rule; see §4.

**Allele spelling.** Three sources, three spellings, and no error if they
disagree: an allele panel written in netMHCpan's *input* form carries no asterisk
(`HLA-A01:01`, because netMHCpan silently returns zero rows if given one), while
its *output* table and typical genotype files carry it (`HLA-A*01:01`). All are
reduced to one form. The script aborts if the prediction table and the panel
share no allele at all, because the natural failure here is a clean page of
zeros rather than a crash.

**Autologous test.** One pass over the typed patients. For each, the genotype is
intersected with the already-computed binding-allele set of every peptide that
patient carries; a non-empty intersection marks the peptide a binder and records
which allele presented it. Peptide-to-patient links are never materialised in
full, so cost is O(peptides).

**Event propagation.** An event survives a branch if at least one of its peptides
does. Downstream funnel rows reuse the definitions in `src/svneo/run.py`
verbatim: `events_rna_supported` admits `rna_tier` in {STRONG, SUGGESTIVE, WEAK},
whereas `events_hc_and_rna` and `events_private_hc_and_rna` admit only
{STRONG, SUGGESTIVE}.

## 3. Limitations

- **`unevaluable` is not `non_binder`.** A peptide contributed only by patients
  with no linkable genotype has no autologous verdict. Such peptides are excluded
  from `n_mhc_autologous`, which therefore **understates** the true count; a
  further class is a patient carrying an allele outside the predicted panel,
  flagged separately in `autologous_verdict`. Read `n_mhc_autologous` as a floor
  whose gap is quantified in `n_unevaluable`, never as an exact value.
- **Coverage is not random.** Typed patients are a subset, and nothing here
  establishes that they are representative, so the unevaluable share is a size,
  not a correction factor.
- **The binder rule is brittle by construction.** Three hard cuts required
  simultaneously means a marginal move in any one flips the verdict. Per-cut
  agreement between two netMHCpan versions stays above 99%, yet binder agreement
  is far lower. Expect the same fragility against any re-run, any predictor and
  any threshold nudge; the rule is kept for comparability, not because it is
  well-conditioned.
- **Prediction, not observation.** netMHCpan output is not evidence of
  presentation. Neither branch says anything about proteasomal processing, TAP
  transport, surface abundance or T-cell recognition.
- **`n_mhc_panel` must never be quoted as a recurrence result.** It exists to
  bound the autologous number from above, and to show how much of a difference is
  the choice of allele set rather than the biology.

## 3b. Margin: how far inside the cut a call sits

The binder rule is three hard cuts required at once, so "passes" is not one
thing: a call at 498.8 nM and a call at 2.5 nM are not the same event, and
treating them alike is how a list ends up resting on coin flips.

`autologous_margin` records, for each of the three cuts, the distance to it as a
fraction of the cut; takes the **smallest** (the constraint actually holding the
call up); then the **largest** across the peptide's presenting alleles, since a
peptide is as robust as its most comfortable allele. 0 is exactly on a boundary.
`autologous_limiting_cut` names the binding constraint, and
`autologous_robustness` bins the margin into `flippable` (<= 10% of a cut),
`marginal` (<= 25%), `solid` (<= 50%) and `robust` (> 50%).

**What this is for.** Nearly all of the calls that disagree between netMHCpan
versions sit close to a threshold, so the margin identifies the population that a
re-run, a version bump or a threshold nudge would move. The `robust` subset is
**not a stricter biological claim** — it is the same criterion minus the calls
that are not reproducible. A peptide-by-peptide downstream analysis built on the
whole list rests partly on unstable calls; built on the robust subset it is
stable, at the cost of some N.

Among marginal calls the limiting cut is spread across all three rather than
concentrated in one, so the fragility is the conjunction, not a single badly
placed threshold.

Propagated to events as `binds_autologous_robust` and
`n_peptides_binding_autologous_robust`.

## 3c. Alternative binder criteria

`--criterion` selects the binder rule. **`three_way` is the default and the only
one comparable with a conventional presentation filter.** The `neo_*` rules apply
netMHCpan 4.2's CEDAR neoepitope head, for which no cut is established; their
thresholds are netMHCpan's generic strong/weak rank conventions borrowed for it,
which is a choice and is labelled as one. `neo_strong` (`%Rank_Neo <= 0.5`) is of
comparable stringency to the three-way rule; `neo_weak` (`<= 2`) is far more
permissive.

Two criteria applied to **identical predictions** agree less well than one
criterion applied across two netMHCpan versions. When reporting a binder count,
which rule produced it matters more than which version did.
`autologous_robustness` predicts stability here too: robust calls survive a
change of predictor head far more often than flippable ones, which is the
argument for the robust subset — it is the part of the list that does not depend
on the arbitrary parts of the criterion.

Alternative criteria write suffixed output files and never overwrite the default.

## 4. netMHCpan 4.1 vs 4.2 — what actually differs

Recorded because the first instinct on seeing low binder agreement between
versions is to suspect the criterion changed. It did not. Sources: the 4.2
release notes, `netMHCpan -h`, its `data/` directory, and the parameter block
netMHCpan 4.1 writes into the header of its own output files.

### The criteria are unchanged

| | 4.1b | 4.2e |
|---|---|---|
| strong-binder rank threshold | `-rth` **0.5** | `-rankS` **0.5** |
| weak-binder rank threshold | `-rlt` **2.0** | `-rankW` **2.0** |
| `-BA` (affinity prediction) | off by default | off by default |

The flags were renamed, the values were not. In any case this repository does not
use them: the cuts are applied downstream to the numeric columns, so the binder
definition is its own and is version-independent.

### Interface differences

| | 4.1b | 4.2e |
|---|---|---|
| default `-l` | `8,9,10,11` | **`9`** (FASTA input only) |
| input types | 0 FASTA, 1 PEPTIDE | 0 FASTA, 1 PEPTIDE, 2 PEPTIDECONT, 3 PEPTIDEMHC, 4 PEPTIDECONTMHC |
| network file | one `synlist.bin` | `synlist_nocontext.bin`, `_context`, `_iedb`, `_cedar` |
| extra modes | — | `-pathogen`, `-neo`, `-context` |

The `-l` change is a latent trap rather than a live one: with peptide input
(`-inptype 1` / `-p`) lengths come from the input file. Switching to FASTA input
under 4.2 would silently yield 9-mers only.

### Model differences — the actual cause of the disagreement

4.1 (Reynisson et al., *NAR* 2020) trains on binding-affinity and eluted-ligand
data with concurrent motif deconvolution. 4.2 (Nilsson et al., *Front Immunol*
2025) adds, on top of that training set, fine-tuning on tens of thousands of
experimentally verified epitopes from IEDB and CEDAR, plus structural features.

The mechanism that moves `%Rank` is worth being explicit about: a rank is not a
property of the peptide, it is a percentile against a **per-allele,
per-scoring-mode background distribution**. 4.2's `data/threshold/` holds a
separate background table for each allele pseudo-sequence in each of six modes
(`ba`, `bacont`, `el`, `elcont`, `eliedb`, `elcedar`); 4.1 shipped one network
file and its corresponding backgrounds. So between versions both the score *and*
the background it is ranked against change. That is why `%Rank` medians shift
several points while per-cut agreement stays above 99%: the reordering is small,
and a hard three-way conjunction amplifies it.

### The mode this analysis does not use

4.2 offers `-neo`, a head fine-tuned on neoepitopes (CEDAR). The default
predictions here are made with plain `-p -a -BA`, i.e. the standard EL+BA method,
with **no** `-neo`, `-pathogen` or `-context`. That is the right choice for
resembling 4.1, but it is a choice: a neoantigen analysis declining a
neoantigen-tuned model in favour of comparability with an older version.

`-neo` **appends** `Score_Neo` and `%Rank_Neo` and moves no other column —
verified by running the bundled test set with and without the flag and confirming
the EL, BA and affinity columns come back identical. `tools/run_netmhcpan.py --neo`
emits them; they are deliberately **not** folded into the binder rule, because no
rank cut for that head is established and inventing one would silently change
what "binder" means between runs.
