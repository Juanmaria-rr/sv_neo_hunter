"""
_neosv_extensions.py — what this pipeline needs from a peptide generator and
NeoSV does not provide.

NeoSV emits only post-MHC-filtered neoantigens and carries no identifier linking
a peptide back to the SV that produced it. Three things are therefore added here,
in this repository's own code, on top of the unmodified upstream library:

    1. SV identifiers, recovered from the VCF and matched back by coordinates
    2. the junction position within the fusion protein, so a peptide can be
       tested for whether it actually crosses the breakpoint
    3. a candidate-peptide table written BEFORE any MHC step

None of this changes peptide generation: the sequences come from NeoSV's own
`generate_neoepitopes`, which this module never touches. Verified — NeoSV with
the frame patch and the reference fork produce byte-identical peptide sets
(`tools/compare_generators.py`).

WHY NOT PATCH THE VENDORED LIBRARY FOR THIS
-------------------------------------------
The vendored copy carries exactly one patch, a three-character frame correction
that is required for the tool to be correct at all with current pyensembl. Every
other need is additive, so it lives here: the smaller the diff against upstream,
the cheaper it is to take a new release, and the clearer it is that peptide
sequences are upstream's and not ours.
"""
from __future__ import annotations

import os

import gzip
import re

#: ALT of a mate-carrying breakend, e.g. `G]chr5:123456]`.
_MATE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


# ============================================================================
# 1. SV identifiers, recovered from the VCF
# ============================================================================

def build_sv_id_map(vcf_path: str) -> dict:
    """Map breakend coordinates to the VCF's own ID column.

    NeoSV's SV objects carry no identifier, so the link from a peptide back to
    its SV call has to be re-established. Coordinates are the only thing both
    sides agree on, and they are exact: the tool does not move a breakend.

    Keyed on the *unordered* coordinate pair, because the two breakends of one
    junction may be visited from either side.
    """
    opener = gzip.open if vcf_path.endswith(".gz") else open
    by_position: dict[tuple, str] = {}
    with opener(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, vid, _, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            chrom = chrom.replace("chr", "")
            mate = _MATE.search(alt)
            if not mate:
                continue
            mate_chrom, mate_pos = mate.group(1).replace("chr", ""), mate.group(2)
            key = tuple(sorted([(chrom, int(pos)), (mate_chrom, int(mate_pos))]))
            by_position.setdefault(key, vid)
    return by_position


def sv_id_for(sv, sv_id_map: dict) -> str:
    """Look up one SV's identifier; 'NA' when it cannot be matched."""
    key = tuple(sorted([(str(sv.chrom1).replace("chr", ""), int(sv.pos1)),
                        (str(sv.chrom2).replace("chr", ""), int(sv.pos2))]))
    return sv_id_map.get(key, "NA")


# ============================================================================
# 2. Where the junction falls in the fusion protein
# ============================================================================

def junction_indices(svfusion) -> tuple[int | None, int | None]:
    """Nucleotide and amino-acid offset of the junction within the fusion.

    The 5' side's coding sequence, plus any inserted bases, is everything that
    precedes the junction; its length is the junction offset. Dividing by three
    gives the residue index.
    """
    try:
        if svfusion.cc_1.part == "5":
            left = svfusion.cc_1.nt_sequence + svfusion.nt_sequence_ins
        else:
            left = svfusion.cc_2.nt_sequence + svfusion.nt_sequence_ins
    except AttributeError:
        return None, None
    return len(left), len(left) // 3


def spans_junction(aa_sequence: str, peptide: str, junction_aa: int | None,
                   junction_nt: int | None = None) -> bool | None:
    """Does this peptide actually cross the breakpoint?

    A sliding window over the fusion protein also produces peptides lying wholly
    on one side of the junction. Those are ordinary peptides of the intact
    reading frame, not SV-derived neoantigens, and they can match a reference
    catalogue while carrying no information about the variant.

    Returns None when the junction offset is unknown, so "not determined" is
    never silently reported as "does not span".
    """
    if junction_aa is None or not peptide or not aa_sequence:
        return None
    # Decided in NUCLEOTIDES when the offset is known. The junction lies between
    # nucleotides junction_nt-1 and junction_nt; a window covering nucleotides
    # [3*start, 3*(start+len)) crosses it iff 3*start < junction_nt <
    # 3*(start+len). In residues this is not the same test: when the junction
    # falls inside a codon, that residue is chimeric (new), and a window starting
    # exactly on it was reported as not spanning by `start < junction_aa`.
    start = aa_sequence.find(peptide)
    while start != -1:
        if junction_nt is not None:
            if 3 * start < junction_nt < 3 * (start + len(peptide)):
                return True
        elif start < junction_aa < start + len(peptide):
            return True
        start = aa_sequence.find(peptide, start + 1)
    return False


# ============================================================================
# 3. The candidate-peptide table, written before any MHC step
# ============================================================================

# ============================================================================
# Transcript preference (vendor patch 005)
# ============================================================================

RULE_NAMES = {0: "MANE_Select", 1: "Ensembl_canonical"}


def transcript_ranks(release: int, cache_dir: str | None) -> dict[str, int]:
    """transcript_id -> 0 (MANE Select) / 1 (Ensembl canonical), from the GTF.

    pyensembl does not load GTF tags, so they are read once from the cached GTF
    and kept beside it as a small TSV. Returns {} if the GTF cannot be found, in
    which case the generator falls back to NeoSV's longest-transcript rule and
    every row says so in `transcript_rule*`.
    """
    import glob, gzip, re
    roots = [cache_dir] if cache_dir else []
    roots.append(os.environ.get("PYENSEMBL_CACHE_DIR", ""))
    gtf = next((g for r in roots if r for g in glob.glob(
        os.path.join(r, "**", f"ensembl{release}", f"*.{release}.gtf.gz"), recursive=True)), None)
    if not gtf:
        return {}
    cache = gtf.replace(".gtf.gz", ".transcript_tags.tsv")
    if not os.path.exists(cache):
        tid_re, tag_re = re.compile(r'transcript_id "([^"]+)"'), re.compile(r'tag "([^"]+)"')
        with gzip.open(gtf, "rt") as src, open(cache, "w") as dst:
            dst.write("transcript_id\trank\n")
            for line in src:
                if "\ttranscript\t" not in line:
                    continue
                tags = set(tag_re.findall(line))
                rank = 0 if "MANE_Select" in tags else 1 if "Ensembl_canonical" in tags else None
                if rank is not None:
                    dst.write(f"{tid_re.search(line).group(1)}\t{rank}\n")
    ranks = {}
    with open(cache) as handle:
        next(handle)
        for line in handle:
            tid, rank = line.rstrip("\n").split("\t")
            ranks[tid] = int(rank)
    return ranks


def breakend_region(transcript, pos) -> str:
    """Where a breakend falls in the transcript the generator used."""
    if transcript is None or pos is None:
        return ""
    pos = int(pos)
    if any(s <= pos <= e for s, e in transcript.coding_sequence_position_ranges):
        return "CDS"
    exons = sorted((x.start, x.end) for x in transcript.exons)
    if any(s <= pos <= e for s, e in exons):
        return "UTR"
    if transcript.start <= pos <= transcript.end:
        return "intron"
    return "outside"


def junction_model(fusion, r1: str, r2: str) -> str:
    """How the mutant transcript is built, and how certain that model is.

    The sequence model is exact when the product is a join of coding sequence
    (both breakends in coding exons) or of whole exons (both in introns). When
    one breakend is exonic and the other intronic, the exon on the intronic side
    has lost its splice site and what the cell transcribes depends on cryptic
    splice-site use, which no sequence rule decides: the generator joins to the
    next exon boundary, and this column says that is an assumption.
    """
    g1, g2 = _attr(fusion, "cc_1", "gene_name"), _attr(fusion, "cc_2", "gene_name")
    if g1 != g2:
        return "inter-gene fusion"
    pair = {r1, r2}
    if pair == {"CDS"}:
        return "coding join"
    if pair == {"intron"}:
        return "whole-exon change (or no change, same intron)"
    if pair == {"CDS", "intron"}:
        return "exon-intron: splice outcome uncertain"
    return f"{r1}-{r2}: splice outcome uncertain"


def transcript_rule(transcript_id, ranks: dict) -> str:
    if not transcript_id:
        return ""
    return RULE_NAMES.get(ranks.get(transcript_id), "longest")


#: Column order of the emitted table. Matches the contract in base.py.
PEPTIDE_COLUMNS = ["prefix", "sv_id", "chrom1", "pos1", "gene1", "transcript_id1",
                   "strand1", "chrom2", "pos2", "gene2", "transcript_id2", "strand2",
                   "svtype", "frame_effect", "junction_nt", "junction_aa",
                   "spans_junction", "neopeptide", "pep_length",
                   "transcript_rule1", "transcript_rule2",
                   "breakend_region1", "breakend_region2", "junction_model"]


def write_all_neopeptides(path: str, fusions: list, sv_id_map: dict,
                          prefix: str, ranks: dict | None = None) -> int:
    """Write every candidate peptide, one row per (fusion, peptide).

    This is the table the whole downstream analysis consumes, and it exists
    because a sequence-identity test is HLA-independent: filtering by predicted
    binding here would discard true matches and would require a licensed
    predictor for a question that does not need one.
    """
    written = 0
    with open(path, "w") as handle:
        handle.write("\t".join(PEPTIDE_COLUMNS) + "\n")
        for fusion in fusions:
            sv = fusion.sv
            junction_nt, junction_aa = junction_indices(fusion)
            base = {
                "prefix": prefix,
                "sv_id": sv_id_for(sv, sv_id_map),
                "chrom1": sv.chrom1, "pos1": sv.pos1,
                "chrom2": sv.chrom2, "pos2": sv.pos2,
                "gene1": _attr(fusion, "cc_1", "gene_name"),
                "transcript_id1": _attr(fusion, "cc_1", "transcript_id"),
                # Transcript strand, emitted so results can be stratified by it.
                "strand1": _attr(fusion, "cc_1", "strand"),
                "gene2": _attr(fusion, "cc_2", "gene_name"),
                "transcript_id2": _attr(fusion, "cc_2", "transcript_id"),
                "strand2": _attr(fusion, "cc_2", "strand"),
                "svtype": getattr(sv, "svtype", "") or _svtype_from_pattern(sv),
                # NeoSV names this `frame_effect`, not `frameshift`. Reading
                # the wrong name returned "" for all 72,351 rows without
                # erroring, so the column existed and was always empty — and
                # the one thing it reports is the reliability warning below.
                "frame_effect": _frame_effect(fusion),
                "junction_nt": junction_nt, "junction_aa": junction_aa,
                # Which rule picked each transcript: MANE_Select,
                # Ensembl_canonical, or longest (genomic span) as the fallback.
                "transcript_rule1": transcript_rule(
                    _attr(fusion, "cc_1", "transcript_id"), ranks or {}),
                "transcript_rule2": transcript_rule(
                    _attr(fusion, "cc_2", "transcript_id"), ranks or {}),
            }
            r1 = breakend_region(getattr(fusion.cc_1, "transcript", None), sv.pos1)
            r2 = breakend_region(getattr(fusion.cc_2, "transcript", None), sv.pos2)
            base.update(breakend_region1=r1, breakend_region2=r2,
                        junction_model=junction_model(fusion, r1, r2))
            for peptide in fusion.neoepitopes:
                row = dict(base)
                row["neopeptide"] = peptide
                row["pep_length"] = len(peptide)
                row["spans_junction"] = spans_junction(
                    fusion.aa_sequence, peptide, junction_aa, junction_nt)
                handle.write("\t".join(
                    "" if row.get(c) is None else str(row.get(c, ""))
                    for c in PEPTIDE_COLUMNS) + "\n")
                written += 1
    return written


def write_annotation(path: str, fusions: list, sv_id_map: dict, prefix: str) -> int:
    """One row per annotated SV, keyed on the same sv_id as the peptide table."""
    columns = ["prefix", "sv_id", "chrom1", "pos1", "gene1", "transcript_id1",
               "strand1", "chrom2", "pos2", "gene2", "transcript_id2", "strand2",
               "svtype"]
    with open(path, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for fusion in fusions:
            sv = fusion.sv
            row = {"prefix": prefix, "sv_id": sv_id_for(sv, sv_id_map),
                   "chrom1": sv.chrom1, "pos1": sv.pos1,
                   "chrom2": sv.chrom2, "pos2": sv.pos2,
                   "gene1": _attr(fusion, "cc_1", "gene_name"),
                   "transcript_id1": _attr(fusion, "cc_1", "transcript_id"),
                   "strand1": _attr(fusion, "cc_1", "strand"),
                   "gene2": _attr(fusion, "cc_2", "gene_name"),
                   "transcript_id2": _attr(fusion, "cc_2", "transcript_id"),
                   "strand2": _attr(fusion, "cc_2", "strand"),
                   "svtype": getattr(sv, "svtype", "") or _svtype_from_pattern(sv)}
            handle.write("\t".join(str(row.get(c, "") or "") for c in columns) + "\n")
    return len(fusions)


# -- small helpers ---------------------------------------------------------

def _frame_effect(fusion) -> str:
    """NeoSV's own verdict on what the fusion does to the reading frame.

    `In-frame`, `Stop-gain`, `Stop-loss`, or `Start-loss`. The last is a
    RELIABILITY WARNING, not a biological category: upstream's own comment says
    that when the 5' CDS is empty or under three residues the start codon is
    lost, the tool falls back to the next ATG it finds, and "such prediction is
    of low reliability, we should annotate it and remove these fusions when
    necessary."

    This pipeline read the attribute as `frameshift`, which does not exist, so
    `getattr` returned the default for every row: the column shipped empty and
    the warning never reached anyone. Peptides that do not span the junction are
    only credible as neoantigens if the frame downstream is genuinely altered,
    and this is the field that says so.
    """
    try:
        return str(fusion.frame_effect)
    except Exception:                       # noqa: BLE001 - upstream may raise
        return ""


def _attr(fusion, side: str, name: str):
    """Read an attribute of one coding-collection side, tolerating absence."""
    collection = getattr(fusion, side, None)
    return getattr(collection, name, None) if collection is not None else None


def _svtype_from_pattern(sv) -> str:
    """NeoSV encodes breakend orientation as patterns 1-4; map to an SV class.

    Intra-chromosomal only — an inter-chromosomal junction is a translocation
    regardless of orientation.
    """
    if str(sv.chrom1) != str(sv.chrom2):
        return "BND"
    return {1: "DEL", 2: "DUP", 3: "INV", 4: "INV"}.get(getattr(sv, "pattern", None), "")
