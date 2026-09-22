"""
rna.py — stages 7-8: expression context, and junction evidence from reads.

This module is the one that two independent workflows disagreed on, so its rules
are stated as rules rather than implied by the code.

    Three events called STRONG by one workflow, over the same BAM:
        event A  BND       very uneven coverage, crossing reads=0  -> NONE
        event B  DEL 2 bp  even coverage,        crossing reads=0  -> test invalid
        event C  DEL 220bp even coverage,        crossing reads>0  -> STRONG
    Three STRONG, zero survivors: two had nothing crossing, and the third is a
    common germline polymorphism by panel-of-normals count, i.e. real but not
    private.

RULE 1  ONLY JUNCTION-CROSSING READS TIER AN EVENT.
        Coverage is context. A breakpoint inside a highly expressed gene has
        thousands of reads whether or not the junction exists.
RULE 2  NEVER max() COVERAGE ACROSS BREAKENDS.
        Reporting max(bp1, bp2) turns a silent locus into a well-covered one
        when its partner is transcribed. Both ends are reported and
        `min` is the summary: an event is in transcribed territory only if BOTH
        ends are.
RULE 3  NO VALID TEST => `UNTESTABLE`, NOT `NONE`.
        `NONE` means tested and negative. A 2 bp deletion cannot produce a CIGAR
        `N` gap at all, so calling it `NONE` fabricates a negative result.
"""
from __future__ import annotations

import os

import pandas as pd

from . import criteria


# ============================================================================
# Stage 7 — expression context
# ============================================================================

def expression(events: pd.DataFrame, anno: pd.DataFrame, isofox_dir: str,
               prefix: str) -> tuple[pd.DataFrame, dict]:
    """Gene TPM, disrupted-isoform TPM, and the background pass-rate gradient.

    The background is not decoration: without it, "18 of 39 events have TPM > 0"
    is uninterpretable, because in a typical transcriptome ~55% of ALL genes
    clear that bar. The gradient says how much information the threshold carries.
    """
    gene_path = os.path.join(isofox_dir, f"{prefix}.gene_data.csv")
    genes = pd.read_csv(gene_path)
    gene_col = next(c for c in genes.columns if c.lower() in ("genename", "gene_name"))
    tpm_col = next(c for c in genes.columns if "tpm" in c.lower())
    gene_tpm = dict(zip(genes[gene_col], genes[tpm_col]))

    all_tpm = genes[tpm_col].dropna()
    background = {f"TPM>{t}": round(100.0 * (all_tpm > t).mean(), 1)
                  for t in criteria.TPM_GRADIENT}
    background["n_genes_quantified"] = int(len(all_tpm))

    isoform_tpm = {}
    tx_path = os.path.join(isofox_dir, f"{prefix}.transcript_data.csv")
    if os.path.exists(tx_path):
        tx = pd.read_csv(tx_path)
        name_col = next((c for c in tx.columns if "transname" in c.lower()), None)
        tx_tpm = next((c for c in tx.columns if "tpm" in c.lower()), None)
        if name_col and tx_tpm:
            isoform_tpm = dict(zip(tx[name_col], tx[tx_tpm]))

    # Join on sv_id explicitly. Using the first column positionally silently
    # matched against `prefix`, so no annotation was ever attached and every
    # event came back with no gene, no transcript and no TPM.
    key = next((c for c in ("sv_id", "sample_sv_id") if c in anno.columns), None)

    rows = []
    for _, event in events.iterrows():
        sv_id = str(event.get("sample_sv_id"))
        match = anno[anno[key].astype(str) == sv_id] if (key and len(anno)) else anno.iloc[0:0]
        record = {"sample_sv_id": sv_id, "gene": event.get("gene")}

        # BOTH breakends' transcripts, not just the first: an inter-chromosomal
        # event has a partner gene whose expression matters just as much.
        tpms = []
        if len(match):
            r = match.iloc[0]
            for side in ("1", "2") if criteria.EXPRESSION_USES_BOTH_BREAKENDS else ("1",):
                gene = r.get(f"gene{side}")
                transcript = r.get(f"transcript_id{side}")
                record[f"gene{side}"] = gene
                record[f"gene{side}_TPM"] = gene_tpm.get(gene)
                record[f"isoform{side}_TPM"] = isoform_tpm.get(transcript)
                value = record[f"isoform{side}_TPM"]
                tpms.append(value if value is not None else record[f"gene{side}_TPM"])

        usable = [t for t in tpms if t is not None]
        # The conservative summary: an event is "expressed" only if the weaker
        # side clears the threshold, matching the min-coverage logic of stage 8.
        record["min_side_TPM"] = min(usable) if usable else None
        record["expressed"] = bool(usable and min(usable) > criteria.TPM_EXPRESSED)
        rows.append(record)
    return pd.DataFrame(rows), background


# ============================================================================
# Stage 7b — independent RNA evidence from the quantifier's own output
# ============================================================================

def isofox_context(events: pd.DataFrame, anno: pd.DataFrame, isofox_dir: str,
                   prefix: str) -> pd.DataFrame:
    """Evidence Isofox already produced, which our read counting does not see.

    Three things, all independent of our own junction test:

    `alt_splice_junc`  novel splice junctions with their fragment support. A
                       breakpoint that coincides with one is transcribed through
                       a rearranged structure, which a CIGAR gap at a fixed
                       offset can miss.
    `pass_fusions`     Isofox's own fusion calls. This is a fusion CALLER's
                       verdict, not our hand-counted supplementary alignments —
                       stronger evidence for inter-chromosomal junctions than
                       anything stage 8 computes, and it was going unused.
    `retained_intron`  a retained intron at the locus, which changes what the
                       transcript actually contains.

    Distances are reported, never thresholded here; `RNA_JUNCTION_MAX_DIST`
    supplies the flag but the base-pair distance travels with it so a reader can
    apply their own cutoff.
    """
    def load(suffix):
        path = os.path.join(isofox_dir, f"{prefix}.{suffix}.csv")
        return pd.read_csv(path, low_memory=False) if os.path.exists(path) \
            else pd.DataFrame()

    sj, fusions, introns = load("alt_splice_junc"), load("pass_fusions"), \
        load("retained_intron")
    key = next((c for c in ("sv_id", "sample_sv_id") if c in anno.columns), None)

    def nearest_sj(chrom, pos):
        """Closest novel splice junction, with its fragment support."""
        if sj.empty or pos is None:
            return None, None, None
        chrom = str(chrom).replace("chr", "")
        sub = sj[sj["Chromosome"].astype(str).str.replace("chr", "", regex=False) == chrom]
        if sub.empty:
            return None, None, None
        distances = pd.concat([(sub["SjStart"] - pos).abs(),
                               (sub["SjEnd"] - pos).abs()], axis=1).min(axis=1)
        best = distances.idxmin()
        return int(distances[best]), int(sub.loc[best, "FragCount"]), \
            str(sub.loc[best, "Type"])

    def fusion_at(chrom1, pos1, chrom2, pos2, tolerance=1000):
        """Does Isofox call a fusion joining these two breakends?"""
        if fusions.empty or pos1 is None or pos2 is None:
            return None, None
        c1 = str(chrom1).replace("chr", "")
        c2 = str(chrom2).replace("chr", "")
        up = fusions["ChrUp"].astype(str).str.replace("chr", "", regex=False)
        down = fusions["ChrDown"].astype(str).str.replace("chr", "", regex=False)
        same = (((up == c1) & (down == c2)) | ((up == c2) & (down == c1)))
        near = ((fusions["PosUp"] - pos1).abs() <= tolerance) | \
               ((fusions["PosUp"] - pos2).abs() <= tolerance)
        hit = fusions[same & near]
        if hit.empty:
            return None, None
        row = hit.iloc[0]
        support = next((row[c] for c in ("TotalFragments", "SplitFrags", "CoverageUp")
                        if c in hit.columns), None)
        return str(row.get("FusionId", "yes")), support

    rows = []
    for _, event in events.iterrows():
        sv_id = str(event.get("sample_sv_id"))
        match = anno[anno[key].astype(str) == sv_id] if (key and len(anno)) else anno.iloc[0:0]
        record = {"sample_sv_id": sv_id}
        if len(match):
            r = match.iloc[0]
            p1 = int(float(r["pos1"])) if pd.notna(r.get("pos1")) else None
            p2 = int(float(r["pos2"])) if pd.notna(r.get("pos2")) else None
            d1, f1, t1 = nearest_sj(r.get("chrom1"), p1)
            d2, f2, t2 = nearest_sj(r.get("chrom2"), p2)
            nearest = min([d for d in (d1, d2) if d is not None], default=None)
            record.update(
                nearest_alt_sj_bp=nearest,
                alt_sj_frags=f1 if (d1 is not None and d1 == nearest) else f2,
                alt_sj_type=t1 if (d1 is not None and d1 == nearest) else t2,
                alt_sj_within_window=(nearest is not None
                                      and nearest <= criteria.RNA_JUNCTION_MAX_DIST))
            fusion_id, support = fusion_at(r.get("chrom1"), p1, r.get("chrom2"), p2)
            record.update(isofox_fusion=fusion_id, isofox_fusion_support=support)
            if not introns.empty and "GeneName" in introns.columns:
                genes = {str(r.get("gene1")), str(r.get("gene2"))}
                hit = introns[introns["GeneName"].astype(str).isin(genes)]
                record["retained_intron"] = len(hit) or None
        rows.append(record)
    return pd.DataFrame(rows)


# ============================================================================
# Stage 8 — junction evidence (the arbiter)
# ============================================================================

def _size_or_none(value):
    """Normalise an event size to int or None.

    NaN is TRUTHY in Python, so `if value:` lets a missing size through and
    `int(float(nan))` then raises. Every size arriving from a DataFrame can be
    NaN — inter-chromosomal events have no size by definition — so the check has
    to be explicit rather than relying on truthiness.
    """
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if number != number else int(number)      # NaN != NaN


def select_test(svtype, chrom1, chrom2, event_size=None,
                insert_len=None) -> tuple[str, str]:
    """Choose the only valid RNA test for this event geometry (RULE 3).

    Returns (test, reason), test in {chimeric, sizegap, insertion, none}.

    Order matters: an insertion-driven event is tested as an insertion even when
    its breakend span is tiny, because the span is not the lesion. An event whose
    span is 3 bp but which inserts 34 bp is an insertion; testing it for a 2 bp
    gap and recording `NONE` is a fabricated negative.
    """
    insert = _size_or_none(insert_len)
    if insert and insert >= criteria.MIN_INSERT_LEN:
        return "insertion", f"insertion-driven (insert_len={insert})"
    if str(chrom1) != str(chrom2) or str(svtype or "").upper() in {"BND", "TRA"}:
        return "chimeric", "inter-chromosomal: needs SA-to-partner / fusion caller"
    size = _size_or_none(event_size)
    if size is None:
        return "none", "intra-chromosomal with no resolvable size"
    if size < criteria.MIN_TESTABLE_GAP_SIZE:
        return "none", (f"{size} bp < {criteria.MIN_TESTABLE_GAP_SIZE} bp: below the "
                        "aligner's minimum intron, no N-gap can exist")
    return "sizegap", f"intra-chromosomal {size} bp: CIGAR N-gap test"


def sa_hits_partner(sa_tag: str, partner_chrom, partner_pos,
                    tolerance: int = None) -> bool:
    """True only if a supplementary alignment lands ON the partner breakend.

    Chromosome identity alone is not sufficient: the reads that produced a false
    STRONG carried SA tags to the right chromosome but 3-156 Mb away.
    """
    limit = criteria.SA_PARTNER_TOLERANCE if tolerance is None else tolerance
    if not sa_tag:
        return False
    want = str(partner_chrom).replace("chr", "")
    for record in str(sa_tag).split(";"):
        parts = record.split(",")
        if len(parts) < 2:
            continue
        chrom, pos = parts[0].replace("chr", ""), parts[1]
        if chrom == want and pos.lstrip("-").isdigit() \
                and abs(int(pos) - int(partner_pos)) <= limit:
            return True
    return False


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTacgtN", "TGCAtgcaN"))[::-1]


def count_at_breakend(bam, chrom, pos, test: str, event_size=None,
                      partner=None, insert_seq=None) -> dict:
    """Evidence at ONE breakend. `junction` is the only field that tiers."""
    contig = str(chrom) if str(chrom).startswith("chr") else f"chr{chrom}"
    out = {"coverage": 0, "softclip": 0, "junction": 0,
           "junction_by_ngap": 0, "junction_by_sa": 0, "junction_by_insert": 0,
           "alignments": 0, "low_mapq": 0}
    # Support is counted per FRAGMENT: one fragment may produce several alignment
    # records (supplementary/chimeric), and counting records inflates the total.
    # The per-mechanism sets are kept separate so a reader can see WHICH kind of
    # evidence a tier rests on, not just how much.
    junction_fragments: set = set()
    by_mechanism = {"ngap": set(), "sa": set(), "insert": set()}
    try:
        reads = bam.fetch(contig, max(0, int(pos) - criteria.COVERAGE_WINDOW),
                          int(pos) + criteria.COVERAGE_WINDOW)
    except (ValueError, KeyError):
        return out

    size = _size_or_none(event_size)
    for read in reads:
        if read.is_unmapped or not read.cigartuples:
            continue
        out["coverage"] += 1
        out["alignments"] += 1
        # Ambiguously placed reads are not evidence for a specific locus.
        if read.mapping_quality < criteria.MIN_READ_MAPQ:
            out["low_mapq"] += 1
            continue

        # Diagnostic only — soft-clips never tier (RULE 1). A clip shows a read
        # ENDS at the breakpoint; it says nothing about where the rest went.
        for op, length in read.cigartuples:
            if op == 4 and length >= criteria.MIN_SOFTCLIP_LEN:
                if min(abs(read.reference_start - pos),
                       abs(read.reference_end - pos)) <= criteria.NGAP_POSITION_TOLERANCE:
                    out["softclip"] += 1
                    break

        if test == "chimeric" and partner and read.has_tag("SA"):
            if sa_hits_partner(read.get_tag("SA"), partner[0], partner[1]):
                junction_fragments.add(read.query_name)
                by_mechanism["sa"].add(read.query_name)

        elif test == "sizegap" and size:
            # A spliced/deleted gap is CIGAR `N` (skip), NOT `D`. Searching for
            # `D` found 0 reads on a 221 bp deletion where `N` found 31. The
            # gene's own introns are also `N` but large and scattered; the
            # discriminators are size AND a fixed position.
            cursor = read.reference_start
            for op, length in read.cigartuples:
                if op == 3 and abs(length - size) <= criteria.NGAP_SIZE_TOLERANCE \
                        and abs(cursor - pos) <= criteria.NGAP_POSITION_TOLERANCE:
                    junction_fragments.add(read.query_name)
                    by_mechanism["ngap"].add(read.query_name)
                    break
                if op in (0, 2, 3, 7, 8):     # reference-consuming ops
                    cursor += length

        elif test == "insertion" and insert_seq:
            seq = read.query_sequence or ""
            if insert_seq in seq or _revcomp(insert_seq) in seq:
                junction_fragments.add(read.query_name)
                by_mechanism["insert"].add(read.query_name)
    out["junction"] = len(junction_fragments)
    for mechanism, names in by_mechanism.items():
        out[f"junction_by_{mechanism}"] = len(names)
    return out


def rna_tier(junction_reads: int, testable: bool, has_rna: bool = True) -> str:
    """Tier from junction reads alone. Coverage is deliberately NOT an argument.

    The workflow this replaces promoted an event on `coverage > background`, with
    zero crossing reads — under which any breakpoint inside an expressed gene
    looks validated.
    """
    if not has_rna:
        return "NA"
    if not testable:
        return "UNTESTABLE"
    if junction_reads >= criteria.STRONG_MIN_JUNCTION_READS:
        return "STRONG"
    if junction_reads >= criteria.SUGGESTIVE_MIN_JUNCTION_READS:
        return "SUGGESTIVE"
    return "WEAK" if junction_reads >= 1 else "NONE"


def evaluate(events: pd.DataFrame, bam_path: str | None) -> pd.DataFrame:
    """Stage 8 over an event table. Missing BAM -> every tier is `NA`."""
    if events.empty:
        return pd.DataFrame()

    bam = None
    if bam_path and os.path.exists(bam_path):
        import pysam
        bam = pysam.AlignmentFile(bam_path, "rb")

    rows = []
    for _, event in events.iterrows():
        test, reason = select_test(event.get("svtype"), event.get("chrom1"),
                                   event.get("chrom2"), event.get("event_size"),
                                   event.get("insert_len"))
        row = {"sample_sv_id": event.get("sample_sv_id"), "gene": event.get("gene"),
               "test": test, "test_reason": reason}

        if bam is None:
            row.update(rna_tier="NA", junction_reads=None)
        elif test == "none":
            row.update(rna_tier="UNTESTABLE", junction_reads=None)
        else:
            inter = str(event.get("chrom1")) != str(event.get("chrom2"))
            partner12 = (event["chrom2"], int(float(event["pos2"]))) if inter else None
            partner21 = (event["chrom1"], int(float(event["pos1"]))) if inter else None
            bp1 = count_at_breakend(bam, event["chrom1"], int(float(event["pos1"])),
                                    test, event.get("event_size"), partner12,
                                    event.get("insert_seq"))
            bp2 = count_at_breakend(bam, event["chrom2"], int(float(event["pos2"])),
                                    test, event.get("event_size"), partner21,
                                    event.get("insert_seq"))
            # The one legitimate max: a crossing read is visible from either end,
            # so the same read must not be counted twice.
            junction = max(bp1["junction"], bp2["junction"])
            row.update(coverage_bp1=bp1["coverage"], coverage_bp2=bp2["coverage"],
                       min_coverage=min(bp1["coverage"], bp2["coverage"]),
                       softclip_bp1=bp1["softclip"], softclip_bp2=bp2["softclip"],
                       # Which mechanism the support rests on, and how much was
                       # discarded for ambiguous placement — a tier backed by
                       # one mechanism at one breakend reads differently from
                       # the same number spread across both.
                       junction_by_ngap=max(bp1["junction_by_ngap"], bp2["junction_by_ngap"]),
                       junction_by_sa=max(bp1["junction_by_sa"], bp2["junction_by_sa"]),
                       junction_by_insert=max(bp1["junction_by_insert"], bp2["junction_by_insert"]),
                       alignments_bp1=bp1["alignments"], alignments_bp2=bp2["alignments"],
                       low_mapq_bp1=bp1["low_mapq"], low_mapq_bp2=bp2["low_mapq"],
                       junction_reads=junction,
                       rna_tier=rna_tier(junction, testable=True))
        rows.append(row)

    if bam is not None:
        bam.close()
    return pd.DataFrame(rows)
