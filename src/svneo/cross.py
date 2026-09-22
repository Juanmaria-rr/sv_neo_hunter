"""
cross.py — stages 3-5: the cross against the reference, sequence QC, and the
collapse from peptide rows to genomic events.

STAGE 3 REPORTS THREE LEVELS, ALWAYS TOGETHER
---------------------------------------------
    level 1  identical peptide      strongest claim, and the narrowest
    level 2  same gene + same SV type   mechanistic recurrence with a different
                                        junction (and therefore a different peptide)
    level 3  breakpoint proximity    reported as a window GRADIENT, never one cutoff

Reporting only level 1 under-claims (a real shared mechanism with a shifted
breakpoint is invisible); reporting only level 2 over-claims (the same gene is
broken in unrelated genomes surprisingly often, especially at fragile sites).

STAGE 5 IS THE MOST CONSEQUENTIAL STEP IN THE WHOLE PIPELINE
------------------------------------------------------------
One SV produces up to ~40 overlapping sliding-window peptides, each able to match
the reference separately. Counting peptide rows overstates the number of
independent genomic findings — observed 21-fold on one 53 bp deletion. Everything
downstream, every headline number and every figure, uses EVENTS.

Two grains are emitted (criteria.REPORT_BOTH_EVENT_GRAINS):
    junction grain  one row per VCF junction         <- the honest unit
    grouped grain   the caller's own grouping id     <- for continuity
They differ: a caller's grouped id merges clustered junctions at one locus and
can double-count a breakend shared by two groups (observed 46 -> 39).
"""
from __future__ import annotations

import gzip

import pandas as pd

from . import criteria


# ============================================================================
# Reference loading and the self-proteome index
# ============================================================================

def load_reference(reference) -> pd.DataFrame:
    """Load the query peptide catalogue and normalise its column names."""
    sep = "," if reference.peptides.endswith(".csv") else "\t"
    frame = pd.read_csv(reference.peptides, sep=sep, dtype=str)
    if reference.peptide_column not in frame.columns:
        raise ValueError(f"peptide column '{reference.peptide_column}' not in "
                         f"{reference.peptides}; columns are {list(frame.columns)}")
    out = frame.rename(columns={reference.peptide_column: "ref_peptide"})
    if reference.gene_column and reference.gene_column in frame.columns:
        out = out.rename(columns={reference.gene_column: "ref_gene"})
    if reference.svtype_column and reference.svtype_column in frame.columns:
        out = out.rename(columns={reference.svtype_column: "ref_svtype"})
    out["ref_peptide"] = out["ref_peptide"].astype(str).str.strip().str.upper()
    return out


def load_proteome(fasta_gz: str) -> str:
    """The proteome as ONE searchable string for the self-peptide test.

    Sequences are joined with a NUL separator so no peptide can match across a
    protein boundary — plain concatenation would invent junction sequences and
    wrongly mark real neopeptides as self.
    """
    sequences, buffer = [], []
    opener = gzip.open if fasta_gz.endswith(".gz") else open
    with opener(fasta_gz, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if buffer:
                    sequences.append("".join(buffer))
                    buffer = []
            else:
                buffer.append(line.strip())
    if buffer:
        sequences.append("".join(buffer))
    return "\0".join(sequences)


# ============================================================================
# Stage 3 — the three levels
# ============================================================================

def level1_identical(peptides: pd.DataFrame, anno: pd.DataFrame,
                     reference: pd.DataFrame) -> pd.DataFrame:
    """Exact peptide-string intersection, annotated with concordance flags.

    `gene_concordant` is the key sanity signal: coincidental short-peptide
    convergence would not preferentially land in the SAME gene, so a high
    concordance fraction among matches argues for genuine shared loci. It is not
    a substitute for the null model (see null_model.py).
    """
    column = criteria.PEPTIDE_COLUMN
    candidates = peptides.copy()

    # A peptide lying entirely on one side of the breakpoint is not SV-derived:
    # it is an ordinary peptide of the intact reading frame, and it can match the
    # reference while carrying no information about the structural variant.
    flag = criteria.SPANS_JUNCTION_COLUMN
    if criteria.REQUIRE_SPANS_JUNCTION and flag in candidates.columns:
        before = len(candidates)
        spans = candidates[flag].astype(str).str.lower().isin(("true", "1"))
        candidates = candidates[spans]
        if before != len(candidates):
            print(f"    spans_junction filter: {before} -> {len(candidates)} "
                  f"peptide rows ({before - len(candidates)} did not cross)")

    candidates["peptide"] = candidates[column].astype(str).str.strip().str.upper()

    merged = candidates.merge(reference, left_on="peptide", right_on="ref_peptide",
                              how="inner")
    if merged.empty:
        return merged

    # Attach the sample-side annotation (gene, svtype) for the SV that produced
    # each peptide, so concordance can be evaluated.
    keys = [c for c in ("sv_id", "sample_sv_id") if c in merged.columns
            and c in anno.columns]
    if keys:
        merged = merged.merge(anno, on=keys[0], how="left", suffixes=("", "_anno"))

    # Concordance is evaluated against BOTH breakends. A junction joins two
    # genes, and the reference peptide may be annotated to either of them: in an
    # observed run, 8 of 64 matches were concordant on gene2 while gene1 differed
    # (both translocations), and testing only gene1 discarded them.
    if "ref_gene" in merged.columns and ("gene1" in merged.columns
                                         or "gene2" in merged.columns):
        reference_gene = merged["ref_gene"].astype(str).str.upper().str.strip()
        concordant = None
        for side in ("gene1", "gene2"):
            if side not in merged.columns:
                continue
            this_side = merged[side].astype(str).str.upper().str.strip()
            match = this_side == reference_gene
            concordant = match if concordant is None else (concordant | match)
            merged[f"concordant_{side}"] = match
        merged["gene_concordant"] = concordant
    else:
        # With no gene column on one side, concordance cannot be evaluated. Do
        # not default it to True: that would silently pass everything through
        # the credibility filter.
        merged["gene_concordant"] = pd.NA

    if "svtype" in merged.columns and "ref_svtype" in merged.columns:
        merged["svtype_concordant"] = (
            merged["svtype"].astype(str).str.upper()
            == merged["ref_svtype"].astype(str).str.upper())
    return merged


def level2_gene_svtype(anno: pd.DataFrame, reference: pd.DataFrame) -> pd.DataFrame:
    """Same gene broken by the same class of SV, independent of peptide sequence."""
    if "ref_gene" not in reference.columns or "gene1" not in anno.columns:
        return pd.DataFrame()
    sample = anno.assign(_gene=anno["gene1"].astype(str).str.upper().str.strip())
    ref = reference.assign(_gene=reference["ref_gene"].astype(str).str.upper().str.strip())
    merged = sample.merge(ref, on="_gene", how="inner")
    if "svtype" in merged.columns and "ref_svtype" in merged.columns:
        merged["svtype_match"] = (merged["svtype"].astype(str).str.upper()
                                  == merged["ref_svtype"].astype(str).str.upper())
    return merged


def level3_proximity(junctions: pd.DataFrame, pool: pd.DataFrame,
                     windows_kb=(1, 10, 50, 100)) -> pd.DataFrame:
    """Distance from each reference breakpoint to the nearest sample breakpoint.

    Reported as a gradient because the count grows smoothly with window size and
    any single cutoff is arbitrary. Both "at least one breakend within window"
    and "both breakends within window" are given: they answer different questions
    and the first is far more permissive.

    COMPARABILITY WARNING, from an error that produced a wrong conclusion once:
    when comparing branches or samples at this level, BOTH sides must use the
    identical reference call. Comparing a curated-subset figure against a
    full-pool figure suggested a filter had *increased* proximity hits.
    """
    if junctions.empty or pool is None or pool.empty:
        return pd.DataFrame()
    windows_kb = windows_kb or criteria.PROXIMITY_WINDOWS_KB

    by_chrom: dict[str, list[int]] = {}
    for _, row in junctions.iterrows():
        for chrom, pos in ((row["chrom1"], row["pos1"]), (row["chrom2"], row["pos2"])):
            by_chrom.setdefault(str(chrom).replace("chr", ""), []).append(int(pos))
    for chrom in by_chrom:
        by_chrom[chrom].sort()

    def nearest(chrom, pos):
        positions = by_chrom.get(str(chrom).replace("chr", ""))
        if not positions:
            return None
        import bisect
        i = bisect.bisect_left(positions, int(pos))
        best = None
        for j in (i - 1, i):
            if 0 <= j < len(positions):
                distance = abs(positions[j] - int(pos))
                best = distance if best is None else min(best, distance)
        return best

    rows = []
    for _, row in pool.iterrows():
        d1 = nearest(row.get("chrom1"), row.get("pos1")) if row.get("pos1") else None
        d2 = nearest(row.get("chrom2"), row.get("pos2")) if row.get("pos2") else None
        record = {"ref_peptide": row.get("ref_peptide"),
                  "dist_bp1": d1, "dist_bp2": d2}
        for kb in windows_kb:
            limit = kb * 1000
            record[f"any_within_{kb}kb"] = any(
                d is not None and d <= limit for d in (d1, d2))
            record[f"both_within_{kb}kb"] = all(
                d is not None and d <= limit for d in (d1, d2))
        rows.append(record)
    return pd.DataFrame(rows)


def junction_recurrence(junctions: pd.DataFrame, pool: pd.DataFrame,
                        anno: pd.DataFrame = None) -> pd.DataFrame:
    """For each of OUR junctions, how close is the nearest patient breakpoint?

    `level3_proximity()` answers the mirror question — for each patient
    breakpoint, is one of ours nearby — which is the right form for summarising a
    cohort but cannot be joined to a per-junction table. This one is keyed on our
    `sample_sv_id`, so every admitted junction can carry an explicit answer to
    "was anything like this seen in patients?".

    Distance is reported in base pairs, not thresholded: a breakpoint 300 bp away
    and one 90 kb away are both "within 100 kb" and mean entirely different
    things. `patient_svtype` and `patient_gene` name what was found there, so a
    match can be judged rather than merely counted.
    """
    if junctions.empty or pool is None or pool.empty:
        return pd.DataFrame()

    import bisect
    needed = {"chrom1", "pos1", "chrom2", "pos2"}
    if not needed.issubset(pool.columns):
        return pd.DataFrame()

    # Patient breakpoints indexed per chromosome, sorted for binary search.
    by_chrom: dict[str, list[tuple]] = {}
    for _, row in pool.iterrows():
        for chrom, pos in ((row["chrom1"], row["pos1"]), (row["chrom2"], row["pos2"])):
            if pd.isna(chrom) or pd.isna(pos):
                continue
            key = str(chrom).replace("chr", "")
            by_chrom.setdefault(key, []).append(
                (int(float(pos)), str(row.get("svtype", "")), str(row.get("gene1", ""))))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    positions = {c: [p for p, _, _ in v] for c, v in by_chrom.items()}

    def nearest(chrom, pos):
        chrom = str(chrom).replace("chr", "")
        entries = by_chrom.get(chrom)
        if not entries or pd.isna(pos):
            return None, None, None
        pos = int(float(pos))
        i = bisect.bisect_left(positions[chrom], pos)
        best = None
        for j in (i - 1, i):
            if 0 <= j < len(entries):
                distance = abs(entries[j][0] - pos)
                if best is None or distance < best[0]:
                    best = (distance, entries[j][1], entries[j][2])
        return best if best else (None, None, None)

    rows = []
    for _, junction in junctions.iterrows():
        d1, type1, gene1 = nearest(junction["chrom1"], junction["pos1"])
        d2, type2, gene2 = nearest(junction["chrom2"], junction["pos2"])
        closest = min([d for d in (d1, d2) if d is not None], default=None)
        rows.append({
            "sample_sv_id": str(junction["sample_sv_id"]),
            "patient_bp_dist_bp": closest,
            "patient_bp_dist_bp1": d1, "patient_bp_dist_bp2": d2,
            "patient_bp_both_within_10kb": (d1 is not None and d2 is not None
                                            and max(d1, d2) <= 10_000),
            "patient_svtype": type1 if (d1 is not None and d1 == closest) else type2,
            "patient_gene": gene1 if (d1 is not None and d1 == closest) else gene2,
        })
    return pd.DataFrame(rows)


# ============================================================================
# Stage 4 — sequence QC
# ============================================================================

def sequence_qc(matches: pd.DataFrame, proteome: str | None) -> pd.DataFrame:
    """Flag low-complexity and self peptides, then apply the credible rule."""
    if matches.empty:
        return matches
    out = matches.copy()
    out["low_complexity"] = out["peptide"].map(criteria.is_low_complexity)
    out["is_self"] = out["peptide"].map(
        lambda p: (p in proteome) if proteome else False)
    if proteome is None:
        out["is_self"] = pd.NA   # not tested is not the same as tested-negative
    out["credible"] = [
        criteria.is_credible(bool(lc), bool(self_) if pd.notna(self_) else False,
                             bool(gc) if pd.notna(gc) else False)
        for lc, self_, gc in zip(out["low_complexity"], out["is_self"],
                                 out["gene_concordant"])]
    return out


# ============================================================================
# Stage 5 — collapse to genomic events
# ============================================================================

def to_events(matches: pd.DataFrame, junctions: pd.DataFrame) -> pd.DataFrame:
    """Collapse credible peptide matches to one row per genomic event.

    Emits the junction grain. The number of peptides behind each event travels
    with it (`n_peptides`) so the collapse is visible rather than hidden.
    """
    if matches.empty:
        return pd.DataFrame()

    credible = matches[matches["credible"].fillna(False)]
    if credible.empty:
        return pd.DataFrame()

    key = criteria.EVENT_DEDUP_COLUMN
    if key not in credible.columns:
        key = next((c for c in ("sv_id", "sample_sv_id", "junction_key")
                    if c in credible.columns), None)
        if key is None:
            raise ValueError("no usable event-deduplication key in the match table")

    grouped = credible.groupby(key).agg(
        n_peptides=("peptide", "nunique"),
        peptides=("peptide", lambda s: ";".join(sorted(set(s))[:10])),
        ref_genes=("ref_gene", lambda s: ";".join(sorted(set(s.dropna().astype(str))))
                   if "ref_gene" in credible.columns else ""),
    ).reset_index().rename(columns={key: "sample_sv_id"})

    if not junctions.empty:
        cols = ["sample_sv_id", "junction_key", "chrom1", "pos1", "chrom2", "pos2",
                "svtype", "span", "event_size", "insert_len", "insert_seq",
                "pon_count", "vf_bp1", "vf_bp2", "qual_bp1", "qual_bp2",
                "segmapq_bp1", "segmapq_bp2"]
        available = [c for c in cols if c in junctions.columns]
        grouped = grouped.merge(junctions[available].astype({"sample_sv_id": str}),
                                on="sample_sv_id", how="left")
    grouped["gene"] = grouped.get("ref_genes")
    return grouped
