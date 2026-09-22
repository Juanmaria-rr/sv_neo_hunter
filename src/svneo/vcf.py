"""
vcf.py — stage 1: read an SV VCF, admit records, pair breakends into junctions.

Everything downstream counts junctions, so this module is where the grain of the
whole analysis is fixed. Three things it gets right that are easy to get wrong:

1. THE BREAKEND SPAN IS NOT THE EVENT SIZE. |pos2 - pos1| is the distance between
   breakends, which equals the lesion only for DUP:
       DEL  -> span - 1 bases are deleted
       INS  -> span is ALWAYS 1 (the breakends are adjacent); the real size is
               the length of the inserted sequence in the ALT
       DUP  -> size = span
       BND  -> undefined
   Using the span made 248 germline insertions look like 1 bp events, and made a
   34 bp insertion look like a 2 bp deletion. Both `span` and the type-aware
   `event_size` are emitted, and downstream code uses `event_size`.

2. INSERTED SEQUENCE MUST SURVIVE TO THE RNA STAGE. `insert_len` and `insert_seq`
   are parsed here because the only valid RNA test for an insertion-driven event
   is searching for that sequence — a gap test can never confirm it.

3. ADMISSION DIFFERS BY VCF KIND. Somatic call sets arrive panel-filtered by the
   caller; germline/parental sets are typically only ANNOTATED with PON_COUNT and
   pass FILTER unchallenged. Applying the same rule to both leaves the parental
   side unfiltered while the somatic side is filtered — an asymmetry that made an
   unfiltered call set look richer in "recurrent" hits than a filtered one.
"""
from __future__ import annotations

import gzip
import re

import pandas as pd

from . import criteria

#: ALT patterns for a mate-carrying breakend, e.g. `G]chr5:123456]`.
_BRACKET = re.compile(r"[\[\]]")
#: The inserted/novel sequence in a breakend ALT is the run of bases outside the
#: bracketed mate coordinate, minus the reference base.
_ALT_BASES = re.compile(r"^([ACGTNacgtn]*)[\[\]]|[\[\]]([ACGTNacgtn]*)$")


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _info(text: str, key: str, cast=str):
    """Extract one INFO field. Returns True for a valueless flag, None if absent."""
    match = re.search(rf"(?:^|;){key}=([^;]+)", text)
    if not match:
        return True if re.search(rf"(?:^|;){key}(?:;|$)", text) else None
    try:
        return cast(match.group(1))
    except (TypeError, ValueError):
        return match.group(1)


def _insert_sequence(alt: str, ref: str) -> str:
    """Novel bases carried by a breakend ALT, excluding the reference anchor."""
    if not _BRACKET.search(alt):
        return alt[len(ref):] if len(alt) > len(ref) else ""
    match = _ALT_BASES.search(alt)
    if not match:
        return ""
    bases = match.group(1) or match.group(2) or ""
    # The anchor base is the reference allele, at whichever end it sits.
    if bases.upper().startswith(ref.upper()):
        return bases[len(ref):]
    if bases.upper().endswith(ref.upper()):
        return bases[:-len(ref)]
    return bases


def read_breakends(path: str) -> pd.DataFrame:
    """Every VCF record as one row, with the fields that carry evidence."""
    rows = []
    with _open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, filt, info = f[:8]
            rows.append({
                "chrom": chrom.replace("chr", ""), "pos": int(pos), "id": vid,
                "ref": ref, "alt": alt, "filter": filt,
                "qual": float(qual) if qual not in (".", "") else None,
                "paired": bool(_BRACKET.search(alt)),
                "svtype": _info(info, "SVTYPE"),
                "mateid": _info(info, "MATEID"),
                "pon_count": _info(info, "PON_COUNT", int),
                "vf": _info(info, "VF", int), "sf": _info(info, "SF", int),
                "df": _info(info, "DF", int), "ref_reads": _info(info, "REF", int),
                "segmapq": _info(info, "SEGMAPQ", int),
                "homseq": _info(info, "HOMSEQ"),
                "imprecise": bool(_info(info, "IMPRECISE")),
                "insert_seq": _insert_sequence(alt, ref),
            })
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["insert_len"] = frame["insert_seq"].str.len()
    return frame


def admit(breakends: pd.DataFrame, vcf_kind: str,
          pon_max: int | None = None,
          admit_panel_filtered: bool = False) -> tuple[pd.DataFrame, dict]:
    """Apply stage-1 admission. Returns the admitted rows and a funnel dict.

    `pon_max=None` disables this pipeline's own PON threshold, which is all that
    is needed on a germline call set. `admit_panel_filtered=True` additionally
    admits records the CALLER already rejected as panel hits (`FILTER=PON`),
    which is what a somatic call set needs — there, `pon_max` never sees those
    records at all. Both are required for a branch that reports the panel
    without filtering on it; see `criteria.SOMATIC_KEEP_FILTERS_REPORT_PON`.

    The funnel records what each rule removed, because "how many were lost
    where" is the part of a cascade that most often turns out to be the finding.
    """
    funnel = {"records": len(breakends)}
    kept = breakends

    if vcf_kind == "somatic":
        # The caller already applied its panel filter here; FILTER carries it.
        keep = (criteria.SOMATIC_KEEP_FILTERS_REPORT_PON if admit_panel_filtered
                else criteria.SOMATIC_KEEP_FILTERS)
        if admit_panel_filtered:
            # Reported, not filtered: state how many the caller had rejected, or
            # the relaxation is invisible in the output.
            funnel["caller_pon_admitted"] = int((kept["filter"] == "PON").sum())
        kept = kept[kept["filter"].isin(keep)]
        funnel["after_filter_pass"] = len(kept)
    else:
        # Germline/parental sets are usually all-PASS by construction, so FILTER
        # alone admits everything; the population filters below do the work.
        funnel["after_filter_pass"] = len(kept)

    kept = kept[kept["paired"]]
    funnel["after_paired_only"] = len(kept)

    if pon_max is not None:
        pon = kept["pon_count"].fillna(criteria.PON_ABSENT_MEANS)
        kept = kept[pon < pon_max]
        funnel[f"after_pon_lt_{pon_max}"] = len(kept)

    return kept.copy(), funnel


def event_size(svtype: str, pos1: int, pos2: int, insert_len: int = 0):
    """Type-aware lesion size. Returns None where size is undefined (BND/SGL).

    This is the function that stops a 34 bp insertion from being reported as a
    2 bp deletion, and stops a 100 bp span from being read as 100 deleted bases.
    """
    kind = (svtype or "").upper()
    span = abs(int(pos2) - int(pos1))
    if kind in ("BND", "TRA", "SGL"):
        return None
    if kind == "INS":
        return int(insert_len) or None
    if kind == "DEL":
        return max(span - 1, 0)
    if kind in ("DUP", "INV"):
        return span
    return span


def pair_junctions(admitted: pd.DataFrame) -> pd.DataFrame:
    """Collapse mated breakends into one row per junction.

    Pairing is by MATEID when present, falling back to the bracketed coordinate
    in ALT. Each junction is emitted once, keyed on the sorted pair of breakend
    ids so a junction cannot be counted twice from its two ends.
    """
    if admitted.empty:
        return pd.DataFrame()

    by_id = {r["id"]: r for _, r in admitted.iterrows()}
    seen, junctions = set(), []

    for _, record in admitted.iterrows():
        mate_id = record["mateid"]
        mate = by_id.get(mate_id) if mate_id else None
        if mate is None:
            # Fall back to the coordinate inside the ALT brackets.
            match = re.search(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", record["alt"])
            if not match:
                continue
            mchrom, mpos = match.group(1).replace("chr", ""), int(match.group(2))
            mate = {"id": f"{mchrom}:{mpos}", "chrom": mchrom, "pos": mpos,
                    "pon_count": record["pon_count"], "vf": record["vf"],
                    "qual": record["qual"], "segmapq": record["segmapq"],
                    "sf": record["sf"], "df": record["df"],
                    "ref_reads": record["ref_reads"], "homseq": record["homseq"],
                    "imprecise": record["imprecise"], "insert_len": 0}

        key = tuple(sorted([str(record["id"]), str(mate["id"])]))
        if key in seen:
            continue
        seen.add(key)

        insert_len = max(int(record.get("insert_len") or 0),
                         int(mate.get("insert_len") or 0))
        junctions.append({
            "junction_key": "|".join(key),
            "sample_sv_id": str(record["id"]),
            "chrom1": record["chrom"], "pos1": int(record["pos"]),
            "chrom2": mate["chrom"], "pos2": int(mate["pos"]),
            "svtype": record["svtype"],
            "span": abs(int(mate["pos"]) - int(record["pos"]))
                    if record["chrom"] == mate["chrom"] else None,
            "event_size": event_size(record["svtype"], record["pos"], mate["pos"],
                                     insert_len)
                          if record["chrom"] == mate["chrom"] else None,
            "insert_len": insert_len,
            "insert_seq": record.get("insert_seq") or "",
            # Event-level PON is the MAX across breakends: if either end is a
            # recurrent panel site, the whole junction is suspect.
            "pon_count": max([v for v in (record["pon_count"], mate.get("pon_count"))
                              if v is not None], default=None),
            # Every evidence field the VCF carries, per breakend. Reporting the
            # verdict without the values it rests on makes a table unusable for
            # judging a call: "not high-confidence" does not say whether the
            # mapping quality, the fragment support, the caller quality or the
            # size floor was the reason.
            "vf_bp1": record["vf"], "vf_bp2": mate.get("vf"),
            "sf_bp1": record["sf"], "sf_bp2": mate.get("sf"),
            "df_bp1": record["df"], "df_bp2": mate.get("df"),
            "ref_bp1": record["ref_reads"], "ref_bp2": mate.get("ref_reads"),
            "qual_bp1": record["qual"], "qual_bp2": mate.get("qual"),
            "segmapq_bp1": record["segmapq"], "segmapq_bp2": mate.get("segmapq"),
            "homseq_bp1": record["homseq"], "homseq_bp2": mate.get("homseq"),
            "imprecise_bp1": record["imprecise"], "imprecise_bp2": mate.get("imprecise"),
            "filter": record["filter"], "alt_bp1": record["alt"],
        })
    return pd.DataFrame(junctions)


def write_admitted_vcf(source_vcf: str, admitted: pd.DataFrame, out_path: str) -> int:
    """Write the admitted records to a plain VCF, preserving the header.

    THIS IS NOT COSMETIC. The peptide generator must see the ADMITTED call set,
    not the raw file: it applies no FILTER or panel logic of its own, so handing
    it the raw VCF builds peptides from PON and INFERRED breakends, wastes the
    generation, and — worse — leaves artefact-derived peptides in the candidate
    universe where they can match a reference peptide by chance.

    Observed on a somatic VCF: 350 raw records (108 PASS, 82 PON, 160 INFERRED),
    of which the generator happily loaded 176. Only the PASS subset should ever
    reach it.

    Output is uncompressed because the generator checks the file *extension* and
    rejects `.vcf.gz` outright.
    """
    # Select by record ID, not by coordinate. Two breakends can sit at the SAME
    # position and belong to different junctions — 16 such positions in one
    # observed germline call set. Filtering on (chrom, pos) then writes out a
    # rejected record because its neighbour passed: a common germline breakend
    # can reach the generator that way and produce peptides, which is precisely
    # what handing the generator the admitted set prevents.
    ids = set(admitted["id"].astype(str))
    usable_ids = "." not in ids and len(ids) == len(admitted)
    keys = set(zip(admitted["chrom"].astype(str), admitted["pos"].astype(int)))

    written = 0
    with _open(source_vcf) as src, open(out_path, "w") as dst:
        for line in src:
            if line.startswith("#"):
                dst.write(line)
                continue
            fields = line.split("\t", 3)
            keep = (fields[2] in ids) if usable_ids else \
                ((fields[0].replace("chr", ""), int(fields[1])) in keys)
            if keep:
                dst.write(line)
                written += 1
    return written


def panel_size_estimate(breakends: pd.DataFrame) -> int:
    """Lower-bound panel size, so PON_COUNT can be reported as a fraction.

    max(PON_COUNT) over the file is a floor on the true panel size: the panel
    cannot be smaller than the largest count observed in it.
    """
    if breakends.empty or breakends["pon_count"].isna().all():
        return 0
    return int(breakends["pon_count"].max())
