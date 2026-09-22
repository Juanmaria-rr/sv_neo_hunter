"""
confidence.py — stage 6: is the SV call itself believable, and is it private?

Two orthogonal questions, deliberately kept apart:

    CONFIDENCE  does the evidence support this junction being real?
                (read support, mapping quality, caller confidence, geometry)
    PRIVACY     is it specific to this sample, or seen in the general population?
                (panel of normals AND population SV frequency)

They are not interchangeable, and the distinction decides real cases: an event can
carry perfect confidence metrics — high fragment support, high mapping quality,
high caller quality, strong expression, junction-crossing RNA reads and a
homozygous-deletion depth profile in DNA — and still be disqualified, because a
high panel-of-normals count means the junction is present in many unrelated
normals. Such an event is real, well-called, well-expressed and a COMMON
POLYMORPHISM.
Confidence says the call is right; privacy says it is not the sample's own.

HC IS REPORTED, NEVER APPLIED AS A FILTER
-----------------------------------------
A non-HC event with orthogonal RNA evidence can be more believable than an HC
event with none. Deleting non-HC events would hide that.

COPY-NUMBER FIELDS ARE EXCLUDED BY DEFAULT
------------------------------------------
`PURPLE_AF` / `PURPLE_CN` / `PURPLE_JCN` depend on the caller's purity/ploidy
fit, which is unreliable for clonal cell lines — an observed fit of purity 0.37
and ploidy 3.3 on a clonal pair that should sit near 1.0. Raw evidence fields are
unaffected. See criteria.TRUST_COPY_NUMBER_FIELDS.
"""
from __future__ import annotations

import os
import subprocess

import pandas as pd

from . import criteria


def annotate_confidence(events: pd.DataFrame) -> pd.DataFrame:
    """Attach per-breakend evidence fields and the event-level HC verdict."""
    if events.empty:
        return events
    out = events.copy()

    def hc(row) -> bool:
        same_chrom = str(row.get("chrom1")) == str(row.get("chrom2"))
        size = row.get("event_size") if same_chrom else None
        breakends = [
            {"SEGMAPQ": row.get("segmapq_bp1"), "VF": row.get("vf_bp1"),
             "QUAL": row.get("qual_bp1"), "size": size},
            {"SEGMAPQ": row.get("segmapq_bp2"), "VF": row.get("vf_bp2"),
             "QUAL": row.get("qual_bp2"), "size": size},
        ]
        return criteria.event_is_hc(breakends)

    out["sv_hc"] = out.apply(hc, axis=1)
    return out


def annotate_privacy(events: pd.DataFrame, panel_size: int = 0,
                     gnomad_af: dict | None = None,
                     pon_in_privacy: bool = True) -> pd.DataFrame:
    """Attach the two recurrence filters and the combined `is_private` verdict.

    `gnomad_af` maps sample_sv_id -> population allele frequency (see
    `annotate_gnomad()`); absent means "not matched in gnomAD", treated as AF 0.

    `pon_in_privacy=False` keeps the panel count as an annotation only, so
    `is_private` rests on population frequency alone. This is the SECOND place
    the panel filter acts: relaxing admission alone still lets `is_private`
    remove the very candidates the relaxed branch was run to see. `pass_pon` is
    still computed and reported either way — what changes is whether it votes.
    """
    if events.empty:
        return events
    out = events.copy()

    out["pass_pon"] = out["pon_count"].map(criteria.passes_pon)
    if panel_size:
        out["pon_fraction"] = out["pon_count"].map(
            lambda c: round(criteria.pon_fraction(c or 0, panel_size), 6))
        out["panel_size_estimate"] = panel_size

    # Every ancestry group is carried into the table. Which one is the right
    # reference depends on the donor's ancestry — knowledge the pipeline does not
    # have — so the columns are emitted and the verdict below is only a default.
    columns = ["gnomad_af", "gnomad_af_popmax", "gnomad_af_popmax_pop"] + \
        [f"gnomad_af_{p}" for p in criteria.GNOMAD_POPULATIONS]
    if gnomad_af:
        mapped = out["sample_sv_id"].astype(str).map(gnomad_af)
        for column in columns:
            out[column] = mapped.map(
                lambda record: record.get(column) if isinstance(record, dict) else None)
    else:
        # Not evaluated is NOT the same as "absent from gnomAD". Recorded as NA
        # so a run without the resource cannot be mistaken for a clean one.
        for column in columns:
            out[column] = pd.NA

    field = criteria.GNOMAD_AF_FIELD
    judged_on = {"popmax": "gnomad_af_popmax", "global": "gnomad_af"}.get(
        field, f"gnomad_af_{field}")
    out["gnomad_af_used"] = out.get(judged_on, pd.Series([pd.NA] * len(out)))
    out["gnomad_af_field"] = judged_on
    out["pass_gnomad"] = out["gnomad_af_used"].map(
        lambda af: True if pd.isna(af) else float(af) < criteria.GNOMAD_MAX_AF)
    out["pass_pon"] = out["pass_pon"].fillna(False).astype(bool)
    out["pon_in_privacy"] = pon_in_privacy
    out["is_private"] = (out["pass_pon"] & out["pass_gnomad"]) if pon_in_privacy \
        else out["pass_gnomad"]

    # The note has to distinguish three different states that would otherwise all
    # read as a bare `is_private` boolean: gnomAD missing, the panel deliberately
    # not voting, and both criteria applied.
    notes = []
    for af, passed_pon in zip(out["gnomad_af"], out["pass_pon"]):
        if not pon_in_privacy:
            notes.append("gnomAD only — panel count reported, not applied"
                         + ("" if passed_pon else "; this event IS in the panel"))
        elif pd.isna(af):
            notes.append("gnomAD not evaluated — PON only")
        else:
            notes.append("")
    out["privacy_note"] = notes
    return out


def annotate_gnomad(events: pd.DataFrame, gnomad_vcf: str,
                    helper: str | None = None) -> dict:
    """Population AF per event by reciprocal overlap against a gnomAD-SV VCF.

    Delegates to the project's existing annotator when available (it already
    implements the reciprocal-overlap join and the local-copy discovery); returns
    an empty mapping when the resource is absent, so the run continues with the
    PON alone and says so in `privacy_note`.
    """
    if not gnomad_vcf or not os.path.exists(gnomad_vcf) or events.empty:
        return {}

    if helper and os.path.exists(helper):
        try:
            subprocess.run(["python3", helper, "--gnomad", gnomad_vcf],
                           check=True, capture_output=True)
        except subprocess.CalledProcessError as error:
            print(f"[confidence] gnomAD helper failed ({error}); continuing with "
                  f"PON only")
            return {}

    return _overlap_join(events, gnomad_vcf)


def _overlap_join(events: pd.DataFrame, gnomad_vcf: str) -> dict:
    """Minimal reciprocal-overlap join, used when no helper is configured.

    Only intra-chromosomal events with a resolvable size can be matched this way;
    inter-chromosomal junctions are left unmatched (not zero) because gnomAD-SV
    does not represent them comparably.
    """
    import gzip
    import re

    wanted = events[(events["chrom1"].astype(str) == events["chrom2"].astype(str))
                    & events["event_size"].notna()]
    if wanted.empty:
        return {}

    # Intervals sorted per chromosome so each gnomAD record can binary-search the
    # events that could overlap it, instead of scanning all of them. With a few
    # dozen events the difference is irrelevant; with thousands it is the
    # difference between minutes and hours.
    import bisect
    by_chrom: dict[str, list[tuple]] = {}
    for _, row in wanted.iterrows():
        start = min(int(row["pos1"]), int(row["pos2"]))
        end = start + int(float(row["event_size"]))
        by_chrom.setdefault(str(row["chrom1"]).replace("chr", ""), []).append(
            (start, end, str(row["sample_sv_id"])))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    starts = {c: [s for s, _, _ in v] for c, v in by_chrom.items()}
    max_len = {c: max((e - s for s, e, _ in v), default=0) for c, v in by_chrom.items()}

    hits: dict[str, float] = {}
    opener = gzip.open if gnomad_vcf.endswith(".gz") else open
    with opener(gnomad_vcf, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            chrom = f[0].replace("chr", "")
            if chrom not in by_chrom:
                continue
            pos = int(f[1])
            info = f[7]
            end_match = re.search(r"(?:^|;)END=(\d+)", info)
            af_match = re.search(r"(?:^|;)AF=([0-9.eE+-]+)", info)
            if not end_match or not af_match:
                continue
            g_start, g_end = pos, int(end_match.group(1))
            # Only events whose start lies within [g_start - longest, g_end] can
            # overlap this record.
            lo = bisect.bisect_left(starts[chrom], g_start - max_len[chrom])
            hi = bisect.bisect_right(starts[chrom], g_end)
            for start, end, sv_id in by_chrom[chrom][lo:hi]:
                overlap = min(end, g_end) - max(start, g_start)
                if overlap <= 0:
                    continue
                if overlap / max(end - start, 1) >= criteria.GNOMAD_RECIPROCAL_OVERLAP \
                        and overlap / max(g_end - g_start, 1) >= criteria.GNOMAD_RECIPROCAL_OVERLAP:
                    # Per-ancestry frequencies are parsed only for records that
                    # actually overlap, so the extra regexes cost nothing on the
                    # ~99.9% of the file that does not match.
                    record = {"gnomad_af": float(af_match.group(1))}
                    for population in criteria.GNOMAD_POPULATIONS:
                        match = re.search(rf"(?:^|;)AF_{population}=([0-9.eE+-]+)", info)
                        if match:
                            record[f"gnomad_af_{population}"] = float(match.group(1))
                    per_pop = [v for k, v in record.items() if k != "gnomad_af"]
                    if per_pop:
                        record["gnomad_af_popmax"] = max(per_pop)
                        record["gnomad_af_popmax_pop"] = max(
                            (k for k in record if k.startswith("gnomad_af_")
                             and k not in ("gnomad_af_popmax", "gnomad_af_popmax_pop")),
                            key=lambda k: record[k]).replace("gnomad_af_", "")
                    # Keep the record with the highest global AF if an event
                    # overlaps several gnomAD entries.
                    if sv_id not in hits or record["gnomad_af"] > hits[sv_id]["gnomad_af"]:
                        hits[sv_id] = record
    return hits
