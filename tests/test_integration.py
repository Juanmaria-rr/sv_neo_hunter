#!/usr/bin/env python3
"""
test_integration.py — does a fresh install reproduce known results?

WHY THIS EXISTS
---------------
The real inputs of this analysis are not public: the reference catalogue is
access-controlled and the sample VCFs and BAMs belong to the group that produced
them. Anyone cloning this repository therefore cannot reproduce our figures, and
has no way to tell a correct installation from a broken one.

This runs the whole chain on a small synthetic dataset shipped in `tests/data/`
and asserts exact values. If it passes, the installation computes what ours
computes. If it fails, something in the environment differs — before any real
data is involved.

WHAT THE FIXTURE EXERCISES
--------------------------
`mini.vcf` is ten breakends chosen to hit the cases that have actually broken
this pipeline:

  sv1  221 bp deletion               the ordinary case, above the gap-test floor
  sv2  3 bp span, 21 inserted bases  insertion-driven: the span is not the lesion
  sv3  inter-chromosomal             no defined size; needs the chimeric test
  sv4  PON_COUNT 2,500               must be rejected by the panel filter
  sv5  shares position 4000015 with sv4b — a rejected record must not be written
       out because its neighbour passed

The peptide fixture has four peptides, and each is removed by a different gate:

  MATCHPEPT  in the catalogue, gene-concordant, high-complexity  -> survives
  NOTINCAT1  absent from the catalogue                           -> never matches
  FFFFFFFFF  in the catalogue but low-complexity                 -> dropped at QC
  OTHERPEPT  in the catalogue, annotated to a different gene     -> discordant

Stage 2 uses the `precomputed` backend, so this test needs neither pyensembl nor
the vendored generator: it verifies OUR logic, and `tools/compare_generators.py`
verifies the generator separately.

    python tests/test_integration.py
"""
from __future__ import annotations

import os
import pathlib
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import pandas as pd  # noqa: E402

from svneo import confidence, criteria, cross, vcf  # noqa: E402

DATA = pathlib.Path(__file__).parent / "data"
FAILURES: list[str] = []


def check(statement: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {statement}" + (f"  ({detail})" if detail and not condition else ""))
    if not condition:
        FAILURES.append(statement)


def test_stage1_admission():
    """Ten breakends -> the panel filter rejects one junction, leaving four."""
    breakends = vcf.read_breakends(str(DATA / "mini.vcf"))
    check("stage 1 reads all 10 breakends", len(breakends) == 10, f"got {len(breakends)}")

    admitted, funnel = vcf.admit(breakends, "germline", criteria.PON_MAX)
    # sv4a/sv4b carry PON_COUNT 2,500 and must go; the other eight remain.
    check("panel filter rejects the PON_COUNT 2,500 pair",
          len(admitted) == 8, f"admitted {len(admitted)}")

    junctions = vcf.pair_junctions(admitted)
    check("four junctions after pairing", len(junctions) == 4, f"got {len(junctions)}")

    sizes = dict(zip(junctions.sample_sv_id, junctions.event_size))
    check("a 221 bp span is a 220 bp deletion", sizes.get("sv1a") == 220,
          f"got {sizes.get('sv1a')}")
    inserts = dict(zip(junctions.sample_sv_id, junctions.insert_len))
    check("the insertion-driven event keeps its 21 inserted bases",
          inserts.get("sv2a") == 21, f"got {inserts.get('sv2a')}")
    check("the inter-chromosomal junction has no defined size",
          pd.isna(sizes.get("sv3a")), f"got {sizes.get('sv3a')}")
    return junctions, admitted


def test_admitted_vcf_excludes_rejected_neighbour(admitted):
    """sv5a sits at the same position as the rejected sv4b."""
    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as handle:
        out = handle.name
    written = vcf.write_admitted_vcf(str(DATA / "mini.vcf"), admitted, out)
    text = pathlib.Path(out).read_text()
    check("the admitted VCF holds exactly the admitted records", written == 8,
          f"wrote {written}")
    check("a rejected breakend is not written out via a shared coordinate",
          "sv4a" not in text and "sv4b" not in text)
    check("its position-sharing neighbour IS written", "sv5a" in text)
    os.unlink(out)


def test_cross_and_qc(junctions):
    """One catalogue match survives QC; a low-complexity one does not."""
    peptides = pd.read_csv(DATA / "mini.all_neopeptides.txt", sep="\t", dtype=str)
    anno = pd.read_csv(DATA / "mini.anno.txt", sep="\t", dtype=str)

    class Reference:
        peptides = str(DATA / "mini_catalogue.tsv")
        peptide_column, gene_column, svtype_column = "neoantigen", "gene1", "svtype"

    reference = cross.load_reference(Reference)
    matches = cross.level1_identical(peptides, anno, reference)
    # Three of the four fixture peptides are in the catalogue; the fourth is not.
    check("three peptides match the catalogue", len(matches) == 3, f"got {len(matches)}")
    # Each of the three is removed by a DIFFERENT gate, which is the point of the
    # fixture: one survives, one is low-complexity, one is gene-discordant.
    discordant = matches[~matches.gene_concordant.fillna(False)]
    check("the catalogue peptide annotated to another gene is discordant",
          set(discordant.peptide) == {"OTHERPEPT"}, f"got {set(discordant.peptide)}")

    qc = cross.sequence_qc(matches, proteome=None)
    check("the low-complexity peptide is flagged",
          bool(qc[qc.peptide == "FFFFFFFFF"].low_complexity.iloc[0]))
    credible = qc[qc.credible.fillna(False)]
    check("exactly one match is credible", len(credible) == 1, f"got {len(credible)}")
    check("and it is the expected peptide",
          set(credible.peptide) == {"MATCHPEPT"}, f"got {set(credible.peptide)}")

    events = cross.to_events(qc, junctions)
    check("one credible genomic event", len(events) == 1, f"got {len(events)}")
    return events


def test_confidence_and_privacy(events):
    """No gnomAD resource: privacy rests on the panel, and says so."""
    annotated = confidence.annotate_confidence(events)
    check("the surviving event is high-confidence",
          bool(annotated.sv_hc.iloc[0]))

    private = confidence.annotate_privacy(annotated, panel_size=3000, gnomad_af={})
    check("an event with no PON_COUNT counts as private",
          bool(private.is_private.iloc[0]))
    check("and the run records that gnomAD was not evaluated",
          "not evaluated" in str(private.privacy_note.iloc[0]))


def test_rna_test_selection():
    """Geometry decides the test; no BAM needed to check the decision."""
    from svneo import rna
    cases = [("DEL", "1", "1", 220, 0, "sizegap"),
             ("DEL", "1", "1", 3, 21, "insertion"),
             ("BND", "1", "7", None, 0, "chimeric"),
             ("DEL", "1", "1", 3, 0, "none")]
    for svtype, c1, c2, size, insert, expected in cases:
        got, _ = rna.select_test(svtype, c1, c2, size, insert)
        check(f"{svtype} size={size} insert={insert} -> {expected}",
              got == expected, f"got {got}")


def main() -> None:
    print("integration test on the synthetic fixture\n")
    junctions, admitted = test_stage1_admission()
    test_admitted_vcf_excludes_rejected_neighbour(admitted)
    events = test_cross_and_qc(junctions)
    test_confidence_and_privacy(events)
    test_rna_test_selection()

    print(f"\n{len(FAILURES)} failure(s)")
    if FAILURES:
        print("\nA failure here means this installation does not compute what the")
        print("reference one does. Investigate before running real data:")
        for failure in FAILURES:
            print(f"  - {failure}")
        sys.exit(1)
    print("This installation reproduces the reference results on the fixture.")


if __name__ == "__main__":
    main()
