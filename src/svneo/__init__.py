"""
svneo — recurrence testing for SV-derived neoantigens.

Asks one question, of any sample set against any reference catalogue:

    do these reference neopeptides reappear in these samples, and if so is the
    recurrence real (a shared genomic event) or an artefact (mapping noise,
    low-complexity sequence, coincidental peptide convergence, or a common
    germline polymorphism)?

The package is the merge of two independently built workflows that disagreed on
the same three events over the same data; docs/PROVENANCE.md records what each
contributed and which rule each disagreement produced.

Stage map (module -> stage):
    vcf         1     admission, breakend pairing, type-aware event size
    generators  2     candidate neopeptide generation (PLUGGABLE:
                      the only stage that touches external tooling;
                      NeoSV-Trace lives in generators/neosv_trace.py
                      and vendor/neosv_trace/)
    cross       3-5   three-level cross, sequence QC, collapse to events
    confidence  6     SV-call confidence (HC) and privacy (PON + population AF)
    rna         7-8   expression context and junction-crossing read evidence
    null_model  -     chance-match null, making samples comparable
    synthesis   9     funnel, lineage attribution, cross-sample tables
    run         -     orchestrator
    criteria    -     every threshold, in one place
"""

__version__ = "0.1.0"

__all__ = ["criteria", "config", "vcf", "generators", "cross", "confidence",
           "rna", "null_model", "synthesis", "run"]
