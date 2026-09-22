"""
generators — stage 2, the only pluggable stage.

    config:  peptide_generator: neosv_trace | precomputed

Stage 2 is the sole point where the pipeline depends on external peptide-
generation code. Everything else — admission, the three-level cross, sequence QC,
event deduplication, SV confidence, privacy, RNA junction evidence, the null
model, attribution — is generator-agnostic and never imports a tool.

Registering a new backend:

    1. add `mytool.py` here, subclassing `PeptideGenerator`
    2. satisfy the contract in `base.py` (`validate_output` enforces it)
    3. add it to REGISTRY below

Nothing else in the codebase changes.

Available backends
------------------
`neosv`         Published NeoSV (MIT) + this repo's candidate-peptide output
                (default). One documented frame patch; see neosv.py.
`neosv_trace`   The NeoSV-Trace fork, vendored. Produces byte-identical peptides
                to `neosv`; kept for reproducing runs made with it.
`precomputed`   Reuse a stored output — reproduce a run without the tool, or
                bring in peptides generated elsewhere.
"""
from __future__ import annotations

from .base import (GenerationResult, GeneratorError, PeptideGenerator,
                   REQUIRED_ANNOTATION_COLUMNS, REQUIRED_PEPTIDE_COLUMNS)
from .neosv import NeoSVGenerator
from .neosv_trace import NeoSVTraceGenerator
from .precomputed import PrecomputedGenerator

REGISTRY = {
    NeoSVGenerator.name: NeoSVGenerator,
    NeoSVTraceGenerator.name: NeoSVTraceGenerator,
    PrecomputedGenerator.name: PrecomputedGenerator,
}

DEFAULT_GENERATOR = NeoSVGenerator.name


def get_generator(name: str | None = None, **kwargs) -> PeptideGenerator:
    """Instantiate a backend by name. Unknown names fail loudly with the list of
    valid ones — a typo must not silently fall back to a default."""
    key = (name or DEFAULT_GENERATOR).strip().lower()
    if key not in REGISTRY:
        raise GeneratorError(
            f"unknown peptide generator '{name}'. Available: "
            f"{sorted(REGISTRY)}. Add a backend in src/svneo/generators/ and "
            f"register it in REGISTRY.")
    return REGISTRY[key](**kwargs)


__all__ = ["get_generator", "REGISTRY", "DEFAULT_GENERATOR", "PeptideGenerator",
           "NeoSVGenerator", "NeoSVTraceGenerator",
           "GenerationResult", "GeneratorError", "REQUIRED_PEPTIDE_COLUMNS",
           "REQUIRED_ANNOTATION_COLUMNS"]
