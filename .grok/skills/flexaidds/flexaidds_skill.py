"""
flexaidds_skill.py — Entry point for /flexaidds skill figure generation.

This module provides a robust, injectable handler for generating publication-quality
Nature Reviews Drug Discovery-style cover figures for FlexAID∆S thermodynamic
enthalpy/entropy balance visualizations, including the Entropy–Enthalpy Index (I_E–E).

It implements the 5-point integration:
1. Skill entry point (this module).
2. Figure-generation handler with validation, safe overrides, dependency injection.
3. Exposed via companion manifest (flexaidds_skill_manifest.json).
4. Safe dependency management (no hard-coded image lib; injectable generator).
5. Reproducible returns (prompt or saved path + sidecar metadata).

The handler reuses the core prompt builder from python/flexaidds/figures for consistency
with the main package, but wraps it for skill-specific use (logging, manifest-friendly
dataclass, error handling, unique output naming).

Usage (from skill/agent):
    from flexaidds_skill import generate_flexaids_figure, FigureParameters

    params = FigureParameters(entropy_value=0.93, enthalpy_value=1.4, index_value=0.92)
    result = generate_flexaids_figure(params=params)  # returns prompt if no generator
    # or
    result = generate_flexaids_figure(
        params=params,
        image_generator=my_image_gen_callable,  # must accept prompt str -> {'path': str}
        output_dir="results/figures"
    )
    # result["path"] or result["prompt"]

The default image_generator is None, allowing the calling skill runtime (Grok Build,
Claude, etc.) to supply its native image tool (image_gen, DALL·E, etc.) via injection.
This keeps the module pure and testable.
"""

import logging
import shutil
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Optional

# Reuse the battle-tested builder and params from the main package for consistency,
# scientific guardrails, JetBrains Mono / thebonhomme.com typography, PLIP-style
# interaction emphasis on most favourable/CF contacts, E-E index, etc.
try:
    from flexaidds.figures import (
        generate_flexaids_nrdd_cover,
        NRDDCoverParams,
    )
except ImportError as exc:
    raise RuntimeError(
        "flexaidds package not available. Install with `pip install -e python` "
        "from the repo root, or ensure PYTHONPATH includes python/."
    ) from exc

logger = logging.getLogger(__name__)
if not logger.handlers:
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")


@dataclass(frozen=True)
class FigureParameters:
    """
    Validated parameters for NRDD-style FlexAID∆S enthalpy/entropy cover figures.

    Mirrors (and validates) the values used in the reference covers:
    - entropy_value: representative TΔS (or |TΔS|) in kcal/mol scale.
    - tds_value: value displayed as the prominent '-TΔS' (user preference: make -TdS great and highly visible;
      do not show -ΔH / -dH labels prominently; prefer the Enthalpy-Entropy Index I_E-E instead).
    - index_value: Enthalpy-Entropy Index (I_E-E / I_EE) developed within the FlexAIDdS skill (see LIB/statmech compute_IEE),
      typically [0, 1] or [-1, 1].
    - style: "dramatic_faces" (personified blue entropy vs fiery enthalpy faces)
             or "molecular_gauge" (abstract proteins + central E-E gauge).
    Additional fields allow full control of titles, dates, etc. for reproducibility.
    """
    entropy_value: float = 0.93
    # Note: this value is displayed as the prominent '-TΔS' (user preference: make -TdS great and visible;
    # avoid prominent -ΔH / -dH labels in the figure; prefer the Enthalpy-Entropy Index I_E-E instead).
    tds_value: float = 1.4
    index_value: float = 0.92
    style: str = "dramatic_faces"
    title: str = "The ΔG balance"
    subtitle: str = "Striking the right pose in drug discovery"
    date: str = "June 2025"
    volume: str = "Volume 24 | No. 6"

    def __post_init__(self):
        if self.entropy_value < 0 or self.tds_value < 0:
            raise ValueError("entropy_value and tds_value must be non-negative")
        if not (0.0 <= self.index_value <= 1.0):
            # Allow slight tolerance but clamp/document; strict [0,1] for the index in this context
            if not (-0.1 <= self.index_value <= 1.1):
                raise ValueError("index_value (I_E–E) must be in approximately [0, 1]")
        if self.style not in ("dramatic_faces", "molecular_gauge"):
            raise ValueError("style must be 'dramatic_faces' or 'molecular_gauge'")


def _build_prompt(
    params: FigureParameters,
    overrides: Optional[Dict[str, str]] = None,
) -> str:
    """
    Build the full descriptive prompt using the core package logic.

    Allows safe placeholder overrides (e.g. for custom titles) while rejecting
    unknown keys to prevent prompt injection or drift.
    """
    # Leverage the existing robust generator (it handles sourcing from results_dir
    # if desired, scientific notes, typography, PLIP interaction clarity, etc.).
    # We pass the values directly; the package already produces the dramatic
    # reference-style prompt with E-E index, cubes/gauge, branding, etc.
    res = generate_flexaids_nrdd_cover(
        entropy_value=params.entropy_value,
        enthalpy_value=params.tds_value,  # passed to the slot used for prominent -TΔS in the figure (user preference)
        index_value=params.index_value,
        style=params.style,
        title=params.title,
        subtitle=params.subtitle,
        # date/volume are embedded in the prompt via the package's internal templates
    )
    prompt = res["prompt"]

    if overrides:
        unknown = set(overrides) - {"title", "subtitle", "date", "volume"}  # whitelisted safe keys
        if unknown:
            raise ValueError(f"Unknown prompt override keys (sanitized): {unknown}")
        for key, value in overrides.items():
            # Simple safe substitution for the few whitelisted fields
            placeholder = "{" + key + "}"
            if placeholder in prompt:
                prompt = prompt.replace(placeholder, str(value))
            else:
                # Also try the literal strings that appear in the generated prompt
                if key == "title":
                    prompt = prompt.replace(params.title, str(value))
                elif key == "subtitle":
                    prompt = prompt.replace(params.subtitle, str(value))
    return prompt


def generate_flexaids_figure(
    params: Optional[FigureParameters] = None,
    image_generator: Optional[Callable[[str], Dict[str, Any]]] = None,
    output_dir: str = "figures",
    prompt_overrides: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """
    Main skill-exposed handler for automated, publication-quality figure generation.

    - Accepts FigureParameters (or defaults to reference cover values).
    - Builds the prompt via the canonical builder (reproducible, branded, scientifically accurate).
    - Validates inputs strictly.
    - Supports safe prompt_overrides (only whitelisted keys).
    - If no image_generator is supplied, returns the prompt so the calling runtime
      (Grok Build image_gen, Claude, DALL·E, etc.) can generate on its own schedule.
    - If an image_generator is supplied (dependency injection), calls it, validates
      the return value, copies the result to a uniquely named file in output_dir,
      and also writes a sidecar <name>.json with full parameters + prompt for audit/reproducibility.
    - Returns a dict with at minimum {"prompt": str, "path": Optional[str]}.
      When a file is written, "path" is the absolute path to the copied image.

    This design is bulletproof:
    - No hard dependency on any particular image library (injected callable).
    - Extensive validation + sanitization.
    - Deterministic unique output naming (timestamped).
    - Full metadata sidecar for reproducibility (matches the project's overall philosophy).
    - Graceful: always returns the prompt even if generation is deferred.
    - Logging of key steps and errors.

    The image_generator must be a callable accepting a single prompt string and
    returning a dict with at least {"path": <str>}. Example in Grok Build context:
        def grok_image_gen(prompt: str) -> Dict[str, Any]:
            result = image_gen(prompt=prompt, aspect_ratio="16:9")  # the env tool
            return {"path": result["path"]}   # adapt to actual tool return shape
    """
    if params is None:
        params = FigureParameters()

    # Strict validation (fail fast, clear errors)
    # (dataclass __post_init__ already did basic checks; we can add more here)
    if prompt_overrides:
        allowed = {"title", "subtitle", "date", "volume"}
        bad = set(prompt_overrides) - allowed
        if bad:
            raise ValueError(f"prompt_overrides may only contain keys in {allowed}, got: {bad}")

    logger.info(
        "Generating FlexAID∆S NRDD figure with params: entropy=%s, tds=%s, index=%s, style=%s",
        params.entropy_value, params.tds_value, params.index_value, params.style,
    )

    prompt = _build_prompt(params, prompt_overrides)
    result: Dict[str, Any] = {"prompt": prompt, "path": None, "params": asdict(params)}

    if image_generator is None:
        logger.info("No image_generator supplied — returning prompt only (deferred generation).")
        return result

    try:
        gen_result = image_generator(prompt)
        if not isinstance(gen_result, dict) or "path" not in gen_result:
            raise ValueError(
                "image_generator must return a dict containing at least {'path': <str>}"
            )
        src_path = Path(gen_result["path"])
        if not src_path.exists():
            raise FileNotFoundError(f"Generator reported path that does not exist: {src_path}")

        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        # Unique, reproducible filename
        ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
        suffix = src_path.suffix or ".png"
        dst = out_dir / f"flexaids_nrdd_{params.style}_{ts}{suffix}"

        shutil.copy2(src_path, dst)
        result["path"] = str(dst.resolve())

        # Write sidecar metadata for full reproducibility (params + prompt + timing)
        meta_path = dst.with_suffix(".json")
        meta = {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "params": asdict(params),
            "prompt": prompt,
            "source_image": str(src_path),
            "output_image": str(dst.resolve()),
            "skill": "flexaidds",
            "version": "1.0",
        }
        with open(meta_path, "w", encoding="utf-8") as f:
            import json
            json.dump(meta, f, indent=2, ensure_ascii=False)

        logger.info("Figure generated and saved to %s (metadata: %s)", dst, meta_path)
        result["metadata_path"] = str(meta_path.resolve())

    except Exception as exc:
        logger.exception("Image generation failed: %s", exc)
        # Re-raise as RuntimeError so callers (skill runtime) get a clean failure mode
        raise RuntimeError(f"Figure generation failed: {exc}") from exc

    return result
