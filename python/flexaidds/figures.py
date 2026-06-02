"""Publication figure & animated cover generation for FlexAID∆S (Grok Imagine integration).

This module is **purely post-processing and promotional**. It never reads or writes
scientific kernels, never changes scores/ranking/poses, and never claims thermodynamic
values are experimental ∆G.

It consumes only the public `load_results()` + reproducibility/audit JSONs already
produced by a successful docking run, selects the best-scoring binding mode (by the
ensemble free energy ledger), and emits ready-to-use prompts + metadata for
imagine_text_to_image / imagine_image_to_video / image_edit.

Gate: figure preparation is intended to run only after Gate 6 (F/S cross-check) has
passed in the thermodynamic audit. The prepare function enforces this by default.

Aesthetic contract (redesigned for NRDD cover + high-end MD viz reference):
- Deep navy gradients (#0a0e14 family) + cyan/teal #22D3EE accents (from site branding).
- Gold #FBBF24 for ∆G, terra #A78BFA for entropy term.
- Exact bottom banner "/flexaids-docking • FlexAID∆S".
- Prominent baked-in equation ∆G=∆H−T∆S with real run values.
- Reproducibility footer (gate status, short git, date, run id).
- Cinematic scientific illustration quality suitable for Nature Reviews Drug Discovery
  cover and X scientific posts, inspired by elegant SwitchCraft-style MD visualizations
  (clean hybrid cartoon/surface, sophisticated lighting, subtle dynamics in animation).
- AI-tool compatible plain-text prompts (Grok, Claude, ChatGPT, etc.).

Primary entry for runs & skill agent flow:
    from flexaidds.figures import prepare_publication_figures
    prepare_publication_figures(Path("results/my_run"), visualize=True, require_gate6=True)

Output layout under <results_dir>/figures/:
    prompt_cover.txt
    prompt_animation.txt
    figure_metadata.json   (full snapshot + paths for audit)
    (optional) base_publication.png   (if use_pymol_base and PyMOL present)
    (later materialized by agent) cover_best_mode.png , animation_6s.mp4

All functions are non-destructive and best-effort.
"""

from __future__ import annotations

import json
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional

from .models import DockingResult
from .results import load_results

# Brand colors (exact match to site/colors_and_type.css and GitHub visual identity)
BRAND = {
    "bg_deep_navy": "#0a0e14",
    "teal": "#22D3EE",      # primary accent, ∆H, CTA
    "gold": "#FBBF24",      # ∆G highlight
    "terra": "#A78BFA",     # entropy / T∆S
    "fg": "#d4dced",
}

# Required banner and footer text (never change without coordinated doc update)
BANNER = "/flexaids-docking • FlexAID∆S"
FOOTER_TEMPLATE = "gate6:{gate} • git:{git} • {date} • {runid} • FlexAIDdS"

# The canonical long prompt template for NRDD-cover + ref MD aesthetic.
# Real values are injected at runtime via .format(). Keep placeholders stable.
TEMPLATE_COVER = """Create a stunning, publication-ready scientific figure suitable for the cover of Nature Reviews Drug Discovery. Cinematic, high-end molecular visualization in the refined aesthetic of premium 2026 protein design molecular dynamics videos (e.g. SwitchCraft-style elegant animations): the best-scoring FlexAID∆S binding mode for {ligand} in {receptor}, shown with exquisite detail.

Hybrid rendering — translucent molecular surface + cartoon ribbons for the receptor protein in deep navy tones ({bg_deep_navy}), ligand rendered in bright cyan/teal sticks/balls ({teal}) with crisp atomic detail and subtle glow, 3-5 key induced-fit side chains highlighted in matching teal with thin H-bond dashes. Subtle blue-to-red entropy heatmap wash on the receptor surface and flexible loops (blue = low configurational entropy/rigid, red = high entropy/flexible regions per the run's Shannon + tENCoM values).

Prominently and elegantly overlay the thermodynamic equation in clean modern typography: ΔG = ΔH − TΔS   with actual values ΔG = {delta_g:.2f} kcal/mol , ΔH = {delta_h:.2f} , −TΔS = {minus_t_ds:.2f}  at T={temperature:.2f} K (gold for ΔG, teal for ΔH, purple for entropy term).

Composition: centered binding interface at a slight 3D angle for depth, soft cinematic volumetric lighting with cyan rim lights, deep navy to near-black gradient background with extremely faint abstract molecular field lines or density, shallow depth-of-field, ultra-clean professional finish, 8K scientific illustration quality, no clutter.

At bottom: sleek integrated banner bar with exact text '{banner}' in refined sans-serif, small subtle icons for code + AI on sides, FlexAID∆S wordmark. Bottom-right or footer in tiny crisp text: reproducibility metadata '{footer}'.

Overall mood: confident, precise, beautiful, suitable for top-tier journal cover and high-engagement X/Twitter scientific thread. Aspect for cover: 3:2 or 16:9 cinematic. High contrast, sharp text overlays baked in. AI-tool compatible plain text prompt."""

TEMPLATE_ANIMATION = """6-second seamless cinematic animation (1080p or 4K, 6s duration, loop-friendly) of the best-scoring FlexAID∆S binding mode for {ligand} in {receptor}.

Smooth slow 360° orbit + gentle dolly around the docked complex (exact same pose and induced-fit geometry as the static cover). Very subtle breathing motion on the 3-5 highlighted induced-fit side chains (low amplitude, physically plausible), ligand micro-fluctuations consistent with the ensemble. Flowing faint entropy color waves (blue low-entropy ↔ red high-entropy) pulsing gently across receptor surface and loops.

Equation block and bottom banner '{banner}' fade in elegantly at t≈1.5s and persist; values and reproducibility footer '{footer}' appear cleanly without jitter. Camera motion, timing, and production quality exactly like the referenced high-end molecular dynamics visualization videos (SwitchCraft aesthetic): fluid, sophisticated lighting, clean dynamic protein representations, no clutter, premium scientific art.

Deep navy gradient background, cyan/teal rim highlights ({teal}), gold accents on ΔG. Ultra-smooth, high dynamic range, 8K detail feel. AI-tool compatible prompt for video generation."""


def _safe_float(v: Any, default: float = 0.0) -> float:
    try:
        if v is None:
            return default
        return float(v)
    except (TypeError, ValueError):
        return default


def check_gate6_passed(results_dir: Path) -> bool:
    """Scan common JSON artifacts under results_dir for a passing Gate 6.

    Looks for:
      - reproducibility.json (provenance.gate_results.gate6_crosscheck.passed)
      - any *_audit*.json or thermo*.json with the same shape
      - direct "gate6_crosscheck": {"passed": true} at top level (legacy)

    Returns False on any error or missing data (conservative).
    """
    results_dir = Path(results_dir)
    candidates = list(results_dir.rglob("*.json"))
    for p in candidates:
        if p.name.endswith((".tar.gz", ".zip")):
            continue
        try:
            data = json.loads(p.read_text(encoding="utf-8", errors="replace"))
        except Exception:
            continue

        # Direct legacy shape
        g6 = data.get("gate6_crosscheck")
        if isinstance(g6, dict) and g6.get("passed") is True:
            return True

        # Provenance shape (thermo audit + run_metadata)
        prov = data.get("provenance") or data.get("gate_results") or {}
        if isinstance(prov, dict):
            gr = prov.get("gate_results") or prov
            g6 = gr.get("gate6_crosscheck") or gr.get("gate6")
            if isinstance(g6, dict) and g6.get("passed") is True:
                return True
            # also accept flat "gate6_crosscheck": {"passed": true} under provenance
            if gr.get("gate6_crosscheck", {}).get("passed") is True:
                return True
    return False


def extract_best_mode_summary(results_dir: Path) -> Dict[str, Any]:
    """Produce a compact dict with everything the prompt builders need.

    Prefers detailed thermodynamics block when present; falls back to scalar
    fields on BindingModeResult. Always includes a best_pose PDB path when
    discoverable.
    """
    results_dir = Path(results_dir).resolve()
    summary: Dict[str, Any] = {
        "results_dir": str(results_dir),
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "gate6_passed": check_gate6_passed(results_dir),
        "ligand": "ligand",
        "receptor": "receptor",
        "best_mode_id": 1,
        "n_poses": 0,
        "delta_g": None,
        "delta_h": None,
        "minus_t_ds": None,
        "temperature": 298.15,
        "git_sha": "unknown",
        "run_id": results_dir.name,
        "best_pose_pdb": None,
    }

    # 1. Try structured load_results (best source)
    try:
        docking: DockingResult = load_results(results_dir)
        if docking.binding_modes:
            # top by free_energy (most negative first, as produced by package)
            modes = sorted(docking.binding_modes, key=lambda m: (m.free_energy or 1e9))
            best = modes[0]
            summary["best_mode_id"] = getattr(best, "mode_id", 1) or getattr(best, "rank", 1) or 1
            summary["n_poses"] = best.n_poses
            summary["temperature"] = best.temperature or summary["temperature"]

            # Prefer detailed thermo if present
            thermo = getattr(best, "thermodynamics", None) or {}
            if isinstance(thermo, dict):
                g = thermo.get("G_total_kcal_mol") or thermo.get("free_energy")
                h = thermo.get("H_eff_kcal_mol") or thermo.get("enthalpy")
                s = thermo.get("minus_T_S_config_kcal_mol") or thermo.get("entropy")
            else:
                g = getattr(best, "free_energy", None)
                h = getattr(best, "enthalpy", None)
                s = getattr(best, "entropy", None)

            summary["delta_g"] = _safe_float(g, 0.0)
            summary["delta_h"] = _safe_float(h, 0.0)
            # convention in ledger: entropy field is often already -TΔS or S; we want the -TΔS term for the equation
            if s is not None and _safe_float(s, 0.0) > 0 and (g is not None and h is not None):
                # rough: if positive and G/H present, treat as -TΔS magnitude
                summary["minus_t_ds"] = -_safe_float(s, 0.0)
            else:
                summary["minus_t_ds"] = _safe_float(s, summary["delta_g"] - summary["delta_h"] if summary["delta_g"] and summary["delta_h"] else 0.0)

            # best pose file
            bp = best.best_pose()
            if bp and getattr(bp, "path", None):
                summary["best_pose_pdb"] = str(bp.path)

            # metadata hints
            md = getattr(best, "metadata", {}) or {}
            if md.get("ligand_name"):
                summary["ligand"] = str(md["ligand_name"])
            if md.get("receptor_name"):
                summary["receptor"] = str(md["receptor_name"])
    except Exception as exc:
        warnings.warn(f"load_results failed for figures: {exc}. Falling back to filename heuristics.", RuntimeWarning)

    # 2. Fallback: scan reproducibility.json for provenance + numbers
    try:
        for j in results_dir.rglob("reproducibility*.json"):
            data = json.loads(j.read_text(errors="replace"))
            prov = data.get("provenance", {})
            if prov.get("flexaidds_git", {}).get("commit"):
                summary["git_sha"] = prov["flexaidds_git"]["commit"][:7]
            # try to pull top mode thermo if present at top level
            if data.get("top_mode"):
                tm = data["top_mode"]
                summary["delta_g"] = _safe_float(tm.get("free_energy") or summary["delta_g"])
                summary["delta_h"] = _safe_float(tm.get("enthalpy") or summary["delta_h"])
                summary["minus_t_ds"] = _safe_float(tm.get("entropy") or summary["minus_t_ds"])
            break
    except Exception:
        pass

    # 3. Last-resort: derive ligand/receptor from dir name or input hints
    if summary["ligand"] == "ligand":
        for p in results_dir.rglob("*.mol2"):
            summary["ligand"] = p.stem
            break
    if summary["receptor"] == "receptor":
        for p in results_dir.rglob("*.pdb"):
            if "receptor" in p.name.lower() or "protein" in p.name.lower():
                summary["receptor"] = p.stem
                break

    # 4. git sha fallback from any .git in ancestor (best effort)
    if summary["git_sha"] == "unknown":
        try:
            from subprocess import run, PIPE
            r = run(["git", "rev-parse", "--short", "HEAD"], cwd=results_dir, capture_output=True, text=True, timeout=2)
            if r.returncode == 0:
                summary["git_sha"] = r.stdout.strip()
        except Exception:
            pass

    return summary


def build_imagine_cover_prompt(summary: Dict[str, Any], *, style: str = "nrdd-cover") -> str:
    """Return the full engineered static cover prompt with real values injected."""
    s = summary
    gate = "PASS" if s.get("gate6_passed") else "UNKNOWN"
    footer = FOOTER_TEMPLATE.format(
        gate=gate,
        git=s.get("git_sha", "unknown")[:7],
        date=s.get("timestamp", "")[:10],
        runid=s.get("run_id", "run"),
    )
    vals = {
        "ligand": s.get("ligand", "ligand"),
        "receptor": s.get("receptor", "receptor"),
        "delta_g": _safe_float(s.get("delta_g"), 0.0),
        "delta_h": _safe_float(s.get("delta_h"), 0.0),
        "minus_t_ds": _safe_float(s.get("minus_t_ds"), 0.0),
        "temperature": _safe_float(s.get("temperature"), 298.15),
        "banner": BANNER,
        "footer": footer,
        **BRAND,
    }
    # style is reserved for future variants; current template is the NRDD one
    return TEMPLATE_COVER.format(**vals)


def build_imagine_animation_prompt(summary: Dict[str, Any], duration_s: float = 6.0) -> str:
    """Return the 6 s animation prompt (motion + same branding)."""
    s = summary
    gate = "PASS" if s.get("gate6_passed") else "UNKNOWN"
    footer = FOOTER_TEMPLATE.format(
        gate=gate,
        git=s.get("git_sha", "unknown")[:7],
        date=s.get("timestamp", "")[:10],
        runid=s.get("run_id", "run"),
    )
    vals = {
        "ligand": s.get("ligand", "ligand"),
        "receptor": s.get("receptor", "receptor"),
        "banner": BANNER,
        "footer": footer,
        "teal": BRAND["teal"],
    }
    base = TEMPLATE_ANIMATION.format(**vals)
    # Inject duration if the template ever uses it (kept for forward compat)
    return base.replace("6-second", f"{duration_s:g}-second")


def prepare_publication_figures(
    results_dir: Path,
    *,
    visualize: bool = True,
    require_gate6: bool = True,
    use_pymol_base: bool = False,
    force: bool = False,
) -> Dict[str, Any]:
    """Prepare prompts + metadata for imagine cover + animation.

    This is the function called by run_flexaidds.sh --visualize and by the
    flexaid-docking skill agent after a Gate-6-successful docking.

    Returns a dict with paths and status. Creates <results_dir>/figures/ only
    when it decides to proceed.
    """
    results_dir = Path(results_dir).resolve()
    out: Dict[str, Any] = {"results_dir": str(results_dir), "figures_dir": None, "proceeded": False}

    if not visualize:
        return out

    gate_ok = check_gate6_passed(results_dir)
    out["gate6_passed"] = gate_ok

    if require_gate6 and not gate_ok and not force:
        out["skipped"] = "gate6"
        warnings.warn(
            "prepare_publication_figures: Gate 6 not passed (or no audit JSON found) and require_gate6=True. "
            "Skipping figure prompt generation. Use force=True or --visualize with explicit confirmation.",
            RuntimeWarning,
        )
        return out

    summary = extract_best_mode_summary(results_dir)

    fig_dir = results_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    out["figures_dir"] = str(fig_dir)

    cover_prompt = build_imagine_cover_prompt(summary)
    anim_prompt = build_imagine_animation_prompt(summary)

    (fig_dir / "prompt_cover.txt").write_text(cover_prompt, encoding="utf-8")
    (fig_dir / "prompt_animation.txt").write_text(anim_prompt, encoding="utf-8")

    meta = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "summary": {k: v for k, v in summary.items() if k != "best_pose_pdb"},  # keep small
        "best_pose_pdb": summary.get("best_pose_pdb"),
        "branding": BRAND,
        "banner": BANNER,
        "required_elements": [
            "bottom banner: " + BANNER,
            "equation ΔG=ΔH−TΔS with injected values",
            "reproducibility footer with gate/git/date",
            "cyan/teal accents + deep navy gradients",
            "entropy blue→red heatmap description",
            "induced-fit side chains",
        ],
        "gate6_passed": gate_ok,
        "use_pymol_base_attempted": use_pymol_base,
    }
    (fig_dir / "figure_metadata.json").write_text(json.dumps(meta, indent=2, sort_keys=True), encoding="utf-8")

    # Optional PyMOL base render (best-effort, never fatal)
    base_path: Optional[str] = None
    if use_pymol_base:
        try:
            from .visualization import setup_publication_view  # type: ignore
            # We cannot easily drive full PyMOL headless here without a loaded complex.
            # For now we write a tiny helper note + a placeholder path convention.
            # Real base render is expected to be done by the caller (PyMOL plugin or
            # explicit `python -c 'import pymol; ...'` script that loads the best pose).
            note = (
                "To produce base_publication.png for image-to-image enhancement:\n"
                "  1. (in PyMOL) flexaids_load_results " + str(results_dir) + "\n"
                "  2. Select best mode, run setup_publication_view with ray_trace=True + output_png\n"
                "  3. Copy the PNG to " + str(fig_dir / "base_publication.png") + "\n"
                "Then re-run prepare with the base present; the imagine prompt will suggest image_to_image."
            )
            (fig_dir / "BASE_RENDER_INSTRUCTIONS.txt").write_text(note, encoding="utf-8")
            candidate = fig_dir / "base_publication.png"
            if candidate.exists():
                base_path = str(candidate)
        except Exception as e:
            warnings.warn(f"PyMOL base render path unavailable: {e}", RuntimeWarning)

    out.update({
        "proceeded": True,
        "cover_prompt_path": str(fig_dir / "prompt_cover.txt"),
        "animation_prompt_path": str(fig_dir / "prompt_animation.txt"),
        "metadata_path": str(fig_dir / "figure_metadata.json"),
        "base_image_path": base_path,
        "summary": summary,
    })
    return out
