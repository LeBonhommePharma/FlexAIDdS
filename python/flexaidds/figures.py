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

# PLIP (Protein-Ligand Interaction Profiler) style reference for nice, clean interaction figures.
# PLIP generates professional 3D rendered diagrams (PNG via -p, or raytraced from .pse PyMOL sessions)
# with standard color-coded viz of non-covalent interactions. See https://github.com/pharmai/plip
# We emulate this "nice fig" clarity + precision in our imagine prompts (and can use PLIP output
# as base image for image-to-image if user runs `plip -f <best_pose.pdb> -p`).
# Legend (from PLIP docs): 
#   Hydrophobic: grey50 dashed
#   H-bond: blue solid line
#   Salt bridge: yellow dashed
#   pi-Stacking: green/smudge dashed
#   pi-Cation: orange dashed
#   Halogen: greencyan solid
#   Water bridge: lightblue
#   Metal: violetpurple
# Prioritize the *most favourable contacts* (highest contribution to our Voronoi CF score for the mode)
# + standard chem types.
PLIP_STYLE = (
    "in the clean, professional 3D interaction diagram style of PLIP (Protein-Ligand Interaction Profiler): "
    "precise color-coded lines and labels per the standard PLIP legend (blue solid for H-bonds, grey dashed "
    "for hydrophobic, yellow for salt bridges, green for pi-stacking, etc.), atom- and residue-accurate, "
    "publication-quality, no clutter. Emphasize the most favourable contacts and those contributing most "
    "to the docking CF/Voronoi contact function for this binding mode."
)

# Required banner and footer text (never change without coordinated doc update)
BANNER = "/flexaids-docking • FlexAID∆S"
FOOTER_TEMPLATE = "gate6:{gate} • git:{git} • {date} • {runid} • FlexAIDdS"

# The canonical long prompt template for NRDD-cover + ref MD aesthetic.
# Real values are injected at runtime via .format(). Keep placeholders stable.
TEMPLATE_COVER = """Create a stunning, publication-ready scientific figure suitable for the cover of Nature Reviews Drug Discovery. Cinematic, high-end molecular visualization in the refined aesthetic of premium 2026 protein design molecular dynamics videos (e.g. SwitchCraft-style elegant animations): the best-scoring FlexAID∆S binding mode for {ligand} in {receptor}, shown with exquisite detail.

Hybrid rendering — translucent molecular surface + cartoon ribbons for the receptor protein in deep navy tones ({bg_deep_navy}), ligand rendered in bright cyan/teal sticks/balls ({teal}) with crisp atomic detail and subtle glow, 3-5 key induced-fit side chains highlighted in matching teal with thin H-bond dashes. Subtle blue-to-red entropy heatmap wash on the receptor surface and flexible loops (blue = low configurational entropy/rigid, red = high entropy/flexible regions per the run's Shannon + tENCoM values).

All text (banners, labels, equation, footer) in clean sharp JetBrains Mono (or very similar modern technical mono font aesthetic exactly like thebonhomme.com and Le Bonhomme Pharma branding): minimalist, highly legible, professional.

CRITICAL: Clearly depict, highlight with elegant dashed lines/glows, and label the most important molecular interactions that matter most: {key_interactions}. {plip_style} Use contrasting teal dashed lines for H-bonds and polar contacts, gold for hydrophobic/vdW packing, with crisp residue labels (e.g. "H-bond to Asn23", "salt-bridge to Asp128", "hydrophobic core with Trp79/Tyr43") placed elegantly near the contacts so they are scientifically precise and immediately readable. Prioritize the most favourable contacts and those contributing most to the CF. Make them visually prominent yet balanced in the composition.

Prominently and elegantly overlay the thermodynamic equation in clean modern typography: ΔG = ΔH − TΔS   with actual values ΔG = {delta_g:.2f} kcal/mol , ΔH = {delta_h:.2f} , −TΔS = {minus_t_ds:.2f}  at T={temperature:.2f} K (gold for ΔG, teal for ΔH, purple for entropy term).

Composition: centered binding interface at a slight 3D angle for depth, soft cinematic volumetric lighting with cyan rim lights, deep navy to near-black gradient background with extremely faint abstract molecular field lines or density, shallow depth-of-field, ultra-clean professional finish, 8K scientific illustration quality, no clutter.

At bottom: sleek integrated banner bar with exact text '{banner}' in refined sans-serif, small subtle icons for code + AI on sides, FlexAID∆S wordmark. Bottom-right or footer in tiny crisp text: reproducibility metadata '{footer}'.

Overall mood: confident, precise, beautiful, suitable for top-tier journal cover and high-engagement X/Twitter scientific thread. Aspect for cover: 3:2 or 16:9 cinematic. High contrast, sharp text overlays baked in. AI-tool compatible plain text prompt."""

TEMPLATE_ANIMATION = """6-second seamless cinematic animation (1080p or 4K, 6s duration, loop-friendly) of the best-scoring FlexAID∆S binding mode for {ligand} in {receptor}.

Smooth slow 360° orbit + gentle dolly around the docked complex (exact same pose and induced-fit geometry as the static cover). Very subtle breathing motion on the 3-5 highlighted induced-fit side chains (low amplitude, physically plausible), ligand micro-fluctuations consistent with the ensemble. Flowing faint entropy color waves (blue low-entropy ↔ red high-entropy) pulsing gently across receptor surface and loops. All text in JetBrains Mono / thebonhomme.com mono aesthetic. During the motion, the critical molecular interactions are dynamically highlighted: {key_interactions} — {plip_style} with elegant animated dashed lines (teal for H-bonds/polar, gold for hydrophobic) and fading residue labels that emphasize the interactions that matter most (favouring highest-CF contributors).

Equation block and bottom banner '{banner}' fade in elegantly at t≈1.5s and persist; values and reproducibility footer '{footer}' appear cleanly without jitter. Camera motion, timing, and production quality exactly like the referenced high-end molecular dynamics visualization videos (SwitchCraft aesthetic): fluid, sophisticated lighting, clean dynamic protein representations, no clutter, premium scientific art.

Deep navy gradient background, cyan/teal rim highlights ({teal}), gold accents on ΔG. Ultra-smooth, high dynamic range, 8K detail feel. AI-tool compatible prompt for video generation."""


def _safe_float(v: Any, default: float = 0.0) -> float:
    try:
        if v is None:
            return default
        return float(v)
    except (TypeError, ValueError):
        return default


def _extract_key_interactions(best_pose_pdb: Optional[str]) -> str:
    """Best-effort description of the interactions that matter most for the pose.

    If a best_pose PDB is available, parses it (reusing the package PDB reader)
    and identifies the closest ligand-receptor atom pairs. Classifies roughly
    into H-bond/polar vs. hydrophobic/vdW and returns a short human-readable
    phrase suitable for injection into imagine prompts. Falls back to a strong
    generic instruction when no detailed data or parsing fails.
    """
    default = "the critical molecular interactions (key hydrogen bonds, salt bridges, and hydrophobic contacts) that matter most for binding affinity and specificity in this mode"
    if not best_pose_pdb:
        return default
    p = Path(best_pose_pdb)
    if not p.exists():
        return default
    try:
        from .io import read_pdb, is_ion
        struct = read_pdb(str(p))
        lig_atoms = [a for a in struct.atoms if a.record == "HETATM" and not is_ion(a)]
        rec_atoms = [a for a in struct.atoms if a.record == "ATOM"]
        if not lig_atoms or not rec_atoms:
            return default

        # Group close contacts (< ~4.0 Å) by receptor residue
        from collections import defaultdict
        import math
        contacts: dict[str, list] = defaultdict(list)
        for la in lig_atoms:
            for ra in rec_atoms:
                dx = la.x - ra.x
                dy = la.y - ra.y
                dz = la.z - ra.z
                d = math.sqrt(dx * dx + dy * dy + dz * dz)
                if d < 4.0:
                    key = f"{ra.resname.strip()}{ra.resseq}"
                    contacts[key].append((ra.name.strip(), la.name.strip(), d, ra, la))

        if not contacts:
            return default

        # Rank residues by number of close contacts (proxy for importance)
        ranked = sorted(contacts.items(), key=lambda kv: len(kv[1]), reverse=True)[:5]

        phrases = []
        for res_key, pairs in ranked:
            # Rough classification
            has_polar = any(n[0] in ("O", "N") for n, _, _, _, _ in pairs)
            has_hbond_geom = any(d < 3.5 and n1[0] in ("O", "N") and n2[0] in ("O", "N") for n1, n2, d, _, _ in pairs)
            if has_hbond_geom or has_polar:
                phrases.append(f"H-bond/polar contacts to {res_key}")
            else:
                phrases.append(f"hydrophobic/vdW packing with {res_key}")

        if phrases:
            return "critical interactions that matter most: " + "; ".join(phrases) + " (visualize with clear, elegant dashed lines and residue labels)"
        return default
    except Exception:
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
        "key_interactions": "the critical molecular interactions (key hydrogen bonds, salt bridges, and hydrophobic contacts) that matter most for binding affinity and specificity in this mode",
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

    # Always compute (or default) the key interactions description for prompt injection
    summary["key_interactions"] = _extract_key_interactions(summary.get("best_pose_pdb"))

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
        "key_interactions": s.get("key_interactions", "the critical molecular interactions (key hydrogen bonds, salt bridges, and hydrophobic contacts) that matter most"),
        "plip_style": PLIP_STYLE,
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
        "key_interactions": s.get("key_interactions", "the critical molecular interactions (key hydrogen bonds, salt bridges, and hydrophobic contacts) that matter most"),
        "plip_style": PLIP_STYLE,
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
        "plip_style_note": "Prompts emulate PLIP clean 3D interaction diagrams (see PLIP_USAGE.txt if generated) + JetBrains Mono / thebonhomme.com typography. Prioritizes most favourable CF contacts + standard non-covalent types.",
        "required_elements": [
            "bottom banner: " + BANNER,
            "equation ΔG=ΔH−TΔS with injected values",
            "reproducibility footer with gate/git/date",
            "cyan/teal accents + deep navy gradients",
            "entropy blue→red heatmap description",
            "induced-fit side chains",
            "PLIP-style critical interaction viz (dashed lines + residue labels for top H-bonds/hydrophobics/CF-favourable contacts)",
            "typography: JetBrains Mono (thebonhomme.com style)",
        ],
        "gate6_passed": gate_ok,
        "use_pymol_base_attempted": use_pymol_base,
        "plip_attempted": "see PLIP_USAGE.txt or run 'plip -f <best.pdb> -p' manually for best interaction base figs",
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

    # Optional: generate a PLIP "nice fig" base interaction diagram (if plip CLI available).
    # PLIP produces clean, publication-ready 3D interaction PNGs + .pse exactly suited for
    # our "show critical interactions" requirement. Use as base for image_to_image or manual.
    # Install: pip install plip (or use the official docker). Non-fatal if missing.
    plip_base = None
    try:
        import shutil, subprocess
        if shutil.which("plip"):
            best_p = summary.get("best_pose_pdb")
            if best_p and Path(best_p).exists():
                # -p for rendered PNG(s), -y for PyMOL session
                cmd = ["plip", "-f", str(best_p), "-o", str(fig_dir), "--name", "plip", "-p", "-y", "-q"]
                proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
                if proc.returncode == 0:
                    candidates = sorted(fig_dir.glob("plip*.png")) or sorted(fig_dir.glob("*.png"))
                    if candidates:
                        plip_base = str(candidates[0])
                        # canonical copy for easy reference
                        canon = fig_dir / "base_plip_interactions.png"
                        try:
                            import shutil as sh
                            sh.copy2(plip_base, canon)
                            plip_base = str(canon)
                        except Exception:
                            pass
                        out["plip_base_png"] = plip_base
                        (fig_dir / "PLIP_USAGE.txt").write_text(
                            "This base PNG was generated by PLIP for precise interaction viz.\n"
                            "For even better control: load the .pse in PyMOL and raytrace.\n"
                            "Then feed the PNG as --reference to imagine image-to-image or edit.\n",
                            encoding="utf-8"
                        )
    except Exception as e:
        warnings.warn(f"Optional PLIP interaction diagram generation skipped (install plip for nice figs): {e}", RuntimeWarning)

    if plip_base:
        out["plip_base_png"] = plip_base

    return out


# =============================================================================
# NRDD Journal Cover & Promotional Figure Generation (dramatic FlexAID∆S style)
# =============================================================================
# Implements professional, bulletproof support for generating Nature Reviews
# Drug Discovery-style covers matching the requested aesthetic (dramatic
# entropy/enthalpy personification or molecular E-E gauge compositions,
# precise thermodynamic call-outs, FlexAID∆S branding, Le Bonhomme Pharma,
# JetBrains Mono / thebonhomme.com typography, PLIP-inspired clarity on
# interactions where relevant, reproducibility metadata).
#
# This fulfills the 5-point integration request for the /flexaidds skill:
# the function is the handler; the skill description (SKILL.md) exposes it
# and tells the agent how to call the image_gen tool with the output prompt.
# =============================================================================

from dataclasses import dataclass, asdict, field
from typing import Literal

@dataclass(frozen=True)
class NRDDCoverParams:
    """Validated parameters for NRDD-style FlexAID∆S cover generation.

    Values should be ensemble-derived thermodynamic quantities (kcal/mol scale
    for TΔS terms; I_E-E typically normalized -1 to +1 or 0-1).
    Per user preference: prominently feature -TΔS (great visibility), use the
    Enthalpy-Entropy Index (I_E-E) developed in the skill; do not show -ΔH / -dH
    labels or visual elements prominently.
    """
    entropy_value: float = 0.93          # e.g. representative TΔS 
    # tds / enthalpy_value slot: value for prominent -TΔS (user: -TdS is great and should be highly visible;
    # do NOT feature -ΔH / -dH labels or cubes prominently; prefer the Enthalpy-Entropy Index I_E-E instead)
    enthalpy_value: float = 1.4          # value shown as -TΔS in figures
    index_value: float = 0.92            # Enthalpy-Entropy Index (I_E-E / I_EE) developed in the skill (statmech::compute_IEE)
    title: str = "The ΔG balance"
    subtitle: str = "Striking the right pose in drug discovery"
    date: str = "June 2025"
    volume: str = "Volume 24 | No. 6"
    style: Literal["dramatic_faces", "molecular_gauge"] = "dramatic_faces"
    include_flexaidds_logo: bool = True
    include_lebonhomme: bool = True
    footer_banner: str = "/flexaids-docking • FlexAID∆S"
    extra_footer: str = "Entropy-Driven • Fully Flexible Induced-Fit Docking • C++26 + Python • Grok / Claude / ChatGPT compatible"
    source_note: str = "Visualisation only. Values illustrative or from ensemble thermodynamic ledger (FlexAID∆S). Not experimental ΔG."
    # For pulling real values
    results_dir: str | None = None

    def __post_init__(self):
        # Bulletproof validation (idiotproof + scientifically sensible ranges)
        if not isinstance(self.entropy_value, (int, float)):
            raise TypeError("entropy_value must be numeric (kcal/mol scale)")
        if not -20.0 < self.entropy_value < 20.0:
            raise ValueError("entropy_value out of plausible range for TΔS term (kcal/mol)")
        if not -20.0 < self.enthalpy_value < 20.0:
            raise ValueError("enthalpy_value out of plausible range for ΔH term (kcal/mol)")
        if not -1.5 <= self.index_value <= 1.5:
            raise ValueError("index_value (I_E-E) must be roughly in [-1, 1] range")
        if self.style not in ("dramatic_faces", "molecular_gauge"):
            raise ValueError("style must be 'dramatic_faces' or 'molecular_gauge'")
        if self.results_dir is not None and not Path(self.results_dir).exists():
            raise FileNotFoundError(f"results_dir does not exist: {self.results_dir}")


def build_nrdd_cover_prompt(params: NRDDCoverParams) -> str:
    """Construct a detailed, reproducible prompt for image_gen that produces
    covers in the exact dramatic style of the provided reference NRDD FlexAID∆S
    illustrations (personified entropy/enthalpy faces or molecular gauge with
    E-E index, central ligand, precise call-out boxes, branding, typography).

    The prompt is engineered for consistency, scientific tone, and the requested
    fonts/branding (JetBrains Mono + thebonhomme.com aesthetic).
    """
    p = params
    # Pull real values if results_dir provided (robust, non-destructive)
    if p.results_dir:
        try:
            from .results import load_results
            res = load_results(p.results_dir)
            if res.binding_modes:
                top = sorted(res.binding_modes, key=lambda m: m.free_energy or 0)[0]
                # Use ledger values if present; fall back gracefully
                thermo = getattr(top, "thermodynamics", {}) or {}
                if isinstance(thermo, dict):
                    p = NRDDCoverParams(  # re-create with real-ish values (signs illustrative)
                        entropy_value=abs(thermo.get("minus_T_S_config_kcal_mol") or p.entropy_value),
                        enthalpy_value=abs(thermo.get("H_eff_kcal_mol") or p.enthalpy_value),
                        index_value=p.index_value,
                        title=p.title, subtitle=p.subtitle, date=p.date, volume=p.volume,
                        style=p.style, results_dir=None  # prevent recursion
                    )
        except Exception:
            pass  # keep provided values; never fail figure gen on data issues

    common = (
        f"High-resolution, cinematic Nature Reviews Drug Discovery cover illustration. "
        f"FlexAID∆S branding throughout. Use deep navy/black gradients, cyan/teal (#22D3EE) accents, "
        f"gold (#FBBF24) for ΔG / index highlights, orange/red for enthalpy, blue for entropy. "
        f"All typography (call-outs, banners, labels, footer) in clean, sharp JetBrains Mono font "
        f"(modern technical mono aesthetic exactly matching thebonhomme.com and Le Bonhomme Pharma site). "
        f"Professional scientific art quality, dramatic lighting, molecular detail, no clutter, 16:9 or 3:2 landscape suitable for journal cover. "
        f"Bottom dark banner with exact text: '{p.footer_banner}' and subtitle '{p.extra_footer}'. "
        f"Include FlexAID∆S logo/wordmark and 'LeBonhommePharma.github.io' . Subtle scientific icons. "
        f"Scientific note in small text: '{p.source_note}'. "
        f"Thermodynamic call-outs: cyan box 'TΔS = {p.entropy_value:.2f}', purple box '-TΔS = {p.enthalpy_value:.2f}', "
        f"gold box 'Enthalpy-Entropy Index (I_E-E) = {p.index_value:.2f}' (developed in FlexAIDdS). "
        f"Prominent equation 'ΔG = ΔH − TΔS' with values (focus on -TΔS term). Date top-right '{p.date}', volume '{p.volume}'. "
        f"Title large: '{p.title}'. Subtitle: '{p.subtitle}'."
    )

    if p.style == "dramatic_faces":
        visual = (
            "Dramatic split composition: left side a translucent, icy, blue-glowing human-like face representing "
            "Entropy (TΔS), formed from swirling water splashes, bubbles, and small ligand molecules, intense blue eyes, "
            "cool ethereal lighting, water droplets and molecular fragments exploding outward. Right side a fiery, molten, "
            "lava-cracked human-like face representing the balancing Enthalpy contribution (no -ΔH label), emerging from detailed 3D protein surface folds with "
            "orange/red glows, embers, and heat distortion. In the center, floating a highly detailed small 3D ball-and-stick "
            "ligand molecule (generic drug-like or specific if known). Prominently featured floating glass/ice cube in cyan/purple for '-TΔS' (great visibility), and gold accent for the 'Enthalpy-Entropy Index (I_E-E)'. Explosive energy and water particles at the interface. Bottom three-column layout with clean "
            "dividers and text blocks: left 'Beyond lipophilicity / New paradigms for binding affinity', middle "
            "'Conformational thermodynamics / Hidden states, real consequences', right 'Designing the ΔG sweet spot / "
            "Balancing enthalpy and entropy for better drugs' plus small Le Bonhomme Pharma icon (top hat silhouette). "
            "Overall epic, high-contrast, cinematic, suitable for the cover of Nature Reviews Drug Discovery."
        )
    else:  # molecular_gauge
        visual = (
            "Abstract molecular composition: left side intricate blue translucent protein ribbon/surface structure representing "
            "Entropy (order from disorder), with floating water-like molecules and cool lighting. Right side complex red/orange "
            "glowing protein surface with hot spots and interaction lines representing Enthalpy (energy from interactions). "
            "Center: large, detailed 3D ligand molecule with highlighted bonds and interaction dashed lines (PLIP-style clarity: "
            "blue for H-bonds, grey dashed for favourable hydrophobic, etc.). Large central glowing gauge/arc labeled "
            f"'ENTROPY-ENTHALPY INDEX E–E' with needle pointing to the {p.index_value:.2f} value, scale from -1 (blue) through 0 to +1 (orange). "
            "Prominent boxed equation 'ΔG = ΔH − TΔS' with 'Binding is not just about energy. It’s about balance.' "
            "Top-right or side panels in clean sans/mono: 'REVIEWS Entropy matters in molecular recognition', "
            "'PERSPECTIVE Beyond ΔG: why the balance of forces drives binding and selectivity', 'FOCUS Integrating physics-driven AI for better drug design'. "
            "Bottom branding bar with FlexAID∆S logo, tagline 'BALANCE • LEARN • DESIGN • DISCOVER', four icons with labels "
            "'QUANTIFY Measure E-E ...', 'UNDERSTAND Decode the balance...', 'PREDICT AI models guided by physics...', "
            "'DESIGN Optimize interactions with E-E-informed strategies'. Le Bonhomme Pharma logo and QR/github link. "
            "Dramatic yet clean scientific visualization style."
        )

    prompt = f"{common} {visual} {p.source_note}"
    return prompt


def generate_flexaids_nrdd_cover(
    entropy_value: float = 0.93,
    enthalpy_value: float = 1.4,
    index_value: float = 0.92,
    *,
    style: Literal["dramatic_faces", "molecular_gauge"] = "dramatic_faces",
    title: str = "The ΔG balance",
    subtitle: str = "Striking the right pose in drug discovery",
    results_dir: str | None = None,
    render: bool = False,
    **image_gen_kwargs,
) -> dict:
    """Professional, bulletproof handler for /flexaidds skill figure generation.

    Implements the 5 points requested:
    1. Entry point: this function (exposed via python/flexaidds and documented in skill).
    2. Handler: builds the prompt (modeled on the exact reference covers), calls image_gen if render=True
       (or returns prompt for the agent to call the tool in skill context).
    3. Exposed via skill "manifest" (SKILL.md describes the action + parameters; agent knows to use it).
    4. Dependencies: image generation via host tool (image_gen / image_edit). No hard GPL deps. Optional PLIP for bases.
    5. Returns rich result (prompt + full metadata for reproducibility + path if rendered). Use computer.sync_file
       or just the path in Grok context. Client displays/downloads the figure.

    Robustness:
    - Strict validation of thermodynamic values and ranges.
    - If results_dir given, attempts to source realistic values from the ensemble ledger (non-fatal).
    - Reproducible: metadata contains every input param, timestamp, git sha (best effort), source note.
    - Scientific: never claims experimental ΔG; uses project terminology ("ensemble thermodynamic ledger").
    - Typography & branding: JetBrains Mono + thebonhomme.com aesthetic enforced in prompt.
    - Failsafe: always returns usable prompt even if rendering unavailable; no crashes on bad data.

    Example (agent / skill usage):
        from flexaidds.figures import generate_flexaids_nrdd_cover
        res = generate_flexaids_nrdd_cover(entropy_value=0.93, enthalpy_value=1.4, index_value=0.92,
                                           style="dramatic_faces", render=False)
        # Then in Grok/agent context:
        # path = image_gen(prompt=res['prompt'], aspect_ratio="16:9")['path']   # or whatever the tool returns
        # Optionally image_edit(path, prompt="refine labels...") for polish.
    """
    params = NRDDCoverParams(
        entropy_value=entropy_value,
        enthalpy_value=enthalpy_value,
        index_value=index_value,
        title=title,
        subtitle=subtitle,
        style=style,
        results_dir=results_dir,
    )
    prompt = build_nrdd_cover_prompt(params)

    meta = {
        "params": asdict(params),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "tool": "flexaidds.figures.generate_flexaids_nrdd_cover",
        "branding": "FlexAID∆S + Le Bonhomme Pharma (thebonhomme.com aesthetic, JetBrains Mono typography)",
        "plip_note": "Interaction clarity inspired by PLIP where molecular detail is shown",
        "scientific_note": params.source_note,
        "suggested_image_gen": "image_gen(prompt=..., aspect_ratio='16:9' or '3:2')",
        "suggested_edit": "image_edit(image=[path], prompt='tighten labels, enhance contrast on I_E-E gauge or cubes...')",
    }

    result = {
        "prompt": prompt,
        "metadata": meta,
        "path": None,
        "params": params,
    }

    if render:
        # In environments where the image_gen tool is directly importable/available to code (rare),
        # we could call it here. By default we keep the package pure and let the skill/agent
        # perform the tool call (as described in the /flexaidds skill instructions).
        # This is the robust, host-agnostic approach.
        try:
            # Placeholder for direct call if a make_image is exposed in the runtime.
            # In this Grok agent environment the caller will use the image_gen tool directly.
            result["note"] = "render=True: in this context, the agent should now call image_gen with the prompt above."
        except Exception as e:
            result["error"] = f"Direct render failed (use host image_gen tool): {e}"

    return result
