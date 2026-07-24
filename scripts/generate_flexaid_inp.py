#!/usr/bin/env python3
"""Generate classic FlexAID CONFIG.inp + ga.inp work trees for arms A / B0 / B.

Primary three-engine claim path (TIER-1 cognate, no native seed):
  A  = FlexAID 2015-era (JCIM 2015) binary  (CF election; TEMPER 0)
  B  = FlexAID master, entropy ranking ON   (TEMPER 21 — operator-optimized)
  C0 = FlexAIDdS (separate runner; not this generator)

Optional / deferred:
  B0 = FlexAID master, entropy OFF          (TEMPER 0 → CF clustering forced)
  B@298 = prior protocol default TEMPER 298 (use --temper 298 to restore)

Inputs per target (from queue):
  <PDB>_apo.pdb, <PDB>_ligand.sdf, cleft spheres (GetCleft or queue site)

ProcessLigand (processligand-py) produces ligand .inp/.ic and typed target PDB.
Sphere files with HETATM records are rewritten as ATOM for classic FlexAID
read_spheres (A/B binaries only accept ATOM in some builds).

P0 canary prep gates (pilot8 failure recovery):
  * Clean apo TARGET: strip HOH/WAT always; metals unless FLEXAIDDS_KEEP_METALS=1
    (see scripts/clean_target_apo.py). Orphan CONECT cleaned.
  * Ligand integrity after ProcessLigand on LIG_ref (heavy count + CONECT bonds).
    INI is only available after FlexAID starts — post-INI preflight:
      python3 scripts/validate_ligand_integrity.py --work <work> --require-ini
  * Native CF oracle (after results exist) — ranking experiments forbidden until pass:
      python3 scripts/native_cf_oracle_gate.py --work <work> --results <out>/<pdb>

Seed note: staged FlexAID A/B binaries seed with time(0); STRTSEED is written
into ga.inp for provenance / future binaries, but may be ignored by current
Mach-O pins. Restarts are separate process invocations.

Usage:
  python3 scripts/generate_flexaid_inp.py \\
    --queue-root "$FLEXAIDDS_QUEUE_ROOT" \\
    --pdb 1GPK --arm B0 \\
    --work-dir /path/to/work/B0/1GPK

  python3 scripts/generate_flexaid_inp.py --queue-root "$Q" --pilot8 --arms A,B0,B
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

# Local script imports (same directory) for P0 canary prep gates.
_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

try:
    from clean_target_apo import clean_apo_file, resolve_keep_flags
except ImportError:  # pragma: no cover
    clean_apo_file = None  # type: ignore[assignment]
    resolve_keep_flags = None  # type: ignore[assignment]

try:
    from ligand_integrity import validate_ligand_integrity, format_result as _fmt_lig
except ImportError:  # pragma: no cover
    validate_ligand_integrity = None  # type: ignore[assignment]
    _fmt_lig = None  # type: ignore[assignment]

# Protocol defaults (three_engine_entropy_comparison.md)
SEED_BASE = 20260714
DEFAULT_POP = 1000
DEFAULT_GEN = 2000  # claim freeze 2026-07-15 (was 6000)
DEFAULT_RESTARTS = 5
DEFAULT_MAXRES = 50  # match FlexAIDdS cluster emit ceiling for fair S3/BCR
MATRIX_NAME = "MC_st0r5.2_6.dat"
# Lab production pin (v119-era claim / three-engine). Re-pinned after df2f36c58
# matrix-value regression (9dc9) inverted CF.com ranking; see ROOTCAUSE 2026-07-23.
MATRIX_MD5_PIN = "72d7c7396702331d96ff12d18f831796"
LIGAND_RESNUM = 9999
ATOM_INDEX = 90000

# PSHARE niching — match LIB/config_defaults.h (sharing_scale=10, sharing_peaks=5)
# and JCIM-era fair logs (scale 10). Pilot8 incorrectly used SHARESCL 0.20 / SHAREPEK 3
# (0.20 is ADAPTKCO k4, not niche scale) → ~50× larger niche radius vs production.
# Override: FLEXAIDDS_GA_SHARESCL / FLEXAIDDS_GA_SHAREPEK / FLEXAIDDS_GA_SHAREALF.
DEFAULT_SHAREALF = 4.0
DEFAULT_SHAREPEK = 5.0
DEFAULT_SHARESCL = 10.0

PILOT8 = ["1G9V", "1GPK", "1MEH", "1P62", "1Q4G", "1R9O", "1T40", "2BYS"]

# P0 canary prep: strip waters/metals for redock TARGET (see clean_target_apo.py).
# Override with FLEXAIDDS_KEEP_HOH=1 / FLEXAIDDS_KEEP_METALS=1 or --keep-hoh/--keep-metals.
# Ligand integrity after ProcessLigand; INI only exists after FlexAID start — see
# scripts/validate_ligand_integrity.py --require-ini for post-INI preflight.

# TEMPER for arm B default = 21 (LP-optimized engine soft-T for FO/ACF emission).
# Soft-β here is FlexAID β=1/T on CF a.u. — not k_B·T kcal; not DatasetRunner Softβ S1.
# Arm B FO + TEMPER ≠ Softβ rescoring of frozen CF ensembles (see softbeta_election_policy.md).
# Override per run: --temper N  (applies to all listed arms that receive entropy).
# Arm B: CLUSTA FO with exactly ONE literature MinPts (Ankerst[10-20]+Sander 2*dim+Ester floor).
# Engine chooses MinPts once in FastOPTICS_cluster.cpp — no multi-scale ladder in CONFIG.
# See docs/implementation/fo_minpts_literature.md and 3dsig_red_pair_protocol.md §2.1.
ARM_SPEC = {
    "A": {
        "bin_arm": "A",
        "temper": 0,
        "clusta": "CF",
        "label": "FlexAID historical pin TEMPER0 CLUSTA CF (must differ SHA from master B)",
    },
    "B0": {
        "bin_arm": "B",
        "temper": 0,
        "clusta": "CF",
        # Twin of A when bin A==B SHA — not an independent control. Prefer arm B for master.
        "label": "master TEMPER0 CLUSTA CF (twin of A if same binary; not claim control)",
        "science_control": False,
    },
    "B": {
        "bin_arm": "B",
        "temper": 21,
        "clusta": "FO",
        "fo_minpts_policy": "single_literature",
        "label": "master TEMPER21 CLUSTA FO single-literature-MinPts",
    },
    # Arm C: FO @ 298K — only after panel native CF oracle PASS (claim gate).
    "C": {
        "bin_arm": "B",
        "temper": 298,
        "clusta": "FO",
        "fo_minpts_policy": "single_literature",
        "label": "master TEMPER298 CLUSTA FO (requires native CF oracle PASS)",
        "requires_oracle_pass": True,
    },
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_pshare_knobs() -> tuple[float, float, float]:
    """Resolve SHAREALF / SHAREPEK / SHARESCL from env or production defaults."""
    def _f(name: str, default: float) -> float:
        raw = os.environ.get(name, "").strip()
        if not raw:
            return default
        try:
            return float(raw)
        except ValueError as e:
            raise ValueError(f"Invalid {name}={raw!r}") from e

    return (
        _f("FLEXAIDDS_GA_SHAREALF", DEFAULT_SHAREALF),
        _f("FLEXAIDDS_GA_SHAREPEK", DEFAULT_SHAREPEK),
        _f("FLEXAIDDS_GA_SHARESCL", DEFAULT_SHARESCL),
    )


def stable_seed(pdb_id: str, restart_i: int, seed_base: int = SEED_BASE) -> int:
    """Protocol: SEED_BASE + stable_hash(pdb_id, restart_i) → positive 31-bit int."""
    h = hashlib.sha256(f"{pdb_id.upper()}:{restart_i}".encode()).hexdigest()
    return int(seed_base + (int(h[:8], 16) % 1_000_000_000))


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve_processligand() -> Path:
    """Locate ProcessLigand binary without baking in machine-specific defaults."""
    env = os.environ.get("FLEXAIDDS_PROCESSLIGAND")
    if env and Path(env).is_file() and os.access(env, os.X_OK):
        return Path(env).resolve()

    try:
        import processligandpy  # type: ignore

        cand = Path(processligandpy.__file__).resolve().parent / "bin" / "ProcessLigand"
        if cand.is_file() and os.access(cand, os.X_OK):
            return cand
    except ImportError:
        pass

    for base in (repo_root() / ".venv-processligand", Path.home() / ".venv-processligand"):
        if base.is_dir():
            matches = list(
                base.glob("lib/python*/site-packages/processligandpy/bin/ProcessLigand")
            )
            if matches and os.access(matches[0], os.X_OK):
                return matches[0].resolve()

    which = shutil.which("ProcessLigand")
    if which:
        return Path(which).resolve()

    raise FileNotFoundError(
        "ProcessLigand not found. Install with: "
        "python3 -m venv .venv-processligand && "
        ".venv-processligand/bin/pip install processligand-py\n"
        "Or set FLEXAIDDS_PROCESSLIGAND to the binary path."
    )


def ensure_babel_env(processligand: Path) -> None:
    if os.environ.get("BABEL_DATADIR"):
        return
    root = processligand.parent.parent
    for cand in (root / "share" / "openbabel" / "3.1.1", root / "bin" / "data"):
        if cand.is_dir():
            os.environ["BABEL_DATADIR"] = str(cand)
            return


def run_processligand(
    processligand: Path,
    input_file: Path,
    out_prefix: Path,
    *,
    target: bool = False,
    atom_index: int = ATOM_INDEX,
    ref: bool = True,
) -> None:
    ensure_babel_env(processligand)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(processligand),
        "-f",
        str(input_file),
        "-o",
        str(out_prefix),
        "-d",
        "-pf",
    ]
    if target:
        cmd.append("-target")
    else:
        cmd.extend(
            ["--atom_index", str(atom_index), "--res_number", str(LIGAND_RESNUM)]
        )
        if ref:
            cmd.append("-ref")

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"ProcessLigand failed (rc={proc.returncode}) for {input_file}\n"
            f"stdout:\n{proc.stdout[-2000:]}\nstderr:\n{proc.stderr[-2000:]}"
        )


def convert_spheres_to_atom(src: Path, dst: Path) -> int:
    """Rewrite HETATM→ATOM so classic FlexAID read_spheres accepts the file."""
    n = 0
    lines_out: List[str] = []
    for line in src.read_text(errors="replace").splitlines():
        if line.startswith("HETATM"):
            lines_out.append("ATOM  " + line[6:])
            n += 1
        elif line.startswith("ATOM  "):
            lines_out.append(line)
            n += 1
        else:
            lines_out.append(line)
    if n == 0:
        raise ValueError(f"no sphere ATOM/HETATM records in {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("\n".join(lines_out) + "\n")
    return n


def _pdb_xyz(line: str) -> Optional[Tuple[float, float, float]]:
    if len(line) < 54:
        return None
    try:
        return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except ValueError:
        return None


def prune_spheres_near_ligand(
    spheres_pdb: Path,
    ligand_ref: Path,
    *,
    max_dist: float = 8.0,
    max_spheres: int = 2500,
) -> Dict[str, Any]:
    """Keep LOCCLF spheres near LIG_ref; cap count for GA budget.

    Oversized native-occupied GetCleft clouds (10k–25k spheres) dilute search
    under R=1 budgets. Prune by distance to any LIG_ref heavy atom, then
    nearest-first cap to ``max_spheres``.
    """
    import math

    lig_xyz: List[Tuple[float, float, float]] = []
    if ligand_ref.is_file():
        for line in ligand_ref.read_text(errors="replace").splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                xyz = _pdb_xyz(line)
                if xyz:
                    lig_xyz.append(xyz)
    if not lig_xyz:
        return {"pruned": False, "reason": "no_ligand_coords", "n_in": 0, "n_out": 0}

    header: List[str] = []
    spheres: List[Tuple[float, str]] = []  # (min_dist, line)
    n_in = 0
    for line in spheres_pdb.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            n_in += 1
            xyz = _pdb_xyz(line)
            if xyz is None:
                continue
            dmin = min(
                math.sqrt(
                    (xyz[0] - lx) ** 2 + (xyz[1] - ly) ** 2 + (xyz[2] - lz) ** 2
                )
                for lx, ly, lz in lig_xyz
            )
            if dmin <= max_dist:
                # normalize to ATOM record for classic FlexAID
                atom_line = ("ATOM  " + line[6:]) if line.startswith("HETATM") else line
                spheres.append((dmin, atom_line))
        else:
            header.append(line)

    spheres.sort(key=lambda t: t[0])
    kept = spheres[: max(1, int(max_spheres))]
    lines_out = header + [ln for _, ln in kept]
    if not kept:
        return {"pruned": False, "reason": "empty_after_prune", "n_in": n_in, "n_out": 0}
    spheres_pdb.write_text("\n".join(lines_out) + "\n")
    return {
        "pruned": True,
        "n_in": n_in,
        "n_out": len(kept),
        "max_dist": max_dist,
        "max_spheres": max_spheres,
        "max_kept_dist": kept[-1][0] if kept else None,
    }


def parse_fledih_count(ligand_inp: Path) -> int:
    text = ligand_inp.read_text(errors="replace")
    return sum(1 for line in text.splitlines() if line.startswith("FLEDIH"))


def resolve_soft_wall_cutoff(explicit: Optional[float] = None) -> float:
    """Soft-core wall cutoff (Å) for classic CONFIG SOFTWA.

    Default 0.40 matches DatasetRunner / top.cpp. Set FLEXAIDDS_SOFT_WALL=0 for
    legacy hard r^-12 (not recommended — CF.com overpacking pathology).
    """
    if explicit is not None:
        return float(explicit)
    env = os.environ.get("FLEXAIDDS_SOFT_WALL", "").strip()
    if env:
        try:
            return float(env)
        except ValueError:
            pass
    return 0.40


def write_config(
    path: Path,
    *,
    target_pdb: Path,
    ligand_inp: Path,
    spheres_pdb: Path,
    matrix_path: Path,
    depspa: Path,
    statep: Path,
    tempop: Path,
    rmsd_ref: Optional[Path],
    temper: int,
    n_flex: int,
    maxres: int = DEFAULT_MAXRES,
    soft_wall_cutoff: Optional[float] = None,
) -> None:
    soft_wa = resolve_soft_wall_cutoff(soft_wall_cutoff)
    lines: List[str] = [
        "# generated by scripts/generate_flexaid_inp.py",
        f"PDBNAM {target_pdb}",
        f"INPLIG {ligand_inp}",
        "COMPLF VCT",
        f"RNGOPT LOCCLF {spheres_pdb}",
        f"OPTIMZ {LIGAND_RESNUM} - -1",
        f"OPTIMZ {LIGAND_RESNUM} - 0",
    ]
    for i in range(1, n_flex + 1):
        lines.append(f"OPTIMZ {LIGAND_RESNUM} - {i}")

    if rmsd_ref is not None:
        lines.append(f"RMSDST {rmsd_ref}")

    lines.extend(
        [
            f"IMATRX {matrix_path}",
            "PERMEA 0.9",
            # SOFTWA: 6-char classic keyword → FA->soft_wall_cutoff (Å)
            f"SOFTWA {soft_wa:.2f}",
            "VARANG 5.0",
            "VARDIH 5.0",
            "VARFLX 10.0",
            "SPACER 0.375",
            "SLVTYP 40",
            "EXCHET",
            "METOPT GA",
            "VCTPLA R",
            "NORMAR",
            "VINDEX",
            f"STATEP {statep}",
            f"TEMPOP {tempop}",
            f"DEPSPA {depspa}",
            f"MAXRES {maxres}",
            f"TEMPER {temper}",
            # FO path: engine runs ONE FastOPTICS pass at literature MinPts (no ladder).
            "CLUSTA CF" if temper <= 0 else "CLUSTA FO",
            "CLRMSD 2.0",
            "ENDINP",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def write_ga(
    path: Path,
    *,
    pop: int,
    gen: int,
    seed: int,
    sharealf: float = DEFAULT_SHAREALF,
    sharepek: float = DEFAULT_SHAREPEK,
    sharescl: float = DEFAULT_SHARESCL,
) -> None:
    """Write ga.inp with production PSHARE niching (JCIM / config_defaults).

    SHARESCL defaults to 10 (engine GB->scale). Pilot8 wrongly used 0.20
    (copy-paste from ADAPTKCO k4) which expands niche radius ~50× vs historical.
    Override via FLEXAIDDS_GA_SHARESCL / FLEXAIDDS_GA_SHAREPEK / FLEXAIDDS_GA_SHAREALF.
    """
    lines = [
        "# generated by scripts/generate_flexaid_inp.py",
        "# STRTSEED may be ignored by staged A/B binaries (time-based srand); kept for provenance",
        f"NUMCHROM {pop}",
        f"NUMGENER {gen}",
        "ADAPTVGA 1",
        "ADAPTKCO 0.95 0.50 0.95 0.20",
        "CROSRATE 0.900",
        "MUTARATE 0.025",
        "POPINIMT RANDOM",
        "FITMODEL PSHARE",
        f"SHAREALF {sharealf:.2f}",
        f"SHAREPEK {sharepek:.2f}",
        f"SHARESCL {sharescl:.2f}",
        "REPMODEL BOOM",
        "BOOMFRAC 1.00",
        # Emit up to DEFAULT_MAXRES cluster heads (parity with FlexAIDdS max_results=50)
        f"PRINTCHR {min(DEFAULT_MAXRES, pop)}",
        "OUTGENER 1",
        f"STRTSEED {seed}",
    ]
    path.write_text("\n".join(lines) + "\n")


def find_sphere_file(queue_root: Path, pdb_id: str) -> Path:
    pdb = pdb_id.upper()
    env_root = os.environ.get("FLEXAIDDS_SPHERE_ROOT")
    candidates: List[Path] = []
    if env_root:
        candidates.extend(Path(env_root).glob(f"{pdb}/**/*_occupied_sph_1.pdb"))
        candidates.extend(Path(env_root).glob(f"{pdb}/*sph*.pdb"))

    q_in = queue_root / "inputs" / "astex_diverse" / pdb
    candidates.extend(q_in.glob("*_sph*.pdb"))
    candidates.extend(q_in.glob("*_cleft*.pdb"))

    icloud = os.environ.get("FLEXAIDDS_ICLOUD")
    if not icloud:
        icloud = str(
            Path.home()
            / "Library"
            / "Mobile Documents"
            / "com~apple~CloudDocs"
            / "FlexAIDdS_benchmarks"
        )
    cav = (
        Path(icloud)
        / "astex_entropy"
        / "native_headtohead_20260707"
        / "cavities"
        / "native"
        / pdb
        / "native_occupied"
        / f"{pdb}_occupied_sph_1.pdb"
    )
    candidates.append(cav)

    for c in candidates:
        if c.is_file():
            return c.resolve()

    site = q_in / f"{pdb}_site.pdb"
    if site.is_file():
        return site.resolve()

    raise FileNotFoundError(
        f"No sphere/site file for {pdb}. Set FLEXAIDDS_SPHERE_ROOT or stage GetCleft spheres."
    )


def prepare_target(
    queue_root: Path,
    pdb_id: str,
    arm: str,
    work_root: Path,
    *,
    pop: int = DEFAULT_POP,
    gen: int = DEFAULT_GEN,
    restarts: int = DEFAULT_RESTARTS,
    seed_base: int = SEED_BASE,
    processligand: Optional[Path] = None,
    force: bool = False,
    temper_override: Optional[int] = None,
    keep_hoh: bool = False,
    keep_metals: bool = False,
    skip_clean_apo: bool = False,
    skip_ligand_gate: bool = False,
    max_bond: float = 3.0,
    sphere_max: int = 0,
    sphere_max_dist: float = 8.0,
) -> Path:
    if arm not in ARM_SPEC:
        raise ValueError(f"unknown arm {arm}; expected one of {list(ARM_SPEC)}")

    spec = ARM_SPEC[arm]
    pdb = pdb_id.upper()
    inp = queue_root / "inputs" / "astex_diverse" / pdb
    apo = inp / f"{pdb}_apo.pdb"
    lig = inp / f"{pdb}_ligand.sdf"
    if not apo.is_file() or not lig.is_file():
        raise FileNotFoundError(f"missing apo/ligand under {inp}")

    matrix = queue_root / "data" / MATRIX_NAME
    if not matrix.is_file():
        raise FileNotFoundError(f"matrix missing: {matrix}")
    got = md5_file(matrix)
    if got != MATRIX_MD5_PIN:
        raise RuntimeError(f"matrix MD5 {got} != pin {MATRIX_MD5_PIN}")

    depspa = (queue_root / "data").resolve()
    work = work_root / arm / pdb
    if work.exists() and force:
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)

    pl = processligand or resolve_processligand()

    lig_work = work / f"{pdb}_ligand.sdf"
    if not lig_work.exists() or force:
        shutil.copy2(lig, lig_work)
    lig_prefix = work / "LIG"
    if not (work / "LIG.inp").is_file() or force:
        run_processligand(pl, lig_work, lig_prefix, target=False, ref=True)

    ligand_inp = work / "LIG.inp"
    ligand_ic = work / "LIG.ic"
    ligand_ref = work / "LIG_ref.pdb"
    if not ligand_inp.is_file():
        raise RuntimeError(f"ProcessLigand did not produce {ligand_inp}")
    if not ligand_ic.is_file():
        raise RuntimeError(f"ProcessLigand did not produce {ligand_ic}")

    # ── P0 ligand integrity gate (post-ProcessLigand, pre-FlexAID) ──────────
    # Self-check LIG_ref heavy atoms + CONECT bonds. INI is only available after
    # FlexAID starts — document post-INI: validate_ligand_integrity.py --require-ini
    if not skip_ligand_gate and validate_ligand_integrity is not None:
        if not ligand_ref.is_file():
            raise RuntimeError(
                f"ProcessLigand did not produce {ligand_ref} "
                "(required for ligand integrity gate)"
            )
        lig_res = validate_ligand_integrity(
            ligand_ref,
            None,
            max_bond=max_bond,
            check_bonds=True,
        )
        gate_path = work / "ligand_integrity_prep.json"
        gate_path.write_text(
            json.dumps(lig_res.as_dict(), indent=2) + "\n", encoding="utf-8"
        )
        if not lig_res.ok:
            detail = _fmt_lig(lig_res) if _fmt_lig else lig_res.messages
            raise RuntimeError(
                f"ligand integrity gate failed for {pdb}:\n{detail}\n"
                f"report: {gate_path}"
            )
        print(
            f"  ligand_integrity OK {pdb} heavy={lig_res.ref_heavy} "
            f"(INI post-flight: scripts/validate_ligand_integrity.py "
            f"--work {work} --require-ini)"
        )
    elif not ligand_ref.is_file():
        # Soft note only when gate skipped
        print(f"  WARN: missing {ligand_ref} (RMSDST will be omitted)", file=sys.stderr)

    # ── Clean apo for redock TARGET (waters always; metals default OFF) ─────
    # Only rewrite apo_work when we will (re)run ProcessLigand -target, so
    # TARGET.inp.pdb stays consistent with the cleaned apo on disk.
    apo_work = work / f"{pdb}_apo.pdb"
    apo_raw = work / f"{pdb}_apo.raw.pdb"
    tgt_prefix = work / "TARGET"
    target_pdb = work / "TARGET.inp.pdb"
    need_target = force or not target_pdb.is_file()

    if resolve_keep_flags is not None:
        env_hoh, env_metals = resolve_keep_flags(
            keep_hoh=True if keep_hoh else None,
            keep_metals=True if keep_metals else None,
        )
        use_hoh = bool(keep_hoh) or bool(env_hoh)
        use_metals = bool(keep_metals) or bool(env_metals)
    else:
        use_hoh, use_metals = keep_hoh, keep_metals

    if need_target:
        if not apo_raw.is_file() or force:
            shutil.copy2(apo, apo_raw)
        if not skip_clean_apo and clean_apo_file is not None:
            clean_report = clean_apo_file(
                apo_raw if apo_raw.is_file() else apo,
                apo_work,
                keep_hoh=use_hoh,
                keep_metals=use_metals,
                ligand_ref=ligand_ref if ligand_ref.is_file() else None,
            )
            (work / "clean_apo_report.json").write_text(
                json.dumps(clean_report.as_dict(), indent=2) + "\n",
                encoding="utf-8",
            )
            print(
                f"  clean_apo {pdb}: waters={clean_report.waters_removed} "
                f"metals={clean_report.metals_removed} "
                f"metals_near_lig={getattr(clean_report, 'metals_kept_near_ligand', 0)} "
                f"atoms {clean_report.atoms_in}->{clean_report.atoms_out} "
                f"(keep_hoh={use_hoh} keep_metals={use_metals})"
            )
        else:
            if not apo_work.exists() or force:
                shutil.copy2(apo, apo_work)
    else:
        if not apo_work.exists():
            # Recover missing apo_work without forcing target rebuild
            if not skip_clean_apo and clean_apo_file is not None:
                if not apo_raw.is_file():
                    shutil.copy2(apo, apo_raw)
                clean_apo_file(
                    apo_raw,
                    apo_work,
                    keep_hoh=use_hoh,
                    keep_metals=use_metals,
                    ligand_ref=ligand_ref if ligand_ref.is_file() else None,
                )
            else:
                shutil.copy2(apo, apo_work)

    if need_target:
        run_processligand(pl, apo_work, tgt_prefix, target=True)
    if not target_pdb.is_file():
        alt = work / f"{pdb}_apo.inp.pdb"
        if alt.is_file():
            shutil.copy2(alt, target_pdb)
        else:
            found = list(work.glob("*.inp.pdb"))
            if found:
                shutil.copy2(found[0], target_pdb)
            else:
                raise RuntimeError(
                    f"ProcessLigand did not produce target .inp.pdb under {work}"
                )

    sph_src = find_sphere_file(queue_root, pdb)
    sph_dst = work / f"{pdb}_spheres.pdb"
    if sph_src.name.endswith("_site.pdb") and "sph" not in sph_src.name.lower():
        (work / "SPHERE_SOURCE_WARNING.txt").write_text(
            f"WARNING: using pocket-atom site file (not GetCleft spheres): {sph_src}\n"
            "Classic FlexAID LOCCLF expects sphere radii in B-factor; results may be invalid.\n"
        )
        shutil.copy2(sph_src, sph_dst)
    else:
        convert_spheres_to_atom(sph_src, sph_dst)
        (work / "sphere_source.txt").write_text(str(sph_src) + "\n")

    # Optional LOCCLF prune (arg or env): oversized native-occupied clouds dilute GA
    sph_max = int(sphere_max or 0)
    if sph_max <= 0:
        sph_max = int(os.environ.get("FLEXAIDDS_SPHERE_MAX", "0") or "0")
    sph_dist = float(sphere_max_dist)
    env_dist = os.environ.get("FLEXAIDDS_SPHERE_MAX_DIST", "").strip()
    if env_dist:
        try:
            sph_dist = float(env_dist)
        except ValueError:
            pass
    if sph_max > 0 and ligand_ref.is_file():
        prune_info = prune_spheres_near_ligand(
            sph_dst,
            ligand_ref,
            max_dist=sph_dist,
            max_spheres=sph_max,
        )
        (work / "sphere_prune.json").write_text(
            json.dumps(prune_info, indent=2) + "\n"
        )
        if prune_info.get("pruned"):
            print(
                f"  sphere prune {pdb}: {prune_info['n_in']} → {prune_info['n_out']} "
                f"(d≤{sph_dist}Å max={sph_max})"
            )

    n_flex = parse_fledih_count(ligand_inp)
    temper = int(spec["temper"]) if temper_override is None else int(temper_override)

    statep = work / "state"
    tempop = work / "tmp"
    statep.mkdir(exist_ok=True)
    tempop.mkdir(exist_ok=True)

    write_config(
        work / "CONFIG.inp",
        target_pdb=target_pdb.resolve(),
        ligand_inp=ligand_inp.resolve(),
        spheres_pdb=sph_dst.resolve(),
        matrix_path=matrix.resolve(),
        depspa=depspa,
        statep=statep.resolve(),
        tempop=tempop.resolve(),
        rmsd_ref=ligand_ref.resolve() if ligand_ref.is_file() else None,
        temper=temper,
        n_flex=n_flex,
        soft_wall_cutoff=resolve_soft_wall_cutoff(),
    )

    sharealf, sharepek, sharescl = resolve_pshare_knobs()
    for r in range(restarts):
        seed = stable_seed(pdb, r, seed_base)
        rdir = work / f"restart_{r}"
        rdir.mkdir(exist_ok=True)
        write_ga(
            rdir / "ga.inp",
            pop=pop,
            gen=gen,
            seed=seed,
            sharealf=sharealf,
            sharepek=sharepek,
            sharescl=sharescl,
        )
        r_state = rdir / "state"
        r_tmp = rdir / "tmp"
        r_state.mkdir(exist_ok=True)
        r_tmp.mkdir(exist_ok=True)
        write_config(
            rdir / "CONFIG.inp",
            target_pdb=target_pdb.resolve(),
            ligand_inp=ligand_inp.resolve(),
            spheres_pdb=sph_dst.resolve(),
            matrix_path=matrix.resolve(),
            depspa=depspa,
            statep=r_state.resolve(),
            tempop=r_tmp.resolve(),
            rmsd_ref=ligand_ref.resolve() if ligand_ref.is_file() else None,
            temper=temper,
            n_flex=n_flex,
            soft_wall_cutoff=resolve_soft_wall_cutoff(),
        )
        (rdir / "seed.txt").write_text(f"{seed}\n")

    clusta = spec.get("clusta") or ("CF" if temper <= 0 else "FO")
    meta = {
        "arm": arm,
        "bin_arm": spec["bin_arm"],
        "label": spec["label"],
        "pdb_id": pdb,
        "temper": temper,
        "clusta": clusta,
        "fo_minpts_policy": spec.get(
            "fo_minpts_policy",
            "single_literature" if clusta == "FO" else "n/a",
        ),
        "pop": pop,
        "gen": gen,
        "restarts": restarts,
        "seed_base": seed_base,
        "sharealf": sharealf,
        "sharepek": sharepek,
        "sharescl": sharescl,
        "n_flex_bonds": n_flex,
        "matrix_md5": got,
        "matrix_path": str(matrix.resolve()),
        "sphere_source": str(sph_src),
        "processligand": str(pl),
        "work_dir": str(work.resolve()),
        "seeds": [stable_seed(pdb, r, seed_base) for r in range(restarts)],
        "cli": "FlexAID CONFIG.inp ga.inp output_prefix",
        "seed_limitation": (
            "Staged A/B FlexAID binaries use time(0) srand; STRTSEED in ga.inp "
            "is provenance-only unless the binary implements it."
        ),
        "clean_apo": not skip_clean_apo,
        "ligand_integrity_gate": not skip_ligand_gate,
        "soft_wall_cutoff": resolve_soft_wall_cutoff(),
        "sphere_max": sph_max if sph_max > 0 else None,
        "sphere_max_dist": sph_dist if sph_max > 0 else None,
        "science_control": bool(spec.get("science_control", True)),
        "requires_oracle_pass": bool(spec.get("requires_oracle_pass", False)),
        "post_ini_preflight": (
            "python3 scripts/validate_ligand_integrity.py "
            f"--work {work.resolve()} --require-ini"
        ),
        "native_cf_oracle": (
            "python3 scripts/native_cf_oracle_gate.py "
            f"--work {work.resolve()} --results <OUT>/{pdb}"
        ),
        "ranking_forbidden_until": (
            "native_cf_oracle_gate.py exit 0 on canary "
            "(CF_native <= best_ga_cf + tol); ranking experiments forbidden until then"
        ),
    }
    (work / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    return work


def load_pdb_list(
    queue_root: Path, pilot8: bool, pdb: Optional[str], list_file: Optional[str]
) -> List[str]:
    if pdb:
        return [pdb.upper()]
    if pilot8:
        return list(PILOT8)
    if list_file:
        return [
            line.strip().upper()
            for line in Path(list_file).read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ]
    d = queue_root / "inputs" / "astex_diverse"
    return sorted(
        p.name
        for p in d.iterdir()
        if p.is_dir() and (p / f"{p.name}_apo.pdb").is_file()
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--queue-root",
        default=os.environ.get("FLEXAIDDS_QUEUE_ROOT", ""),
        help="Three-engine queue root",
    )
    ap.add_argument(
        "--work-root", default="", help="Work tree root (default: $QUEUE/work)"
    )
    ap.add_argument("--pdb", default="", help="Single PDB id")
    ap.add_argument("--pilot8", action="store_true", help="Generate pilot8 panel")
    ap.add_argument("--list-file", default="", help="File with one PDB id per line")
    ap.add_argument("--arms", default="A,B0,B", help="Comma-separated arms")
    ap.add_argument("--pop", type=int, default=DEFAULT_POP)
    ap.add_argument("--gen", type=int, default=DEFAULT_GEN)
    ap.add_argument("--restarts", type=int, default=DEFAULT_RESTARTS)
    ap.add_argument("--seed-base", type=int, default=SEED_BASE)
    ap.add_argument(
        "--temper",
        type=int,
        default=None,
        help="Override TEMPER for all prepared arms (e.g. 21 for optimized entropy B)",
    )
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="Print plan only")
    ap.add_argument(
        "--keep-hoh",
        action="store_true",
        help="Keep waters in TARGET apo (default strip; or FLEXAIDDS_KEEP_HOH=1)",
    )
    ap.add_argument(
        "--keep-metals",
        action="store_true",
        help="Keep metals in TARGET apo (default strip; or FLEXAIDDS_KEEP_METALS=1)",
    )
    ap.add_argument(
        "--skip-clean-apo",
        action="store_true",
        help="Do not strip waters/metals from apo before ProcessLigand -target",
    )
    ap.add_argument(
        "--skip-ligand-gate",
        action="store_true",
        help="Skip LIG_ref integrity gate after ProcessLigand",
    )
    ap.add_argument(
        "--max-bond",
        type=float,
        default=3.0,
        help="Max CONECT bond (Å) for ligand integrity gate (default 3.0)",
    )
    ap.add_argument(
        "--sphere-max",
        type=int,
        default=int(os.environ.get("FLEXAIDDS_SPHERE_MAX", "0") or "0"),
        help="Prune LOCCLF to this many spheres near LIG_ref (0=off; env FLEXAIDDS_SPHERE_MAX)",
    )
    ap.add_argument(
        "--sphere-max-dist",
        type=float,
        default=float(os.environ.get("FLEXAIDDS_SPHERE_MAX_DIST", "8.0") or "8.0"),
        help="Max distance (Å) from LIG_ref for sphere keep (default 8; env FLEXAIDDS_SPHERE_MAX_DIST)",
    )
    args = ap.parse_args(argv)

    if not args.queue_root:
        print("ERROR: --queue-root or FLEXAIDDS_QUEUE_ROOT required", file=sys.stderr)
        return 2
    q = Path(args.queue_root).expanduser().resolve()
    if not q.is_dir():
        print(f"ERROR: queue root missing: {q}", file=sys.stderr)
        return 2

    work_root = (
        Path(args.work_root).expanduser().resolve()
        if args.work_root
        else (q / "work")
    )
    arms = [a.strip() for a in args.arms.split(",") if a.strip()]
    pdbs = load_pdb_list(q, args.pilot8, args.pdb or None, args.list_file or None)

    print(f"queue={q}")
    print(f"work_root={work_root}")
    print(f"arms={arms} pdbs={len(pdbs)} pop={args.pop} gen={args.gen} R={args.restarts}")

    if args.dry_run:
        try:
            pl = resolve_processligand()
            print(f"ProcessLigand={pl}")
        except FileNotFoundError as e:
            print(f"ProcessLigand MISSING: {e}")
        for arm in arms:
            for pdb in pdbs:
                try:
                    sph = find_sphere_file(q, pdb)
                except FileNotFoundError as e:
                    sph = f"MISSING: {e}"
                print(f"  would prepare arm={arm} pdb={pdb} sphere={sph}")
        return 0

    pl = resolve_processligand()
    print(f"ProcessLigand={pl}")
    n_ok = 0
    for arm in arms:
        for pdb in pdbs:
            try:
                w = prepare_target(
                    q,
                    pdb,
                    arm,
                    work_root,
                    pop=args.pop,
                    gen=args.gen,
                    temper_override=args.temper,
                    restarts=args.restarts,
                    seed_base=args.seed_base,
                    processligand=pl,
                    force=args.force,
                    keep_hoh=args.keep_hoh,
                    keep_metals=args.keep_metals,
                    skip_clean_apo=args.skip_clean_apo,
                    skip_ligand_gate=args.skip_ligand_gate,
                    max_bond=args.max_bond,
                    sphere_max=args.sphere_max,
                    sphere_max_dist=args.sphere_max_dist,
                )
                print(f"OK {arm}/{pdb} → {w}")
                n_ok += 1
            except Exception as exc:
                print(f"FAIL {arm}/{pdb}: {exc}", file=sys.stderr)
    print(f"done: {n_ok}/{len(arms)*len(pdbs)}")
    return 0 if n_ok == len(arms) * len(pdbs) else 1


if __name__ == "__main__":
    sys.exit(main())
