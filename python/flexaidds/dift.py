"""DiFT — Discrete Fourier Transform torsional parametrization for FlexAID∆S.

Pure-Python (NumPy) mirror of the C++ ``dift`` engine (``LIB/DiFT/``). Produces
numerically identical results to the compiled core; the C++ path exists only as
a speed option. This module is always importable — it has no dependency on the
compiled ``_core`` extension.

Method
------
A dihedral's torsional potential ``V(φ)`` and its Fourier spectrum ``{Aₙ, ωₙ}``
are conjugate representations of the same object. One FFT yields two payoffs:

* **Energy**  — a truncated cosine series ``V(φ) = mean + Σ Aₙ cos(nφ − ωₙ)``.
* **Entropy** — the same spectrum gives the 1-D torsional partition function
  ``z = ⟨exp(−βV)⟩`` and hence a rigorous per-bond torsional entropy
  ``S_tors = (⟨V⟩ − F)/T`` — the ΔS contribution an entropy-driven docking
  engine should actually use, instead of a heuristic rotatable-bond count.

Truncation is **Shannon-collapse**: the spectral entropy ``H_spec`` of the
normalized power spectrum sets the effective mode count ``N_eff = exp(H_spec)``;
the ``⌈N_eff⌉`` highest-power terms are retained. No user threshold.

Reference
---------
Flores-Trujillo, Rodríguez-Segura, Amador-Bedolla & Domínguez (2026).
"Fast Fourier Transform Enables Automated Parametrization of Complex Dihedral
Potentials in All-Atom and Coarse-Grained Force Fields."
J. Chem. Inf. Model.  DOI: 10.1021/acs.jcim.6c00123

SPDX-License-Identifier: Apache-2.0
Copyright 2026 Le Bonhomme Pharma
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Sequence

import numpy as np

__all__ = [
    "kB_kcal",
    "FourierTerm",
    "TorsionalPotential",
    "TorsionalThermo",
    "DiFTEngine",
    "spectral_entropy",
    "RotatableBondTorsion",
    "TorsionalScore",
    "score_torsional",
    "make_bond_torsion",
]

# Physical constant — matches statmech / the C++ dift namespace.
kB_kcal = 0.001987206  # kcal mol⁻¹ K⁻¹

_TWO_PI = 2.0 * math.pi


def _wrap_pi(a: float) -> float:
    """Wrap an angle into (−π, π]."""
    a = math.fmod(a + math.pi, _TWO_PI)
    if a < 0.0:
        a += _TWO_PI
    return a - math.pi


# ─────────────────────────────────────────────────────────────────────────────
# Data structures
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class FourierTerm:
    """A single cosine term: ``amplitude · cos(multiplicity·φ − phase)``."""

    multiplicity: int = 0      # n — integer frequency
    amplitude: float = 0.0     # Aₙ — force constant (kcal/mol), ≥ 0
    phase: float = 0.0         # ωₙ — phase shift (radians), in (−π, π]
    power: float = 0.0         # Aₙ² — spectral power

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (f"<FourierTerm n={self.multiplicity} "
                f"A={self.amplitude:.4g} ω={self.phase:.4g}>")


@dataclass
class TorsionalPotential:
    """Analytical torsional potential ``V(φ) = mean + Σ Aₙ cos(nφ − ωₙ)``."""

    terms: List[FourierTerm] = field(default_factory=list)
    mean: float = 0.0               # A₀ — DC offset (kcal/mol)
    v_min: float = 0.0              # global minimum of V(φ) (kcal/mol)
    n_samples: int = 0              # M — input grid resolution
    r_squared: float = 0.0          # fit quality vs the target profile
    spectral_entropy: float = 0.0   # H_spec (nats) of the FULL spectrum
    effective_modes: float = 0.0    # N_eff = exp(H_spec)
    refinement_iters: int = 0       # QM–MM refinement iterations used

    def evaluate(self, phi: float) -> float:
        """Absolute potential at angle ``phi`` (radians)."""
        v = self.mean
        for t in self.terms:
            v += t.amplitude * math.cos(t.multiplicity * phi - t.phase)
        return v

    def relative(self, phi: float) -> float:
        """Potential relative to its global minimum — the form used to score."""
        return self.evaluate(phi) - self.v_min

    def sample(self, n: int) -> np.ndarray:
        """Sample the analytical potential on ``n`` points over [0, 2π)."""
        if n < 1:
            return np.empty(0, dtype=float)
        phi = np.arange(n, dtype=float) * (_TWO_PI / n)
        out = np.full(n, self.mean, dtype=float)
        for t in self.terms:
            out += t.amplitude * np.cos(t.multiplicity * phi - t.phase)
        return out

    @property
    def n_terms(self) -> int:
        return len(self.terms)


@dataclass
class TorsionalThermo:
    """Per-bond torsional thermodynamics (excess, vs a free rotor).

    A free rotor (``V ≡ 0``) gives ``entropy = 0``; a confined bond gives
    ``entropy < 0`` (an entropy loss) and ``minus_TS > 0`` (a ΔG penalty).
    """

    temperature: float = 300.0
    partition_function: float = 1.0   # z = ⟨exp(−βV)⟩
    free_energy: float = 0.0          # F = −kT ln z      (kcal/mol)
    mean_energy: float = 0.0          # ⟨V⟩               (kcal/mol)
    entropy: float = 0.0              # S_tors            (kcal mol⁻¹ K⁻¹)
    minus_TS: float = 0.0             # −T·S_tors         (kcal/mol)


# ─────────────────────────────────────────────────────────────────────────────
# Spectral Shannon entropy
# ─────────────────────────────────────────────────────────────────────────────

def spectral_entropy(spectrum: Sequence[FourierTerm]) -> float:
    """Shannon entropy ``H_spec`` (nats) of a power spectrum.

    Normalized powers ``pₙ = Aₙ²/ΣAₘ²`` form a distribution; ``H = −Σ pₙ ln pₙ``.
    ``exp(H)`` is the effective number of contributing frequencies.
    """
    total = sum(t.power for t in spectrum)
    if total <= 0.0:
        return 0.0
    h = 0.0
    for t in spectrum:
        p = t.power / total
        if p > 0.0:
            h -= p * math.log(p)
    return h


# ─────────────────────────────────────────────────────────────────────────────
# Core engine
# ─────────────────────────────────────────────────────────────────────────────

class DiFTEngine:
    """DiFT torsional parametrization engine (pure-Python / NumPy)."""

    def __init__(self, temperature_K: float = 300.0) -> None:
        self.set_temperature(temperature_K)

    def set_temperature(self, t_k: float) -> None:
        self._t = t_k if t_k > 0.0 else 300.0
        self._beta = 1.0 / (kB_kcal * self._t)

    @property
    def temperature(self) -> float:
        return self._t

    # ── forward transform ────────────────────────────────────────────────────
    def transform(self, profile: Sequence[float]) -> tuple[List[FourierTerm], float]:
        """Forward DiFT of an M-point torsional profile over [0, 2π).

        Returns ``(spectrum, mean)`` where ``spectrum`` lists terms n = 1…⌊M/2⌋.
        """
        x = np.asarray(profile, dtype=float)
        m = x.size
        if m < 2:
            raise ValueError("DiFT.transform: profile needs ≥ 2 samples")

        # rfft → X_0 … X_{M//2}; sign convention X_n = Σ x_k exp(−i2πnk/M).
        x_full = np.fft.rfft(x)
        inv_m = 1.0 / m
        mean = x_full[0].real * inv_m

        n_max = m // 2
        spectrum: List[FourierTerm] = []
        for n in range(1, n_max + 1):
            nyquist = (m % 2 == 0) and (n == n_max)
            scale = inv_m if nyquist else (2.0 * inv_m)
            amp = scale * abs(x_full[n])
            phase = -math.atan2(x_full[n].imag, x_full[n].real)
            if amp < 0.0:
                amp, phase = -amp, _wrap_pi(phase + math.pi)
            else:
                phase = _wrap_pi(phase)
            spectrum.append(FourierTerm(multiplicity=n, amplitude=amp,
                                        phase=phase, power=amp * amp))
        return spectrum, mean

    # ── parametrize with Shannon-collapse truncation ─────────────────────────
    def parametrize(self, profile: Sequence[float],
                    max_multiplicity: int = 0) -> TorsionalPotential:
        """Forward transform + Shannon-collapse spectral truncation."""
        spectrum, mean = self.transform(profile)

        h_spec = spectral_entropy(spectrum)
        n_eff = math.exp(h_spec)

        if max_multiplicity > 0:
            spectrum = [t for t in spectrum if t.multiplicity <= max_multiplicity]

        keep = max(1, min(len(spectrum), math.ceil(n_eff)))
        spectrum.sort(key=lambda t: t.power, reverse=True)
        spectrum = spectrum[:keep]
        spectrum.sort(key=lambda t: t.multiplicity)

        pot = TorsionalPotential(
            terms=spectrum,
            mean=mean,
            n_samples=len(profile),
            spectral_entropy=h_spec,
            effective_modes=n_eff,
        )
        model = pot.sample(len(profile))
        pot.r_squared = self.r_squared(np.asarray(profile, dtype=float), model)
        pot.v_min = float(np.min(pot.sample(720)))
        return pot

    # ── iterative QM–MM refinement (paper eq. 18) ────────────────────────────
    def refine(self, qm: Sequence[float], mm_initial: Sequence[float],
               lambda_: float = 0.5, r2_target: float = 0.98,
               max_iter: int = 50, max_multiplicity: int = 6) -> TorsionalPotential:
        """Damped, FFT band-limited refinement: ``V_{i+1} = V_i + λ·D_i``."""
        qm_a = np.asarray(qm, dtype=float)
        mm_a = np.asarray(mm_initial, dtype=float)
        m = qm_a.size
        if m < 2 or mm_a.size != m:
            raise ValueError("DiFT.refine: qm and mm_initial must share a "
                             "grid of ≥ 2 samples")

        correction = np.zeros(m, dtype=float)
        pot = TorsionalPotential()
        iters = 0
        r2 = 0.0
        for it in range(1, max_iter + 1):
            iters = it
            mm_current = mm_a + correction
            correction = correction + lambda_ * (qm_a - mm_current)
            pot = self.parametrize(correction, max_multiplicity)
            correction = pot.sample(m)
            r2 = self.r_squared(qm_a, mm_a + correction)
            if r2 >= r2_target:
                break

        pot.r_squared = r2
        pot.refinement_iters = iters
        return pot

    # ── Boltzmann inversion of a CG dihedral histogram (paper eq. 20) ────────
    def boltzmann_invert(self, histogram: Sequence[float]) -> np.ndarray:
        """Invert a dihedral-angle histogram into an energy profile.

        ``E(φ) = −kT ln p(φ)``, shifted so the minimum is 0. No Jacobian:
        dihedral angles are uniformly distributed for a free rotor.
        """
        h = np.asarray(histogram, dtype=float)
        if h.size < 2:
            raise ValueError("DiFT.boltzmann_invert: need ≥ 2 bins")
        if np.any(h < 0.0):
            raise ValueError("DiFT.boltzmann_invert: negative count")
        total = float(h.sum())
        if total <= 0.0:
            raise ValueError("DiFT.boltzmann_invert: empty histogram")

        kt = kB_kcal * self._t
        energy = np.full(h.size, np.nan, dtype=float)
        nonzero = h > 0.0
        energy[nonzero] = -kt * np.log(h[nonzero] / total)
        e_max = float(np.nanmax(energy))
        energy[~nonzero] = e_max          # empty bins → capped well wall
        energy -= float(np.nanmin(energy))  # shift minimum to 0
        return energy

    # ── per-bond torsional thermodynamics ────────────────────────────────────
    def thermodynamics(self, pot: TorsionalPotential) -> TorsionalThermo:
        """Excess torsional thermodynamics from a parametrized potential."""
        n_fine = 1440
        phi = np.arange(n_fine, dtype=float) * (_TWO_PI / n_fine)
        v = np.full(n_fine, pot.mean, dtype=float)
        for t in pot.terms:
            v += t.amplitude * np.cos(t.multiplicity * phi - t.phase)
        v -= pot.v_min                       # energies relative to the minimum

        w = np.exp(-self._beta * v)
        z_sum = float(w.sum())
        z = z_sum / n_fine                   # ⟨exp(−βV)⟩; 1 for a free rotor
        mean_e = float((w * v).sum() / z_sum) if z_sum > 0.0 else 0.0
        f = -kB_kcal * self._t * math.log(z)

        th = TorsionalThermo(temperature=self._t, partition_function=z,
                             free_energy=f, mean_energy=mean_e)
        th.entropy = (mean_e - f) / self._t
        th.minus_TS = -self._t * th.entropy
        return th

    # ── helpers ──────────────────────────────────────────────────────────────
    @staticmethod
    def circular_mean(angles: Sequence[float]) -> float:
        """Directional mean of phase offsets — ``atan2(Σ sin, Σ cos)``."""
        a = np.asarray(angles, dtype=float)
        s = float(np.sin(a).sum())
        c = float(np.cos(a).sum())
        if s == 0.0 and c == 0.0:
            return 0.0
        return math.atan2(s, c)

    @staticmethod
    def r_squared(observed: np.ndarray, model: np.ndarray) -> float:
        """Coefficient of determination R² between two equal-length profiles."""
        obs = np.asarray(observed, dtype=float)
        mod = np.asarray(model, dtype=float)
        if obs.size == 0 or mod.size != obs.size:
            return 0.0
        ss_res = float(np.sum((obs - mod) ** 2))
        ss_tot = float(np.sum((obs - obs.mean()) ** 2))
        if ss_tot <= 0.0:
            return 1.0 if ss_res <= 1e-18 else 0.0
        return 1.0 - ss_res / ss_tot


# ─────────────────────────────────────────────────────────────────────────────
# GA scoring adapter — mirrors LIB/DiFT/DiFTGAAdapter.h
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class RotatableBondTorsion:
    """A ligand rotatable bond and its DiFT-parametrized potential."""

    potential: TorsionalPotential
    gene_index: int = -1


@dataclass
class TorsionalScore:
    """Decomposed torsional contribution to the binding free energy."""

    energy: float = 0.0     # Σ V_tors,b relative to each well minimum
    minus_TS: float = 0.0   # Σ −T·S_tors,b — torsional confinement penalty
    n_bonds: int = 0

    def total(self) -> float:
        """Free-energy contribution the GA adds to the CF score."""
        return self.energy + self.minus_TS


def score_torsional(bonds: Sequence[RotatableBondTorsion],
                    dihedral_angles_rad: Sequence[float],
                    temperature_K: float = 300.0) -> TorsionalScore:
    """Score a pose's torsional state — one FFT-derived potential per bond."""
    score = TorsionalScore()
    if len(bonds) != len(dihedral_angles_rad):
        return score  # caller contract violated → zero contribution
    engine = DiFTEngine(temperature_K)
    for bond, phi in zip(bonds, dihedral_angles_rad):
        score.energy += bond.potential.relative(phi)
        score.minus_TS += engine.thermodynamics(bond.potential).minus_TS
        score.n_bonds += 1
    return score


def make_bond_torsion(profile: Sequence[float], gene_index: int,
                      temperature_K: float = 300.0,
                      max_multiplicity: int = 6) -> RotatableBondTorsion:
    """Parametrize a raw torsional profile into a :class:`RotatableBondTorsion`."""
    engine = DiFTEngine(temperature_K)
    return RotatableBondTorsion(
        potential=engine.parametrize(profile, max_multiplicity),
        gene_index=gene_index,
    )
