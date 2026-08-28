"""Terminal colour for the ``flexaidds`` CLI, matching the design system.

This is the Python mirror of ``LIB/TuiColor.h``. Same seven key colors, same
bindings — a key color is never reassigned to a different quantity, because
that binding *is* the system::

    mint       #45E0A8   ΔH · enthalpy · brand primary
    violet     #8B5CF6   ΔS · configurational entropy
    tangerine  #FF9300   ΔG · free energy · results
    firetruck  #F5232B   T  · temperature · and failure
    aqua       #00A2FF   ΔS_vib · vibrational · tENCoM
    strawberry #FF2F92   receptor · pocket · section eyebrows · warnings
    magnesium  #DCDCE4   apo baseline · reference

The SERIES ramp is ordered by energy along the binding coordinate, not by
hue, so a staged run reads as a reaction path rather than a log.

Why nothing here can corrupt a file
-----------------------------------
Every accessor returns ``""`` unless stdout is an interactive terminal, so
``python -m flexaidds results/ > summary.txt`` and any pipe produce exactly
the bytes they produced before this module existed. ``--json`` additionally
calls :func:`set_enabled(False)` — machine-readable output is never coloured
even on a TTY, because a JSON consumer on the other end of a terminal is
still a JSON consumer.

Honoured, in precedence order:

* ``FLEXAIDDS_NO_COLOR`` / ``NO_COLOR`` set  → off (https://no-color.org)
* ``TERM=dumb`` or unset                     → off
* ``FORCE_COLOR`` / ``CLICOLOR_FORCE`` set   → on, even when piped
* stdout is not a TTY                        → off
"""

from __future__ import annotations

import os
import sys
from typing import Optional

__all__ = [
    "enabled", "set_enabled", "reset", "bold", "dim",
    "mint", "violet", "tangerine", "firetruck", "aqua", "strawberry",
    "magnesium", "fg", "muted", "failtext",
    "dH", "dS", "dG", "T", "dSvib", "pocket", "apo",
    "series", "brand", "rule", "stage", "equation",
    "ok", "warn", "fail", "converged", "tag", "kv",
]

# (r, g, b, xterm-256 fallback)
_MINT       = (0x45, 0xE0, 0xA8, 79)
_VIOLET     = (0x8B, 0x5C, 0xF6, 99)
_TANGERINE  = (0xFF, 0x93, 0x00, 208)
_FIRETRUCK  = (0xF5, 0x23, 0x2B, 196)
_AQUA       = (0x00, 0xA2, 0xFF, 39)
_STRAWBERRY = (0xFF, 0x2F, 0x92, 198)
_MAGNESIUM  = (0xDC, 0xDC, 0xE4, 253)
_FG         = (0xE4, 0xE3, 0xF5, 189)
_MUTED      = (0x8D, 0x8C, 0xB0, 103)
_FAILTEXT   = (0xFF, 0x6B, 0x6B, 203)

_state: Optional[str] = None  # None = undetected, "" = off, "true"/"256" = on


def _detect() -> str:
    if os.environ.get("FLEXAIDDS_NO_COLOR") or "NO_COLOR" in os.environ:
        return ""
    term = os.environ.get("TERM")
    if not term or term == "dumb":
        return ""
    forced = bool(os.environ.get("FORCE_COLOR") or os.environ.get("CLICOLOR_FORCE"))
    if not forced:
        try:
            if not sys.stdout.isatty():
                return ""
        except Exception:
            # A stdout with no isatty (captured, replaced) is not a terminal.
            return ""
    ct = os.environ.get("COLORTERM", "")
    if "truecolor" in ct or "24bit" in ct:
        return "true"
    return "256"


def _mode() -> str:
    global _state
    if _state is None:
        _state = _detect()
    return _state


def enabled() -> bool:
    return _mode() != ""


def set_enabled(on: bool) -> None:
    """Force colour off (``--json``, ``--no-color``) or back on."""
    global _state
    _state = (_detect() or "256") if on else ""


def _seq(hue) -> str:
    m = _mode()
    if not m:
        return ""
    r, g, b, x = hue
    return f"\033[38;2;{r};{g};{b}m" if m == "true" else f"\033[38;5;{x}m"


def reset() -> str: return "\033[0m" if enabled() else ""
def bold() -> str:  return "\033[1m" if enabled() else ""
def dim() -> str:   return "\033[2m" if enabled() else ""

def mint() -> str:       return _seq(_MINT)
def violet() -> str:     return _seq(_VIOLET)
def tangerine() -> str:  return _seq(_TANGERINE)
def firetruck() -> str:  return _seq(_FIRETRUCK)
def aqua() -> str:       return _seq(_AQUA)
def strawberry() -> str: return _seq(_STRAWBERRY)
def magnesium() -> str:  return _seq(_MAGNESIUM)
def fg() -> str:         return _seq(_FG)
def muted() -> str:      return _seq(_MUTED)
def failtext() -> str:   return _seq(_FAILTEXT)

# Quantity aliases — prefer these at call sites so the colour→quantity
# binding stays legible in the source itself.
def dH() -> str:    return mint()
def dS() -> str:    return violet()
def dG() -> str:    return tangerine()
def T() -> str:     return firetruck()
def dSvib() -> str: return aqua()
def pocket() -> str: return strawberry()
def apo() -> str:   return magnesium()


def series(step: int) -> str:
    """Series ramp, ordered by energy along the binding coordinate.

    Firetruck is deliberately absent: it is a scalar and a failure signal,
    never a data class.
    """
    return {
        1: magnesium(),   # apo baseline
        2: violet(),      # unbound · ΔS dominates
        3: strawberry(),  # first pocket contact
        4: aqua(),        # rigidification · ΔS_vib
        5: mint(),        # contacts formed · ΔH
        6: tangerine(),   # converged · ΔG
    }.get(step, muted())


def brand() -> str:
    """``FlexAID∆S`` — mint mark, tangerine ∆S, as everywhere else."""
    return f"{bold()}{mint()}FlexAID{reset()}{tangerine()}∆S{reset()}"


def rule(title: str, width: int = 58) -> str:
    pad = max(0, width - len(title))
    return f"{violet()}──{reset()} {bold()}{title}{reset()} {violet()}{'-' * pad}{reset()}"


def stage(n: int, total: int, label: str, detail: str = "") -> str:
    """``[n/total] label            detail`` — one step of the binding
    coordinate, coloured by its place on the series ramp."""
    return (f"  {muted()}[{n}/{total}]{reset()} "
            f"{series(n)}{label}{reset()}{' ' * max(0, 28 - len(label))} "
            f"{muted()}{detail}{reset()}")


def equation(temperature_K: Optional[float] = None) -> str:
    """The identity, read term by term in colour."""
    s = (f"  {dG()}ΔG{reset()} = {dH()}ΔH{reset()} − {T()}T{reset()}{dS()}ΔS{reset()}"
         f" − {T()}T{reset()}{dSvib()}ΔS_vib{reset()}")
    if temperature_K is not None:
        s += f"        {T()}T = {temperature_K} K{reset()}"
    return s


def ok(msg: str) -> str:        return f"  {mint()}●{reset()} {fg()}{msg}{reset()}"
def warn(msg: str) -> str:      return f"  {strawberry()}● WARN{reset()} {fg()}{msg}{reset()}"
def fail(msg: str) -> str:      return f"  {failtext()}● FAIL{reset()} {fg()}{msg}{reset()}"
def converged(msg: str) -> str: return f"  {tangerine()}● CONVERGED{reset()} {muted()}{msg}{reset()}"


def tag(name: str, hue: Optional[str] = None) -> str:
    """A ``[TAG]`` eyebrow. Strawberry by default, matching the site."""
    return f"{hue or strawberry()}[{name}]{reset()}"


def kv(key: str, value, hue: Optional[str] = None) -> str:
    """``key=value`` with the value carrying its quantity's colour."""
    return f"{muted()}{key}={reset()}{hue or fg()}{value}{reset()}"
