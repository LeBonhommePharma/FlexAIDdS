"""Shared FlexAID∆S PyMOL plugin session state.

Both the legacy ``flexaids_load`` path (``visualization``) and the
read-only ``flexaids_load_results`` path (``results_adapter``) write into
this module so entropy heatmaps, animation, ITC plots, and the GUI always
see a consistent loaded result set.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclass
class PluginSession:
    """Mutable plugin session for the currently loaded docking results."""

    result: Any = None  # Optional[DockingResult]
    objects: Dict[int, List[str]] = field(default_factory=dict)
    prefix: str = "flexaids"
    temperature_K: float = 300.0
    output_dir: Optional[Path] = None
    # Legacy visualization bookkeeping: mode_name -> _ModeRecord
    mode_records: Dict[str, Any] = field(default_factory=dict)

    def clear(self) -> None:
        """Drop all loaded result state (does not delete PyMOL objects)."""
        self.result = None
        self.objects.clear()
        self.mode_records.clear()
        self.output_dir = None
        self.prefix = "flexaids"
        self.temperature_K = 300.0

    @property
    def n_modes(self) -> int:
        if self.result is None:
            return 0
        return getattr(self.result, "n_modes", len(getattr(self.result, "binding_modes", [])))

    def get_mode(self, mode_id: int) -> Any:
        """Return BindingModeResult for *mode_id*, or None."""
        if self.result is None:
            return None
        mid = int(mode_id)
        for mode in self.result.binding_modes:
            if mode.mode_id == mid:
                return mode
        return None

    def object_names_for(self, mode_id: int) -> List[str]:
        return list(self.objects.get(int(mode_id), []))


# Process-wide singleton used by all plugin modules.
SESSION = PluginSession()


def object_name(prefix: str, mode_id: int, pose_rank: int) -> str:
    """Canonical PyMOL object name for one pose."""
    return f"{prefix}_mode{int(mode_id)}_pose{int(pose_rank)}"


def group_name(prefix: str, mode_id: int) -> str:
    """Canonical PyMOL group name for one binding mode."""
    return f"{prefix}_mode{int(mode_id)}"


def mode_label(mode_id: int) -> str:
    """Legacy mode name used by flexaids_show_ensemble / flexaids_thermo."""
    return f"mode{int(mode_id)}"
