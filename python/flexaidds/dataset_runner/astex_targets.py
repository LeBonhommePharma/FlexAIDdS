"""Astex benchmark target tables — parity with LIB/DatasetRunner.cpp."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

# Verdonk / Gaudreault expanded 85-structure Astex Diverse set (DatasetRunner.cpp).
ASTEX_DIVERSE_CODES: Tuple[str, ...] = (
    "1G9V", "1GM8", "1GPK", "1HNN", "1HP0", "1HQ2", "1IA1", "1IGJ",
    "1J3J", "1JD0", "1JJE", "1K3U", "1KE5", "1KZK", "1L2S", "1L7F",
    "1LPZ", "1M2Z", "1MEH", "1MQ6", "1N1M", "1N2J", "1N2V", "1N46",
    "1NAV", "1OF1", "1OF6", "1OPK", "1OQ5", "1OWE", "1P2Y", "1P62",
    "1PMN", "1Q1G", "1Q41", "1Q4G", "1R1H", "1R55", "1R58", "1R9O",
    "1S19", "1S3V", "1SG0", "1SJ0", "1SQ5", "1T40", "1T46", "1T9B",
    "1TT1", "1TW6", "1TZ8", "1U1C", "1U4D", "1UML", "1UNL", "1UOU",
    "1V0P", "1V48", "1V4S", "1VCJ", "1W1P", "1W2G", "1X8X", "1XM6",
    "1XOZ", "1Y6B", "1Y6R", "1YGC", "1YQY", "1YV3", "1YVF", "1YWR",
    "1Z95", "2BM2", "2BR1", "2BSM", "2BYS", "2C3I", "2CET", "2CGR",
    "2D3U", "2GBP", "2HB1", "2HR7", "2J62",
)


@dataclass(frozen=True)
class AstexNonNativeFamily:
    """One protein family in the Astex non-native cross-docking set."""

    name: str
    native_pdb: str
    alternative_pdbs: Tuple[str, ...]


# Representative tier-2 YAML families + native PDB codes (DatasetRunner.cpp).
ASTEX_NONNATIVE_FAMILIES: Tuple[AstexNonNativeFamily, ...] = (
    AstexNonNativeFamily("ACE", "1G9V", ("1EVE", "1GQR", "1QTI", "2ACE")),
    AstexNonNativeFamily("CA2", "1V4S", ("1A42", "1AM6", "2CA2", "2CBA")),
    AstexNonNativeFamily("CDK2", "1KE5", ("1AQ1", "1B38", "1HCK", "1KE6")),
    AstexNonNativeFamily("ER", "1SJ0", ("1ERE", "1A52", "2POG", "1NDE")),
    AstexNonNativeFamily("HIVPR", "1HQ2", ("1A30", "1HVR", "1HXW", "1HSG")),
    AstexNonNativeFamily("HSP90", "1UY6", ("1BYQ", "1UY8", "2FWZ", "2XAB")),
    AstexNonNativeFamily("PPARg", "1K74", ("1FM6", "1PRG", "2P4Y", "2PRG")),
    AstexNonNativeFamily("PTP1B", "1Q1G", ("1BZH", "1G1F", "1OEM", "2B07")),
    AstexNonNativeFamily("REN", "1R9O", ("1RNE", "2REN", "1BIL", "3G6Z")),
    AstexNonNativeFamily("THR", "1TT1", ("1PPB", "1A2C", "1ABJ", "2CF8")),
)

ASTEX_NONNATIVE_BY_NAME: Dict[str, AstexNonNativeFamily] = {
    f.name: f for f in ASTEX_NONNATIVE_FAMILIES
}


def lookup_nonnative_family(name: str) -> Optional[AstexNonNativeFamily]:
    """Case-insensitive family lookup (ACE, ace, etc.)."""
    return ASTEX_NONNATIVE_BY_NAME.get(name) or ASTEX_NONNATIVE_BY_NAME.get(name.upper())


def parse_crossdock_entry_id(entry_id: str) -> Optional[Tuple[str, str]]:
    """Parse catalog entry ids like ``1G9V_1EVE`` → (native_pdb, receptor_pdb)."""
    if "_" not in entry_id:
        return None
    native, receptor = entry_id.split("_", 1)
    if len(native) != 4 or len(receptor) != 4:
        return None
    return native.upper(), receptor.upper()