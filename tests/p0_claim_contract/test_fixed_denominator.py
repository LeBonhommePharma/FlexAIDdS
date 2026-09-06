"""P0 claim-contract fixtures — fail-closed proof of the aggregator invariants.

Run: python tests/p0_claim_contract/test_fixed_denominator.py
Exit 0 = all invariants hold; nonzero = contract broken.

Covers Codex P0 kill/promote gates:
  1. Fixed 85-target denominator: missing/dropped targets count as FAILURES,
     never removed from the denominator (anti-inflation).
  2. Same-pose hash receipts: a hash mismatch fails claim admission (closed).
  3. Seed/native-seed gates fail closed on missing columns.
"""
import sys, importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "agg", ROOT / "scripts" / "aggregate_claim_metrics.py"
)
agg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(agg)

PIN = "9dc93717dfed0698006d88dd6a9627bc"
failures = []


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        failures.append(name)


def strict_row(pid, sha="a" * 64):
    """A fully claim-ready, strict-success row."""
    return {
        "pdb_id": pid, "seed_echo": "0", "native_pose_seeded": "0",
        "matrix_md5": PIN, "protocol_claim_eligible": "1", "claim_ready": "1",
        "rmsd_to_crystal": "1.0", "pb_pass": "1", "success_pb": "1",
        "success_rmsd": "1", "score_pose_consistent": "1", "score_pose_delta": "0",
        "posebusters_input_sha256": "b" * 64,
        "pb_backend": "bust_cli", "tencom_status": "ok", "eigen_status": "ok",
        "pb_ran": "1", "pb_n_pass": "27", "pb_n_fail": "0", "pb_n_checks": "27",
        "eigen_n_modes": "1", "elected_H_vib": "-1.5",
        "pose_sha256": sha, "rmsd_pose_sha256": sha,
        "posebusters_pose_sha256": sha, "tencom_pose_sha256": sha,
    }


manifest_codes, _ = agg.load_target_manifest()
print(f"manifest loaded: {len(manifest_codes) if manifest_codes else 0} targets")
check("manifest present with 85 targets", bool(manifest_codes) and len(manifest_codes) == 85)

# --- Invariant 1: fixed denominator, missing targets = failures --------------
# Only 10 of 85 present, all strict successes. Naive would report 10/10=100%.
# Correct: 10/85 = 11.76%.
some = manifest_codes[:10]
rows = [strict_row(pid) for pid in some]
rep = agg.aggregate_rows(rows, PIN, "test", fixed_denominator=True)
check("denominator is 85 not 10", rep["N_denominator"] == 85)
check("STRICT rate = 10/85 (not 10/10)",
      abs(rep["metrics"]["STRICT"]["rate"] - 10/85) < 1e-9)
check("N_missing = 75", rep["N_missing_from_manifest"] == 75)
check("headline N = 85", rep["headline"]["N"] == 85)

# --- Invariant 2: hash mismatch fails closed ---------------------------------
bad = strict_row("1G9V")
bad["posebusters_pose_sha256"] = "b" * 64  # different pose validated by PB
rep2 = agg.aggregate_rows([bad], PIN, "test", fixed_denominator=True)
check("hash-mismatch row is NOT claim-eligible", rep2["N_claim"] == 0)
check("hash-mismatch STRICT n = 0", rep2["metrics"]["STRICT"]["n"] == 0)

# --- Invariant 3: seed gates fail closed on missing columns ------------------
noseed = {"pdb_id": "1G9V", "matrix_md5": PIN, "claim_ready": "1",
          "rmsd_to_crystal": "1.0"}  # no seed_echo / native_pose_seeded columns
ok, reasons = agg.is_claim_eligible(noseed, PIN)
check("missing seed_echo/native_pose_seeded fails closed", not ok)
check("  reasons name both gates",
      any("seed_echo" in r for r in reasons) and
      any("native_pose_seeded" in r for r in reasons))

# --- Invariant 4: extra rows can't shrink denominator below 85 ---------------
extra = [strict_row(pid) for pid in manifest_codes] + [strict_row("9XXX")]
rep3 = agg.aggregate_rows(extra, PIN, "test", fixed_denominator=True)
check("denominator stays 85 with an off-manifest extra row",
      rep3["N_denominator"] == 85)

# --- Invariant 5: off-manifest success must not inflate STRICT numerator ------
check("STRICT n stays 85 (off-manifest extra does not inflate numerator)",
      rep3["metrics"]["STRICT"]["n"] == 85)
check("STRICT ids exclude off-manifest 9XXX",
      "9XXX" not in {str(x).upper() for x in rep3["metrics"]["STRICT"]["ids"]})

print(f"\n{'ALL PASS' if not failures else f'{len(failures)} FAILURES: {failures}'}")

# Only exit the interpreter when run as a script. Under pytest this module is
# imported during collection, and a module-level sys.exit() aborts the whole
# session with INTERNALERROR before any test runs -- see the guard below.
if __name__ == "__main__":
    sys.exit(1 if failures else 0)


def test_fixed_denominator_invariants():
    """Expose the script's checks to pytest so a failure is reported, not exited."""
    assert not failures, f"{len(failures)} failures: {failures}"
