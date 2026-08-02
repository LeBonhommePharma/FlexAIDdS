"""P0 PoseBusters fixtures — prove real PB 0.6.5 discriminates good vs broken geometry.

The 0/85 failure mode is a harness that schema-rejects EVERY pose (so nothing can
pass) OR silently passes everything. This test proves PB 0.6.5, as installed,
returns PASS on a valid molecule and FAIL on a physically broken one.

Run: python tests/p0_claim_contract/test_posebusters_fixtures.py
"""
import sys

# rdkit and posebusters are optional at collection time. Importing them
# unguarded aborts the ENTIRE pytest session with a collection error, not just
# this module. As a script this still raises normally; under pytest the module
# is skipped when the dependency is absent.
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from posebusters import PoseBusters
except ImportError as _exc:  # pragma: no cover - environment-dependent
    if __name__ == "__main__":
        raise
    import pytest
    pytest.skip(
        f"optional dependency unavailable: {_exc.name}", allow_module_level=True
    )

failures = []
def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        failures.append(name)

# --- clean molecule: RDKit-embedded, MMFF-optimized (valid geometry) ---------
m = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1C(=O)NCC(=O)O"))  # hippuric-acid-like
AllChem.EmbedMolecule(m, randomSeed=1)
AllChem.MMFFOptimizeMolecule(m)

# --- broken molecule: same graph, two atoms forced to collide -----------------
mb = Chem.Mol(m)
conf = mb.GetConformer()
p0 = conf.GetAtomPosition(0)
conf.SetAtomPosition(1, p0)  # atom 1 coincident with atom 0 -> bond-length/clash fail

pb = PoseBusters(config="mol")
df_ok = pb.bust([m], None, None)
df_bad = pb.bust([mb], None, None)

# every check column is a bool; a valid mol passes (nearly) all, broken fails some
ok_pass = int(df_ok.sum(axis=1).iloc[0])
ok_total = df_ok.shape[1]
bad_pass = int(df_bad.sum(axis=1).iloc[0])
print(f"  clean pose: {ok_pass}/{ok_total} checks pass")
print(f"  broken pose: {bad_pass}/{ok_total} checks pass")

check("PB returns a real per-check table (not empty/schema-reject)", ok_total >= 8)
check("clean molecule passes ALL checks", ok_pass == ok_total)
check("broken molecule FAILS at least one check", bad_pass < ok_total)
check("broken fails specifically on geometry (bond length or clash)",
      not bool(df_bad.get("bond_lengths", [True]).iloc[0]) if "bond_lengths" in df_bad.columns
      else bad_pass < ok_total)

# report which checks the broken mol failed (diagnostic)
if bad_pass < ok_total:
    failed_checks = [c for c in df_bad.columns if not bool(df_bad[c].iloc[0])]
    print(f"  broken pose failed: {failed_checks}")

print(f"\n{'ALL PASS' if not failures else f'{len(failures)} FAILURES: {failures}'}")

# Only exit the interpreter when run as a script. Under pytest this module is
# imported during collection, and a module-level sys.exit() aborts the whole
# session with INTERNALERROR before any test runs -- see the guard below.
if __name__ == "__main__":
    sys.exit(1 if failures else 0)


def test_posebusters_fixtures():
    """Expose the script's checks to pytest so a failure is reported, not exited."""
    assert not failures, f"{len(failures)} failures: {failures}"
