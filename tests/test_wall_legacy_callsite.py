"""Exercise vcfunction's actual wall dispatch without linking the docking engine.

The extracted block consumes the precomputed Ewall, which may already contain
the upstream softcore correction. Keeping that value is part of legacy parity.
"""

from pathlib import Path
import os
import shlex
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def wall_dispatch(tmp_path_factory):
    source = (ROOT / "LIB" / "vcfunction.cpp").read_text()
    start = "// ── PHYSICAL WALL (default since the CF steric-asymptotics fix)"
    end = "// Structural-presence proof. A gate that is set but never reaches a"
    assert source.count(start) == source.count(end) == 1
    block = source.split(start, 1)[1].split(end, 1)[0]
    # Discard the remainder of the opening comment, retain the actual dispatch.
    block = block.split("\n", 1)[1]
    compiler = shlex.split(os.environ.get("CXX", "c++"))
    if not compiler or not shutil.which(compiler[0]):
        pytest.skip("A C++ compiler is required for the wall call-site test")
    work = tmp_path_factory.mktemp("wall_legacy_callsite")
    harness = work / "wall_dispatch.cpp"
    harness.write_text(
        '#include "soft_wall.h"\n'
        "#include <cstdlib>\n#include <iomanip>\n#include <iostream>\n"
        "int main(int argc, char** argv) {\n"
        "  if (argc != 11) return 2;\n"
        "  const double d = std::strtod(argv[1], nullptr);\n"
        "  const double cr = std::strtod(argv[2], nullptr);\n"
        "  struct Config { float soft_wall_cutoff; };\n"
        "  const Config config{std::strtof(argv[3], nullptr)};\n"
        "  const Config* FA = &config;\n"
        "  const double Ewall = std::strtod(argv[4], nullptr);\n"
        "  const bool wal_legacy_capped = std::atoi(argv[5]);\n"
        "  const bool wal_flex_contact = std::atoi(argv[6]);\n"
        "  const double wal_cap_flex = std::strtod(argv[7], nullptr);\n"
        "  const bool wal_coercive = std::atoi(argv[8]);\n"
        "  const double wal_stiff = std::strtod(argv[9], nullptr);\n"
        "  const bool wal_c1 = std::atoi(argv[10]);\n"
        "  double Ewall_fitness;\n"
        + block
        + '\n  std::cout << std::setprecision(17) << Ewall_fitness << "\\n";\n'
        "}\n"
    )
    executable = work / "wall_dispatch"
    subprocess.run(
        [*compiler, "-std=c++17", "-O2", f"-I{ROOT / 'LIB'}", str(harness),
         "-o", str(executable)],
        check=True, capture_output=True, text=True, timeout=60,
    )

    def evaluate(*, cutoff, precomputed_wall, flex=False, cap=0,
                 coercive=False, legacy=True, c1=False):
        result = subprocess.run(
            [str(executable), "1.5", "3.5", str(cutoff), str(precomputed_wall),
             str(int(legacy)), str(int(flex)), str(cap), str(int(coercive)),
             "50", str(int(c1))],
            check=True, capture_output=True, text=True, timeout=10,
        )
        return float(result.stdout)

    return evaluate


@pytest.mark.parametrize(
    "flex,cap,coercive,expected",
    [(True, 100, False, 100), (True, 250, False, 200),
     (True, 0, False, 200), (False, 100, False, 50),
     (False, 100, True, 200)],
)
def test_legacy_soft_wall_cap(wall_dispatch, flex, cap, coercive, expected):
    # At overlap 2 A, the uncapped historical quadratic is 50 * 2**2 = 200.
    assert wall_dispatch(
        cutoff=0.4, precomputed_wall=999, flex=flex, cap=cap, coercive=coercive,
    ) == pytest.approx(expected)


@pytest.mark.parametrize("legacy", [False, True])
@pytest.mark.parametrize(
    "precomputed_wall,flex,cap,expected",
    [(200, True, 100, 100), (200, True, 0, 50), (200, False, 100, 50),
     (7.25, False, 0, 7.25), (75, True, 100, 75)],
)
def test_zero_cutoff_preserves_precomputed_wall(
    wall_dispatch, legacy, precomputed_wall, flex, cap, expected,
):
    # Distinct supplied values catch accidental recomputation of raw r^-12:
    # upstream softcore corrections must survive the zero-cutoff path.
    assert wall_dispatch(
        cutoff=0, precomputed_wall=precomputed_wall, flex=flex, cap=cap,
        legacy=legacy,
    ) == expected


@pytest.mark.parametrize("c1", [False, True])
def test_physical_wall_still_applies_flex_cap(wall_dispatch, c1):
    assert wall_dispatch(
        cutoff=0.4, precomputed_wall=999, flex=True, cap=100,
        legacy=False, c1=c1,
    ) == 100
