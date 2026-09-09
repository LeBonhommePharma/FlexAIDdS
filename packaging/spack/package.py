# Spack recipe for FlexAID∆S (not yet in builtin Spack; copy into a repo overlay).
# SPDX-License-Identifier: Apache-2.0
from spack.package import *


class Flexaidds(CMakePackage):
    """Entropy-aware molecular docking engine (native FlexAIDdS + tENCoM)."""

    homepage = "https://github.com/LeBonhommePharma/FlexAIDdS"
    url = "https://github.com/LeBonhommePharma/FlexAIDdS/archive/refs/tags/v2.0.3.tar.gz"
    git = "https://github.com/LeBonhommePharma/FlexAIDdS.git"

    maintainers("LeBonhommePharma")
    license("Apache-2.0")

    version("main", branch="main")
    version(
        "2.0.3",
        sha256="6c8442fc672a127db354ff3b6e08a2252e8c921372d902d062ecbf4296aef186",
    )

    variant("python", default=False, description="Install flexaidds Python package")
    variant("metal", default=False, description="Metal GPU (macOS only)")

    depends_on("cmake@3.28:", type="build")
    depends_on("ninja", type="build")
    depends_on("eigen@3.4:")
    depends_on("python@3.9:", when="+python", type=("build", "run"))

    def cmake_args(self):
        args = [
            self.define("BUILD_TESTING", False),
            self.define("BUILD_PYTHON_BINDINGS", False),
            self.define("FLEXAIDS_USE_CUDA", False),
            self.define_from_variant("FLEXAIDS_USE_METAL", "metal"),
            self.define("FLEXAIDS_GIT_COMMIT_OVERRIDE", str(self.spec.version)),
        ]
        return args
