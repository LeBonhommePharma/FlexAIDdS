class Flexaidds < Formula
  desc "Entropy-driven molecular docking engine (FlexAID + ΔS thermodynamic analysis)"
  homepage "https://github.com/LeBonhommePharma/FlexAIDdS"
  url "https://github.com/LeBonhommePharma/FlexAIDdS/archive/refs/tags/v2.0.0.tar.gz"
  sha256 "31488777f361bb02cb29df7a5e65e62cea15434e451bb574a14b10e702ae21e3"
  license "Apache-2.0"
  head "https://github.com/LeBonhommePharma/FlexAIDdS.git", branch: "master"

  livecheck do
    url :stable
    strategy :github_latest
  end

  depends_on "cmake" => :build
  depends_on "ninja" => :build
  depends_on "eigen"
  depends_on "libomp" if OS.mac?

  def install
    # Build the main tools. Focus on the fast optimized FlexAIDdS binary,
    # the tENCoM tool, and legacy FlexAID for compatibility.
    args = std_cmake_args + %W[
      -GNinja
      -DCMAKE_BUILD_TYPE=Release
      -DBUILD_TESTING=OFF
      -DBUILD_PYTHON_BINDINGS=OFF
      -DFLEXAIDS_USE_CUDA=OFF
      -DFLEXAIDS_USE_METAL=#{OS.mac? ? "ON" : "OFF"}
      -DFLEXAIDS_USE_OPENMP=ON
      -DBUILD_FLEXAIDDS_FAST=ON
      -DENABLE_TENCOM_TOOL=ON
    ]

    if OS.mac?
      # Help CMake find Homebrew's libomp (keg-only)
      libomp = formula_opt_prefix("libomp")
      args += %W[
        -DOpenMP_C_FLAGS=-Xpreprocessor\ -fopenmp\ -I#{libomp}/include
        -DOpenMP_C_LIB_NAMES=omp
        -DOpenMP_CXX_FLAGS=-Xpreprocessor\ -fopenmp\ -I#{libomp}/include
        -DOpenMP_CXX_LIB_NAMES=omp
        -DOpenMP_omp_LIBRARY=#{libomp}/lib/libomp.dylib
      ]
    end

    system "cmake", "-S", ".", "-B", "build", *args
    system "cmake", "--build", "build", "--parallel"

    # The build system stages runtime data files (MC_*.dat, AMINO.def, etc.)
    # next to the built binaries via POST_BUILD rules. Install them together
    # so relative-path lookup in the executables succeeds out of the box.
    bin.install "build/FlexAIDdS" if File.exist?("build/FlexAIDdS")
    bin.install "build/tENCoM" if File.exist?("build/tENCoM")
    bin.install "build/FlexAID" if File.exist?("build/FlexAID")

    # Install co-located data files (staged by the build).
    # These are looked up relative to the executable directory.
    bin.install Dir["build/MC_*.dat"]
    bin.install "build/AMINO.def" if File.exist?("build/AMINO.def")
    bin.install "build/NUCLEOTIDES.def" if File.exist?("build/NUCLEOTIDES.def")
    bin.install Dir["build/FA_matrix*.dat"]

    # Optional useful tools
    bin.install "build/flexaids_process_ligand" if File.exist?("build/flexaids_process_ligand")
    bin.install "build/tencom_entropy_diff" if File.exist?("build/tencom_entropy_diff")

    # NOTE: The Python package "flexaidds" is best installed via pip from GitHub
    # until the first public PyPI release (see docs/INSTALLATION.md).
    # This formula focuses on the high-performance native CLI tools.
  end

  def caveats
    <<~EOS
      The FlexAIDdS native tools (FlexAIDdS, tENCoM, FlexAID) have been installed.

      Runtime data files (MC matrices, AMINO.def, etc.) are installed alongside the
      binaries so they are found automatically.

      For the Python package (analysis, results loader, thermodynamics).
      Not yet on public PyPI — install from GitHub:
        pip install "git+https://github.com/LeBonhommePharma/FlexAIDdS.git#subdirectory=python"

      To upgrade later:
        python -m flexaidds --self-update
        # or re-run the git+https pip install above

      To rebuild the latest development tip of this formula:
        brew reinstall --HEAD flexaidds

      Example usage:
        FlexAIDdS receptor.pdb ligand.mol2
        tENCoM --help
    EOS
  end

  test do
    assert_path_exists bin/"FlexAIDdS"
    system bin/"FlexAIDdS", "--help"
    assert_path_exists bin/"tENCoM"
  end
end
