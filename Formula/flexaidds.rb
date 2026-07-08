class Flexaidds < Formula
  desc "Entropy-driven molecular docking engine (FlexAID + ΔS thermodynamic analysis)"
  homepage "https://github.com/LeBonhommePharma/FlexAIDdS"
  license "Apache-2.0"
  head "https://github.com/LeBonhommePharma/FlexAIDdS.git", branch: "master"
  # Stable releases will be added here with proper tarball + sha256 when cutting releases.
  # For now the formula is optimized for `brew install --HEAD ...` and direct formula URLs.

  depends_on "cmake" => :build
  depends_on "ninja" => :build
  depends_on "eigen"
  depends_on "libomp" if OS.mac?
  depends_on "python@3.11" => :recommended # for the flexaidds Python package (recommended install separately via pip too)

  def install
    # Build the main tools. Focus on the fast optimized FlexAIDdS binary,
    # the tENCoM tool, and legacy FlexAID for compatibility.
    args = std_cmake_args + %W[
      -G Ninja
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
      libomp = Formula["libomp"].opt_prefix
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
    bin.install "build/FlexAIDdS"
    bin.install "build/tENCoM"
    bin.install "build/FlexAID" if File.exist?("build/FlexAID")

    # Install co-located data files (staged by the build)
    # These are looked up relative to the executable directory.
    bin.install Dir["build/MC_*.dat"]
    bin.install "build/AMINO.def" if File.exist?("build/AMINO.def")
    bin.install "build/NUCLEOTIDES.def" if File.exist?("build/NUCLEOTIDES.def")
    bin.install Dir["build/FA_matrix*.dat"]

    # Optional useful tools
    bin.install "build/flexaids_process_ligand" if File.exist?("build/flexaids_process_ligand")
    bin.install "build/tencom_entropy_diff" if File.exist?("build/tencom_entropy_diff")

    # Note: The Python package "flexaidds" is best installed via pip:
    #   pip install flexaidds
    # This formula focuses on the high-performance native CLI tools.
  end

  def caveats
    <<~EOS
      The FlexAIDdS native tools (FlexAIDdS, tENCoM, FlexAID) have been installed.

      Runtime data files (MC matrices, AMINO.def, etc.) are installed alongside the
      binaries so they are found automatically.

      For the Python package (analysis, results loader, thermodynamics):
        pip install flexaidds

      To use the full Python + native integration, install the Python package
      after this formula.

      Example usage:
        FlexAIDdS receptor.pdb ligand.mol2
        tENCoM --help
    EOS
  end

  test do
    # Basic smoke test - the binaries should exist and report something
    assert_predicate bin/"FlexAIDdS", :exist?
    # Most of these tools accept -h or --help; fall back to existence + basic run
    system "#{bin}/FlexAIDdS", "--help" if (bin/"FlexAIDdS").exist?
    assert_predicate bin/"tENCoM", :exist?
  end
end