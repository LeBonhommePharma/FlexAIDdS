class Flexaidds < Formula
  desc "Entropy-driven molecular docking engine (FlexAID + ΔS thermodynamic analysis)"
  homepage "https://github.com/LeBonhommePharma/FlexAIDdS"
  # v2.0.2 ships production MC_st0r5.2_6.dat (md5 9dc93717…) at repo root and
  # WRK/. v2.0.0 only had an outdated WRK matrix (md5 204b75ef…) that broke typing.
  url "https://github.com/LeBonhommePharma/FlexAIDdS/archive/refs/tags/v2.0.2.tar.gz"
  sha256 "11109a3eb6cac4185cee10390be7bf09be2520e89f13064a524e984be1366cc0"
  license "Apache-2.0"
  # HEAD carries the flexaid_core Metal OBJCXX link fix (stable v2.0.2 does not).
  # Default branch is main (Git 3.0 / repo rename); brew install -s --HEAD uses this.
  head "https://github.com/LeBonhommePharma/FlexAIDdS.git", branch: "main"

  livecheck do
    url :stable
    strategy :github_latest
  end

  # Optional Metal GPU path (macOS only). Off by default — see install notes.
  # Requires a working Metal toolchain (Xcode + MetalToolchain). Stable v2.0.2
  # lacks the flexaid_core OBJCXX membership fix; use --HEAD --with-metal until
  # the next release tag that includes that CMake change (PR #260).
  option "with-metal", "Build with Metal GPU acceleration (macOS; needs Metal toolchain)"

  depends_on "cmake" => :build
  depends_on "ninja" => :build
  depends_on "eigen"
  depends_on "libomp" if OS.mac?

  def install
    # Metal OFF by default for CLT-only SDKs without a working metalc; use
    # --with-metal when Xcode/Metal toolchain is present. On sources that include
    # the flexaid_core Metal fix, OBJCXX bridges + frameworks are linked via
    # flexaid_core so every consumer (FlexAID, FlexAIDdS, cavity tools) links.
    metal = build.with?("metal") ? "ON" : "OFF"

    # Stable v2.0.2 tarball still attaches Metal .mm only to some executables;
    # targets that link flexaid_core alone (e.g. cavity_detect_cli) then fail with
    # undefined metal_eval_* / metal_rmsd::*. Prefer HEAD for Metal until the next
    # release ships the core membership fix.
    if build.with?("metal") && !build.head?
      odie <<~EOS
        --with-metal requires source that links Metal bridges via flexaid_core.
        Stable v2.0.2 still fails that link (undefined metal_eval_* from gaboom /
        hardware_detect / FOPTICS when building auxiliary targets).

        Until the next release tag, install Metal from HEAD. Homebrew 6 accepts
        --HEAD on install (not reinstall):
          brew uninstall flexaidds 2>/dev/null
          brew install -s --HEAD --with-metal lebonhommepharma/flexaidds/flexaidds

        Default (CPU + OpenMP, no Metal) still works on stable:
          brew reinstall lebonhommepharma/flexaidds/flexaidds
      EOS
    end

    args = std_cmake_args + %W[
      -GNinja
      -DCMAKE_BUILD_TYPE=Release
      -DBUILD_TESTING=OFF
      -DBUILD_PYTHON_BINDINGS=OFF
      -DFLEXAIDS_USE_CUDA=OFF
      -DFLEXAIDS_USE_METAL=#{metal}
      -DFLEXAIDS_USE_OPENMP=ON
      -DBUILD_FLEXAIDDS_FAST=ON
      -DENABLE_TENCOM_TOOL=ON
      -DENABLE_CAVITY_DETECT_CLI=OFF
      -DENABLE_BENCHMARK_DATASETS=OFF
      -DENABLE_DUAL_ASSEMBLY_TOOL=OFF
      -DENABLE_DIFT_TOOL=OFF
    ]

    if OS.mac?
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

    # Real binaries + data live in libexec so PATH symlinks under /opt/homebrew/bin
    # do not break relative data lookup. Wrappers set FLEXAIDDS_DATA_DIR.
    libexec.mkpath
    (libexec/"bin").mkpath
    (libexec/"share").mkpath

    %w[FlexAIDdS tENCoM FlexAID tencom_entropy_diff flexaids_process_ligand].each do |exe|
      %W[build/#{exe} build_lto/#{exe}].each do |candidate|
        next unless File.exist?(candidate)

        (libexec/"bin").install candidate
        break
      end
    end

    # Install data BOTH next to the real binary (base_path lookup) and under
    # share/ (FLEXAIDDS_DATA_DIR for engines that honor it after the path fix).
    install_runtime_data!(libexec/"bin")
    install_runtime_data!(libexec/"share")

    # Prefer root matrix over WRK/ when both exist (v2.0.2+ ships production both).
    write_wrappers!
  end

  def install_runtime_data!(dest)
    dest = Pathname(dest)
    dest.mkpath
    def_names = %w[AMINO.def NUCLEOTIDES.def]
    # Prefer repo-root / build-staged matrices over WRK/ (legacy outdated copy).
    data_candidates = %w[build build_lto . data]
    wrk_fallback = %w[WRK]

    copy_data = lambda do |dirs|
      dirs.each do |dir|
        Dir["#{dir}/MC_*.dat"].each do |f|
          target = dest/File.basename(f)
          cp f, target unless target.exist?
        end
        Dir["#{dir}/FA_matrix*.dat"].each do |f|
          target = dest/File.basename(f)
          cp f, target unless target.exist?
        end
        def_names.each do |name|
          path = Pathname("#{dir}/#{name}")
          next unless path.exist?

          target = dest/name
          cp path, target unless target.exist?
        end
      end
    end

    copy_data.call(data_candidates)
    # Only use WRK if we still lack the primary matrix.
    copy_data.call(wrk_fallback) if Dir["#{dest}/MC_*.dat"].empty?
  end

  def write_wrappers!
    data_dir = libexec/"share"
    %w[FlexAIDdS tENCoM FlexAID tencom_entropy_diff flexaids_process_ligand].each do |exe|
      real = libexec/"bin"/exe
      next unless real.exist?

      (bin/exe).write <<~EOS
        #!/bin/bash
        # Homebrew wrapper: force data dir next to installed runtime files and
        # exec the real binary (avoids argv[0] resolving to /opt/homebrew/bin).
        export FLEXAIDDS_DATA_DIR="${FLEXAIDDS_DATA_DIR:-#{data_dir}}"
        exec "#{real}" "$@"
      EOS
      chmod 0755, bin/exe
    end
  end

  def caveats
    <<~EOS
      This formula installs the *native* docking tools only:
        FlexAIDdS, FlexAID, tENCoM, tencom_entropy_diff

      It does *not* install the Python analysis package. Those are separate:
        # Python CLI + load_results / StatMech (GitHub until public PyPI):
        pip install "git+https://github.com/LeBonhommePharma/FlexAIDdS.git#subdirectory=python"
        # After the first PyPI release: pip install flexaidds
        # Then: flexaidds --help   or   python -m flexaidds --help

      Wrappers under #{bin} set FLEXAIDDS_DATA_DIR=#{libexec}/share so PATH
      symlinks still find MC matrices and AMINO.def.

      Stable v2.0.2+ includes the production docking matrix (atom typing works
      out of the box). Upgrade from broken v2.0.0 installs with:
        brew update && brew reinstall lebonhommepharma/flexaidds/flexaidds

      Install / reinstall path (Homebrew 6+ requires a real tap; raw URL installs
      are rejected):
        brew tap lebonhommepharma/flexaidds https://github.com/LeBonhommePharma/FlexAIDdS
        brew install lebonhommepharma/flexaidds/flexaidds

      Default brew build uses CPU + OpenMP (no Metal) and is the portable path.

      Metal GPU (macOS + Xcode Metal toolchain). Stable v2.0.2 cannot link Metal
      into every flexaid_core consumer; use HEAD until the next release.
      Homebrew 6: --HEAD is an install flag (reinstall --HEAD is invalid):
        brew uninstall flexaidds 2>/dev/null
        brew install -s --HEAD --with-metal lebonhommepharma/flexaidds/flexaidds
      After a post-v2.0.2 tag with the flexaid_core Metal fix:
        brew install --with-metal lebonhommepharma/flexaidds/flexaidds

      Example:
        FlexAIDdS receptor.pdb ligand.sdf --rigid -o /tmp/out
        tENCoM --help
    EOS
  end

  test do
    assert_path_exists bin/"FlexAIDdS"
    assert_path_exists libexec/"bin"/"FlexAIDdS"
    assert_path_exists libexec/"share"/"AMINO.def"
    assert Dir["#{libexec}/share/MC_*.dat"].any?, "expected MC_*.dat in libexec/share"

    # Wrapper must advertise the Cellar data dir, not /opt/homebrew/bin.
    output = shell_output("#{bin}/FlexAIDdS --help")
    assert_match "base path", output

    system bin/"tENCoM", "--help"
  end
end
