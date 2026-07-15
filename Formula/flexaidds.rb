class Flexaidds < Formula
  desc "Entropy-driven molecular docking engine (FlexAID + ΔS thermodynamic analysis)"
  homepage "https://github.com/LeBonhommePharma/FlexAIDdS"
  # v2.0.2 ships production MC_st0r5.2_6.dat (md5 9dc93717…) at repo root and
  # WRK/. v2.0.0 only had an outdated WRK matrix (md5 204b75ef…) that broke typing.
  url "https://github.com/LeBonhommePharma/FlexAIDdS/archive/refs/tags/v2.0.2.tar.gz"
  sha256 "11109a3eb6cac4185cee10390be7bf09be2520e89f13064a524e984be1366cc0"
  license "Apache-2.0"
  head "https://github.com/LeBonhommePharma/FlexAIDdS.git", branch: "master"

  livecheck do
    url :stable
    strategy :github_latest
  end

  # Optional Metal GPU path (macOS only). Off by default — see install notes.
  option "with-metal", "Build with Metal GPU acceleration (macOS; may fail on pure CLT SDKs)"

  depends_on "cmake" => :build
  depends_on "ninja" => :build
  depends_on "eigen"
  depends_on "libomp" if OS.mac?

  def install
    # Metal OFF by default: HEAD can fail to link metal_eval/metal_rmsd under
    # pure Command Line Tools SDKs. CPU + OpenMP is the portable brew path.
    metal = build.with?("metal") ? "ON" : "OFF"

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
      The FlexAIDdS native tools (FlexAIDdS, tENCoM, FlexAID) have been installed.

      Wrappers under #{bin} set FLEXAIDDS_DATA_DIR=#{libexec}/share so PATH
      symlinks still find MC matrices and AMINO.def.

      Stable v2.0.2+ includes the production docking matrix (atom typing works
      out of the box). Upgrade from broken v2.0.0 installs with:
        brew update && brew reinstall flexaidds

      Python package:
        pip install flexaidds
        # or from git: pip install "git+https://github.com/LeBonhommePharma/FlexAIDdS.git#subdirectory=python"

      Default brew build uses CPU + OpenMP (no Metal).
      Metal: brew install --with-metal flexaidds

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
