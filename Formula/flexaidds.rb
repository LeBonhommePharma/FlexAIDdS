class Flexaidds < Formula
  desc "Entropy-driven molecular docking engine (FlexAID + ΔS thermodynamic analysis)"
  homepage "https://github.com/LeBonhommePharma/FlexAIDdS"
  # v2.0.3 includes flexaid_core Metal OBJCXX membership (PR #260) so stable
  # --with-metal links. Still ships production MC_st0r5.2_6.dat (md5 9dc93717…).
  url "https://github.com/LeBonhommePharma/FlexAIDdS/archive/refs/tags/v2.0.3.tar.gz"
  sha256 "6c8442fc672a127db354ff3b6e08a2252e8c921372d902d062ecbf4296aef186"
  license "Apache-2.0"
  # HEAD always tracks main only. Never pin ephemeral fix/* branches — the
  # Homebrew tap clones this monorepo; a deleted branch breaks every reinstall
  # (fatal: couldn't find remote ref refs/heads/fix/…). Metal is a CMake option
  # (FLEXAIDS_USE_METAL), not a git branch.
  head "https://github.com/LeBonhommePharma/FlexAIDdS.git", branch: "main"

  livecheck do
    url :stable
    strategy :github_latest
  end

  # Optional Metal GPU path (macOS only). Off by default — see install notes.
  # Requires a working Metal toolchain (Xcode + MetalToolchain). Stable v2.0.3+
  # links Metal bridges via flexaid_core (PR #260).
  option "with-metal", "Build with Metal GPU acceleration (macOS; needs Metal toolchain)"

  depends_on "cmake" => :build
  depends_on "ninja" => :build
  depends_on "eigen"
  depends_on "libomp" if OS.mac?

  def install
    # Metal OFF by default for CLT-only SDKs without a working metalc; use
    # --with-metal when Xcode/Metal toolchain is present. v2.0.3+ attaches
    # OBJCXX bridges + frameworks PUBLIC on flexaid_core so every consumer
    # (FlexAID, FlexAIDdS, cavity tools) links cleanly.
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
      -DENABLE_CAVITY_DETECT_CLI=OFF
      -DENABLE_BENCHMARK_DATASETS=OFF
      -DENABLE_DUAL_ASSEMBLY_TOOL=OFF
      -DENABLE_DIFT_TOOL=OFF
    ]

    # Provenance. A `brew install` builds from an extracted tarball with no
    # .git, so CMake's git probe finds nothing. Two routes cover it:
    #   * stable — the release tarball carries .git_archival.txt, populated by
    #     `export-subst` at `git archive` time, and CMake reads the commit from
    #     it (source_provenance=archive).
    #   * --HEAD — Homebrew clones the repo, so the git probe works directly.
    # Tarballs cut BEFORE .git_archival.txt was added (v2.0.3 and earlier) have
    # neither; those builds stamp commit=unknown dirty=2, which is honest but
    # untraceable. Pass the tag as an explicit override so even an old tarball
    # is identifiable.
    unless build.head?
      args << "-DFLEXAIDS_GIT_COMMIT_OVERRIDE=v#{version}"
    end

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

    # One build directory, one producer.
    #
    # This selected the binary here and the provenance JSON in a separate loop
    # below, each re-deciding build/ vs build_lto/ from scratch. Two independent
    # decisions over one fact can disagree: CMake writes the provenance file to
    # CMAKE_BINARY_DIR, so it describes one configure of one directory. A stale
    # build/ holding a binary but no JSON, next to a build_lto/ holding a JSON,
    # installed a provenance record for a build that was not the installed
    # build -- and said nothing. The per-exe inner loop made it worse: FlexAIDdS
    # could come from build/ and tENCoM from build_lto/, so a single JSON could
    # not have described both even in principle.
    #
    # That is the defect this branch exists to remove, in the branch's own code:
    # a claim asserted about something that was never checked. Resolve the
    # directory once; take everything from it.
    build_dir = %w[build build_lto].find { |dir| File.exist?("#{dir}/FlexAIDdS") }
    odie "no FlexAIDdS binary in build/ or build_lto/ after a successful build" if build_dir.nil?
    ohai "installing from #{build_dir}/"

    %w[FlexAIDdS tENCoM FlexAID tencom_entropy_diff flexaids_process_ligand].each do |exe|
      candidate = "#{build_dir}/#{exe}"
      next unless File.exist?(candidate)

      (libexec/"bin").install candidate
    end

    # Install data BOTH next to the real binary (base_path lookup) and under
    # share/ (FLEXAIDDS_DATA_DIR for engines that honor it after the path fix).
    install_runtime_data!(libexec/"bin")
    install_runtime_data!(libexec/"share")

    # Make the installed engine identifiable without running a dock. The commit
    # otherwise only surfaces in the REMARK block of emitted pose files, so a
    # freshly installed binary cannot be told apart from any other build.
    # Bound to build_dir, not re-decided. If that directory produced a binary
    # but no provenance file, the install is unidentifiable; fail rather than
    # borrow a plausible-looking JSON from a directory that produced nothing
    # that is being installed.
    prov_src = "#{build_dir}/flexaidds-build-provenance.json"
    unless File.exist?(prov_src)
      odie "#{build_dir}/ produced FlexAIDdS but no flexaidds-build-provenance.json; " \
           "installing would put an unidentifiable binary on PATH"
    end
    (libexec/"share").install prov_src

    (bin/"flexaidds-buildinfo").write <<~EOS
      #!/bin/bash
      # Report what the installed FlexAIDdS actually is.
      set -euo pipefail
      prov="#{libexec}/share/flexaidds-build-provenance.json"
      if [ -f "$prov" ]; then
        cat "$prov"
      else
        echo '{"error":"no provenance file installed"}'
        exit 1
      fi
      # The previous "second, independent source" grepped the binary for
      # ^REMARK FLEXAID.commit= . That pattern matches the printf FORMAT string
      # in LIB/BindingMode.cpp:790, so on a real binary it returns
      #   REMARK FLEXAID.commit=%s FLEXAID.dirty=%d FLEXAID.seed=%llu
      # -- a literal %s, never a commit. It was not a second source, and
      # nothing compared it to the first, so a disagreement could not surface.
      #
      # The commit IS in the binary, as the separate string constant that %s
      # consumes, but it is an unanchored token: its presence can be tested,
      # its identity cannot be read back. So this checks the one direction that
      # is sound. Absence is proof of disagreement; presence is consistency,
      # not verification, and is labelled as such.
      json_commit=$(grep -o '"git_commit"[^,]*' "$prov" | head -1 | cut -d'"' -f4 || true)
      echo "--- identity cross-check ---"
      echo "provenance_json_commit=${json_commit:-missing}"
      if [ -z "${json_commit:-}" ] || [ "$json_commit" = "unknown" ]; then
        echo "status=unverifiable (provenance file carries no commit to look for)"
        exit 0
      fi
      if /usr/bin/strings "#{libexec}/bin/FlexAIDdS" 2>/dev/null | grep -qxF "$json_commit"; then
        echo "status=consistent (commit from the provenance file is present in the binary)"
        exit 0
      fi
      echo "status=MISMATCH"
      echo "The installed provenance file names a commit that does not appear in" >&2
      echo "the installed binary. They describe different builds; do not trust" >&2
      echo "results attributed to this install." >&2
      exit 1
    EOS
    chmod 0755, bin/"flexaidds-buildinfo"

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

      Stable v2.0.3+ includes the production docking matrix (atom typing works
      out of the box). Upgrade from broken v2.0.0 / stale HEAD installs with:
        brew update && brew reinstall lebonhommepharma/flexaidds/flexaidds

      Install / reinstall path (Homebrew 6+ requires a real tap; raw URL installs
      are rejected). The tap is this monorepo — keep the tap checkout on main:
        brew tap lebonhommepharma/flexaidds https://github.com/LeBonhommePharma/FlexAIDdS
        # Prefer formula-scoped trust when HOMEBREW_REQUIRE_TAP_TRUST is set
        # (https://docs.brew.sh/Tap-Trust). Do not use HOMEBREW_NO_REQUIRE_TAP_TRUST.
        brew trust --formula lebonhommepharma/flexaidds/flexaidds
        brew install lebonhommepharma/flexaidds/flexaidds

      Default brew build uses CPU + OpenMP (no Metal) and is the portable path.
      Formula head tracks main only (never ephemeral fix/* branches).

      If reinstall fails with "couldn't find remote ref" after an old --HEAD
      install that pinned a deleted branch (e.g. fix/homebrew-metal-link):
        cd "$(brew --repository lebonhommepharma/flexaidds)"
        git fetch --prune && git checkout main && git reset --hard origin/main
        brew uninstall --force lebonhommepharma/flexaidds/flexaidds
        brew install --build-from-source lebonhommepharma/flexaidds/flexaidds

      Metal GPU (macOS + Xcode Metal toolchain). Stable v2.0.3+ links Metal
      bridges via flexaid_core (PR #260); no special git branch required:
        brew install --build-from-source --with-metal lebonhommepharma/flexaidds/flexaidds
      Or after tap update if already installed:
        brew reinstall --build-from-source --with-metal lebonhommepharma/flexaidds/flexaidds

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

    # The installed engine must be identifiable. A build that cannot be tied to
    # a source revision silently invalidates any comparison it takes part in.
    assert_path_exists libexec/"share"/"flexaidds-build-provenance.json"
    prov = JSON.parse((libexec/"share"/"flexaidds-build-provenance.json").read)
    refute_empty prov["git_commit"].to_s, "provenance file has an empty commit"
    assert_match(/^(git|archive|override)$/, prov["source_provenance"],
                 "installed build has unrecoverable provenance: #{prov}")

    # This asserted that buildinfo's output contained the JSON's commit -- but
    # buildinfo begins by printing that same JSON, so the assertion compared a
    # string to itself and could not fail. Assert against the cross-check line,
    # which is derived from the binary.
    info = shell_output("#{bin}/flexaidds-buildinfo")
    assert_match "provenance_json_commit=#{prov['git_commit']}", info
    refute_match "status=MISMATCH", info
  end
end
