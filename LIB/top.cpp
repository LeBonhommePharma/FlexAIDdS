#include "TuiColor.h"
#include "version_info.h"
#include "gaboom.h"
#include "top_helpers.h"
#include "fileio.h"
#include "flexaid_exception.h"
#include "Vcontacts.h"
#include "config_parser.h"
#include "config_defaults.h"
#include "Mol2Reader.h"
#include "SdfReader.h"
#include "CifReader.h"
#include "CleftDetector.h"
#include "site_confine.h"
#include "statmech.h"
#include "TargetServer.h"  // P1: for grand canonical TargetServer context in cluster paths
#include "ProcessLigand/ProcessLigand.h"
#include "ProcessLigand/CoordBuilder.h"
#include "LibrarySplitter.h"
#include "ReferenceEntropy.h"
#include "assign_formal_charges.h"
#include "atom_typing_256.h"
#include "RefLigSeed.h"
#include "CoarseScreen.h"
#include "TwoStageScreen.h"
#include "GISTEvaluator.h"
#include "ParallelDock.h"
#include "ParallelCampaign.h"
#include "GAContext.h"
// FLEXAIDDS_CMAES_INCLUDE_BEGIN
#include "cmaes_search.h"
// FLEXAIDDS_CMAES_INCLUDE_END
#include "MIFGrid.h"
#include "CavityDetect/SpatialGrid.h"
#include "native_score.h"
#include "rescore_pool.h"
#include "hbond_potential.h"
#include "RngSeed.h"
#include "ensemble_pipeline.h"
#include "ProtocolConfig.h"
#include "shell_exec.h"
#include "UnifiedHardwareDispatch.h"
#include "flexaidds_flags.h"
#if defined(FLEXAIDDS_ENABLE_REDOCK)
#include "DatasetRunner.h"
#include "DatasetRunnerProvenance.h"
#include "RunReceipt.h"
#if defined(__APPLE__)
#include <mach-o/dyld.h>
#endif
#endif

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <string>
#include <vector>
#include <filesystem>
#include <system_error>
#include <unistd.h>
#ifndef _WIN32
#include <sys/wait.h>
#endif

// ── Tier 2 typing: SYBYL name → canonical VCT matrix index ──────────────────
// ProcessLigand (BonMol) perceives full hybridisation/aromaticity but stores
// the result in its own internal SYBYL numbering (see SybylTyper.h). For the
// scorer, atom.type must be a *canonical* VCT row index matching
// nrgrank_matrix.h / assign_radii_types.cpp and Mol2Reader's canonical mapping:
//   1=C.1  2=C.2  3=C.3  4=C.AR  5=C.CAT
//   6=N.1  7=N.2  8=N.3  9=N.4  10=N.AR 11=N.AM 12=N.PL3
//  13=O.2 14=O.3 15=O.CO2 16=O.AR
//  17=S.2 18=S.3 19=S.O  20=S.O2 21=S.AR
//  22=P.3 23=F   24=CL   25=BR   26=I   27=SE
//  28=MG 29=SR 30=CU 31=MN 32=HG 33=CD 34=NI 35=ZN 36=CA 37=FE 38=CO.OH
//  39=DUMMY 40=SOLVENT
// We bridge via the canonical SYBYL string name produced by
// bonmol::sybyl::sybyl_type_name(), keeping a single source of truth for the
// string→canonical mapping (identical to Mol2Reader::sybyl_to_flexaid_type).
//
// ── Row liveness in MC_st0r5.2_6.dat, and every substitution ────────────────
// A row with no non-zero entries makes the atom *invisible* to the contact
// scorer: it still occupies volume (wall/clash) but contributes nothing to CF.
// Any perceived type landing on a dead row must therefore be aliased onto the
// nearest live row that reproduces its chemistry. Dead/near-dead rows:
//
//   row  8 N.3    0 entries → alias N.3  → 11 (N.am)  sp3 amine
//   row 16 O.ar   0 entries → alias O.ar → 14 (O.3)   furan / oxazole O
//   row 20 S.O2   0 entries → alias S.O2 → 19 (S.O)   sulfone / sulfonamide S
//   row 21 S.ar   0 entries → alias S.ar → 18 (S.3)   thiophene / thiazole S
//   row 27 Se     0 entries → alias Se   → 18 (S.3)   selenomethionine
//   row 26 I      3 entries → alias I    → 25 (Br)    heavy halogen
//   row  1 C.1   10 entries → alias C.1  →  2 (C.2)   judgement call, not dead:
//                 sp C is rare in the training set and row 1 carries extreme,
//                 poorly-sampled values (1-13 = -198.3, the matrix maximum).
//   rows 32 Hg, 38 Co.oh, 39 DUMMY: 0 entries; 29 Sr, 33 Cd: 1 entry. Left as
//                 is — these are steric-only by construction (DUMMY) or too
//                 rare to be worth a surrogate.
//
// Rows deliberately NOT substituted: 7 (N.2, 13 entries) and 6 (N.1, 5
// entries) are live; coercing them to N.ar/N.am discards real chemistry.
//
// An alias changes only the *scoring* row. Where it would also discard H-bond
// chemistry, the pre-substitution row is recorded in atom_struct::sybyl_orig
// and the geometry consumers dispatch on that instead — currently N.3, whose
// sp3 virtual-H geometry differs from N.am's planar amide bisector. The O.ar /
// S.ar / Se aliases need no such record: their geometry recipes and implicit-H
// counts are gated on heavy_bonds<=1, and a ring or chain heteroatom taking
// those aliases always has 2 heavy neighbours, so donor status is unchanged.
//
// Any change here MUST be mirrored in Mol2Reader::sybyl_to_flexaid_type,
// SdfReader::element_to_flexaid_type and read_coor.cpp:
// canonical_vct_type_for_element, or the same element will score differently
// depending on which file format it arrived in and which side of the complex
// it sits on.

// ════════════════════════════════════════════════════════════════════════════
//  Strategy A — on-disk VCT grid cache  (env: FLEXAIDDS_GRID_CACHE_DIR)
// ════════════════════════════════════════════════════════════════════════════
// For non-native cross-docking the same receptor is paired with dozens of
// ligands; the receptor cleft grid (detect_cleft → generate_grid → calc_cleftic
// → SITE-CONFINE) is purely receptor-determined yet currently rebuilt from
// scratch for every pair.  When FLEXAIDDS_GRID_CACHE_DIR is set we key a binary
// snapshot of the finalized (site-confined) grid on the receptor + site PDB
// content and the grid spacing/permeability, so the first ligand of a receptor
// builds-and-writes the grid and every subsequent ligand loads it from disk.
//
// The key is an inline FNV-1a 64-bit hash — no external crypto dependency, and
// collision risk is negligible for the ~85 receptors in a benchmark run.  The
// feature is gated entirely on the env var: when it is absent the helpers below
// are never invoked and behaviour is bit-for-bit identical to before.
namespace gridcache {

static constexpr uint32_t kMagic   = 0x56435400u; // "VCT\0"
static constexpr uint32_t kVersion = 1u;

// FNV-1a (64-bit) — fold a byte buffer into a running hash.
static inline uint64_t fnv1a_update(uint64_t h, const void* data, size_t n) {
	const unsigned char* p = static_cast<const unsigned char*>(data);
	for (size_t i = 0; i < n; ++i) {
		h ^= static_cast<uint64_t>(p[i]);
		h *= 0x100000001b3ull;
	}
	return h;
}

// Fold an entire file's content into the running hash.  Returns false when the
// file cannot be opened.
static bool fnv1a_file(uint64_t& h, const char* path) {
	FILE* fp = fopen(path, "rb");
	if (!fp) return false;
	unsigned char buf[65536];
	size_t got;
	while ((got = fread(buf, 1, sizeof(buf), fp)) > 0)
		h = fnv1a_update(h, buf, got);
	fclose(fp);
	return true;
}

// Resolve the cache file path for this (receptor, site, spacer, permea) tuple,
// or an empty string when caching is disabled or the receptor is unreadable.
static std::string cache_path(const char* receptor_file,
                              const char* site_file,
                              float spacer, float permea) {
	const char* dir = std::getenv("FLEXAIDDS_GRID_CACHE_DIR");
	if (!dir || dir[0] == '\0') return std::string();
	if (!receptor_file || receptor_file[0] == '\0') return std::string();

	uint64_t h = 0xcbf29ce484222325ull; // FNV-1a offset basis
	if (!fnv1a_file(h, receptor_file)) return std::string(); // no receptor → no key
	// Site PDB is optional (AUTO mode has none); fold it only when present so
	// AUTO and oracle runs of the same receptor never share a key.
	if (site_file && site_file[0] != '\0')
		fnv1a_file(h, site_file); // absent/unreadable site contributes nothing
	h = fnv1a_update(h, &spacer, sizeof(spacer));
	h = fnv1a_update(h, &permea, sizeof(permea));

	char name[32];
	snprintf(name, sizeof(name), "%016llx.vct", static_cast<unsigned long long>(h));
	return std::string(dir) + "/" + name;
}

// Load num_grd + gridpoint[] from a cache file.  On a validated hit *out_grid is
// malloc'd (caller owns / frees with free()) and *out_num is set.  gridpoint is
// a fixed-layout POD (ints + floats), so the raw struct array is portable across
// runs of the same binary.
static bool load(const std::string& path, gridpoint** out_grid, int* out_num) {
	FILE* fp = fopen(path.c_str(), "rb");
	if (!fp) return false;
	uint32_t magic = 0, version = 0;
	int32_t  num = 0;
	bool ok = (fread(&magic,   sizeof(magic),   1, fp) == 1) &&
	          (fread(&version, sizeof(version), 1, fp) == 1) &&
	          (fread(&num,     sizeof(num),     1, fp) == 1) &&
	          magic == kMagic && version == kVersion &&
	          num > 0 && num < (1 << 24);
	gridpoint* grid = nullptr;
	if (ok) {
		grid = static_cast<gridpoint*>(malloc(static_cast<size_t>(num) * sizeof(gridpoint)));
		ok = grid && (fread(grid, sizeof(gridpoint), static_cast<size_t>(num), fp)
		              == static_cast<size_t>(num));
	}
	fclose(fp);
	if (!ok) { if (grid) free(grid); return false; }
	*out_grid = grid;
	*out_num  = num;
	return true;
}

// Atomically persist num_grd + gridpoint[] (temp file + rename) so concurrent
// workers building the same receptor never observe a half-written cache file.
static void save(const std::string& path, const gridpoint* grid, int num) {
	if (num <= 0 || grid == nullptr) return;
	std::string tmp = path + ".tmp";
	FILE* fp = fopen(tmp.c_str(), "wb");
	if (!fp) return;
	uint32_t magic = kMagic, version = kVersion;
	int32_t  n = num;
	bool ok = (fwrite(&magic,   sizeof(magic),   1, fp) == 1) &&
	          (fwrite(&version, sizeof(version), 1, fp) == 1) &&
	          (fwrite(&n,       sizeof(n),       1, fp) == 1) &&
	          (fwrite(grid, sizeof(gridpoint), static_cast<size_t>(num), fp)
	           == static_cast<size_t>(num));
	fclose(fp);
	std::error_code ec;
	if (ok) {
		std::filesystem::rename(tmp, path, ec);
		if (ec) std::filesystem::remove(tmp, ec);
	} else {
		std::filesystem::remove(tmp, ec);
	}
}

} // namespace gridcache

// ════════════════════════════════════════════════════════════════════════════
//  Strategy B — multi-ligand batch dispatch  (env: FLEXAIDDS_LIGAND_BATCH)
// ════════════════════════════════════════════════════════════════════════════
// Re-exec this same executable once per ligand in a fresh process.  A clean
// process per ligand deliberately avoids resetting the deeply-shared atom /
// residue / FA / VC / chromosome state that main() allocates once and frees only
// at exit — an in-process loop would have to unwind all of it correctly between
// ligands.  The receptor VCT grid is still built only once per receptor because
// each child consults the Strategy A on-disk cache (FLEXAIDDS_GRID_CACHE_DIR):
// the first ligand writes the grid, the rest load it.  FLEXAIDDS_LIGAND_BATCH is
// cleared in the child to prevent infinite recursion.
#ifndef _WIN32
static int batch_exec_child(const char* exe, const std::vector<std::string>& args) {
	pid_t pid = fork();
	if (pid < 0) return -1;
	if (pid == 0) {
		unsetenv("FLEXAIDDS_LIGAND_BATCH"); // child docks a single ligand
		std::vector<char*> cargv;
		cargv.reserve(args.size() + 1);
		for (const auto& a : args) cargv.push_back(const_cast<char*>(a.c_str()));
		cargv.push_back(nullptr);
		execvp(exe, cargv.data()); // PATH-search if exe has no '/', else literal path
		_exit(127);                // exec failed
	}
	int status = 0;
	while (waitpid(pid, &status, 0) < 0 && errno == EINTR) { /* retry */ }
	if (WIFEXITED(status)) return WEXITSTATUS(status);
	return -1;
}
#endif

// Resolve a fully-perceived BonMol atom to its canonical VCT index.
// BonMol's SybylTyper follows the classic SYBYL convention of typing aromatic
// ring heteroatoms structurally: thiophene/thiadiazole S becomes S.3 and
// furan/oxazole O becomes O.2/O.3. The canonical VCT energy matrix, however,
// carries dedicated aromatic-heteroatom rows (C.AR=4, N.AR=10, O.AR=16,
// S.AR=21). When the atom sits in a perceived aromatic ring, route it to the
// matching .ar row so the scorer sees aromatic chemistry rather than the
// sp2/sp3 row. (Aromatic C and N already resolve to C.ar/N.ar via the SYBYL
// name, but we override uniformly for clarity and to cover both heteroatoms.)
// out_name (optional) receives the human-readable SYBYL name used for logging.
static int bonmol_atom_to_canonical_vct(const bonmol::Atom& a, const char** out_name) {
	if (a.is_aromatic) {
		switch (a.element) {
			case bonmol::Element::C: if (out_name) *out_name = "C.ar"; return 4;
			case bonmol::Element::N: if (out_name) *out_name = "N.ar"; return 10;
			// O.ar/S.ar rows (16/21) are all-zero in the matrix; route to the
			// live ether/thioether surrogates. See sybyl_name_to_canonical_vct.
			// Without this, BonMol's aromatic perception actively *downgraded*
			// SdfReader's already-correct 14/18 to a zero-scoring row.
			case bonmol::Element::O: if (out_name) *out_name = "O.ar"; return 14;
			case bonmol::Element::S: if (out_name) *out_name = "S.ar"; return 18;
			default: break;
		}
	}
	const char* sname = bonmol::sybyl::sybyl_type_name(a.sybyl_type);
	if (out_name) *out_name = sname;
	return sybyl_name_to_canonical_vct(sname);
}

// ── Idiotproof file role detection ──────────────────────────────────────────

/// Download RCSB entry, extract cognate ligand, strip apo receptor.
/// Returns true and fills paths on success. Cache under FLEXAIDDS_REDOCK_CACHE
/// or ~/.flexaidds/benchmarks/redock/<PDBID>/.
/// Built only into the FlexAIDdS target (FLEXAIDDS_ENABLE_REDOCK).
static bool prepare_redock_from_rcsb(const std::string& pdb_id,
                                     std::string& receptor_path,
                                     std::string& ligand_path,
                                     std::string& out_pdb_id) {
#if !defined(FLEXAIDDS_ENABLE_REDOCK)
	(void)pdb_id;
	(void)receptor_path;
	(void)ligand_path;
	(void)out_pdb_id;
	fprintf(stderr,
	        "ERROR: --redock is only available in the FlexAIDdS binary "
	        "(this build is FlexAID without DatasetRunner cognate prep).\n");
	return false;
#else
	std::string upper = pdb_id;
	std::transform(upper.begin(), upper.end(), upper.begin(),
	               [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
	if (!is_valid_pdb_id(upper)) {
		fprintf(stderr,
		        "ERROR: --redock requires a valid RCSB PDB ID (4–8 alphanumeric), got '%s'\n",
		        pdb_id.c_str());
		return false;
	}

	std::string cache_dir;
	if (const char* env = std::getenv("FLEXAIDDS_REDOCK_CACHE")) {
		if (env[0] != '\0') cache_dir = env;
	}
	dataset::DatasetRunner runner(cache_dir);
	printf("[REDOCK] Preparing cognate redock for %s (RCSB download + ligand extract + apo strip)\n",
	       upper.c_str());
	dataset::DatasetEntry entry = runner.prepare_pdb_entry(upper, "redock");
	if (entry.receptor_path.empty() || !std::filesystem::exists(entry.receptor_path)) {
		fprintf(stderr, "ERROR: --redock %s failed to prepare apo receptor from RCSB\n",
		        upper.c_str());
		return false;
	}
	if (entry.ligand_path.empty() || !std::filesystem::exists(entry.ligand_path)) {
		fprintf(stderr,
		        "ERROR: --redock %s failed to extract a cognate ligand "
		        "(no non-water/ion HETATM after cofactor filters)\n",
		        upper.c_str());
		return false;
	}
	receptor_path = entry.receptor_path;
	ligand_path = entry.ligand_path;
	out_pdb_id = entry.pdb_id.empty() ? upper : entry.pdb_id;
	printf("[REDOCK] PDB %s\n", out_pdb_id.c_str());
	printf("[REDOCK] Apo receptor: %s\n", receptor_path.c_str());
	printf("[REDOCK] Cognate ligand: %s\n", ligand_path.c_str());
	printf("[REDOCK] docking_mode: self_docking (cognate ligand → native holo-derived apo)\n");
	return true;
#endif
}

#if defined(FLEXAIDDS_ENABLE_REDOCK)
/// Write RUN_RECEIPT.json for a `--redock` run whose `-o` names an existing
/// directory. Before this, `--redock` wrote no receipt at all, so a redock run
/// carried no engine identity.
///
/// `-o` is a file *prefix* everywhere else in this file: `-o run/1STP` yields
/// `run/1STP.rrg`, `run/1STP_INI.pdb` and so on, and that is unchanged. This is
/// a strictly additive second behaviour, fired only when the prefix happens to
/// name a directory that already exists, in which case the receipt is written
/// inside it. No docking artifact moves.
///
/// Conforms to the engine dialect in docs/run-uniformity/RUN_RECEIPT_CONTRACT.md.
/// This is the *third* writer into a format that already has two dialects, so it
/// calls flexaids::write_run_receipt rather than emitting JSON by hand (§7.1).
/// That keeps two contract details automatic instead of re-derived here: the
/// two-booleans-two-encodings rule (§2.1 — seed_elitism as 1/0 while
/// oracle_site_dir_set is true/false, inconsistent and load-bearing), and the
/// fixed 6-decimal float formatting (§2.2).
///
/// Best-effort by design: a read-only or full output directory must still yield
/// a successful docking run, so every failure path warns and returns.
static void write_redock_run_receipt(const std::string& output_dir,
                                     const std::string& pdb_id,
                                     double temperature_K,
                                     int pop,
                                     int gen) {
	const flexaids::ProtocolConfig proto = flexaids::ProtocolConfig::from_env();

	// In --redock there is no separate engine process: this binary both drives
	// the run and does the docking, so binary_* and runner_* are the same file.
	// DatasetRunner distinguishes them because it spawns a child; recording
	// them as equal here is the accurate statement, not a shortcut.
	std::string self_path;
#if defined(__APPLE__)
	{
		char buf[4096];
		std::uint32_t size = static_cast<std::uint32_t>(sizeof(buf));
		if (_NSGetExecutablePath(buf, &size) == 0) self_path = buf;
	}
#elif defined(__linux__)
	{
		std::error_code lec;
		const auto self = std::filesystem::read_symlink("/proc/self/exe", lec);
		if (!lec) self_path = self.string();
	}
#endif

	const std::string matrix_path =
	    dataset::resolve_scoring_matrix_path(proto.data_dir, self_path);

	flexaids::RunReceiptInput receipt;
	receipt.run_id       = pdb_id.empty() ? std::string("redock") : pdb_id;
	receipt.started_utc  = flexaids::utc_now_iso8601();
	receipt.output       = output_dir;
	receipt.dataset      = receipt.run_id;
	receipt.mode         = "defined-cleft-redock";
	receipt.temperature_K = temperature_K;
	receipt.pop          = pop;
	receipt.gen          = gen;
	receipt.restarts     = std::max(1, proto.restarts);
	receipt.seed_base    = proto.seed_base;

	// seed_elitism is DERIVED, never copied — §3 / DatasetRunner.cpp:5341-5347.
	// Only ORACLE_CEILING forces it true; DEFINED_CLEFT_REDOCK, AUTONOMOUS and
	// UNSET all force false. --redock is a defined-cleft redock and is never the
	// oracle-ceiling benchmark, so the derived value is false. Copying
	// proto.seed_elitism straight through would publish the same key with a
	// different meaning: silent, and visible later only as an inexplicable
	// cross-arm difference.
	receipt.seed_elitism = false;

	receipt.matrix_path   = matrix_path;
	receipt.matrix_md5    = dataset::provenance_file_md5(matrix_path);
	receipt.matrix_sha256 = dataset::provenance_file_sha256(matrix_path);
	receipt.binary_path   = self_path;
	receipt.binary_sha256 = dataset::provenance_file_sha256(self_path);
	receipt.runner_path   = self_path;
	receipt.runner_sha256 = receipt.binary_sha256;

	// Build-stamped commit rather than shelling out to git — §4 and §7.5.
	// Drivers cd into the arm directory before launching, so `git rev-parse HEAD`
	// runs outside any checkout and comes back empty; it is empty in the sampled
	// production receipt. CMakeLists.txt:54-61 stamps this at configure time and
	// run_t13_twotarget.sh recovers the same value with `strings`, refusing to
	// launch an unstamped binary. Note it is `rev-parse --short`, so this field
	// is a short SHA where DatasetRunner's — when it succeeds — is full length.
#if defined(FLEXAIDS_GIT_COMMIT)
	receipt.git_commit = FLEXAIDS_GIT_COMMIT;
#endif

	receipt.oracle_site_dir     = proto.oracle_site_dir;
	receipt.oracle_site_dir_set = !proto.oracle_site_dir.empty();
	receipt.protocol            = proto;

	// provenance.json is deliberately NOT written; §2.6 requires the second
	// writer to decide rather than default. --redock has never written a receipt
	// of any kind, so no existing tool can be looking for a provenance.json from
	// this path and nothing regresses by its absence. DatasetRunner still passes
	// true and is untouched.
	bool ok = false;
	try {
		ok = flexaids::write_run_receipt(output_dir, receipt,
		                                 /*also_write_provenance_json=*/false);
	} catch (...) {
		ok = false;
	}

	if (ok) {
		fprintf(stderr, "[RECEIPT] wrote %s/RUN_RECEIPT.json\n", output_dir.c_str());
	} else {
		// Warn, never abort. A read-only or full output directory currently
		// yields a successful docking run and must continue to.
		fprintf(stderr,
		        "[WARN] could not write RUN_RECEIPT.json to %s — docking continues\n",
		        output_dir.c_str());
	}
}
#endif  // FLEXAIDDS_ENABLE_REDOCK

static void print_usage(const char* progname) {
	tui::brand(); printf(" %s—%s Entropy-driven molecular docking\n\n", tui::muted(), tui::reset());
	printf("Usage:\n");
	printf("  %s <receptor> <ligand> [options]\n\n", progname);
	printf("  Files can be in any order. FlexAIDdS auto-detects which is\n");
	printf("  the receptor and which is the ligand from file content.\n\n");
	printf("  Receptor: .pdb, .cif, .mmcif (protein/nucleic acid)\n");
	printf("  Ligand:   .mol2, .sdf, .mol, .pdb (small molecule)\n");
	printf("            or a SMILES string directly on the command line\n\n");
	printf("  %s --redock <PDBid> [options]\n", progname);
	printf("      Download RCSB entry <PDBid>, extract the cognate ligand,\n");
	printf("      strip it from the target (apo), and redock automatically.\n\n");
	printf("  %s --legacy <config.inp> <ga.inp> <output_prefix>\n\n", progname);
	printf("Options:\n");
	printf("  -c, --config <file.json>   JSON config (overrides defaults)\n");
	printf("  -o, --output <prefix>      Output prefix (default: flexaid_out;\n");
	printf("                             default <PDBid>_redock with --redock)\n");
	printf("  --backend <cpu|metal|webgpu>  GPU compute backend for CF scoring (default: cpu)\n");
	printf("  --rigid                    Fast rigid-body screening\n");
	printf("  --screen                   Coarse-grained cube screening (Stage 1)\n");
	printf("  --screen-top-n <N>         Return top N from coarse screen (default: 100)\n");
	printf("  --screen-dock              Stage 2 hook on top-N (default OFF; surrogate\n");
	printf("                             unless a real GA callback is registered)\n");
	printf("  --screen-target-mol2 <f>   MOL2 target for cube screen / prefilter\n");
	printf("  --parallel-dock            Grid-decomposed parallel docking (ParallelDock)\n");
	printf("  --parallel-dock-regions <N> Number of spatial regions (default: 128)\n");
	printf("  --campaign                 Parallel virtual screening campaign mode\n");
	printf("  --coarse-prefilter         Campaign: cube-screen library, dock top-N only\n");
	printf("  --coarse-prefilter-top-n N Top-N for --coarse-prefilter (default: 100)\n");
	printf("  --folded                   Skip NATURaL chain growth\n");
	printf("  --legacy                   Legacy 3-file input mode\n");
	printf("  --redock <PDBid>           Cognate redock from RCSB PDB ID\n");
	printf("  --benchmark <set>          Run benchmark dataset (astex, casf2016, etc.)\n");
	printf("  -h, --help                 Show this help\n");
	printf("  --version                  Build identity as key=value lines\n\n");
	printf("Library input (virtual screening):\n");
	printf("  Ligand can be a multi-molecule SDF, a SMILES file (.smi),\n");
	printf("  or a directory of MOL2/SDF files. Each ligand is docked\n");
	printf("  independently against the receptor.\n\n");
	printf("Multi-model receptor (NMR / cryo-EM / MD):\n");
	printf("  PDB with MODEL/ENDMDL records or multi-model CIF.\n");
	printf("  Each model is used as a separate receptor conformer.\n");
	printf("  Results are combined via Boltzmann ensemble consensus\n");
	printf("  with reference entropy correction.\n\n");
	printf("Examples:\n");
	printf("  %s receptor.pdb ligand.mol2\n", progname);
	printf("  %s ligand.sdf receptor.pdb          # order doesn't matter\n", progname);
	printf("  %s receptor.pdb 'c1ccccc1' --rigid  # SMILES input\n", progname);
	printf("  %s --redock 1STP                   # cognate biotin redock from RCSB\n", progname);
	printf("  %s --redock 1GPK -o results/1gpk   # redock with custom output prefix\n", progname);
	printf("  %s protein.pdb drug.sdf -c config.json -o results\n", progname);
	printf("  %s receptor.pdb library.sdf          # multi-molecule SDF\n", progname);
	printf("  %s receptor.pdb ligands/             # directory of files\n", progname);
	printf("  %s receptor.pdb compounds.smi        # SMILES file\n", progname);
	printf("  %s nmr_ensemble.pdb ligand.mol2      # NMR ensemble\n\n", progname);
	printf("Defaults: T=300K, full flexibility, Voronoi contacts, intramolecular ON.\n");
	printf("Redock cache: $FLEXAIDDS_REDOCK_CACHE or ~/.flexaidds/benchmarks/redock/<PDBid>/\n");
}

int main(int argc, char **argv){
	// ── --version ────────────────────────────────────────────────────────
	// Deliberately the FIRST statement in main(): ahead of the try block,
	// ahead of flexaids_rng::init_from_env(), ahead of the FA allocation and
	// the base-path / data-directory resolution below.
	//
	// This placement is a correctness requirement, not a style choice.  The
	// existing --help scan sits ~370 lines further down, and by the time
	// control reaches it main() has already written
	//     base path is '<...>'
	//     auto-detected data directory: '<...>'
	// to STDOUT.  Handling --version there would put two non-key=value lines
	// ahead of the stamp and break the one-key=value-per-line contract that
	// makes the output parseable -- silently, for any reader that does not
	// happen to skip them.
	//
	// Being first also satisfies the other two requirements as a consequence
	// rather than by separate effort: nothing has opened a file yet, so
	// --version works with no input files present; and nothing has been
	// allocated or seeded, so it cannot perturb a run.
	//
	// Scanned across all of argv rather than tested at argv[1] only, matching
	// the --help loop, so `FlexAIDdS <receptor> <ligand> --version` answers
	// instead of starting a dock.
	for (int a = 1; a < argc; ++a) {
		if (strcmp(argv[a], "--version") == 0) {
			flexaids::version::print_build_identity();
			return 0;
		}
	}

  try {
	flexaids_rng::init_from_env();
	int   i,j;

	char remark[MAX_REMARK];
	char tmpremark[MAX_REMARK];
	char dockinp[MAX_PATH__];
	char gainp[MAX_PATH__];
	#ifdef _WIN32
	char *pch;                               // for finding base path
	#endif
	// Initialized so the output prefix is always a valid C string. Every
	// assignment below (legacy argv[4], legacy_files[2], output_prefix) is
	// CONDITIONAL, so an unrecognized invocation used to leave these
	// uninitialized. A `end_strfile ? ... : "flexaid"` guard downstream could
	// never catch that — the address of an array is never null — so the
	// fallback has to live here, in the value.
	char end_strfile[MAX_PATH__] = "flexaid";
	char tmp_end_strfile[MAX_PATH__] = "";

	int memchrom=0;

	time_t sta_timer,end_timer;
	struct tm *sta,*end;
	int sta_val[3],end_val[3];
	long ct; // computational time

	atom *atoms = NULL;
	resid *residue = NULL;
	resid *res_ptr = NULL;
	cfstr cf;
	cfstr* cf_ptr = NULL;
	rot* rotamer = NULL;
	chromosome* chrom = NULL;
	chromosome* chrom_snapshot = NULL;
	genlim* gene_lim = NULL;
	gridpoint* cleftgrid = NULL;

	//flexaid global variables
	FA_Global* FA = NULL;
	GB_Global* GB = NULL;
	VC_Global* VC = NULL;

	try {
		FA = new FA_Global{};
		GB = new GB_Global{};
		VC = new VC_Global{};
	} catch (const std::bad_alloc&) {
		fprintf(stderr,"ERROR: Could not allocate memory for FA || GB || VC\n");
		Terminate(2);
	}
	// Honour FLEXAIDDS_FLAGS=… overlay + mutual-exclusion losers before any
	// later getenv() in scoring/search. Does not remove knobs — losers are
	// unset in the environment so legacy call sites follow the winner.
	flexaidds::flags::apply_to_environ();
	GB->metal_batch_n = 2;  // N=2: safe batch size verified
	// MIF/RefLig/GridPrio non-zero defaults (pointers already NULL via value-init)
	FA->mif_temperature = 300.0f;
	FA->grid_prio_percent = 100.0f;
	FA->reflig_seed_fraction = 0.25f;
	FA->reflig_k_nearest = 10;
	FA->reflig_hetatm_fallback = 1;
	FA->autoflex_enabled = 1;  // auto-flex key binding residues by default
	FA->autoflex_max = 5;

	// calloc (not malloc): FLEXAIDDS_CONTACTS_EPOCH mode never memsets this
	// buffer between vcfunction() calls (see vcfunction.cpp), so it must start
	// all-zero once here to match the epoch-0-means-"never touched" invariant.
	// Legacy mode still memsets it every call, so the zero-fill costs nothing
	// extra there (one-time, at startup, vs. the malloc it replaces).
	// CONTACTS_BUFFER_SIZE (not MAX_ATOM_NUMBER): the trailing slot carries the
	// epoch counter — see flexaid.h.
	FA->contacts = (int*)calloc(CONTACTS_BUFFER_SIZE,sizeof(int));
	if(FA->contacts == NULL){
		fprintf(stderr,"ERROR: Could not allocate memory for contacts\n");
		Terminate(2);
	}

	VC->ptorder = (ptindex*)malloc(MAX_PT*sizeof(ptindex));
	VC->centerpt = (vertex*)malloc(MAX_PT*sizeof(vertex));
	VC->poly = (vertex*)malloc(MAX_POLY*sizeof(vertex));
	VC->cont = (plane*)malloc(MAX_PT*sizeof(plane));
	VC->vedge = (edgevector*)malloc(MAX_POLY*sizeof(edgevector));

	if(!VC->ptorder || !VC->centerpt || !VC->poly ||
	   !VC->cont || !VC->vedge){
		fprintf(stderr,"ERROR: Could not allocate memory for ptorder || centerpt || poly || cont || vedge\n");
		Terminate(2);
	}

	VC->recalc = 1;

	// set minimal default values
	FA->MIN_NUM_ATOM = 1000;
	FA->MIN_NUM_RESIDUE = 250;
	FA->MIN_ROTAMER_LIBRARY_SIZE = 155;
	FA->MIN_ROTAMER = 1;
	FA->MIN_FLEX_BONDS = 5;
	FA->MIN_CLEFTGRID_POINTS = 250;
	FA->MIN_PAR = 6;
	FA->MIN_FLEX_RESIDUE = 5;
	FA->MIN_NORMAL_GRID_POINTS = 250;
	FA->MIN_OPTRES = 1;
	FA->MIN_CONSTRAINTS = 1;

	FA->vindex = 0;
	FA->rotout = 0;
	FA->num_optres = 0;
	FA->nflexbonds = 0;
	FA->normal_grid = NULL;
	FA->supernode = 0;
	FA->eigenvector = NULL;
	FA->psFlexDEENode = NULL;
	FA->FlexDEE_Nodes = 0;
	FA->dee_clash = 0.5;
	FA->intrafraction = 1.0;
	FA->cluster_rmsd = 2.0f;
	FA->use_super_cluster = false;
	FA->rotamer_permeability = 0.8;
	FA->temperature = 0;
	FA->beta = 0.0;
	// Classic FlexAID entropy ranking is the product default when T>0
	// (see config_parser / docs/classic_entropy_ranking.md). Explicit zero
	// here so non-JSON paths never inherit garbage for the emission gate.
	FA->force_cf_rank_emission = false;

	FA->force_interaction=0;
	FA->interaction_factor=5.0;
	// Classic CONFIG path never set soft_wall_cutoff (stayed 0 → hard r^-12 only).
	// Default matches JSON/DatasetRunner (0.40 Å) so CF.wal competes with CF.com;
	// override with SOFTWA 0.0 in CONFIG.inp for pure legacy hard wall.
	FA->soft_wall_cutoff = 0.40f;
	// PoseBust physical-realism clash penalty (opt-in, both legacy + JSON paths).
	// UNCAPPED severity-scaled term added to CF.pb_clash in vcfunction.cpp; sums
	// into the GA fitness (ic2cf.cpp get_apparent_cf_evalue). Fixes the arm-A 0%
	// where capped CF.wal could not overcome unbounded CF.com overpacking.
	// Env-set so the classic --legacy .inp path (no CONFIG keyword) can enable it.
	FA->pb_clash_weight = 0.0;      // OFF by default
	FA->pb_clash_exponent = 3.0;    // steep tail, smooth onset
	FA->pb_clash_ratio = 0.75;      // PoseBusters intermolecular clash default (fraction of summed vdW radii)
	if (const char* e = std::getenv("FLEXAIDDS_PB_CLASH_WEIGHT")) { FA->pb_clash_weight = atof(e); }
	if (const char* e = std::getenv("FLEXAIDDS_PB_CLASH_EXP"))    { double v = atof(e); if (v > 0.0) FA->pb_clash_exponent = v; }
	if (const char* e = std::getenv("FLEXAIDDS_PB_CLASH_RATIO"))  { double v = atof(e); if (v > 0.0) FA->pb_clash_ratio = v; }
	// PoseBust pocket-presence penalty (opt-in). Soft quadratic ramp on the ligand
	// centroid's distance to the nearest receptor heavy atom, so a pose that drifts
	// out of the pocket is steered against during the GA search rather than merely
	// reported afterwards by the post-dock `bust` gate.
	FA->pb_pocket_weight = 0.0;     // OFF by default
	FA->pb_pocket_radius = 6.0;     // free radius (A) — no penalty inside it
	if (const char* e = std::getenv("FLEXAIDDS_PB_POCKET_WEIGHT")) { FA->pb_pocket_weight = atof(e); }
	if (const char* e = std::getenv("FLEXAIDDS_PB_POCKET_RADIUS")) { double v = atof(e); if (v > 0.0) FA->pb_pocket_radius = v; }
	// VCT com normalization (Lever 2): divide CF.com by contact count (intensive
	// score) instead of the extensive sum, rescaled by VCT_NREF=100. Fixes the
	// arm-A overpacking where a dense non-native pose out-scores native purely by
	// burying more surface (more contacts). Read from env so the classic --legacy
	// .inp path (no CONFIG "scoring" JSON block) can enable it; the JSON path sets
	// it in config_parser.cpp:69. Default OFF (extensive) — unchanged legacy behavior.
	FA->vct_normalize_contacts = std::getenv("FLEXAIDDS_VCT_NORM") ? 1 : 0;
	FA->atm_cnt=0;
	FA->atm_cnt_real=0;
	FA->res_cnt=0;
	FA->nors=0;

	FA->htpmode=false;
	FA->nrg_suite=0;
	FA->nrg_suite_timeout=60;
	FA->translational=0;
	FA->refstructure=0;
	FA->omit_buried=0;
	FA->assume_folded=0;
	FA->natural_deltaG=0.0;
	FA->is_protein=1;

	FA->delta_angstron=0.25;
	FA->delta_angle=5.0;
	FA->delta_dihedral=5.0;
	FA->delta_flexible=10.0;
	FA->delta_index=1.0;
	FA->max_results=10;
	FA->deelig_flex = 0;
	FA->resligand = NULL;
	FA->useacs = 0;
	FA->acsweight = 1.0;

	GB->outgen=0;
	GB->entropy_weight=0.5;
	GB->entropy_interval=0;
	GB->use_shannon=0;
	FA->num_grd=0;
	FA->exclude_het=0;
	FA->remove_water=1;
	FA->normalize_area=0;

	FA->recalci=0;
	FA->skipped=0;
	FA->clashed=0;

	FA->spacer_length=0.375;
	FA->opt_grid=0;

	FA->pbloops=1;
	FA->bloops=2;

	FA->rotobs=0;
	FA->contributions=NULL;
	FA->output_scored_only=0;
	FA->score_ligand_only=0;
	FA->permeability=1.0;
	FA->intramolecular=1;
	FA->solventterm=0.0f;
	FA->use_elec=0;
	FA->use_memetic=0;
	FA->memetic_armed_at_gen=0;
	FA->dielectric=4.0f;

	FA->use_gist=0;
	FA->gist_dg_file[0]='\0';
	FA->gist_dens_file[0]='\0';
	FA->gist_weight=1.0f;
	FA->gist_dg_cutoff=1.0f;
	FA->gist_rho_cutoff=4.8f;
	FA->gist_divisor=2.0f;
	FA->gist_evaluator=NULL;

	FA->use_hbond=0;
	FA->use_hbond_search=0;
	FA->use_hbond_rank=0;
	FA->hbond_rank_rescore=0;
	FA->hbond_weight=-2.5;
	// v116: tuneable via env — routed through ProtocolConfig (typed adapter).
	FA->hbond_weight = flexaids::ProtocolConfig::from_env().hbond_weight;
	FA->hbond_optimal_dist=2.8;
	FA->hbond_optimal_angle=180.0;
	FA->hbond_sigma_dist=0.4;
	FA->hbond_sigma_angle=30.0;
	FA->hbond_salt_bridge_weight=-5.0;

	FA->use_metal_coord=0;
	FA->metal_coord_weight=1.0;
	FA->metal_coord_sigma=0.45;
	FA->metal_coord_cn_weight=0.5;
	FA->tencom_weight=0.0f;

	FA->useflexdee=0;
	FA->num_constraints=0;

	FA->npar=0;

	FA->mov[0] = NULL;
	FA->mov[1] = NULL;
	strncpy(FA->clustering_algorithm,"CF",sizeof(FA->clustering_algorithm)-1); FA->clustering_algorithm[sizeof(FA->clustering_algorithm)-1]='\0';
	strncpy(FA->vcontacts_self_consistency,"MAX",sizeof(FA->vcontacts_self_consistency)-1); FA->vcontacts_self_consistency[sizeof(FA->vcontacts_self_consistency)-1]='\0';
	FA->vcontacts_planedef = 'X';

	// ── Determine base path from executable location ──
	// Prefer the real path of argv[0] so Homebrew/PATH symlinks still resolve
	// to Cellar/.../bin (where runtime data is installed), not /opt/homebrew/bin.
#ifndef _WIN32
	{
		// Resolve argv[0] with a glibc-allocated buffer (second arg NULL) rather
		// than a fixed MAX_PATH__ stack buffer. _FORTIFY_SOURCE rejects any
		// realpath() destination smaller than PATH_MAX (4096): __realpath_chk
		// calls __chk_fail() when resolvedlen < PATH_MAX, before resolving —
		// so a MAX_PATH__ (255) buffer aborts ("buffer overflow detected")
		// UNCONDITIONALLY on fortified glibc builds, every invocation, any
		// path length. Passing NULL makes glibc size the buffer itself.
		// Do not reintroduce a fixed-size destination here. Falls back to argv[0].
		char* rp = realpath(argv[0], NULL);
		const char* src = (rp != NULL) ? rp : argv[0];
		// strrchr(const char*) returns const char* — keep a local const pointer
		// instead of assigning into the legacy char* pch variable.
		const char* slash = strrchr(src, '/');
		if (slash != NULL) {
			size_t n = (size_t)(slash - src);
			if (n >= MAX_PATH__) n = MAX_PATH__ - 1;
			memcpy(FA->base_path, src, n);
			FA->base_path[n] = '\0';
		} else {
			strncpy(FA->base_path, ".", MAX_PATH__ - 1);
			FA->base_path[MAX_PATH__ - 1] = '\0';
		}
		free(rp);  // free(NULL) is a no-op
	}
#else
	pch = strrchr(argv[0], '\\');
	if (pch == NULL) {
		pch = strrchr(argv[0], '/');
	}
	if (pch != NULL) {
		for (i = 0; i < (int)(pch - argv[0]); i++) {
			FA->base_path[i] = argv[0][i];
			FA->base_path[i + 1] = '\0';
		}
	} else {
		strncpy(FA->base_path, ".", MAX_PATH__ - 1);
		FA->base_path[MAX_PATH__ - 1] = '\0';
	}
#endif //_WIN32

	printf("base path is '%s'\n", FA->base_path);

		// ── Auto-detect data directory ──────────────────────────────────────
		// Priority:
		//   1. FLEXAIDDS_DATA_DIR env (Homebrew wrappers set this)
		//   2. base_path itself (MC_st0r5.2_6.dat next to the binary)
		//   3. base_path/../WRK (dev layout: binary in build/, data in WRK/)
		//   4. base_path/../share/flexaidds (Homebrew share layout)
		// The --data-dir flag (parsed below) overrides this auto-detection.
		{
			char probe[MAX_PATH__];
			int found = 0;
			const flexaids::ProtocolConfig proto = flexaids::ProtocolConfig::from_env();
			if (!proto.data_dir.empty()) {
				snprintf(probe, MAX_PATH__, "%s/MC_st0r5.2_6.dat", proto.data_dir.c_str());
				FILE* fp = fopen(probe, "r");
				if (fp) {
					fclose(fp);
					strncpy(FA->dependencies_path, proto.data_dir.c_str(), MAX_PATH__ - 1);
					FA->dependencies_path[MAX_PATH__ - 1] = '\0';
					printf("data directory from FLEXAIDDS_DATA_DIR: '%s'\n", FA->dependencies_path);
					found = 1;
				}
			}
			if (!found) {
				snprintf(probe, MAX_PATH__, "%s/MC_st0r5.2_6.dat", FA->base_path);
				FILE* fp = fopen(probe, "r");
				if (fp) {
					fclose(fp);
					found = 1; // data co-located with binary
				}
			}
			if (!found) {
				const char* suffixes[] = { "/../WRK", "/../share/flexaidds", "/share/flexaidds" };
				for (size_t si = 0; si < sizeof(suffixes) / sizeof(suffixes[0]); ++si) {
					snprintf(probe, MAX_PATH__, "%s%s/MC_st0r5.2_6.dat", FA->base_path, suffixes[si]);
					FILE* fp = fopen(probe, "r");
					if (fp) {
						fclose(fp);
						snprintf(FA->dependencies_path, MAX_PATH__, "%s%s", FA->base_path, suffixes[si]);
						printf("auto-detected data directory: '%s'\n", FA->dependencies_path);
						found = 1;
						break;
					}
				}
			}
		}


	// ── CLI argument parsing ──────────────────────────────────────────────
	bool legacy_mode = false;
	bool use_rigid = false;
	bool use_folded = false;
	bool use_screen = false;
	bool use_screen_dock = false;
	int  screen_top_n = 100;
	std::string screen_receptor_path;  // populated in auto-detect path for --screen
	std::string screen_ligand_path;
	std::string screen_target_mol2;
	bool use_parallel_dock = false;
	int  parallel_dock_regions = 128;
	bool use_campaign = false;
	bool use_coarse_prefilter = false;
	int  coarse_prefilter_top_n = 100;
	std::string config_path;
	std::string output_prefix = "flexaid_out";
	std::string cached_grid_path;  // Strategy A: .rrg grid cache path from "grid_file" JSON key
	double user_conc_M = 1.0;  // P3: per-run concentration for grand canonical (default 1M)

	if (argc < 2) {
		print_usage(argv[0]);
		Terminate(1);
	}

	// Check for --help
	for (int a = 1; a < argc; a++) {
		if (strcmp(argv[a], "-h") == 0 || strcmp(argv[a], "--help") == 0) {
			print_usage(argv[0]);
			Terminate(0);
		}
	}

	// ── Idiotproof argument parsing ──────────────────────────────────────────
	// Accepts files in any order. Auto-detects receptor vs ligand.
	// Handles: PDB receptor, MOL2/SDF/PDB ligand, SMILES, JSON config.

	// Check for --benchmark mode
	if (strcmp(argv[1], "--benchmark") == 0) {
		if (argc < 3) {
			fprintf(stderr, "ERROR: --benchmark requires a dataset name\n");
			fprintf(stderr, "  Available: astex, astex_nonnative, hap2, casf2016, posebusters,\n");
			fprintf(stderr, "             dude, bindingdb_itc, sampl6, sampl7, pdbbind, all\n");
			fprintf(stderr, "  Also: doi:<DOI>, pdb_list:<file>\n");
			Terminate(1);
		}
		// Forward to benchmark_datasets via argv exec (no shell).
		// Single child only -- never dual-launch.
		std::vector<std::string> bench_argv;
		bench_argv.reserve(static_cast<size_t>(argc));
		bench_argv.emplace_back("benchmark_datasets");
		for (int a = 1; a < argc; a++) {
			if (!flexaids::shell_exec::is_safe_exec_path(argv[a])) {
				fprintf(stderr,
				        "ERROR: unsafe argument for --benchmark (NUL/newline/control)\n");
				Terminate(1);
			}
			bench_argv.emplace_back(argv[a]);
		}
		printf("Launching benchmark runner (argv exec, no shell)\n");
		int ret = flexaids::shell_exec::run_argv(bench_argv);
		Terminate(ret < 0 ? 127 : ret);
	}

	// Check for --legacy mode first
	if (strcmp(argv[1], "--legacy") == 0) {
		if (argc < 5) {
			fprintf(stderr, "ERROR: --legacy requires 3 arguments: <config.inp> <ga.inp> <output_prefix>\n");
			Terminate(1);
		}
		legacy_mode = true;
		strncpy(dockinp, argv[2], MAX_PATH__-1); dockinp[MAX_PATH__-1]='\0';
		strncpy(gainp, argv[3], MAX_PATH__-1); gainp[MAX_PATH__-1]='\0';
		strncpy(end_strfile, argv[4], MAX_PATH__-1); end_strfile[MAX_PATH__-1]='\0';
		strncpy(FA->rrgfile, end_strfile, MAX_PATH__-1); FA->rrgfile[MAX_PATH__-1]='\0';
	}
	else {
		// ── Auto-detect mode: scan ALL arguments, classify each ──
		std::string receptor_path;
		std::string ligand_path;
		std::string redock_pdb_id;
		// Hoisted out of the --redock block below so the receipt writer can name
		// the run by its prepared (upper-cased) PDB id. Declaration only — the
		// value and every use of it are unchanged.
		std::string redock_prepared_id;
		bool user_set_output = false;
		std::vector<std::string> legacy_files;

		for (int a = 1; a < argc; a++) {
			std::string arg(argv[a]);

			// Skip flags and their values
			if (arg == "-c" || arg == "--config") {
				if (a + 1 < argc) config_path = argv[++a];
				continue;
			}
			if (arg == "-o" || arg == "--output") {
				if (a + 1 < argc) {
					output_prefix = argv[++a];
					user_set_output = true;
				}
				continue;
			}
			if (arg == "--data-dir") {
				if (a + 1 < argc) {
					strncpy(FA->dependencies_path, argv[++a], MAX_PATH__-1);
					FA->dependencies_path[MAX_PATH__-1] = '\0';
				} else {
					fprintf(stderr, "ERROR: --data-dir requires a directory path\n");
					Terminate(1);
				}
				continue;
			}
			if (arg == "--redock") {
				if (a + 1 >= argc) {
					fprintf(stderr, "ERROR: --redock requires a RCSB PDB ID (e.g. --redock 1STP)\n");
					Terminate(1);
				}
				if (!redock_pdb_id.empty()) {
					fprintf(stderr, "ERROR: --redock specified more than once\n");
					Terminate(1);
				}
				redock_pdb_id = argv[++a];
				continue;
			}
			if (arg == "--backend") {
				if (a + 1 < argc) {
					std::string be = argv[++a];
					hw::Backend requested;
					if (be == "cpu") {
						requested = hw::Backend::SCALAR;
					} else if (be == "metal") {
						requested = hw::Backend::METAL;
					} else if (be == "webgpu") {
						requested = hw::Backend::WEBGPU;
					} else {
						fprintf(stderr, "ERROR: --backend must be one of cpu|metal|webgpu (got '%s')\n", be.c_str());
						Terminate(1);
						continue;
					}
					if (requested != hw::Backend::SCALAR &&
					    !hw::UnifiedHardwareDispatch::instance().is_available(requested)) {
						fprintf(stderr, "WARNING: --backend %s requested but not available on this build/host; "
						                "falling back to CPU.\n", be.c_str());
						requested = hw::Backend::SCALAR;
					}
					hw::UnifiedHardwareDispatch::instance().set_override(requested);
				} else {
					fprintf(stderr, "ERROR: --backend requires an argument (cpu|metal|webgpu)\n");
					Terminate(1);
				}
				continue;
			}
			if (arg == "--rigid")  { use_rigid = true;  continue; }
			if (arg == "--screen") { use_screen = true; continue; }
			if (arg == "--screen-dock") { use_screen_dock = true; continue; }
			if (arg == "--screen-top-n") {
				if (a + 1 < argc) screen_top_n = std::atoi(argv[++a]);
				continue;
			}
			if (arg == "--screen-target-mol2") {
				if (a + 1 < argc) screen_target_mol2 = argv[++a];
				continue;
			}
			if (arg == "--parallel-dock") { use_parallel_dock = true; continue; }
			if (arg == "--parallel-dock-regions") {
				if (a + 1 < argc) parallel_dock_regions = std::atoi(argv[++a]);
				continue;
			}
			if (arg == "--campaign") { use_campaign = true; continue; }
			if (arg == "--coarse-prefilter") { use_coarse_prefilter = true; continue; }
			if (arg == "--coarse-prefilter-top-n") {
				if (a + 1 < argc) coarse_prefilter_top_n = std::atoi(argv[++a]);
				continue;
			}
			if (arg == "--folded") { use_folded = true; continue; }
			if (arg == "--conc" || arg == "--concentration") {
				if (a + 1 < argc) user_conc_M = std::atof(argv[++a]);
				continue;
			}
			if (arg == "-h" || arg == "--help") { print_usage(argv[0]); Terminate(0); }

			// Classify this positional argument
			std::string role = detect_file_role(arg);

			if (role == "receptor") {
				if (receptor_path.empty()) {
					receptor_path = arg;
				} else {
					fprintf(stderr, "WARNING: Multiple receptor files detected.\n");
					fprintf(stderr, "  Using: %s\n  Ignoring: %s\n",
					        receptor_path.c_str(), arg.c_str());
				}
			} else if (role == "ligand" || role == "smiles") {
				if (ligand_path.empty()) {
					ligand_path = arg;
				} else {
					fprintf(stderr, "WARNING: Multiple ligand inputs detected.\n");
					fprintf(stderr, "  Using: %s\n  Ignoring: %s\n",
					        ligand_path.c_str(), arg.c_str());
				}
			} else if (role == "config") {
				config_path = arg;
			} else if (role == "legacy") {
				legacy_files.push_back(arg);
			} else {
				// Unknown file — try to be helpful
				if (std::filesystem::exists(arg)) {
					fprintf(stderr, "WARNING: Cannot determine role of '%s'.\n", arg.c_str());
					fprintf(stderr, "  Supported: .pdb/.cif/.mmcif (receptor), .mol2/.sdf/.mol (ligand), .json (config)\n");
				} else {
					fprintf(stderr, "ERROR: File not found and not valid SMILES: '%s'\n", arg.c_str());
					print_usage(argv[0]);
					Terminate(1);
				}
			}
		}

		// ── Cognate redock from RCSB (--redock PDBid) ─────────────────────
		// Downloads the deposit, extracts the cognate ligand, writes an apo
		// receptor with that ligand stripped, then continues in direct mode.
		if (!redock_pdb_id.empty()) {
			if (!receptor_path.empty() || !ligand_path.empty()) {
				fprintf(stderr,
				        "ERROR: --redock cannot be combined with explicit receptor/ligand files\n");
				Terminate(1);
			}
			if (!legacy_files.empty()) {
				fprintf(stderr, "ERROR: --redock cannot be combined with legacy .inp inputs\n");
				Terminate(1);
			}
			std::string& prepared_id = redock_prepared_id;
			if (!prepare_redock_from_rcsb(redock_pdb_id, receptor_path, ligand_path, prepared_id)) {
				Terminate(1);
			}
			if (!user_set_output) {
				output_prefix = prepared_id + "_redock";
			}
		}

		// Legacy auto-detect: if we got legacy .inp files instead of PDB/MOL2
		if (!legacy_files.empty() && receptor_path.empty() && ligand_path.empty()) {
			if (legacy_files.size() >= 2) {
				legacy_mode = true;
				strncpy(dockinp, legacy_files[0].c_str(), MAX_PATH__-1); dockinp[MAX_PATH__-1]='\0';
				strncpy(gainp, legacy_files[1].c_str(), MAX_PATH__-1); gainp[MAX_PATH__-1]='\0';
				if (legacy_files.size() >= 3) {
					strncpy(end_strfile, legacy_files[2].c_str(), MAX_PATH__-1);
				} else {
					strncpy(end_strfile, output_prefix.c_str(), MAX_PATH__-1);
				}
				end_strfile[MAX_PATH__-1]='\0';
				strncpy(FA->rrgfile, end_strfile, MAX_PATH__-1); FA->rrgfile[MAX_PATH__-1]='\0';
			} else {
				fprintf(stderr, "ERROR: Legacy mode requires at least 2 .inp files.\n");
				print_usage(argv[0]);
				Terminate(1);
			}
		}

		// Validate we have what we need for direct mode
		if (!legacy_mode) {
			if (receptor_path.empty()) {
				fprintf(stderr, "ERROR: No receptor file detected.\n");
				fprintf(stderr, "  Provide a .pdb or .cif file containing a protein or nucleic acid,\n");
				fprintf(stderr, "  or use --redock <PDBid> for automatic cognate redocking from RCSB.\n\n");
				print_usage(argv[0]);
				Terminate(1);
			}
			if (ligand_path.empty()) {
				fprintf(stderr, "ERROR: No ligand input detected.\n");
				fprintf(stderr, "  Provide a .mol2, .sdf, .mol, or .pdb ligand file,\n");
				fprintf(stderr, "  pass a SMILES string directly, or use --redock <PDBid>.\n\n");
				print_usage(argv[0]);
				Terminate(1);
			}

			printf("Receptor: %s\n", receptor_path.c_str());
			printf("Ligand:   %s\n", ligand_path.c_str());
			if (use_screen) {
				screen_receptor_path = receptor_path;
				screen_ligand_path   = ligand_path;
			}
		}

		// ── Apply config ──
		if (!legacy_mode) {
			json::Value config = load_config(config_path);
			if (use_rigid) config = json::merge(config, flexaid_rigid_overrides());
			if (use_folded) {
				using V = json::Value;
				using O = json::Object;
				config = json::merge(config, V(O{{"advanced", V(O{{"assume_folded", V(true)}})}}));
			}
			// Snapshot protocol once for apply_config (no mid-apply getenv dual path).
			const flexaids::ProtocolConfig apply_proto = flexaids::ProtocolConfig::from_env();
			apply_config(config, FA, GB, &apply_proto);
			if (GB->seed == 0) {
				std::uint64_t env = 0;
				if (flexaids_rng::has_master_seed()) {
					const auto ms = static_cast<unsigned int>(flexaids_rng::master_seed() & 0x7fffffffu);
					GB->seed = ms ? static_cast<int>(ms) : 1;
				} else if (flexaids_rng::env_seed(env)) {
					const auto es = static_cast<unsigned int>(env & 0x7fffffffu);
					GB->seed = es ? static_cast<int>(es) : 1;
					flexaids_rng::set_master_seed(static_cast<std::uint64_t>(GB->seed));
				} else {
					const std::string key = !redock_pdb_id.empty()
					    ? (std::string("redock:") + redock_pdb_id)
					    : (receptor_path + "|" + ligand_path);
					GB->seed = flexaids_rng::deterministic_seed_from_key(key.c_str());
					flexaids_rng::set_master_seed(static_cast<std::uint64_t>(GB->seed));
					printf("[SEED] ga.seed=0; assigned deterministic seed %d from '%s'\n",
					       GB->seed, key.c_str());
				}
			} else {
				flexaids_rng::set_master_seed(static_cast<std::uint64_t>(GB->seed));
			}

			// ── Strategy A: extract grid cache path from JSON "grid_file" key ──────
			if (config.contains("grid_file")) {
				cached_grid_path = config["grid_file"].as_string("");
				if (!cached_grid_path.empty())
					fprintf(stderr, "[GRID-CACHE] grid_file=%s\n", cached_grid_path.c_str());
			}

			printf("%sFlexAIDdS config:%s %sT=%uK%s, ligand_flex=%s, intramolecular=%s, "
			       "scoring=%s, intermolecular_clash_ratio=%.3f\n",
				tui::muted(), tui::reset(), tui::T(), FA->temperature, tui::reset(),
				FA->deelig_flex ? "ON" : "OFF",
				FA->intramolecular ? "ON" : "OFF",
				FA->complf,
				FA->intermolecular_clash_ratio);
		}


		// Set output prefix
		strncpy(end_strfile, output_prefix.c_str(), MAX_PATH__ - 1);
		end_strfile[MAX_PATH__ - 1] = '\0';
		strncpy(FA->rrgfile, end_strfile, MAX_PATH__-1); FA->rrgfile[MAX_PATH__-1]='\0';

		// ── RUN_RECEIPT.json when --redock's -o names an existing directory ──
		// Additive only. `-o` stays a file prefix; nothing above or below this
		// block changes. Three conditions must all hold, so every pre-existing
		// invocation reaches this point and does nothing: --redock was used, -o
		// was given explicitly, and that path already exists as a directory.
		// Written here — after apply_config, before any docking — to match the
		// engine's own semantics: the receipt is a statement of intent, complete
		// before a single pose exists (RUN_RECEIPT_CONTRACT.md §5).
#if defined(FLEXAIDDS_ENABLE_REDOCK)
		if (!redock_pdb_id.empty() && user_set_output) {
			std::error_code rec_ec;
			if (std::filesystem::is_directory(output_prefix, rec_ec) && !rec_ec) {
				write_redock_run_receipt(
				    output_prefix,
				    redock_prepared_id.empty() ? redock_pdb_id : redock_prepared_id,
				    static_cast<double>(FA->temperature),
				    GB->num_chrom,
				    GB->max_generations);
			}
		}
#endif

		// GA input not used in direct mode
		dockinp[0] = '\0';
		gainp[0] = '\0';

		// ── Strategy B — multi-ligand batch dispatch ──────────────────────────
		// When FLEXAIDDS_LIGAND_BATCH=<dir> is set, dock every *.sdf in that
		// directory against this receptor by re-exec'ing this binary once per
		// ligand (see batch_exec_child above for the rationale).  The single
		// ligand argument parsed above is ignored — the batch directory is the
		// ligand source.  Each child writes to <output_dir>/<sdf-stem>* and
		// reuses the receptor grid through the Strategy A cache.  Absent env var
		// → this block is skipped and the normal single-ligand path runs.
#ifndef _WIN32
		{
			const char* batch_env = std::getenv("FLEXAIDDS_LIGAND_BATCH");
			if (batch_env && batch_env[0] != '\0' &&
			    std::filesystem::is_directory(batch_env)) {

				// Enumerate *.sdf, sorted lexicographically for determinism.
				std::vector<std::string> sdfs;
				for (const auto& de : std::filesystem::directory_iterator(batch_env)) {
					if (!de.is_regular_file()) continue;
					std::string ext = de.path().extension().string();
					for (char& c : ext)
						c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
					if (ext == ".sdf") sdfs.push_back(de.path().string());
				}
				std::sort(sdfs.begin(), sdfs.end());

				if (sdfs.empty()) {
					fprintf(stderr, "[BATCH] no .sdf files in %s — nothing to dock\n",
					        batch_env);
					Terminate(0);
				}

				// Output directory derived from the existing -o prefix.
				std::filesystem::path out_base(output_prefix);
				const std::string out_dir = out_base.has_parent_path()
				                          ? out_base.parent_path().string()
				                          : std::string(".");

				const int M = static_cast<int>(sdfs.size());
				fprintf(stderr,
				        "[BATCH] %d ligand(s) vs receptor %s "
				        "(grid built once, shared via FLEXAIDDS_GRID_CACHE_DIR)\n",
				        M, receptor_path.c_str());

				int n_ok = 0;
				for (int li = 0; li < M; ++li) {
					const std::string stem =
					    std::filesystem::path(sdfs[li]).stem().string();
					const std::string child_prefix = out_dir + "/" + stem;
					fprintf(stderr, "[BATCH] ligand %d/%d: %s\n",
					        li + 1, M, stem.c_str());

					std::vector<std::string> cargs;
					cargs.push_back(argv[0]);
					cargs.push_back(receptor_path);
					cargs.push_back(sdfs[li]);
					cargs.push_back("-o");
					cargs.push_back(child_prefix);
					if (!config_path.empty()) {
						cargs.push_back("-c");
						cargs.push_back(config_path);
					}
					if (FA->dependencies_path[0] != '\0') {
						cargs.push_back("--data-dir");
						cargs.push_back(FA->dependencies_path);
					}

					const int rc = batch_exec_child(argv[0], cargs);
					if (rc == 0) {
						++n_ok;
					} else {
						fprintf(stderr,
						        "[BATCH] ligand %s exited with code %d\n",
						        stem.c_str(), rc);
					}
				}
				fprintf(stderr, "%s[BATCH]%s complete: %d/%d ligand(s) succeeded\n", tui::err::strawberry(), tui::err::reset(),
				        n_ok, M);
				Terminate(0);
			}
		}
#endif // _WIN32

		// Direct loading pipeline — use auto-detected paths
		const char* receptor_file = receptor_path.c_str();
		const char* ligand_file   = ligand_path.c_str();

		// ── 1. Interaction matrix ──
		{
			char emat[MAX_PATH__];
			if (!strcmp(FA->dependencies_path, "")) {
				strcpy(emat, FA->base_path);
			} else {
				strcpy(emat, FA->dependencies_path);
			}
#ifdef _WIN32
			strcat(emat, "\\MC_st0r5.2_6.dat");
#else
			strcat(emat, "/MC_st0r5.2_6.dat");
#endif
			printf("interaction matrix is <%s>\n", emat);
			read_emat(FA, emat);
		}

		// ── 2. Check if target is RNA ──
		if (rna_structure(const_cast<char*>(receptor_file))) {
			printf("target molecule is a RNA structure\n");
			FA->is_protein = 0;
		}

		// ── 3. Definition of types ──
		char deftyp[MAX_PATH__];
		{
			if (!strcmp(FA->dependencies_path, "")) {
				strcpy(deftyp, FA->base_path);
			} else {
				strcpy(deftyp, FA->dependencies_path);
			}
			if (FA->is_protein) {
#ifdef _WIN32
				strcat(deftyp, "\\AMINO.def");
#else
				strcat(deftyp, "/AMINO.def");
#endif
			} else {
#ifdef _WIN32
				strcat(deftyp, "\\NUCLEOTIDES.def");
#else
				strcat(deftyp, "/NUCLEOTIDES.def");
#endif
			}
			printf("definition of types is <%s>\n", deftyp);
		}

		// ── 4. Read receptor PDB ──
		{
			// Create temporary cleaned PDB
			char tmpprotname[MAX_PATH__];
			int random_num = static_cast<int>(std::random_device{}() % 900000 + 100000);
			const std::string tmpdir = std::filesystem::temp_directory_path().string();
			snprintf(tmpprotname, MAX_PATH__, "%s/flexaid_receptor_%d.pdb",
			         tmpdir.c_str(), random_num);

			modify_pdb(const_cast<char*>(receptor_file), tmpprotname, FA->exclude_het, FA->remove_water, FA->is_protein,
			           FA->keep_ions, FA->keep_structural_waters, FA->structural_water_bfactor_max);
			read_pdb(FA, &atoms, &residue, tmpprotname);
			remove(tmpprotname);
		}

		residue[FA->res_cnt].latm[0] = FA->atm_cnt;
		for (int k = 1; k <= FA->res_cnt; k++) {
			FA->atm_cnt_real += residue[k].latm[0] - residue[k].fatm[0] + 1;
		}

		calc_center(FA, atoms, residue);

		if (FA->is_protein) {
			residue_conect(FA, atoms, residue, deftyp);
		}
		assign_types(FA, atoms, residue, deftyp);

		// ── 5. Read ligand (auto-detect: SMILES / SDF / MOL2 / PDB) ──
		// ProcessLigand handles format detection, validation, ring
		// perception, aromaticity, SYBYL typing, and failsafe fallback.
		{
			int lig_ok = 0;
			std::string lig_input(ligand_file);

			// Auto-detect: is this a file path or a SMILES string?
			bool is_file = std::filesystem::exists(lig_input);
			bool is_smiles = false;

			if (!is_file) {
				// Not a file — treat as SMILES string if it contains
				// typical SMILES characters and no whitespace/path separators
				bool has_path_chars = (lig_input.find('/') != std::string::npos ||
				                      lig_input.find('\\') != std::string::npos);
				if (!has_path_chars && !lig_input.empty()) {
					is_smiles = true;
					printf("Ligand input detected as SMILES: %s\n", ligand_file);
				} else {
					fprintf(stderr, "ERROR: Ligand file not found: %s\n", ligand_file);
					Terminate(2);
				}
			}

			if (is_smiles) {
				// SMILES → ProcessLigand pipeline → BonMol
				// Note: SMILES provides topology only (no 3D coords).
				// ProcessLigand validates, perceives rings, assigns types.
				// For docking, 3D coordinates are required — user should
				// provide SDF/MOL2 with coords, or use external conformer
				// generation (RDKit Python, OpenBabel) first.
				bonmol::ProcessOptions opts;
				opts.input  = lig_input;
				opts.format = bonmol::InputFormat::SMILES;

				bonmol::ProcessLigand pl;
				auto result = pl.run(opts);

				if (!result.success) {
					fprintf(stderr, "ERROR: ProcessLigand failed for SMILES '%s': %s\n",
					        ligand_file, result.error.c_str());
					Terminate(2);
				}

				// Build 3D coordinates from topology
				printf("Building 3D coordinates from SMILES topology...\n");
				bonmol::CoordBuilderOptions cb_opts;
				if (!bonmol::build_3d_coords(result.mol, cb_opts)) {
					fprintf(stderr, "ERROR: Failed to generate 3D coordinates from SMILES.\n");
					Terminate(2);
				}

				printf("ProcessLigand (SMILES): %d atoms, %d rings (%d aromatic), "
				       "%d rotatable bonds, MW=%.1f\n",
				       result.num_heavy_atoms, result.num_rings,
				       result.num_arom_rings, result.num_rot_bonds,
				       result.molecular_weight);

				// Write temporary MOL2 from BonMol and read back through standard path
				char tmp_mol2[MAX_PATH__];
				const std::string tmpdir = std::filesystem::temp_directory_path().string();
				snprintf(tmp_mol2, MAX_PATH__, "%s/flexaid_smiles_%d.mol2",
				         tmpdir.c_str(), static_cast<int>(std::random_device{}() % 900000 + 100000));

				{
					FILE* fp = fopen(tmp_mol2, "w");
					if (!fp) {
						fprintf(stderr, "ERROR: Cannot write temp MOL2 for SMILES ligand.\n");
						Terminate(2);
					}

					const auto& m = result.mol;
					int na = m.num_atoms();
					int nb = m.num_bonds();

					fprintf(fp, "@<TRIPOS>MOLECULE\nLIG\n%d %d 1 0 0\nSMALL\nNO_CHARGES\n\n", na, nb);

					// Map SYBYL type codes to strings
					auto sybyl_str = [](int t) -> const char* {
						switch (t) {
							case 1:  return "C.3";
							case 2:  return "C.2";
							case 3:  return "C.ar";
							case 4:  return "N.3";
							case 5:  return "N.2";
							case 6:  return "N.ar";
							case 7:  return "N.am";
							case 8:  return "N.pl3";
							case 9:  return "N.4";
							case 10: return "O.3";
							case 11: return "O.2";
							case 12: return "O.co2";
							case 13: return "F";
							case 14: return "Cl";
							case 15: return "Br";
							case 16: return "S.3";
							case 17: return "S.2";
							case 18: return "S.O";
							case 19: return "S.O2";
							case 20: return "P.3";
							case 21: return "I";
							case 22: return "H";
							default: return "Du";
						}
					};

					// Element symbol from enum
					auto elem_str = [](bonmol::Element e) -> const char* {
						switch (e) {
							case bonmol::Element::H:  return "H";
							case bonmol::Element::C:  return "C";
							case bonmol::Element::N:  return "N";
							case bonmol::Element::O:  return "O";
							case bonmol::Element::F:  return "F";
							case bonmol::Element::P:  return "P";
							case bonmol::Element::S:  return "S";
							case bonmol::Element::Cl: return "Cl";
							case bonmol::Element::Br: return "Br";
							case bonmol::Element::I:  return "I";
							default: return "X";
						}
					};

					fprintf(fp, "@<TRIPOS>ATOM\n");
					for (int i = 0; i < na; i++) {
						fprintf(fp, "%6d %4s %9.4f %9.4f %9.4f %6s %3d LIG %8.4f\n",
						        i + 1,
						        elem_str(m.atoms[i].element),
						        m.coords(0, i), m.coords(1, i), m.coords(2, i),
						        sybyl_str(m.atoms[i].sybyl_type),
						        1,
						        m.atoms[i].partial_charge);
					}

					fprintf(fp, "@<TRIPOS>BOND\n");
					for (int i = 0; i < nb; i++) {
						const char* bt = "1";
						switch (m.bonds[i].order) {
							case bonmol::BondOrder::SINGLE:   bt = "1"; break;
							case bonmol::BondOrder::DOUBLE:   bt = "2"; break;
							case bonmol::BondOrder::TRIPLE:   bt = "3"; break;
							case bonmol::BondOrder::AROMATIC: bt = "ar"; break;
						}
						fprintf(fp, "%6d %5d %5d %s\n",
						        i + 1, m.bonds[i].atom_i + 1, m.bonds[i].atom_j + 1, bt);
					}

					fclose(fp);
				}

				// Read the generated MOL2 through the standard FlexAID reader
				printf("read ligand MOL2 (from SMILES) <%s>\n", tmp_mol2);
				lig_ok = read_mol2_ligand(FA, &atoms, &residue, tmp_mol2);
				remove(tmp_mol2);

				if (!lig_ok) {
					fprintf(stderr, "ERROR: Failed to process SMILES-derived ligand.\n");
					Terminate(2);
				}
			} else {
				// File input — detect format from extension
				const char* ext = strrchr(ligand_file, '.');
				bool is_sdf = false;
				if (ext) {
					is_sdf = (strcmp(ext, ".sdf") == 0 || strcmp(ext, ".SDF") == 0 ||
					          strcmp(ext, ".mol") == 0 || strcmp(ext, ".MOL") == 0);
				}

				// Run ProcessLigand for validation + typing enrichment
				// (failsafe: if ProcessLigand fails, fall back to raw readers)
				bonmol::ProcessOptions opts;
				opts.input  = lig_input;
				opts.format = is_sdf ? bonmol::InputFormat::SDF
				                     : bonmol::InputFormat::MOL2;

				bonmol::ProcessLigand pl;
				auto result = pl.run(opts);

				if (result.success) {
					printf("ProcessLigand: %d atoms, %d rings (%d aromatic), "
					       "%d rotatable bonds, MW=%.1f\n",
					       result.num_heavy_atoms, result.num_rings,
					       result.num_arom_rings, result.num_rot_bonds,
					       result.molecular_weight);
				} else {
					printf("ProcessLigand info: %s (continuing with raw reader)\n",
					       result.error.c_str());
				}

				// Always use the existing readers for FlexAID atom/resid population
				// (ProcessLigand enrichment is diagnostic; the readers do the
				// actual struct population that gaboom.cpp expects)
				if (is_sdf) {
					printf("read ligand SDF <%s>\n", ligand_file);
					lig_ok = read_sdf_ligand(FA, &atoms, &residue, ligand_file);
				} else {
					printf("read ligand MOL2 <%s>\n", ligand_file);
					lig_ok = read_mol2_ligand(FA, &atoms, &residue, ligand_file);
				}

				if (!lig_ok) {
					fprintf(stderr, "ERROR: Failed to read ligand file: %s\n", ligand_file);
					Terminate(2);
				}

				// ── Tier 2: enrich SDF ligand atom types with ProcessLigand ───
				// read_sdf_ligand() above assigned element-only generic types
				// (C.3/N.3/O.3/S.3) because SDF/MOL files carry no hybridisation
				// or aromaticity. ProcessLigand perceives rings, aromaticity and
				// hybridisation, so for the organic C/N/O/S atoms we replace the
				// generic type with the perceived canonical VCT index.
				//
				// Scope (deliberately narrow to avoid regressions):
				//   • SDF only. MOL2 atoms keep the file's native SYBYL types,
				//     which read_mol2_ligand() already maps to canonical VCT
				//     indices and which are authoritative (Tripos standard);
				//     re-perceiving them adds risk with no benefit.
				//   • C/N/O/S only. bonmol::SybylTyper defaults elements it does
				//     not handle (Se, metals, …) to C.3, so overriding those
				//     would corrupt the element-based types the reader assigns
				//     correctly.
				// BonMol atom order matches the reader — both parse the source
				// file sequentially (H included) — so BonMol atom i maps to
				// ligand atom fatm[0]+i. assign_radii_types() (called next) skips
				// ligand residues, so these overrides persist into the GA.
				if (is_sdf && result.success && !result.mol.atoms.empty()) {
					int lig_res = FA->res_cnt;
					int fa     = residue[lig_res].fatm[0];
					int la     = residue[lig_res].latm[0];
					int n_lig  = la - fa + 1;
					int n_bm   = result.mol.num_atoms();

					if (n_bm != n_lig) {
						printf("[TYPING] WARNING: ProcessLigand atom count (%d) != "
						       "reader ligand atom count (%d); cannot map by index — "
						       "keeping fallback types from the reader.\n",
						       n_bm, n_lig);
					} else {
						int upgraded = 0;
						for (int i = 0; i < n_bm; ++i) {
							const bonmol::Atom& bm = result.mol.atoms[i];
							// Only organic C/N/O/S gain information from SYBYL
							// perception; every other element is already canonical
							// from the element-based reader (and BonMol would
							// default unhandled elements such as Se/metals to C.3).
							if (bm.element != bonmol::Element::C &&
							    bm.element != bonmol::Element::N &&
							    bm.element != bonmol::Element::O &&
							    bm.element != bonmol::Element::S)
								continue;
							const char* sname = nullptr;
							int canon  = bonmol_atom_to_canonical_vct(bm, &sname);
							int fa_idx = fa + i;
							int old_t  = atoms[fa_idx].type;
							// read_sdf_ligand() now perceives C.ar/N.ar/O.co2
							// directly from the SDF connection table (MDL bond
							// type 4 = aromatic; COO single+double pattern). Those
							// three come from explicit bond orders and are more
							// reliable than BonMol's ring perception, which can
							// miss aromaticity (reports "0 aromatic") and would
							// otherwise silently downgrade C.ar->C.3 / N.ar->N.am
							// / O.co2->O.3. Preserve the reader's topology type
							// in that case instead of letting BonMol overwrite it.
							const bool reader_perceived_hybrid =
							    (old_t == 4  /*C.ar */ ||
							     old_t == 10 /*N.ar */ ||
							     old_t == 15 /*O.co2*/);
							// Leave the reader's fallback in place for DUMMY (e.g. H),
							// otherwise upgrade to the SYBYL-based canonical type.
							if (canon != FA_TYPE_DUMMY && canon != old_t &&
							    !reader_perceived_hybrid) {
								printf("[TYPING] atom %d (%s): element-only type %d "
								       "-> SYBYL %s -> canonical type %d\n",
								       i, atoms[fa_idx].name, old_t, sname, canon);
								atoms[fa_idx].type = canon;
								++upgraded;
							} else if (reader_perceived_hybrid && canon != old_t &&
							           canon != FA_TYPE_DUMMY) {
								printf("[TYPING] atom %d (%s): keeping reader "
								       "topology type %d over BonMol %s (canonical "
								       "%d) — SDF bond orders are authoritative for "
								       "aromatic/carboxylate perception\n",
								       i, atoms[fa_idx].name, old_t, sname, canon);
							}
						}
						printf("ProcessLigand typing applied: %d/%d ligand atoms "
						       "upgraded to SYBYL-based canonical VCT types\n",
						       upgraded, n_bm);
					}
				} else if (!is_sdf) {
					printf("MOL2 ligand: keeping native SYBYL types from "
					       "read_mol2_ligand() (authoritative)\n");
				} else {
					printf("ProcessLigand typing not applied — using fallback "
					       "types from the reader\n");
				}
			}
		}

		// ── 6. Assign radii and types ──
		assign_radii_types(FA, atoms, residue);
		printf("radii are now assigned\n");

		// Recalculate the global frame with the ligand loaded, matching
		// the dedicated direct-input pipeline.
		calc_center(FA, atoms, residue);

		// ── 6a. Assign formal charges to receptor atoms from PDB ──
		// PDB files carry no charge data — this assigns AMBER ff14SB
		// partial charges to titratable side-chains and formal charges
		// to metal ions, enabling Coulomb electrostatics and salt bridge
		// detection. Skips atoms that already have charges (MOL2/PTM).
		formal_charges::assign_formal_charges(FA, atoms, residue);

		// ── 6c. Populate type256 for ALL atoms (v57: donor/acceptor roles) ───
		// Before this call type256 = 0 for every atom because atom256::encode()
		// was never wired into the main docking binary (it only existed in the
		// standalone ProcessLigand prep tool).  With type256 = 0:
		//   • atom256::get_hbond_donor(0)    → false  (bit 7 = 0)
		//   • atom256::get_hbond_acceptor(0) → false  (bit 6 = 0)
		// Result: hbond_potential.h early-returns 0.0 for every pair →
		//   E_hb = 0 in every pose → cf.hbond = 0 in every PDB output.
		// Fix: encode_from_sybyl() maps the SYBYL type (1–40), partial charge
		// (atoms[i].charge, now populated by assign_formal_charges above), and
		// effective H count into the 8-bit layout [D:1][A:1][base:6] that
		// hbond_potential.h reads. Most benchmark structures are heavy-atom
		// only, so donor roles use explicit bonded H when present and a
		// conservative implicit-H estimate otherwise.
		// Called after assign_formal_charges() so receptor charges are final.
		{
			auto is_hydrogen_atom = [](const atom& a) {
				return a.element[0] == 'H' ||
				       (a.element[0] == ' ' && a.element[1] == 'H') ||
				       a.name[0] == 'H';
			};
			auto bonded_hydrogen_count = [&](int atom_idx) {
				int n_h = 0;
				for (int b = 1; b <= atoms[atom_idx].bond[0] && b <= 6; ++b) {
					int nb = atoms[atom_idx].bond[b];
					if (nb >= 0 && is_hydrogen_atom(atoms[nb])) ++n_h;
				}
				return n_h;
			};
			auto heavy_neighbor_count = [&](int atom_idx) {
				int n_heavy = 0;
				for (int b = 1; b <= atoms[atom_idx].bond[0] && b <= 6; ++b) {
					int nb = atoms[atom_idx].bond[b];
					if (nb >= 0 && !is_hydrogen_atom(atoms[nb])) ++n_heavy;
				}
				return n_heavy;
			};
			auto conservative_implicit_h_count = [&](int atom_idx, int explicit_h, int res_k) {
				if (explicit_h > 0) return 0;
				const atom& a = atoms[atom_idx];
				const int heavy_bonds = heavy_neighbor_count(atom_idx);
				switch (a.type) {
					case 7:  // N.2
						return heavy_bonds <= 1 ? 1 : 0;
					case 8: { // N.3
						const int valence = (a.charge >= 0.3f) ? 4 : 3;
						const int h = valence - heavy_bonds;
						return h > 0 ? h : 0;
					}
								case 11: // N.am — restored: virtual-H (VHG_AMIDE) provides
					         // planar angular discrimination; no longer suppressed.
					         // PRO N is tertiary (no labile H): return 0.
					    if (strcmp(residue[res_k].name, "PRO") == 0) return 0;
					    return heavy_bonds <= 2 ? 1 : 0;
					case 12: { // N.pl3
						const int h = 3 - heavy_bonds;
						return h > 0 ? h : 0;
					}
					case 14: // O.3
						return heavy_bonds <= 1 ? 1 : 0;
					case 18: // S.3
						return heavy_bonds <= 1 ? 1 : 0;
					case 10: {  // N.ar: topology-based pyrrole-donor vs pyridine-acceptor discrimination.
						// 2-hop BFS from the two heavy neighbors detects a closing X–Y bond that
						// completes atom_idx–nbA–X–Y–nbB–atom_idx (5-membered ring).
						// 5-ring → pyrrole/indole/benzimidazole NH (donor, return 1).
						// No 5-ring → 6-ring or larger → pyridine-like (acceptor, return 0).
						// For diazoles: explicit-H on partner (R1) or partner with ≥3 heavy bonds
						// (R2, covers purine N7/N9) disambiguates; charge fallback for MOL2 Gasteiger.
						if (heavy_bonds >= 3) return 0;  // bridgehead/tertiary: no NH
						if (heavy_bonds != 2) return 0;  // degenerate edge case
						int nbA = -1, nbB = -1;
						for (int b = 1; b <= atoms[atom_idx].bond[0] && b <= 6; ++b) {
							int nb = atoms[atom_idx].bond[b];
							if (nb < 0 || is_hydrogen_atom(atoms[nb])) continue;
							if (nbA < 0) nbA = nb; else { nbB = nb; break; }
						}
						if (nbA < 0 || nbB < 0) return (a.charge > -0.15f) ? 1 : 0; // topology unavailable
						// 2-hop BFS: for each X (heavy nbr of nbA) and Y (heavy nbr of nbB),
						// if X bonds Y → closes 5-ring: atom_idx–nbA–X–Y–nbB–atom_idx.
						for (int bA = 1; bA <= atoms[nbA].bond[0] && bA <= 6; ++bA) {
							int X = atoms[nbA].bond[bA];
							if (X < 0 || X == atom_idx || is_hydrogen_atom(atoms[X])) continue;
							for (int bB = 1; bB <= atoms[nbB].bond[0] && bB <= 6; ++bB) {
								int Y = atoms[nbB].bond[bB];
								if (Y < 0 || Y == atom_idx || is_hydrogen_atom(atoms[Y])) continue;
								for (int bXY = 1; bXY <= atoms[X].bond[0] && bXY <= 6; ++bXY) {
									if (atoms[X].bond[bXY] != Y) continue;
									// 5-ring confirmed: atom_idx–nbA–X–Y–nbB–atom_idx.
									// Scan ring atoms {nbA,X,Y,nbB} for another N.ar (diazole check).
									const int ring5[4] = {nbA, X, Y, nbB};
									int other_nar = -1;
									for (int ri = 0; ri < 4; ++ri)
									    if (atoms[ring5[ri]].type == 10) { other_nar = ring5[ri]; break; }
									if (other_nar < 0) return 1;  // sole N.ar in 5-ring → pyrrole-type donor
									// Diazole: 2 N.ar in same 5-ring (imidazole/pyrazole/purine half).
									// R1: partner has explicit H → it IS the pyrrole-NH → we are acceptor.
									if (bonded_hydrogen_count(other_nar) > 0) return 0;
									// R2: partner has ≥3 heavy bonds (e.g. purine N9 with glycosidic bond)
									//     → it occupies the substituted bridgehead → we are also not pyrrole-NH.
									if (heavy_neighbor_count(other_nar) >= 3) return 0;
									// Charge fallback: pyrrole-NH ≈ 0 to −0.10; pyridine-like ≈ −0.20 to −0.45.
									// −0.15 threshold cleanly separates MOL2/Gasteiger. For SDF (charge=0)
									// R1 (partner explicit-H) should have fired; this handles MOL2.
									return (a.charge > -0.15f) ? 1 : 0;
								}
							}
						}
						return 0;  // no 5-ring found → 6-ring or larger → pyridine-like → no NH
					}
					default:
						return 0;
				}
			};
			for (int k = 1; k <= FA->res_cnt; k++) {
				for (int i = residue[k].fatm[0]; i <= residue[k].latm[0]; i++) {
					if (atoms[i].type > 0) {
						const int explicit_h = bonded_hydrogen_count(i);
						const int n_hydrogens = explicit_h +
							conservative_implicit_h_count(i, explicit_h, k);
						// Heavy-atom substitution evidence for amine/alcohol roles.
						// Only trusted when the atom actually carries a bond list:
						// PDB receptor atoms can arrive with bond[0]==0, and a
						// fabricated heavy count of 0 would read as a primary amine.
						// With known=false the classifier reproduces its previous
						// verdict exactly.
						atom256::HbondTopology topo;
						topo.n_heavy_neighbors = heavy_neighbor_count(i);
						topo.known             = (atoms[i].bond[0] > 0);
						atoms[i].type256 = atom256::encode_from_sybyl(
							atoms[i].type,   // SYBYL type 1–40
							atoms[i].charge, // partial charge (MOL2 or AMBER ff14SB)
							n_hydrogens,     // explicit + conservative implicit H
							topo             // heavy-atom substitution evidence
						);
						// Virtual-H geometry recipe: stores heavy-neighbor indices so
						// hbond_potential.h reconstructs H direction from live coords
						// at every scoring call. N.am uses VHG_AMIDE (planar bisector).
						hbond::assign_virtual_h_geometry(atoms, i, explicit_h, heavy_neighbor_count(i),
							strcmp(residue[k].name, "PRO") == 0);
					}
				}
			}
			printf("[vH] type256 + virtual-H populated (N.am=VHG_AMIDE, N.3=SP3, O.3/S.3=HYDROXYL)\n");

			// ── vH assignment diagnostics ─────────────────────────────────────
			// Activated by env FLEXAIDS_VH_DEBUG=1. Dumps per-atom vH recipe and
			// bond[] population. Confirms whether receptor atoms (PDB, bond[]=0)
			// are getting VHG_NONE instead of VHG_AMIDE — the suspected root cause
			// of cf_native(1JD0)=-1.23 (should be ~-23 from backbone amide H-bonds).
			if (std::getenv("FLEXAIDS_VH_DEBUG")) {
				const char* kind_names[] = {"NONE","AMIDE","SP2_1NBR","SP2_2NBR",
				                            "SP3_2NBR","SP3_1NBR","HYDROXYL"};
				int n_am_amide=0, n_am_none=0, n_active=0, n_none_donor=0;
				for (int k=1; k<=FA->res_cnt; k++) {
					for (int i=residue[k].fatm[0]; i<=residue[k].latm[0]; i++) {
						if (atoms[i].type <= 0) continue;
						const bool is_donor = atom256::get_hbond_donor(atoms[i].type256);
						const int bond_cnt = atoms[i].bond[0];
						const uint8_t kind  = atoms[i].vH_kind;
						// Classify N.am specifically
						if (atoms[i].type == 11) {
							if (kind == hbond::VHG_AMIDE) ++n_am_amide; else ++n_am_none;
						}
						if (is_donor) {
							if (kind != hbond::VHG_NONE) ++n_active; else ++n_none_donor;
						}
						// Print every donor atom's assignment
						if (is_donor) {
							printf("[vHdbg] atom[%4d] res%-4d %-5s type=%-2d bond_cnt=%-2d "
							       "vH_kind=%-10s vH_n=%d nbr=[%d,%d]%s\n",
							       i, k, atoms[i].name, atoms[i].type, bond_cnt,
							       (kind<7?kind_names[kind]:"?"), atoms[i].vH_n,
							       atoms[i].vH_nbr[0], atoms[i].vH_nbr[1],
							       (bond_cnt==0 && atoms[i].type==11)
							           ? "  <-- N.am NO BONDS (PDB receptor?)" : "");
						}
					}
				}
				printf("[vHdbg] ── Summary ──────────────────────────────────────────\n");
				printf("[vHdbg]   N.am with VHG_AMIDE : %d\n", n_am_amide);
				printf("[vHdbg]   N.am with VHG_NONE  : %d  (<-- if >0, receptor bond[] empty)\n", n_am_none);
				printf("[vHdbg]   Donors with active vH: %d\n", n_active);
				printf("[vHdbg]   Donors with VHG_NONE : %d  (use 0.3 fallback)\n", n_none_donor);
			}
		}

		// ── 6b. Set up GPA and IC origin for MOL2/SDF ligand ──
		// generate_grid() requires residue[last].gpa to be non-NULL and
		// atoms[gpa[0]].dis/ang/dih to be computed (normally done by
		// read_lig for legacy .inp/.ic format). For direct-mode ligands,
			// Direct readers normally provide a bonded, topology-derived GPA
			// triad. The first-three fallback is retained only for legacy readers.
		{
			int lig_res = FA->res_cnt;
			if (residue[lig_res].gpa == NULL) {
				int fa = residue[lig_res].fatm[0];
				int la = residue[lig_res].latm[0];
				int n_lig = la - fa + 1;

				// FA->ori was set to the receptor centroid by calc_center() above.
				// Do NOT overwrite with ligand centroid here: the cleftgrid IC
				// (calc_cleftic) are encoded relative to this receptor-center ori,
				// and buildcc uses ori as the GPA1/GPA2 grandparent reference.
				// Overwriting with the ligand centroid breaks that reference frame
				// when gene[0] translates GPA0 far from the ligand starting position.
				printf("the receptor center of coordinates is: %8.3f %8.3f %8.3f\n",
				       FA->ori[0], FA->ori[1], FA->ori[2]);

				// Allocate gpa (3 global-positioning atoms)
				residue[lig_res].gpa = (int*)malloc(3 * sizeof(int));
				if (!residue[lig_res].gpa) {
					fprintf(stderr, "ERROR: malloc for residue.gpa\n");
					Terminate(2);
				}
				residue[lig_res].gpa[0] = fa;
				residue[lig_res].gpa[1] = (n_lig > 1) ? fa + 1 : fa;
				residue[lig_res].gpa[2] = (n_lig > 2) ? fa + 2 : fa;

				// Compute IC for GPA atom relative to FA->ori
				buildic_point(FA, atoms[fa].coor,
				              &atoms[fa].dis, &atoms[fa].ang, &atoms[fa].dih);
			}
		}

		// ── 7. Binding site detection (oracle LOCCLF or SURFNET AUTO) ──
		{
			strcpy(FA->rngopt, "locclf");

			// ── Strategy A: lazy grid cache load ──────────────────────────────────
			// When DatasetRunner passes "grid_file" in the JSON config (pointing to
			// <out_prefix>.rrg from a prior same-receptor run), load the fully
			// pruned + IC-transformed cleftgrid directly, bypassing SURFNET +
			// generate_grid + site-confine + MIF.  The gridpoint struct stores
			// pre-computed IC values (dis/ang/dih) relative to FA->ori, which is
			// set during receptor loading (before this block) and is identical for
			// the same receptor PDB.  No calc_cleftic() call needed after load.
			bool grid_loaded_from_cache = false;
			static constexpr uint32_t RRG_MAGIC   = 0x56435400U;
			static constexpr uint32_t RRG_VERSION = 1U;

			// ── Multi-cleft restoration hook ────────────────────────────────────
			// FLEXAIDDS_CLEFT_SPHERE_FILE is an opt-in direct-mode override for
			// process-per-cleft docking.  Each child process receives one ranked
			// Get_Cleft/FlexAID sphere file and builds its grid from that cleft
			// only, reproducing the old independent-GA-per-major-cleft behavior.
			const flexaids::ProtocolConfig proto_grid = flexaids::ProtocolConfig::from_env();
			const char* explicit_cleft_env =
			    proto_grid.cleft_sphere_file.empty() ? nullptr
			                                        : proto_grid.cleft_sphere_file.c_str();
			const bool explicit_cleft_requested =
			    explicit_cleft_env && explicit_cleft_env[0] != '\0';
			if (explicit_cleft_requested &&
			    !std::filesystem::exists(explicit_cleft_env)) {
				fprintf(stderr,
				        "ERROR: FLEXAIDDS_CLEFT_SPHERE_FILE does not exist: %s\n",
				        explicit_cleft_env);
				Terminate(2);
			}

			// ── Probe FLEXAIDDS_GRID_CACHE_DIR (content-hash keyed) ─────────────
			// Used by Strategy B batch children: the first child builds and saves;
			// subsequent children of the same receptor load and skip the grid build.
			// GRID_CACHE_DIR stays raw getenv (infra path, not a science knob).
			const char* oracle_site_for_cache =
			    proto_grid.oracle_site.empty() ? nullptr : proto_grid.oracle_site.c_str();
			const std::string gc_vct_path = explicit_cleft_requested
			    ? std::string()
			    : gridcache::cache_path(receptor_file, oracle_site_for_cache,
			                            FA->spacer_length, FA->permeability);
			if (!explicit_cleft_requested && !gc_vct_path.empty()) {
				if (gridcache::load(gc_vct_path, &cleftgrid, &FA->num_grd)) {
					grid_loaded_from_cache = true;
					fprintf(stderr,
					        "[GRID-CACHE] Loaded %d pts from GRID_CACHE_DIR cache: %s\n",
					        FA->num_grd, gc_vct_path.c_str());
				}
			}

			if (!explicit_cleft_requested && !cached_grid_path.empty() &&
			    std::filesystem::exists(cached_grid_path)) {
				FILE* gf = fopen(cached_grid_path.c_str(), "rb");
				if (gf) {
					uint32_t magic = 0, version = 0;
					int32_t  ng    = 0;
					if (fread(&magic,   sizeof magic,   1, gf) == 1 &&
					    fread(&version, sizeof version, 1, gf) == 1 &&
					    fread(&ng,      sizeof ng,      1, gf) == 1 &&
					    magic == RRG_MAGIC && version == RRG_VERSION && ng > 0) {
						gridpoint* cached_gp = static_cast<gridpoint*>(
						    malloc(static_cast<size_t>(ng) * sizeof(gridpoint)));
						if (cached_gp &&
						    static_cast<int32_t>(
						        fread(cached_gp, sizeof(gridpoint), static_cast<size_t>(ng), gf))
						        == ng) {
							cleftgrid           = cached_gp;
							FA->num_grd         = static_cast<int>(ng);
							grid_loaded_from_cache = true;
							fprintf(stderr, "[GRID-CACHE] Loaded %d grid points from: %s\n",
							        ng, cached_grid_path.c_str());
						} else {
							free(cached_gp);
							fprintf(stderr,
							        "[GRID-CACHE] WARN: incomplete read from %s — regenerating\n",
							        cached_grid_path.c_str());
						}
					} else {
						fprintf(stderr,
						        "[GRID-CACHE] WARN: bad header in %s (magic=%08X ver=%u) — regenerating\n",
						        cached_grid_path.c_str(), magic, version);
					}
					fclose(gf);
				} else {
					fprintf(stderr,
					        "[GRID-CACHE] WARN: cannot open %s — regenerating\n",
					        cached_grid_path.c_str());
				}
			}

			auto initialize_direct_mif = [&](bool allow_grid_prune) {
				if (!(FA->mif_enabled || FA->grid_prio_percent < 100.0f)) return;

				std::vector<atom> protein_atoms(atoms, atoms + FA->atm_cnt_real);
				cavity_detect::SpatialGrid sg;
				sg.build(protein_atoms);
				auto field = mif::compute_mif(cleftgrid, FA->num_grd,
				                              atoms, FA->atm_cnt_real, sg);

				free(FA->mif_energies);
				free(FA->mif_sorted);
				free(FA->mif_cdf);
				FA->mif_count = static_cast<int>(field.sorted_indices.size());
				FA->mif_energies = static_cast<float*>(
				    malloc(field.energies.size() * sizeof(float)));
				FA->mif_sorted = static_cast<int*>(
				    malloc(field.sorted_indices.size() * sizeof(int)));
				std::copy_n(field.energies.data(), field.energies.size(),
				            FA->mif_energies);
				std::copy_n(field.sorted_indices.data(), field.sorted_indices.size(),
				            FA->mif_sorted);

				mif::build_sampling_cdf(field, FA->mif_temperature);
				FA->mif_cdf = static_cast<double*>(
				    malloc(field.cdf.size() * sizeof(double)));
				std::copy_n(field.cdf.data(), field.cdf.size(), FA->mif_cdf);

				if (allow_grid_prune && FA->grid_prio_percent < 100.0f) {
					auto kept = mif::prioritize_grid(field, FA->grid_prio_percent);
					gridpoint* new_grid = nullptr;
					int new_count = mif::rebuild_cleftgrid(
					    cleftgrid, FA->num_grd, kept, &new_grid);
					if (new_grid && new_count > 0) {
						int old_count = FA->num_grd;
						free(cleftgrid);
						cleftgrid = new_grid;
						FA->num_grd = new_count;
						calc_cleftic(FA, cleftgrid);
						printf("GRIDPRIO: kept %d/%d grid points (top %.0f%%)\n",
						       new_count - 1, old_count - 1,
						       FA->grid_prio_percent);
					}
				}

				printf("MIF: computed for %d grid points (T=%.0fK%s)\n",
				       FA->mif_count, FA->mif_temperature,
				       grid_loaded_from_cache ? ", cached grid" : "");
			};

			if (!grid_loaded_from_cache) {
			// Oracle mode: FLEXAIDDS_ORACLE_SITE (ProtocolConfig) points to a binding
			// site PDB.  Parse for centroid only — SURFNET void-space detection
			// always runs (probes placed in void between atoms, not on atoms).
			// Oracle centroid guides cleft selection and site-confinement.
			const char* oracle_site_env =
			    proto_grid.oracle_site.empty() ? nullptr : proto_grid.oracle_site.c_str();
			bool using_oracle = false;
			float oracle_cx = 0.0f, oracle_cy = 0.0f, oracle_cz = 0.0f;
			sphere* spheres = NULL;
			bool using_explicit_cleft = false;

			if (explicit_cleft_requested) {
				std::vector<char> cleft_file(explicit_cleft_env,
				                             explicit_cleft_env + strlen(explicit_cleft_env) + 1);
				spheres = read_spheres(cleft_file.data());
				using_explicit_cleft = true;
				fprintf(stderr, "[MULTI-CLEFT] using explicit cleft sphere file: %s\n",
				        explicit_cleft_env);
			} else if (oracle_site_env && oracle_site_env[0] != '\0' &&
			           std::filesystem::exists(oracle_site_env)) {
				int n = 0;
				FILE* fp_oracle = fopen(oracle_site_env, "r");
				if (fp_oracle) {
					char line[256];
					while (fgets(line, sizeof(line), fp_oracle)) {
						if (strncmp(line, "ATOM", 4) != 0 && strncmp(line, "HETATM", 6) != 0) continue;
						if ((int)strlen(line) < 54) continue;
						float x, y, z;
						if (sscanf(line + 30, "%8f%8f%8f", &x, &y, &z) != 3) continue;
						oracle_cx += x; oracle_cy += y; oracle_cz += z; ++n;
					}
					fclose(fp_oracle);
				}
				if (n > 0) {
					oracle_cx /= n; oracle_cy /= n; oracle_cz /= n;
					using_oracle = true;
					fprintf(stderr, "[ORACLE] centroid: %.2f %.2f %.2f (%d atoms)\n",
					        oracle_cx, oracle_cy, oracle_cz, n);
					fprintf(stderr, "[ORACLE] using SURFNET void-space + oracle centroid guidance\n");
				} else {
					fprintf(stderr, "[WARN] Oracle site load failed, falling back to AUTO\n");
				}
			}

			// --redock: the extracted cognate ligand is already loaded at crystal
			// coordinates. Use that centroid as the oracle so SURFNET + SITE-CONFINE
			// target the known pocket (FLEXAIDDS_ORACLE_SITE not required).
			if (!using_oracle && !using_explicit_cleft && !redock_pdb_id.empty()) {
				const int lig_res = FA->res_cnt;
				if (lig_res >= 1 && residue[lig_res].fatm && residue[lig_res].latm) {
					const int fa = residue[lig_res].fatm[0];
					const int la = residue[lig_res].latm[0];
					if (la >= fa) {
						double sx = 0, sy = 0, sz = 0;
						int nn = 0;
						for (int a = fa; a <= la; ++a) {
							sx += atoms[a].coor[0];
							sy += atoms[a].coor[1];
							sz += atoms[a].coor[2];
							++nn;
						}
						if (nn > 0) {
							oracle_cx = static_cast<float>(sx / nn);
							oracle_cy = static_cast<float>(sy / nn);
							oracle_cz = static_cast<float>(sz / nn);
							using_oracle = true;
							printf("[REDOCK] oracle centroid from cognate ligand: "
							       "%.2f %.2f %.2f (%d atoms)\n",
							       oracle_cx, oracle_cy, oracle_cz, nn);
						}
					}
				}
			}

			// Always run SURFNET void-space detection (probes placed in void between atoms)
			// unless a child process was given one explicit cleft sphere file.
			// In oracle mode, spatially pre-filter atoms to the oracle centroid sphere so
			// that multimeric receptors (e.g. 1OF6 octamer: 20826 atoms, 8 chains) do not
			// cause O(N^3) blowup.  15 A covers the binding pocket + one shell of framework
			// atoms needed to correctly cap probe spheres at the cavity boundary.
			if (!using_explicit_cleft) {
				printf("SURFNET binding-site detection (CleftDetector) ...\n");
				CleftDetectorParams cleft_params = default_cleft_params();
				if (using_oracle) {
					cleft_params.oracle_center[0] = oracle_cx;
					cleft_params.oracle_center[1] = oracle_cy;
					cleft_params.oracle_center[2] = oracle_cz;
					cleft_params.oracle_radius     = 15.0f;  // Å — covers pocket + 1 shell
				}
				spheres = detect_cleft(atoms, residue, FA->atm_cnt_real, FA->res_cnt,
				                       cleft_params);
			}
			if (spheres == NULL) {
				fprintf(stderr, "ERROR: cleft detection found no cavities.\n");
				Terminate(2);
			}

			cleftgrid = generate_grid(FA, spheres, atoms, residue);
			calc_cleftic(FA, cleftgrid);

			// Layer 2 reproducibility: capture cleft centroid/extent BEFORE free
			// so explicit multi-cleft runs confine Ω to the sphere support
			// (not whole-protein translation).
			ensemble::CleftCentroid cleft_geom{};
			const bool have_cleft_geom =
				ensemble::cleft_centroid_extent(spheres, &cleft_geom);

			// Free sphere linked list (oracle or SURFNET)
			while (spheres) { sphere* p = spheres->prev; free(spheres); spheres = p; }

			// ── Confine search to the cognate (reference-ligand) site ──────
			// Skipped in oracle mode — the binding site PDB already defines
			// the precise search space without further confinement needed.
			// Re-docking benchmark: the binding site is known.  Restrict the
			// auto-detected cleftgrid to grid points within (ligand_radius +
			// margin) of the cognate centroid; otherwise the GA can settle in a
			// wrong cavity (1IGJ 74 Å, 1GM8 32 Å, 1GPK 6.8 Å off-centre).
			// Index 0 (reflig reference conformation) is always preserved by
			// mif::rebuild_cleftgrid.
			if (!using_explicit_cleft) {
				const float kSiteMargin = 0.0f;            // Å beyond ligand extent
				int lig_res = FA->res_cnt;
				int fa = residue[lig_res].fatm[0];
				int la = residue[lig_res].latm[0];
				if (la >= fa && FA->num_grd > 1) {
					double cx=0, cy=0, cz=0; int nn=0;
					for (int a = fa; a <= la; ++a) {
						cx += atoms[a].coor[0];
						cy += atoms[a].coor[1];
						cz += atoms[a].coor[2]; ++nn;
					}
					cx/=nn; cy/=nn; cz/=nn;
					bool using_explicit_reflig_site = false;
					std::vector<reflig::RefLigAtom> site_atoms;
					if (!using_oracle && FA->reflig_file[0] != '\0' &&
					    std::filesystem::exists(FA->reflig_file)) {
						site_atoms = reflig::parse_reflig(FA->reflig_file);
						if (!site_atoms.empty()) {
							float site_centroid[3];
							reflig::compute_centroid(site_atoms, site_centroid);
							cx = site_centroid[0];
							cy = site_centroid[1];
							cz = site_centroid[2];
							using_explicit_reflig_site = true;
							printf("SITE-CONFINE: reference ligand centroid override %.2f %.2f %.2f from %s\n",
							       cx, cy, cz, FA->reflig_file);
						} else {
							printf("SITE-CONFINE: WARNING reference ligand file had no atoms: %s\n",
							       FA->reflig_file);
						}
					}
					// Oracle mode: override centroid with oracle-derived position
					// (rmax2 still computed from ligand extent for a sensible rcut_initial)
					if (using_oracle) {
						cx = oracle_cx; cy = oracle_cy; cz = oracle_cz;
						printf("SITE-CONFINE: oracle centroid override %.2f %.2f %.2f\n", cx, cy, cz);
					}
					double rmax2 = 0.0;
					if (using_explicit_reflig_site) {
						for (const auto& ra : site_atoms) {
							double dx=ra.x-cx, dy=ra.y-cy, dz=ra.z-cz;
							double d2 = dx*dx+dy*dy+dz*dz;
							if (d2 > rmax2) rmax2 = d2;
						}
					} else {
						for (int a = fa; a <= la; ++a) {
							double dx=atoms[a].coor[0]-cx, dy=atoms[a].coor[1]-cy, dz=atoms[a].coor[2]-cz;
							double d2 = dx*dx+dy*dy+dz*dz;
							if (d2 > rmax2) rmax2 = d2;
						}
					}
					// Expanding-radius confinement: always confine to the cognate
					// site.  Start at (ligand_extent + margin) and grow rcut in 2 Å
					// steps (up to 30 Å) until at least MIN_SITE_GRID points fall
					// inside the sphere.  The previous behaviour kept a single static
					// radius and, when it yielded < MIN_SITE_GRID points, silently
					// fell back to the full ~123k-point grid — letting the GA roam the
					// entire protein surface.
					const int    MIN_SITE_GRID  = 500;   // floor: target site density (engine min ~250, 2x for safety)
					const double rcut_initial   = std::sqrt(rmax2) + kSiteMargin;
					const double rcut_max       = 8.0;   // Å — hard ceiling on expansion (v13: oracle-style ~6 A pocket + 2 A margin)
					const double rcut_step       = 2.0;  // Å — expansion increment
					double rcut = rcut_initial;
					std::vector<int> keep;
					keep.reserve(FA->num_grd);
					for (;;) {
						keep.clear();
						const double rcut2 = rcut*rcut;
						for (int i = 1; i < FA->num_grd; ++i) {  // i=0 = reflig ref conf
							double dx=cleftgrid[i].coor[0]-cx, dy=cleftgrid[i].coor[1]-cy, dz=cleftgrid[i].coor[2]-cz;
							if (dx*dx+dy*dy+dz*dz <= rcut2) keep.push_back(i);
						}
						if ((int)keep.size() >= MIN_SITE_GRID || rcut >= rcut_max) break;
						rcut += rcut_step;
					}
					// Crash guard: rcut_max=8 Å may leave < MIN_SITE_GRID pts for sparse pockets.
					// Retry up to 12 Å before falling through to full-grid fallback.
					if ((int)keep.size() < MIN_SITE_GRID && rcut >= rcut_max) {
						const double rcut_extended = 12.0;
						while (rcut < rcut_extended) {
							rcut += rcut_step;
							keep.clear();
							const double rcut2 = rcut * rcut;
							for (int i = 1; i < FA->num_grd; ++i) {
								double dx = cleftgrid[i].coor[0]-cx,
								       dy = cleftgrid[i].coor[1]-cy,
								       dz = cleftgrid[i].coor[2]-cz;
								if (dx*dx+dy*dy+dz*dz <= rcut2) keep.push_back(i);
							}
							if ((int)keep.size() >= MIN_SITE_GRID) break;
						}
						printf("SITE-CONFINE: sparse-pocket retry to %.1f A -> %d pts\n",
						       rcut, (int)keep.size());
					}
					if (flexaids::site_confine_should_rebuild((int)keep.size(), FA->num_grd - 1)) {
						gridpoint* confined = nullptr;
						int new_count = mif::rebuild_cleftgrid(cleftgrid, FA->num_grd, keep, &confined);
						if (confined && new_count > 0) {
							int old_count = FA->num_grd;
							free(cleftgrid);
							cleftgrid = confined;
							FA->num_grd = new_count;
							calc_cleftic(FA, cleftgrid);
							printf("SITE-CONFINE: %d pts within %.1f A of cognate centroid (expanded from %.1f A, %d->%d grid pts)\n",
							       new_count - 1, rcut, rcut_initial, old_count - 1, new_count - 1);
							if ((int)keep.size() < MIN_SITE_GRID) {
								printf("SITE-CONFINE: WARNING confined to only %d pts (< MIN_SITE_GRID=%d) after expanding to %.1f A — keeping pocket grid (not full-grid fallback)\n",
								       new_count - 1, MIN_SITE_GRID, rcut);
							}
						}
					} else {
						printf("SITE-CONFINE: WARNING full-grid fallback after expanding to %.1f A (keep=%d, total=%d)\n",
						       rcut, (int)keep.size(), FA->num_grd - 1);
					}
				}
			} else if (have_cleft_geom && FA->num_grd > 1) {
				// Explicit multi-cleft: confine translation to sphere support
				// (reproducible Ω_cleft; not whole-protein).
				const double cx = cleft_geom.cx, cy = cleft_geom.cy, cz = cleft_geom.cz;
				const int    MIN_SITE_GRID = 500;
				const double rcut_initial  = std::max(cleft_geom.extent_A, 4.0);
				const double rcut_max      = 12.0;
				const double rcut_step     = 2.0;
				double rcut = rcut_initial;
				std::vector<int> keep;
				keep.reserve(static_cast<size_t>(FA->num_grd));
				for (;;) {
					keep.clear();
					const double rcut2 = rcut * rcut;
					for (int i = 1; i < FA->num_grd; ++i) {
						const double dx = cleftgrid[i].coor[0] - cx;
						const double dy = cleftgrid[i].coor[1] - cy;
						const double dz = cleftgrid[i].coor[2] - cz;
						if (dx * dx + dy * dy + dz * dz <= rcut2) keep.push_back(i);
					}
					if ((int)keep.size() >= MIN_SITE_GRID || rcut >= rcut_max) break;
					rcut += rcut_step;
				}
				printf("SITE-CONFINE: cleft-centroid %.2f %.2f %.2f extent=%.1f A "
				       "n_spheres=%d rcut=%.1f keep=%d/%d\n",
				       cx, cy, cz, cleft_geom.extent_A, cleft_geom.n_spheres,
				       rcut, (int)keep.size(), FA->num_grd - 1);
				if (flexaids::site_confine_should_rebuild((int)keep.size(), FA->num_grd - 1)) {
					gridpoint* confined = nullptr;
					int new_count = mif::rebuild_cleftgrid(
						cleftgrid, FA->num_grd, keep, &confined);
					if (confined && new_count > 0) {
						const int old_count = FA->num_grd;
						free(cleftgrid);
						cleftgrid = confined;
						FA->num_grd = new_count;
						calc_cleftic(FA, cleftgrid);
						printf("SITE-CONFINE: %d pts within %.1f A of cleft centroid "
						       "(%d->%d grid pts)\n",
						       new_count - 1, rcut, old_count - 1, new_count - 1);
					}
				} else {
					printf("SITE-CONFINE: WARNING explicit-cleft full-grid fallback "
					       "(keep=%d, total=%d)\n",
					       (int)keep.size(), FA->num_grd - 1);
				}
			} else {
				printf("SITE-CONFINE: skipped for explicit multi-cleft (no sphere geometry)\n");
			}

			// Cavity-only MIF is a legitimate known-site prior; it contains no
			// reference-ligand coordinates.
			initialize_direct_mif(true);

			// Free spheres linked list
			while (spheres != NULL) {
				sphere* prev = spheres->prev;
				free(spheres);
				spheres = prev;
			}

			// ── Strategy A: save newly generated grid ────────────────────────────
			// 1. To FLEXAIDDS_GRID_CACHE_DIR/<hash>.vct (atomic temp+rename)
			//    so Strategy B batch children of the same receptor find it.
			// 2. To <output_prefix>.rrg so DatasetRunner's "grid_file" mechanism
			//    can pass it to subsequent same-receptor entries.
			if (!gc_vct_path.empty())
				gridcache::save(gc_vct_path, cleftgrid, FA->num_grd);

			if (cleftgrid && FA->num_grd > 0 && !output_prefix.empty()) {
				std::string rrg_path = output_prefix + ".rrg";
				FILE* gf_save = fopen(rrg_path.c_str(), "wb");
				if (gf_save) {
					int32_t ng_save = static_cast<int32_t>(FA->num_grd);
					fwrite(&RRG_MAGIC,   sizeof RRG_MAGIC,   1, gf_save);
					fwrite(&RRG_VERSION, sizeof RRG_VERSION, 1, gf_save);
					fwrite(&ng_save,     sizeof ng_save,     1, gf_save);
					size_t written = fwrite(cleftgrid, sizeof(gridpoint),
					                        static_cast<size_t>(ng_save), gf_save);
					fclose(gf_save);
					if (static_cast<int32_t>(written) == ng_save) {
						fprintf(stderr, "[GRID-CACHE] Saved %d grid points → %s\n",
						        ng_save, rrg_path.c_str());
					} else {
						fprintf(stderr,
						        "[GRID-CACHE] WARN: partial write (%zu/%d pts) — removing %s\n",
						        written, ng_save, rrg_path.c_str());
						std::filesystem::remove(rrg_path);
					}
				} else {
					fprintf(stderr, "[GRID-CACHE] WARN: cannot write %s\n", rrg_path.c_str());
				}
			}
			} // end if (!grid_loaded_from_cache)
			if (grid_loaded_from_cache) {
				// The grid cache stores geometry, not MIF sampling arrays.
				initialize_direct_mif(false);
			}
		}

			// Direct-mode GA seeding is optional.  A loaded ligand must not become
			// an implicit spatial reference when HETATM fallback is disabled.
			if (FA->reflig_file[0] != '\0' || FA->reflig_hetatm_fallback) {
			int lig_res = FA->res_cnt;
			int fa = residue[lig_res].fatm[0];
			int la = residue[lig_res].latm[0];
			int n_lig = la - fa + 1;
			int anchor = (residue[lig_res].gpa != NULL)
			           ? residue[lig_res].gpa[0]
			           : fa;

			if (n_lig > 0 && cleftgrid && FA->num_grd > 1) {
				float seed[3] = {
					atoms[anchor].coor[0],
					atoms[anchor].coor[1],
					atoms[anchor].coor[2]
				};
				bool seed_from_explicit_reflig = false;
				if (FA->reflig_file[0] != '\0' &&
				    std::filesystem::exists(FA->reflig_file)) {
					auto reflig_atoms = reflig::parse_reflig(FA->reflig_file);
					if (!reflig_atoms.empty()) {
						reflig::compute_centroid(reflig_atoms, seed);
						seed_from_explicit_reflig = true;
					}
				}

				auto nearest = reflig::find_nearest_grid_points(
				    seed, cleftgrid, FA->num_grd, FA->reflig_k_nearest);
				free(FA->reflig_nearest_grid);
				FA->reflig_nearest_count = static_cast<int>(nearest.size());
				FA->reflig_nearest_grid = nullptr;
				if (!nearest.empty()) {
					FA->reflig_nearest_grid = static_cast<int*>(
					    malloc(nearest.size() * sizeof(int)));
					if (FA->reflig_nearest_grid) {
						std::copy_n(nearest.data(), nearest.size(),
						            FA->reflig_nearest_grid);
						printf("REFLIG: direct-mode %s seed (%.1f, %.1f, %.1f), %d nearest points\n",
						       seed_from_explicit_reflig ? "reference" : "loaded-ligand",
						       seed[0], seed[1], seed[2],
						       FA->reflig_nearest_count);
					} else {
						FA->reflig_nearest_count = 0;
						fprintf(stderr, "WARNING: direct-mode reflig seed allocation failed.\n");
					}
				}
			}
		}

		printf("Direct loading: receptor/ligand structures loaded, cleft detected\n");
	}

	//printf("END FILE:<%s>\n",end_strfile);
	//PAUSE;

	/*
	  if(IS_BIG_ENDIAN())
	  printf("platform is big-endian\n");
	  else
	  printf("platform is little-endian\n");    
	*/

	wif083(FA); // initialization of FA->sphere[]
	
	///////////////////////////////////////////////////////////////////////////////
	// memory allocations for param structures
  
	//printf("memory allocation for opt_par\n");

	FA->map_par = (optmap*)malloc(FA->MIN_PAR*sizeof(optmap));
	FA->opt_par = (double*)malloc(FA->MIN_PAR*sizeof(double));
	FA->del_opt_par = (double*)malloc(FA->MIN_PAR*sizeof(double));
	FA->min_opt_par = (double*)malloc(FA->MIN_PAR*sizeof(double));
	FA->max_opt_par = (double*)malloc(FA->MIN_PAR*sizeof(double));
	FA->map_opt_par = (int*)malloc(FA->MIN_PAR*sizeof(int));

	if(!FA->map_par || !FA->opt_par ||
	   !FA->del_opt_par || !FA->min_opt_par || 
	   !FA->max_opt_par || !FA->map_opt_par)
	{
		fprintf(stderr,"ERROR: memory allocation error for opt_par\n");
		Terminate(2);
	}

	memset(FA->map_par,0,FA->MIN_PAR*sizeof(optmap));
	memset(FA->opt_par,0,FA->MIN_PAR*sizeof(double));
	memset(FA->del_opt_par,0,FA->MIN_PAR*sizeof(double));
	memset(FA->min_opt_par,0,FA->MIN_PAR*sizeof(double));
	memset(FA->max_opt_par,0,FA->MIN_PAR*sizeof(double));

	FA->map_par_flexbond_first_index = -1;
	FA->map_par_flexbond_first = NULL;
	FA->map_par_flexbond_last = NULL;
	
	FA->map_par_sidechain_first_index = -1;
	FA->map_par_sidechain_first = NULL;
	FA->map_par_sidechain_last = NULL;
	
	/////////////////////////////////////////////////////////////////////////////////

	if (legacy_mode) {
		printf("Reading input (%s)...\n",dockinp);
		read_input(FA,&atoms,&residue,&rotamer,&cleftgrid,dockinp);
	} else {
		// Direct mode: set up IC bounds and optimization parameters
		// (receptor, ligand, and cleft grid were already loaded above)
		ic_bounds(FA, FA->rngopt);
		if (FA->reflig_nearest_count > 0) {
			// In direct redocking, preserve the native ligand anchor.
			FA->index_min = 0.0;
		}

		int opt[2];
		char chain = ' ';

		// Translation: grid-index gene (typ=-1), picks anchor point from cleft grid
		opt[0] = FA->resligand->number;
		opt[1] = -1;
		add2_optimiz_vec(FA, atoms, residue, opt, chain, "");

		// Rotation: 3 Euler-angle genes (ang + dih + dih of GPA atoms)
		opt[1] = 0;
		add2_optimiz_vec(FA, atoms, residue, opt, chain, "");

		// Ligand torsional flexibility: one dihedral gene per perceived
		// rotatable bond (resligand->fdih, populated by the SDF/MOL2 reader
		// after ring/terminal-bond perception). Legacy mode enumerates these
		// from explicit OPTIMIZ config lines; direct mode previously skipped
		// them entirely, so every ligand docked as a rigid body (npar=4) and
		// could not reach its bound conformer. Gated on FA->intramolecular so
		// the torsional DoF and the internal-clash scoring it requires are
		// always switched on together — a flexible chromosome without
		// intramolecular scoring lets the ligand fold through itself. A
		// genuinely rigid ligand (fdih==0) makes this loop a no-op, preserving
		// the 4-gene rigid search; force_rigid (config intramolecular=false)
		// pins that behaviour for ablation regardless of fdih.
		if (FA->intramolecular && FA->resligand != NULL) {
			for (int b = 1; b <= FA->resligand->fdih; ++b) {
				opt[1] = b;
				add2_optimiz_vec(FA, atoms, residue, opt, chain, "");
			}
		}

		// Side-chain and normal-mode extensions
		add2_optimiz_vec(FA, atoms, residue, opt, chain, "SC");
		add2_optimiz_vec(FA, atoms, residue, opt, chain, "NM");

		if (FA->translational && FA->num_grd == 1) {
			fprintf(stderr, "%sERROR:%s the binding-site has no anchor points\n", tui::err::failtext(), tui::err::reset());
			Terminate(2);
		}

		update_optres(atoms, residue, FA->atm_cnt, FA->optres, FA->num_optres);

		printf("Direct loading complete: %d atoms, %d residues, %d grid points, %d params\n",
		       FA->atm_cnt, FA->res_cnt, FA->num_grd, FA->npar);
	}

	// memory allocation and initialization of VC struct
	if (strcmp(FA->complf,"VCT")==0)
	{
		VC->planedef = FA->vcontacts_planedef;
		
		// Vcontacts memory allocations...
		// ca_rec can be reallocated
		VC->Calc = (atomsas*)malloc(FA->atm_cnt_real*sizeof(atomsas));
		VC->Calclist = (int*)malloc(FA->atm_cnt_real*sizeof(int));
		VC->ca_index = (int*)malloc(FA->atm_cnt_real*sizeof(int));
		VC->seed = (int*)malloc(3*FA->atm_cnt_real*sizeof(int));
		VC->contlist = (contactlist*)malloc(10000*sizeof(contactlist));
    
		// initialize contact atom index
		VC->ca_recsize = 5*FA->atm_cnt_real;
		VC->ca_rec = (ca_struct*)malloc(VC->ca_recsize*sizeof(ca_struct));
		
		if(!VC->ca_rec) {
			fprintf(stderr,"ERROR: memory allocation error for ca_rec\n"); 
			Terminate(2);
		}
		
		if((!VC->Calc) || (!VC->ca_index) || 
		   (!VC->seed) || (!VC->contlist) || (!VC->Calclist)) {
			fprintf(stderr, "ERROR: memory allocation error for (Calc or Calclist or ca_index or seed or contlist)\n");
			Terminate(2);
		}

		for(i=0;i<FA->atm_cnt_real;i++){
			VC->Calc[i].atom = NULL;
			VC->Calc[i].residue = NULL;
			VC->Calc[i].exposed = true;
		}

		if(FA->omit_buried){
			printf("calcuting SAS of non-scorable atoms...\n");
			Vcontacts(FA,atoms,residue,VC,NULL,true);
			
			//FILE* surffile = fopen("surfpdb.pdb", "w");

			int n_buried = 0;
			for(i=0;i<FA->atm_cnt_real;i++){
				if(!VC->Calc[i].score){
					double radoA = VC->Calc[i].atom->radius + Rw;
					double SAS = 4.0*PI*radoA*radoA;
			
					int currindex = VC->ca_index[i];
					while(currindex != -1) {
						double area = VC->ca_rec[currindex].area;
						SAS -= area;
						currindex = VC->ca_rec[currindex].prev;
					}

					if(SAS <= 0.0){
						VC->Calc[i].exposed = false;
						n_buried++;
					}

					/*
					//ATOM    135  CG2 ILE A  30      26.592   6.245  -4.544  1.00 21.36           3
					fprintf(surffile, "ATOM  %5d  XX  XXX A%4d    %8.3f%8.3f%8.3f  1.00  1.00           %2s\n",
							VC->Calc[i].atom->number,VC->Calc[i].residue->number,
							VC->Calc[i].atom->coor[0],VC->Calc[i].atom->coor[1],VC->Calc[i].atom->coor[2],
							VC->Calc[i].exposed? "C ": "O ");
					*/
				}			
			}
			printf("%d atoms set as buried\n", n_buried);
			//fclose(surffile);

			for(i=0;i<FA->atm_cnt_real;i++){
				if(VC->Calc[i].score){VC->Calc[i].atom = NULL;}
			}
			//getchar();
		}
	}  
	
	// ── GIST evaluator initialization ──
	if(FA->use_gist && FA->gist_dg_file[0] != '\0' && FA->gist_dens_file[0] != '\0'){
		GISTEvaluator* gist = new GISTEvaluator();
		gist->delta_G_cutoff = FA->gist_dg_cutoff;
		gist->rho_cutoff     = FA->gist_rho_cutoff;
		gist->divisor        = FA->gist_divisor;
		gist->weight         = FA->gist_weight;
		if(gist->load_dx(FA->gist_dg_file, FA->gist_dens_file)){
			FA->gist_evaluator = gist;
			printf("GIST water displacement scoring enabled\n");
		}else{
			fprintf(stderr,"WARNING: GIST grid loading failed, disabling GIST scoring\n");
			delete gist;
			FA->use_gist = 0;
		}
	}

	if(FA->use_hbond){
		printf("Directional H-bond scoring enabled (weight=%.2f, search=%s, rank=%s)\n",
		       FA->hbond_weight,
		       FA->use_hbond_search ? "on" : "off",
		       FA->use_hbond_rank ? "on" : "off");
	}

	FA->deelig_root_node = new struct deelig_node_struct;
	FA->deelig_root_node->parent = NULL;

	FA->contributions = (float*)malloc(FA->ntypes*FA->ntypes*sizeof(float));
	if(!FA->contributions){
		fprintf(stderr,"ERROR: memory allocation error for contributions\n");
		Terminate(2);
	}
	
	//printf("Create rebuild list...\n");
	create_rebuild_list(FA,atoms,residue);

	//printf("atm_cnt=%d\tres_cnt=%d\n",FA->atm_cnt,FA->res_cnt);
	//printf("npar=%d\n",FA->npar);
	//cf=ic2cf(FA,VC,atoms,residue,cleftgrid,FA->npar,FA->opt_par);
	//for(i=0;i<FA->npar;i++){printf("[%8.3f]",FA->opt_par[i]);}
	//printf("=%8.5f\n",cf);

	// ── DUMP_POP audit refstructure (FLEXAIDDS_DUMP_POP + FLEXAIDDS_RMSDST) ──
// Enable calc_rmsd / .rrd / .pop.tsv without seeding the GA (coor[] unchanged).
// No-op unless FLEXAIDDS_DUMP_POP is truthy. Safe for AUTONOMOUS audit docks.
(void)load_dump_pop_refstructure(FA, atoms, residue);

// ── Native-pose CF diagnostic (FLEXAIDDS_SCORE_NATIVE=1) ─────────────────
	// Score the crystal/reference pose with the CF before the GA runs.
	// Answers: "Is the scorer broken for this ligand?" (cf<<0 means scorer works).
	// DatasetRunner sets FLEXAIDDS_SCORE_NATIVE=1 + FLEXAIDDS_RMSDST=<crystal.sdf>
	// per entry; parses [NATIVE_CF] from stderr.log alongside the GA results.
	// Prints one line to stderr, then returns — the GA continues normally.
	// FLEXAIDDS_NATIVE_ONLY=1: exit immediately after [NATIVE_CF] (smoke-test / CI).
	// DatasetRunner uses FLEXAIDDS_SCORE_NATIVE=1 alone so the GA still runs and
	// produces pose files that DatasetRunner parses alongside [NATIVE_CF].
	{
		const flexaids::ProtocolConfig native_proto = flexaids::ProtocolConfig::from_env();
		const bool do_native = native_proto.score_native || native_proto.native_only;
		if (do_native) {
			score_native_pose(FA, VC, atoms, residue, cleftgrid);
			if (native_proto.native_only) {
				std::exit(0);  // native-only mode: bail before GA (does not write pose files)
			}
		}
	}

	// ── Offline pool rescoring (FLEXAIDDS_RESCORE_POOL=<dir>) ───────────────
	// Score pre-emitted pose pools with the exact production CF. Full-complex
	// coordinates are restored per pose (serial-mapped), so optimisable DoF on
	// BOTH sides — receptor side chains and ligand torsions — are honoured
	// exactly as docked. The process exits after scoring; the GA never runs.
	// See LIB/rescore_pool.h for env vars and file conventions.
	if (const char* rescore_pool_dir = std::getenv("FLEXAIDDS_RESCORE_POOL")) {
		if (rescore_pool_dir[0] != '\0') {
			rescore_pool_mode(FA, VC, atoms, residue, cleftgrid);
			std::exit(0);  // rescore mode: scoring only, no docking
		}
	}

	//-----------------------------------------------------------------------------------
	snprintf(tmp_end_strfile, MAX_PATH__, "%s_INI.pdb", end_strfile);
	size_t remark_len = 0; remark[0] = '\0'; safe_remark_cat(remark, "REMARK initial structure\n", &remark_len);

	// Should execute cf-vcfunction instead to avoid rotamer change for INI conf.
	cf=ic2cf(FA,VC,atoms,residue,cleftgrid,FA->npar,FA->opt_par);
	VC->recalc = 0;

	for(i=0;i<FA->npar;i++){printf("[%8.3f]",FA->opt_par[i]);}
	printf("=%8.5f\n", get_cf_evalue(&cf, FA));
	//getchar();
  
	snprintf(tmpremark,MAX_REMARK,"REMARK CF=%8.5f\n", get_cf_evalue(&cf, FA));
	safe_remark_cat(remark,tmpremark,&remark_len);
	snprintf(tmpremark,MAX_REMARK,"REMARK CF.app=%8.5f\n", get_apparent_cf_evalue(&cf));
	safe_remark_cat(remark,tmpremark,&remark_len);

	for(i=0;i<FA->num_optres;i++){
    
		res_ptr = &residue[FA->optres[i].rnum];
		cf_ptr = &FA->optres[i].cf;

		snprintf(tmpremark,MAX_REMARK,"REMARK optimizable residue %s %c %d\n",
			res_ptr->name,res_ptr->chn,res_ptr->number);
		safe_remark_cat(remark,tmpremark,&remark_len);

		snprintf(tmpremark,MAX_REMARK,"REMARK CF.com=%8.5f\n",cf_ptr->com);
		safe_remark_cat(remark,tmpremark,&remark_len);
		snprintf(tmpremark,MAX_REMARK,"REMARK CF.sas=%8.5f\n",cf_ptr->sas);
		safe_remark_cat(remark,tmpremark,&remark_len);
		snprintf(tmpremark,MAX_REMARK,"REMARK CF.wal=%8.5f\n",cf_ptr->wal);
		safe_remark_cat(remark,tmpremark,&remark_len);
		snprintf(tmpremark,MAX_REMARK,"REMARK CF.con=%8.5f\n",cf_ptr->con);
		safe_remark_cat(remark,tmpremark,&remark_len);
		snprintf(tmpremark,MAX_REMARK,"REMARK CF.gist=%8.5f\n",cf_ptr->gist);
		safe_remark_cat(remark,tmpremark,&remark_len);
		snprintf(tmpremark,MAX_REMARK,"REMARK CF.hbond=%8.5f\n",cf_ptr->hbond);
		safe_remark_cat(remark,tmpremark,&remark_len);
		snprintf(tmpremark,MAX_REMARK,"REMARK CF.elec=%8.5f\n",cf_ptr->elec);
		safe_remark_cat(remark,tmpremark,&remark_len);
		snprintf(tmpremark,MAX_REMARK,"REMARK CF.gist_desolv=%8.5f\n",cf_ptr->gist_desolv);
		safe_remark_cat(remark,tmpremark,&remark_len);
		snprintf(tmpremark,MAX_REMARK,
			"REMARK CF.elec_gist_con_status = gated_inert_on_claim_path (use_elec default off; GIST hard-disabled; con constraints-only; CF.gist unused vs gist_desolv)\n");
		safe_remark_cat(remark,tmpremark,&remark_len);
		snprintf(tmpremark,MAX_REMARK,"REMARK Residue has an overall SAS of %.3f\n",cf_ptr->totsas);
		safe_remark_cat(remark,tmpremark,&remark_len);
		
	}
	
	for(i=0;i<FA->npar;i++){
		snprintf(tmpremark,MAX_REMARK,"REMARK [%8.3f]\n",FA->opt_par[i]);
		safe_remark_cat(remark,tmpremark,&remark_len);
	}
	snprintf(tmpremark,MAX_REMARK,"REMARK inputs: %s & %s\n",dockinp,gainp);
	safe_remark_cat(remark,tmpremark,&remark_len);
	
	if (FA->htpmode == false) {write_pdb(FA,atoms,residue,tmp_end_strfile,remark);}

	//printf("wrote initial PDB structure on %s\n",tmp_end_strfile);
	//-----------------------------------------------------------------------------------


	/* PRINTS ALL ACCEPTED ROTAMER LIST
	   for (i=0;i<FA->nflxsc;i++){
	   resnum=FA->flex_res[i].rnum;
	   for (j=1;j<=FA->res_cnt;j++){
	   if (residue[j].number==resnum){
	   printf("ROTAMERS RESIDUE %s%d%c\n-----------------\n",
	   residue[j].name,residue[j].number,residue[j].chn);
	   for (k=0;k<residue[j].trot+1;k++){
	   firstatm=residue[j].fatm[k];
	   lastatm=residue[j].latm[k];
	   printf("Rotamer[%d]\tFATM=%d\tLATM=%d\n",residue[j].rotid[k],firstatm,lastatm);
	   printf("COOR=");
	   for (l=0;l<3;l++){
	   printf("[%1.3f] ",atoms[lastatm].coor[l]);
	   }
	   printf("\n");
	   }
	   }
	   }
	   PAUSE;
	   }
	*/
  
	if(strcmp(FA->metopt,"GA") == 0)
	{
		////////////////////////////////
		////// Genetic Algorithm ///////
		////////////////////////////////

		// calculate time
		sta_timer=time(NULL);
		sta=localtime(&sta_timer);
		sta_val[0]=sta->tm_sec;
		sta_val[1]=sta->tm_min;
		sta_val[2]=sta->tm_hour;

		int n_chrom_snapshot = 0;

		// P1: create local TargetServer for this GA run (P1 wiring; default 1M conc, later from config)
		std::unique_ptr<target::TargetServer> local_ts;
		target::TargetServer* active_ts = nullptr;
		// TEMPER 0 is a SUPPORTED legacy configuration, not an error: read_input.cpp
		// detects it and deliberately forces clustering to CF, printing "does not
		// allow the consideration of conformational entropy". The grand-canonical
		// machinery is therefore already excluded by the user's own config.
		//
		// TargetServer holds a GrandPartitionFunction by value, whose constructor
		// throws on temperature <= 0. Constructing it unconditionally turned that
		// supported config into a hard fatal ("Fatal error: Temperature must be
		// positive") before the GA could start. Guard the construction instead of
		// relaxing the positive-temperature invariant in the thermodynamic classes,
		// which is correct and should stay.
		//
		// Leaving active_ts null is safe on this path: cluster() accepts the pointer
		// but never dereferences it, and the post-run [GRAND] reporting block below
		// is already wrapped in `if (active_ts)`.
		if (FA->temperature > 0) {
			target::TargetConfig tcfg;
			tcfg.temperature_K = static_cast<double>(FA->temperature);
			tcfg.default_conc_M = user_conc_M; // P3: from --conc or default
			local_ts = std::make_unique<target::TargetServer>(tcfg);
			active_ts = local_ts.get();
		} else {
			printf("[GRAND] disabled: TEMPER=0 (legacy CF-only mode)\n");
		}

		if (use_parallel_dock) {
			// ── ParallelDock: grid-decomposed parallel GA instances ──
			printf("%s=== ParallelDock mode: %d spatial regions ===%s\n", tui::violet(), parallel_dock_regions, tui::reset());

			ParallelDockConfig pdcfg;
			pdcfg.target_regions = parallel_dock_regions;
			// Pass parent gene_lim so each region inherits correct flexible-bond
			// and rotation gene limits (genes 1..N-1).  Previously these were
			// uninitialized, producing physically nonsensical conformations.
			ParallelDockManager pdm(FA, GB, VC, atoms, residue, cleftgrid, pdcfg,
			                        gene_lim, GB->num_genes, active_ts, "parallel-dock");
			pdm.decompose();
			pdm.run(ic2cf);
			auto global_thermo = pdm.aggregate();

			printf("ParallelDock: claim_validity=proxy_only F_like=%.4f minus_T_S_like=%.4f [CF proxy scale]\n",
			       global_thermo.free_energy, -FA->temperature * global_thermo.entropy);
			printf("ParallelDock: %zu regions completed\n", pdm.region_results().size());

			// ── Extract best chromosome from best region -> chrom[0] ──────────
			// Previously: summed all regions' num_snapshots into n_chrom_snapshot
			// but chrom[] still held the pre-GA initialisation state, so clustering
			// received garbage and wrote garbage PDB output.
			// Fix: allocate one chromosome, extract the globally best from the best
			// region (remap gene[0] to global grid index), put in chrom[0] so the
			// standard clustering/write_rrd path produces a valid pose file.
			//
			// NOTE: in the parallel dock path GA() is never called for the main
			// chrom pointer, so chrom == NULL here.  Allocate one chromosome.
			int pd_num_genes = FA->npar;

			chrom = (chromosome*)calloc(1, sizeof(chromosome));
			if (chrom) {
				chrom[0].genes = (gene*)calloc(pd_num_genes, sizeof(gene));
			}
			memchrom = 1;

			int pd_global_grd_idx = -1;

			if (chrom && chrom[0].genes &&
			    pdm.get_best_chromosome(chrom[0], pd_global_grd_idx)) {

				chrom[0].genes[0].to_ic = (double)pd_global_grd_idx;
				n_chrom_snapshot = 1;
			} else {
				fprintf(stderr, "ParallelDock: failed to extract best chromosome -- "
				        "no output will be written\n");
				n_chrom_snapshot = 0;
			}

			if (n_chrom_snapshot > 0) {
				// ── Wire the best chromosome into the downstream pipeline ────────
				// The regular GA path sets chrom_snapshot, gene_lim, GB->num_genes
				// automatically. For the parallel dock path we must do it manually.

				// chrom_snapshot: 1 slot only — parallel dock produces exactly 1 best
				// chromosome. Allocating num_chrom * max_generations = 3,500,000 slots
				// wastes ~350 MB and causes the cleanup loop to read 3.5M entries
				// even though only slot 0 has valid genes. The fix: allocate 1, then
				// patch GB->num_chrom/max_generations to 1 so cleanup iterates once.
				const int pd_snap_size = 1;
				chrom_snapshot = (chromosome*)calloc(pd_snap_size, sizeof(chromosome));
				if (chrom_snapshot) {
					chrom_snapshot[0] = chrom[0];  // shallow copy struct fields
					chrom_snapshot[0].genes = (gene*)calloc(pd_num_genes, sizeof(gene));
					if (chrom_snapshot[0].genes) {
						memcpy(chrom_snapshot[0].genes, chrom[0].genes,
						       pd_num_genes * sizeof(gene));
					}

				} else {
					fprintf(stderr, "ParallelDock: chrom_snapshot alloc failed\n");
					n_chrom_snapshot = 0;
				}

				// gene_lim: allocate using global FA parameters (GLOBAL grid limits).
				// Gene[0] in chrom_snapshot[0] already holds the global grid index,
				// so these limits are correct for downstream scoring + clustering.
				GB->num_genes = FA->npar;
				gene_lim = (genlim*)malloc(FA->npar * sizeof(genlim));
				if (gene_lim) {
					set_gene_lim(FA, GB, gene_lim);

				} else {
					fprintf(stderr, "ParallelDock: gene_lim alloc failed\n");
					n_chrom_snapshot = 0;
				}

				// Patch GB counters so cleanup loop iterates 1 slot (not 3,500,000).
				// Downstream code uses n_chrom_snapshot directly; this is safe to do now.
				GB->num_chrom = 1;
				GB->max_generations = 1;
			}
		} else if (use_campaign) {
			// ── ParallelCampaign: multi-ligand virtual screening ──
			printf("%s=== Campaign mode: parallel virtual screening ===%s\n", tui::violet(), tui::reset());

			auto ccfg = campaign::auto_configure(
				"", "",  // paths already loaded in FA globals
				config_path,
				output_prefix.empty() ? "campaign" : output_prefix,
				use_rigid, use_folded
			);
			ccfg.default_conc_M = user_conc_M;  // P3: forward --conc for grand canonical in campaign path
			ccfg.coarse_prefilter = use_coarse_prefilter;
			ccfg.coarse_prefilter_top_n = coarse_prefilter_top_n;
			ccfg.coarse_target_mol2 = screen_target_mol2;
			{
				const flexaids::ProtocolConfig camp_proto = flexaids::ProtocolConfig::from_env();
				ccfg.coarse_cleft_pdb = camp_proto.oracle_site.empty()
					? camp_proto.cleft_sphere_file : camp_proto.oracle_site;
			}
			auto summary = campaign::run_campaign(ccfg,
				[](int done, int total, const campaign::LigandResult& lr) {
					printf("\r  [%d/%d] %s: score_proxy=%.2f (%.1fs)",
					       done, total, lr.name.c_str(), lr.dG_corrected, lr.dock_time_sec);
					fflush(stdout);
				}
			);
			printf("\n%sCampaign complete:%s %d/%d successful, %.0f ligands/hour\n", tui::mint(), tui::reset(),
			       summary.successful, summary.total_ligands, summary.throughput_per_hour);
			n_chrom_snapshot = summary.successful;
		} else if (use_screen) {
			// ── CoarseScreen / optional Stage-2 hook (default: Stage 1 only) ──
			printf("%s=== CoarseScreen mode: cube screening (top %d)%s ===%s\n",
			       tui::violet(), screen_top_n, use_screen_dock ? " + Stage-2 hook" : "", tui::reset());

			const std::string target_mol2 = !screen_target_mol2.empty()
				? screen_target_mol2 : screen_receptor_path;
			auto screen_target = nrgrank::parse_target_mol2(target_mol2);
			if (screen_target.empty()) {
				fprintf(stderr, "[SCREEN] ERROR: could not parse target atoms from %s\n",
				        target_mol2.c_str());
				fprintf(stderr, "[SCREEN]   --screen requires a receptor in MOL2 format "
				        "(or --screen-target-mol2).\n");
				Terminate(1);
			}

			const flexaids::ProtocolConfig screen_proto = flexaids::ProtocolConfig::from_env();
			std::string site_pdb_path = screen_proto.oracle_site;
			if (site_pdb_path.empty()) site_pdb_path = screen_proto.cleft_sphere_file;
			auto screen_spheres = site_pdb_path.empty()
			    ? std::vector<nrgrank::BindingSiteSphere>{}
			    : nrgrank::parse_binding_site_pdb(site_pdb_path);
			if (screen_spheres.empty()) {
				fprintf(stderr, "[SCREEN] WARNING: no binding site spheres loaded; "
				        "set FLEXAIDDS_ORACLE_SITE or FLEXAIDDS_CLEFT_SPHERE_FILE.\n");
			}

			std::vector<nrgrank::ScreenLigand> screen_ligands;
			{
				const std::string ext = [&]() {
					auto p = screen_ligand_path.rfind('.');
					return p != std::string::npos
					    ? screen_ligand_path.substr(p)
					    : std::string{};
				}();
				if (ext == ".sdf" || ext == ".mol")
					screen_ligands = nrgrank::CoarseScreener::load_ligands_sdf(screen_ligand_path);
				else
					screen_ligands = nrgrank::CoarseScreener::load_ligands_mol2(screen_ligand_path);
			}
			if (screen_ligands.empty()) {
				fprintf(stderr, "[SCREEN] ERROR: no ligands loaded from %s\n",
				        screen_ligand_path.c_str());
				Terminate(1);
			}

			nrgrank::TwoStageScreener ts;
			nrgrank::TwoStageConfig tcfg;
			tcfg.top_n = screen_top_n;
			tcfg.coarse.top_n = screen_top_n;
			tcfg.write_coarse_csv = true;
			tcfg.output_dir = output_prefix + "_screen";
			tcfg.verbose = true;
			tcfg.stage2_kind_label = "surrogate";
			ts.set_config(tcfg);
			ts.prepare_target(screen_target, screen_spheres);
			if (!ts.coarse_screener().is_prepared()) {
				fprintf(stderr, "[SCREEN] ERROR: target preparation failed.\n");
				Terminate(1);
			}
			printf("[SCREEN] Target prepared: %d anchors; screening %d ligands...\n",
			       static_cast<int>(ts.coarse_screener().num_anchors()),
			       static_cast<int>(screen_ligands.size()));

			if (use_screen_dock) {
				// Composition placeholder — same honesty as ParallelCampaign's
				// surrogate_model_dock_score. Not a Voronoi CF and not ΔG.
				ts.set_full_dock_callback(
					[](const nrgrank::ScreenLigand& lig,
					   const nrgrank::ScreenResult&) {
						return static_cast<float>(
							-0.04 * static_cast<double>(lig.atoms.size()));
					});
			}

			auto staged = ts.run(screen_ligands);
			nrgrank::TwoStageScreener::write_unified_csv(
				tcfg.output_dir + "/unified.csv", staged);
			nrgrank::TwoStageScreener::write_screen_receipt(
				tcfg.output_dir,
				static_cast<int>(screen_ligands.size()),
				screen_top_n,
				use_screen_dock);

			printf("\n%-6s  %-10s  %s\n", "Rank", "Score", "Name");
			printf("------  ----------  ----\n");
			for (int i = 0; i < static_cast<int>(staged.size()); ++i) {
				printf("%-6d  %-10.3f  %s\n", i + 1,
				       staged[static_cast<size_t>(i)].coarse_result.score,
				       staged[static_cast<size_t>(i)].coarse_result.name.c_str());
			}

			n_chrom_snapshot = static_cast<int>(staged.size());
		} else {
			// ── Standard single search run (GA default; CMA-ES opt-in) ──
			// FLEXAIDDS_CMAES_BEGIN
			// Opt-in CMA-ES search backend. Env FLEXAIDDS_SEARCH=cmaes (or CMAES)
			// swaps the operator; scoring path (ic2cf / 5 seam fns) is unchanged.
			// Eval budget: λ×gens ≡ pop×gens (claim 1000×2000 = 2e6). See CMAES_INTEGRATION.md.
			// Allocation contract matches GA(): sets chrom, chrom_snapshot, gene_lim,
			// memchrom so existing top.cpp free / ranking loops remain valid.
			// ic2cf.cpp / gaboom.cpp are never modified by this branch.
			const char* flexaidds_search = std::getenv("FLEXAIDDS_SEARCH");
			const bool use_cmaes = flexaidds_search &&
				(std::strcmp(flexaidds_search, "cmaes") == 0 ||
				 std::strcmp(flexaidds_search, "CMAES") == 0);
			if (use_cmaes) {
				// Minimal GA-input plumbing so GB->num_chrom / max_generations / seed
				// match the claim budget path (read_gainputs when a .inp is present).
				GB->num_genes = FA->npar;
				if (GB->num_genes <= 0) {
					fprintf(stderr, "ERROR: CMA-ES: no parameters to optimize (FA->npar=0).\n");
					Terminate(1);
				}
				if (gainp[0] != 0) {
					int geninterval = 0, popszpartition = 0;
					read_gainputs(FA, GB, &geninterval, &popszpartition, gainp);
				}
				if (GB->num_chrom <= 0) GB->num_chrom = 1000;
				if (GB->max_generations <= 0) GB->max_generations = 2000;

				gene_lim = (genlim*)malloc(static_cast<size_t>(GB->num_genes) * sizeof(genlim));
				if (!gene_lim) {
					fprintf(stderr, "ERROR: CMA-ES: gene_lim allocation failed.\n");
					Terminate(1);
				}

				CmaesConfig cma_cfg;
				cma_cfg.population = GB->num_chrom > 0 ? GB->num_chrom : 1000;
				cma_cfg.max_evals = 0;
				if (const char* me = std::getenv("FLEXAIDDS_CMAES_MAX_EVALS")) {
					if (me[0] != 0) {
						char* endp = nullptr;
						const long long v = std::strtoll(me, &endp, 10);
						if (endp != me && v > 0)
							cma_cfg.max_evals = static_cast<std::int64_t>(v);
					}
				}
				if (cma_cfg.max_evals <= 0) {
					const long long budget =
						static_cast<long long>(cma_cfg.population) *
						static_cast<long long>(GB->max_generations);
					cma_cfg.max_evals = budget > 0 ? budget : 2000000LL;
				}
				cma_cfg.seed = static_cast<std::uint32_t>(GB->seed != 0 ? GB->seed : 1);
				cma_cfg.enable_entropy_trace = true;
				if (!output_prefix.empty()) {
					cma_cfg.write_trace = output_prefix + "_cmaes_entropy.csv";
				}

				CmaesResult cma_res;
				std::vector<EntropyTraceSample> cma_trace;
				const int cma_rc = cmaes_run_dock(
					FA, GB, VC, gene_lim, atoms, residue, cleftgrid,
					ic2cf, cma_cfg, &cma_res, &cma_trace);
				if (cma_rc != 0) {
					fprintf(stderr, "ERROR: cmaes_run_dock failed (rc=%d status=%d)\n",
					        cma_rc, cma_res.status);
					n_chrom_snapshot = 0;
				} else {
					if (!cma_cfg.write_trace.empty() && !cma_trace.empty()) {
						cmaes_write_trace_csv(cma_cfg.write_trace, cma_trace);
					}
					// Archive size K — each chrom[i].genes / chrom_snapshot[i].genes
					// is a separate malloc so top.cpp free loops stay correct.
					const int arch_n = static_cast<int>(cma_res.archive_genes.size());
					const int K = std::max(1, std::min(cma_cfg.archive_size,
						arch_n > 0 ? arch_n : 1));
					memchrom = K;
					chrom = (chromosome*)calloc(static_cast<size_t>(K), sizeof(chromosome));
					chrom_snapshot = (chromosome*)calloc(static_cast<size_t>(K), sizeof(chromosome));
					// Temporary contiguous storage for cmaes_fill_chromosomes only;
					// immediately re-homed into per-slot malloc buffers.
					gene* tmp_storage = (gene*)calloc(
						static_cast<size_t>(K) * static_cast<size_t>(GB->num_genes),
						sizeof(gene));
					if (!chrom || !chrom_snapshot || !tmp_storage) {
						fprintf(stderr, "ERROR: CMA-ES: chromosome/snapshot allocation failed.\n");
						n_chrom_snapshot = 0;
					} else {
						const int filled = cmaes_fill_chromosomes(
							cma_res, GB->num_genes, chrom, K, tmp_storage);
						n_chrom_snapshot = 0;
						for (int i = 0; i < filled; ++i) {
							gene* cgenes = (gene*)malloc(
								static_cast<size_t>(GB->num_genes) * sizeof(gene));
							gene* sgenes = (gene*)malloc(
								static_cast<size_t>(GB->num_genes) * sizeof(gene));
							if (!cgenes || !sgenes) {
								fprintf(stderr, "ERROR: CMA-ES: gene buffer allocation failed.\n");
								free(cgenes);
								free(sgenes);
								break;
							}
							std::memcpy(cgenes, chrom[i].genes,
								static_cast<size_t>(GB->num_genes) * sizeof(gene));
							std::memcpy(sgenes, chrom[i].genes,
								static_cast<size_t>(GB->num_genes) * sizeof(gene));
							chrom[i].genes = cgenes;
							chrom_snapshot[i] = chrom[i];
							chrom_snapshot[i].genes = sgenes;
							++n_chrom_snapshot;
						}
						free(tmp_storage);
						// Patch counters so free(chrom_snapshot) walks K slots, not
						// num_chrom * max_generations (ParallelDock uses the same fix).
						GB->num_chrom = std::max(1, n_chrom_snapshot);
						GB->max_generations = 1;
						printf("[SEARCH] backend=cmaes evals=%d n_snap=%d best_cf=%.6f "
						       "lambda=%d max_evals=%lld\n",
						       cma_res.n_evals, n_chrom_snapshot, cma_res.best_cf,
						       cma_cfg.population,
						       static_cast<long long>(cma_cfg.max_evals));
						if (!cma_cfg.write_trace.empty()) {
							printf("[SEARCH] cmaes entropy trace: %s\n",
							       cma_cfg.write_trace.c_str());
						}
					}
				}
			} else {
				GAContext ga_ctx;
				n_chrom_snapshot = GA(FA,GB,VC,&chrom,&chrom_snapshot,&gene_lim,atoms,residue,&cleftgrid,gainp,&memchrom,ic2cf, &ga_ctx);
			}
			// FLEXAIDDS_CMAES_END
		}
    
		if(n_chrom_snapshot > 0){

			end_timer=time(NULL);
			end=localtime(&end_timer);
			end_val[0]=end->tm_sec;
			end_val[1]=end->tm_min;
			end_val[2]=end->tm_hour;
      
			printf("GA:Start time =%0d:%0d:%0d\n",sta_val[2],sta_val[1],sta_val[0]);
			printf("GA:End time   =%0d:%0d:%0d\n",end_val[2],end_val[1],end_val[0]);
      
			ct=0;
			if (sta_val[0]>end_val[0]){
				end_val[1]--;
				end_val[0]+=60;
			}
			if (sta_val[1]>end_val[1]){
				end_val[2]--;
				end_val[1]+=60;
			}
			ct+=((end_val[0]-sta_val[0])+(end_val[1]-sta_val[1])*60);
      
			printf("GA Computational time %ld sec (%4.2f min)\n",ct,(double)ct/60.0);
      
			printf("atoms recalculated=%d\n",FA->recalci);
			printf("individuals skipped=%d\n",FA->skipped);
			printf("individuals clashed=%d\n",FA->clashed);
			
			////////////////////////////////
			//////       END         ///////
			////////////////////////////////
      
			/******************************************************************/

				// Independently rebuild and re-score every retained chromosome on the
				// serial master state. Clustering, thermodynamics, and PDB emission must
				// never inherit a stale score from an OpenMP worker or a prior geometry.
				if (n_chrom_snapshot > 0) {
					const bool rank_hbond = FA->use_hbond && FA->use_hbond_rank &&
					                        !FA->use_hbond_search;
					double max_score_delta = 0.0;
					int inconsistent_scores = 0;
					printf("Exact-pose re-scoring %d retained chromosomes%s...\n",
					       n_chrom_snapshot,
					       rank_hbond ? " with H-bond rank term" : "");
					for (int si = 0; si < n_chrom_snapshot; ++si) {
						if (FA->ring_flex_active) {
							for (int s = 0; s < FA->ring_n_sugars && s < MAX_RING_FLEX; ++s)
								FA->ring_cur_phases[s] = chrom_snapshot[si].ring_phases[s];
						}
						const double search_score = chrom_snapshot[si].evalue;
						FA->hbond_rank_rescore = rank_hbond ? 1 : 0;
						cfstr exact_cf = eval_chromosome(
						    FA, GB, VC, gene_lim, atoms, residue, cleftgrid,
						    chrom_snapshot[si].genes, ic2cf);
						const double exact_score = get_cf_evalue(&exact_cf, FA);
						const double delta = std::abs(exact_score - search_score);
						max_score_delta = std::max(max_score_delta, delta);
						if (!std::isfinite(exact_score) || delta > 1e-4) {
							++inconsistent_scores;
							if (inconsistent_scores <= 20) {
								printf("WARNING: retained chromosome %d search CF=%.8f "
								       "exact-pose CF=%.8f delta=%.8f\n",
								       si, search_score, exact_score, delta);
							}
						}
						// Fail closed: even a clash penalty is the score of the exact pose.
						// Retaining the old value would reintroduce score/geometry mismatch.
						chrom_snapshot[si].cf = exact_cf;
						chrom_snapshot[si].evalue = exact_score;
						chrom_snapshot[si].app_evalue = get_apparent_cf_evalue(&exact_cf);
					}
					FA->hbond_rank_rescore = 0;
					QuickSort(chrom_snapshot, 0, n_chrom_snapshot - 1, true);
					printf("Exact-pose score audit: inconsistent=%d/%d max_delta=%.8f\n",
					       inconsistent_scores, n_chrom_snapshot, max_score_delta);
				}

				// ── Post-GA ensemble thermodynamic summary ──
				if (FA->temperature > 0 && n_chrom_snapshot > 0) {
					statmech::StatMechEngine post_engine(
						static_cast<double>(FA->temperature),
						statmech::make_contact_function_optimizer_provenance());
					for (int si = 0; si < n_chrom_snapshot; si++) {
						post_engine.add_sample(chrom_snapshot[si].evalue);
					}
					auto post_thermo = post_engine.compute();
					printf("\n%s======= Post-GA CF-proxy ensemble diagnostics (%sT parameter=%uK%s) =======%s\n",
					       tui::tangerine(), tui::T(), FA->temperature, tui::tangerine(), tui::reset());
					printf("  claim_validity = proxy_only\n");
					printf("  F-like proxy   = %10.4f [legacy transform]\n", post_thermo.free_energy);
					printf("  Mean CF       = %10.4f [CF units]\n", post_thermo.mean_energy);
					printf("  S-like value  = %10.6f [proxy scale/K]\n", post_thermo.entropy);
					printf("  -T*S-like     = %10.4f [proxy scale]\n", -static_cast<double>(FA->temperature) * post_thermo.entropy);
					printf("  Heat capacity  = %10.4f\n", post_thermo.heat_capacity);
					printf("  CF std dev     = %10.4f [CF units]\n", post_thermo.std_energy);
					printf("  Ensemble size  = %d\n", n_chrom_snapshot);
					printf("%s========================================================%s\n\n", tui::tangerine(), tui::reset());
				}

			printf("clustering all individuals in GA...");
			fflush(stdout);

			printf("n_chrom_snapshot=%d\n", n_chrom_snapshot);
			fflush(stdout);

			if( strcmp(FA->clustering_algorithm,"FO") == 0 )
			{
				printf("using the Fast OPTICS (FO) density based clustering algorithm.\n");
				FastOPTICS_cluster(FA,GB,VC,chrom_snapshot,gene_lim,atoms,residue,cleftgrid,n_chrom_snapshot,end_strfile,tmp_end_strfile,dockinp,gainp, active_ts, "ga-ligand");
			}
			else if( strcmp(FA->clustering_algorithm,"DP") == 0 )
			{
				printf("using the Density Peak (DP) based clustering algorithm.\n");
				DensityPeak_cluster(FA,GB,VC,chrom_snapshot,gene_lim,atoms,residue,cleftgrid,n_chrom_snapshot,end_strfile,tmp_end_strfile,dockinp,gainp, active_ts, "ga-ligand");
			}
			else
			{
				printf("using the Complementarity Function (CF) based clustering algorithm.\n");
				cluster(FA,GB,VC,chrom_snapshot,gene_lim,atoms,residue,cleftgrid,n_chrom_snapshot,end_strfile,tmp_end_strfile,dockinp,gainp, active_ts, "ga-ligand");
			}

			// P1/P5: augment output with grand canonical info if ts active.
			// LigandRank fields: name, log_Z, dG, p_bound (see GrandPartitionFunction.h).
			// p_bound is CF-proxy occupancy (p_bind_like), not calibrated ΔG occupancy.
			if (active_ts) {
				printf("[GRAND] sessions=%d claim_validity=proxy_only\n", active_ts->completed_sessions());
				auto ranks = active_ts->rank_ligands();
				for (const auto& r : ranks) {
					printf("[GRAND] %s: log_Z=%.6g p_bind_like=%.6g p_bind=%.6g dG=%.6g claim_validity=proxy_only\n",
					       r.name.c_str(), r.log_Z, r.p_bound, r.p_bound, r.dG);
				}
				// P5: write sidecar .grand.txt for richer output (Xi, p_bind_like per ligand)
				char grandfile[512];
				// end_strfile is a stack array, so it is never null; it is
				// initialized to "flexaid" at declaration, which is where the
				// fallback belongs. -Wpointer-bool-conversion flagged the old
				// ternary as always-true, and it was.
				snprintf(grandfile, sizeof(grandfile), "%s.grand.txt", end_strfile);
				FILE* gf = fopen(grandfile, "w");
				if (gf) {
					fprintf(gf, "# Grand canonical summary (P3/P5) claim_validity=proxy_only (CF-proxy Z; p_bind_like)\n");
					for (const auto& r : ranks) {
						fprintf(gf, "ligand=%s log_Z=%.6g p_bind_like=%.6g p_bind=%.6g dG=%.6g claim_validity=proxy_only\n",
						        r.name.c_str(), r.log_Z, r.p_bound, r.p_bound, r.dG);
					}
					fclose(gf);
				}
				// P5: also emit as REMARK GRAND for parsers (per plan)
				printf("REMARK GRAND_SESSIONS %d claim_validity=proxy_only\n", active_ts->completed_sessions());
				for (const auto& r : ranks) {
					printf("REMARK GRAND %s log_Z=%.6g p_bind_like=%.6g p_bind=%.6g dG=%.6g claim_validity=proxy_only\n",
					       r.name.c_str(), r.log_Z, r.p_bound, r.p_bound, r.dG);
				}
			}

			//////////////////////////////////////////
			// Looking at cleftgrid chrom's density //
			//////////////////////////////////////////
// 			int* gridcount;
// 			gridcount = (int*) malloc(FA->MIN_CLEFTGRID_POINTS * sizeof(int));
// 			if(!gridcount)
// 			{
// 				fprintf(stderr, "ERROR: memory allocation error for gridcount\n");
// 				Terminate(2);
// 			}
//             for(i = 0; i < FA->MIN_CLEFTGRID_POINTS; ++i)
//             {
//                 gridcount[i] = 0;
//                 cleftgrid[i].number = 0;
//             }
// 			for(i = 0; i < n_chrom_snapshot; ++i)
// 			{
// 				gridcount[(unsigned int)chrom_snapshot[i].genes[0].to_ic]++;
//                 cleftgrid[(unsigned int)chrom_snapshot[i].genes[0].to_ic].number++;
// 			}
//             std::sort(&gridcount[0],&gridcount[FA->MIN_CLEFTGRID_POINTS-1]);
// 		    for(i = 0, j = 0; j < FA->MIN_CLEFTGRID_POINTS; ++j)
// 		    {
//                  if(gridcount[j] > 0) {printf("%d: %d\n",j,gridcount[j]); ++i;}
// //                if(cleftgrid[j].number > 0) { printf("%d: %d\n",j,cleftgrid[j].number); ++i;}
// 		    }
//             printf("there is a total of :\n\t%d occupied grid points\n\t%d empty grid points\n\t%d total grid points\n",i, FA->MIN_CLEFTGRID_POINTS-i, FA->MIN_CLEFTGRID_POINTS);
//             // Grid Density Count (free-ing) 
//             if(gridcount != NULL) free(gridcount);
		}
	}
    
	//////////////////////////////////////////
	// free up memory allocated using malloc//
	//////////////////////////////////////////
	printf("free-ing up memory\n");
	
	// Genes properties
	if(gene_lim != NULL) free(gene_lim);
	
	// Chromosomes
	if(chrom != NULL){
		for(i=0;i<memchrom;++i){
			if(chrom[i].genes != NULL) free(chrom[i].genes);
		}
		free(chrom);
	}
	
	if(chrom_snapshot != NULL){
		for(i=0;i<(GB->num_chrom*GB->max_generations);++i){
			if(chrom_snapshot[i].genes != NULL) free(chrom_snapshot[i].genes);
		}
		free(chrom_snapshot);
	}
	
	// Vcontacts
	if(VC->Calc != NULL) {
		free(VC->Calc);
		free(VC->Calclist);
		free(VC->ca_index);
		free(VC->seed);
		free(VC->contlist);
	}

	// Cleft Grid
	if(cleftgrid != NULL) free(cleftgrid);

	// Atoms
	if(atoms != NULL) {
	  
		for(i=0;i<FA->MIN_NUM_ATOM;i++){

			if(atoms[i].cons != NULL) { free(atoms[i].cons); }
			if(atoms[i].coor_ref != NULL) { free(atoms[i].coor_ref); }

			if(atoms[i].eigen != NULL){
				for(j=0;j<FA->normal_modes;j++)
					if(atoms[i].eigen[j] != NULL)
						free(atoms[i].eigen[j]);
				free(atoms[i].eigen);
			}
		}

		free(atoms);
		
	}
  
	free(FA->num_atm);
	
	// de-allocate flat arrays (A2) and energy_matrix
	for(int _k = 0; _k < FA->ntypes * FA->ntypes; ++_k)
		if(FA->energy_matrix[_k].flat_x) free(FA->energy_matrix[_k].flat_x);
	free(FA->energy_matrix);

	// Constraints
	if(FA->constraints != NULL) free(FA->constraints);

	// Residues
	if(residue != NULL) {
		// From 0, not 1: slot 0 is allocated by read_pdb.cpp:44-45 and its
		// fatm/latm went unreleased for as long as the loop started at 1.
		// free_resid guards the natm derivation on the trio, so the slot that
		// has only fatm/latm needs no special case here -- the reason that
		// guard exists is documented at LIB/free_resid.cpp:22.
		for(i=0;i<=FA->res_cnt;i++){
			free_resid(&residue[i]);
		}

		free(residue);
	}

	// Mov (buildlist)
	for(i=0;i<2;i++){ if(FA->mov[i] != NULL) free(FA->mov[i]); }

	// Optimizable residues
	if(FA->optres != NULL) free(FA->optres);

	// Rotamers
	if(rotamer != NULL) free(rotamer);
  
	// Flexible Residues
	if (FA->flex_res != NULL) {
		for(i=0;i<FA->MIN_FLEX_RESIDUE;i++){
			if(FA->flex_res[i].close != NULL) {
				free(FA->flex_res[i].close);
			}
		}
		free(FA->flex_res);
	}

	// eigen vectors
	if(FA->eigenvector != NULL){
		for(i=0;i<3*FA->MIN_NUM_ATOM;i++)
			if(FA->eigenvector[i] != NULL) 
				free(FA->eigenvector[i]);
		
		free(FA->eigenvector);
	}

	// normal grid
	if(FA->normal_grid != NULL){
		for(i=0;i<FA->MIN_NORMAL_GRID_POINTS;i++)
			if(FA->normal_grid[i] != NULL)
				free(FA->normal_grid[i]);
		
		free(FA->normal_grid);
	}
  
	// Param
	if(FA->map_par != NULL) free(FA->map_par);
	if(FA->opt_par != NULL) free(FA->opt_par);
	if(FA->del_opt_par != NULL) free(FA->del_opt_par);
	if(FA->min_opt_par != NULL) free(FA->min_opt_par);
	if(FA->max_opt_par != NULL) free(FA->max_opt_par);

	/*
	// RMSD
	for(i=0;i<FA->num_het;i++){
	if(FA->res_rmsd[i].fatm != NULL) free(FA->res_rmsd[i].fatm);
	if(FA->res_rmsd[i].latm != NULL) free(FA->res_rmsd[i].latm);  // was .fatm: double-free + leaked .latm
	}
	//if(FA->atoms_rmsd != NULL) free(FA->atoms_rmsd);
	*/

  
	// FlexDEE Nodes
	if ( FA->psFlexDEENode != NULL ) {
		FA->psFlexDEENode = FA->psFlexDEENode->last;

		while( FA->psFlexDEENode->prev != NULL ) {
			free(FA->psFlexDEENode->rotlist);
			FA->psFlexDEENode = FA->psFlexDEENode->prev;
			free(FA->psFlexDEENode->next);
		}
    
		free(FA->psFlexDEENode->rotlist);
		free(FA->psFlexDEENode);

	}

	if(VC != NULL){
		free(VC->ptorder);
		free(VC->centerpt);
		free(VC->poly);
		free(VC->cont);
		free(VC->vedge);
		free(VC->ca_rec);
		delete VC;
	}

	if(GB != NULL) { delete GB; }

	if(FA != NULL) {
		free(FA->contacts);
		free(FA->mif_energies);
		free(FA->mif_sorted);
		free(FA->mif_cdf);
		free(FA->reflig_nearest_grid);
		free(FA->coarse_seeds_grid);
		free(FA->coarse_seeds_genes);
		delete FA->thermo_engine;
		FA->thermo_engine = nullptr;
		delete FA;
	}


	//////////////////////////////////////////
	/////////////////  END   /////////////////
	//////////////////////////////////////////

	printf("Done.\n");

	return (0);
  } catch (const FlexAIDException& e) {
	if (e.exit_code() == 0) return 0;
	fprintf(stderr, "%sFlexAID Error:%s %s\n", tui::err::failtext(), tui::err::reset(), e.what());
	return e.exit_code();
  } catch (const std::exception& e) {
	fprintf(stderr, "Fatal error: %s\n", e.what());
	return 1;
  }
}
