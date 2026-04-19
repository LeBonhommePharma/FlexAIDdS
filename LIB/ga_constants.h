#ifndef GA_CONSTANTS_H
#define GA_CONSTANTS_H

// Named constants for hardcoded numeric literals in gaboom.cpp.
// Extracted for clarity; values are unchanged from the original code.

// ── Sleep & polling ─────────────────────────────────────────────────
constexpr int GA_SLEEP_MS = 25;                // milliseconds between state file checks
constexpr int GA_STATE_CHECK_INTERVAL = 1;     // seconds between pause/abort/stop checks

// ── Grid partitioning defaults ──────────────────────────────────────
constexpr int GA_DEFAULT_GEN_INTERVAL = 50;    // generation interval for grid repartitioning
constexpr int GA_DEFAULT_POP_PARTITION = 100;  // population partition size for grid operations

// ── Default GA parameter values (when reading from config) ──────────
constexpr int GA_DEFAULT_NUM_PRINT = 10;       // default number of chromosomes to print
constexpr int GA_DEFAULT_PRINT_INT = 1;        // default print interval (every generation)
constexpr int GA_DEFAULT_SEED = 0;             // 0 = use current time as seed

// ── Entropy convergence defaults ────────────────────────────────────
constexpr int GA_DEFAULT_ENTROPY_CHECK_INTERVAL = 10;   // check every N generations
constexpr int GA_DEFAULT_ENTROPY_WINDOW = 5;            // plateau detection window size
constexpr double GA_DEFAULT_ENTROPY_REL_THRESHOLD = 0.01; // relative convergence threshold

// ── Thermodynamic analysis ──────────────────────────────────────────
constexpr double GA_DEFAULT_TEMPERATURE_K = 300.0; // fallback temperature when FA->temperature == 0

// ── TurboQuant ensemble compression ─────────────────────────────────
constexpr int GA_TQENS_MIN_SNAPSHOTS = 64;     // minimum snapshots for TurboQuant compression
constexpr int GA_TQENS_BITS = 3;               // bits/coordinate (97% fidelity, 10.7x compression)
constexpr int GA_TQENS_ENERGY_DIM = 4;         // energy descriptor dimensions (com, wal, sas, elec)

// ── FastOPTICS super-clustering ─────────────────────────────────────
constexpr int GA_FOPTICS_MIN_POINTS = 4;       // minimum minPts for FastOPTICS
constexpr int GA_FOPTICS_DIVISOR = 20;         // minPts = max(MIN_POINTS, n / DIVISOR)

// ── Fitness statistics ──────────────────────────────────────────────
constexpr double GA_FITNESS_DENOM_FLOOR = 1e-15; // minimum denominator to prevent division by zero

// ── Logging ─────────────────────────────────────────────────────────
constexpr int GA_SMFREE_LOG_INTERVAL = 50;     // log SMFREE thermodynamics every N generations

// ── Duplicate gene tolerance ────────────────────────────────────────
constexpr double GA_GENE_MATCH_TOLERANCE = 0.1; // genes within this delta are considered identical

// ── DEE ligand conformer filtering ──────────────────────────────────
constexpr int GA_DEELIG_SENTINEL = -1000;       // "no dihedral assigned" sentinel in deelig_list[]
constexpr int GA_MAX_DEELIG_DIHEDRALS = 100;    // max flexible dihedrals for DEE conformer lists

// ── TurboQuant contact matrix (TQCM) ──────────────────────────────
constexpr double GA_TQCM_SAMPLE_AREA = 0.5;     // energy matrix spline sample point
constexpr int GA_TQCM_BIT_WIDTH = 2;            // bits/entry for QuantizedContactMatrix (16x compression)
constexpr int GA_TQCM_MAX_SPOT_CHECKS = 1000;   // max type-pair validations after TQCM build

// ── tENCoM thresholds ───────────────────────────────────────────────
constexpr int GA_TENCM_MIN_RESIDUES = 6;         // min residues for torsional ENM analysis

// ── Contact list workspace ──────────────────────────────────────────
constexpr int GA_CONTLIST_SIZE = 10000;           // per-thread Vcontacts contact list workspace size

// ── FlexDEE rotamer comparison ──────────────────────────────────────
constexpr int GA_MAX_FLEXDEE_PARAMS = 100;        // max rotamer parameters in cmp_chrom2rotlist

// ── NATURaL co-translational defaults ───────────────────────────────
constexpr double GA_NATURAL_DEFAULT_TEMP = 310.0; // default NATURaL temperature (body temp, K)

// ── InStreamClustering (online medoid clustering during GA) ──────────
constexpr int   GA_INSTREAM_INTERVAL        = 100;   // merge elites every N generations
constexpr int   GA_INSTREAM_TOP_K           = 50;    // top K individuals to merge per batch
constexpr int   GA_INSTREAM_MAX_MEDOIDS     = 5000;  // upper bound on medoid count
constexpr float GA_INSTREAM_RMSD_THRESHOLD  = 2.0f;  // RMSD threshold for cluster membership (Angstrom)

#endif // GA_CONSTANTS_H
