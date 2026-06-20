// tests/stubs.cpp — minimal stubs for functions referenced by LIB sources
// but not needed by unit tests. Prevents linker errors when compiling
// individual source files in isolation.

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <utility>
#include "../LIB/gaboom.h"

// ---------------------------------------------------------------------------
// fileio.cpp
// ---------------------------------------------------------------------------
void Terminate(int status) {
    std::fprintf(stderr, "Terminate(%d) called in test stub\n", status);
    std::exit(status);
}
int OpenFile_B(char*, const char*, FILE**) { return 0; }
void CloseFile_B(FILE**, const char*) {}

// ---------------------------------------------------------------------------
// ic2cf.cpp, vcfunction.cpp (separate from Vcontacts.cpp)
// ---------------------------------------------------------------------------
cfstr ic2cf(FA_Global*, VC_Global*, atom*, resid*, gridpoint*, int, double*) {
    cfstr cf{}; return cf;
}
double vcfunction(FA_Global*, VC_Global*, atom*, resid*,
                  std::vector<std::pair<int,int>>&, bool*) {
    return 0.0;
}

// ---------------------------------------------------------------------------
// write_pdb.cpp, calc_rmsd.cpp, etc. — referenced by BindingMode.cpp
// ---------------------------------------------------------------------------
int write_pdb(FA_Global*, atom*, resid*, char[], char[]) { return 0; }
int write_MODEL_pdb(bool, bool, int, FA_Global*, atom*, resid*, char[], char[]) { return 0; }
float calc_rmsd(FA_Global*, atom*, resid*, gridpoint*, int, const double*, bool) { return 0.0f; }
double get_cf_evalue(cfstr* cf, FA_Global* /*FA*/) {
    return cf ? cf->com + cf->wal + cf->sas + cf->con + cf->elec : 0.0;
}
double get_apparent_cf_evalue(cfstr* cf) { return cf ? cf->com + cf->wal + cf->sas + cf->elec : 0.0; }
void write_contributions(FA_Global*, FILE*, bool) {}

// ---------------------------------------------------------------------------
// geometry/build functions — referenced by gaboom.cpp
// ---------------------------------------------------------------------------
void buildcc(FA_Global*, atom*, int, int[]) {}
void buildic(FA_Global*, atom*, resid*, int) {}
void build_rotamers(FA_Global*, atom**, resid*, rot*) {}
void bondedlist(atom*, int, int, int*, int*, int*) {}
void update_bonded(resid*, int, int, int*, int*) {}
void shortest_path(resid*, int, atom*) {}
void assign_shortflex(resid*, int, int, atom*) {}
int check_clash(FA_Global*, atom*, resid*, int, int, int[]) { return 0; }
void create_rebuild_list(FA_Global*, atom*, resid*) {}

// ---------------------------------------------------------------------------
// dee_pivot.cpp
// ---------------------------------------------------------------------------
int dee_pivot(psFlexDEE_Node, psFlexDEE_Node*, int, int, int, int, int) { return 0; }

// ---------------------------------------------------------------------------
// cluster.cpp, DensityPeak_Cluster.cpp, FastOPTICS_cluster.cpp
// ---------------------------------------------------------------------------
void cluster(FA_Global*, GB_Global*, VC_Global*, chromosome*, genlim*,
             atom*, resid*, gridpoint*, int, char*, char*, char*, char*) {}
void DensityPeak_cluster(FA_Global*, GB_Global*, VC_Global*, chromosome*, genlim*,
                         atom*, resid*, gridpoint*, int, char*, char*, char*, char*) {}
void FastOPTICS_cluster(FA_Global*, GB_Global*, VC_Global*, chromosome*, genlim*,
                        atom*, resid*, gridpoint*, int, char*, char*, char*, char*) {}
int write_DensityPeak_rrd(FA_Global*, GB_Global*, const chromosome*, const genlim*,
                          atom*, resid*, gridpoint*, ClusterChrom*, DPcluster*,
                          float*, char[]) { return 0; }

// ---------------------------------------------------------------------------
// vcfunction.cpp — get_yval (energy matrix lookup)
// ---------------------------------------------------------------------------
double get_yval(energy_matrix*, double) { return 0.0; }

// ---------------------------------------------------------------------------
// partition_grid.cpp, slice_grid.cpp, write_grid/generate_grid
// ---------------------------------------------------------------------------
void partition_grid(FA_Global*, chromosome*, genlim*, atom*, resid*, gridpoint**, int, int) {}
void slice_grid(FA_Global*, genlim*, atom*, resid*, gridpoint**) {}
void write_grid(FA_Global*, const gridpoint*, char[]) {}

// ---------------------------------------------------------------------------
// InStreamClustering — stubs replaced by real LIB/InStreamClustering.cpp
// ---------------------------------------------------------------------------
#if 0
#include "../LIB/InStreamClustering.h"
namespace flexaids {
InStreamCluster::InStreamCluster(float, int, int) {}
void InStreamCluster::merge_elites(const float*, const double*, int, int, int) {}
std::vector<ClusterMedoid> InStreamCluster::finalize() { return {}; }
int InStreamCluster::cluster_count() const { return 0; }
void InStreamCluster::reset() {}
void InStreamCluster::merge_elites_from_chrom(const void*, int, int, int, int) {}
float InStreamCluster::gene_rmsd(const float*, const float*, int) { return 0.0f; }
void InStreamCluster::merge_closest_pair(int) {}
} // namespace flexaids
#endif
