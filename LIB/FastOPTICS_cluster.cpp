#include "FOPTICS.h"
#include "fast_optics.hpp"
#include "MinibatchSampler.h"
#include "TargetServer.h"
#include <set>
#include <string>

void FastOPTICS_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp, target::TargetServer* ts, const std::string& ligand_name)
{
    // Adaptive minPoints based on conformational diversity.
    // Compute the ratio of distinct CF values to total snapshots.
    // When the landscape is flat (few distinct CFs), minPoints must be small
    // to find clusters.  When diverse, larger minPoints gives robust clusters.
    int minPoints;
    {
        std::set<double> distinct_cf;
        for (int i = 0; i < nChrom; ++i)
            distinct_cf.insert(chrom[i].evalue);
        double diversity_ratio = static_cast<double>(distinct_cf.size()) / static_cast<double>(nChrom);

        // Map diversity [0,1] → minPoints fraction [0.5%, 5%] of population
        //   diversity=0.001 (flat)  → minPoints ≈ 0.5% × nChrom (very small)
        //   diversity=0.50  (rich)  → minPoints ≈ 2.8% × nChrom
        //   diversity=1.00  (max)   → minPoints ≈ 5% × nChrom
        double fraction = 0.005 + 0.045 * diversity_ratio;
        minPoints = std::max(5, std::min(50, static_cast<int>(fraction * nChrom)));
    }

    // Optional super-cluster pre-filter using lightweight FastOPTICS.
    // Identifies the dominant energy basin and compacts filtered poses
    // to the front of the chrom array so downstream OPTICS runs operate
    // on a cleaner, smaller ensemble (~40% faster Shannon entropy collapse).
    if (FA->use_super_cluster && nChrom > 4) {
        std::vector<fast_optics::Point> energy_pts(nChrom);
        for (int i = 0; i < nChrom; ++i)
            energy_pts[i].coords = { chrom[i].evalue };

        fast_optics::FastOPTICS sc_optics(energy_pts, std::max(4, nChrom / 20));
        auto sc_indices = sc_optics.extractSuperCluster(fast_optics::ClusterMode::SUPER_CLUSTER_ONLY);

        if (!sc_indices.empty() && sc_indices.size() < static_cast<size_t>(nChrom)) {
            // Mark which chromosomes belong to the super-cluster
            std::vector<bool> in_sc(nChrom, false);
            for (size_t idx : sc_indices)
                in_sc[idx] = true;

            // Compact: swap super-cluster members to front of array
            int write_pos = 0;
            for (int i = 0; i < nChrom; ++i) {
                if (in_sc[i]) {
                    if (i != write_pos)
                        std::swap(chrom[write_pos], chrom[i]);
                    ++write_pos;
                }
            }

            printf("--- SuperCluster pre-filter: %zu / %d poses in dominant basin ---\n",
                   sc_indices.size(), nChrom);
            nChrom = static_cast<int>(sc_indices.size());
        }
    }

    // ── Minibatch pre-filter (farthest-point sampling) ─────────────────────
    // When the population exceeds a configurable threshold, reduce to ~5K
    // diverse representatives before clustering.  This turns O(N^2) clustering
    // into O(k*N + k^2) where k << N, giving ~100x speedup for large ensembles.
    {
        const int MINIBATCH_THRESHOLD = 10000;  // only activate above this size
        const int MINIBATCH_TARGET    = 5000;   // target representative count

        if (nChrom > MINIBATCH_THRESHOLD) {
            // Build coordinate cache (same pattern as cluster.cpp)
            const int nAtoms_mb = residue[atoms[FA->map_par[0].atm].ofres].latm[0]
                                - residue[atoms[FA->map_par[0].atm].ofres].fatm[0] + 1;
            const int stride_mb = nAtoms_mb * 3;

            minibatch::CoordCache coord_cache;
            coord_cache.n_chrom = nChrom;
            coord_cache.stride  = stride_mb;
            coord_cache.data.resize(static_cast<std::size_t>(nChrom) * stride_mb);

            for (int c = 0; c < nChrom; ++c) {
                if (c + 1 < nChrom) {
                    calc_rmsd_chrom(FA, GB, chrom, gene_lim, atoms, residue, cleftgrid,
                                    GB->num_genes, c, c + 1,
                                    &coord_cache.data[static_cast<std::size_t>(c) * stride_mb],
                                    &coord_cache.data[static_cast<std::size_t>(c + 1) * stride_mb], false);
                } else {
                    calc_rmsd_chrom(FA, GB, chrom, gene_lim, atoms, residue, cleftgrid,
                                    GB->num_genes, c, c,
                                    &coord_cache.data[static_cast<std::size_t>(c) * stride_mb], NULL, false);
                }
            }

            // Collect energies
            std::vector<double> energies(static_cast<std::size_t>(nChrom));
            for (int i = 0; i < nChrom; ++i)
                energies[i] = chrom[i].app_evalue;

            // Run farthest-point sampling
            auto sample = minibatch::MinibatchSampler::farthest_point_sample(
                coord_cache, energies.data(), nAtoms_mb, MINIBATCH_TARGET, /*verbose=*/true);

            if (sample.n_selected > 0 && sample.n_selected < nChrom) {
                // Compact selected chromosomes to front of array
                const auto& sel = sample.selected_indices;
                std::vector<bool> is_selected(static_cast<std::size_t>(nChrom), false);
                for (int idx : sel)
                    is_selected[idx] = true;

                int write_pos = 0;
                for (int i = 0; i < nChrom; ++i) {
                    if (is_selected[i]) {
                        if (i != write_pos)
                            std::swap(chrom[write_pos], chrom[i]);
                        ++write_pos;
                    }
                }
                nChrom = sample.n_selected;
            }
        }
    }

    // BindingPopulation() : BindingPopulation constructor *non-overridable*
    BindingPopulation Population1(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
    BindingPopulation Population2(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
    BindingPopulation Population3(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
 //    BindingPopulation::BindingPopulation Population4(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
	// BindingPopulation::BindingPopulation Population5(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
    
    // FastOPTICS() : calling FastOPTICS constructors
    FastOPTICS Algo1(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population1, minPoints);
    minPoints = std::floor(minPoints * 1.5);
    FastOPTICS Algo2(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population2, minPoints);
    minPoints = std::floor(minPoints * 1.5);
    FastOPTICS Algo3(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population3, minPoints);
    // minPoints = std::floor(minPoints * 1.5);
    // FastOPTICS::FastOPTICS Algo4(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population4, minPoints);
    // minPoints = std::floor(minPoints * 1.5);
    // FastOPTICS::FastOPTICS Algo5(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population5, minPoints);
    
    // 	1. Partition Sets using Random Vectorial Projections
    // 	2. Calculate Neighborhood
    // 	3. Calculate reachability distance
    // 	4. Compute the Ordering of Points To Identify Cluster Structure (OPTICS)
    // 	5. Populate BindingPopulation::Population after analyzing OPTICS
    Algo1.Execute_FastOPTICS(end_strfile, tmp_end_strfile);
    Algo2.Execute_FastOPTICS(end_strfile, tmp_end_strfile);
    Algo3.Execute_FastOPTICS(end_strfile, tmp_end_strfile);
    // Algo4.Execute_FastOPTICS(end_strfile, tmp_end_strfile);
    // Algo5.Execute_FastOPTICS(end_strfile, tmp_end_strfile);

    // Algo1.output_OPTICS(end_strfile, tmp_end_strfile);
    // Algo2.output_OPTICS(end_strfile, tmp_end_strfile);
    // Algo3.output_OPTICS(end_strfile, tmp_end_strfile);
    // Algo4.output_OPTICS(end_strfile, tmp_end_strfile);
    // Algo5.output_OPTICS(end_strfile, tmp_end_strfile);

    // output the 3D poses ordered with Fast OPTICS (done only once for the purpose as the order should not change)
    // Algo1.output_3d_OPTICS_ordering(end_strfile, tmp_end_strfile);
    // Algo2.output_3d_OPTICS_ordering(end_strfile, tmp_end_strfile);
    // Algo3.output_3d_OPTICS_ordering(end_strfile, tmp_end_strfile);
    // Algo4.output_3d_OPTICS_ordering(end_strfile, tmp_end_strfile);
    // Algo5.output_3d_OPTICS_ordering(end_strfile, tmp_end_strfile);
    
    std::cout << "Size of Population 1 is " << Population1.get_Population_size() << " Binding Modes." << std::endl;
    std::cout << "Size of Population 2 is " << Population2.get_Population_size() << " Binding Modes." << std::endl;
    std::cout << "Size of Population 3 is " << Population3.get_Population_size() << " Binding Modes." << std::endl;
    // std::cout << "Size of Population 4 is " << Population4.get_Population_size() << " Binding Modes." << std::endl;
    // std::cout << "Size of Population 5 is " << Population5.get_Population_size() << " Binding Modes." << std::endl;
    
    // output FA->max_result BindingModes
    Population1.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp, Algo1.get_minPoints());
    Population2.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp, Algo2.get_minPoints());
    Population3.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp, Algo3.get_minPoints());
    // Population4.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp, Algo4.get_minPoints());
    // Population5.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp, Algo5.get_minPoints());

    // P1 cluster hook: if TargetServer provided, register real log_Z from BindingPopulation ensemble (uses total_energy for CCBM)
    if (ts && !ligand_name.empty()) {
        auto sess = ts->create_session(ligand_name);
        sess.completed = true;
        sess.n_poses = Population1.get_Population_size();
        sess.log_Z = Population1.get_log_Z();  // P1 accessor: full ensemble, not dG approx
        // TODO (later chunks): sess.best_energy, best_center, conformer_populations from modes
        ts->register_result(sess);
    }

    printf("-- end of FastOPTICS_cluster --\n");
}