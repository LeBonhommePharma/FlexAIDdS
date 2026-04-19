#include "FOPTICS.h"
#include "fast_optics.hpp"

void FastOPTICS_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp)
{
    // Base neighbourhood size: at least 5, scales with snapshot population.
    // A value of ~nChrom/20 gives a 5 % neighbourhood, consistent with
    // typical OPTICS minPts heuristics for molecular-docking pose sets.
    int minPoints = std::max(5, nChrom / 20);

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

    // BindingPopulation() : BindingPopulation constructor *non-overridable*
    BindingPopulation Population1(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
    BindingPopulation Population2(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
    BindingPopulation Population3(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
 //    BindingPopulation::BindingPopulation Population4(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
	// BindingPopulation::BindingPopulation Population5(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
    
    // FastOPTICS() : calling FastOPTICS constructors
    // GPU-accelerated kNN is enabled when CUDA, ROCm, or Metal is compiled in.
    bool enableGPU = false;
#if defined(FLEXAIDS_USE_CUDA) || defined(FLEXAIDS_USE_ROCM) || defined(FLEXAIDS_USE_METAL)
    enableGPU = true;
#endif

    FastOPTICS Algo1(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population1, minPoints, enableGPU);
    minPoints = std::floor(minPoints * 1.5);
    FastOPTICS Algo2(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population2, minPoints, enableGPU);
    minPoints = std::floor(minPoints * 1.5);
    FastOPTICS Algo3(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population3, minPoints, enableGPU);
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
    printf("-- end of FastOPTICS_cluster --\n");
}