#pragma once

#include "gaboom.h"
#include <algorithm>
#include <numeric>
#include <vector>

namespace flexaids {

// The generation loop uses this same capture/restore path before and after
// reproduction. A cached score belongs to the complete chromosome, including
// its ring-pucker arrays; genes are the only externally owned chromosome data.
class GaEliteBuffer {
public:
    GaEliteBuffer(int count, int gene_count)
        : gene_count_(gene_count), elites_(count),
          genes_(static_cast<std::size_t>(count) * gene_count) {
        for (int i = 0; i < count; ++i)
            elites_[i].genes = genes_.data() + static_cast<std::size_t>(i) * gene_count;
    }
    GaEliteBuffer(const GaEliteBuffer&) = delete;
    GaEliteBuffer& operator=(const GaEliteBuffer&) = delete;

    void capture(const chromosome* population, int population_count) {
        if (elites_.empty()) return;
        std::vector<int> indices(population_count);
        std::iota(indices.begin(), indices.end(), 0);
        std::partial_sort(indices.begin(), indices.begin() + elites_.size(), indices.end(),
            [&](int a, int b) { return population[a].evalue < population[b].evalue; });
        for (std::size_t i = 0; i < elites_.size(); ++i)
            copy_chrom(&elites_[i], &population[indices[i]], gene_count_);
    }

    void restore(chromosome* population, int population_count) const {
        if (elites_.empty()) return;
        std::vector<int> indices(population_count);
        std::iota(indices.begin(), indices.end(), 0);
        std::partial_sort(indices.begin(), indices.begin() + elites_.size(), indices.end(),
            [&](int a, int b) { return population[a].evalue > population[b].evalue; });
        for (std::size_t i = 0; i < elites_.size(); ++i) {
            chromosome& destination = population[indices[i]];
            copy_chrom(&destination, &elites_[i], gene_count_);
            destination.fitnes = 0.0;
            destination.boltzmann_weight = 0.0;
            destination.free_energy = 0.0;
            // Preserve the captured validity flag: never certify an unscored
            // chromosome merely because it passed through elitism.
        }
    }

private:
    int gene_count_;
    std::vector<chromosome> elites_;
    std::vector<gene> genes_;
};
}  // namespace flexaids
