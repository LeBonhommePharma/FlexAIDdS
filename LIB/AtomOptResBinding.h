#pragma once

#include "flexaid.h"

#include <span>
#include <unordered_map>
#include <vector>

namespace flexaids {

// Capture ownership once from the master arrays. Selective atom refresh leaves
// a mixture of master and workspace pointers, so rebinding must never infer an
// index by subtracting a retained workspace pointer from the master allocation.
class AtomOptResBinding {
public:
    AtomOptResBinding(std::span<const atom> master_atoms,
                      std::span<const OptRes> master_optres)
        : indices_(master_atoms.size(), no_optres), optres_count_(master_optres.size()) {
        std::unordered_map<const OptRes*, std::size_t> owned_indices;
        owned_indices.reserve(master_optres.size());
        for (std::size_t i = 0; i < master_optres.size(); ++i)
            owned_indices.emplace(&master_optres[i], i);

        for (std::size_t i = 0; i < master_atoms.size(); ++i) {
            if (master_atoms[i].optres == nullptr) continue;
            const auto found = owned_indices.find(master_atoms[i].optres);
            if (found == owned_indices.end())
                throw FlexAIDException("Atom OptRes pointer is not owned by the master OptRes array");
            indices_[i] = found->second;
        }
    }

    // Read-only mapping shared by all evaluation threads. The arrays supplied
    // here belong to one workspace; repeated calls are deliberately idempotent.
    void bind(std::span<atom> workspace_atoms, std::span<OptRes> workspace_optres) const {
        if (workspace_atoms.size() != indices_.size() ||
            workspace_optres.size() != optres_count_)
            throw FlexAIDException("Atom/OptRes workspace sizes changed after binding capture");
        for (std::size_t i = 0; i < indices_.size(); ++i) {
            const auto index = indices_[i];
            workspace_atoms[i].optres = index == no_optres ? nullptr : &workspace_optres[index];
        }
    }

private:
    static constexpr std::size_t no_optres = static_cast<std::size_t>(-1);
    std::vector<std::size_t> indices_;
    std::size_t optres_count_;
};

}  // namespace flexaids
