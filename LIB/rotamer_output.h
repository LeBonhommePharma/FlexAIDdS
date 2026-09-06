#pragma once

#include "flexaid.h"

#include <span>

namespace flexaids {

// Atom indices survive growth of the atom allocation. Resolve them here only
// after the rotamer has been built and accepted, never before a possible realloc.
inline void write_rotamer_model(FILE* output, std::span<const atom> atoms,
                                std::span<const resid> residues, int cb_index,
                                std::span<const int> built_indices, int model) {
    if (cb_index <= 0 || static_cast<std::size_t>(cb_index) >= atoms.size())
        throw FlexAIDException("Cannot write rotamer without a valid CB atom index");
    const atom& cb = atoms[cb_index];
    if (cb.ofres <= 0 || static_cast<std::size_t>(cb.ofres) >= residues.size())
        throw FlexAIDException("Rotamer CB atom has an invalid residue index");
    for (const int index : built_indices) {
        if (index <= 0 || static_cast<std::size_t>(index) >= atoms.size())
            throw FlexAIDException("Rotamer output contains an invalid built atom index");
    }
    const resid& res = residues[cb.ofres];
    fprintf(output, "MODEL       %2d\n", model);
    fprintf(output, "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f\n",
            cb.number, cb.name, res.name, res.chn, res.number,
            cb.coor[0], cb.coor[1], cb.coor[2]);
    for (const int index : built_indices) {
        const atom& at = atoms[index];
        fprintf(output, "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                at.number, at.name, res.name, res.chn, res.number,
                at.coor[0], at.coor[1], at.coor[2]);
    }
    fprintf(output, "ENDMDL\n");
}

}  // namespace flexaids
