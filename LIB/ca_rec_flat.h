// ca_rec_flat.h — materialize ca_rec prev-chain in walk order (default OFF)
//
// FLEXAIDDS_CA_REC_FLAT=1: walk a packed index array filled in the same order
// as `curr = ca_index[i]; while (curr != -1) curr = ca_rec[curr].prev`.
// That is the live vcfunction order. Unset → callers keep the pointer chase.
//
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "EnvFlags.h"
#include "Vcontacts.h"

#include <vector>

namespace flexaids {

inline bool ca_rec_flat_enabled() noexcept
{
    return env_bool("FLEXAIDDS_CA_REC_FLAT", false);
}

/// Uncapped prev-chain length (legacy walk). `curr` follows prev to -1.
inline int ca_rec_chain_len(const int* ca_index, const ca_struct* ca_rec, int atom)
{
    if (!ca_index || !ca_rec || atom < 0) return 0;
    int n = 0;
    int curr = ca_index[atom];
    while (curr != -1) {
        ++n;
        curr = ca_rec[curr].prev;
    }
    return n;
}

/// Fill `out[0..)` with ca_rec indices in walk order. Returns count.
/// If the chain is longer than `cap`, returns `-(needed)` so callers cannot
/// silently treat a prefix as the full walk.
inline int flatten_ca_rec(const int* ca_index, const ca_struct* ca_rec,
                          int atom, int* out, int cap)
{
    if (!ca_index || !ca_rec || !out || cap <= 0 || atom < 0) return 0;
    const int need = ca_rec_chain_len(ca_index, ca_rec, atom);
    if (need > cap) return -need;
    int n = 0;
    int curr = ca_index[atom];
    while (curr != -1 && n < cap) {
        out[n++] = curr;
        curr = ca_rec[curr].prev;
    }
    return n;
}

/// Flatten the full chain; never truncates. Returns count.
inline int flatten_ca_rec_all(const int* ca_index, const ca_struct* ca_rec,
                              int atom, std::vector<int>& out)
{
    out.clear();
    const int need = ca_rec_chain_len(ca_index, ca_rec, atom);
    if (need <= 0) return 0;
    out.resize(static_cast<std::size_t>(need));
    const int n = flatten_ca_rec(ca_index, ca_rec, atom, out.data(), need);
    if (n < 0) {
        out.clear();
        return 0;
    }
    return n;
}

/// Next ca_rec index after `curr`. Used by every continue/advance in
/// vcfunction so a skip cannot leave flat_k stale and revisit a node.
/// When `use_flat`, `flat_k` is the next unread slot in `flat_idx`.
inline int ca_rec_next(bool use_flat, int curr, const ca_struct* ca_rec,
                       const int* flat_idx, int nflat, int& flat_k)
{
    if (use_flat) {
        if (!flat_idx || flat_k >= nflat) return -1;
        return flat_idx[flat_k++];
    }
    if (!ca_rec || curr < 0) return -1;
    return ca_rec[curr].prev;
}

}  // namespace flexaids
