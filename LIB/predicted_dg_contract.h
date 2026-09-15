// predicted_dg_contract.h — DatasetRunner predicted_dG provenance flags
//
// CSV column `predicted_dG` is a historical name. It is NEVER experimental
// binding free energy ΔG_bind. When the Post-GA StatMech ledger is present the
// value is an ensemble free-energy *estimate* F = −kT ln Z; otherwise it is a
// CF/contact-function (or parsed-dG) fallback.
//
// Apache-2.0 © 2026 Le Bonhomme Pharma
#pragma once

namespace dataset {

struct PredictedDGContract {
    float predicted_dG{0.0f};
    bool  has_free_energy{false};
    bool  cf_fallback{true};
};

/// Assign predicted_dG plus machine-readable provenance.
/// has_free_energy is true only for Helmholtz F from the StatMech ledger.
/// cf_fallback is the complement: CF or parsed-dG stand-in, not experimental ΔG.
inline PredictedDGContract make_predicted_dg(bool have_free_energy,
                                             float free_energy_F,
                                             float best_dG,
                                             float best_cf) noexcept
{
    PredictedDGContract c;
    c.has_free_energy = have_free_energy;
    c.cf_fallback = !have_free_energy;
    c.predicted_dG = have_free_energy ? free_energy_F
                    : (best_dG != 0.0f) ? best_dG
                                        : best_cf;
    return c;
}

}  // namespace dataset
