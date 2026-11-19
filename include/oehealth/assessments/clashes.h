#ifndef OEHEALTH_ASSESSMENTS_CLASHES_H
#define OEHEALTH_ASSESSMENTS_CLASHES_H

#include "oehealth/assessment.h"

namespace OEHEALTH {

/// \brief Overlap thresholds and the polar contact distances, in angstroms.
struct ClashOptions {
    double warning_overlap = 0.4;
    double error_overlap = 1.0;
    double polar_contact_distance = 2.9; ///< Contact distance for N/O pairs.
    // 1.9 is the 2.9 heavy-atom rule less a 1.0 donor-H bond. It sits within
    // the limit Probe (the MolProbity clash engine) applies: donor-H/acceptor
    // pairs may overlap by 0.6 ("HBregular") on its own radii (polar H 1.05,
    // O 1.40, N 1.55), giving 1.85 for H...O and 2.00 for H...N.
    double hbond_contact_distance =
        1.9; ///< Contact distance for a polar hydrogen and an N/O acceptor.
};

/// \brief Finds pairs of atoms closer than their contact distance.
///
/// The contact distance is polar_contact_distance when both atoms are N or O,
/// hbond_contact_distance when one atom is a polar hydrogen and the other is
/// N or O, and the sum of Bondi radii otherwise. A polar hydrogen is a
/// hydrogen bonded to at least one N or O. Atoms within three bonds of each
/// other, atoms with atomic number 0, and metals are excluded. The atoms are
/// assessed as given; the result records whether hydrogens were present.
///
/// A polar hydrogen and its hydrogen-bond acceptor are therefore not reported
/// as a clash at ordinary hydrogen-bond distances. A hydrogen bonded only to
/// carbon, or not bonded at all, keeps its Bondi radius.
class ClashAssessment : public Assessment {
public:
    /// \throws std::invalid_argument for NaN, infinite, or negative values,
    ///         a zero polar_contact_distance or hbond_contact_distance, or
    ///         error_overlap < warning_overlap.
    explicit ClashAssessment(const ClashOptions& options = ClashOptions());

    const ClashOptions& GetOptions() const;
    std::string GetName() const override;
    std::string GetDescription() const override;
    AssessmentResult Assess(const StructureContext& ctx) const override;

private:
    ClashOptions options_;
};

} // namespace OEHEALTH

#endif // OEHEALTH_ASSESSMENTS_CLASHES_H
