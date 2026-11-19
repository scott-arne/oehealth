#ifndef OEHEALTH_ASSESSMENTS_CLASHES_H
#define OEHEALTH_ASSESSMENTS_CLASHES_H

#include "oehealth/assessment.h"

namespace OEHEALTH {

/// \brief Overlap thresholds and the polar contact distance, in angstroms.
struct ClashOptions {
    double warning_overlap = 0.4;
    double error_overlap = 1.0;
    double polar_contact_distance = 2.9; ///< Contact distance for N/O pairs.
};

/// \brief Finds pairs of atoms closer than their contact distance.
///
/// The contact distance is the sum of Bondi radii, or polar_contact_distance
/// when both atoms are N or O. Atoms within three bonds of each other,
/// atoms with atomic number 0, and metals are excluded. The atoms are
/// assessed as given; the result records whether hydrogens were present.
///
/// With explicit hydrogens, a polar hydrogen and its hydrogen-bond acceptor are
/// compared by Bondi radii (1.20 + 1.52 angstroms), not polar_contact_distance,
/// so ordinary hydrogen bonds are reported as clashes. Assess heavy-atom
/// structures when that matters.
class ClashAssessment : public Assessment {
public:
    /// \throws std::invalid_argument for NaN, infinite, or negative values,
    ///         a zero polar_contact_distance, or error_overlap < warning_overlap.
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
