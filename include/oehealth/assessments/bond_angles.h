#ifndef OEHEALTH_ASSESSMENTS_BOND_ANGLES_H
#define OEHEALTH_ASSESSMENTS_BOND_ANGLES_H

#include "oehealth/assessment.h"

namespace OEHEALTH {

/// \brief Thresholds on |Z| for bond angle outliers.
struct BondAngleOptions {
    double warning_z = 4.0;
    double error_z = 10.0;
};

/// \brief Compares heavy-atom bond angles of standard amino acids with
/// Engh & Huber (2001) targets.
///
/// Z = (observed - ideal) / sigma, in degrees, for every intra-residue and
/// peptide link angle. |Z| >= warning_z gives a WARNING finding and
/// |Z| >= error_z an ERROR finding, code "bond_angle_outlier".
class BondAngleAssessment : public Assessment {
public:
    /// \throws std::invalid_argument for NaN, infinite, or negative
    ///         thresholds, or error_z < warning_z.
    explicit BondAngleAssessment(const BondAngleOptions& options = BondAngleOptions());

    const BondAngleOptions& GetOptions() const;
    std::string GetName() const override;
    std::string GetDescription() const override;
    AssessmentResult Assess(const StructureContext& ctx) const override;

private:
    BondAngleOptions options_;
};

} // namespace OEHEALTH

#endif // OEHEALTH_ASSESSMENTS_BOND_ANGLES_H
