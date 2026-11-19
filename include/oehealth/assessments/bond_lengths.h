#ifndef OEHEALTH_ASSESSMENTS_BOND_LENGTHS_H
#define OEHEALTH_ASSESSMENTS_BOND_LENGTHS_H

#include "oehealth/assessment.h"

namespace OEHEALTH {

/// \brief Thresholds on |Z| for bond length outliers.
struct BondLengthOptions {
    double warning_z = 4.0;
    double error_z = 10.0;
};

/// \brief Compares heavy-atom bond lengths of standard amino acids with
/// Engh & Huber (2001) targets.
///
/// Z = (observed - ideal) / sigma for every intra-residue, peptide link,
/// and disulfide bond. |Z| >= warning_z gives a WARNING finding and
/// |Z| >= error_z an ERROR finding, code "bond_length_outlier".
class BondLengthAssessment : public Assessment {
public:
    /// \throws std::invalid_argument for NaN, infinite, or negative
    ///         thresholds, or error_z < warning_z.
    explicit BondLengthAssessment(const BondLengthOptions& options = BondLengthOptions());

    const BondLengthOptions& GetOptions() const;
    std::string GetName() const override;
    std::string GetDescription() const override;
    AssessmentResult Assess(const StructureContext& ctx) const override;

private:
    BondLengthOptions options_;
};

} // namespace OEHEALTH

#endif // OEHEALTH_ASSESSMENTS_BOND_LENGTHS_H
