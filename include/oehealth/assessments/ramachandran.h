#ifndef OEHEALTH_ASSESSMENTS_RAMACHANDRAN_H
#define OEHEALTH_ASSESSMENTS_RAMACHANDRAN_H

#include "oehealth/assessment.h"

namespace OEHEALTH {

/// \brief Classifies standard amino-acid backbone torsions with OEBio's
/// Ramachandran analysis.
///
/// Residues without a CA atom, and residues OEBio cannot classify (termini
/// and chain breaks), are not assessed. Each Outlier is a WARNING finding,
/// code "ramachandran_outlier".
class RamachandranAssessment : public Assessment {
public:
    RamachandranAssessment();

    std::string GetName() const override;
    std::string GetDescription() const override;
    AssessmentResult Assess(const StructureContext& ctx) const override;
};

} // namespace OEHEALTH

#endif // OEHEALTH_ASSESSMENTS_RAMACHANDRAN_H
