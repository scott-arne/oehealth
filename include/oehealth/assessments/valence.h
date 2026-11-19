#ifndef OEHEALTH_ASSESSMENTS_VALENCE_H
#define OEHEALTH_ASSESSMENTS_VALENCE_H

#include "oehealth/assessment.h"

namespace OEHEALTH {

/// \brief Reports atoms whose valence OEChem considers invalid.
///
/// Runs OECheckAtomValence on every atom of an element in OEChem's valence
/// model (H, B, C, N, O, F, Si, P, S, Cl, Se, Br, I); other elements, including
/// metals and atomic number 0, are skipped because OEChem accepts them
/// unchecked and prints a warning. Each failure is an ERROR finding,
/// code "invalid_valence". Implicit hydrogen counts must be assigned, as
/// OpenEye's readers do by default: with zeroed counts, normal atoms are
/// reported as under-valent.
class ValenceAssessment : public Assessment {
public:
    ValenceAssessment();

    std::string GetName() const override;
    std::string GetDescription() const override;
    AssessmentResult Assess(const StructureContext& ctx) const override;
};

} // namespace OEHEALTH

#endif // OEHEALTH_ASSESSMENTS_VALENCE_H
