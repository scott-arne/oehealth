#ifndef OEHEALTH_ASSESSMENT_H
#define OEHEALTH_ASSESSMENT_H

#include "oehealth/context.h"
#include "oehealth/result.h"

#include <oechem.h>

#include <string>

namespace OEHEALTH {

/// \brief Tolerance for inclusive threshold comparisons.
///
/// A value meets a threshold t when value >= t - THRESHOLD_TOLERANCE, so that
/// floating-point rounding does not decide boundary cases (1.7 + 1.7 - 3.0
/// evaluates to 0.39999..., just under a 0.4 threshold).
constexpr double THRESHOLD_TOLERANCE = 1e-6;

/// \brief A single check, or a named set of checks, over one structure.
///
/// Assessments are immutable after construction, so one instance may be
/// shared between suites and threads.
class Assessment {
public:
    virtual ~Assessment() = default;

    /// \brief Stable identifier, e.g. "bond_lengths".
    virtual std::string GetName() const = 0;

    /// \brief One-sentence description.
    virtual std::string GetDescription() const = 0;

    /// \brief Assess the structure behind a context.
    virtual AssessmentResult Assess(const StructureContext& ctx) const = 0;

    /// \brief Validate a molecule, build a context, and assess it.
    ///
    /// A molecule with no atoms is passed through, and each assessment
    /// reports NOT_APPLICABLE. For a multi-conformer molecule the active
    /// conformer is assessed.
    ///
    /// \throws std::invalid_argument if the molecule has atoms but its
    ///         dimension is not 3.
    AssessmentResult Run(const OEChem::OEMolBase& mol) const;
};

} // namespace OEHEALTH

#endif // OEHEALTH_ASSESSMENT_H
