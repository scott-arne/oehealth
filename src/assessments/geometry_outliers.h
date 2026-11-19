#ifndef OEHEALTH_SRC_ASSESSMENTS_GEOMETRY_OUTLIERS_H
#define OEHEALTH_SRC_ASSESSMENTS_GEOMETRY_OUTLIERS_H

#include "geometry/targets.h"
#include "oehealth/result.h"

#include <oechem.h>

#include <string>
#include <vector>

namespace OEHEALTH {

/// \brief What differs between the bond length and bond angle assessments.
struct GeometryOutlierSettings {
    std::string assessment; ///< Name of the producing assessment.
    std::string code;       ///< Finding code.
    std::string unit;       ///< Unit used in messages: "A" or "deg".
    int decimals = 3;       ///< Decimals for observed and ideal values in messages.
    double warning_z = 4.0;
    double error_z = 10.0;
};

/// \brief Measure each term, compute Z, and build the result shared by
/// bond_lengths and bond_angles.
///
/// Two-atom terms are measured as distances and three-atom terms as angles.
AssessmentResult assess_geometry_terms(const OEChem::OEMolBase& mol,
                                       const std::vector<geometry::GeometryTerm>& terms,
                                       const GeometryOutlierSettings& settings);

} // namespace OEHEALTH

#endif // OEHEALTH_SRC_ASSESSMENTS_GEOMETRY_OUTLIERS_H
