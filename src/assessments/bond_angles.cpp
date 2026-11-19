#include "oehealth/assessments/bond_angles.h"

#include "assessments/geometry_outliers.h"
#include "common.h"
#include "geometry/targets.h"

namespace OEHEALTH {

BondAngleAssessment::BondAngleAssessment(const BondAngleOptions& options) : options_(options) {
    validate_thresholds("bond angle", options_.warning_z, options_.error_z);
}

const BondAngleOptions& BondAngleAssessment::GetOptions() const { return options_; }

std::string BondAngleAssessment::GetName() const { return "bond_angles"; }

std::string BondAngleAssessment::GetDescription() const {
    return "Compares heavy-atom bond angles of standard amino acids with Engh & Huber (2001) "
           "targets.";
}

AssessmentResult BondAngleAssessment::Assess(const StructureContext& ctx) const {
    GeometryOutlierSettings settings;
    settings.assessment = GetName();
    settings.code = "bond_angle_outlier";
    settings.unit = "deg";
    settings.decimals = 1;
    settings.warning_z = options_.warning_z;
    settings.error_z = options_.error_z;
    return assess_geometry_terms(ctx.GetMolecule(), geometry::collect_angle_terms(ctx), settings);
}

} // namespace OEHEALTH
