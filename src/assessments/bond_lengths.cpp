#include "oehealth/assessments/bond_lengths.h"

#include "assessments/geometry_outliers.h"
#include "common.h"
#include "geometry/targets.h"

namespace OEHEALTH {

BondLengthAssessment::BondLengthAssessment(const BondLengthOptions& options) : options_(options) {
    validate_thresholds("bond length", options_.warning_z, options_.error_z);
}

const BondLengthOptions& BondLengthAssessment::GetOptions() const { return options_; }

std::string BondLengthAssessment::GetName() const { return "bond_lengths"; }

std::string BondLengthAssessment::GetDescription() const {
    return "Compares heavy-atom bond lengths of standard amino acids with Engh & Huber (2001) "
           "targets.";
}

AssessmentResult BondLengthAssessment::Assess(const StructureContext& ctx) const {
    GeometryOutlierSettings settings;
    settings.assessment = GetName();
    settings.code = "bond_length_outlier";
    settings.unit = "A";
    settings.decimals = 3;
    settings.warning_z = options_.warning_z;
    settings.error_z = options_.error_z;
    return assess_geometry_terms(ctx.GetMolecule(), geometry::collect_bond_terms(ctx), settings);
}

} // namespace OEHEALTH
