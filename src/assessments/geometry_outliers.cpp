#include "assessments/geometry_outliers.h"

#include "common.h"
#include "geometry/measure.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace OEHEALTH {

AssessmentResult assess_geometry_terms(const OEChem::OEMolBase& mol,
                                       const std::vector<geometry::GeometryTerm>& terms,
                                       const GeometryOutlierSettings& settings) {
    AssessmentResult result;
    result.assessment = settings.assessment;
    if (terms.empty()) {
        result.status = Status::NOT_APPLICABLE;
        result.message = "no standard amino-acid terms to assess";
        result.metrics = {{"n_assessed", 0.0}, {"n_outliers", 0.0}};
        return result;
    }

    double sum_z2 = 0.0;
    double max_abs_z = 0.0;
    for (const geometry::GeometryTerm& term : terms) {
        const auto& atoms = term.atoms;
        const double observed = atoms.size() == 2
                                    ? geometry::distance(mol, atoms[0], atoms[1])
                                    : geometry::angle_degrees(mol, atoms[0], atoms[1], atoms[2]);
        const double z = (observed - term.ideal) / term.sigma;
        const double abs_z = std::fabs(z);
        sum_z2 += z * z;
        max_abs_z = std::max(max_abs_z, abs_z);

        Severity severity = Severity::INFO;
        if (meets_threshold(abs_z, settings.error_z)) {
            severity = Severity::ERROR;
        } else if (meets_threshold(abs_z, settings.warning_z)) {
            severity = Severity::WARNING;
        } else {
            continue;
        }

        Finding finding;
        finding.code = settings.code;
        finding.severity = severity;
        for (const OEChem::OEAtomBase* atom : atoms) {
            finding.atoms.push_back(make_atom_ref(atom));
        }
        finding.observed = observed;
        finding.expected = term.ideal;
        finding.deviation = z;
        finding.message = atom_path(finding.atoms) + " " +
                          format_fixed(observed, settings.decimals) + " " + settings.unit +
                          " (ideal " + format_fixed(term.ideal, settings.decimals) + ", Z " +
                          format_fixed(z, 1) + ")";
        result.findings.push_back(std::move(finding));
    }

    const auto n_assessed = static_cast<double>(terms.size());
    result.metrics = {
        {"n_assessed", n_assessed},
        {"n_outliers", static_cast<double>(result.findings.size())},
        {"rmsz", std::sqrt(sum_z2 / n_assessed)},
        {"max_abs_z", max_abs_z},
    };
    result.status = status_from_findings(result.findings);
    return result;
}

} // namespace OEHEALTH
