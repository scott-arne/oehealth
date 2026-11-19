#include "oehealth/assessments/ramachandran.h"

#include "common.h"

#include <oebio.h>

#include <utility>

namespace OEHEALTH {

namespace {

constexpr double RADIANS_TO_DEGREES = 180.0 / 3.14159265358979323846;

} // namespace

RamachandranAssessment::RamachandranAssessment() = default;

std::string RamachandranAssessment::GetName() const { return "ramachandran"; }

std::string RamachandranAssessment::GetDescription() const {
    return "Classifies backbone phi/psi torsions of standard amino acids as favored, allowed, or "
           "outlier.";
}

AssessmentResult RamachandranAssessment::Assess(const StructureContext& ctx) const {
    AssessmentResult result;
    result.assessment = GetName();

    std::size_t n_favored = 0;
    std::size_t n_allowed = 0;
    std::size_t n_outliers = 0;
    for (const ResidueInfo& residue : ctx.GetResidues()) {
        if (!residue.is_standard_amino_acid) {
            continue;
        }
        const OEChem::OEAtomBase* alpha = residue.FindAtom("CA");
        if (alpha == nullptr) {
            continue;
        }
        // The (mol, residue) constructor gives the same results but scans the whole
        // molecule on every call, which made large structures quadratic.
        const OEBio::OERamachandranAnalysis analysis(alpha);
        const unsigned int category = analysis.GetRamaCategory();
        if (category == OEBio::OERamaCategory::Favored) {
            ++n_favored;
        } else if (category == OEBio::OERamaCategory::Allowed) {
            ++n_allowed;
        } else if (category == OEBio::OERamaCategory::Outlier) {
            ++n_outliers;
            const double phi = analysis.GetPhi() * RADIANS_TO_DEGREES;
            const double psi = analysis.GetPsi() * RADIANS_TO_DEGREES;
            Finding finding;
            finding.code = "ramachandran_outlier";
            finding.severity = Severity::WARNING;
            for (const char* name : {"N", "CA", "C"}) {
                if (const OEChem::OEAtomBase* atom = residue.FindAtom(name)) {
                    finding.atoms.push_back(make_atom_ref(atom));
                }
            }
            finding.observed = analysis.GetRamaScore();
            finding.details = {
                {"phi", phi},
                {"psi", psi},
                {"score", analysis.GetRamaScore()},
                {"rama_type", static_cast<double>(analysis.GetRamaType())},
            };
            finding.message = residue_label(make_atom_ref(alpha)) + " phi " + format_fixed(phi, 1) +
                              " psi " + format_fixed(psi, 1) + " (" +
                              OEBio::OEGetRamachandranTypeName(analysis.GetRamaType()) +
                              ") is a Ramachandran outlier";
            result.findings.push_back(std::move(finding));
        }
    }

    const std::size_t n_assessed = n_favored + n_allowed + n_outliers;
    if (n_assessed == 0) {
        result.status = Status::NOT_APPLICABLE;
        result.message = "no standard amino-acid residues could be classified";
        result.metrics = {{"n_assessed", 0.0}};
        return result;
    }
    const auto total = static_cast<double>(n_assessed);
    result.metrics = {
        {"n_assessed", total},
        {"n_favored", static_cast<double>(n_favored)},
        {"n_allowed", static_cast<double>(n_allowed)},
        {"n_outliers", static_cast<double>(n_outliers)},
        {"pct_favored", 100.0 * static_cast<double>(n_favored) / total},
        {"pct_outliers", 100.0 * static_cast<double>(n_outliers) / total},
    };
    result.status = status_from_findings(result.findings);
    return result;
}

} // namespace OEHEALTH
