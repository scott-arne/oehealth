#include "oehealth/assessments/valence.h"

#include "common.h"

#include <utility>

namespace OEHEALTH {

namespace {

/// OEChem's valence model covers only these elements. OECheckAtomValence
/// accepts every other element unchecked (it returns true); all but atomic
/// number 0 also print an "Unknown element in TypeCheck" warning through
/// OEThrow, so skipping them changes no result and keeps metalloproteins from
/// filling stderr.
bool has_valence_model(unsigned int atomic_num) {
    using namespace OEChem::OEElemNo;
    switch (atomic_num) {
        case H:
        case B:
        case C:
        case N:
        case O:
        case F:
        case Si:
        case P:
        case S:
        case Cl:
        case Se:
        case Br:
        case I:
            return true;
        default:
            return false;
    }
}

} // namespace

ValenceAssessment::ValenceAssessment() = default;

std::string ValenceAssessment::GetName() const { return "valence"; }

std::string ValenceAssessment::GetDescription() const {
    return "Reports atoms whose valence is invalid for their element and formal charge.";
}

AssessmentResult ValenceAssessment::Assess(const StructureContext& ctx) const {
    const OEChem::OEMolBase& mol = ctx.GetMolecule();
    AssessmentResult result;
    result.assessment = GetName();
    if (mol.NumAtoms() == 0) {
        result.status = Status::NOT_APPLICABLE;
        result.message = "the molecule has no atoms";
        result.metrics = {{"n_invalid", 0.0}};
        return result;
    }

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (!has_valence_model(atom->GetAtomicNum()) || OEChem::OECheckAtomValence(atom)) {
            continue;
        }
        Finding finding;
        finding.code = "invalid_valence";
        finding.severity = Severity::ERROR;
        finding.atoms = {make_atom_ref(atom)};
        finding.observed = static_cast<double>(atom->GetValence());
        finding.details = {
            {"formal_charge", static_cast<double>(atom->GetFormalCharge())},
            {"implicit_h_count", static_cast<double>(atom->GetImplicitHCount())},
            {"degree", static_cast<double>(atom->GetDegree())},
        };
        finding.message = atom_path(finding.atoms) + " " +
                          OEChem::OEGetAtomicSymbol(atom->GetAtomicNum()) + " with formal charge " +
                          std::to_string(atom->GetFormalCharge()) + " has invalid valence " +
                          std::to_string(atom->GetValence());
        result.findings.push_back(std::move(finding));
    }

    result.metrics = {{"n_invalid", static_cast<double>(result.findings.size())}};
    result.status = status_from_findings(result.findings);
    return result;
}

} // namespace OEHEALTH
