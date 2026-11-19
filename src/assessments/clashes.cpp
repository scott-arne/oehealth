#include "oehealth/assessments/clashes.h"

#include "common.h"
#include "geometry/measure.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace OEHEALTH {

namespace {

// Margin added to the neighbor-search cutoff so pairs exactly at the warning
// threshold are candidates whether or not OENearestNbrs treats its cutoff as
// inclusive. The overlap test decides which candidates are findings.
constexpr double SEARCH_MARGIN = 0.01;

bool is_polar(const OEChem::OEAtomBase* atom) {
    const unsigned int element = atom->GetAtomicNum();
    return element == OEChem::OEElemNo::N || element == OEChem::OEElemNo::O;
}

// A hydrogen bonded to N or O. An unbonded hydrogen is not polar.
bool is_polar_hydrogen(const OEChem::OEAtomBase* atom) {
    if (atom->GetAtomicNum() != OEChem::OEElemNo::H) {
        return false;
    }
    for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom->GetAtoms(); nbr; ++nbr) {
        if (is_polar(nbr)) {
            return true;
        }
    }
    return false;
}

bool within_three_bonds(const OEChem::OEAtomBase* a, const OEChem::OEAtomBase* b) {
    for (OESystem::OEIter<OEChem::OEAtomBase> first = a->GetAtoms(); first; ++first) {
        if (first->GetIdx() == b->GetIdx()) {
            return true;
        }
        for (OESystem::OEIter<OEChem::OEAtomBase> second = first->GetAtoms(); second; ++second) {
            if (second->GetIdx() == a->GetIdx()) {
                continue;
            }
            if (second->GetIdx() == b->GetIdx()) {
                return true;
            }
            for (OESystem::OEIter<OEChem::OEAtomBase> third = second->GetAtoms(); third; ++third) {
                if (third->GetIdx() == b->GetIdx()) {
                    return true;
                }
            }
        }
    }
    return false;
}

struct Clash {
    const OEChem::OEAtomBase* first;
    const OEChem::OEAtomBase* second;
    double distance;
    double contact;
    double overlap;
    Severity severity;
};

} // namespace

ClashAssessment::ClashAssessment(const ClashOptions& options) : options_(options) {
    validate_thresholds("clash overlap", options_.warning_overlap, options_.error_overlap);
    if (!std::isfinite(options_.polar_contact_distance) || options_.polar_contact_distance <= 0.0) {
        throw std::invalid_argument("polar_contact_distance must be finite and positive");
    }
    if (!std::isfinite(options_.hbond_contact_distance) || options_.hbond_contact_distance <= 0.0) {
        throw std::invalid_argument("hbond_contact_distance must be finite and positive");
    }
}

const ClashOptions& ClashAssessment::GetOptions() const { return options_; }

std::string ClashAssessment::GetName() const { return "clashes"; }

std::string ClashAssessment::GetDescription() const {
    return "Finds non-bonded atom pairs that overlap by more than the clash threshold.";
}

AssessmentResult ClashAssessment::Assess(const StructureContext& ctx) const {
    const OEChem::OEMolBase& mol = ctx.GetMolecule();
    AssessmentResult result;
    result.assessment = GetName();
    const double hydrogens_present = ctx.HasHydrogens() ? 1.0 : 0.0;

    std::vector<double> radius(mol.GetMaxAtomIdx(), 0.0);
    std::vector<bool> in_scope(mol.GetMaxAtomIdx(), false);
    std::vector<bool> polar_hydrogen(mol.GetMaxAtomIdx(), false);
    std::size_t n_in_scope = 0;
    double max_radius = 0.0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        const unsigned int element = atom->GetAtomicNum();
        if (element == 0 || atom->IsMetal()) {
            continue;
        }
        in_scope[atom->GetIdx()] = true;
        polar_hydrogen[atom->GetIdx()] = is_polar_hydrogen(atom);
        radius[atom->GetIdx()] = OEChem::OEGetBondiVdWRadius(element);
        max_radius = std::max(max_radius, radius[atom->GetIdx()]);
        ++n_in_scope;
    }

    if (n_in_scope == 0) {
        result.status = Status::NOT_APPLICABLE;
        result.message = "no atoms in scope for clash detection";
        result.metrics = {{"hydrogens_present", hydrogens_present}};
        return result;
    }

    std::vector<Clash> clashes;
    const double max_contact = std::max(
        {2.0 * max_radius, options_.polar_contact_distance, options_.hbond_contact_distance});
    const double cutoff = max_contact - options_.warning_overlap + SEARCH_MARGIN;
    if (cutoff > 0.0) {
        for (OESystem::OEIter<OEChem::OENbrs> pair = OEChem::OEGetNearestNbrs(mol, cutoff); pair;
             ++pair) {
            const OEChem::OEAtomBase* a = pair->GetBgn();
            const OEChem::OEAtomBase* b = pair->GetEnd();
            if (!in_scope[a->GetIdx()] || !in_scope[b->GetIdx()] || within_three_bonds(a, b)) {
                continue;
            }
            if (a->GetIdx() > b->GetIdx()) {
                std::swap(a, b);
            }
            double contact = radius[a->GetIdx()] + radius[b->GetIdx()];
            if (is_polar(a) && is_polar(b)) {
                contact = options_.polar_contact_distance;
            } else if ((polar_hydrogen[a->GetIdx()] && is_polar(b)) ||
                       (polar_hydrogen[b->GetIdx()] && is_polar(a))) {
                contact = options_.hbond_contact_distance;
            }
            const double d = geometry::distance(mol, a, b);
            const double overlap = contact - d;
            if (meets_threshold(overlap, options_.error_overlap)) {
                clashes.push_back({a, b, d, contact, overlap, Severity::ERROR});
            } else if (meets_threshold(overlap, options_.warning_overlap)) {
                clashes.push_back({a, b, d, contact, overlap, Severity::WARNING});
            }
        }
    }

    std::sort(clashes.begin(), clashes.end(), [](const Clash& x, const Clash& y) {
        return std::make_pair(x.first->GetIdx(), x.second->GetIdx()) <
               std::make_pair(y.first->GetIdx(), y.second->GetIdx());
    });

    double max_overlap = 0.0;
    for (const Clash& clash : clashes) {
        max_overlap = std::max(max_overlap, clash.overlap);
        Finding finding;
        finding.code = "clash";
        finding.severity = clash.severity;
        finding.atoms = {make_atom_ref(clash.first), make_atom_ref(clash.second)};
        finding.observed = clash.distance;
        finding.expected = clash.contact;
        finding.deviation = clash.overlap;
        finding.message = atom_path(finding.atoms) + " " + format_fixed(clash.distance, 3) +
                          " A apart (contact " + format_fixed(clash.contact, 3) + ", overlap " +
                          format_fixed(clash.overlap, 3) + ")";
        result.findings.push_back(std::move(finding));
    }

    const auto n_clashes = static_cast<double>(result.findings.size());
    result.metrics = {
        {"n_clashes", n_clashes},
        {"clashscore", n_clashes * 1000.0 / static_cast<double>(n_in_scope)},
        {"max_overlap", max_overlap},
        {"hydrogens_present", hydrogens_present},
    };
    result.status = status_from_findings(result.findings);
    return result;
}

} // namespace OEHEALTH
