#include "common.h"

#include "oehealth/assessment.h"

#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace OEHEALTH {

std::string trim(const std::string& text) {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

bool is_standard_amino_acid(const std::string& residue_name) {
    static const std::array<const char*, 20> STANDARD = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
    for (const char* name : STANDARD) {
        if (residue_name == name) {
            return true;
        }
    }
    return false;
}

std::string trimmed_atom_name(const OEChem::OEAtomBase* atom) { return trim(atom->GetName()); }

std::string residue_char(char value) {
    if (value == ' ' || value == '\0') {
        return "";
    }
    return {value};
}

std::string chain_id_of(const OEChem::OEResidue& residue) { return trim(residue.GetExtChainID()); }

AtomRef make_atom_ref(const OEChem::OEAtomBase* atom) {
    const OEChem::OEResidue& residue = OEChem::OEAtomGetResidue(atom);
    AtomRef ref;
    ref.atom_idx = atom->GetIdx();
    ref.atom_name = trimmed_atom_name(atom);
    ref.residue_name = trim(residue.GetName());
    ref.chain_id = chain_id_of(residue);
    ref.residue_number = residue.GetResidueNumber();
    ref.insert_code = residue_char(residue.GetInsertCode());
    return ref;
}

std::string residue_label(const AtomRef& atom) {
    std::string label;
    if (!atom.chain_id.empty()) {
        label = atom.chain_id + " ";
    }
    label += atom.residue_name + " " + std::to_string(atom.residue_number) + atom.insert_code;
    return label;
}

namespace {

bool same_residue(const AtomRef& a, const AtomRef& b) {
    return a.chain_id == b.chain_id && a.residue_name == b.residue_name &&
           a.residue_number == b.residue_number && a.insert_code == b.insert_code;
}

} // namespace

std::string atom_path(const std::vector<AtomRef>& atoms) {
    std::string path;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        if (i > 0) {
            path += '-';
        }
        if (i == 0 || !same_residue(atoms[i - 1], atoms[i])) {
            path += residue_label(atoms[i]) + " ";
        }
        path += atoms[i].atom_name;
    }
    return path;
}

std::string format_fixed(double value, int decimals) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(decimals) << value;
    return stream.str();
}

Status status_from_findings(const std::vector<Finding>& findings) {
    Status status = Status::PASS;
    for (const Finding& finding : findings) {
        if (finding.severity == Severity::ERROR) {
            return Status::FAIL;
        }
        if (finding.severity == Severity::WARNING) {
            status = Status::WARN;
        }
    }
    return status;
}

bool meets_threshold(double value, double threshold) {
    return value >= threshold - THRESHOLD_TOLERANCE;
}

void validate_thresholds(const std::string& what, double warning, double error) {
    if (!std::isfinite(warning) || !std::isfinite(error)) {
        throw std::invalid_argument(what + " thresholds must be finite");
    }
    if (warning < 0.0 || error < 0.0) {
        throw std::invalid_argument(what + " thresholds must not be negative");
    }
    if (error < warning) {
        throw std::invalid_argument(what +
                                    " error threshold must not be less than its warning threshold");
    }
}

} // namespace OEHEALTH
