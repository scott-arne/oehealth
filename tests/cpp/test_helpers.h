#ifndef OEHEALTH_TESTS_TEST_HELPERS_H
#define OEHEALTH_TESTS_TEST_HELPERS_H

#include "oehealth/result.h"

#include <oechem.h>

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace OEHEALTH::test {

/// Read tests/data/1crn.pdb (crambin) with OpenEye's default PDB reader.
inline void read_crambin(OEChem::OEMolBase& mol) {
    const std::string path = std::string(OEHEALTH_TEST_DATA_DIR) + "/1crn.pdb";
    OEChem::oemolistream stream;
    if (!stream.open(path) || !OEChem::OEReadMolecule(stream, mol)) {
        throw std::runtime_error("cannot read " + path);
    }
}

inline std::string trimmed(const std::string& text) {
    const auto first = text.find_first_not_of(' ');
    if (first == std::string::npos) {
        return "";
    }
    return text.substr(first, text.find_last_not_of(' ') - first + 1);
}

/// Find an atom by residue number and trimmed atom name, or return nullptr.
inline OEChem::OEAtomBase* find_atom(const OEChem::OEMolBase& mol, int residue_number,
                                     const std::string& name) {
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (OEChem::OEAtomGetResidue(atom).GetResidueNumber() == residue_number &&
            trimmed(atom->GetName()) == name) {
            return atom;
        }
    }
    return nullptr;
}

/// Add an atom with PDB-style residue information and coordinates.
inline OEChem::OEAtomBase* add_atom(OEChem::OEMolBase& mol, unsigned int element,
                                    const std::string& name, const std::string& residue_name,
                                    int residue_number, double x, double y, double z) {
    OEChem::OEAtomBase* atom = mol.NewAtom(element);
    atom->SetName(name.c_str());
    OEChem::OEResidue residue;
    residue.SetName(residue_name.c_str());
    residue.SetChainID('A');
    residue.SetResidueNumber(residue_number);
    OEChem::OEAtomSetResidue(atom, residue);
    const std::array<double, 3> xyz{x, y, z};
    mol.SetCoords(atom, xyz.data());
    mol.SetDimension(3);
    return atom;
}

inline std::array<double, 3> coords(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* atom) {
    std::array<double, 3> xyz{0.0, 0.0, 0.0};
    mol.GetCoords(atom, xyz.data());
    return xyz;
}

inline void set_coords(OEChem::OEMolBase& mol, OEChem::OEAtomBase* atom,
                       const std::array<double, 3>& xyz) {
    mol.SetCoords(atom, xyz.data());
}

/// Return a metric's value; throws if the result has no such metric.
inline double metric(const AssessmentResult& result, const std::string& name) {
    for (const Metric& m : result.metrics) {
        if (m.name == name) {
            return m.value;
        }
    }
    throw std::runtime_error("no metric named " + name);
}

/// Return whether any of a finding's atoms has the given atom index.
inline bool mentions_atom(const Finding& finding, const OEChem::OEAtomBase* atom) {
    for (const AtomRef& ref : finding.atoms) {
        if (ref.atom_idx == atom->GetIdx()) {
            return true;
        }
    }
    return false;
}

} // namespace OEHEALTH::test

#endif // OEHEALTH_TESTS_TEST_HELPERS_H
