#include "oehealth/context.h"

#include "common.h"

#include <map>
#include <memory>
#include <tuple>

namespace OEHEALTH {

const OEChem::OEAtomBase* ResidueInfo::FindAtom(const std::string& atom_name) const {
    const std::string wanted = trim(atom_name);
    for (const OEChem::OEAtomBase* atom : atoms) {
        if (trimmed_atom_name(atom) == wanted) {
            return atom;
        }
    }
    return nullptr;
}

StructureContext::StructureContext(const OEChem::OEMolBase& mol) : mol_(&mol) {
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->IsHydrogen()) {
            has_hydrogens_ = true;
            break;
        }
    }
}

const OEChem::OEMolBase& StructureContext::GetMolecule() const { return *mol_; }

bool StructureContext::HasHydrogens() const { return has_hydrogens_; }

const std::vector<ResidueInfo>& StructureContext::GetResidues() const {
    if (residues_) {
        return *residues_;
    }

    // Visit atoms in index order so residues are ordered by their first atom.
    std::vector<const OEChem::OEAtomBase*> by_index(mol_->GetMaxAtomIdx(), nullptr);
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol_->GetAtoms(); atom; ++atom) {
        by_index[atom->GetIdx()] = atom;
    }

    // Atoms of one residue are usually contiguous, but not always, so group by
    // OEChem residue identity. The key narrows the OESameResidue comparisons.
    using Key = std::tuple<std::string, int, char, std::string>;
    std::map<Key, std::vector<std::size_t>> groups_by_key;
    std::vector<OEChem::OEResidue> group_residues;
    auto residues = std::make_unique<std::vector<ResidueInfo>>();

    for (const OEChem::OEAtomBase* atom : by_index) {
        if (atom == nullptr) {
            continue;
        }
        const OEChem::OEResidue& residue = OEChem::OEAtomGetResidue(atom);
        const std::string chain_id = chain_id_of(residue);
        const Key key{chain_id, residue.GetResidueNumber(), residue.GetInsertCode(),
                      residue.GetName()};
        std::vector<std::size_t>& candidates = groups_by_key[key];
        std::size_t group = residues->size();
        for (std::size_t candidate : candidates) {
            if (OEChem::OESameResidue(group_residues[candidate], residue)) {
                group = candidate;
                break;
            }
        }
        if (group == residues->size()) {
            ResidueInfo info;
            info.name = trim(residue.GetName());
            info.chain_id = chain_id;
            info.number = residue.GetResidueNumber();
            info.insert_code = residue_char(residue.GetInsertCode());
            info.is_standard_amino_acid = is_standard_amino_acid(info.name);
            residues->push_back(info);
            group_residues.push_back(residue);
            candidates.push_back(group);
        }
        (*residues)[group].atoms.push_back(atom);
    }

    residues_ = std::move(residues);
    return *residues_;
}

} // namespace OEHEALTH
