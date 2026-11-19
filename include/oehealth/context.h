#ifndef OEHEALTH_CONTEXT_H
#define OEHEALTH_CONTEXT_H

#include <oechem.h>

#include <memory>
#include <string>
#include <vector>

namespace OEHEALTH {

/// \brief One residue of the structure, built from OEChem residue data.
struct ResidueInfo {
    std::string name; ///< Residue name, e.g. "ALA".
    std::string chain_id;
    int number = 0;
    std::string insert_code;                      ///< Empty when the residue has none.
    bool is_standard_amino_acid = false;          ///< One of the 20 standard residues.
    std::vector<const OEChem::OEAtomBase*> atoms; ///< In molecule order.

    /// \brief Find an atom by name, compared after trimming whitespace.
    /// \returns The atom, or nullptr when the residue has no such atom.
    const OEChem::OEAtomBase* FindAtom(const std::string& atom_name) const;
};

/// \brief Shared, lazily computed views of one molecule for assessments.
///
/// The context borrows the molecule: the caller keeps it alive and unchanged
/// for the context's lifetime. Lazily computed members are cached on first
/// use, so a single context must not be used from several threads at once.
class StructureContext {
public:
    explicit StructureContext(const OEChem::OEMolBase& mol);

    const OEChem::OEMolBase& GetMolecule() const;

    /// \brief Residues ordered by the index of their first atom (file order for PDB input).
    const std::vector<ResidueInfo>& GetResidues() const;

    /// \brief Whether the molecule contains explicit hydrogen atoms.
    bool HasHydrogens() const;

private:
    const OEChem::OEMolBase* mol_;
    bool has_hydrogens_ = false;
    mutable std::unique_ptr<std::vector<ResidueInfo>> residues_;
};

} // namespace OEHEALTH

#endif // OEHEALTH_CONTEXT_H
