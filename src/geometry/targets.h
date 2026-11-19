#ifndef OEHEALTH_SRC_GEOMETRY_TARGETS_H
#define OEHEALTH_SRC_GEOMETRY_TARGETS_H

#include "data/engh_huber_2001.h"
#include "oehealth/context.h"

#include <oechem.h>

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace OEHEALTH::geometry {

/// \brief One bond or angle to assess, with its Engh & Huber (2001) target.
struct GeometryTerm {
    /// Two atoms for a bond; three for an angle, with the vertex second.
    std::vector<const OEChem::OEAtomBase*> atoms;
    double ideal = 0.0; ///< Angstroms or degrees.
    double sigma = 0.0;
    /// The residue variant of the row (".", "SS", "SH", "cis", "trans",
    /// "HISD", "HISE", "HISH") or, for link terms, the link variant name
    /// ("DEFAULT", "GLY", "PRO", "PRO_CIS", "PRO_TRANS", "DISULFIDE").
    std::string variant;
};

/// \brief Every heavy-atom bond term in the structure: intra-residue,
/// peptide link, and disulfide. Terms with a missing atom are left out.
std::vector<GeometryTerm> collect_bond_terms(const StructureContext& ctx);

/// \brief Every heavy-atom angle term in the structure: intra-residue and
/// peptide link. Terms with a missing atom are left out.
std::vector<GeometryTerm> collect_angle_terms(const StructureContext& ctx);

/// \brief The variant used for a residue's intra-residue rows: "SS" or "SH"
/// for CYS, "cis" or "trans" for PRO, "HISD", "HISE", or "HISH" for HIS,
/// and "." for every other residue.
std::string residue_variant(const StructureContext& ctx, std::size_t residue_index);

/// \brief Whether the peptide bond preceding residue j is cis (|omega| < 90).
///
/// False when residue j has no peptide predecessor or a CA atom is missing.
bool preceding_peptide_is_cis(const StructureContext& ctx, std::size_t residue_index);

/// \brief The HIS variant implied by explicit hydrogens: "HISH" for both HD1
/// and HE2, "HISD" for HD1 only, "HISE" for HE2 only, "" for neither.
std::string_view his_variant_from_hydrogens(bool has_hd1, bool has_he2);

/// \brief The HIS variant with the lowest RMSZ; ties go to HISE, then HISD,
/// then HISH.
std::string lowest_rmsz_his_variant(double rmsz_hise, double rmsz_hisd, double rmsz_hish);

/// \brief RMSZ of one HIS residue against one variant's bond and angle rows.
///
/// Only rows with that exact variant are used, and only those whose atoms
/// are present. Returns infinity when no such row can be measured.
double his_variant_rmsz(const StructureContext& ctx, std::size_t residue_index,
                        const std::string& variant);

/// \brief The name of a link variant, e.g. "PRO_CIS".
std::string_view link_variant_name(data::LinkVariant variant);

} // namespace OEHEALTH::geometry

#endif // OEHEALTH_SRC_GEOMETRY_TARGETS_H
