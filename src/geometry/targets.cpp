#include "geometry/targets.h"

#include "geometry/measure.h"

#include <cmath>
#include <cstring>
#include <limits>
#include <optional>
#include <utility>

namespace OEHEALTH::geometry {

namespace {

using data::LinkVariant;

// SG atoms of two CYS residues this close, in angstroms, form a disulfide.
constexpr double DISULFIDE_CUTOFF = 2.5;

bool is_residue(const ResidueInfo& residue, const char* name) { return residue.name == name; }

bool row_applies(const char* row_variant, const std::string& variant) {
    return std::strcmp(row_variant, ".") == 0 || variant == row_variant;
}

bool has_disulfide(const StructureContext& ctx, std::size_t residue_index) {
    const auto& residues = ctx.GetResidues();
    const OEChem::OEAtomBase* sulfur = residues[residue_index].FindAtom("SG");
    if (sulfur == nullptr) {
        return false;
    }
    for (std::size_t other = 0; other < residues.size(); ++other) {
        if (other == residue_index || !is_residue(residues[other], "CYS")) {
            continue;
        }
        const OEChem::OEAtomBase* partner = residues[other].FindAtom("SG");
        if (partner != nullptr &&
            distance(ctx.GetMolecule(), sulfur, partner) <= DISULFIDE_CUTOFF) {
            return true;
        }
    }
    return false;
}

std::string his_variant(const StructureContext& ctx, std::size_t residue_index) {
    const ResidueInfo& residue = ctx.GetResidues()[residue_index];
    const std::string_view from_hydrogens = his_variant_from_hydrogens(
        residue.FindAtom("HD1") != nullptr, residue.FindAtom("HE2") != nullptr);
    if (!from_hydrogens.empty()) {
        return std::string(from_hydrogens);
    }
    return lowest_rmsz_his_variant(his_variant_rmsz(ctx, residue_index, "HISE"),
                                   his_variant_rmsz(ctx, residue_index, "HISD"),
                                   his_variant_rmsz(ctx, residue_index, "HISH"));
}

const data::LinkBond* find_link_bond(const char* atom1, const char* atom2, LinkVariant variant) {
    for (const data::LinkBond& row : data::link_bonds()) {
        if (std::strcmp(row.atom1, atom1) == 0 && std::strcmp(row.atom2, atom2) == 0 &&
            row.link_variant == variant) {
            return &row;
        }
    }
    return nullptr;
}

const data::LinkAngle* find_link_angle(const char* atom1, const char* atom2, const char* atom3,
                                       LinkVariant variant) {
    for (const data::LinkAngle& row : data::link_angles()) {
        if (std::strcmp(row.atom1, atom1) == 0 && std::strcmp(row.atom2, atom2) == 0 &&
            std::strcmp(row.atom3, atom3) == 0 && row.link_variant == variant) {
            return &row;
        }
    }
    return nullptr;
}

LinkVariant pro_link_variant(const StructureContext& ctx, std::size_t j) {
    return preceding_peptide_is_cis(ctx, j) ? LinkVariant::PRO_CIS : LinkVariant::PRO_TRANS;
}

void add_term(std::vector<GeometryTerm>& terms, std::vector<const OEChem::OEAtomBase*> atoms,
              double ideal, double sigma, std::string_view variant) {
    for (const OEChem::OEAtomBase* atom : atoms) {
        if (atom == nullptr) {
            return;
        }
    }
    terms.push_back(GeometryTerm{std::move(atoms), ideal, sigma, std::string(variant)});
}

void add_residue_bonds(const StructureContext& ctx, std::size_t i,
                       std::vector<GeometryTerm>& terms) {
    const ResidueInfo& residue = ctx.GetResidues()[i];
    const std::string variant = residue_variant(ctx, i);
    for (const data::ResidueBond& row : data::residue_bonds()) {
        if (residue.name == row.residue && row_applies(row.variant, variant)) {
            add_term(terms, {residue.FindAtom(row.atom1), residue.FindAtom(row.atom2)}, row.ideal,
                     row.sigma, row.variant);
        }
    }
}

void add_residue_angles(const StructureContext& ctx, std::size_t i,
                        std::vector<GeometryTerm>& terms) {
    const ResidueInfo& residue = ctx.GetResidues()[i];
    const std::string variant = residue_variant(ctx, i);
    for (const data::ResidueAngle& row : data::residue_angles()) {
        if (residue.name == row.residue && row_applies(row.variant, variant)) {
            add_term(terms,
                     {residue.FindAtom(row.atom1), residue.FindAtom(row.atom2),
                      residue.FindAtom(row.atom3)},
                     row.ideal, row.sigma, row.variant);
        }
    }
}

void add_link_bond(const StructureContext& ctx, std::size_t i, std::size_t j,
                   std::vector<GeometryTerm>& terms) {
    const ResidueInfo& first = ctx.GetResidues()[i];
    const ResidueInfo& second = ctx.GetResidues()[j];
    // Assumption from the design: the GLY C-N row is keyed on residue j.
    LinkVariant variant = LinkVariant::DEFAULT;
    if (is_residue(second, "PRO")) {
        variant = LinkVariant::PRO;
    } else if (is_residue(second, "GLY")) {
        variant = LinkVariant::GLY;
    }
    if (const data::LinkBond* row = find_link_bond("C", "N", variant)) {
        add_term(terms, {first.FindAtom("C"), second.FindAtom("N")}, row->ideal, row->sigma,
                 link_variant_name(variant));
    }
}

void add_link_angles(const StructureContext& ctx, std::size_t i, std::size_t j,
                     std::vector<GeometryTerm>& terms) {
    const ResidueInfo& first = ctx.GetResidues()[i];
    const ResidueInfo& second = ctx.GetResidues()[j];
    const bool second_is_pro = is_residue(second, "PRO");

    // CA(i)-C(i)-N(j) and O(i)-C(i)-N(j): the varying atom belongs to residue i.
    LinkVariant before = LinkVariant::DEFAULT;
    if (second_is_pro) {
        before = LinkVariant::PRO;
    } else if (is_residue(first, "GLY")) {
        before = LinkVariant::GLY;
    }
    for (const char* atom : {"CA", "O"}) {
        if (const data::LinkAngle* row = find_link_angle(atom, "C", "N", before)) {
            add_term(terms, {first.FindAtom(atom), first.FindAtom("C"), second.FindAtom("N")},
                     row->ideal, row->sigma, link_variant_name(before));
        }
    }

    // C(i)-N(j)-CA(j), and C(i)-N(j)-CD(j) for PRO: the varying atom belongs to residue j.
    LinkVariant after = LinkVariant::DEFAULT;
    if (second_is_pro) {
        after = pro_link_variant(ctx, j);
    } else if (is_residue(second, "GLY")) {
        after = LinkVariant::GLY;
    }
    if (const data::LinkAngle* row = find_link_angle("C", "N", "CA", after)) {
        add_term(terms, {first.FindAtom("C"), second.FindAtom("N"), second.FindAtom("CA")},
                 row->ideal, row->sigma, link_variant_name(after));
    }
    if (second_is_pro) {
        if (const data::LinkAngle* row = find_link_angle("C", "N", "CD", after)) {
            add_term(terms, {first.FindAtom("C"), second.FindAtom("N"), second.FindAtom("CD")},
                     row->ideal, row->sigma, link_variant_name(after));
        }
    }
}

void add_disulfide_bonds(const StructureContext& ctx, std::vector<GeometryTerm>& terms) {
    const auto& residues = ctx.GetResidues();
    const data::LinkBond* row = find_link_bond("SG", "SG", LinkVariant::DISULFIDE);
    if (row == nullptr) {
        return;
    }
    for (std::size_t a = 0; a < residues.size(); ++a) {
        if (!is_residue(residues[a], "CYS")) {
            continue;
        }
        const OEChem::OEAtomBase* first = residues[a].FindAtom("SG");
        for (std::size_t b = a + 1; b < residues.size() && first != nullptr; ++b) {
            if (!is_residue(residues[b], "CYS")) {
                continue;
            }
            const OEChem::OEAtomBase* second = residues[b].FindAtom("SG");
            if (second != nullptr &&
                distance(ctx.GetMolecule(), first, second) <= DISULFIDE_CUTOFF) {
                add_term(terms, {first, second}, row->ideal, row->sigma,
                         link_variant_name(LinkVariant::DISULFIDE));
            }
        }
    }
}

} // namespace

std::string residue_variant(const StructureContext& ctx, std::size_t residue_index) {
    const ResidueInfo& residue = ctx.GetResidues()[residue_index];
    if (is_residue(residue, "CYS")) {
        return has_disulfide(ctx, residue_index) ? "SS" : "SH";
    }
    if (is_residue(residue, "PRO")) {
        return preceding_peptide_is_cis(ctx, residue_index) ? "cis" : "trans";
    }
    if (is_residue(residue, "HIS")) {
        return his_variant(ctx, residue_index);
    }
    return ".";
}

bool preceding_peptide_is_cis(const StructureContext& ctx, std::size_t residue_index) {
    const auto predecessor = peptide_predecessor(ctx, residue_index);
    if (!predecessor) {
        return false;
    }
    const ResidueInfo& first = ctx.GetResidues()[*predecessor];
    const ResidueInfo& second = ctx.GetResidues()[residue_index];
    const OEChem::OEAtomBase* ca_i = first.FindAtom("CA");
    const OEChem::OEAtomBase* c_i = first.FindAtom("C");
    const OEChem::OEAtomBase* n_j = second.FindAtom("N");
    const OEChem::OEAtomBase* ca_j = second.FindAtom("CA");
    if (ca_i == nullptr || ca_j == nullptr) {
        return false;
    }
    const double omega = torsion_degrees(ctx.GetMolecule(), ca_i, c_i, n_j, ca_j);
    return std::fabs(omega) < 90.0;
}

std::string_view his_variant_from_hydrogens(bool has_hd1, bool has_he2) {
    if (has_hd1 && has_he2) {
        return "HISH";
    }
    if (has_hd1) {
        return "HISD";
    }
    if (has_he2) {
        return "HISE";
    }
    return "";
}

std::string lowest_rmsz_his_variant(double rmsz_hise, double rmsz_hisd, double rmsz_hish) {
    std::string best = "HISE";
    double best_rmsz = rmsz_hise;
    if (rmsz_hisd < best_rmsz) {
        best = "HISD";
        best_rmsz = rmsz_hisd;
    }
    if (rmsz_hish < best_rmsz) {
        best = "HISH";
    }
    return best;
}

double his_variant_rmsz(const StructureContext& ctx, std::size_t residue_index,
                        const std::string& variant) {
    const OEChem::OEMolBase& mol = ctx.GetMolecule();
    const ResidueInfo& residue = ctx.GetResidues()[residue_index];
    double sum_z2 = 0.0;
    std::size_t count = 0;
    for (const data::ResidueBond& row : data::residue_bonds()) {
        if (residue.name != row.residue || variant != row.variant) {
            continue;
        }
        const OEChem::OEAtomBase* a = residue.FindAtom(row.atom1);
        const OEChem::OEAtomBase* b = residue.FindAtom(row.atom2);
        if (a != nullptr && b != nullptr) {
            const double z = (distance(mol, a, b) - row.ideal) / row.sigma;
            sum_z2 += z * z;
            ++count;
        }
    }
    for (const data::ResidueAngle& row : data::residue_angles()) {
        if (residue.name != row.residue || variant != row.variant) {
            continue;
        }
        const OEChem::OEAtomBase* a = residue.FindAtom(row.atom1);
        const OEChem::OEAtomBase* b = residue.FindAtom(row.atom2);
        const OEChem::OEAtomBase* c = residue.FindAtom(row.atom3);
        if (a != nullptr && b != nullptr && c != nullptr) {
            const double z = (angle_degrees(mol, a, b, c) - row.ideal) / row.sigma;
            sum_z2 += z * z;
            ++count;
        }
    }
    if (count == 0) {
        return std::numeric_limits<double>::infinity();
    }
    return std::sqrt(sum_z2 / static_cast<double>(count));
}

std::string_view link_variant_name(data::LinkVariant variant) {
    switch (variant) {
        case LinkVariant::DEFAULT:
            return "DEFAULT";
        case LinkVariant::GLY:
            return "GLY";
        case LinkVariant::PRO:
            return "PRO";
        case LinkVariant::PRO_CIS:
            return "PRO_CIS";
        case LinkVariant::PRO_TRANS:
            return "PRO_TRANS";
        case LinkVariant::PRO_ALL:
            return "PRO_ALL";
        case LinkVariant::DISULFIDE:
            return "DISULFIDE";
    }
    return "UNKNOWN";
}

std::vector<GeometryTerm> collect_bond_terms(const StructureContext& ctx) {
    std::vector<GeometryTerm> terms;
    const auto& residues = ctx.GetResidues();
    for (std::size_t i = 0; i < residues.size(); ++i) {
        if (!residues[i].is_standard_amino_acid) {
            continue;
        }
        add_residue_bonds(ctx, i, terms);
        if (const auto j = peptide_successor(ctx, i)) {
            add_link_bond(ctx, i, *j, terms);
        }
    }
    add_disulfide_bonds(ctx, terms);
    return terms;
}

std::vector<GeometryTerm> collect_angle_terms(const StructureContext& ctx) {
    std::vector<GeometryTerm> terms;
    const auto& residues = ctx.GetResidues();
    for (std::size_t i = 0; i < residues.size(); ++i) {
        if (!residues[i].is_standard_amino_acid) {
            continue;
        }
        add_residue_angles(ctx, i, terms);
        if (const auto j = peptide_successor(ctx, i)) {
            add_link_angles(ctx, i, *j, terms);
        }
    }
    return terms;
}

} // namespace OEHEALTH::geometry
