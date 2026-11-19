#include <gtest/gtest.h>

#include "geometry/measure.h"
#include "geometry/targets.h"
#include "test_helpers.h"

#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace OEHEALTH::test {

namespace {

// Index of a residue number in a crambin context (residues are numbered 1-46 in order).
std::size_t index_of(int residue_number) { return static_cast<std::size_t>(residue_number - 1); }

bool term_has_atoms(const geometry::GeometryTerm& term,
                    const std::vector<const OEChem::OEAtomBase*>& atoms) {
    if (term.atoms.size() != atoms.size()) {
        return false;
    }
    for (std::size_t k = 0; k < atoms.size(); ++k) {
        if (atoms[k] == nullptr || term.atoms[k]->GetIdx() != atoms[k]->GetIdx()) {
            return false;
        }
    }
    return true;
}

const geometry::GeometryTerm* find_term(const std::vector<geometry::GeometryTerm>& terms,
                                        const std::vector<const OEChem::OEAtomBase*>& atoms) {
    for (const auto& term : terms) {
        if (term_has_atoms(term, atoms)) {
            return &term;
        }
    }
    return nullptr;
}

std::size_t count_c_n_cd_terms(const std::vector<geometry::GeometryTerm>& terms) {
    std::size_t count = 0;
    for (const auto& term : terms) {
        if (term.atoms.size() == 3 && trimmed(term.atoms[0]->GetName()) == "C" &&
            trimmed(term.atoms[1]->GetName()) == "N" && trimmed(term.atoms[2]->GetName()) == "CD") {
            ++count;
        }
    }
    return count;
}

// Heavy atoms of HIS from the ideal coordinates in standard_geometry.cif.
void add_ideal_histidine(OEChem::OEMolBase& mol) {
    struct IdealAtom {
        const char* name;
        unsigned int element;
        double x, y, z;
    };
    const std::vector<IdealAtom> atoms = {
        {"N", OEChem::OEElemNo::N, -0.040, -1.210, 0.053},
        {"CA", OEChem::OEElemNo::C, 1.172, -1.709, 0.652},
        {"C", OEChem::OEElemNo::C, 1.083, -3.207, 0.905},
        {"O", OEChem::OEElemNo::O, 0.040, -3.770, 1.222},
        {"CB", OEChem::OEElemNo::C, 1.484, -0.975, 1.962},
        {"CG", OEChem::OEElemNo::C, 2.940, -1.060, 2.353},
        {"ND1", OEChem::OEElemNo::N, 3.380, -2.075, 3.129},
        {"CD2", OEChem::OEElemNo::C, 3.960, -0.251, 2.046},
        {"CE1", OEChem::OEElemNo::C, 4.693, -1.908, 3.317},
        {"NE2", OEChem::OEElemNo::N, 5.058, -0.801, 2.662},
    };
    for (const IdealAtom& atom : atoms) {
        add_atom(mol, atom.element, atom.name, "HIS", 1, atom.x, atom.y, atom.z);
    }
}

} // namespace

TEST(Targets, CrambinCysteinesAreDisulfideBonded) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const StructureContext ctx(mol);
    for (int number : {3, 4, 16, 26, 32, 40}) {
        EXPECT_EQ(geometry::residue_variant(ctx, index_of(number)), "SS") << number;
    }
}

TEST(Targets, CysteineWithoutPartnerIsSH) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    mol.DeleteAtom(find_atom(mol, 40, "SG"));
    const StructureContext ctx(mol);
    EXPECT_EQ(geometry::residue_variant(ctx, index_of(3)), "SH");
    EXPECT_EQ(geometry::residue_variant(ctx, index_of(40)), "SH");
    EXPECT_EQ(geometry::residue_variant(ctx, index_of(4)), "SS");
}

TEST(Targets, CrambinProlinesAreTrans) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const StructureContext ctx(mol);
    for (int number : {5, 19, 22, 36, 41}) {
        EXPECT_FALSE(geometry::preceding_peptide_is_cis(ctx, index_of(number))) << number;
        EXPECT_EQ(geometry::residue_variant(ctx, index_of(number)), "trans") << number;
    }
}

TEST(Targets, ProlineIsCisWhenOmegaIsNearZero) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    // The C40-N41 bond is outside every disulfide ring, so it can be rotated.
    ASSERT_TRUE(OEChem::OESetTorsion(mol, find_atom(mol, 40, "CA"), find_atom(mol, 40, "C"),
                                     find_atom(mol, 41, "N"), find_atom(mol, 41, "CA"), 0.0));
    const StructureContext ctx(mol);
    EXPECT_TRUE(geometry::preceding_peptide_is_cis(ctx, index_of(41)));
    EXPECT_EQ(geometry::residue_variant(ctx, index_of(41)), "cis");

    const auto angles = geometry::collect_angle_terms(ctx);
    const geometry::GeometryTerm* c_n_cd = find_term(
        angles, {find_atom(mol, 40, "C"), find_atom(mol, 41, "N"), find_atom(mol, 41, "CD")});
    ASSERT_NE(c_n_cd, nullptr);
    EXPECT_EQ(c_n_cd->variant, "PRO_CIS");
    EXPECT_DOUBLE_EQ(c_n_cd->ideal, 120.6);
    const geometry::GeometryTerm* c_n_ca = find_term(
        angles, {find_atom(mol, 40, "C"), find_atom(mol, 41, "N"), find_atom(mol, 41, "CA")});
    ASSERT_NE(c_n_ca, nullptr);
    EXPECT_EQ(c_n_ca->variant, "PRO_CIS");
}

TEST(Targets, HistidineVariantFromHydrogens) {
    EXPECT_EQ(geometry::his_variant_from_hydrogens(true, true), "HISH");
    EXPECT_EQ(geometry::his_variant_from_hydrogens(true, false), "HISD");
    EXPECT_EQ(geometry::his_variant_from_hydrogens(false, true), "HISE");
    EXPECT_EQ(geometry::his_variant_from_hydrogens(false, false), "");
}

TEST(Targets, HistidineLowestRmszWithTieOrder) {
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_EQ(geometry::lowest_rmsz_his_variant(1.0, 2.0, 3.0), "HISE");
    EXPECT_EQ(geometry::lowest_rmsz_his_variant(2.0, 1.0, 3.0), "HISD");
    EXPECT_EQ(geometry::lowest_rmsz_his_variant(3.0, 2.0, 1.0), "HISH");
    EXPECT_EQ(geometry::lowest_rmsz_his_variant(1.0, 1.0, 1.0), "HISE");
    EXPECT_EQ(geometry::lowest_rmsz_his_variant(2.0, 1.0, 1.0), "HISD");
    EXPECT_EQ(geometry::lowest_rmsz_his_variant(inf, inf, inf), "HISE");
}

TEST(Targets, HistidineUsesExplicitHydrogens) {
    for (const auto& [hydrogens, expected] :
         std::vector<std::pair<std::vector<std::string>, std::string>>{
             {{"HD1"}, "HISD"}, {{"HE2"}, "HISE"}, {{"HD1", "HE2"}, "HISH"}}) {
        OEChem::OEGraphMol mol;
        add_ideal_histidine(mol);
        for (const std::string& name : hydrogens) {
            add_atom(mol, OEChem::OEElemNo::H, name, "HIS", 1, 9.0, 9.0, 9.0);
        }
        const StructureContext ctx(mol);
        EXPECT_EQ(geometry::residue_variant(ctx, 0), expected);
    }
}

TEST(Targets, HistidineWithoutHydrogensUsesBestFit) {
    OEChem::OEGraphMol mol;
    add_ideal_histidine(mol);
    const StructureContext ctx(mol);
    const std::string expected = geometry::lowest_rmsz_his_variant(
        geometry::his_variant_rmsz(ctx, 0, "HISE"), geometry::his_variant_rmsz(ctx, 0, "HISD"),
        geometry::his_variant_rmsz(ctx, 0, "HISH"));
    EXPECT_EQ(geometry::residue_variant(ctx, 0), expected);
    EXPECT_LT(geometry::his_variant_rmsz(ctx, 0, "HISE"), std::numeric_limits<double>::infinity());

    // The chosen variant's rows, and no other variant's, become terms.
    for (const auto& term : geometry::collect_bond_terms(ctx)) {
        EXPECT_TRUE(term.variant == "." || term.variant == expected) << term.variant;
    }
}

TEST(Targets, CrambinLinkTerms) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const StructureContext ctx(mol);
    const auto bonds = geometry::collect_bond_terms(ctx);
    const auto angles = geometry::collect_angle_terms(ctx);

    std::size_t peptide_bonds = 0;
    std::size_t disulfides = 0;
    for (const auto& term : bonds) {
        if (term.variant == "DISULFIDE") {
            ++disulfides;
        } else if (term.variant == "DEFAULT" || term.variant == "GLY" || term.variant == "PRO") {
            ++peptide_bonds;
        }
    }
    EXPECT_EQ(peptide_bonds, 45U);
    EXPECT_EQ(disulfides, 3U);
    EXPECT_EQ(count_c_n_cd_terms(angles), 5U);

    // C(4)-N(5): residue j is PRO.
    const auto* cys_pro = find_term(bonds, {find_atom(mol, 4, "C"), find_atom(mol, 5, "N")});
    ASSERT_NE(cys_pro, nullptr);
    EXPECT_EQ(cys_pro->variant, "PRO");
    EXPECT_DOUBLE_EQ(cys_pro->ideal, 1.338);
    // C(19)-N(20): residue j is GLY, so the C-N bond uses the GLY row.
    const auto* pro_gly = find_term(bonds, {find_atom(mol, 19, "C"), find_atom(mol, 20, "N")});
    ASSERT_NE(pro_gly, nullptr);
    EXPECT_EQ(pro_gly->variant, "GLY");
    // CA(20)-C(20)-N(21): residue i is GLY.
    const auto* gly_thr = find_term(
        angles, {find_atom(mol, 20, "CA"), find_atom(mol, 20, "C"), find_atom(mol, 21, "N")});
    ASSERT_NE(gly_thr, nullptr);
    EXPECT_EQ(gly_thr->variant, "GLY");
    EXPECT_DOUBLE_EQ(gly_thr->ideal, 116.2);
    // C(20)-N(21)-CA(21): residue j is THR, so DEFAULT.
    const auto* default_angle = find_term(
        angles, {find_atom(mol, 20, "C"), find_atom(mol, 21, "N"), find_atom(mol, 21, "CA")});
    ASSERT_NE(default_angle, nullptr);
    EXPECT_EQ(default_angle->variant, "DEFAULT");
}

TEST(Targets, ProlineWithoutPredecessorHasNoCNCDTerm) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    std::vector<OEChem::OEAtomBase*> residue4;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (OEChem::OEAtomGetResidue(atom).GetResidueNumber() == 4) {
            residue4.push_back(atom);
        }
    }
    for (OEChem::OEAtomBase* atom : residue4) {
        mol.DeleteAtom(atom);
    }
    const StructureContext ctx(mol);
    EXPECT_EQ(count_c_n_cd_terms(geometry::collect_angle_terms(ctx)), 4U);
}

} // namespace OEHEALTH::test
