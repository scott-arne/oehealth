#include <gtest/gtest.h>

#include "geometry/measure.h"
#include "geometry/targets.h"
#include "test_helpers.h"

#include <array>
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

// Count three-atom terms named first-center-last, where the first two atoms belong to the same
// residue as each other and the last to the next residue (CA-C-N, O-C-N) or the reverse (C-N-CA).
std::size_t count_angle_family(const std::vector<geometry::GeometryTerm>& terms,
                               const std::string& first, const std::string& center,
                               const std::string& last) {
    std::size_t count = 0;
    for (const auto& term : terms) {
        if (term.atoms.size() == 3 && trimmed(term.atoms[0]->GetName()) == first &&
            trimmed(term.atoms[1]->GetName()) == center &&
            trimmed(term.atoms[2]->GetName()) == last) {
            ++count;
        }
    }
    return count;
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
    // RMSZ of the ideal coordinates against each variant's 14 ring rows, computed independently
    // in Python from the table: HISE 1.2989, HISD 1.2833, HISH 1.9774, so HISD wins narrowly.
    const std::string expected = "HISD";
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
    // Each link-angle family is emitted once per peptide link. The residue-internal CA-C-N-like
    // names cannot collide: no residue-internal angle runs CA-C-N, O-C-N or C-N-CA.
    EXPECT_EQ(count_angle_family(angles, "CA", "C", "N"), peptide_bonds);
    EXPECT_EQ(count_angle_family(angles, "O", "C", "N"), peptide_bonds);
    EXPECT_EQ(count_angle_family(angles, "C", "N", "CA"), peptide_bonds);

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

TEST(Targets, GlyProLinkAnglesFollowDangle) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    // Crambin has no Gly-Pro pair, so relabel residue 38 (ALA, after GLY 37) as PRO.
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        OEChem::OEResidue residue = OEChem::OEAtomGetResidue(atom);
        if (residue.GetResidueNumber() == 38) {
            residue.SetName("PRO");
            OEChem::OEAtomSetResidue(atom, residue);
        }
    }
    const StructureContext ctx(mol);
    const auto angles = geometry::collect_angle_terms(ctx);

    // CA-C-N prefers GLY in residue i over PRO in residue j.
    const auto* ca_c_n = find_term(
        angles, {find_atom(mol, 37, "CA"), find_atom(mol, 37, "C"), find_atom(mol, 38, "N")});
    ASSERT_NE(ca_c_n, nullptr);
    EXPECT_EQ(ca_c_n->variant, "GLY");
    EXPECT_DOUBLE_EQ(ca_c_n->ideal, 116.2);
    // O-C-N prefers PRO in residue j over GLY in residue i.
    const auto* o_c_n = find_term(
        angles, {find_atom(mol, 37, "O"), find_atom(mol, 37, "C"), find_atom(mol, 38, "N")});
    ASSERT_NE(o_c_n, nullptr);
    EXPECT_EQ(o_c_n->variant, "PRO");
    EXPECT_DOUBLE_EQ(o_c_n->ideal, 121.1);
}

TEST(Targets, VariantRowsBecomeTerms) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const StructureContext ctx(mol);
    const auto bonds = geometry::collect_bond_terms(ctx);
    const auto angles = geometry::collect_angle_terms(ctx);

    const auto* ss = find_term(bonds, {find_atom(mol, 3, "CB"), find_atom(mol, 3, "SG")});
    ASSERT_NE(ss, nullptr);
    EXPECT_EQ(ss->variant, "SS");
    EXPECT_DOUBLE_EQ(ss->ideal, 1.818);

    const auto* n_cd = find_term(bonds, {find_atom(mol, 5, "N"), find_atom(mol, 5, "CD")});
    ASSERT_NE(n_cd, nullptr);
    EXPECT_EQ(n_cd->variant, "trans");
    EXPECT_DOUBLE_EQ(n_cd->ideal, 1.474);

    const auto* c_n_cd = find_term(
        angles, {find_atom(mol, 4, "C"), find_atom(mol, 5, "N"), find_atom(mol, 5, "CD")});
    ASSERT_NE(c_n_cd, nullptr);
    EXPECT_EQ(c_n_cd->variant, "PRO_TRANS");
    EXPECT_DOUBLE_EQ(c_n_cd->ideal, 128.4);
    const auto* c_n_ca = find_term(
        angles, {find_atom(mol, 4, "C"), find_atom(mol, 5, "N"), find_atom(mol, 5, "CA")});
    ASSERT_NE(c_n_ca, nullptr);
    EXPECT_EQ(c_n_ca->variant, "PRO_TRANS");
    EXPECT_DOUBLE_EQ(c_n_ca->ideal, 119.3);
}

TEST(Targets, FreeCysteineBondTermIsSH) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    mol.DeleteAtom(find_atom(mol, 40, "SG"));
    const StructureContext ctx(mol);
    const auto bonds = geometry::collect_bond_terms(ctx);
    const auto* sh = find_term(bonds, {find_atom(mol, 3, "CB"), find_atom(mol, 3, "SG")});
    ASSERT_NE(sh, nullptr);
    EXPECT_EQ(sh->variant, "SH");
    EXPECT_DOUBLE_EQ(sh->ideal, 1.812);
}

TEST(Targets, DisulfideTermsFollowResidueOrderNotCoordinateOrder) {
    // Residues 1-4 carry only SG atoms. The x-sorted order is 2, 4, 1, 3, so the
    // pair (2, 4) is found before (1, 3) and must still come out second.
    OEChem::OEGraphMol mol;
    const std::array<double, 4> x_values = {10.0, 0.0, 11.5, 1.5};
    for (int i = 0; i < 4; ++i) {
        add_atom(mol, OEChem::OEElemNo::S, "SG", "CYS", i + 1, x_values[i], 0.0, 0.0);
    }
    const StructureContext ctx(mol);
    std::vector<std::array<int, 2>> pairs;
    for (const auto& term : geometry::collect_bond_terms(ctx)) {
        if (term.variant == "DISULFIDE") {
            pairs.push_back({OEChem::OEAtomGetResidue(term.atoms[0]).GetResidueNumber(),
                             OEChem::OEAtomGetResidue(term.atoms[1]).GetResidueNumber()});
        }
    }
    const std::vector<std::array<int, 2>> expected{{1, 3}, {2, 4}};
    EXPECT_EQ(pairs, expected);
}

TEST(Targets, DisulfideSweepSkipsNonFiniteCoordinates) {
    // A NaN x would break the strict weak ordering of the x-sorted sweep.
    OEChem::OEGraphMol mol;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    add_atom(mol, OEChem::OEElemNo::S, "SG", "CYS", 1, nan, 0.0, 0.0);
    add_atom(mol, OEChem::OEElemNo::S, "SG", "CYS", 2, 0.0, nan, 0.0);
    add_atom(mol, OEChem::OEElemNo::S, "SG", "CYS", 3, 20.0, 0.0, 0.0);
    add_atom(mol, OEChem::OEElemNo::S, "SG", "CYS", 4, 21.5, 0.0, 0.0);
    const StructureContext ctx(mol);
    EXPECT_EQ(geometry::residue_variant(ctx, 0), "SH");
    EXPECT_EQ(geometry::residue_variant(ctx, 1), "SH");
    EXPECT_EQ(geometry::residue_variant(ctx, 2), "SS");
    EXPECT_EQ(geometry::residue_variant(ctx, 3), "SS");
}

TEST(Measure, TorsionOfPlanarTransIsPositive180) {
    // This zigzag made atan2 return -180.0 (y == -0.0, x < 0) before the result was normalized.
    OEChem::OEGraphMol mol;
    const auto* a = add_atom(mol, OEChem::OEElemNo::C, "C1", "UNK", 1, 1.0, 1.0, 0.0);
    const auto* b = add_atom(mol, OEChem::OEElemNo::C, "C2", "UNK", 1, 0.0, 0.0, 0.0);
    const auto* c = add_atom(mol, OEChem::OEElemNo::C, "C3", "UNK", 1, 1.5, 0.0, 0.0);
    const auto* d = add_atom(mol, OEChem::OEElemNo::C, "C4", "UNK", 1, 2.5, -1.0, 0.0);
    EXPECT_EQ(geometry::torsion_degrees(mol, a, b, c, d), 180.0);
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
