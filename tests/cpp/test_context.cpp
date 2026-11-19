#include <gtest/gtest.h>

#include "common.h"
#include "geometry/measure.h"
#include "oehealth/context.h"
#include "test_helpers.h"

#include <cmath>
#include <string>
#include <vector>

namespace OEHEALTH::test {

namespace {

constexpr double DEGREES_PER_RADIAN = 180.0 / 3.14159265358979323846;

std::size_t count_peptide_links(const StructureContext& ctx) {
    std::size_t links = 0;
    for (std::size_t i = 0; i < ctx.GetResidues().size(); ++i) {
        if (geometry::peptide_successor(ctx, i)) {
            ++links;
        }
    }
    return links;
}

void set_chain(OEChem::OEAtomBase* atom, const std::string& chain_id) {
    OEChem::OEResidue residue = OEChem::OEAtomGetResidue(atom);
    residue.SetExtChainID(chain_id);
    OEChem::OEAtomSetResidue(atom, residue);
}

} // namespace

TEST(StructureContext, CrambinHas46StandardResiduesInFileOrder) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const StructureContext ctx(mol);
    const auto& residues = ctx.GetResidues();
    ASSERT_EQ(residues.size(), 46U);
    for (std::size_t i = 0; i < residues.size(); ++i) {
        EXPECT_EQ(residues[i].number, static_cast<int>(i) + 1);
        EXPECT_TRUE(residues[i].is_standard_amino_acid) << residues[i].name;
        EXPECT_EQ(residues[i].chain_id, "A");
        EXPECT_EQ(residues[i].insert_code, "");
    }
    EXPECT_EQ(residues.front().name, "THR");
    EXPECT_EQ(residues.back().name, "ASN");
    std::size_t atom_count = 0;
    for (const ResidueInfo& residue : residues) {
        atom_count += residue.atoms.size();
    }
    EXPECT_EQ(atom_count, 327U);
}

TEST(StructureContext, FindAtomComparesTrimmedNames) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const StructureContext ctx(mol);
    const ResidueInfo& first = ctx.GetResidues().front();
    EXPECT_NE(first.FindAtom("CA"), nullptr);
    EXPECT_EQ(first.FindAtom(" CA "), first.FindAtom("CA"));
    EXPECT_EQ(first.FindAtom("XX"), nullptr);
}

TEST(StructureContext, ReportsExplicitHydrogens) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    EXPECT_FALSE(StructureContext(mol).HasHydrogens());
    add_atom(mol, OEChem::OEElemNo::H, "H", "HOH", 900, 50.0, 50.0, 50.0);
    EXPECT_TRUE(StructureContext(mol).HasHydrogens());
}

TEST(StructureContext, KeepsExtendedChainIds) {
    // mmCIF chain ids can be longer than one character. GetChainID() returns
    // only the first, which would make chains "AA" and "AB" look the same.
    OEChem::OEGraphMol mol;
    const std::vector<OEChem::OEAtomBase*> first = {
        add_atom(mol, OEChem::OEElemNo::N, "N", "ALA", 1, 0.0, 0.0, 0.0),
        add_atom(mol, OEChem::OEElemNo::C, "CA", "ALA", 1, 1.46, 0.0, 0.0),
        add_atom(mol, OEChem::OEElemNo::C, "C", "ALA", 1, 2.0, 1.4, 0.0)};
    const std::vector<OEChem::OEAtomBase*> second = {
        add_atom(mol, OEChem::OEElemNo::N, "N", "ALA", 2, 3.3, 1.4, 0.0),
        add_atom(mol, OEChem::OEElemNo::C, "CA", "ALA", 2, 4.0, 2.6, 0.0),
        add_atom(mol, OEChem::OEElemNo::C, "C", "ALA", 2, 5.5, 2.6, 0.0)};
    for (OEChem::OEAtomBase* atom : first) {
        set_chain(atom, "AA");
    }
    for (OEChem::OEAtomBase* atom : second) {
        set_chain(atom, "AB");
    }
    const StructureContext ctx(mol);
    ASSERT_EQ(ctx.GetResidues().size(), 2U);
    EXPECT_EQ(ctx.GetResidues()[0].chain_id, "AA");
    EXPECT_EQ(ctx.GetResidues()[1].chain_id, "AB");
    EXPECT_EQ(make_atom_ref(second[0]).chain_id, "AB");
    // C(1)-N(2) is 1.3 A, but the residues are in different chains.
    EXPECT_FALSE(geometry::peptide_successor(ctx, 0));
}

TEST(StructureContext, EmptyMoleculeHasNoResidues) {
    const OEChem::OEGraphMol mol;
    const StructureContext ctx(mol);
    EXPECT_TRUE(ctx.GetResidues().empty());
    EXPECT_FALSE(ctx.HasHydrogens());
}

TEST(PeptideLinks, CrambinHas45) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const StructureContext ctx(mol);
    EXPECT_EQ(count_peptide_links(ctx), 45U);
    EXPECT_EQ(geometry::peptide_successor(ctx, 0), std::optional<std::size_t>(1));
    EXPECT_FALSE(geometry::peptide_successor(ctx, 45));
    EXPECT_EQ(geometry::peptide_predecessor(ctx, 1), std::optional<std::size_t>(0));
    EXPECT_FALSE(geometry::peptide_predecessor(ctx, 0));
}

TEST(PeptideLinks, DeletingNitrogenBreaksTheLink) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    mol.DeleteAtom(find_atom(mol, 10, "N"));
    const StructureContext ctx(mol);
    // Residue 9 is at index 8 and residue 10 at index 9.
    EXPECT_FALSE(geometry::peptide_successor(ctx, 8));
    EXPECT_FALSE(geometry::peptide_predecessor(ctx, 9));
    EXPECT_EQ(count_peptide_links(ctx), 44U);
}

TEST(Measure, MatchesOpenEyeGeometry) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const OEChem::OEAtomBase* c9 = find_atom(mol, 9, "C");
    const OEChem::OEAtomBase* n10 = find_atom(mol, 10, "N");
    const OEChem::OEAtomBase* ca10 = find_atom(mol, 10, "CA");
    const OEChem::OEAtomBase* c10 = find_atom(mol, 10, "C");
    EXPECT_NEAR(geometry::distance(mol, n10, ca10), OEChem::OEGetDistance(mol, n10, ca10), 1e-5);
    EXPECT_NEAR(geometry::angle_degrees(mol, n10, ca10, c10),
                OEChem::OEGetAngle(mol, n10, ca10, c10) * DEGREES_PER_RADIAN, 1e-3);
    EXPECT_NEAR(geometry::torsion_degrees(mol, c9, n10, ca10, c10),
                OEChem::OEGetTorsion(mol, c9, n10, ca10, c10) * DEGREES_PER_RADIAN, 1e-3);
}

TEST(Measure, UsesDoublePrecisionCoordinates) {
    OEChem::OEMol mol(OEChem::OEMCMolType::DoubleCartesian);
    const double length = 1.520 + 4.0 * 0.021;
    const OEChem::OEAtomBase* a = add_atom(mol, OEChem::OEElemNo::C, "CA", "ALA", 1, 0.0, 0.0, 0.0);
    const OEChem::OEAtomBase* b =
        add_atom(mol, OEChem::OEElemNo::C, "CB", "ALA", 1, length, 0.0, 0.0);
    EXPECT_DOUBLE_EQ(geometry::distance(mol, a, b), length);
}

} // namespace OEHEALTH::test
