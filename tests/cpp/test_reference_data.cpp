#include <gtest/gtest.h>

#include "data/engh_huber_2001.h"

#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace OEHEALTH::test {

namespace {

using data::LinkVariant;

bool equal(const char* a, const char* b) { return std::strcmp(a, b) == 0; }

const data::ResidueBond* residue_bond(const char* residue, const char* a1, const char* a2,
                                      const char* variant) {
    for (const auto& row : data::residue_bonds()) {
        if (equal(row.residue, residue) && equal(row.atom1, a1) && equal(row.atom2, a2) &&
            equal(row.variant, variant)) {
            return &row;
        }
    }
    return nullptr;
}

const data::ResidueAngle* residue_angle(const char* residue, const char* a1, const char* a2,
                                        const char* a3, const char* variant) {
    for (const auto& row : data::residue_angles()) {
        if (equal(row.residue, residue) && equal(row.atom1, a1) && equal(row.atom2, a2) &&
            equal(row.atom3, a3) && equal(row.variant, variant)) {
            return &row;
        }
    }
    return nullptr;
}

int count_link_bonds(const char* a1, const char* a2, LinkVariant variant) {
    int count = 0;
    for (const auto& row : data::link_bonds()) {
        count += static_cast<int>(equal(row.atom1, a1) && equal(row.atom2, a2) &&
                                  row.link_variant == variant);
    }
    return count;
}

int count_link_angles(const char* a1, const char* a2, const char* a3, LinkVariant variant) {
    int count = 0;
    for (const auto& row : data::link_angles()) {
        count += static_cast<int>(equal(row.atom1, a1) && equal(row.atom2, a2) &&
                                  equal(row.atom3, a3) && row.link_variant == variant);
    }
    return count;
}

const data::LinkAngle* link_angle(const char* a1, const char* a2, const char* a3,
                                  LinkVariant variant) {
    for (const auto& row : data::link_angles()) {
        if (equal(row.atom1, a1) && equal(row.atom2, a2) && equal(row.atom3, a3) &&
            row.link_variant == variant) {
            return &row;
        }
    }
    return nullptr;
}

void expect_bond(const char* residue, const char* a1, const char* a2, double ideal, double sigma) {
    const data::ResidueBond* row = residue_bond(residue, a1, a2, ".");
    ASSERT_NE(row, nullptr) << residue << " " << a1 << "-" << a2;
    EXPECT_DOUBLE_EQ(row->ideal, ideal);
    EXPECT_DOUBLE_EQ(row->sigma, sigma);
}

void expect_angle(const char* residue, const char* a1, const char* a2, const char* a3,
                  const char* variant, double ideal, double sigma) {
    const data::ResidueAngle* row = residue_angle(residue, a1, a2, a3, variant);
    ASSERT_NE(row, nullptr) << residue << " " << a1 << "-" << a2 << "-" << a3;
    EXPECT_DOUBLE_EQ(row->ideal, ideal);
    EXPECT_DOUBLE_EQ(row->sigma, sigma);
}

// PDB heavy-atom names of each standard residue, OXT included.
const std::map<std::string, std::vector<std::string>>& heavy_atoms() {
    static const std::map<std::string, std::vector<std::string>> atoms = {
        {"ALA", {"N", "CA", "C", "O", "CB", "OXT"}},
        {"ARG", {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2", "OXT"}},
        {"ASN", {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2", "OXT"}},
        {"ASP", {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2", "OXT"}},
        {"CYS", {"N", "CA", "C", "O", "CB", "SG", "OXT"}},
        {"GLN", {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2", "OXT"}},
        {"GLU", {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2", "OXT"}},
        {"GLY", {"N", "CA", "C", "O", "OXT"}},
        {"HIS", {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2", "OXT"}},
        {"ILE", {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1", "OXT"}},
        {"LEU", {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "OXT"}},
        {"LYS", {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ", "OXT"}},
        {"MET", {"N", "CA", "C", "O", "CB", "CG", "SD", "CE", "OXT"}},
        {"PHE", {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OXT"}},
        {"PRO", {"N", "CA", "C", "O", "CB", "CG", "CD", "OXT"}},
        {"SER", {"N", "CA", "C", "O", "CB", "OG", "OXT"}},
        {"THR", {"N", "CA", "C", "O", "CB", "OG1", "CG2", "OXT"}},
        {"TRP",
         {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2",
          "OXT"}},
        {"TYR", {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "OXT"}},
        {"VAL", {"N", "CA", "C", "O", "CB", "CG1", "CG2", "OXT"}},
    };
    return atoms;
}

} // namespace

TEST(ReferenceData, SpotValuesMatchEnghHuber2001) {
    expect_bond("ALA", "N", "CA", 1.459, 0.020);
    expect_bond("ALA", "CA", "C", 1.525, 0.026);
    expect_bond("ALA", "C", "O", 1.229, 0.019);
    expect_bond("ALA", "CA", "CB", 1.520, 0.021);
    expect_angle("ALA", "N", "CA", "C", ".", 111.0, 2.7);
    expect_angle("ALA", "CA", "C", "O", ".", 120.1, 2.1);
    expect_angle("GLY", "N", "CA", "C", ".", 113.1, 2.5);
    expect_angle("PRO", "CA", "N", "CD", "trans", 111.7, 1.4);
    expect_angle("PRO", "CA", "N", "CD", "cis", 111.5, 1.4);

    ASSERT_EQ(count_link_bonds("C", "N", LinkVariant::DEFAULT), 1);
    for (const auto& row : data::link_bonds()) {
        if (equal(row.atom1, "C") && row.link_variant == LinkVariant::DEFAULT) {
            EXPECT_DOUBLE_EQ(row.ideal, 1.336);
            EXPECT_DOUBLE_EQ(row.sigma, 0.023);
        }
    }
    const data::LinkAngle* ca_c_n = link_angle("CA", "C", "N", LinkVariant::DEFAULT);
    ASSERT_NE(ca_c_n, nullptr);
    EXPECT_DOUBLE_EQ(ca_c_n->ideal, 117.2);
    EXPECT_DOUBLE_EQ(ca_c_n->sigma, 2.2);
    const data::LinkAngle* o_c_n = link_angle("O", "C", "N", LinkVariant::DEFAULT);
    ASSERT_NE(o_c_n, nullptr);
    EXPECT_DOUBLE_EQ(o_c_n->ideal, 122.7);
    EXPECT_DOUBLE_EQ(o_c_n->sigma, 1.6);
    const data::LinkAngle* c_n_ca = link_angle("C", "N", "CA", LinkVariant::DEFAULT);
    ASSERT_NE(c_n_ca, nullptr);
    EXPECT_DOUBLE_EQ(c_n_ca->ideal, 121.7);
    EXPECT_DOUBLE_EQ(c_n_ca->sigma, 2.5);
}

TEST(ReferenceData, ResidueBondsConnectExactlyTheHeavyAtoms) {
    for (const auto& [residue, names] : heavy_atoms()) {
        const std::set<std::string> expected(names.begin(), names.end());
        std::map<std::string, std::set<std::string>> graph;
        for (const auto& row : data::residue_bonds()) {
            if (residue == row.residue) {
                graph[row.atom1].insert(row.atom2);
                graph[row.atom2].insert(row.atom1);
            }
        }
        std::set<std::string> nodes;
        for (const auto& entry : graph) {
            nodes.insert(entry.first);
        }
        EXPECT_EQ(nodes, expected) << residue;

        // Walk the graph from N; every heavy atom must be reachable.
        std::set<std::string> reached{"N"};
        std::vector<std::string> frontier{"N"};
        while (!frontier.empty()) {
            const std::string atom = frontier.back();
            frontier.pop_back();
            for (const std::string& next : graph[atom]) {
                if (reached.insert(next).second) {
                    frontier.push_back(next);
                }
            }
        }
        EXPECT_EQ(reached, expected) << residue;
    }
}

TEST(ReferenceData, ResidueAnglesUseOnlyHeavyAtoms) {
    for (const auto& row : data::residue_angles()) {
        const auto found = heavy_atoms().find(row.residue);
        ASSERT_NE(found, heavy_atoms().end()) << row.residue;
        const std::set<std::string> names(found->second.begin(), found->second.end());
        EXPECT_TRUE(names.count(row.atom1) == 1 && names.count(row.atom2) == 1 &&
                    names.count(row.atom3) == 1)
            << row.residue << " " << row.atom1 << "-" << row.atom2 << "-" << row.atom3;
    }
}

TEST(ReferenceData, EveryPeptideTermHasOneRowPerVariant) {
    for (LinkVariant variant : {LinkVariant::DEFAULT, LinkVariant::GLY, LinkVariant::PRO}) {
        EXPECT_EQ(count_link_bonds("C", "N", variant), 1);
        EXPECT_EQ(count_link_angles("CA", "C", "N", variant), 1);
        EXPECT_EQ(count_link_angles("O", "C", "N", variant), 1);
    }
    for (LinkVariant variant :
         {LinkVariant::DEFAULT, LinkVariant::GLY, LinkVariant::PRO_CIS, LinkVariant::PRO_TRANS}) {
        EXPECT_EQ(count_link_angles("C", "N", "CA", variant), 1);
    }
    for (LinkVariant variant : {LinkVariant::PRO_CIS, LinkVariant::PRO_TRANS}) {
        EXPECT_EQ(count_link_angles("C", "N", "CD", variant), 1);
    }
    EXPECT_EQ(count_link_bonds("SG", "SG", LinkVariant::DISULFIDE), 1);
}

TEST(ReferenceData, ProlineHasNoIntraResidueCNCDAngle) {
    for (const auto& row : data::residue_angles()) {
        EXPECT_FALSE(equal(row.residue, "PRO") && equal(row.atom1, "C") && equal(row.atom2, "N") &&
                     equal(row.atom3, "CD"));
    }
}

TEST(ReferenceData, RecordsSourceVersion) {
    EXPECT_STREQ(data::source_version(), "Fri Feb 22 17:25:15 GMT 2013");
}

} // namespace OEHEALTH::test
