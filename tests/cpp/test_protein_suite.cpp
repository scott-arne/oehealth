#include <gtest/gtest.h>

#include "oehealth/oehealth.h"
#include "test_helpers.h"

#include <string>
#include <vector>

namespace OEHEALTH::test {

TEST(Version, MacrosDefined) {
    EXPECT_GE(OEHEALTH_VERSION_MAJOR, 0);
    EXPECT_GE(OEHEALTH_VERSION_MINOR, 0);
    EXPECT_GE(OEHEALTH_VERSION_PATCH, 0);
}

TEST(ProteinGeometrySuite, HasTheFiveAssessmentsInOrder) {
    const auto suite = make_protein_geometry_suite();
    EXPECT_EQ(suite->GetName(), "protein_geometry");
    EXPECT_FALSE(suite->GetDescription().empty());
    std::vector<std::string> names;
    for (const auto& child : suite->GetChildren()) {
        names.push_back(child->GetName());
    }
    EXPECT_EQ(names, (std::vector<std::string>{"bond_lengths", "bond_angles", "clashes", "valence",
                                               "ramachandran"}));
}

TEST(ProteinGeometrySuite, CrambinHasNoErrors) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const AssessmentResult result = make_protein_geometry_suite()->Run(mol);
    ASSERT_EQ(result.children.size(), 5U);
    for (const AssessmentResult& child : result.children) {
        EXPECT_NE(child.status, Status::ERROR) << child.assessment << ": " << child.message;
        EXPECT_NE(child.status, Status::FAIL) << child.assessment;
        for (const Finding& finding : child.findings) {
            EXPECT_NE(finding.severity, Severity::ERROR) << finding.message;
        }
    }
}

TEST(ProteinGeometrySuite, EmptyMoleculeIsNotApplicableEverywhere) {
    const OEChem::OEGraphMol mol;
    const AssessmentResult result = make_protein_geometry_suite()->Run(mol);
    EXPECT_EQ(result.status, Status::NOT_APPLICABLE);
    ASSERT_EQ(result.children.size(), 5U);
    for (const AssessmentResult& child : result.children) {
        EXPECT_EQ(child.status, Status::NOT_APPLICABLE) << child.assessment;
    }
}

} // namespace OEHEALTH::test
