#include <gtest/gtest.h>

#include "oehealth/assessments/clashes.h"
#include "test_helpers.h"

#include <limits>
#include <stdexcept>

namespace OEHEALTH::test {

namespace {

// Two unbonded atoms on the x axis; coordinates like these are exact in floats.
AssessmentResult assess_pair(unsigned int first, unsigned int second, double separation) {
    OEChem::OEGraphMol mol;
    add_atom(mol, first, "A1", "LIG", 1, 0.0, 0.0, 0.0);
    add_atom(mol, second, "A2", "LIG", 2, separation, 0.0, 0.0);
    return ClashAssessment().Run(mol);
}

} // namespace

TEST(ClashAssessment, DefaultsAndName) {
    const ClashAssessment assessment;
    EXPECT_EQ(assessment.GetName(), "clashes");
    EXPECT_FALSE(assessment.GetDescription().empty());
    EXPECT_DOUBLE_EQ(assessment.GetOptions().warning_overlap, 0.4);
    EXPECT_DOUBLE_EQ(assessment.GetOptions().error_overlap, 1.0);
    EXPECT_DOUBLE_EQ(assessment.GetOptions().polar_contact_distance, 2.9);
}

TEST(ClashAssessment, RejectsInvalidOptions) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW(ClashAssessment(ClashOptions{nan, 1.0, 2.9}), std::invalid_argument);
    EXPECT_THROW(ClashAssessment(ClashOptions{0.4, inf, 2.9}), std::invalid_argument);
    EXPECT_THROW(ClashAssessment(ClashOptions{-0.1, 1.0, 2.9}), std::invalid_argument);
    EXPECT_THROW(ClashAssessment(ClashOptions{1.0, 0.4, 2.9}), std::invalid_argument);
    EXPECT_THROW(ClashAssessment(ClashOptions{0.4, 1.0, 0.0}), std::invalid_argument);
    EXPECT_THROW(ClashAssessment(ClashOptions{0.4, 1.0, -2.9}), std::invalid_argument);
    EXPECT_THROW(ClashAssessment(ClashOptions{0.4, 1.0, nan}), std::invalid_argument);
    EXPECT_NO_THROW(ClashAssessment(ClashOptions{0.4, 0.4, 2.9}));
}

TEST(ClashAssessment, UnbondedCarbonsTwoAngstromsApartAreAnError) {
    const AssessmentResult result = assess_pair(OEChem::OEElemNo::C, OEChem::OEElemNo::C, 2.0);
    ASSERT_EQ(result.findings.size(), 1U);
    const Finding& clash = result.findings[0];
    EXPECT_EQ(clash.code, "clash");
    EXPECT_EQ(clash.severity, Severity::ERROR);
    EXPECT_EQ(clash.atoms.size(), 2U);
    EXPECT_NEAR(clash.observed, 2.0, 1e-9);
    EXPECT_NEAR(clash.expected, 3.4, 1e-9);
    EXPECT_NEAR(clash.deviation, 1.4, 1e-9);
    EXPECT_EQ(result.status, Status::FAIL);
    ASSERT_EQ(result.metrics.size(), 4U);
    EXPECT_EQ(result.metrics[0].name, "n_clashes");
    EXPECT_EQ(result.metrics[1].name, "clashscore");
    EXPECT_EQ(result.metrics[2].name, "max_overlap");
    EXPECT_EQ(result.metrics[3].name, "hydrogens_present");
    EXPECT_EQ(metric(result, "n_clashes"), 1.0);
    EXPECT_DOUBLE_EQ(metric(result, "clashscore"), 500.0);
    EXPECT_NEAR(metric(result, "max_overlap"), 1.4, 1e-9);
    EXPECT_EQ(metric(result, "hydrogens_present"), 0.0);
}

TEST(ClashAssessment, ExactBoundariesAreInclusive) {
    const auto at_warning = assess_pair(OEChem::OEElemNo::C, OEChem::OEElemNo::C, 3.0);
    ASSERT_EQ(at_warning.findings.size(), 1U);
    EXPECT_EQ(at_warning.findings[0].severity, Severity::WARNING);

    const auto at_error = assess_pair(OEChem::OEElemNo::C, OEChem::OEElemNo::C, 2.4);
    ASSERT_EQ(at_error.findings.size(), 1U);
    EXPECT_EQ(at_error.findings[0].severity, Severity::ERROR);

    const auto beyond = assess_pair(OEChem::OEElemNo::C, OEChem::OEElemNo::C, 3.01);
    EXPECT_TRUE(beyond.findings.empty());
    EXPECT_EQ(beyond.status, Status::PASS);
    EXPECT_EQ(metric(beyond, "max_overlap"), 0.0);
}

TEST(ClashAssessment, PolarPairsUsePolarContactDistance) {
    const auto hydrogen_bond = assess_pair(OEChem::OEElemNo::N, OEChem::OEElemNo::O, 2.8);
    EXPECT_TRUE(hydrogen_bond.findings.empty());

    const auto close = assess_pair(OEChem::OEElemNo::N, OEChem::OEElemNo::O, 2.4);
    ASSERT_EQ(close.findings.size(), 1U);
    EXPECT_EQ(close.findings[0].severity, Severity::WARNING);
    EXPECT_NEAR(close.findings[0].expected, 2.9, 1e-9);
    EXPECT_NEAR(close.findings[0].deviation, 0.5, 1e-6);
}

TEST(ClashAssessment, MetalsAreExcluded) {
    const auto result = assess_pair(OEChem::OEElemNo::Zn, OEChem::OEElemNo::N, 2.1);
    EXPECT_TRUE(result.findings.empty());
}

TEST(ClashAssessment, AtomsWithinThreeBondsAreExcluded) {
    OEChem::OEGraphMol mol;
    OEChem::OEAtomBase* c1 = add_atom(mol, OEChem::OEElemNo::C, "C1", "LIG", 1, 0.0, 0.0, 0.0);
    OEChem::OEAtomBase* c2 = add_atom(mol, OEChem::OEElemNo::C, "C2", "LIG", 1, 1.5, 0.0, 0.0);
    OEChem::OEAtomBase* c3 = add_atom(mol, OEChem::OEElemNo::C, "C3", "LIG", 1, 1.5, 1.5, 0.0);
    OEChem::OEAtomBase* c4 = add_atom(mol, OEChem::OEElemNo::C, "C4", "LIG", 1, 0.0, 2.0, 0.0);
    mol.NewBond(c1, c2, 1);
    mol.NewBond(c2, c3, 1);
    mol.NewBond(c3, c4, 1);
    // C1 and C4 are 2.0 A apart, which would be an ERROR clash if not 1-4.
    const AssessmentResult result = ClashAssessment().Run(mol);
    EXPECT_TRUE(result.findings.empty());
}

TEST(ClashAssessment, RecordsWhetherHydrogensArePresent) {
    OEChem::OEGraphMol mol;
    add_atom(mol, OEChem::OEElemNo::C, "C1", "LIG", 1, 0.0, 0.0, 0.0);
    add_atom(mol, OEChem::OEElemNo::H, "H1", "LIG", 1, 10.0, 0.0, 0.0);
    EXPECT_EQ(metric(ClashAssessment().Run(mol), "hydrogens_present"), 1.0);
}

TEST(ClashAssessment, CrambinHasNoErrorClashes) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const AssessmentResult result = ClashAssessment().Run(mol);
    for (const Finding& finding : result.findings) {
        EXPECT_NE(finding.severity, Severity::ERROR) << finding.message;
    }
    EXPECT_NE(result.status, Status::FAIL);
    EXPECT_EQ(metric(result, "hydrogens_present"), 0.0);
}

TEST(ClashAssessment, NothingInScopeIsNotApplicable) {
    const OEChem::OEGraphMol empty;
    const AssessmentResult none = ClashAssessment().Run(empty);
    EXPECT_EQ(none.status, Status::NOT_APPLICABLE);
    EXPECT_FALSE(none.message.empty());
    EXPECT_EQ(metric(none, "hydrogens_present"), 0.0);

    OEChem::OEGraphMol metals;
    add_atom(metals, OEChem::OEElemNo::Zn, "ZN", "ZN", 1, 0.0, 0.0, 0.0);
    EXPECT_EQ(ClashAssessment().Run(metals).status, Status::NOT_APPLICABLE);
}

} // namespace OEHEALTH::test
