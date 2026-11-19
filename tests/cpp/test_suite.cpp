#include <gtest/gtest.h>

#include "oehealth/assessment.h"
#include "oehealth/suite.h"
#include "test_helpers.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace OEHEALTH::test {

namespace {

// Returns a fixed status, so suite behavior can be tested in isolation.
class FixedAssessment : public Assessment {
public:
    FixedAssessment(std::string name, Status status) : name_(std::move(name)), status_(status) {}
    std::string GetName() const override { return name_; }
    std::string GetDescription() const override { return "Returns a fixed status."; }
    AssessmentResult Assess(const StructureContext& /*ctx*/) const override {
        AssessmentResult result;
        result.assessment = name_;
        result.status = status_;
        return result;
    }

private:
    std::string name_;
    Status status_;
};

class ThrowingAssessment : public Assessment {
public:
    std::string GetName() const override { return "thrower"; }
    std::string GetDescription() const override { return "Always throws."; }
    AssessmentResult Assess(const StructureContext& /*ctx*/) const override {
        throw std::runtime_error("boom");
    }
};

// Records the number of atoms it saw, to check what Run passes through.
class CountingAssessment : public Assessment {
public:
    std::string GetName() const override { return "counter"; }
    std::string GetDescription() const override { return "Counts atoms."; }
    AssessmentResult Assess(const StructureContext& ctx) const override {
        AssessmentResult result;
        result.assessment = GetName();
        result.status = ctx.GetMolecule().NumAtoms() == 0 ? Status::NOT_APPLICABLE : Status::PASS;
        result.metrics = {{"n_atoms", static_cast<double>(ctx.GetMolecule().NumAtoms())}};
        return result;
    }
};

AssessmentResult assess(const Assessment& assessment) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    return assessment.Run(mol);
}

} // namespace

TEST(AssessmentSuite, ReturnsConstructorValues) {
    auto child = std::make_shared<FixedAssessment>("a", Status::PASS);
    const AssessmentSuite suite("my_suite", "My description.", {child});
    EXPECT_EQ(suite.GetName(), "my_suite");
    EXPECT_EQ(suite.GetDescription(), "My description.");
    ASSERT_EQ(suite.GetChildren().size(), 1U);
    EXPECT_EQ(suite.GetChildren()[0], child);
}

TEST(AssessmentSuite, RunsChildrenInOrderAndTakesWorstStatus) {
    const AssessmentSuite suite("s", "",
                                {std::make_shared<FixedAssessment>("a", Status::PASS),
                                 std::make_shared<FixedAssessment>("b", Status::FAIL),
                                 std::make_shared<FixedAssessment>("c", Status::WARN)});
    const AssessmentResult result = assess(suite);
    EXPECT_EQ(result.assessment, "s");
    EXPECT_EQ(result.status, Status::FAIL);
    EXPECT_TRUE(result.message.empty());
    EXPECT_TRUE(result.metrics.empty());
    EXPECT_TRUE(result.findings.empty());
    ASSERT_EQ(result.children.size(), 3U);
    EXPECT_EQ(result.children[0].assessment, "a");
    EXPECT_EQ(result.children[1].assessment, "b");
    EXPECT_EQ(result.children[2].assessment, "c");
}

TEST(AssessmentSuite, RecordsThrowingChildAndRunsTheRest) {
    const AssessmentSuite suite("s", "",
                                {std::make_shared<ThrowingAssessment>(),
                                 std::make_shared<FixedAssessment>("after", Status::PASS)});
    const AssessmentResult result = assess(suite);
    ASSERT_EQ(result.children.size(), 2U);
    const AssessmentResult& failed = result.children[0];
    EXPECT_EQ(failed.assessment, "thrower");
    EXPECT_EQ(failed.status, Status::ERROR);
    EXPECT_EQ(failed.message, "boom");
    EXPECT_TRUE(failed.metrics.empty());
    EXPECT_TRUE(failed.findings.empty());
    EXPECT_TRUE(failed.children.empty());
    EXPECT_EQ(result.children[1].status, Status::PASS);
    EXPECT_EQ(result.status, Status::ERROR);
}

TEST(AssessmentSuite, EmptySuiteIsNotApplicable) {
    const AssessmentSuite suite("empty", "", {});
    const AssessmentResult result = assess(suite);
    EXPECT_EQ(result.status, Status::NOT_APPLICABLE);
    EXPECT_FALSE(result.message.empty());
    EXPECT_TRUE(result.children.empty());
}

TEST(AssessmentSuite, SuitesNest) {
    auto inner = std::make_shared<AssessmentSuite>(
        "inner", "",
        std::vector<std::shared_ptr<Assessment>>{
            std::make_shared<FixedAssessment>("x", Status::WARN)});
    const AssessmentSuite outer("outer", "", {inner});
    const AssessmentResult result = assess(outer);
    EXPECT_EQ(result.status, Status::WARN);
    ASSERT_EQ(result.children.size(), 1U);
    ASSERT_EQ(result.children[0].children.size(), 1U);
    EXPECT_EQ(result.children[0].children[0].assessment, "x");
}

TEST(AssessmentSuite, RejectsNullChildAndEmptyName) {
    EXPECT_THROW(AssessmentSuite("s", "", {nullptr}), std::invalid_argument);
    EXPECT_THROW(AssessmentSuite("", "", {}), std::invalid_argument);
}

TEST(AssessmentRun, Rejects2DMoleculeWithAtoms) {
    OEChem::OEGraphMol mol;
    OEChem::OESmilesToMol(mol, "CCO");
    ASSERT_NE(mol.GetDimension(), 3U);
    EXPECT_THROW(CountingAssessment().Run(mol), std::invalid_argument);
}

TEST(AssessmentRun, PassesEmptyMoleculeThrough) {
    const OEChem::OEGraphMol mol;
    const AssessmentResult result = CountingAssessment().Run(mol);
    EXPECT_EQ(result.status, Status::NOT_APPLICABLE);
    EXPECT_EQ(metric(result, "n_atoms"), 0.0);
}

} // namespace OEHEALTH::test
