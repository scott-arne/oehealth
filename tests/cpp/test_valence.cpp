#include <gtest/gtest.h>

#include "oehealth/assessments/valence.h"
#include "test_helpers.h"

namespace OEHEALTH::test {

namespace {

/// Redirects OEThrow to a string for the lifetime of the object.
class ScopedOEThrowCapture {
public:
    ScopedOEThrowCapture() { OESystem::OEThrow.SetOutputStream(stream_); }
    ~ScopedOEThrowCapture() { OESystem::OEThrow.SetOutputStream(OEPlatform::oeerr); }
    ScopedOEThrowCapture(const ScopedOEThrowCapture&) = delete;
    ScopedOEThrowCapture& operator=(const ScopedOEThrowCapture&) = delete;

    std::string Text() { return stream_.str(); }

private:
    OEPlatform::oeosstream stream_;
};

} // namespace

TEST(ValenceAssessment, Name) {
    const ValenceAssessment assessment;
    EXPECT_EQ(assessment.GetName(), "valence");
    EXPECT_FALSE(assessment.GetDescription().empty());
}

TEST(ValenceAssessment, PentavalentCarbonIsReported) {
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(OEChem::OESmilesToMol(mol, "C(C)(C)(C)(C)C"));
    const StructureContext ctx(mol);
    const AssessmentResult result = ValenceAssessment().Assess(ctx);
    ASSERT_EQ(result.findings.size(), 1U);
    const Finding& finding = result.findings[0];
    EXPECT_EQ(finding.code, "invalid_valence");
    EXPECT_EQ(finding.severity, Severity::ERROR);
    ASSERT_EQ(finding.atoms.size(), 1U);
    EXPECT_EQ(finding.observed, 5.0);
    ASSERT_EQ(finding.details.size(), 3U);
    EXPECT_EQ(finding.details[0].name, "formal_charge");
    EXPECT_EQ(finding.details[1].name, "implicit_h_count");
    EXPECT_EQ(finding.details[2].name, "degree");
    EXPECT_EQ(finding.details[2].value, 5.0);
    EXPECT_EQ(result.status, Status::FAIL);
    EXPECT_EQ(metric(result, "n_invalid"), 1.0);
}

TEST(ValenceAssessment, CrambinAsReadIsValid) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const AssessmentResult result = ValenceAssessment().Run(mol);
    EXPECT_TRUE(result.findings.empty());
    EXPECT_EQ(result.status, Status::PASS);
    EXPECT_EQ(metric(result, "n_invalid"), 0.0);
}

TEST(ValenceAssessment, ZeroedImplicitHydrogensAreReported) {
    // Documents the requirement that implicit hydrogen counts be assigned.
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        atom->SetImplicitHCount(0);
    }
    EXPECT_GT(ValenceAssessment().Run(mol).findings.size(), 0U);
}

TEST(ValenceAssessment, ElementsOutsideTheValenceModelAreSkippedSilently) {
    // OEChem accepts these unchecked but prints a warning for each.
    ScopedOEThrowCapture capture;
    OEChem::OEGraphMol mol;
    add_atom(mol, OEChem::OEElemNo::Zn, "ZN", "ZN", 1, 0.0, 0.0, 0.0);
    add_atom(mol, OEChem::OEElemNo::Xe, "XE", "XE", 2, 10.0, 0.0, 0.0);
    const AssessmentResult result = ValenceAssessment().Run(mol);
    EXPECT_TRUE(result.findings.empty());
    EXPECT_EQ(result.status, Status::PASS);
    EXPECT_EQ(capture.Text(), "");
}

TEST(ValenceAssessment, EmptyMoleculeIsNotApplicable) {
    const OEChem::OEGraphMol mol;
    const AssessmentResult result = ValenceAssessment().Run(mol);
    EXPECT_EQ(result.status, Status::NOT_APPLICABLE);
    EXPECT_FALSE(result.message.empty());
}

} // namespace OEHEALTH::test
