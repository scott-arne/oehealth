#include <gtest/gtest.h>

#include "common.h"
#include "oehealth/assessment.h"
#include "oehealth/result.h"
#include "test_helpers.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace OEHEALTH::test {

namespace {

Finding finding_with(Severity severity) {
    Finding finding;
    finding.code = "test";
    finding.severity = severity;
    return finding;
}

} // namespace

TEST(WorstStatus, FollowsDeclarationOrder) {
    EXPECT_EQ(worst_status(Status::NOT_APPLICABLE, Status::PASS), Status::PASS);
    EXPECT_EQ(worst_status(Status::PASS, Status::WARN), Status::WARN);
    EXPECT_EQ(worst_status(Status::FAIL, Status::WARN), Status::FAIL);
    EXPECT_EQ(worst_status(Status::FAIL, Status::ERROR), Status::ERROR);
    EXPECT_EQ(worst_status(Status::ERROR, Status::NOT_APPLICABLE), Status::ERROR);
    EXPECT_EQ(worst_status(Status::PASS, Status::PASS), Status::PASS);
}

TEST(StatusFromFindings, DerivesFromSeverities) {
    EXPECT_EQ(status_from_findings({}), Status::PASS);
    EXPECT_EQ(status_from_findings({finding_with(Severity::INFO)}), Status::PASS);
    EXPECT_EQ(status_from_findings({finding_with(Severity::INFO), finding_with(Severity::WARNING)}),
              Status::WARN);
    EXPECT_EQ(
        status_from_findings({finding_with(Severity::WARNING), finding_with(Severity::ERROR)}),
        Status::FAIL);
}

TEST(MeetsThreshold, IsInclusiveWithinTolerance) {
    EXPECT_TRUE(meets_threshold(0.4, 0.4));
    EXPECT_TRUE(meets_threshold(1.7 + 1.7 - 3.0, 0.4));
    EXPECT_TRUE(meets_threshold(0.4 - 0.5 * THRESHOLD_TOLERANCE, 0.4));
    EXPECT_FALSE(meets_threshold(0.4 - 2.0 * THRESHOLD_TOLERANCE, 0.4));
    EXPECT_FALSE(meets_threshold(0.39, 0.4));
}

TEST(ValidateThresholds, RejectsInvalidPairs) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW(validate_thresholds("x", nan, 1.0), std::invalid_argument);
    EXPECT_THROW(validate_thresholds("x", 1.0, inf), std::invalid_argument);
    EXPECT_THROW(validate_thresholds("x", -1.0, 1.0), std::invalid_argument);
    EXPECT_THROW(validate_thresholds("x", 2.0, 1.0), std::invalid_argument);
    EXPECT_NO_THROW(validate_thresholds("x", 1.0, 1.0));
    EXPECT_NO_THROW(validate_thresholds("x", 0.0, 1.0));
}

TEST(Finding, NumericFieldsDefaultToNaN) {
    const Finding finding;
    EXPECT_TRUE(std::isnan(finding.observed));
    EXPECT_TRUE(std::isnan(finding.expected));
    EXPECT_TRUE(std::isnan(finding.deviation));
}

TEST(AtomRef, CopiesCrambinIdentifiers) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const AtomRef ref = make_atom_ref(find_atom(mol, 5, "CA"));
    EXPECT_EQ(ref.atom_name, "CA");
    EXPECT_EQ(ref.residue_name, "PRO");
    EXPECT_EQ(ref.residue_number, 5);
    EXPECT_EQ(ref.chain_id, "A");
    EXPECT_EQ(ref.insert_code, "");
    EXPECT_EQ(ref.atom_idx, find_atom(mol, 5, "CA")->GetIdx());
}

TEST(AtomPath, RepeatsResidueLabelWhenResidueChanges) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const std::vector<AtomRef> same = {make_atom_ref(find_atom(mol, 5, "CA")),
                                       make_atom_ref(find_atom(mol, 5, "CB"))};
    EXPECT_EQ(atom_path(same), "A PRO 5 CA-CB");
    const std::vector<AtomRef> across = {make_atom_ref(find_atom(mol, 4, "C")),
                                         make_atom_ref(find_atom(mol, 5, "N"))};
    EXPECT_EQ(atom_path(across), "A CYS 4 C-A PRO 5 N");
}

} // namespace OEHEALTH::test
