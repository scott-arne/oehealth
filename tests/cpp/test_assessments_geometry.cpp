#include <gtest/gtest.h>

#include "oehealth/assessments/bond_angles.h"
#include "oehealth/assessments/bond_lengths.h"
#include "test_helpers.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace OEHEALTH::test {

namespace {

// Ideal SER CA-CB-OG angle and its sigma in Engh & Huber (2001), degrees.
constexpr double SER_CA_CB_OG = 111.2;
constexpr double SER_CA_CB_OG_SIGMA = 2.7;
// Ideal ALA CA-CB bond length and its sigma, angstroms.
constexpr double ALA_CA_CB = 1.520;
constexpr double ALA_CA_CB_SIGMA = 0.021;
constexpr double PI = 3.14159265358979323846;

std::vector<const Finding*> findings_mentioning(const AssessmentResult& result,
                                                const OEChem::OEAtomBase* atom) {
    std::vector<const Finding*> found;
    for (const Finding& finding : result.findings) {
        if (mentions_atom(finding, atom)) {
            found.push_back(&finding);
        }
    }
    return found;
}

bool has_error_finding(const AssessmentResult& result) {
    for (const Finding& finding : result.findings) {
        if (finding.severity == Severity::ERROR) {
            return true;
        }
    }
    return false;
}

// Move `moving` along the anchor->moving direction to sit `length` from the anchor.
void set_bond_length(OEChem::OEMolBase& mol, const OEChem::OEAtomBase* anchor,
                     OEChem::OEAtomBase* moving, double length) {
    const auto a = coords(mol, anchor);
    const auto m = coords(mol, moving);
    const std::array<double, 3> v{m[0] - a[0], m[1] - a[1], m[2] - a[2]};
    const double norm = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    set_coords(
        mol, moving,
        {a[0] + v[0] / norm * length, a[1] + v[1] / norm * length, a[2] + v[2] / norm * length});
}

// Rotate `moving` about `vertex`, in the plane of the three atoms, so the
// end-vertex-moving angle becomes `degrees` while |vertex-moving| is kept.
void set_bond_angle(OEChem::OEMolBase& mol, const OEChem::OEAtomBase* end,
                    const OEChem::OEAtomBase* vertex, OEChem::OEAtomBase* moving, double degrees) {
    const auto e = coords(mol, end);
    const auto v = coords(mol, vertex);
    const auto m = coords(mol, moving);
    std::array<double, 3> u{e[0] - v[0], e[1] - v[1], e[2] - v[2]};
    const double u_norm = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    for (double& c : u) {
        c /= u_norm;
    }
    const std::array<double, 3> w{m[0] - v[0], m[1] - v[1], m[2] - v[2]};
    const double length = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    const double along = w[0] * u[0] + w[1] * u[1] + w[2] * u[2];
    std::array<double, 3> p{w[0] - along * u[0], w[1] - along * u[1], w[2] - along * u[2]};
    const double p_norm = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    for (double& c : p) {
        c /= p_norm;
    }
    const double theta = degrees * PI / 180.0;
    set_coords(mol, moving,
               {v[0] + length * (std::cos(theta) * u[0] + std::sin(theta) * p[0]),
                v[1] + length * (std::cos(theta) * u[1] + std::sin(theta) * p[1]),
                v[2] + length * (std::cos(theta) * u[2] + std::sin(theta) * p[2])});
}

// An ALA fragment with only CA and CB, stored in double precision, so the
// CA-CB distance is exactly the requested length.
void make_ala_fragment(OEChem::OEMol& mol, double length) {
    add_atom(mol, OEChem::OEElemNo::C, "CA", "ALA", 1, 0.0, 0.0, 0.0);
    add_atom(mol, OEChem::OEElemNo::C, "CB", "ALA", 1, length, 0.0, 0.0);
}

} // namespace

TEST(BondLengthAssessment, DefaultsAndName) {
    const BondLengthAssessment assessment;
    EXPECT_EQ(assessment.GetName(), "bond_lengths");
    EXPECT_FALSE(assessment.GetDescription().empty());
    EXPECT_DOUBLE_EQ(assessment.GetOptions().warning_z, 4.0);
    EXPECT_DOUBLE_EQ(assessment.GetOptions().error_z, 10.0);
}

TEST(BondAngleAssessment, DefaultsAndName) {
    const BondAngleAssessment assessment;
    EXPECT_EQ(assessment.GetName(), "bond_angles");
    EXPECT_FALSE(assessment.GetDescription().empty());
    EXPECT_DOUBLE_EQ(assessment.GetOptions().warning_z, 4.0);
    EXPECT_DOUBLE_EQ(assessment.GetOptions().error_z, 10.0);
}

TEST(GeometryOptions, RejectInvalidThresholds) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW(BondLengthAssessment(BondLengthOptions{nan, 10.0}), std::invalid_argument);
    EXPECT_THROW(BondLengthAssessment(BondLengthOptions{4.0, inf}), std::invalid_argument);
    EXPECT_THROW(BondLengthAssessment(BondLengthOptions{-1.0, 10.0}), std::invalid_argument);
    EXPECT_THROW(BondLengthAssessment(BondLengthOptions{5.0, 4.0}), std::invalid_argument);
    EXPECT_NO_THROW(BondLengthAssessment(BondLengthOptions{4.0, 4.0}));
    EXPECT_THROW(BondAngleAssessment(BondAngleOptions{nan, 10.0}), std::invalid_argument);
    EXPECT_THROW(BondAngleAssessment(BondAngleOptions{4.0, inf}), std::invalid_argument);
    EXPECT_THROW(BondAngleAssessment(BondAngleOptions{-1.0, 10.0}), std::invalid_argument);
    EXPECT_THROW(BondAngleAssessment(BondAngleOptions{5.0, 4.0}), std::invalid_argument);
    EXPECT_NO_THROW(BondAngleAssessment(BondAngleOptions{4.0, 4.0}));
}

TEST(BondLengthAssessment, CrambinHasNoErrors) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const AssessmentResult result = BondLengthAssessment().Run(mol);
    EXPECT_EQ(result.assessment, "bond_lengths");
    EXPECT_FALSE(has_error_finding(result));
    ASSERT_EQ(result.metrics.size(), 4U);
    EXPECT_EQ(result.metrics[0].name, "n_assessed");
    EXPECT_EQ(result.metrics[1].name, "n_outliers");
    EXPECT_EQ(result.metrics[2].name, "rmsz");
    EXPECT_EQ(result.metrics[3].name, "max_abs_z");
    EXPECT_GT(metric(result, "n_assessed"), 300.0);
    EXPECT_EQ(metric(result, "n_outliers"), static_cast<double>(result.findings.size()));
}

TEST(BondAngleAssessment, CrambinHasNoErrors) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const AssessmentResult result = BondAngleAssessment().Run(mol);
    EXPECT_EQ(result.assessment, "bond_angles");
    EXPECT_FALSE(has_error_finding(result));
    EXPECT_GT(metric(result, "n_assessed"), 400.0);
}

TEST(BondLengthAssessment, StretchedBondIsTheOnlyNewOutlier) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const OEChem::OEAtomBase* ca = find_atom(mol, 27, "CA");
    OEChem::OEAtomBase* cb = find_atom(mol, 27, "CB");
    const AssessmentResult before = BondLengthAssessment().Run(mol);
    ASSERT_TRUE(findings_mentioning(before, cb).empty()) << "pick another ALA";

    // CB is bonded only to CA among heavy atoms, so this changes one bond.
    set_bond_length(mol, ca, cb, ALA_CA_CB + 6.0 * ALA_CA_CB_SIGMA);
    const AssessmentResult after = BondLengthAssessment().Run(mol);
    const auto found = findings_mentioning(after, cb);
    ASSERT_EQ(found.size(), 1U);
    const Finding& finding = *found[0];
    EXPECT_EQ(finding.code, "bond_length_outlier");
    EXPECT_EQ(finding.severity, Severity::WARNING);
    ASSERT_EQ(finding.atoms.size(), 2U);
    EXPECT_EQ(finding.atoms[0].atom_name, "CA");
    EXPECT_EQ(finding.atoms[1].atom_name, "CB");
    EXPECT_NEAR(finding.deviation, 6.0, 1e-3);
    EXPECT_NEAR(finding.expected, ALA_CA_CB, 1e-12);
    EXPECT_EQ(finding.message.rfind("A ALA 27 CA-CB", 0), 0U) << finding.message;
    EXPECT_EQ(metric(after, "n_outliers"), metric(before, "n_outliers") + 1.0);
    EXPECT_EQ(after.status, Status::WARN);
}

TEST(BondLengthAssessment, BondExactlyAtWarningThresholdIsAWarning) {
    OEChem::OEMol mol(OEChem::OEMCMolType::DoubleCartesian);
    make_ala_fragment(mol, ALA_CA_CB + 4.0 * ALA_CA_CB_SIGMA);
    const AssessmentResult result = BondLengthAssessment().Run(mol);
    EXPECT_EQ(metric(result, "n_assessed"), 1.0);
    ASSERT_EQ(result.findings.size(), 1U);
    EXPECT_EQ(result.findings[0].severity, Severity::WARNING);
    EXPECT_NEAR(result.findings[0].deviation, 4.0, 1e-9);
}

TEST(BondLengthAssessment, ValueAtSharedThresholdIsAnError) {
    OEChem::OEMol mol(OEChem::OEMCMolType::DoubleCartesian);
    make_ala_fragment(mol, ALA_CA_CB + 4.0 * ALA_CA_CB_SIGMA);
    const AssessmentResult result = BondLengthAssessment(BondLengthOptions{4.0, 4.0}).Run(mol);
    ASSERT_EQ(result.findings.size(), 1U);
    EXPECT_EQ(result.findings[0].severity, Severity::ERROR);
    EXPECT_EQ(result.status, Status::FAIL);
}

TEST(BondAngleAssessment, DistortedAngleIsTheOnlyNewOutlier) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const OEChem::OEAtomBase* ca = find_atom(mol, 11, "CA");
    const OEChem::OEAtomBase* cb = find_atom(mol, 11, "CB");
    OEChem::OEAtomBase* og = find_atom(mol, 11, "OG");
    const AssessmentResult before = BondAngleAssessment().Run(mol);
    ASSERT_TRUE(findings_mentioning(before, og).empty()) << "pick another SER";

    // OG is bonded only to CB, so the only angle that changes is CA-CB-OG.
    set_bond_angle(mol, ca, cb, og, SER_CA_CB_OG + 6.0 * SER_CA_CB_OG_SIGMA);
    const AssessmentResult after = BondAngleAssessment().Run(mol);
    const auto found = findings_mentioning(after, og);
    ASSERT_EQ(found.size(), 1U);
    const Finding& finding = *found[0];
    EXPECT_EQ(finding.code, "bond_angle_outlier");
    EXPECT_EQ(finding.severity, Severity::WARNING);
    ASSERT_EQ(finding.atoms.size(), 3U);
    EXPECT_EQ(finding.atoms[0].atom_name, "CA");
    EXPECT_EQ(finding.atoms[1].atom_name, "CB");
    EXPECT_EQ(finding.atoms[2].atom_name, "OG");
    EXPECT_NEAR(finding.deviation, 6.0, 1e-3);
    EXPECT_EQ(metric(after, "n_outliers"), metric(before, "n_outliers") + 1.0);
}

TEST(GeometryAssessments, EmptyAndNonProteinMoleculesAreNotApplicable) {
    const OEChem::OEGraphMol empty;
    EXPECT_EQ(BondLengthAssessment().Run(empty).status, Status::NOT_APPLICABLE);
    EXPECT_EQ(BondAngleAssessment().Run(empty).status, Status::NOT_APPLICABLE);

    OEChem::OEGraphMol ethanol;
    OEChem::OESmilesToMol(ethanol, "CCO");
    ethanol.SetDimension(3);
    const AssessmentResult result = BondLengthAssessment().Run(ethanol);
    EXPECT_EQ(result.status, Status::NOT_APPLICABLE);
    EXPECT_FALSE(result.message.empty());
    EXPECT_TRUE(result.findings.empty());
    EXPECT_EQ(metric(result, "n_assessed"), 0.0);
}

} // namespace OEHEALTH::test
