#include <gtest/gtest.h>

#include "oehealth/assessments/bond_lengths.h"
#include "oehealth/json.h"
#include "oehealth/oehealth.h"
#include "test_helpers.h"

#include <nlohmann/json.hpp>

#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace OEHEALTH::test {

namespace {

using Json = nlohmann::ordered_json;

AssessmentResult sample_result() {
    AtomRef ca;
    ca.atom_idx = 10;
    ca.atom_name = "CA";
    ca.residue_name = "ALA";
    ca.chain_id = "A";
    ca.residue_number = 5;

    Finding finding;
    finding.code = "bond_length_outlier";
    finding.severity = Severity::WARNING;
    finding.message = "A ALA 5 CA-CB";
    finding.atoms = {ca};
    finding.observed = 1.617;
    finding.expected = 1.520;
    finding.deviation = std::numeric_limits<double>::quiet_NaN();
    finding.details = {{"phi", -60.0}, {"psi", -45.0}};

    AssessmentResult child;
    child.assessment = "child";
    child.status = Status::PASS;

    AssessmentResult result;
    result.assessment = "bond_lengths";
    result.status = Status::WARN;
    result.metrics = {{"n_assessed", 412.0}, {"n_outliers", 1.0}, {"rmsz", 0.83}};
    result.findings = {finding};
    result.children = {child};
    return result;
}

std::vector<std::string> keys(const Json& object) {
    std::vector<std::string> names;
    for (const auto& item : object.items()) {
        names.push_back(item.key());
    }
    return names;
}

} // namespace

TEST(ToJson, RootLayoutAndOrder) {
    const Json root = Json::parse(to_json(sample_result()));
    EXPECT_EQ(keys(root),
              (std::vector<std::string>{"schema_version", "oehealth_version", "assessment",
                                        "status", "message", "metrics", "findings", "children"}));
    EXPECT_EQ(root["schema_version"], 1);
    EXPECT_EQ(root["oehealth_version"], std::to_string(OEHEALTH_VERSION_MAJOR) + "." +
                                            std::to_string(OEHEALTH_VERSION_MINOR) + "." +
                                            std::to_string(OEHEALTH_VERSION_PATCH));
    EXPECT_EQ(root["status"], "WARN");
    EXPECT_EQ(keys(root["metrics"]),
              (std::vector<std::string>{"n_assessed", "n_outliers", "rmsz"}));
    EXPECT_DOUBLE_EQ(root["metrics"]["rmsz"].get<double>(), 0.83);
}

TEST(ToJson, FindingLayout) {
    const Json finding = Json::parse(to_json(sample_result()))["findings"][0];
    EXPECT_EQ(keys(finding),
              (std::vector<std::string>{"code", "severity", "message", "atoms", "observed",
                                        "expected", "deviation", "details"}));
    EXPECT_EQ(finding["severity"], "WARNING");
    EXPECT_EQ(keys(finding["atoms"][0]),
              (std::vector<std::string>{"atom_idx", "atom_name", "residue_name", "residue_number",
                                        "insert_code", "chain_id"}));
    EXPECT_EQ(finding["atoms"][0]["atom_idx"], 10);
    EXPECT_EQ(finding["atoms"][0]["insert_code"], "");
    EXPECT_TRUE(finding["deviation"].is_null());
    EXPECT_EQ(keys(finding["details"]), (std::vector<std::string>{"phi", "psi"}));
}

TEST(ToJson, ChildrenOmitVersionKeys) {
    const Json child = Json::parse(to_json(sample_result()))["children"][0];
    EXPECT_FALSE(child.contains("schema_version"));
    EXPECT_FALSE(child.contains("oehealth_version"));
    EXPECT_EQ(child["assessment"], "child");
    EXPECT_EQ(child["status"], "PASS");
}

TEST(ToJson, DuplicateNamesThrow) {
    AssessmentResult duplicate_metric = sample_result();
    duplicate_metric.metrics.push_back({"rmsz", 1.0});
    EXPECT_THROW(to_json(duplicate_metric), std::invalid_argument);

    AssessmentResult duplicate_detail = sample_result();
    duplicate_detail.findings[0].details.push_back({"phi", 1.0});
    EXPECT_THROW(to_json(duplicate_detail), std::invalid_argument);
}

TEST(ToJson, NegativeIndentIsCompact) {
    const std::string compact = to_json(sample_result(), -1);
    EXPECT_EQ(compact.find('\n'), std::string::npos);
    EXPECT_NE(to_json(sample_result()).find('\n'), std::string::npos);
}

TEST(ToJson, RealResultRoundTrips) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const Json root = Json::parse(to_json(BondLengthAssessment().Run(mol)));
    EXPECT_EQ(root["assessment"], "bond_lengths");
    EXPECT_TRUE(root["metrics"].contains("n_assessed"));
}

} // namespace OEHEALTH::test
