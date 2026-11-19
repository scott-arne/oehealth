#include "oehealth/json.h"

#include "oehealth/oehealth.h"

#include <nlohmann/json.hpp>

#include <cmath>
#include <set>
#include <stdexcept>
#include <string>

namespace OEHEALTH {

namespace {

using Json = nlohmann::ordered_json;

const char* severity_name(Severity severity) {
    switch (severity) {
        case Severity::INFO:
            return "INFO";
        case Severity::WARNING:
            return "WARNING";
        case Severity::ERROR:
            return "ERROR";
    }
    return "UNKNOWN";
}

const char* status_name(Status status) {
    switch (status) {
        case Status::NOT_APPLICABLE:
            return "NOT_APPLICABLE";
        case Status::PASS:
            return "PASS";
        case Status::WARN:
            return "WARN";
        case Status::FAIL:
            return "FAIL";
        case Status::ERROR:
            return "ERROR";
    }
    return "UNKNOWN";
}

Json number(double value) {
    if (std::isnan(value)) {
        return nullptr;
    }
    return value;
}

Json metrics_object(const std::vector<Metric>& metrics, const std::string& where) {
    Json object = Json::object();
    std::set<std::string> seen;
    for (const Metric& metric : metrics) {
        if (!seen.insert(metric.name).second) {
            throw std::invalid_argument("duplicate name '" + metric.name + "' in " + where);
        }
        object[metric.name] = number(metric.value);
    }
    return object;
}

Json atom_object(const AtomRef& atom) {
    Json object = Json::object();
    object["atom_idx"] = atom.atom_idx;
    object["atom_name"] = atom.atom_name;
    object["residue_name"] = atom.residue_name;
    object["residue_number"] = atom.residue_number;
    object["insert_code"] = atom.insert_code;
    object["chain_id"] = atom.chain_id;
    return object;
}

Json finding_object(const Finding& finding) {
    Json object = Json::object();
    object["code"] = finding.code;
    object["severity"] = severity_name(finding.severity);
    object["message"] = finding.message;
    Json atoms = Json::array();
    for (const AtomRef& atom : finding.atoms) {
        atoms.push_back(atom_object(atom));
    }
    object["atoms"] = atoms;
    object["observed"] = number(finding.observed);
    object["expected"] = number(finding.expected);
    object["deviation"] = number(finding.deviation);
    object["details"] =
        metrics_object(finding.details, "details of finding '" + finding.code + "'");
    return object;
}

Json result_object(const AssessmentResult& result, bool is_root) {
    Json object = Json::object();
    if (is_root) {
        object["schema_version"] = JSON_SCHEMA_VERSION;
        object["oehealth_version"] = std::to_string(OEHEALTH_VERSION_MAJOR) + "." +
                                     std::to_string(OEHEALTH_VERSION_MINOR) + "." +
                                     std::to_string(OEHEALTH_VERSION_PATCH);
    }
    object["assessment"] = result.assessment;
    object["status"] = status_name(result.status);
    object["message"] = result.message;
    object["metrics"] = metrics_object(result.metrics, "metrics of '" + result.assessment + "'");
    Json findings = Json::array();
    for (const Finding& finding : result.findings) {
        findings.push_back(finding_object(finding));
    }
    object["findings"] = findings;
    Json children = Json::array();
    for (const AssessmentResult& child : result.children) {
        children.push_back(result_object(child, false));
    }
    object["children"] = children;
    return object;
}

} // namespace

std::string to_json(const AssessmentResult& result, int indent) {
    return result_object(result, true).dump(indent);
}

} // namespace OEHEALTH
