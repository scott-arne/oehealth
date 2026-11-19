#include "oehealth/suite.h"

#include <exception>
#include <stdexcept>
#include <utility>

namespace OEHEALTH {

AssessmentSuite::AssessmentSuite(std::string name, std::string description,
                                 std::vector<std::shared_ptr<Assessment>> children)
    : name_(std::move(name)), description_(std::move(description)), children_(std::move(children)) {
    if (name_.empty()) {
        throw std::invalid_argument("an assessment suite needs a name");
    }
    for (const auto& child : children_) {
        if (!child) {
            throw std::invalid_argument("assessment suite '" + name_ + "' has a null child");
        }
    }
}

std::string AssessmentSuite::GetName() const { return name_; }

std::string AssessmentSuite::GetDescription() const { return description_; }

const std::vector<std::shared_ptr<Assessment>>& AssessmentSuite::GetChildren() const {
    return children_;
}

AssessmentResult AssessmentSuite::Assess(const StructureContext& ctx) const {
    AssessmentResult result;
    result.assessment = name_;
    result.status = Status::NOT_APPLICABLE;
    if (children_.empty()) {
        result.message = "the suite has no assessments";
        return result;
    }
    for (const auto& child : children_) {
        AssessmentResult child_result;
        try {
            child_result = child->Assess(ctx);
        } catch (const std::exception& error) {
            child_result = AssessmentResult();
            child_result.assessment = child->GetName();
            child_result.status = Status::ERROR;
            child_result.message = error.what();
        }
        result.status = worst_status(result.status, child_result.status);
        result.children.push_back(std::move(child_result));
    }
    if (result.status == Status::NOT_APPLICABLE) {
        result.message = "no assessment in the suite applied";
    }
    return result;
}

} // namespace OEHEALTH
