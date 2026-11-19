#ifndef OEHEALTH_SUITE_H
#define OEHEALTH_SUITE_H

#include "oehealth/assessment.h"

#include <memory>
#include <string>
#include <vector>

namespace OEHEALTH {

/// \brief A named, ordered set of assessments that is itself an assessment.
///
/// Children are fixed at construction, so a suite cannot contain itself.
class AssessmentSuite : public Assessment {
public:
    /// \throws std::invalid_argument if the name is empty or any child is null.
    AssessmentSuite(std::string name, std::string description,
                    std::vector<std::shared_ptr<Assessment>> children);

    std::string GetName() const override;
    std::string GetDescription() const override;
    const std::vector<std::shared_ptr<Assessment>>& GetChildren() const;

    /// \brief Run every child in order against the same context.
    ///
    /// A child that throws a std::exception is recorded as a child result
    /// with status ERROR and the exception message; the remaining children
    /// still run. The suite's status is the worst child status, or
    /// NOT_APPLICABLE when it has no children.
    AssessmentResult Assess(const StructureContext& ctx) const override;

private:
    std::string name_;
    std::string description_;
    std::vector<std::shared_ptr<Assessment>> children_;
};

} // namespace OEHEALTH

#endif // OEHEALTH_SUITE_H
