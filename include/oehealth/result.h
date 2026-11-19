#ifndef OEHEALTH_RESULT_H
#define OEHEALTH_RESULT_H

#include <limits>
#include <string>
#include <vector>

namespace OEHEALTH {

// windows.h, reached through OpenEye's oetime.h, defines ERROR. Save and remove
// it around the enums so this header compiles regardless of include order, then
// restore it so including the header never changes the caller's macros.
#pragma push_macro("ERROR")
#undef ERROR

/// \brief How serious a single finding is.
enum class Severity { INFO, WARNING, ERROR };

/// \brief Overall outcome of an assessment, from best to worst.
///
/// ERROR is reserved for a suite child that threw; it is never derived from
/// findings.
enum class Status { NOT_APPLICABLE, PASS, WARN, FAIL, ERROR };

#pragma pop_macro("ERROR")

/// \brief Identifies an atom by copied values, so results outlive the molecule.
struct AtomRef {
    unsigned int atom_idx = 0;
    std::string atom_name;
    std::string residue_name;
    std::string chain_id;
    int residue_number = 0;
    std::string insert_code; ///< Empty when the residue has none.
};

/// \brief A named number. Names are unique within the vector that holds them.
struct Metric {
    std::string name;
    double value = 0.0;
};

/// \brief One problem, or one noteworthy observation, at specific atoms.
///
/// ``observed``, ``expected``, and ``deviation`` are NaN when not meaningful
/// for the finding's code.
struct Finding {
    std::string code; ///< Stable identifier, e.g. "bond_length_outlier".
    Severity severity = Severity::INFO;
    std::string message; ///< Human-readable; names the residue and atoms.
    std::vector<AtomRef> atoms;
    double observed = std::numeric_limits<double>::quiet_NaN();
    double expected = std::numeric_limits<double>::quiet_NaN();
    double deviation = std::numeric_limits<double>::quiet_NaN();
    std::vector<Metric> details; ///< Extra numbers, e.g. phi and psi.
};

/// \brief The outcome of running one assessment or suite.
struct AssessmentResult {
    std::string assessment; ///< GetName() of the producer.
    Status status = Status::NOT_APPLICABLE;
    std::string message; ///< Set for ERROR and NOT_APPLICABLE.
    std::vector<Metric> metrics;
    std::vector<Finding> findings;
    std::vector<AssessmentResult> children; ///< Filled only by suites.
};

/// \brief Return the worse of two statuses, by declaration order.
Status worst_status(Status a, Status b);

} // namespace OEHEALTH

#endif // OEHEALTH_RESULT_H
