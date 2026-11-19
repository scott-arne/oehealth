#ifndef OEHEALTH_JSON_H
#define OEHEALTH_JSON_H

#include "oehealth/result.h"

#include <string>

namespace OEHEALTH {

/// \brief Version of the JSON layout written by to_json.
constexpr int JSON_SCHEMA_VERSION = 1;

/// \brief Serialize a result, and its children, to JSON.
///
/// The root object carries ``schema_version`` and ``oehealth_version``;
/// children omit them. Metrics and details are objects that keep their
/// order, NaN is written as null, and enums are written by name.
///
/// \param indent Spaces per level; a negative value gives compact output.
/// \throws std::invalid_argument if a metric or detail name repeats within
///         one vector.
std::string to_json(const AssessmentResult& result, int indent = 2);

} // namespace OEHEALTH

#endif // OEHEALTH_JSON_H
