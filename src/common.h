#ifndef OEHEALTH_SRC_COMMON_H
#define OEHEALTH_SRC_COMMON_H

#include "oehealth/result.h"

#include <oechem.h>

#include <string>
#include <vector>

namespace OEHEALTH {

/// \brief Return text with leading and trailing whitespace removed.
std::string trim(const std::string& text);

/// \brief Return whether a residue name is one of the 20 standard amino acids.
bool is_standard_amino_acid(const std::string& residue_name);

/// \brief Return an atom's name with surrounding whitespace removed.
std::string trimmed_atom_name(const OEChem::OEAtomBase* atom);

/// \brief Convert an OEChem insertion code character to a string.
///
/// A blank (space or NUL) becomes the empty string.
std::string residue_char(char value);

/// \brief Return a residue's full chain id, trimmed; a blank id becomes "".
///
/// Uses GetExtChainID(): GetChainID() returns only the first character, so
/// mmCIF chains such as "AA" and "AB" would collide.
std::string chain_id_of(const OEChem::OEResidue& residue);

/// \brief Copy an atom's identifiers into an AtomRef.
AtomRef make_atom_ref(const OEChem::OEAtomBase* atom);

/// \brief Return "chain name number[insert]", e.g. "A ALA 5" or "B GLY 52A".
///
/// The chain is omitted when empty.
std::string residue_label(const AtomRef& atom);

/// \brief Describe atoms for a message, e.g. "A ALA 5 CA-CB" or "A ALA 5 C-A GLY 6 N".
///
/// The residue label is repeated whenever the residue changes along the path.
std::string atom_path(const std::vector<AtomRef>& atoms);

/// \brief Format a number with a fixed number of decimals.
std::string format_fixed(double value, int decimals);

/// \brief Derive an assessment status from its findings.
///
/// Any ERROR finding gives FAIL, otherwise any WARNING gives WARN, otherwise
/// PASS. INFO findings never affect status.
Status status_from_findings(const std::vector<Finding>& findings);

/// \brief Whether a value meets an inclusive threshold, within THRESHOLD_TOLERANCE.
bool meets_threshold(double value, double threshold);

/// \brief Throw std::invalid_argument unless a warning/error threshold pair is valid.
///
/// Both must be finite and non-negative, and error must not be less than
/// warning. ``what`` names the option pair in the exception message.
void validate_thresholds(const std::string& what, double warning, double error);

} // namespace OEHEALTH

#endif // OEHEALTH_SRC_COMMON_H
