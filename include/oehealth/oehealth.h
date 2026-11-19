#ifndef OEHEALTH_H
#define OEHEALTH_H

// Version information
#define OEHEALTH_VERSION_MAJOR 0
#define OEHEALTH_VERSION_MINOR 1
#define OEHEALTH_VERSION_PATCH 0

#include "oehealth/assessment.h"
#include "oehealth/assessments/bond_angles.h"
#include "oehealth/assessments/bond_lengths.h"
#include "oehealth/assessments/clashes.h"
#include "oehealth/assessments/ramachandran.h"
#include "oehealth/assessments/valence.h"
#include "oehealth/context.h"
#include "oehealth/json.h"
#include "oehealth/result.h"
#include "oehealth/suite.h"

#include <oechem.h>

namespace OEHEALTH {

/// \brief Calculate the molecular weight of a molecule.
///
/// Passes the molecule natively from Python to C++ via SWIG typemaps,
/// then delegates to OEChem's OECalculateMolecularWeight.
///
/// \param mol Reference to an OEMolBase object.
/// \returns Molecular weight in Daltons.
double calculate_molecular_weight(const OEChem::OEMolBase& mol);

} // namespace OEHEALTH

#endif // OEHEALTH_H
