#include "oehealth/oehealth.h"

#include <oechem.h>

namespace OEHEALTH {

double calculate_molecular_weight(const OEChem::OEMolBase& mol) {
    return OEChem::OECalculateMolecularWeight(mol);
}

} // namespace OEHEALTH
