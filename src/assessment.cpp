#include "oehealth/assessment.h"

#include <stdexcept>

namespace OEHEALTH {

AssessmentResult Assessment::Run(const OEChem::OEMolBase& mol) const {
    if (mol.NumAtoms() > 0 && mol.GetDimension() != 3) {
        throw std::invalid_argument(GetName() + ": the molecule has no 3D coordinates");
    }
    const StructureContext ctx(mol);
    return Assess(ctx);
}

} // namespace OEHEALTH
