#ifndef OEHEALTH_SRC_GEOMETRY_MEASURE_H
#define OEHEALTH_SRC_GEOMETRY_MEASURE_H

#include "oehealth/context.h"

#include <oechem.h>

#include <array>
#include <cstddef>
#include <optional>

namespace OEHEALTH::geometry {

/// \brief C(i)-N(j) distances at or above this, in angstroms, are chain breaks.
constexpr double PEPTIDE_BOND_CUTOFF = 2.5;

/// \brief Coordinates of an atom in the active conformer, in double precision.
///
/// OEGetDistance and friends compute in single precision; measuring from
/// double coordinates keeps threshold boundaries exact for molecules that
/// store doubles (OEMCMolType::DoubleCartesian).
std::array<double, 3> coordinates(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* atom);

/// \brief Distance between two atoms, in angstroms.
double distance(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* a,
                const OEChem::OEAtomBase* b);

/// \brief Angle a-b-c with b as the vertex, in degrees.
double angle_degrees(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* a,
                     const OEChem::OEAtomBase* b, const OEChem::OEAtomBase* c);

/// \brief Torsion a-b-c-d in degrees, in (-180, 180], IUPAC sign convention.
double torsion_degrees(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* a,
                       const OEChem::OEAtomBase* b, const OEChem::OEAtomBase* c,
                       const OEChem::OEAtomBase* d);

/// \brief The peptide successor of residue i, if it has one.
///
/// Residue j is the successor when it is the next entry in the residue
/// table, both are standard amino acids with the same chain id, C(i) and
/// N(j) are present, and their distance is under PEPTIDE_BOND_CUTOFF.
std::optional<std::size_t> peptide_successor(const StructureContext& ctx, std::size_t i);

/// \brief The peptide predecessor of residue j, if it has one.
std::optional<std::size_t> peptide_predecessor(const StructureContext& ctx, std::size_t j);

} // namespace OEHEALTH::geometry

#endif // OEHEALTH_SRC_GEOMETRY_MEASURE_H
