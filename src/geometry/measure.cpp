#include "geometry/measure.h"

#include <cmath>

namespace OEHEALTH::geometry {

namespace {

using Vec3 = std::array<double, 3>;

constexpr double PI = 3.14159265358979323846;
constexpr double RADIANS_TO_DEGREES = 180.0 / PI;

Vec3 subtract(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }

double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

} // namespace

std::array<double, 3> coordinates(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* atom) {
    Vec3 xyz{0.0, 0.0, 0.0};
    mol.GetCoords(atom, xyz.data());
    return xyz;
}

double distance(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* a,
                const OEChem::OEAtomBase* b) {
    return norm(subtract(coordinates(mol, a), coordinates(mol, b)));
}

double angle_degrees(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* a,
                     const OEChem::OEAtomBase* b, const OEChem::OEAtomBase* c) {
    const Vec3 vertex = coordinates(mol, b);
    const Vec3 u = subtract(coordinates(mol, a), vertex);
    const Vec3 v = subtract(coordinates(mol, c), vertex);
    // atan2 of |u x v| and u.v is accurate near 0 and 180 degrees, unlike acos.
    return std::atan2(norm(cross(u, v)), dot(u, v)) * RADIANS_TO_DEGREES;
}

double torsion_degrees(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* a,
                       const OEChem::OEAtomBase* b, const OEChem::OEAtomBase* c,
                       const OEChem::OEAtomBase* d) {
    const Vec3 pb = coordinates(mol, b);
    const Vec3 pc = coordinates(mol, c);
    const Vec3 b1 = subtract(pb, coordinates(mol, a));
    const Vec3 b2 = subtract(pc, pb);
    const Vec3 b3 = subtract(coordinates(mol, d), pc);
    const double y = norm(b2) * dot(b1, cross(b2, b3));
    const double x = dot(cross(b1, b2), cross(b2, b3));
    const double degrees = std::atan2(y, x) * RADIANS_TO_DEGREES;
    // atan2 returns -pi for y == -0.0 and x < 0 (planar trans); the range is (-180, 180].
    return degrees <= -180.0 ? 180.0 : degrees;
}

std::optional<std::size_t> peptide_successor(const StructureContext& ctx, std::size_t i) {
    const auto& residues = ctx.GetResidues();
    const std::size_t j = i + 1;
    if (j >= residues.size()) {
        return std::nullopt;
    }
    const ResidueInfo& first = residues[i];
    const ResidueInfo& second = residues[j];
    if (!first.is_standard_amino_acid || !second.is_standard_amino_acid ||
        first.chain_id != second.chain_id) {
        return std::nullopt;
    }
    const OEChem::OEAtomBase* carbon = first.FindAtom("C");
    const OEChem::OEAtomBase* nitrogen = second.FindAtom("N");
    if (carbon == nullptr || nitrogen == nullptr) {
        return std::nullopt;
    }
    if (distance(ctx.GetMolecule(), carbon, nitrogen) >= PEPTIDE_BOND_CUTOFF) {
        return std::nullopt;
    }
    return j;
}

std::optional<std::size_t> peptide_predecessor(const StructureContext& ctx, std::size_t j) {
    if (j == 0) {
        return std::nullopt;
    }
    const auto successor = peptide_successor(ctx, j - 1);
    if (successor && *successor == j) {
        return j - 1;
    }
    return std::nullopt;
}

} // namespace OEHEALTH::geometry
