#ifndef OEHEALTH_DATA_ENGH_HUBER_2001_H
#define OEHEALTH_DATA_ENGH_HUBER_2001_H

#include <cstddef>

namespace OEHEALTH::data {

/// \brief Which peptide or disulfide link row applies to a pair of residues.
///
/// The source file names these by link id; the generator normalizes both
/// spellings of the default row ("NOT GLY PRO" and "NOT PRO GLY") to DEFAULT.
enum class LinkVariant { DEFAULT, GLY, PRO, PRO_CIS, PRO_TRANS, PRO_ALL, DISULFIDE };

/// \brief An intra-residue bond target. Lengths are in angstroms.
///
/// ``variant`` is "." when the row applies to every form of the residue, and
/// otherwise one of "SS", "SH" (CYS), "cis", "trans" (PRO), or "HISD",
/// "HISE", "HISH" (HIS).
struct ResidueBond {
    const char* residue;
    const char* atom1;
    const char* atom2;
    const char* variant;
    double ideal;
    double sigma;
};

/// \brief An intra-residue angle target. Angles are in degrees; atom2 is the vertex.
struct ResidueAngle {
    const char* residue;
    const char* atom1;
    const char* atom2;
    const char* atom3;
    const char* variant;
    double ideal;
    double sigma;
};

/// \brief A link bond target: atom1 is in residue i, atom2 in residue j.
struct LinkBond {
    const char* atom1;
    const char* atom2;
    LinkVariant link_variant;
    double ideal;
    double sigma;
};

/// \brief A link angle target, in degrees, with atom2 as the vertex.
///
/// For peptide links, atoms named before the C-N pair belong to residue i
/// and atoms after it to residue j.
struct LinkAngle {
    const char* atom1;
    const char* atom2;
    const char* atom3;
    LinkVariant link_variant;
    double ideal;
    double sigma;
};

/// \brief A read-only view of one generated table.
template <typename Row>
struct RowSpan {
    const Row* first;
    std::size_t count;

    const Row* begin() const { return first; }
    const Row* end() const { return first + count; }
    std::size_t size() const { return count; }
};

RowSpan<ResidueBond> residue_bonds();
RowSpan<ResidueAngle> residue_angles();
RowSpan<LinkBond> link_bonds();
RowSpan<LinkAngle> link_angles();

/// \brief The ``_standard_geometry.version`` string of the source file.
const char* source_version();

} // namespace OEHEALTH::data

#endif // OEHEALTH_DATA_ENGH_HUBER_2001_H
