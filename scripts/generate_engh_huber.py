"""Generate the Engh & Huber (2001) geometry tables from ``standard_geometry.cif``.

The input is the NDB/wwPDB file distributed with CCP4 as
``lib/data/standard_geometry.cif``. The output is a C++ source file holding
``constexpr`` tables of the heavy-atom bond and angle targets for the 20
standard amino acids, plus the peptide and disulfide link targets.

Usage::

    python3 scripts/generate_engh_huber.py standard_geometry.cif src/data/engh_huber_2001.cpp

The output is already clang-format clean; the tables are kept one row per
line inside ``clang-format off`` regions.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

STANDARD_RESIDUES = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    }
)  # fmt: skip

# The source spells the default peptide link both ways, so both map to DEFAULT.
LINK_VARIANTS = {
    "NOT GLY PRO": "DEFAULT",
    "NOT PRO GLY": "DEFAULT",
    "GLY": "GLY",
    "PRO": "PRO",
    "PRO-CIS": "PRO_CIS",
    "PRO-TRANS": "PRO_TRANS",
    "PRO-ALL": "PRO_ALL",
    "CYS": "DISULFIDE",
}

HYDROGEN_SYMBOLS = frozenset({"H", "D"})

# A CIF token is a quoted string or a run of non-blank characters. Atom names
# such as O3' contain a quote character but are not quoted strings.
TOKEN_PATTERN = re.compile(r"\"[^\"]*\"|'[^']*'|\S+")


class GeneratorError(Exception):
    """Raised when the input cannot be converted faithfully."""


@dataclass(frozen=True)
class Loop:
    """One CIF loop: the field names of its category and its rows."""

    fields: list[str]
    rows: list[dict[str, str]]


@dataclass(frozen=True)
class CifData:
    """The parts of the CIF file the generator uses."""

    items: dict[str, str]
    loops: dict[str, Loop]
    malformed: set[str]


@dataclass(frozen=True)
class ResidueRow:
    """An intra-residue bond (two atoms) or angle (three atoms) target."""

    residue: str
    atoms: tuple[str, ...]
    variant: str
    ideal: str
    sigma: str


@dataclass(frozen=True)
class LinkRow:
    """An inter-residue bond or angle target with its normalized link variant."""

    atoms: tuple[str, ...]
    link_variant: str
    ideal: str
    sigma: str


def tokenize(line: str) -> list[str]:
    """Split one CIF line into tokens, dropping quotes and trailing comments."""
    tokens: list[str] = []
    for match in TOKEN_PATTERN.finditer(line):
        token = match.group(0)
        if token.startswith("#"):
            break
        if len(token) >= 2 and token[0] == token[-1] and token[0] in "\"'":
            token = token[1:-1]
        tokens.append(token)
    return tokens


def split_tag(tag: str) -> tuple[str, str]:
    """Split ``_category.field`` into ``("_category", "field")``."""
    category, _, field = tag.partition(".")
    if not field:
        raise GeneratorError(f"malformed tag {tag!r}")
    return category, field


def read_cif(text: str) -> CifData:
    """Parse the items and loops of a single-block CIF file."""
    items: dict[str, str] = {}
    loops: dict[str, Loop] = {}
    malformed: set[str] = set()
    category = ""
    fields: list[str] = []
    values: list[str] = []
    in_loop = False

    def close_loop() -> None:
        if not in_loop:
            return
        width = len(fields)
        # Some loops the generator never reads have more values per row than
        # fields; record them and fail only if one of them is needed.
        if len(values) % width != 0:
            malformed.add(category)
            return
        rows = [
            dict(zip(fields, values[start : start + width], strict=True))
            for start in range(0, len(values), width)
        ]
        loops[category] = Loop(fields=list(fields), rows=rows)

    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("data_"):
            continue
        if stripped == "loop_":
            close_loop()
            in_loop = True
            category = ""
            fields = []
            values = []
            continue
        if stripped.startswith("_"):
            tokens = tokenize(stripped)
            tag_category, field = split_tag(tokens[0])
            if in_loop and not values:
                if category and tag_category != category:
                    raise GeneratorError(f"loop mixes {category} and {tag_category}")
                category = tag_category
                fields.append(field)
                continue
            close_loop()
            in_loop = False
            if len(tokens) != 2:
                raise GeneratorError(f"item {tokens[0]} has no single value")
            items[tokens[0]] = tokens[1]
            continue
        if not in_loop:
            raise GeneratorError(f"unexpected data line: {stripped!r}")
        values.extend(tokenize(stripped))
    close_loop()
    return CifData(items=items, loops=loops, malformed=malformed)


def require_loop(cif: CifData, category: str) -> Loop:
    """Return a loop or fail with the missing category's name."""
    if category in cif.malformed:
        raise GeneratorError(f"loop {category} has a partial row")
    loop = cif.loops.get(category)
    if loop is None:
        raise GeneratorError(f"missing loop {category}")
    return loop


def has_peptide_pair(atoms: tuple[str, ...]) -> bool:
    """Return whether C and N are consecutive anywhere in the atom path."""
    return any({first, second} == {"C", "N"} for first, second in pairwise(atoms))


def hydrogen_atoms(cif: CifData) -> set[tuple[str, str]]:
    """Return ``(residue, atom)`` pairs whose element is hydrogen."""
    loop = require_loop(cif, "_chem_comp_atom")
    return {
        (row["comp_id"], row["atom_id"])
        for row in loop.rows
        if row["type_symbol"] in HYDROGEN_SYMBOLS
    }


def residue_rows(
    cif: CifData, category: str, atom_fields: tuple[str, ...], value_field: str
) -> list[ResidueRow]:
    """Extract heavy-atom rows for the standard residues from one loop."""
    hydrogens = hydrogen_atoms(cif)
    loop = require_loop(cif, category)
    rows: list[ResidueRow] = []
    for row in loop.rows:
        residue = row["comp_id"]
        if residue not in STANDARD_RESIDUES:
            continue
        atoms = tuple(row[field] for field in atom_fields)
        if any((residue, atom) in hydrogens for atom in atoms):
            continue
        # Rows with consecutive C and N describe the bond to the neighboring
        # residue (PRO C-N-CD); that term comes from the link table instead.
        if has_peptide_pair(atoms):
            continue
        rows.append(
            ResidueRow(
                residue=residue,
                atoms=atoms,
                variant=row["description"],
                ideal=row[value_field],
                sigma=row[value_field + "_esd"],
            )
        )
    return rows


def link_rows(
    cif: CifData, category: str, atom_fields: tuple[str, ...], value_field: str
) -> list[LinkRow]:
    """Extract peptide and disulfide link rows, skipping nucleic-acid links."""
    loop = require_loop(cif, category)
    rows: list[LinkRow] = []
    for row in loop.rows:
        atoms = tuple(row[field] for field in atom_fields)
        is_peptide = has_peptide_pair(atoms)
        is_disulfide = all(atom == "SG" for atom in atoms)
        if not is_peptide and not is_disulfide:
            continue
        link_id = row["link_id"]
        variant = LINK_VARIANTS.get(link_id)
        if variant is None or (variant == "DISULFIDE") != is_disulfide:
            raise GeneratorError(f"unknown link id {link_id!r} on {category} row {'-'.join(atoms)}")
        rows.append(
            LinkRow(
                atoms=atoms,
                link_variant=variant,
                ideal=row[value_field],
                sigma=row[value_field + "_esd"],
            )
        )
    return rows


def quote(text: str) -> str:
    """Return a C++ string literal for an atom, residue, or variant name."""
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def emit_table(type_name: str, array_name: str, entries: list[str]) -> list[str]:
    """Emit one ``constexpr std::array`` definition, one row per line."""
    # clang-format would pack several rows onto each line, which makes the
    # tables hard to compare with the source file.
    lines = [
        "// clang-format off",
        f"constexpr std::array<{type_name}, {len(entries)}> {array_name}{{{{",
    ]
    lines.extend(f"    {entry}," for entry in entries)
    lines += ["}};", "// clang-format on", ""]
    return lines


def accessor(row_type: str, function: str, array: str) -> str:
    """Emit the accessor that returns a view of one table."""
    return f"RowSpan<{row_type}> {function}() {{ return {{{array}.data(), {array}.size()}}; }}"


def render(
    version: str,
    source_name: str,
    bonds: list[ResidueRow],
    angles: list[ResidueRow],
    link_bonds: list[LinkRow],
    link_angles: list[LinkRow],
) -> str:
    """Render the generated C++ source file."""

    def residue_entry(row: ResidueRow) -> str:
        names = ", ".join(quote(name) for name in (row.residue, *row.atoms, row.variant))
        return f"{{{names}, {row.ideal}, {row.sigma}}}"

    def link_entry(row: LinkRow) -> str:
        names = ", ".join(quote(name) for name in row.atoms)
        return f"{{{names}, LinkVariant::{row.link_variant}, {row.ideal}, {row.sigma}}}"

    lines = [
        "// Engh & Huber (2001) ideal bond lengths and angles for the standard amino acids.",
        "//",
        "// Generated by scripts/generate_engh_huber.py. Do not edit by hand.",
        "//",
        f'// Source: {source_name}, _standard_geometry.version "{version}"',
        "// Reference: Engh, R.A. & Huber, R. (2001). International Tables for",
        "// Crystallography, Vol. F, pp. 382-392. Dordrecht: Kluwer Academic Publishers.",
        "// Command:",
        f"//   python3 scripts/generate_engh_huber.py {source_name} src/data/engh_huber_2001.cpp",
        "",
        '#include "data/engh_huber_2001.h"',
        "",
        "#include <array>",
        "",
        "namespace OEHEALTH::data {",
        "namespace {",
        "",
    ]
    lines += emit_table("ResidueBond", "RESIDUE_BONDS", [residue_entry(r) for r in bonds])
    lines += emit_table("ResidueAngle", "RESIDUE_ANGLES", [residue_entry(r) for r in angles])
    lines += emit_table("LinkBond", "LINK_BONDS", [link_entry(r) for r in link_bonds])
    lines += emit_table("LinkAngle", "LINK_ANGLES", [link_entry(r) for r in link_angles])
    lines += [
        "} // namespace",
        "",
        accessor("ResidueBond", "residue_bonds", "RESIDUE_BONDS"),
        "",
        accessor("ResidueAngle", "residue_angles", "RESIDUE_ANGLES"),
        "",
        accessor("LinkBond", "link_bonds", "LINK_BONDS"),
        "",
        accessor("LinkAngle", "link_angles", "LINK_ANGLES"),
        "",
        f"const char* source_version() {{ return {quote(version)}; }}",
        "",
        "} // namespace OEHEALTH::data",
        "",
    ]
    return "\n".join(lines)


def generate(source: Path) -> str:
    """Read ``standard_geometry.cif`` and return the generated C++ source."""
    cif = read_cif(source.read_text(encoding="utf-8"))
    version = cif.items.get("_standard_geometry.version")
    if version is None:
        raise GeneratorError("missing _standard_geometry.version")
    return render(
        version=version,
        source_name=source.name,
        bonds=residue_rows(cif, "_chem_comp_bond", ("atom_id_1", "atom_id_2"), "value_dist"),
        angles=residue_rows(
            cif, "_chem_comp_angle", ("atom_id_1", "atom_id_2", "atom_id_3"), "value_angle"
        ),
        link_bonds=link_rows(cif, "_chem_link_bond", ("atom_id_1", "atom_id_2"), "value_dist"),
        link_angles=link_rows(
            cif, "_chem_link_angle", ("atom_id_1", "atom_id_2", "atom_id_3"), "value_angle"
        ),
    )


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description="Generate the Engh & Huber (2001) tables from standard_geometry.cif."
    )
    parser.add_argument("source", type=Path, help="path to standard_geometry.cif")
    parser.add_argument("output", type=Path, help="path of the C++ file to write")
    args = parser.parse_args(argv)
    try:
        text = generate(args.source)
    except GeneratorError as error:
        print(f"generate_engh_huber: {error}", file=sys.stderr)
        return 1
    args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
