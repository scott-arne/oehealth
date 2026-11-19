"""Tests for scripts/generate_engh_huber.py, run on a committed excerpt."""

from __future__ import annotations

import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "generate_engh_huber.py"
EXCERPT = ROOT / "tests" / "data" / "standard_geometry_excerpt.cif"
COMMITTED_TABLE = ROOT / "src" / "data" / "engh_huber_2001.cpp"
SHA256_LINE = re.compile(r"^// Source SHA-256: ([0-9a-f]{64})$", re.MULTILINE)

# Excerpt rows that the error-path tests modify.
PEPTIDE_BOND_ROW = ' C    N        1.336    0.023    "NOT GLY PRO"'
DISULFIDE_BOND_ROW = 'SG   SG        2.033    0.016    "CYS"'


def run_generator(source: Path, output: Path) -> subprocess.CompletedProcess[str]:
    """Run the generator as a separate process, as a developer would."""
    return subprocess.run(
        [sys.executable, str(SCRIPT), str(source), str(output)],
        capture_output=True,
        text=True,
        check=False,
    )


def generate_excerpt(tmp_path: Path) -> str:
    """Generate C++ from the excerpt and return its text."""
    output = tmp_path / "engh_huber_2001.cpp"
    result = run_generator(EXCERPT, output)
    assert result.returncode == 0, result.stderr
    return output.read_text(encoding="utf-8")


def test_emits_heavy_atom_residue_rows(tmp_path: Path) -> None:
    text = generate_excerpt(tmp_path)
    assert '{"ALA", "N", "CA", ".", 1.459, 0.020},' in text
    assert '{"ALA", "CA", "CB", ".", 1.520, 0.021},' in text
    assert '{"PRO", "N", "CD", "trans", 1.474, 0.014},' in text
    assert '{"ALA", "N", "CA", "C", ".", 111.0, 2.7},' in text
    assert '{"PRO", "CA", "N", "CD", "cis", 111.5, 1.4},' in text
    assert "constexpr std::array<ResidueBond, 4> RESIDUE_BONDS{{" in text
    assert "constexpr std::array<ResidueAngle, 3> RESIDUE_ANGLES{{" in text


def test_skips_hydrogens_and_nonstandard_residues(tmp_path: Path) -> None:
    text = generate_excerpt(tmp_path)
    assert '"H"' not in text
    assert '"DA"' not in text


def test_drops_intra_residue_rows_with_consecutive_c_and_n(tmp_path: Path) -> None:
    text = generate_excerpt(tmp_path)
    assert '{"PRO", "C", "N", "CD"' not in text


def test_normalizes_link_ids(tmp_path: Path) -> None:
    text = generate_excerpt(tmp_path)
    assert '{"C", "N", LinkVariant::DEFAULT, 1.336, 0.023},' in text
    # "NOT PRO GLY" is the source's other spelling of the default row.
    assert '{"O", "C", "N", LinkVariant::DEFAULT, 122.7, 1.6},' in text
    assert '{"C", "N", "CD", LinkVariant::PRO_TRANS, 128.4, 2.1},' in text
    assert '{"SG", "SG", LinkVariant::DISULFIDE, 2.033, 0.016},' in text


def test_skips_nucleic_acid_links(tmp_path: Path) -> None:
    text = generate_excerpt(tmp_path)
    assert "O3'" not in text
    assert "constexpr std::array<LinkBond, 2> LINK_BONDS{{" in text
    assert "constexpr std::array<LinkAngle, 2> LINK_ANGLES{{" in text


def test_records_source_version_and_command(tmp_path: Path) -> None:
    text = generate_excerpt(tmp_path)
    assert 'return "Fri Feb 22 17:25:15 GMT 2013";' in text
    assert "standard_geometry_excerpt.cif" in text
    assert "Engh, R.A. & Huber, R. (2001)" in text


def test_unknown_link_id_is_an_error(tmp_path: Path) -> None:
    source = tmp_path / "bad.cif"
    source.write_text(
        EXCERPT.read_text(encoding="utf-8").replace("PRO-TRANS", "PRO-OTHER"),
        encoding="utf-8",
    )
    output = tmp_path / "out.cpp"
    result = run_generator(source, output)
    assert result.returncode == 1
    assert "PRO-OTHER" in result.stderr
    assert "C-N-CD" in result.stderr
    assert not output.exists()


def generate_modified_excerpt(
    tmp_path: Path, old: str, new: str
) -> tuple[subprocess.CompletedProcess[str], Path]:
    """Run the generator on an excerpt copy with one text replacement."""
    text = EXCERPT.read_text(encoding="utf-8")
    assert old in text
    source = tmp_path / "modified.cif"
    source.write_text(text.replace(old, new), encoding="utf-8")
    output = tmp_path / "out.cpp"
    return run_generator(source, output), output


def test_partial_row_in_an_unread_loop_is_ignored(tmp_path: Path) -> None:
    result, output = generate_modified_excerpt(tmp_path, "ALA  A   1  16", "ALA  A   1")
    assert result.returncode == 0, result.stderr
    assert output.exists()


def test_partial_row_in_a_needed_loop_is_an_error(tmp_path: Path) -> None:
    result, output = generate_modified_excerpt(
        tmp_path, PEPTIDE_BOND_ROW, PEPTIDE_BOND_ROW + " extra"
    )
    assert result.returncode == 1
    assert "partial row" in result.stderr
    assert "_chem_link_bond" in result.stderr
    assert not output.exists()


def test_unknown_link_id_on_a_disulfide_row_is_an_error(tmp_path: Path) -> None:
    result, output = generate_modified_excerpt(
        tmp_path, DISULFIDE_BOND_ROW, DISULFIDE_BOND_ROW.replace('"CYS"', '"CYS-OTHER"')
    )
    assert result.returncode == 1
    assert "CYS-OTHER" in result.stderr
    assert "SG-SG" in result.stderr
    assert not output.exists()


def test_peptide_link_id_on_a_disulfide_row_is_an_error(tmp_path: Path) -> None:
    result, output = generate_modified_excerpt(
        tmp_path, DISULFIDE_BOND_ROW, DISULFIDE_BOND_ROW.replace('"CYS"', '"GLY"')
    )
    assert result.returncode == 1
    assert "'GLY'" in result.stderr
    assert "SG-SG" in result.stderr
    assert not output.exists()


def test_disulfide_link_id_on_a_peptide_row_is_an_error(tmp_path: Path) -> None:
    result, output = generate_modified_excerpt(
        tmp_path, PEPTIDE_BOND_ROW, PEPTIDE_BOND_ROW.replace('"NOT GLY PRO"', '"CYS"')
    )
    assert result.returncode == 1
    assert "'CYS'" in result.stderr
    assert "C-N" in result.stderr
    assert not output.exists()


def test_records_source_checksum(tmp_path: Path) -> None:
    text = generate_excerpt(tmp_path)
    digest = hashlib.sha256(EXCERPT.read_bytes()).hexdigest()
    assert f"// Source SHA-256: {digest}\n" in text


def test_committed_table_records_a_checksum() -> None:
    text = COMMITTED_TABLE.read_text(encoding="utf-8")
    assert SHA256_LINE.search(text) is not None


def find_full_source() -> Path | None:
    """Locate the full standard_geometry.cif, or return None if it is not installed."""
    override = os.environ.get("OEHEALTH_STANDARD_GEOMETRY_CIF")
    if override:
        candidate = Path(override)
        # An explicit setting that is wrong is a configuration error, not a reason to skip.
        if not candidate.is_file():
            pytest.fail(f"OEHEALTH_STANDARD_GEOMETRY_CIF names a missing file: {candidate}")
    elif os.environ.get("CCP4"):
        candidate = Path(os.environ["CCP4"]) / "lib" / "data" / "standard_geometry.cif"
    else:
        return None
    return candidate if candidate.is_file() else None


def test_regeneration_reproduces_the_committed_table(tmp_path: Path) -> None:
    source = find_full_source()
    if source is None:
        pytest.skip(
            "full standard_geometry.cif not found; set OEHEALTH_STANDARD_GEOMETRY_CIF "
            "(or CCP4) to check the committed table"
        )
    committed = COMMITTED_TABLE.read_bytes()
    match = SHA256_LINE.search(committed.decode("utf-8"))
    assert match is not None
    actual = hashlib.sha256(source.read_bytes()).hexdigest()
    if actual != match.group(1):
        pytest.skip(
            f"{source.name} has SHA-256 {actual}, but the committed table was generated "
            f"from {match.group(1)}; a different file proves nothing about this table"
        )
    # The generator writes the source file name into its header, so a correct file under another
    # name (for example a versioned copy) would otherwise differ from the committed table.
    named_source = tmp_path / "standard_geometry.cif"
    shutil.copyfile(source, named_source)
    output = tmp_path / "engh_huber_2001.cpp"
    result = run_generator(named_source, output)
    assert result.returncode == 0, result.stderr
    assert output.read_bytes() == committed
