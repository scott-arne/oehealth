"""Tests for scripts/generate_engh_huber.py, run on a committed excerpt."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "generate_engh_huber.py"
EXCERPT = ROOT / "tests" / "data" / "standard_geometry_excerpt.cif"


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
