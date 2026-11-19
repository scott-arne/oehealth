"""Shared fixtures and configuration for oehealth Python tests."""

from pathlib import Path

import pytest

pytest.importorskip("openeye.oechem", reason="OpenEye Toolkits not installed")

DATA_DIR = Path(__file__).resolve().parents[1] / "data"


@pytest.fixture
def crambin():
    """Read PDB entry 1CRN (crambin, 46 residues, no hydrogens)."""
    from openeye import oechem

    mol = oechem.OEGraphMol()
    stream = oechem.oemolistream(str(DATA_DIR / "1crn.pdb"))
    assert oechem.OEReadMolecule(stream, mol)
    return mol


@pytest.fixture
def ethanol_2d():
    """Create ethanol from SMILES, which has atoms but no 3D coordinates."""
    from openeye import oechem

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "CCO")
    return mol
