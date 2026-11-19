# Chain identifiers use OEChem's extended chain id

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 19)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**Chain identifiers use OEChem's extended chain id.** `chain_id_of`
reads `OEResidue::GetExtChainID()`, not the one-character
`GetChainID()`, so multi-character mmCIF chains stay distinct in residue
grouping, peptide-neighbor checks, and `AtomRef`. Rejected: the
one-character chain id (merges chains such as `AA` and `AB` in large
mmCIF assemblies).

## What implementation changed

Implemented with the residue table in commit 6b0b8d5. `chain_id_of` in `src/common.cpp` returns the trimmed `GetExtChainID()`. `StructureContext` keys residues on it, `peptide_successor` in `src/geometry/measure.cpp` compares it to reject cross-chain neighbors, and `make_atom_ref` copies it into `AtomRef::chain_id`, which the JSON output writes as `chain_id`.
