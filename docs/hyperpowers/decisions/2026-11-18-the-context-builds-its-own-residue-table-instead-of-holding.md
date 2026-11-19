# The context builds its own residue table instead of holding an `OEHierView`

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 3)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**The context builds its own residue table instead of holding an
`OEHierView`.** `OEBio::OEHierView` can only be constructed from a
non-const `OEMolBase&`, and its iteration methods are non-const, so it
cannot be built from the borrowed `const OEMolBase&`. Rejected: copying
the molecule to build a view (doubles memory and complicates atom index
mapping); `const_cast` (breaks the promise not to modify the input).

## What implementation changed

No change. `StructureContext::GetResidues` builds the residue table lazily from the borrowed `const OEMolBase&`, grouping atoms by `OEAtomGetResidue` and `OESameResidue` and ordering residues by first atom index.
