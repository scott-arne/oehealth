# Findings refer to atoms by copied identifiers, not OpenEye pointers

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 10)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**Findings refer to atoms by copied identifiers, not OpenEye pointers.**
Results then outlive the molecule and serialize without lookups.

## What implementation changed

No change. `make_atom_ref` in `src/common.cpp` copies the identifiers into an `AtomRef`.
