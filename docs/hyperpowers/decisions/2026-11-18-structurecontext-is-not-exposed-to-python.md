# `StructureContext` is not exposed to Python

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 11)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**`StructureContext` is not exposed to Python.** It borrows the molecule,
and Python could collect the molecule while the context lives. Python
calls `Run(mol)` instead.

## What implementation changed

No change.
