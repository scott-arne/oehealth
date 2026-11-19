# Clashes are assessed on the atoms as given

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 6)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**Clashes are assessed on the atoms as given.** The input is never
modified, and the result records whether hydrogens were present.
Rejected: adding hydrogens to a copy (slower, depends on protonation, and
adds a failure mode to a simple check); requiring hydrogens (most raw PDB
input would get no result).

## What implementation changed

No change.
