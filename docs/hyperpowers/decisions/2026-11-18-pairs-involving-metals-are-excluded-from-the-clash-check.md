# Pairs involving metals are excluded from the clash check

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 7)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**Pairs involving metals are excluded from the clash check.** Bondi radii
make normal metal coordination (Zn–N near 2.1 Å) look like a 0.8 Å
overlap. Metal coordination belongs in its own future assessment.

## What implementation changed

No change. `ClashAssessment` skips atoms for which `OEIsMetal` is true, along with atoms of atomic number 0.
