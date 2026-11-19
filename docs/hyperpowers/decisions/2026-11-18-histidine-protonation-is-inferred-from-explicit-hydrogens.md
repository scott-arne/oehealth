# Histidine protonation is inferred from explicit hydrogens, and otherwise chosen by best fit

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 13)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**Histidine protonation is inferred from explicit hydrogens, and
otherwise chosen by best fit.** Without hydrogens the protonation state is
unknown, and the variants differ by more than one sigma for some bonds.
Rejected: always using one variant (produces false outliers for the other
states).

## What implementation changed

No change. `src/geometry/targets.cpp` infers the HIS variant from HD1 and HE2 when present and otherwise picks the lowest-RMSZ variant, with ties going to HISE, then HISD, then HISH.
