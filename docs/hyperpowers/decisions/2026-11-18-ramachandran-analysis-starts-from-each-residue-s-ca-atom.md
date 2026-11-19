# Ramachandran analysis starts from each residue's CA atom

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 20)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**Ramachandran analysis starts from each residue's CA atom.**
`RamachandranAssessment` builds `OEBio::OERamachandranAnalysis` from the
residue's CA atom. Both constructors agree for structures without
alternate locations; a residue with alternate CA atoms is classified
from its first CA, as the bond and angle checks use the first atom of
each name. Rejected: the `(mol, residue)` constructor (it scans the
molecule for every residue, so time grows quadratically: 25 s at 20,000
residues against 0.5 s; it also leaves residues with alternate CA atoms
unassessed when alternate locations are kept).

## What implementation changed

Changed in commit b8abcc4. `RamachandranAssessment::Assess` constructs `OERamachandranAnalysis(alpha)` from the residue's CA instead of `OERamachandranAnalysis(mol, residue)`; on crambin copies, ramachandran time went from 24.96 s to 0.52 s at 20,010 residues (debug build). `test_ramachandran.cpp` counts Favored, Allowed and Outlier directly with the `(mol, residue)` constructor and compares each metric, and an alternate-location test checks that a residue with a second CA is still assessed.
