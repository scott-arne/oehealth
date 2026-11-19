# The valence check assumes implicit hydrogen counts are assigned

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 12)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**The valence check assumes implicit hydrogen counts are assigned.**
`OECheckAtomValence` flags under-valent atoms as well as over-valent ones.
Probe on PDB entry 1CRN read with OpenEye's default PDB reader (no
explicit hydrogens, implicit counts assigned by the reader): 0 failures.
With implicit counts zeroed: 204 failures. The check reports every
failure, and the documentation states the requirement. Rejected:
reporting only over-valent atoms (would hide real under-valence in
prepared structures).

## What implementation changed

The valence check now runs `OECheckAtomValence` only on the elements OEChem's valence model covers (H, B, C, N, O, F, Si, P, S, Cl, Se, Br, I). Other elements, such as metals and noble gases, are skipped because OEChem accepts them unchecked and prints an "Unknown element in TypeCheck" warning through OEThrow. A table test in `test_valence.cpp` covers atomic numbers 1 to 118 and checks that nothing is written through OEThrow. The implicit-hydrogen requirement is unchanged.
