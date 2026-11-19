# Peptide-link variants follow MolProbity's Dangle

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 18)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**Peptide-link variants follow MolProbity's Dangle.** Engh & Huber
(2001) does not say which residue owns the GLY and PRO peptide-link
rows. Dangle's parameter table (`EnghHuber_IntlTblsF_1999.txt` in
`rlabduke/javadev`) resolves this against the 1991 paper: the C–N bond
and the C–N–CA angle key on the residue that owns N; CA–C–N keys on GLY
in residue i first, then PRO in residue j; O–C–N keys on PRO in residue
j first, then GLY in residue i. oehealth uses the same rules. Rejected:
PRO precedence for both CA–C–N and O–C–N in a Gly–Pro pair (the first
design, which differs from the reference implementation for CA–C–N).

## What implementation changed

Added in commit fbc1cbf. `src/geometry/targets.cpp` selects link variants by the Dangle rules. C-N and C-N-CA key on the residue that owns N, CA-C-N prefers GLY in residue i and then PRO in j, and O-C-N prefers PRO in j and then GLY in i. `test_targets.cpp` and `test_reference_data.cpp` check that every peptide term resolves to exactly one row for each variant it can take.
