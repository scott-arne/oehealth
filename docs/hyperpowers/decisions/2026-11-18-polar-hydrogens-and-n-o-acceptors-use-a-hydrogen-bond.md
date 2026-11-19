# Polar hydrogens and N/O acceptors use a hydrogen-bond contact distance

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 15)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**Polar hydrogens and N/O acceptors use a hydrogen-bond contact
distance.** A pair of a hydrogen bonded to N or O and an N or O atom
uses `hbond_contact_distance` (1.9 Å) instead of Bondi radii, so
hydrogen bonds in protonated structures are not reported as clashes. 1.9
Å is `polar_contact_distance` less a 1.0 Å donor–H bond, and lies within
the limit Probe applies (it allows donor-H/acceptor pairs 0.6 Å of
overlap on its own radii: 1.85 Å for H···O, 2.00 Å for H···N). Rejected:
Bondi radii for these pairs (1.20 + 1.52 Å flags every normal hydrogen
bond; crambin after `OEPlaceHydrogens` gave 60 warnings, 17 with this
rule); Probe's per-element radii and overlap allowance (adds a second
radius set to a check built on Bondi radii).

## What implementation changed

Added after the plan finished, in commit 7421251. `ClashOptions` gained `hbond_contact_distance` (default 1.9 Å) and its validation. `ClashAssessment` uses it for a polar hydrogen paired with an N or O atom, and the N/O rule takes precedence. The `OENearestNbrs` cutoff now accounts for both polar distances. `test_clashes.cpp` covers crambin after `OEPlaceHydrogens` and a hydrogen-bond-distance pair.
