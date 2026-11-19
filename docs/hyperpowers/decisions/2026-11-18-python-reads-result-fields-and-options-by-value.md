# Python reads result fields and options by value

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 17)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**Python reads result fields and options by value.** SWIG's `naturalvar`
feature makes member access return copies, and `OEHEALTH_RETURN_COPY`
makes `GetOptions()` return a copy, so a value taken from a temporary
result or assessment stays valid after Python collects the temporary.
Every new assessment with an options struct adds its own
`OEHEALTH_RETURN_COPY` line. Rejected: returning references (a reference
into a collected temporary dangles and crashes).

## What implementation changed

Implemented with the Python bindings in commit 41b3dca and recorded as a decision after the plan finished. `swig/oehealth.i` sets `%feature("naturalvar", "1")` and defines `OEHEALTH_RETURN_COPY`, applied to `BondLengthOptions`, `BondAngleOptions` and `ClashOptions`. The ownership tests in `tests/python/test_oehealth.py` check owning vector fields, a copied `GetOptions()` result, and `GetChildren()` returning owning proxies.
