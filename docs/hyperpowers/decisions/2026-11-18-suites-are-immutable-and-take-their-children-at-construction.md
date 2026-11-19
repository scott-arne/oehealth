# Suites are immutable and take their children at construction

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 9)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**Suites are immutable and take their children at construction.** Suites
hold `std::shared_ptr<Assessment>`, which SWIG handles cleanly from
Python. With no `Add` method, a suite cannot contain itself or form a
cycle, and every assessment, suites included, is immutable after
construction, so sharing children between suites and threads is safe.
Rejected: a mutable `Add` (needs cycle detection and a freeze step);
`std::unique_ptr` ownership (fragile to transfer from Python).

## What implementation changed

No change.
