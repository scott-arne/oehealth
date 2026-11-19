# Public headers guard against the Windows `ERROR` macro themselves

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 16)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**Public headers guard against the Windows `ERROR` macro themselves.**
On MSVC, OpenEye's `oetime.h` includes `winsock.h`, which brings in
`wingdi.h` and its `#define ERROR 0`, breaking `Severity::ERROR` and
`Status::ERROR`. `result.h` saves the macro with `#pragma push_macro`,
removes it around the enums, and restores it, and `NOGDI` is defined
privately on the library, the SWIG module, and the tests. Rejected:
renaming the enumerators (breaks the documented API and the JSON
strings); a public `NOGDI` (removes GDI from every downstream target that
links oehealth).

## What implementation changed

Added in commit e43c17a. `include/oehealth/result.h` saves `ERROR` with `#pragma push_macro`, undefines it around `Severity` and `Status`, and restores it. `NOGDI` is defined privately on the `oehealth` target, the SWIG module and the test target (`$<$<PLATFORM_ID:Windows>:NOGDI>`); commit 8922c5a had first defined it PUBLIC on `oehealth`. `tests/cpp/test_error_macro.cpp` checks that the macro survives the include and that the enumerators remain usable.
