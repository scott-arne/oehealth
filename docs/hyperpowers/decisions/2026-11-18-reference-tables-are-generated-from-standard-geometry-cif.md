# Reference tables are generated from `standard_geometry.cif`, not transcribed by hand

- Status: accepted
- Spec: docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md (decision 5)
- Recorded: 2026-11-19; notes from the branch diff at finish

## Decision

**Reference tables are generated from `standard_geometry.cif`, not
transcribed by hand.** The file (NDB/wwPDB, version "Fri Feb 22 17:25:15
GMT 2013", distributed with CCP4 9 as `lib/data/standard_geometry.cif`)
states that its amino-acid values come from Engh & Huber (2001), and spot
checks against the published values match. A committed Python script
converts it to C++; the generated file is committed so builds and tests
never need CCP4. Rejected: manual transcription (error-prone and
unverifiable); the current wwPDB validation targets (they no longer use
Engh & Huber; for example N–CA is 1.447 ± 0.012 Å there versus 1.459 ±
0.020 Å in Engh & Huber).

## What implementation changed

The generated `src/data/engh_huber_2001.cpp` header now records a `// Source SHA-256:` line, and `scripts/generate_engh_huber.py` writes it. `tests/python/test_generate_engh_huber.py` checks that the line is present and regenerates the table byte for byte from the full CIF when it is available. The regeneration test skips when the CIF is absent or its digest differs, and fails when `OEHEALTH_STANDARD_GEOMETRY_CIF` names a missing file.
