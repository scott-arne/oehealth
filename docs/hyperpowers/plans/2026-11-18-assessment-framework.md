# Structure Health Assessment Framework Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use hyperpowers:subagent-driven-development (recommended) or hyperpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. A correction added during execution that alters one of the spec's decisions names it: `Correction (<date>), decision <N> (<title>):`.

**Spec:** `docs/hyperpowers/specs/2026-11-18-assessment-framework-design.md`

**Goal:** Add a C++ framework for assessing structure health (result model, structure context, assessment interface, suites, JSON) with five protein checks (bond lengths, bond angles, clashes, valence, Ramachandran), expose it to Python through SWIG, and remove the template's `calculate_molecular_weight`.

**Architecture:** An abstract `OEHEALTH::Assessment` assesses a `StructureContext`, which borrows the molecule and caches a residue table. Every assessment returns the same `AssessmentResult` type, and `AssessmentSuite` is itself an assessment that runs its children in order and records a child that throws as an ERROR result. Bond and angle targets come from Engh & Huber (2001) tables that a committed script generates from CCP4's `standard_geometry.cif`.

**Tech Stack:** C++17 with CMake presets, OpenEye Toolkits 2026.1 (OEChem, OEBio), GoogleTest 1.15.2, nlohmann/json 3.11.3, SWIG 4 through cmake-openeye 1.2.2, Python 3.12+ with pytest, ruff and mypy, clang-format and clang-tidy.

## Global Constraints

- C++17, CMake as already configured; new sources go in `OEHEALTH_SOURCES`.
- clang-format: a committed `.clang-format`; all C++ sources formatted with
  it. The generated `engh_huber_2001.cpp` is formatted too.
- clang-tidy: a committed `.clang-tidy` enabling bugprone, modernize, and
  performance checks; new code is clean against it.
- Python: `[tool.ruff]` and `[tool.mypy]` sections in `pyproject.toml`
  covering `python/`, `scripts/generate_engh_huber.py`, and `tests/python/`;
  that code passes both. The SWIG-generated module and the existing
  `scripts/build_python.py` are excluded.
- Every behavior change comes with GoogleTest and, where it reaches Python,
  pytest coverage.
- Version strings are managed by vrzn; none are hand-edited.
- **Merge-expected files:** `CMakeLists.txt`, `tests/cpp/CMakeLists.txt`

## Grounding

- `CMakeLists.txt:52-70`: the library target: the `OEHEALTH_SOURCES` list, `target_include_directories`, and `target_link_libraries` with `PUBLIC` sections. New sources and links follow this layout.
- `CMakeLists.txt:30-38`: the FetchContent pattern (`FetchContent_Declare` with `GIT_REPOSITORY`, `GIT_TAG`, `GIT_SHALLOW TRUE`, then `FetchContent_MakeAvailable`) that the nlohmann/json block copies.
- `tests/cpp/CMakeLists.txt:1-22`: the test executable, the GoogleTest fetch, and `gtest_discover_tests`.
- `tests/cpp/test_oehealth.cpp:1-40`: GoogleTest file layout: `<gtest/gtest.h>` first, then project headers, then `<oechem.h>`, with tests inside the `OEHEALTH::test` namespace. New tests keep the layout but write `namespace OEHEALTH::test {`, the C++17 form that clang-tidy's `modernize-concat-nested-namespaces` requires.
- `src/oehealth.cpp:1-11`: implementation file layout: own header first, then OpenEye headers, code inside `namespace OEHEALTH`, closed with `} // namespace OEHEALTH`.
- `include/oehealth/oehealth.h:1-24`: header guard naming (`OEHEALTH_H`), the vrzn-managed version macros (lines 5-7), and Doxygen `///` comments with `\brief`, `\param` and `\returns`.
- `swig/oehealth.i:355-370`: the STL typemap includes and the `%exception` block that turns C++ exceptions into Python exceptions.
- `swig/oehealth.i:379-386`: the wrapped API section, which Task 12 replaces.
- `swig/CMakeLists.txt:1-14`: `openeye_add_swig_module` with `SWIG_FLAGS` and `EXPECTED_LIB_VARS`.
- `pyproject.toml:63-68`: `expected-missing-libs` for wheel repair.
- `python/oehealth/__init__.py:559-567`: the extension import list and `__all__`.
- `tests/python/conftest.py:1-25`: fixtures skip when OpenEye is missing (`pytest.importorskip`) and import `openeye.oechem` inside the fixture.
- `tests/python/test_oehealth.py:1-45`: pytest classes with docstrings; lines 44-45 are the version assertions that vrzn manages, which Task 12 keeps character for character.
- Error handling in C++: none: the template throws nothing. The spec defines `std::invalid_argument` for invalid input, mapped to `ValueError` in the SWIG `%exception` block.
- Code generation: none: the repository has no code generator. `scripts/build_python.py` builds wheels and is unrelated.

## How to run things

Run every command from the repository root.

- Configure once per build directory: `cmake --preset debug`. A fresh build directory downloads cmake-openeye, GoogleTest and nlohmann/json from GitHub.
- Build: `cmake --build build-debug -j8`. `OEHEALTH_BUILD_PYTHON` is on by default, so this also regenerates and builds the SWIG module into `python/oehealth/` (ignored by git).
- One C++ test group: `./build-debug/tests/cpp/oehealth_tests --gtest_filter='Suite.*'`. All C++ tests: `ctest --test-dir build-debug --output-on-failure`.
- Python tests: `PYTHONPATH=python python -m pytest tests/python -q`, where `python` is an interpreter with the OpenEye toolkits and pytest installed.
- clang-format: `clang-format --dry-run --Werror <files>`.
- clang-tidy: `clang-tidy -p build-debug --quiet --extra-arg=-isysroot --extra-arg="$(xcrun --show-sdk-path)" <.cpp files>`. Pass the file names directly; zsh does not split an unquoted variable into words. A clean run prints only "N warnings generated" lines, or nothing.
- Python lint: `ruff check .`, `ruff format --check scripts tests/python`, and `mypy --python-executable "$(command -v python)"`.

## Design notes beyond the spec

These choices were made while writing the code and are part of this plan:

- Three private helper units not named in the spec's file layout: `src/common.{h,cpp}` (trimming, residue labels, `AtomRef` construction, status derivation, threshold checks), `src/geometry/measure.{h,cpp}` (double-precision distance, angle, torsion, and peptide neighbors), and `src/assessments/geometry_outliers.{h,cpp}` (the logic shared by bond lengths and bond angles).
- Measurements read coordinates with `OEMolBase::GetCoords(atom, double*)`, because `OEGetDistance` computes in single precision. Tests of exact thresholds build `OEMol(OEMCMolType::DoubleCartesian)` molecules.
- `json.h` declares `constexpr int JSON_SCHEMA_VERSION = 1`, the value written as `schema_version`.
- Chain ids come from `OEResidue::GetExtChainID()`, trimmed (`chain_id_of` in `src/common.h`). `GetChainID()` returns only the first character, so mmCIF chains such as "AA" and "AB" would get the same id and could be joined by a false peptide link; `StructureContext.KeepsExtendedChainIds` covers this.
- The numeric fields of `Finding` default to NaN, so a finding that never sets them serializes them as `null`.
- Metrics that stay defined on NOT_APPLICABLE: `n_assessed` = 0 and `n_outliers` = 0 (bond lengths, bond angles), `n_assessed` = 0 (Ramachandran), `hydrogens_present` (clashes), `n_invalid` = 0 (valence).
- The Ramachandran outlier test sets residue 44's phi and psi to 0°, which OEBio classifies as Outlier.
- `ruff format --check` covers `scripts` and `tests/python`. `python/oehealth/__init__.py` comes from the cmake-openeye template and is not ruff-formatted; reformatting it would be an unrelated change. `ruff check .` covers it.
- Tasks 11 and 12 split the template removal: Task 11 keeps `calculate_molecular_weight` so the SWIG module, which still wraps it, keeps building; Task 12 removes it together with the SWIG changes.
- The spec, this plan and the decision records are committed; `docs/hyperpowers/.gitignore` keeps everything else in that directory out of git.

Assumption: the GLY row for the peptide C–N bond is keyed on residue j (the residue that owns N), validate via comparison with the Engh & Huber (2001) table text by the project owner (Scott Arne Johnson), before the branch is merged. Task 5 proceeds on the assumption because the spec bounds its effect (the GLY and DEFAULT C–N targets differ by 0.010 Å, about half a sigma) and the choice sits in one condition: `is_residue(second, "GLY")` in `add_link_bond` in `src/geometry/targets.cpp`, checked by the C(19)–N(20) expectation in `Targets.CrambinLinkTerms`. A different answer changes that condition and that expectation.

Assumption: redistributing the Engh & Huber (2001) values in this MIT-licensed repository is acceptable because they are the published numbers, validate via review by the project owner before the first public release. No task depends on the answer.

---

### Task 1: Tooling, build scaffolding, and the crambin fixture

**Risk tier:** standard — multi-file build configuration (CMake, SWIG CMake, pyproject) that every later task builds on.

**Depends on:** none

**Files:**
- Create: `.clang-format`
- Create: `.clang-tidy`
- Create: `tests/data/1crn.pdb` (downloaded)
- Modify: `CMakeLists.txt:50-70`
- Modify: `tests/cpp/CMakeLists.txt:11-19`
- Modify: `swig/CMakeLists.txt:9-12`
- Modify: `pyproject.toml:65-68`, and append ruff and mypy sections at the end

**Interfaces:**
- Consumes: nothing.
- Produces: `src/` as a private include directory of the `oehealth` library and of `oehealth_tests` (so `#include "common.h"` and `#include "geometry/measure.h"` resolve); `tests/cpp/` as an include directory of `oehealth_tests` (for `test_helpers.h`); the compile definition `OEHEALTH_TEST_DATA_DIR`, the absolute path of `tests/data` as a string literal; `OpenEye::OEBio` linked `PUBLIC`; `nlohmann_json::nlohmann_json` linked privately by the library and directly by the tests; `tests/data/1crn.pdb`.

**Mirror:** `CMakeLists.txt:30-38`, the FetchContent block layout.

This task changes configuration only, so it has no failing test of its own. Its check is that the template's existing C++ and Python tests still build and pass with the new configuration.

- [ ] **Step 1: Create `.clang-format`**

`.clang-format`:

```yaml
---
Language: Cpp
BasedOnStyle: Google
IndentWidth: 4
ColumnLimit: 100
AccessModifierOffset: -4
SpacesBeforeTrailingComments: 1
PointerAlignment: Left
DerivePointerAlignment: false
IncludeBlocks: Preserve
...
```

- [ ] **Step 2: Create `.clang-tidy`**

`.clang-tidy`:

```yaml
---
# bugprone, modernize, and performance checks, minus the ones that conflict
# with project conventions:
#   modernize-use-trailing-return-type  project style uses leading return types
#   modernize-macro-to-enum             version macros are managed by vrzn
#   modernize-use-nodiscard             not used in this codebase
#   bugprone-easily-swappable-parameters  flags every (double, double) signature
#   performance-enum-size               SWIG and JSON code expect plain int enums
Checks: >
  -*,
  bugprone-*,
  modernize-*,
  performance-*,
  -modernize-use-trailing-return-type,
  -modernize-macro-to-enum,
  -modernize-use-nodiscard,
  -bugprone-easily-swappable-parameters,
  -performance-enum-size
WarningsAsErrors: ''
HeaderFilterRegex: '.*/(include/oehealth|src|tests/cpp)/.*'
FormatStyle: file
...
```

- [ ] **Step 3: Update the root `CMakeLists.txt`**

Add the nlohmann/json fetch after `find_package(OpenEye REQUIRED)`, the private `src` include directory, and the OEBio and json links. Leave `OEHEALTH_SOURCES` as it is. Apply this diff with `git apply`:

```diff
diff --git a/CMakeLists.txt b/CMakeLists.txt
--- a/CMakeLists.txt
+++ b/CMakeLists.txt
@@ -49,6 +49,15 @@ endif()
 
 find_package(OpenEye REQUIRED)
 
+# nlohmann/json is header-only and used only by src/json.cpp.
+FetchContent_Declare(
+    nlohmann_json
+    GIT_REPOSITORY https://github.com/nlohmann/json.git
+    GIT_TAG v3.11.3
+    GIT_SHALLOW TRUE
+)
+FetchContent_MakeAvailable(nlohmann_json)
+
 # Library
 set(OEHEALTH_SOURCES
     src/oehealth.cpp
@@ -62,11 +71,19 @@ target_include_directories(oehealth
     PUBLIC
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include>
+    PRIVATE
+        ${CMAKE_CURRENT_SOURCE_DIR}/src
 )
 
+# OERamachandranAnalysis lives in OEBio, which depends on OEGrid. The
+# BUILD_INTERFACE wrapper keeps the header-only json target out of the
+# installed export set.
 target_link_libraries(oehealth
     PUBLIC
+        OpenEye::OEBio
         OpenEye::OEChem
+    PRIVATE
+        $<BUILD_INTERFACE:nlohmann_json::nlohmann_json>
 )
 
 # Tests
```

- [ ] **Step 4: Update `tests/cpp/CMakeLists.txt`**

```diff
diff --git a/tests/cpp/CMakeLists.txt b/tests/cpp/CMakeLists.txt
--- a/tests/cpp/CMakeLists.txt
+++ b/tests/cpp/CMakeLists.txt
@@ -12,9 +12,23 @@ add_executable(oehealth_tests
     test_oehealth.cpp
 )
 
+# Tests reach private headers under src/ (reference tables, variant
+# selection, helpers) as well as the public API.
+target_include_directories(oehealth_tests
+    PRIVATE
+        ${PROJECT_SOURCE_DIR}/src
+        ${CMAKE_CURRENT_SOURCE_DIR}
+)
+
+target_compile_definitions(oehealth_tests
+    PRIVATE
+        OEHEALTH_TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/tests/data"
+)
+
 target_link_libraries(oehealth_tests
     PRIVATE
         oehealth
+        nlohmann_json::nlohmann_json
         GTest::gtest_main
 )
 
```

- [ ] **Step 5: Update `swig/CMakeLists.txt`**

SWIG needs `-I include` to find the headers that Task 12 `%include`s, and the wheel check needs the two new OpenEye libraries.

```diff
diff --git a/swig/CMakeLists.txt b/swig/CMakeLists.txt
--- a/swig/CMakeLists.txt
+++ b/swig/CMakeLists.txt
@@ -6,9 +6,9 @@ openeye_add_swig_module(
     LINK_LIBS oehealth
     PYTHON_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/python/oehealth
     STABLE_ABI ${OEHEALTH_USE_STABLE_ABI}
-    SWIG_FLAGS "-doxygen"
+    SWIG_FLAGS "-doxygen" "-I${CMAKE_SOURCE_DIR}/include"
     EXPECTED_LIB_VARS
         OECHEM_LIBRARY OESYSTEM_LIBRARY OEPLATFORM_LIBRARY
-        OEMATH_LIBRARY OEZSTD_LIBRARY
+        OEMATH_LIBRARY OEZSTD_LIBRARY OEBIO_LIBRARY OEGRID_LIBRARY
     INIT_PY ${CMAKE_SOURCE_DIR}/python/oehealth/__init__.py
 )
```

- [ ] **Step 6: Update `pyproject.toml`**

mypy's `files` list leaves out `scripts/generate_engh_huber.py` for now, because mypy exits with an error for a listed file that does not exist. Task 4 adds it.

```diff
diff --git a/pyproject.toml b/pyproject.toml
--- a/pyproject.toml
+++ b/pyproject.toml
@@ -64,7 +64,7 @@ rpath-strategy = "platform"  # "platform" = @loader_path/../openeye/libs/{platfo
 # (they are provided at runtime by openeye-toolkits, not bundled in the wheel)
 expected-missing-libs = [
     "liboechem", "liboesystem", "liboeplatform",
-    "liboemath", "libzstd",
+    "liboemath", "libzstd", "liboebio", "liboegrid",
 ]
 
 # Additional CMake defines passed during wheel build (beyond the standard ones)
@@ -87,3 +87,52 @@ exclude = [
     "**/*.egg-info",
 ]
 
+
+# =============================================================================
+# Ruff (lint) and mypy (type check)
+# =============================================================================
+# Both cover the package, the reference-table generator, and the Python tests.
+# The SWIG-generated module, the generated build info, and the shared
+# scripts/build_python.py are excluded.
+[tool.ruff]
+line-length = 100
+target-version = "py312"
+include = [
+    "python/oehealth/**/*.py",
+    "scripts/generate_engh_huber.py",
+    "tests/python/**/*.py",
+]
+extend-exclude = [
+    "python/oehealth/oehealth.py",
+    "python/oehealth/_build_info.py",
+    "scripts/build_python.py",
+]
+force-exclude = true
+
+[tool.ruff.lint]
+select = ["E", "F", "W", "I", "UP", "B", "RUF"]
+
+[tool.ruff.lint.isort]
+known-first-party = ["oehealth"]
+
+[tool.ruff.lint.per-file-ignores]
+# The extension import must follow the runtime-library setup calls, and the
+# template's warnings keep their default stack level.
+"python/oehealth/__init__.py" = ["E402", "B028"]
+
+[tool.mypy]
+python_version = "3.12"
+mypy_path = "python"
+files = ["python/oehealth", "tests/python"]
+exclude = ['^python/oehealth/oehealth\.py$', '^python/oehealth/_build_info\.py$']
+warn_unused_ignores = true
+warn_redundant_casts = true
+
+[[tool.mypy.overrides]]
+module = ["oehealth.oehealth", "oehealth._build_info", "openeye", "openeye.*"]
+ignore_missing_imports = true
+follow_imports = "skip"
+
+[[tool.mypy.overrides]]
+module = ["generate_engh_huber", "test_generate_engh_huber"]
+disallow_untyped_defs = true
```

- [ ] **Step 7: Download the crambin fixture**

```bash
mkdir -p tests/data
curl -sSfL -o tests/data/1crn.pdb https://files.rcsb.org/download/1CRN.pdb
shasum -a 256 tests/data/1crn.pdb
```

Expected: `42199a30a0701864a2a5cc76cd7f35cc544cd0e65fbcf63e03c166543249b811  tests/data/1crn.pdb`. If the checksum differs, stop and report it: later tests count crambin's atoms and residues and compare against the wwPDB validation report.

- [ ] **Step 8: Configure, build, and run the existing tests**

```bash
cmake --preset debug
cmake --build build-debug -j8
ctest --test-dir build-debug --output-on-failure
PYTHONPATH=python python -m pytest tests/python -q
```

Expected: the build succeeds; `100% tests passed out of 4` (the template's tests); pytest reports `5 passed`. The configure step prints a CMake author warning from cmake-openeye's `OpenEyeSWIG.cmake` about `_oehealth.py`; it predates this work and needs no change.

- [ ] **Step 9: Run the Python linters**

```bash
ruff check .
mypy --python-executable "$(command -v python)"
```

Expected: `All checks passed!` and `Success: no issues found in 3 source files`. Do not run `ruff format --check tests/python` yet: the template's `test_oehealth.py` is not ruff-formatted, and Task 12 replaces it.

- [ ] **Step 10: Commit**

```bash
git add .clang-format .clang-tidy CMakeLists.txt tests/cpp/CMakeLists.txt swig/CMakeLists.txt pyproject.toml tests/data/1crn.pdb
git commit -m "Add lint configuration, OEBio and nlohmann/json linking, and the crambin test fixture"
```

---

### Task 2: Result model, structure context, assessment base, and helpers

**Risk tier:** standard — the core types and helpers that every assessment uses, across thirteen new files.

**Depends on:** Task 1

**Files:**
- Create: `include/oehealth/result.h`, `include/oehealth/context.h`, `include/oehealth/assessment.h`
- Create: `src/result.cpp`, `src/context.cpp`, `src/assessment.cpp`
- Create: `src/common.h`, `src/common.cpp`
- Create: `src/geometry/measure.h`, `src/geometry/measure.cpp`
- Test: `tests/cpp/test_helpers.h`, `tests/cpp/test_result.cpp`, `tests/cpp/test_context.cpp`
- Modify: `CMakeLists.txt` (`OEHEALTH_SOURCES`), `tests/cpp/CMakeLists.txt` (`add_executable(oehealth_tests ...)`)

**Interfaces:**
- Consumes: `OEHEALTH_TEST_DATA_DIR` and the `src/` include directory from Task 1.
- Produces (public, namespace `OEHEALTH`):
  - `enum class Severity { INFO, WARNING, ERROR }` and `enum class Status { NOT_APPLICABLE, PASS, WARN, FAIL, ERROR }`.
  - `struct AtomRef { unsigned int atom_idx; std::string atom_name, residue_name, chain_id; int residue_number; std::string insert_code; }`.
  - `struct Metric { std::string name; double value; }`.
  - `struct Finding { std::string code; Severity severity; std::string message; std::vector<AtomRef> atoms; double observed, expected, deviation /* NaN by default */; std::vector<Metric> details; }`.
  - `struct AssessmentResult { std::string assessment; Status status; std::string message; std::vector<Metric> metrics; std::vector<Finding> findings; std::vector<AssessmentResult> children; }`.
  - `Status worst_status(Status a, Status b)`.
  - `struct ResidueInfo { std::string name, chain_id; int number; std::string insert_code; bool is_standard_amino_acid; std::vector<const OEChem::OEAtomBase*> atoms; const OEChem::OEAtomBase* FindAtom(const std::string& atom_name) const; }`.
  - `class StructureContext { explicit StructureContext(const OEChem::OEMolBase&); const OEChem::OEMolBase& GetMolecule() const; const std::vector<ResidueInfo>& GetResidues() const; bool HasHydrogens() const; }`.
  - `constexpr double THRESHOLD_TOLERANCE = 1e-6`.
  - `class Assessment { virtual std::string GetName() const = 0; virtual std::string GetDescription() const = 0; virtual AssessmentResult Assess(const StructureContext&) const = 0; AssessmentResult Run(const OEChem::OEMolBase&) const; }`. `Run` throws `std::invalid_argument("<name>: the molecule has no 3D coordinates")` for a molecule with atoms whose dimension is not 3.
- Produces (private, `src/common.h`, namespace `OEHEALTH`): `std::string trim(const std::string&)`, `bool is_standard_amino_acid(const std::string&)`, `std::string trimmed_atom_name(const OEChem::OEAtomBase*)`, `std::string residue_char(char)` (insertion codes), `std::string chain_id_of(const OEChem::OEResidue&)` (the full chain id from `GetExtChainID()`, trimmed), `AtomRef make_atom_ref(const OEChem::OEAtomBase*)`, `std::string residue_label(const AtomRef&)`, `std::string atom_path(const std::vector<AtomRef>&)`, `std::string format_fixed(double value, int decimals)`, `Status status_from_findings(const std::vector<Finding>&)`, `bool meets_threshold(double value, double threshold)`, `void validate_thresholds(const std::string& what, double warning, double error)`.
- Produces (private, `src/geometry/measure.h`, namespace `OEHEALTH::geometry`): `constexpr double PEPTIDE_BOND_CUTOFF = 2.5`, `std::array<double, 3> coordinates(const OEChem::OEMolBase&, const OEChem::OEAtomBase*)`, `double distance(mol, a, b)`, `double angle_degrees(mol, a, b, c)`, `double torsion_degrees(mol, a, b, c, d)`, `std::optional<std::size_t> peptide_successor(const StructureContext&, std::size_t i)`, `std::optional<std::size_t> peptide_predecessor(const StructureContext&, std::size_t j)`.
- Produces (tests, `tests/cpp/test_helpers.h`, namespace `OEHEALTH::test`): `read_crambin(OEMolBase&)`, `trimmed(const std::string&)`, `find_atom(mol, residue_number, name)`, `add_atom(mol, element, name, residue_name, residue_number, x, y, z)`, `coords(mol, atom)`, `set_coords(mol, atom, xyz)`, `metric(result, name)`, `mentions_atom(finding, atom)`.

**Mirror:** `src/oehealth.cpp:1-11` for implementation file layout; `tests/cpp/test_oehealth.cpp:1-40` for test layout.

- [ ] **Step 1: Write the test helpers**

`tests/cpp/test_helpers.h`:

```cpp
#ifndef OEHEALTH_TESTS_TEST_HELPERS_H
#define OEHEALTH_TESTS_TEST_HELPERS_H

#include "oehealth/result.h"

#include <oechem.h>

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace OEHEALTH::test {

/// Read tests/data/1crn.pdb (crambin) with OpenEye's default PDB reader.
inline void read_crambin(OEChem::OEMolBase& mol) {
    const std::string path = std::string(OEHEALTH_TEST_DATA_DIR) + "/1crn.pdb";
    OEChem::oemolistream stream;
    if (!stream.open(path) || !OEChem::OEReadMolecule(stream, mol)) {
        throw std::runtime_error("cannot read " + path);
    }
}

inline std::string trimmed(const std::string& text) {
    const auto first = text.find_first_not_of(' ');
    if (first == std::string::npos) {
        return "";
    }
    return text.substr(first, text.find_last_not_of(' ') - first + 1);
}

/// Find an atom by residue number and trimmed atom name, or return nullptr.
inline OEChem::OEAtomBase* find_atom(const OEChem::OEMolBase& mol, int residue_number,
                                     const std::string& name) {
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (OEChem::OEAtomGetResidue(atom).GetResidueNumber() == residue_number &&
            trimmed(atom->GetName()) == name) {
            return atom;
        }
    }
    return nullptr;
}

/// Add an atom with PDB-style residue information and coordinates.
inline OEChem::OEAtomBase* add_atom(OEChem::OEMolBase& mol, unsigned int element,
                                    const std::string& name, const std::string& residue_name,
                                    int residue_number, double x, double y, double z) {
    OEChem::OEAtomBase* atom = mol.NewAtom(element);
    atom->SetName(name.c_str());
    OEChem::OEResidue residue;
    residue.SetName(residue_name.c_str());
    residue.SetChainID('A');
    residue.SetResidueNumber(residue_number);
    OEChem::OEAtomSetResidue(atom, residue);
    const std::array<double, 3> xyz{x, y, z};
    mol.SetCoords(atom, xyz.data());
    mol.SetDimension(3);
    return atom;
}

inline std::array<double, 3> coords(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* atom) {
    std::array<double, 3> xyz{0.0, 0.0, 0.0};
    mol.GetCoords(atom, xyz.data());
    return xyz;
}

inline void set_coords(OEChem::OEMolBase& mol, OEChem::OEAtomBase* atom,
                       const std::array<double, 3>& xyz) {
    mol.SetCoords(atom, xyz.data());
}

/// Return a metric's value; throws if the result has no such metric.
inline double metric(const AssessmentResult& result, const std::string& name) {
    for (const Metric& m : result.metrics) {
        if (m.name == name) {
            return m.value;
        }
    }
    throw std::runtime_error("no metric named " + name);
}

/// Return whether any of a finding's atoms has the given atom index.
inline bool mentions_atom(const Finding& finding, const OEChem::OEAtomBase* atom) {
    for (const AtomRef& ref : finding.atoms) {
        if (ref.atom_idx == atom->GetIdx()) {
            return true;
        }
    }
    return false;
}

} // namespace OEHEALTH::test

#endif // OEHEALTH_TESTS_TEST_HELPERS_H
```

- [ ] **Step 2: Write the failing result and helper tests**

`tests/cpp/test_result.cpp`:

```cpp
#include <gtest/gtest.h>

#include "common.h"
#include "oehealth/assessment.h"
#include "oehealth/result.h"
#include "test_helpers.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace OEHEALTH::test {

namespace {

Finding finding_with(Severity severity) {
    Finding finding;
    finding.code = "test";
    finding.severity = severity;
    return finding;
}

} // namespace

TEST(WorstStatus, FollowsDeclarationOrder) {
    EXPECT_EQ(worst_status(Status::NOT_APPLICABLE, Status::PASS), Status::PASS);
    EXPECT_EQ(worst_status(Status::PASS, Status::WARN), Status::WARN);
    EXPECT_EQ(worst_status(Status::FAIL, Status::WARN), Status::FAIL);
    EXPECT_EQ(worst_status(Status::FAIL, Status::ERROR), Status::ERROR);
    EXPECT_EQ(worst_status(Status::ERROR, Status::NOT_APPLICABLE), Status::ERROR);
    EXPECT_EQ(worst_status(Status::PASS, Status::PASS), Status::PASS);
}

TEST(StatusFromFindings, DerivesFromSeverities) {
    EXPECT_EQ(status_from_findings({}), Status::PASS);
    EXPECT_EQ(status_from_findings({finding_with(Severity::INFO)}), Status::PASS);
    EXPECT_EQ(status_from_findings({finding_with(Severity::INFO), finding_with(Severity::WARNING)}),
              Status::WARN);
    EXPECT_EQ(
        status_from_findings({finding_with(Severity::WARNING), finding_with(Severity::ERROR)}),
        Status::FAIL);
}

TEST(MeetsThreshold, IsInclusiveWithinTolerance) {
    EXPECT_TRUE(meets_threshold(0.4, 0.4));
    EXPECT_TRUE(meets_threshold(1.7 + 1.7 - 3.0, 0.4));
    EXPECT_TRUE(meets_threshold(0.4 - 0.5 * THRESHOLD_TOLERANCE, 0.4));
    EXPECT_FALSE(meets_threshold(0.4 - 2.0 * THRESHOLD_TOLERANCE, 0.4));
    EXPECT_FALSE(meets_threshold(0.39, 0.4));
}

TEST(ValidateThresholds, RejectsInvalidPairs) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW(validate_thresholds("x", nan, 1.0), std::invalid_argument);
    EXPECT_THROW(validate_thresholds("x", 1.0, inf), std::invalid_argument);
    EXPECT_THROW(validate_thresholds("x", -1.0, 1.0), std::invalid_argument);
    EXPECT_THROW(validate_thresholds("x", 2.0, 1.0), std::invalid_argument);
    EXPECT_NO_THROW(validate_thresholds("x", 1.0, 1.0));
    EXPECT_NO_THROW(validate_thresholds("x", 0.0, 1.0));
}

TEST(Finding, NumericFieldsDefaultToNaN) {
    const Finding finding;
    EXPECT_TRUE(std::isnan(finding.observed));
    EXPECT_TRUE(std::isnan(finding.expected));
    EXPECT_TRUE(std::isnan(finding.deviation));
}

TEST(AtomRef, CopiesCrambinIdentifiers) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const AtomRef ref = make_atom_ref(find_atom(mol, 5, "CA"));
    EXPECT_EQ(ref.atom_name, "CA");
    EXPECT_EQ(ref.residue_name, "PRO");
    EXPECT_EQ(ref.residue_number, 5);
    EXPECT_EQ(ref.chain_id, "A");
    EXPECT_EQ(ref.insert_code, "");
    EXPECT_EQ(ref.atom_idx, find_atom(mol, 5, "CA")->GetIdx());
}

TEST(AtomPath, RepeatsResidueLabelWhenResidueChanges) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const std::vector<AtomRef> same = {make_atom_ref(find_atom(mol, 5, "CA")),
                                       make_atom_ref(find_atom(mol, 5, "CB"))};
    EXPECT_EQ(atom_path(same), "A PRO 5 CA-CB");
    const std::vector<AtomRef> across = {make_atom_ref(find_atom(mol, 4, "C")),
                                         make_atom_ref(find_atom(mol, 5, "N"))};
    EXPECT_EQ(atom_path(across), "A CYS 4 C-A PRO 5 N");
}

} // namespace OEHEALTH::test
```

- [ ] **Step 3: Write the failing context and measurement tests**

`tests/cpp/test_context.cpp`:

```cpp
#include <gtest/gtest.h>

#include "common.h"
#include "geometry/measure.h"
#include "oehealth/context.h"
#include "test_helpers.h"

#include <cmath>
#include <string>
#include <vector>

namespace OEHEALTH::test {

namespace {

constexpr double DEGREES_PER_RADIAN = 180.0 / 3.14159265358979323846;

std::size_t count_peptide_links(const StructureContext& ctx) {
    std::size_t links = 0;
    for (std::size_t i = 0; i < ctx.GetResidues().size(); ++i) {
        if (geometry::peptide_successor(ctx, i)) {
            ++links;
        }
    }
    return links;
}

void set_chain(OEChem::OEAtomBase* atom, const std::string& chain_id) {
    OEChem::OEResidue residue = OEChem::OEAtomGetResidue(atom);
    residue.SetExtChainID(chain_id);
    OEChem::OEAtomSetResidue(atom, residue);
}

} // namespace

TEST(StructureContext, CrambinHas46StandardResiduesInFileOrder) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const StructureContext ctx(mol);
    const auto& residues = ctx.GetResidues();
    ASSERT_EQ(residues.size(), 46U);
    for (std::size_t i = 0; i < residues.size(); ++i) {
        EXPECT_EQ(residues[i].number, static_cast<int>(i) + 1);
        EXPECT_TRUE(residues[i].is_standard_amino_acid) << residues[i].name;
        EXPECT_EQ(residues[i].chain_id, "A");
        EXPECT_EQ(residues[i].insert_code, "");
    }
    EXPECT_EQ(residues.front().name, "THR");
    EXPECT_EQ(residues.back().name, "ASN");
    std::size_t atom_count = 0;
    for (const ResidueInfo& residue : residues) {
        atom_count += residue.atoms.size();
    }
    EXPECT_EQ(atom_count, 327U);
}

TEST(StructureContext, FindAtomComparesTrimmedNames) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const StructureContext ctx(mol);
    const ResidueInfo& first = ctx.GetResidues().front();
    EXPECT_NE(first.FindAtom("CA"), nullptr);
    EXPECT_EQ(first.FindAtom(" CA "), first.FindAtom("CA"));
    EXPECT_EQ(first.FindAtom("XX"), nullptr);
}

TEST(StructureContext, ReportsExplicitHydrogens) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    EXPECT_FALSE(StructureContext(mol).HasHydrogens());
    add_atom(mol, OEChem::OEElemNo::H, "H", "HOH", 900, 50.0, 50.0, 50.0);
    EXPECT_TRUE(StructureContext(mol).HasHydrogens());
}

TEST(StructureContext, KeepsExtendedChainIds) {
    // mmCIF chain ids can be longer than one character. GetChainID() returns
    // only the first, which would make chains "AA" and "AB" look the same.
    OEChem::OEGraphMol mol;
    const std::vector<OEChem::OEAtomBase*> first = {
        add_atom(mol, OEChem::OEElemNo::N, "N", "ALA", 1, 0.0, 0.0, 0.0),
        add_atom(mol, OEChem::OEElemNo::C, "CA", "ALA", 1, 1.46, 0.0, 0.0),
        add_atom(mol, OEChem::OEElemNo::C, "C", "ALA", 1, 2.0, 1.4, 0.0)};
    const std::vector<OEChem::OEAtomBase*> second = {
        add_atom(mol, OEChem::OEElemNo::N, "N", "ALA", 2, 3.3, 1.4, 0.0),
        add_atom(mol, OEChem::OEElemNo::C, "CA", "ALA", 2, 4.0, 2.6, 0.0),
        add_atom(mol, OEChem::OEElemNo::C, "C", "ALA", 2, 5.5, 2.6, 0.0)};
    for (OEChem::OEAtomBase* atom : first) {
        set_chain(atom, "AA");
    }
    for (OEChem::OEAtomBase* atom : second) {
        set_chain(atom, "AB");
    }
    const StructureContext ctx(mol);
    ASSERT_EQ(ctx.GetResidues().size(), 2U);
    EXPECT_EQ(ctx.GetResidues()[0].chain_id, "AA");
    EXPECT_EQ(ctx.GetResidues()[1].chain_id, "AB");
    EXPECT_EQ(make_atom_ref(second[0]).chain_id, "AB");
    // C(1)-N(2) is 1.3 A, but the residues are in different chains.
    EXPECT_FALSE(geometry::peptide_successor(ctx, 0));
}

TEST(StructureContext, EmptyMoleculeHasNoResidues) {
    const OEChem::OEGraphMol mol;
    const StructureContext ctx(mol);
    EXPECT_TRUE(ctx.GetResidues().empty());
    EXPECT_FALSE(ctx.HasHydrogens());
}

TEST(PeptideLinks, CrambinHas45) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const StructureContext ctx(mol);
    EXPECT_EQ(count_peptide_links(ctx), 45U);
    EXPECT_EQ(geometry::peptide_successor(ctx, 0), std::optional<std::size_t>(1));
    EXPECT_FALSE(geometry::peptide_successor(ctx, 45));
    EXPECT_EQ(geometry::peptide_predecessor(ctx, 1), std::optional<std::size_t>(0));
    EXPECT_FALSE(geometry::peptide_predecessor(ctx, 0));
}

TEST(PeptideLinks, DeletingNitrogenBreaksTheLink) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    mol.DeleteAtom(find_atom(mol, 10, "N"));
    const StructureContext ctx(mol);
    // Residue 9 is at index 8 and residue 10 at index 9.
    EXPECT_FALSE(geometry::peptide_successor(ctx, 8));
    EXPECT_FALSE(geometry::peptide_predecessor(ctx, 9));
    EXPECT_EQ(count_peptide_links(ctx), 44U);
}

TEST(Measure, MatchesOpenEyeGeometry) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const OEChem::OEAtomBase* c9 = find_atom(mol, 9, "C");
    const OEChem::OEAtomBase* n10 = find_atom(mol, 10, "N");
    const OEChem::OEAtomBase* ca10 = find_atom(mol, 10, "CA");
    const OEChem::OEAtomBase* c10 = find_atom(mol, 10, "C");
    EXPECT_NEAR(geometry::distance(mol, n10, ca10), OEChem::OEGetDistance(mol, n10, ca10), 1e-5);
    EXPECT_NEAR(geometry::angle_degrees(mol, n10, ca10, c10),
                OEChem::OEGetAngle(mol, n10, ca10, c10) * DEGREES_PER_RADIAN, 1e-3);
    EXPECT_NEAR(geometry::torsion_degrees(mol, c9, n10, ca10, c10),
                OEChem::OEGetTorsion(mol, c9, n10, ca10, c10) * DEGREES_PER_RADIAN, 1e-3);
}

TEST(Measure, UsesDoublePrecisionCoordinates) {
    OEChem::OEMol mol(OEChem::OEMCMolType::DoubleCartesian);
    const double length = 1.520 + 4.0 * 0.021;
    const OEChem::OEAtomBase* a = add_atom(mol, OEChem::OEElemNo::C, "CA", "ALA", 1, 0.0, 0.0, 0.0);
    const OEChem::OEAtomBase* b =
        add_atom(mol, OEChem::OEElemNo::C, "CB", "ALA", 1, length, 0.0, 0.0);
    EXPECT_DOUBLE_EQ(geometry::distance(mol, a, b), length);
}

} // namespace OEHEALTH::test
```

- [ ] **Step 4: Register the tests and confirm they fail**

Add these entries to `add_executable(oehealth_tests ...)` in `tests/cpp/CMakeLists.txt`, keeping the list sorted:

```cmake
    test_context.cpp
    test_result.cpp
```

Run: `cmake --build build-debug -j8`
Expected: FAIL to compile with `'oehealth/result.h' file not found` (or `'common.h' file not found`).

- [ ] **Step 5: Write the result model**

`include/oehealth/result.h`:

```cpp
#ifndef OEHEALTH_RESULT_H
#define OEHEALTH_RESULT_H

#include <limits>
#include <string>
#include <vector>

namespace OEHEALTH {

/// \brief How serious a single finding is.
enum class Severity { INFO, WARNING, ERROR };

/// \brief Overall outcome of an assessment, from best to worst.
///
/// ERROR is reserved for a suite child that threw; it is never derived from
/// findings.
enum class Status { NOT_APPLICABLE, PASS, WARN, FAIL, ERROR };

/// \brief Identifies an atom by copied values, so results outlive the molecule.
struct AtomRef {
    unsigned int atom_idx = 0;
    std::string atom_name;
    std::string residue_name;
    std::string chain_id;
    int residue_number = 0;
    std::string insert_code; ///< Empty when the residue has none.
};

/// \brief A named number. Names are unique within the vector that holds them.
struct Metric {
    std::string name;
    double value = 0.0;
};

/// \brief One problem, or one noteworthy observation, at specific atoms.
///
/// ``observed``, ``expected``, and ``deviation`` are NaN when not meaningful
/// for the finding's code.
struct Finding {
    std::string code; ///< Stable identifier, e.g. "bond_length_outlier".
    Severity severity = Severity::INFO;
    std::string message; ///< Human-readable; names the residue and atoms.
    std::vector<AtomRef> atoms;
    double observed = std::numeric_limits<double>::quiet_NaN();
    double expected = std::numeric_limits<double>::quiet_NaN();
    double deviation = std::numeric_limits<double>::quiet_NaN();
    std::vector<Metric> details; ///< Extra numbers, e.g. phi and psi.
};

/// \brief The outcome of running one assessment or suite.
struct AssessmentResult {
    std::string assessment; ///< GetName() of the producer.
    Status status = Status::NOT_APPLICABLE;
    std::string message; ///< Set for ERROR and NOT_APPLICABLE.
    std::vector<Metric> metrics;
    std::vector<Finding> findings;
    std::vector<AssessmentResult> children; ///< Filled only by suites.
};

/// \brief Return the worse of two statuses, by declaration order.
Status worst_status(Status a, Status b);

} // namespace OEHEALTH

#endif // OEHEALTH_RESULT_H
```

`src/result.cpp`:

```cpp
#include "oehealth/result.h"

namespace OEHEALTH {

Status worst_status(Status a, Status b) {
    return static_cast<int>(a) >= static_cast<int>(b) ? a : b;
}

} // namespace OEHEALTH
```

- [ ] **Step 6: Write the structure context**

`include/oehealth/context.h`:

```cpp
#ifndef OEHEALTH_CONTEXT_H
#define OEHEALTH_CONTEXT_H

#include <oechem.h>

#include <memory>
#include <string>
#include <vector>

namespace OEHEALTH {

/// \brief One residue of the structure, built from OEChem residue data.
struct ResidueInfo {
    std::string name; ///< Residue name, e.g. "ALA".
    std::string chain_id;
    int number = 0;
    std::string insert_code;                      ///< Empty when the residue has none.
    bool is_standard_amino_acid = false;          ///< One of the 20 standard residues.
    std::vector<const OEChem::OEAtomBase*> atoms; ///< In molecule order.

    /// \brief Find an atom by name, compared after trimming whitespace.
    /// \returns The atom, or nullptr when the residue has no such atom.
    const OEChem::OEAtomBase* FindAtom(const std::string& atom_name) const;
};

/// \brief Shared, lazily computed views of one molecule for assessments.
///
/// The context borrows the molecule: the caller keeps it alive and unchanged
/// for the context's lifetime. Lazily computed members are cached on first
/// use, so a single context must not be used from several threads at once.
class StructureContext {
public:
    explicit StructureContext(const OEChem::OEMolBase& mol);

    const OEChem::OEMolBase& GetMolecule() const;

    /// \brief Residues ordered by the index of their first atom (file order for PDB input).
    const std::vector<ResidueInfo>& GetResidues() const;

    /// \brief Whether the molecule contains explicit hydrogen atoms.
    bool HasHydrogens() const;

private:
    const OEChem::OEMolBase* mol_;
    bool has_hydrogens_ = false;
    mutable std::unique_ptr<std::vector<ResidueInfo>> residues_;
};

} // namespace OEHEALTH

#endif // OEHEALTH_CONTEXT_H
```

`src/context.cpp`:

```cpp
#include "oehealth/context.h"

#include "common.h"

#include <map>
#include <memory>
#include <tuple>

namespace OEHEALTH {

const OEChem::OEAtomBase* ResidueInfo::FindAtom(const std::string& atom_name) const {
    const std::string wanted = trim(atom_name);
    for (const OEChem::OEAtomBase* atom : atoms) {
        if (trimmed_atom_name(atom) == wanted) {
            return atom;
        }
    }
    return nullptr;
}

StructureContext::StructureContext(const OEChem::OEMolBase& mol) : mol_(&mol) {
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->IsHydrogen()) {
            has_hydrogens_ = true;
            break;
        }
    }
}

const OEChem::OEMolBase& StructureContext::GetMolecule() const { return *mol_; }

bool StructureContext::HasHydrogens() const { return has_hydrogens_; }

const std::vector<ResidueInfo>& StructureContext::GetResidues() const {
    if (residues_) {
        return *residues_;
    }

    // Visit atoms in index order so residues are ordered by their first atom.
    std::vector<const OEChem::OEAtomBase*> by_index(mol_->GetMaxAtomIdx(), nullptr);
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol_->GetAtoms(); atom; ++atom) {
        by_index[atom->GetIdx()] = atom;
    }

    // Atoms of one residue are usually contiguous, but not always, so group by
    // OEChem residue identity. The key narrows the OESameResidue comparisons.
    using Key = std::tuple<std::string, int, char, std::string>;
    std::map<Key, std::vector<std::size_t>> groups_by_key;
    std::vector<OEChem::OEResidue> group_residues;
    auto residues = std::make_unique<std::vector<ResidueInfo>>();

    for (const OEChem::OEAtomBase* atom : by_index) {
        if (atom == nullptr) {
            continue;
        }
        const OEChem::OEResidue& residue = OEChem::OEAtomGetResidue(atom);
        const std::string chain_id = chain_id_of(residue);
        const Key key{chain_id, residue.GetResidueNumber(), residue.GetInsertCode(),
                      residue.GetName()};
        std::vector<std::size_t>& candidates = groups_by_key[key];
        std::size_t group = residues->size();
        for (std::size_t candidate : candidates) {
            if (OEChem::OESameResidue(group_residues[candidate], residue)) {
                group = candidate;
                break;
            }
        }
        if (group == residues->size()) {
            ResidueInfo info;
            info.name = trim(residue.GetName());
            info.chain_id = chain_id;
            info.number = residue.GetResidueNumber();
            info.insert_code = residue_char(residue.GetInsertCode());
            info.is_standard_amino_acid = is_standard_amino_acid(info.name);
            residues->push_back(info);
            group_residues.push_back(residue);
            candidates.push_back(group);
        }
        (*residues)[group].atoms.push_back(atom);
    }

    residues_ = std::move(residues);
    return *residues_;
}

} // namespace OEHEALTH
```

- [ ] **Step 7: Write the assessment base**

`include/oehealth/assessment.h`:

```cpp
#ifndef OEHEALTH_ASSESSMENT_H
#define OEHEALTH_ASSESSMENT_H

#include "oehealth/context.h"
#include "oehealth/result.h"

#include <oechem.h>

#include <string>

namespace OEHEALTH {

/// \brief Tolerance for inclusive threshold comparisons.
///
/// A value meets a threshold t when value >= t - THRESHOLD_TOLERANCE, so that
/// floating-point rounding does not decide boundary cases (1.7 + 1.7 - 3.0
/// evaluates to 0.39999..., just under a 0.4 threshold).
constexpr double THRESHOLD_TOLERANCE = 1e-6;

/// \brief A single check, or a named set of checks, over one structure.
///
/// Assessments are immutable after construction, so one instance may be
/// shared between suites and threads.
class Assessment {
public:
    virtual ~Assessment() = default;

    /// \brief Stable identifier, e.g. "bond_lengths".
    virtual std::string GetName() const = 0;

    /// \brief One-sentence description.
    virtual std::string GetDescription() const = 0;

    /// \brief Assess the structure behind a context.
    virtual AssessmentResult Assess(const StructureContext& ctx) const = 0;

    /// \brief Validate a molecule, build a context, and assess it.
    ///
    /// A molecule with no atoms is passed through, and each assessment
    /// reports NOT_APPLICABLE. For a multi-conformer molecule the active
    /// conformer is assessed.
    ///
    /// \throws std::invalid_argument if the molecule has atoms but its
    ///         dimension is not 3.
    AssessmentResult Run(const OEChem::OEMolBase& mol) const;
};

} // namespace OEHEALTH

#endif // OEHEALTH_ASSESSMENT_H
```

`src/assessment.cpp`:

```cpp
#include "oehealth/assessment.h"

#include <stdexcept>

namespace OEHEALTH {

AssessmentResult Assessment::Run(const OEChem::OEMolBase& mol) const {
    if (mol.NumAtoms() > 0 && mol.GetDimension() != 3) {
        throw std::invalid_argument(GetName() + ": the molecule has no 3D coordinates");
    }
    const StructureContext ctx(mol);
    return Assess(ctx);
}

} // namespace OEHEALTH
```

- [ ] **Step 8: Write the common helpers**

`src/common.h`:

```cpp
#ifndef OEHEALTH_SRC_COMMON_H
#define OEHEALTH_SRC_COMMON_H

#include "oehealth/result.h"

#include <oechem.h>

#include <string>
#include <vector>

namespace OEHEALTH {

/// \brief Return text with leading and trailing whitespace removed.
std::string trim(const std::string& text);

/// \brief Return whether a residue name is one of the 20 standard amino acids.
bool is_standard_amino_acid(const std::string& residue_name);

/// \brief Return an atom's name with surrounding whitespace removed.
std::string trimmed_atom_name(const OEChem::OEAtomBase* atom);

/// \brief Convert an OEChem insertion code character to a string.
///
/// A blank (space or NUL) becomes the empty string.
std::string residue_char(char value);

/// \brief Return a residue's full chain id, trimmed; a blank id becomes "".
///
/// Uses GetExtChainID(): GetChainID() returns only the first character, so
/// mmCIF chains such as "AA" and "AB" would collide.
std::string chain_id_of(const OEChem::OEResidue& residue);

/// \brief Copy an atom's identifiers into an AtomRef.
AtomRef make_atom_ref(const OEChem::OEAtomBase* atom);

/// \brief Return "chain name number[insert]", e.g. "A ALA 5" or "B GLY 52A".
///
/// The chain is omitted when empty.
std::string residue_label(const AtomRef& atom);

/// \brief Describe atoms for a message, e.g. "A ALA 5 CA-CB" or "A ALA 5 C-A GLY 6 N".
///
/// The residue label is repeated whenever the residue changes along the path.
std::string atom_path(const std::vector<AtomRef>& atoms);

/// \brief Format a number with a fixed number of decimals.
std::string format_fixed(double value, int decimals);

/// \brief Derive an assessment status from its findings.
///
/// Any ERROR finding gives FAIL, otherwise any WARNING gives WARN, otherwise
/// PASS. INFO findings never affect status.
Status status_from_findings(const std::vector<Finding>& findings);

/// \brief Whether a value meets an inclusive threshold, within THRESHOLD_TOLERANCE.
bool meets_threshold(double value, double threshold);

/// \brief Throw std::invalid_argument unless a warning/error threshold pair is valid.
///
/// Both must be finite and non-negative, and error must not be less than
/// warning. ``what`` names the option pair in the exception message.
void validate_thresholds(const std::string& what, double warning, double error);

} // namespace OEHEALTH

#endif // OEHEALTH_SRC_COMMON_H
```

`src/common.cpp`:

```cpp
#include "common.h"

#include "oehealth/assessment.h"

#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace OEHEALTH {

std::string trim(const std::string& text) {
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

bool is_standard_amino_acid(const std::string& residue_name) {
    static const std::array<const char*, 20> STANDARD = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
    for (const char* name : STANDARD) {
        if (residue_name == name) {
            return true;
        }
    }
    return false;
}

std::string trimmed_atom_name(const OEChem::OEAtomBase* atom) { return trim(atom->GetName()); }

std::string residue_char(char value) {
    if (value == ' ' || value == '\0') {
        return "";
    }
    return {value};
}

std::string chain_id_of(const OEChem::OEResidue& residue) { return trim(residue.GetExtChainID()); }

AtomRef make_atom_ref(const OEChem::OEAtomBase* atom) {
    const OEChem::OEResidue& residue = OEChem::OEAtomGetResidue(atom);
    AtomRef ref;
    ref.atom_idx = atom->GetIdx();
    ref.atom_name = trimmed_atom_name(atom);
    ref.residue_name = trim(residue.GetName());
    ref.chain_id = chain_id_of(residue);
    ref.residue_number = residue.GetResidueNumber();
    ref.insert_code = residue_char(residue.GetInsertCode());
    return ref;
}

std::string residue_label(const AtomRef& atom) {
    std::string label;
    if (!atom.chain_id.empty()) {
        label = atom.chain_id + " ";
    }
    label += atom.residue_name + " " + std::to_string(atom.residue_number) + atom.insert_code;
    return label;
}

namespace {

bool same_residue(const AtomRef& a, const AtomRef& b) {
    return a.chain_id == b.chain_id && a.residue_name == b.residue_name &&
           a.residue_number == b.residue_number && a.insert_code == b.insert_code;
}

} // namespace

std::string atom_path(const std::vector<AtomRef>& atoms) {
    std::string path;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        if (i > 0) {
            path += '-';
        }
        if (i == 0 || !same_residue(atoms[i - 1], atoms[i])) {
            path += residue_label(atoms[i]) + " ";
        }
        path += atoms[i].atom_name;
    }
    return path;
}

std::string format_fixed(double value, int decimals) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(decimals) << value;
    return stream.str();
}

Status status_from_findings(const std::vector<Finding>& findings) {
    Status status = Status::PASS;
    for (const Finding& finding : findings) {
        if (finding.severity == Severity::ERROR) {
            return Status::FAIL;
        }
        if (finding.severity == Severity::WARNING) {
            status = Status::WARN;
        }
    }
    return status;
}

bool meets_threshold(double value, double threshold) {
    return value >= threshold - THRESHOLD_TOLERANCE;
}

void validate_thresholds(const std::string& what, double warning, double error) {
    if (!std::isfinite(warning) || !std::isfinite(error)) {
        throw std::invalid_argument(what + " thresholds must be finite");
    }
    if (warning < 0.0 || error < 0.0) {
        throw std::invalid_argument(what + " thresholds must not be negative");
    }
    if (error < warning) {
        throw std::invalid_argument(what +
                                    " error threshold must not be less than its warning threshold");
    }
}

} // namespace OEHEALTH
```

- [ ] **Step 9: Write the measurement helpers**

```bash
mkdir -p src/geometry
```

`src/geometry/measure.h`:

```cpp
#ifndef OEHEALTH_SRC_GEOMETRY_MEASURE_H
#define OEHEALTH_SRC_GEOMETRY_MEASURE_H

#include "oehealth/context.h"

#include <oechem.h>

#include <array>
#include <cstddef>
#include <optional>

namespace OEHEALTH::geometry {

/// \brief C(i)-N(j) distances at or above this, in angstroms, are chain breaks.
constexpr double PEPTIDE_BOND_CUTOFF = 2.5;

/// \brief Coordinates of an atom in the active conformer, in double precision.
///
/// OEGetDistance and friends compute in single precision; measuring from
/// double coordinates keeps threshold boundaries exact for molecules that
/// store doubles (OEMCMolType::DoubleCartesian).
std::array<double, 3> coordinates(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* atom);

/// \brief Distance between two atoms, in angstroms.
double distance(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* a,
                const OEChem::OEAtomBase* b);

/// \brief Angle a-b-c with b as the vertex, in degrees.
double angle_degrees(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* a,
                     const OEChem::OEAtomBase* b, const OEChem::OEAtomBase* c);

/// \brief Torsion a-b-c-d in degrees, in (-180, 180], IUPAC sign convention.
double torsion_degrees(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* a,
                       const OEChem::OEAtomBase* b, const OEChem::OEAtomBase* c,
                       const OEChem::OEAtomBase* d);

/// \brief The peptide successor of residue i, if it has one.
///
/// Residue j is the successor when it is the next entry in the residue
/// table, both are standard amino acids with the same chain id, C(i) and
/// N(j) are present, and their distance is under PEPTIDE_BOND_CUTOFF.
std::optional<std::size_t> peptide_successor(const StructureContext& ctx, std::size_t i);

/// \brief The peptide predecessor of residue j, if it has one.
std::optional<std::size_t> peptide_predecessor(const StructureContext& ctx, std::size_t j);

} // namespace OEHEALTH::geometry

#endif // OEHEALTH_SRC_GEOMETRY_MEASURE_H
```

`src/geometry/measure.cpp`:

```cpp
#include "geometry/measure.h"

#include <cmath>

namespace OEHEALTH::geometry {

namespace {

using Vec3 = std::array<double, 3>;

constexpr double PI = 3.14159265358979323846;
constexpr double RADIANS_TO_DEGREES = 180.0 / PI;

Vec3 subtract(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }

double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

} // namespace

std::array<double, 3> coordinates(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* atom) {
    Vec3 xyz{0.0, 0.0, 0.0};
    mol.GetCoords(atom, xyz.data());
    return xyz;
}

double distance(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* a,
                const OEChem::OEAtomBase* b) {
    return norm(subtract(coordinates(mol, a), coordinates(mol, b)));
}

double angle_degrees(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* a,
                     const OEChem::OEAtomBase* b, const OEChem::OEAtomBase* c) {
    const Vec3 vertex = coordinates(mol, b);
    const Vec3 u = subtract(coordinates(mol, a), vertex);
    const Vec3 v = subtract(coordinates(mol, c), vertex);
    // atan2 of |u x v| and u.v is accurate near 0 and 180 degrees, unlike acos.
    return std::atan2(norm(cross(u, v)), dot(u, v)) * RADIANS_TO_DEGREES;
}

double torsion_degrees(const OEChem::OEMolBase& mol, const OEChem::OEAtomBase* a,
                       const OEChem::OEAtomBase* b, const OEChem::OEAtomBase* c,
                       const OEChem::OEAtomBase* d) {
    const Vec3 pb = coordinates(mol, b);
    const Vec3 pc = coordinates(mol, c);
    const Vec3 b1 = subtract(pb, coordinates(mol, a));
    const Vec3 b2 = subtract(pc, pb);
    const Vec3 b3 = subtract(coordinates(mol, d), pc);
    const double y = norm(b2) * dot(b1, cross(b2, b3));
    const double x = dot(cross(b1, b2), cross(b2, b3));
    return std::atan2(y, x) * RADIANS_TO_DEGREES;
}

std::optional<std::size_t> peptide_successor(const StructureContext& ctx, std::size_t i) {
    const auto& residues = ctx.GetResidues();
    const std::size_t j = i + 1;
    if (j >= residues.size()) {
        return std::nullopt;
    }
    const ResidueInfo& first = residues[i];
    const ResidueInfo& second = residues[j];
    if (!first.is_standard_amino_acid || !second.is_standard_amino_acid ||
        first.chain_id != second.chain_id) {
        return std::nullopt;
    }
    const OEChem::OEAtomBase* carbon = first.FindAtom("C");
    const OEChem::OEAtomBase* nitrogen = second.FindAtom("N");
    if (carbon == nullptr || nitrogen == nullptr) {
        return std::nullopt;
    }
    if (distance(ctx.GetMolecule(), carbon, nitrogen) >= PEPTIDE_BOND_CUTOFF) {
        return std::nullopt;
    }
    return j;
}

std::optional<std::size_t> peptide_predecessor(const StructureContext& ctx, std::size_t j) {
    if (j == 0) {
        return std::nullopt;
    }
    const auto successor = peptide_successor(ctx, j - 1);
    if (successor && *successor == j) {
        return j - 1;
    }
    return std::nullopt;
}

} // namespace OEHEALTH::geometry
```

- [ ] **Step 10: Register the sources**

Add these entries to `set(OEHEALTH_SOURCES ...)` in `CMakeLists.txt`, keeping the list sorted (`src/oehealth.cpp` stays until Task 12):

```cmake
    src/assessment.cpp
    src/common.cpp
    src/context.cpp
    src/geometry/measure.cpp
    src/result.cpp
```

- [ ] **Step 11: Build and run the tests**

```bash
cmake --build build-debug -j8
./build-debug/tests/cpp/oehealth_tests --gtest_filter='WorstStatus.*:StatusFromFindings.*:MeetsThreshold.*:ValidateThresholds.*:Finding.*:AtomRef.*:AtomPath.*:StructureContext.*:PeptideLinks.*:Measure.*'
ctest --test-dir build-debug --output-on-failure
```

Expected: `[  PASSED  ] 16 tests.` from the filtered run, and `100% tests passed out of 20` from ctest.

- [ ] **Step 12: Format and lint**

```bash
clang-format --dry-run --Werror include/oehealth/result.h include/oehealth/context.h include/oehealth/assessment.h src/result.cpp src/context.cpp src/assessment.cpp src/common.h src/common.cpp src/geometry/measure.h src/geometry/measure.cpp tests/cpp/test_helpers.h tests/cpp/test_result.cpp tests/cpp/test_context.cpp
clang-tidy -p build-debug --quiet --extra-arg=-isysroot --extra-arg="$(xcrun --show-sdk-path)" src/result.cpp src/context.cpp src/assessment.cpp src/common.cpp src/geometry/measure.cpp tests/cpp/test_result.cpp tests/cpp/test_context.cpp
```

Expected: no output from clang-format; no findings from clang-tidy.

- [ ] **Step 13: Commit**

```bash
git add include/oehealth/result.h include/oehealth/context.h include/oehealth/assessment.h src/result.cpp src/context.cpp src/assessment.cpp src/common.h src/common.cpp src/geometry/measure.h src/geometry/measure.cpp tests/cpp/test_helpers.h tests/cpp/test_result.cpp tests/cpp/test_context.cpp CMakeLists.txt tests/cpp/CMakeLists.txt
git commit -m "Add the assessment result model, structure context, and assessment base"
```

---

### Task 3: Assessment suites

**Risk tier:** standard — new public class with exception-isolation behavior, three files plus CMake registration.

**Depends on:** Task 2

**Files:**
- Create: `include/oehealth/suite.h`, `src/suite.cpp`
- Test: `tests/cpp/test_suite.cpp`
- Modify: `CMakeLists.txt` (`OEHEALTH_SOURCES`), `tests/cpp/CMakeLists.txt`

**Interfaces:**
- Consumes: `Assessment`, `AssessmentResult`, `Status`, `worst_status`, `StructureContext` from Task 2.
- Produces: `class AssessmentSuite : public Assessment { AssessmentSuite(std::string name, std::string description, std::vector<std::shared_ptr<Assessment>> children); std::string GetName() const override; std::string GetDescription() const override; const std::vector<std::shared_ptr<Assessment>>& GetChildren() const; AssessmentResult Assess(const StructureContext&) const override; }`. The constructor throws `std::invalid_argument` for an empty name or a null child. Task 11 adds `make_protein_geometry_suite()` to these same files.

**Mirror:** `src/oehealth.cpp:1-11`, implementation file layout.

- [ ] **Step 1: Write the failing tests**

`tests/cpp/test_suite.cpp`:

```cpp
#include <gtest/gtest.h>

#include "oehealth/assessment.h"
#include "oehealth/suite.h"
#include "test_helpers.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace OEHEALTH::test {

namespace {

// Returns a fixed status, so suite behavior can be tested in isolation.
class FixedAssessment : public Assessment {
public:
    FixedAssessment(std::string name, Status status) : name_(std::move(name)), status_(status) {}
    std::string GetName() const override { return name_; }
    std::string GetDescription() const override { return "Returns a fixed status."; }
    AssessmentResult Assess(const StructureContext& /*ctx*/) const override {
        AssessmentResult result;
        result.assessment = name_;
        result.status = status_;
        return result;
    }

private:
    std::string name_;
    Status status_;
};

class ThrowingAssessment : public Assessment {
public:
    std::string GetName() const override { return "thrower"; }
    std::string GetDescription() const override { return "Always throws."; }
    AssessmentResult Assess(const StructureContext& /*ctx*/) const override {
        throw std::runtime_error("boom");
    }
};

// Records the number of atoms it saw, to check what Run passes through.
class CountingAssessment : public Assessment {
public:
    std::string GetName() const override { return "counter"; }
    std::string GetDescription() const override { return "Counts atoms."; }
    AssessmentResult Assess(const StructureContext& ctx) const override {
        AssessmentResult result;
        result.assessment = GetName();
        result.status = ctx.GetMolecule().NumAtoms() == 0 ? Status::NOT_APPLICABLE : Status::PASS;
        result.metrics = {{"n_atoms", static_cast<double>(ctx.GetMolecule().NumAtoms())}};
        return result;
    }
};

AssessmentResult assess(const Assessment& assessment) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    return assessment.Run(mol);
}

} // namespace

TEST(AssessmentSuite, ReturnsConstructorValues) {
    auto child = std::make_shared<FixedAssessment>("a", Status::PASS);
    const AssessmentSuite suite("my_suite", "My description.", {child});
    EXPECT_EQ(suite.GetName(), "my_suite");
    EXPECT_EQ(suite.GetDescription(), "My description.");
    ASSERT_EQ(suite.GetChildren().size(), 1U);
    EXPECT_EQ(suite.GetChildren()[0], child);
}

TEST(AssessmentSuite, RunsChildrenInOrderAndTakesWorstStatus) {
    const AssessmentSuite suite("s", "",
                                {std::make_shared<FixedAssessment>("a", Status::PASS),
                                 std::make_shared<FixedAssessment>("b", Status::FAIL),
                                 std::make_shared<FixedAssessment>("c", Status::WARN)});
    const AssessmentResult result = assess(suite);
    EXPECT_EQ(result.assessment, "s");
    EXPECT_EQ(result.status, Status::FAIL);
    EXPECT_TRUE(result.message.empty());
    EXPECT_TRUE(result.metrics.empty());
    EXPECT_TRUE(result.findings.empty());
    ASSERT_EQ(result.children.size(), 3U);
    EXPECT_EQ(result.children[0].assessment, "a");
    EXPECT_EQ(result.children[1].assessment, "b");
    EXPECT_EQ(result.children[2].assessment, "c");
}

TEST(AssessmentSuite, RecordsThrowingChildAndRunsTheRest) {
    const AssessmentSuite suite("s", "",
                                {std::make_shared<ThrowingAssessment>(),
                                 std::make_shared<FixedAssessment>("after", Status::PASS)});
    const AssessmentResult result = assess(suite);
    ASSERT_EQ(result.children.size(), 2U);
    const AssessmentResult& failed = result.children[0];
    EXPECT_EQ(failed.assessment, "thrower");
    EXPECT_EQ(failed.status, Status::ERROR);
    EXPECT_EQ(failed.message, "boom");
    EXPECT_TRUE(failed.metrics.empty());
    EXPECT_TRUE(failed.findings.empty());
    EXPECT_TRUE(failed.children.empty());
    EXPECT_EQ(result.children[1].status, Status::PASS);
    EXPECT_EQ(result.status, Status::ERROR);
}

TEST(AssessmentSuite, EmptySuiteIsNotApplicable) {
    const AssessmentSuite suite("empty", "", {});
    const AssessmentResult result = assess(suite);
    EXPECT_EQ(result.status, Status::NOT_APPLICABLE);
    EXPECT_FALSE(result.message.empty());
    EXPECT_TRUE(result.children.empty());
}

TEST(AssessmentSuite, SuitesNest) {
    auto inner = std::make_shared<AssessmentSuite>(
        "inner", "",
        std::vector<std::shared_ptr<Assessment>>{
            std::make_shared<FixedAssessment>("x", Status::WARN)});
    const AssessmentSuite outer("outer", "", {inner});
    const AssessmentResult result = assess(outer);
    EXPECT_EQ(result.status, Status::WARN);
    ASSERT_EQ(result.children.size(), 1U);
    ASSERT_EQ(result.children[0].children.size(), 1U);
    EXPECT_EQ(result.children[0].children[0].assessment, "x");
}

TEST(AssessmentSuite, RejectsNullChildAndEmptyName) {
    EXPECT_THROW(AssessmentSuite("s", "", {nullptr}), std::invalid_argument);
    EXPECT_THROW(AssessmentSuite("", "", {}), std::invalid_argument);
}

TEST(AssessmentRun, Rejects2DMoleculeWithAtoms) {
    OEChem::OEGraphMol mol;
    OEChem::OESmilesToMol(mol, "CCO");
    ASSERT_NE(mol.GetDimension(), 3U);
    EXPECT_THROW(CountingAssessment().Run(mol), std::invalid_argument);
}

TEST(AssessmentRun, PassesEmptyMoleculeThrough) {
    const OEChem::OEGraphMol mol;
    const AssessmentResult result = CountingAssessment().Run(mol);
    EXPECT_EQ(result.status, Status::NOT_APPLICABLE);
    EXPECT_EQ(metric(result, "n_atoms"), 0.0);
}

} // namespace OEHEALTH::test
```

- [ ] **Step 2: Register the test and confirm it fails**

Add to `add_executable(oehealth_tests ...)` in `tests/cpp/CMakeLists.txt`, keeping the list sorted:

```cmake
    test_suite.cpp
```

Run: `cmake --build build-debug -j8`
Expected: FAIL to compile with `'oehealth/suite.h' file not found`.

- [ ] **Step 3: Write the suite**

`include/oehealth/suite.h`:

```cpp
#ifndef OEHEALTH_SUITE_H
#define OEHEALTH_SUITE_H

#include "oehealth/assessment.h"

#include <memory>
#include <string>
#include <vector>

namespace OEHEALTH {

/// \brief A named, ordered set of assessments that is itself an assessment.
///
/// Children are fixed at construction, so a suite cannot contain itself.
class AssessmentSuite : public Assessment {
public:
    /// \throws std::invalid_argument if the name is empty or any child is null.
    AssessmentSuite(std::string name, std::string description,
                    std::vector<std::shared_ptr<Assessment>> children);

    std::string GetName() const override;
    std::string GetDescription() const override;
    const std::vector<std::shared_ptr<Assessment>>& GetChildren() const;

    /// \brief Run every child in order against the same context.
    ///
    /// A child that throws a std::exception is recorded as a child result
    /// with status ERROR and the exception message; the remaining children
    /// still run. The suite's status is the worst child status, or
    /// NOT_APPLICABLE when it has no children.
    AssessmentResult Assess(const StructureContext& ctx) const override;

private:
    std::string name_;
    std::string description_;
    std::vector<std::shared_ptr<Assessment>> children_;
};

} // namespace OEHEALTH

#endif // OEHEALTH_SUITE_H
```

`src/suite.cpp`:

```cpp
#include "oehealth/suite.h"

#include <exception>
#include <stdexcept>
#include <utility>

namespace OEHEALTH {

AssessmentSuite::AssessmentSuite(std::string name, std::string description,
                                 std::vector<std::shared_ptr<Assessment>> children)
    : name_(std::move(name)), description_(std::move(description)), children_(std::move(children)) {
    if (name_.empty()) {
        throw std::invalid_argument("an assessment suite needs a name");
    }
    for (const auto& child : children_) {
        if (!child) {
            throw std::invalid_argument("assessment suite '" + name_ + "' has a null child");
        }
    }
}

std::string AssessmentSuite::GetName() const { return name_; }

std::string AssessmentSuite::GetDescription() const { return description_; }

const std::vector<std::shared_ptr<Assessment>>& AssessmentSuite::GetChildren() const {
    return children_;
}

AssessmentResult AssessmentSuite::Assess(const StructureContext& ctx) const {
    AssessmentResult result;
    result.assessment = name_;
    if (children_.empty()) {
        result.status = Status::NOT_APPLICABLE;
        result.message = "the suite has no assessments";
        return result;
    }
    result.status = Status::NOT_APPLICABLE;
    for (const auto& child : children_) {
        AssessmentResult child_result;
        try {
            child_result = child->Assess(ctx);
        } catch (const std::exception& error) {
            child_result = AssessmentResult();
            child_result.assessment = child->GetName();
            child_result.status = Status::ERROR;
            child_result.message = error.what();
        }
        result.status = worst_status(result.status, child_result.status);
        result.children.push_back(std::move(child_result));
    }
    return result;
}

} // namespace OEHEALTH
```

- [ ] **Step 4: Register the source**

Add to `set(OEHEALTH_SOURCES ...)` in `CMakeLists.txt`, keeping the list sorted:

```cmake
    src/suite.cpp
```

- [ ] **Step 5: Build and run the tests**

```bash
cmake --build build-debug -j8
./build-debug/tests/cpp/oehealth_tests --gtest_filter='AssessmentSuite.*:AssessmentRun.*'
ctest --test-dir build-debug --output-on-failure
```

Expected: `[  PASSED  ] 8 tests.` from the filtered run; ctest passes every test (28 after Tasks 1-3).

- [ ] **Step 6: Format and lint**

```bash
clang-format --dry-run --Werror include/oehealth/suite.h src/suite.cpp tests/cpp/test_suite.cpp
clang-tidy -p build-debug --quiet --extra-arg=-isysroot --extra-arg="$(xcrun --show-sdk-path)" src/suite.cpp tests/cpp/test_suite.cpp
```

Expected: no output from clang-format; no findings from clang-tidy.

- [ ] **Step 7: Commit**

```bash
git add include/oehealth/suite.h src/suite.cpp tests/cpp/test_suite.cpp CMakeLists.txt tests/cpp/CMakeLists.txt
git commit -m "Add assessment suites"
```

---

### Task 4: Engh & Huber (2001) reference tables and their generator

**Risk tier:** high — a new script whose generated tables every geometry target trusts without further review, plus tests in two languages.

**Depends on:** Task 1

**Files:**
- Create: `scripts/generate_engh_huber.py`
- Create: `tests/data/standard_geometry_excerpt.cif`
- Create: `src/data/engh_huber_2001.h`
- Create: `src/data/engh_huber_2001.cpp` (generated by the script; do not type it in)
- Test: `tests/python/test_generate_engh_huber.py`, `tests/cpp/test_reference_data.cpp`
- Modify: `pyproject.toml` (mypy `files`), `CMakeLists.txt` (`OEHEALTH_SOURCES`), `tests/cpp/CMakeLists.txt`

**Interfaces:**
- Consumes: the ruff and mypy configuration from Task 1.
- Produces (private, `src/data/engh_huber_2001.h`, namespace `OEHEALTH::data`): `enum class LinkVariant { DEFAULT, GLY, PRO, PRO_CIS, PRO_TRANS, PRO_ALL, DISULFIDE }`; row structs `ResidueBond { residue, atom1, atom2, variant, ideal, sigma }`, `ResidueAngle { residue, atom1, atom2, atom3, variant, ideal, sigma }`, `LinkBond { atom1, atom2, link_variant, ideal, sigma }`, `LinkAngle { atom1, atom2, atom3, link_variant, ideal, sigma }` (names are `const char*`, values `double`); `template <typename Row> struct RowSpan` with `begin()`, `end()`, `size()`; accessors `RowSpan<ResidueBond> residue_bonds()`, `RowSpan<ResidueAngle> residue_angles()`, `RowSpan<LinkBond> link_bonds()`, `RowSpan<LinkAngle> link_angles()`, `const char* source_version()`.
- Produces (script): `python3 scripts/generate_engh_huber.py <standard_geometry.cif> <output.cpp>`; exits non-zero and names the row for an unknown link id on a C–N or SG–SG row.

The source file is CCP4 9's copy of `standard_geometry.cif`, at `$CCP4/lib/data/standard_geometry.cif` (SHA-256 `2e08ceb7998546dbe36fe552fdb243744737fe1fb06f25bc0bbf729e72f0b124`, `_standard_geometry.version "Fri Feb 22 17:25:15 GMT 2013"`). If it is missing, stop and report; do not write the tables by hand.

- [ ] **Step 1: Add the CIF excerpt the generator tests run on**

`tests/data/standard_geometry_excerpt.cif`:

```text
# Excerpt of standard_geometry.cif (NDB/wwPDB, distributed with CCP4) used to
# test scripts/generate_engh_huber.py. Values are copied from the full file.
data_standard_geometry
_standard_geometry.version "Fri Feb 22 17:25:15 GMT 2013"
#
# This loop has more values per row than fields, as in the full file. The
# generator never reads it and must not reject the file because of it.
loop_
_chem_comp.comp_id
_chem_comp.one_letter_code
  ALA  A   1  16
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
ALA   N      N
ALA   CA     C
ALA   C      C
ALA   O      O
ALA   CB     C
ALA   H      H
PRO   N      N
PRO   CA     C
PRO   CD     C
PRO   H2     H
DA    "O5'"  O
DA    "C5'"  C
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.description
ALA    N     CA        1.459    0.020 .
ALA    CA    CB        1.520    0.021 .
ALA    N     H         0.860    0.020 .
# .
PRO    N     CD        1.474    0.014 "trans"
PRO    N     CD        1.474    0.014 "cis"
DA     "O5'" "C5'"     1.440    0.016 .
#
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
_chem_comp_angle.description
 ALA      N      CA     C       111.0      2.7  .
 PRO      CA     N      CD      111.7      1.4  "trans"
 PRO      CA     N      CD      111.5      1.4  "cis"
 PRO      C      N      CD      128.4      2.1  "trans"
 PRO      C      N      CD      120.6      2.2  "cis"
#
loop_
_chem_link_bond.atom_id_1
_chem_link_bond.atom_id_2
_chem_link_bond.value_dist
_chem_link_bond.value_dist_esd
_chem_link_bond.link_id
 C    N        1.336    0.023    "NOT GLY PRO"
 O3'  P        1.607    0.012    "NUC-ACID-ALL"
SG   SG        2.033    0.016    "CYS"
#
loop_
_chem_link_angle.atom_id_1
_chem_link_angle.atom_id_2
_chem_link_angle.atom_id_3
_chem_link_angle.value_angle
_chem_link_angle.value_angle_esd
_chem_link_angle.link_id
     O    C    N      122.7      1.6    "NOT PRO GLY"
     O3'  P    O5'    104.0      1.9    "NA s"
     C    N    CD     128.4      2.1    PRO-TRANS
#
```

- [ ] **Step 2: Write the failing generator tests**

`tests/python/test_generate_engh_huber.py`:

```python
"""Tests for scripts/generate_engh_huber.py, run on a committed excerpt."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "generate_engh_huber.py"
EXCERPT = ROOT / "tests" / "data" / "standard_geometry_excerpt.cif"


def run_generator(source: Path, output: Path) -> subprocess.CompletedProcess[str]:
    """Run the generator as a separate process, as a developer would."""
    return subprocess.run(
        [sys.executable, str(SCRIPT), str(source), str(output)],
        capture_output=True,
        text=True,
        check=False,
    )


def generate_excerpt(tmp_path: Path) -> str:
    """Generate C++ from the excerpt and return its text."""
    output = tmp_path / "engh_huber_2001.cpp"
    result = run_generator(EXCERPT, output)
    assert result.returncode == 0, result.stderr
    return output.read_text(encoding="utf-8")


def test_emits_heavy_atom_residue_rows(tmp_path: Path) -> None:
    text = generate_excerpt(tmp_path)
    assert '{"ALA", "N", "CA", ".", 1.459, 0.020},' in text
    assert '{"ALA", "CA", "CB", ".", 1.520, 0.021},' in text
    assert '{"PRO", "N", "CD", "trans", 1.474, 0.014},' in text
    assert '{"ALA", "N", "CA", "C", ".", 111.0, 2.7},' in text
    assert '{"PRO", "CA", "N", "CD", "cis", 111.5, 1.4},' in text
    assert "constexpr std::array<ResidueBond, 4> RESIDUE_BONDS{{" in text
    assert "constexpr std::array<ResidueAngle, 3> RESIDUE_ANGLES{{" in text


def test_skips_hydrogens_and_nonstandard_residues(tmp_path: Path) -> None:
    text = generate_excerpt(tmp_path)
    assert '"H"' not in text
    assert '"DA"' not in text


def test_drops_intra_residue_rows_with_consecutive_c_and_n(tmp_path: Path) -> None:
    text = generate_excerpt(tmp_path)
    assert '{"PRO", "C", "N", "CD"' not in text


def test_normalizes_link_ids(tmp_path: Path) -> None:
    text = generate_excerpt(tmp_path)
    assert '{"C", "N", LinkVariant::DEFAULT, 1.336, 0.023},' in text
    # "NOT PRO GLY" is the source's other spelling of the default row.
    assert '{"O", "C", "N", LinkVariant::DEFAULT, 122.7, 1.6},' in text
    assert '{"C", "N", "CD", LinkVariant::PRO_TRANS, 128.4, 2.1},' in text
    assert '{"SG", "SG", LinkVariant::DISULFIDE, 2.033, 0.016},' in text


def test_skips_nucleic_acid_links(tmp_path: Path) -> None:
    text = generate_excerpt(tmp_path)
    assert "O3'" not in text
    assert "constexpr std::array<LinkBond, 2> LINK_BONDS{{" in text
    assert "constexpr std::array<LinkAngle, 2> LINK_ANGLES{{" in text


def test_records_source_version_and_command(tmp_path: Path) -> None:
    text = generate_excerpt(tmp_path)
    assert 'return "Fri Feb 22 17:25:15 GMT 2013";' in text
    assert "standard_geometry_excerpt.cif" in text
    assert "Engh, R.A. & Huber, R. (2001)" in text


def test_unknown_link_id_is_an_error(tmp_path: Path) -> None:
    source = tmp_path / "bad.cif"
    source.write_text(
        EXCERPT.read_text(encoding="utf-8").replace("PRO-TRANS", "PRO-OTHER"),
        encoding="utf-8",
    )
    output = tmp_path / "out.cpp"
    result = run_generator(source, output)
    assert result.returncode == 1
    assert "PRO-OTHER" in result.stderr
    assert "C-N-CD" in result.stderr
    assert not output.exists()
```

- [ ] **Step 3: Run the generator tests to confirm they fail**

Run: `PYTHONPATH=python python -m pytest tests/python/test_generate_engh_huber.py -q`
Expected: FAIL; the generator exits with code 2 (`can't open file '.../scripts/generate_engh_huber.py'`), so the return-code assertions fail.

- [ ] **Step 4: Write the generator**

`scripts/generate_engh_huber.py`:

```python
"""Generate the Engh & Huber (2001) geometry tables from ``standard_geometry.cif``.

The input is the NDB/wwPDB file distributed with CCP4 as
``lib/data/standard_geometry.cif``. The output is a C++ source file holding
``constexpr`` tables of the heavy-atom bond and angle targets for the 20
standard amino acids, plus the peptide and disulfide link targets.

Usage::

    python3 scripts/generate_engh_huber.py standard_geometry.cif src/data/engh_huber_2001.cpp

The output is already clang-format clean; the tables are kept one row per
line inside ``clang-format off`` regions.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path

STANDARD_RESIDUES = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    }
)  # fmt: skip

# The source spells the default peptide link both ways, so both map to DEFAULT.
LINK_VARIANTS = {
    "NOT GLY PRO": "DEFAULT",
    "NOT PRO GLY": "DEFAULT",
    "GLY": "GLY",
    "PRO": "PRO",
    "PRO-CIS": "PRO_CIS",
    "PRO-TRANS": "PRO_TRANS",
    "PRO-ALL": "PRO_ALL",
    "CYS": "DISULFIDE",
}

HYDROGEN_SYMBOLS = frozenset({"H", "D"})

# A CIF token is a quoted string or a run of non-blank characters. Atom names
# such as O3' contain a quote character but are not quoted strings.
TOKEN_PATTERN = re.compile(r"\"[^\"]*\"|'[^']*'|\S+")


class GeneratorError(Exception):
    """Raised when the input cannot be converted faithfully."""


@dataclass(frozen=True)
class Loop:
    """One CIF loop: the field names of its category and its rows."""

    fields: list[str]
    rows: list[dict[str, str]]


@dataclass(frozen=True)
class CifData:
    """The parts of the CIF file the generator uses."""

    items: dict[str, str]
    loops: dict[str, Loop]
    malformed: set[str]


@dataclass(frozen=True)
class ResidueRow:
    """An intra-residue bond (two atoms) or angle (three atoms) target."""

    residue: str
    atoms: tuple[str, ...]
    variant: str
    ideal: str
    sigma: str


@dataclass(frozen=True)
class LinkRow:
    """An inter-residue bond or angle target with its normalized link variant."""

    atoms: tuple[str, ...]
    link_variant: str
    ideal: str
    sigma: str


def tokenize(line: str) -> list[str]:
    """Split one CIF line into tokens, dropping quotes and trailing comments."""
    tokens: list[str] = []
    for match in TOKEN_PATTERN.finditer(line):
        token = match.group(0)
        if token.startswith("#"):
            break
        if len(token) >= 2 and token[0] == token[-1] and token[0] in "\"'":
            token = token[1:-1]
        tokens.append(token)
    return tokens


def split_tag(tag: str) -> tuple[str, str]:
    """Split ``_category.field`` into ``("_category", "field")``."""
    category, _, field = tag.partition(".")
    if not field:
        raise GeneratorError(f"malformed tag {tag!r}")
    return category, field


def read_cif(text: str) -> CifData:
    """Parse the items and loops of a single-block CIF file."""
    items: dict[str, str] = {}
    loops: dict[str, Loop] = {}
    malformed: set[str] = set()
    category = ""
    fields: list[str] = []
    values: list[str] = []
    in_loop = False

    def close_loop() -> None:
        if not in_loop:
            return
        width = len(fields)
        # Some loops the generator never reads have more values per row than
        # fields; record them and fail only if one of them is needed.
        if len(values) % width != 0:
            malformed.add(category)
            return
        rows = [
            dict(zip(fields, values[start : start + width], strict=True))
            for start in range(0, len(values), width)
        ]
        loops[category] = Loop(fields=list(fields), rows=rows)

    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("data_"):
            continue
        if stripped == "loop_":
            close_loop()
            in_loop = True
            category = ""
            fields = []
            values = []
            continue
        if stripped.startswith("_"):
            tokens = tokenize(stripped)
            tag_category, field = split_tag(tokens[0])
            if in_loop and not values:
                if category and tag_category != category:
                    raise GeneratorError(f"loop mixes {category} and {tag_category}")
                category = tag_category
                fields.append(field)
                continue
            close_loop()
            in_loop = False
            if len(tokens) != 2:
                raise GeneratorError(f"item {tokens[0]} has no single value")
            items[tokens[0]] = tokens[1]
            continue
        if not in_loop:
            raise GeneratorError(f"unexpected data line: {stripped!r}")
        values.extend(tokenize(stripped))
    close_loop()
    return CifData(items=items, loops=loops, malformed=malformed)


def require_loop(cif: CifData, category: str) -> Loop:
    """Return a loop or fail with the missing category's name."""
    if category in cif.malformed:
        raise GeneratorError(f"loop {category} has a partial row")
    loop = cif.loops.get(category)
    if loop is None:
        raise GeneratorError(f"missing loop {category}")
    return loop


def has_peptide_pair(atoms: tuple[str, ...]) -> bool:
    """Return whether C and N are consecutive anywhere in the atom path."""
    return any({first, second} == {"C", "N"} for first, second in pairwise(atoms))


def hydrogen_atoms(cif: CifData) -> set[tuple[str, str]]:
    """Return ``(residue, atom)`` pairs whose element is hydrogen."""
    loop = require_loop(cif, "_chem_comp_atom")
    return {
        (row["comp_id"], row["atom_id"])
        for row in loop.rows
        if row["type_symbol"] in HYDROGEN_SYMBOLS
    }


def residue_rows(
    cif: CifData, category: str, atom_fields: tuple[str, ...], value_field: str
) -> list[ResidueRow]:
    """Extract heavy-atom rows for the standard residues from one loop."""
    hydrogens = hydrogen_atoms(cif)
    loop = require_loop(cif, category)
    rows: list[ResidueRow] = []
    for row in loop.rows:
        residue = row["comp_id"]
        if residue not in STANDARD_RESIDUES:
            continue
        atoms = tuple(row[field] for field in atom_fields)
        if any((residue, atom) in hydrogens for atom in atoms):
            continue
        # Rows with consecutive C and N describe the bond to the neighboring
        # residue (PRO C-N-CD); that term comes from the link table instead.
        if has_peptide_pair(atoms):
            continue
        rows.append(
            ResidueRow(
                residue=residue,
                atoms=atoms,
                variant=row["description"],
                ideal=row[value_field],
                sigma=row[value_field + "_esd"],
            )
        )
    return rows


def link_rows(
    cif: CifData, category: str, atom_fields: tuple[str, ...], value_field: str
) -> list[LinkRow]:
    """Extract peptide and disulfide link rows, skipping nucleic-acid links."""
    loop = require_loop(cif, category)
    rows: list[LinkRow] = []
    for row in loop.rows:
        atoms = tuple(row[field] for field in atom_fields)
        is_peptide = has_peptide_pair(atoms)
        is_disulfide = all(atom == "SG" for atom in atoms)
        if not is_peptide and not is_disulfide:
            continue
        link_id = row["link_id"]
        variant = LINK_VARIANTS.get(link_id)
        if variant is None or (variant == "DISULFIDE") != is_disulfide:
            raise GeneratorError(f"unknown link id {link_id!r} on {category} row {'-'.join(atoms)}")
        rows.append(
            LinkRow(
                atoms=atoms,
                link_variant=variant,
                ideal=row[value_field],
                sigma=row[value_field + "_esd"],
            )
        )
    return rows


def quote(text: str) -> str:
    """Return a C++ string literal for an atom, residue, or variant name."""
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def emit_table(type_name: str, array_name: str, entries: list[str]) -> list[str]:
    """Emit one ``constexpr std::array`` definition, one row per line."""
    # clang-format would pack several rows onto each line, which makes the
    # tables hard to compare with the source file.
    lines = [
        "// clang-format off",
        f"constexpr std::array<{type_name}, {len(entries)}> {array_name}{{{{",
    ]
    lines.extend(f"    {entry}," for entry in entries)
    lines += ["}};", "// clang-format on", ""]
    return lines


def accessor(row_type: str, function: str, array: str) -> str:
    """Emit the accessor that returns a view of one table."""
    return f"RowSpan<{row_type}> {function}() {{ return {{{array}.data(), {array}.size()}}; }}"


def render(
    version: str,
    source_name: str,
    bonds: list[ResidueRow],
    angles: list[ResidueRow],
    link_bonds: list[LinkRow],
    link_angles: list[LinkRow],
) -> str:
    """Render the generated C++ source file."""

    def residue_entry(row: ResidueRow) -> str:
        names = ", ".join(quote(name) for name in (row.residue, *row.atoms, row.variant))
        return f"{{{names}, {row.ideal}, {row.sigma}}}"

    def link_entry(row: LinkRow) -> str:
        names = ", ".join(quote(name) for name in row.atoms)
        return f"{{{names}, LinkVariant::{row.link_variant}, {row.ideal}, {row.sigma}}}"

    lines = [
        "// Engh & Huber (2001) ideal bond lengths and angles for the standard amino acids.",
        "//",
        "// Generated by scripts/generate_engh_huber.py. Do not edit by hand.",
        "//",
        f'// Source: {source_name}, _standard_geometry.version "{version}"',
        "// Reference: Engh, R.A. & Huber, R. (2001). International Tables for",
        "// Crystallography, Vol. F, pp. 382-392. Dordrecht: Kluwer Academic Publishers.",
        "// Command:",
        f"//   python3 scripts/generate_engh_huber.py {source_name} src/data/engh_huber_2001.cpp",
        "",
        '#include "data/engh_huber_2001.h"',
        "",
        "#include <array>",
        "",
        "namespace OEHEALTH::data {",
        "namespace {",
        "",
    ]
    lines += emit_table("ResidueBond", "RESIDUE_BONDS", [residue_entry(r) for r in bonds])
    lines += emit_table("ResidueAngle", "RESIDUE_ANGLES", [residue_entry(r) for r in angles])
    lines += emit_table("LinkBond", "LINK_BONDS", [link_entry(r) for r in link_bonds])
    lines += emit_table("LinkAngle", "LINK_ANGLES", [link_entry(r) for r in link_angles])
    lines += [
        "} // namespace",
        "",
        accessor("ResidueBond", "residue_bonds", "RESIDUE_BONDS"),
        "",
        accessor("ResidueAngle", "residue_angles", "RESIDUE_ANGLES"),
        "",
        accessor("LinkBond", "link_bonds", "LINK_BONDS"),
        "",
        accessor("LinkAngle", "link_angles", "LINK_ANGLES"),
        "",
        f"const char* source_version() {{ return {quote(version)}; }}",
        "",
        "} // namespace OEHEALTH::data",
        "",
    ]
    return "\n".join(lines)


def generate(source: Path) -> str:
    """Read ``standard_geometry.cif`` and return the generated C++ source."""
    cif = read_cif(source.read_text(encoding="utf-8"))
    version = cif.items.get("_standard_geometry.version")
    if version is None:
        raise GeneratorError("missing _standard_geometry.version")
    return render(
        version=version,
        source_name=source.name,
        bonds=residue_rows(cif, "_chem_comp_bond", ("atom_id_1", "atom_id_2"), "value_dist"),
        angles=residue_rows(
            cif, "_chem_comp_angle", ("atom_id_1", "atom_id_2", "atom_id_3"), "value_angle"
        ),
        link_bonds=link_rows(cif, "_chem_link_bond", ("atom_id_1", "atom_id_2"), "value_dist"),
        link_angles=link_rows(
            cif, "_chem_link_angle", ("atom_id_1", "atom_id_2", "atom_id_3"), "value_angle"
        ),
    )


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description="Generate the Engh & Huber (2001) tables from standard_geometry.cif."
    )
    parser.add_argument("source", type=Path, help="path to standard_geometry.cif")
    parser.add_argument("output", type=Path, help="path of the C++ file to write")
    args = parser.parse_args(argv)
    try:
        text = generate(args.source)
    except GeneratorError as error:
        print(f"generate_engh_huber: {error}", file=sys.stderr)
        return 1
    args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 5: Run the generator tests to confirm they pass**

Run: `PYTHONPATH=python python -m pytest tests/python/test_generate_engh_huber.py -q`
Expected: `7 passed`.

- [ ] **Step 6: Add the generator to mypy's file list**

```diff
diff --git a/pyproject.toml b/pyproject.toml
--- a/pyproject.toml
+++ b/pyproject.toml
@@ -123,7 +123,7 @@ known-first-party = ["oehealth"]
 [tool.mypy]
 python_version = "3.12"
 mypy_path = "python"
-files = ["python/oehealth", "tests/python"]
+files = ["python/oehealth", "scripts/generate_engh_huber.py", "tests/python"]
 exclude = ['^python/oehealth/oehealth\.py$', '^python/oehealth/_build_info\.py$']
 warn_unused_ignores = true
 warn_redundant_casts = true
```

- [ ] **Step 7: Write the failing reference data tests**

`tests/cpp/test_reference_data.cpp`:

```cpp
#include <gtest/gtest.h>

#include "data/engh_huber_2001.h"

#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace OEHEALTH::test {

namespace {

using data::LinkVariant;

bool equal(const char* a, const char* b) { return std::strcmp(a, b) == 0; }

const data::ResidueBond* residue_bond(const char* residue, const char* a1, const char* a2,
                                      const char* variant) {
    for (const auto& row : data::residue_bonds()) {
        if (equal(row.residue, residue) && equal(row.atom1, a1) && equal(row.atom2, a2) &&
            equal(row.variant, variant)) {
            return &row;
        }
    }
    return nullptr;
}

const data::ResidueAngle* residue_angle(const char* residue, const char* a1, const char* a2,
                                        const char* a3, const char* variant) {
    for (const auto& row : data::residue_angles()) {
        if (equal(row.residue, residue) && equal(row.atom1, a1) && equal(row.atom2, a2) &&
            equal(row.atom3, a3) && equal(row.variant, variant)) {
            return &row;
        }
    }
    return nullptr;
}

int count_link_bonds(const char* a1, const char* a2, LinkVariant variant) {
    int count = 0;
    for (const auto& row : data::link_bonds()) {
        count += static_cast<int>(equal(row.atom1, a1) && equal(row.atom2, a2) &&
                                  row.link_variant == variant);
    }
    return count;
}

int count_link_angles(const char* a1, const char* a2, const char* a3, LinkVariant variant) {
    int count = 0;
    for (const auto& row : data::link_angles()) {
        count += static_cast<int>(equal(row.atom1, a1) && equal(row.atom2, a2) &&
                                  equal(row.atom3, a3) && row.link_variant == variant);
    }
    return count;
}

const data::LinkAngle* link_angle(const char* a1, const char* a2, const char* a3,
                                  LinkVariant variant) {
    for (const auto& row : data::link_angles()) {
        if (equal(row.atom1, a1) && equal(row.atom2, a2) && equal(row.atom3, a3) &&
            row.link_variant == variant) {
            return &row;
        }
    }
    return nullptr;
}

void expect_bond(const char* residue, const char* a1, const char* a2, double ideal, double sigma) {
    const data::ResidueBond* row = residue_bond(residue, a1, a2, ".");
    ASSERT_NE(row, nullptr) << residue << " " << a1 << "-" << a2;
    EXPECT_DOUBLE_EQ(row->ideal, ideal);
    EXPECT_DOUBLE_EQ(row->sigma, sigma);
}

void expect_angle(const char* residue, const char* a1, const char* a2, const char* a3,
                  const char* variant, double ideal, double sigma) {
    const data::ResidueAngle* row = residue_angle(residue, a1, a2, a3, variant);
    ASSERT_NE(row, nullptr) << residue << " " << a1 << "-" << a2 << "-" << a3;
    EXPECT_DOUBLE_EQ(row->ideal, ideal);
    EXPECT_DOUBLE_EQ(row->sigma, sigma);
}

// PDB heavy-atom names of each standard residue, OXT included.
const std::map<std::string, std::vector<std::string>>& heavy_atoms() {
    static const std::map<std::string, std::vector<std::string>> atoms = {
        {"ALA", {"N", "CA", "C", "O", "CB", "OXT"}},
        {"ARG", {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2", "OXT"}},
        {"ASN", {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2", "OXT"}},
        {"ASP", {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2", "OXT"}},
        {"CYS", {"N", "CA", "C", "O", "CB", "SG", "OXT"}},
        {"GLN", {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2", "OXT"}},
        {"GLU", {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2", "OXT"}},
        {"GLY", {"N", "CA", "C", "O", "OXT"}},
        {"HIS", {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2", "OXT"}},
        {"ILE", {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1", "OXT"}},
        {"LEU", {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "OXT"}},
        {"LYS", {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ", "OXT"}},
        {"MET", {"N", "CA", "C", "O", "CB", "CG", "SD", "CE", "OXT"}},
        {"PHE", {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OXT"}},
        {"PRO", {"N", "CA", "C", "O", "CB", "CG", "CD", "OXT"}},
        {"SER", {"N", "CA", "C", "O", "CB", "OG", "OXT"}},
        {"THR", {"N", "CA", "C", "O", "CB", "OG1", "CG2", "OXT"}},
        {"TRP",
         {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2",
          "OXT"}},
        {"TYR", {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "OXT"}},
        {"VAL", {"N", "CA", "C", "O", "CB", "CG1", "CG2", "OXT"}},
    };
    return atoms;
}

} // namespace

TEST(ReferenceData, SpotValuesMatchEnghHuber2001) {
    expect_bond("ALA", "N", "CA", 1.459, 0.020);
    expect_bond("ALA", "CA", "C", 1.525, 0.026);
    expect_bond("ALA", "C", "O", 1.229, 0.019);
    expect_bond("ALA", "CA", "CB", 1.520, 0.021);
    expect_angle("ALA", "N", "CA", "C", ".", 111.0, 2.7);
    expect_angle("ALA", "CA", "C", "O", ".", 120.1, 2.1);
    expect_angle("GLY", "N", "CA", "C", ".", 113.1, 2.5);
    expect_angle("PRO", "CA", "N", "CD", "trans", 111.7, 1.4);
    expect_angle("PRO", "CA", "N", "CD", "cis", 111.5, 1.4);

    ASSERT_EQ(count_link_bonds("C", "N", LinkVariant::DEFAULT), 1);
    for (const auto& row : data::link_bonds()) {
        if (equal(row.atom1, "C") && row.link_variant == LinkVariant::DEFAULT) {
            EXPECT_DOUBLE_EQ(row.ideal, 1.336);
            EXPECT_DOUBLE_EQ(row.sigma, 0.023);
        }
    }
    const data::LinkAngle* ca_c_n = link_angle("CA", "C", "N", LinkVariant::DEFAULT);
    ASSERT_NE(ca_c_n, nullptr);
    EXPECT_DOUBLE_EQ(ca_c_n->ideal, 117.2);
    EXPECT_DOUBLE_EQ(ca_c_n->sigma, 2.2);
    const data::LinkAngle* o_c_n = link_angle("O", "C", "N", LinkVariant::DEFAULT);
    ASSERT_NE(o_c_n, nullptr);
    EXPECT_DOUBLE_EQ(o_c_n->ideal, 122.7);
    EXPECT_DOUBLE_EQ(o_c_n->sigma, 1.6);
    const data::LinkAngle* c_n_ca = link_angle("C", "N", "CA", LinkVariant::DEFAULT);
    ASSERT_NE(c_n_ca, nullptr);
    EXPECT_DOUBLE_EQ(c_n_ca->ideal, 121.7);
    EXPECT_DOUBLE_EQ(c_n_ca->sigma, 2.5);
}

TEST(ReferenceData, ResidueBondsConnectExactlyTheHeavyAtoms) {
    for (const auto& [residue, names] : heavy_atoms()) {
        const std::set<std::string> expected(names.begin(), names.end());
        std::map<std::string, std::set<std::string>> graph;
        for (const auto& row : data::residue_bonds()) {
            if (residue == row.residue) {
                graph[row.atom1].insert(row.atom2);
                graph[row.atom2].insert(row.atom1);
            }
        }
        std::set<std::string> nodes;
        for (const auto& entry : graph) {
            nodes.insert(entry.first);
        }
        EXPECT_EQ(nodes, expected) << residue;

        // Walk the graph from N; every heavy atom must be reachable.
        std::set<std::string> reached{"N"};
        std::vector<std::string> frontier{"N"};
        while (!frontier.empty()) {
            const std::string atom = frontier.back();
            frontier.pop_back();
            for (const std::string& next : graph[atom]) {
                if (reached.insert(next).second) {
                    frontier.push_back(next);
                }
            }
        }
        EXPECT_EQ(reached, expected) << residue;
    }
}

TEST(ReferenceData, ResidueAnglesUseOnlyHeavyAtoms) {
    for (const auto& row : data::residue_angles()) {
        const auto found = heavy_atoms().find(row.residue);
        ASSERT_NE(found, heavy_atoms().end()) << row.residue;
        const std::set<std::string> names(found->second.begin(), found->second.end());
        EXPECT_TRUE(names.count(row.atom1) == 1 && names.count(row.atom2) == 1 &&
                    names.count(row.atom3) == 1)
            << row.residue << " " << row.atom1 << "-" << row.atom2 << "-" << row.atom3;
    }
}

TEST(ReferenceData, EveryPeptideTermHasOneRowPerVariant) {
    for (LinkVariant variant : {LinkVariant::DEFAULT, LinkVariant::GLY, LinkVariant::PRO}) {
        EXPECT_EQ(count_link_bonds("C", "N", variant), 1);
        EXPECT_EQ(count_link_angles("CA", "C", "N", variant), 1);
        EXPECT_EQ(count_link_angles("O", "C", "N", variant), 1);
    }
    for (LinkVariant variant :
         {LinkVariant::DEFAULT, LinkVariant::GLY, LinkVariant::PRO_CIS, LinkVariant::PRO_TRANS}) {
        EXPECT_EQ(count_link_angles("C", "N", "CA", variant), 1);
    }
    for (LinkVariant variant : {LinkVariant::PRO_CIS, LinkVariant::PRO_TRANS}) {
        EXPECT_EQ(count_link_angles("C", "N", "CD", variant), 1);
    }
    EXPECT_EQ(count_link_bonds("SG", "SG", LinkVariant::DISULFIDE), 1);
}

TEST(ReferenceData, ProlineHasNoIntraResidueCNCDAngle) {
    for (const auto& row : data::residue_angles()) {
        EXPECT_FALSE(equal(row.residue, "PRO") && equal(row.atom1, "C") && equal(row.atom2, "N") &&
                     equal(row.atom3, "CD"));
    }
}

TEST(ReferenceData, RecordsSourceVersion) {
    EXPECT_STREQ(data::source_version(), "Fri Feb 22 17:25:15 GMT 2013");
}

} // namespace OEHEALTH::test
```

Add to `add_executable(oehealth_tests ...)` in `tests/cpp/CMakeLists.txt`, keeping the list sorted:

```cmake
    test_reference_data.cpp
```

Run: `cmake --build build-debug -j8`
Expected: FAIL to compile with `'data/engh_huber_2001.h' file not found`.

- [ ] **Step 8: Write the table declarations**

```bash
mkdir -p src/data
```

`src/data/engh_huber_2001.h`:

```cpp
#ifndef OEHEALTH_DATA_ENGH_HUBER_2001_H
#define OEHEALTH_DATA_ENGH_HUBER_2001_H

#include <cstddef>

namespace OEHEALTH::data {

/// \brief Which peptide or disulfide link row applies to a pair of residues.
///
/// The source file names these by link id; the generator normalizes both
/// spellings of the default row ("NOT GLY PRO" and "NOT PRO GLY") to DEFAULT.
enum class LinkVariant { DEFAULT, GLY, PRO, PRO_CIS, PRO_TRANS, PRO_ALL, DISULFIDE };

/// \brief An intra-residue bond target. Lengths are in angstroms.
///
/// ``variant`` is "." when the row applies to every form of the residue, and
/// otherwise one of "SS", "SH" (CYS), "cis", "trans" (PRO), or "HISD",
/// "HISE", "HISH" (HIS).
struct ResidueBond {
    const char* residue;
    const char* atom1;
    const char* atom2;
    const char* variant;
    double ideal;
    double sigma;
};

/// \brief An intra-residue angle target. Angles are in degrees; atom2 is the vertex.
struct ResidueAngle {
    const char* residue;
    const char* atom1;
    const char* atom2;
    const char* atom3;
    const char* variant;
    double ideal;
    double sigma;
};

/// \brief A link bond target: atom1 is in residue i, atom2 in residue j.
struct LinkBond {
    const char* atom1;
    const char* atom2;
    LinkVariant link_variant;
    double ideal;
    double sigma;
};

/// \brief A link angle target, in degrees, with atom2 as the vertex.
///
/// For peptide links, atoms named before the C-N pair belong to residue i
/// and atoms after it to residue j.
struct LinkAngle {
    const char* atom1;
    const char* atom2;
    const char* atom3;
    LinkVariant link_variant;
    double ideal;
    double sigma;
};

/// \brief A read-only view of one generated table.
template <typename Row>
struct RowSpan {
    const Row* first;
    std::size_t count;

    const Row* begin() const { return first; }
    const Row* end() const { return first + count; }
    std::size_t size() const { return count; }
};

RowSpan<ResidueBond> residue_bonds();
RowSpan<ResidueAngle> residue_angles();
RowSpan<LinkBond> link_bonds();
RowSpan<LinkAngle> link_angles();

/// \brief The ``_standard_geometry.version`` string of the source file.
const char* source_version();

} // namespace OEHEALTH::data

#endif // OEHEALTH_DATA_ENGH_HUBER_2001_H
```

- [ ] **Step 9: Generate the tables**

```bash
python3 scripts/generate_engh_huber.py $CCP4/lib/data/standard_geometry.cif src/data/engh_huber_2001.cpp
grep -n 'std::array<' src/data/engh_huber_2001.cpp
shasum -a 256 src/data/engh_huber_2001.cpp
```

Expected: the four arrays are `std::array<ResidueBond, 191>`, `std::array<ResidueAngle, 199>`, `std::array<LinkBond, 4>` and `std::array<LinkAngle, 13>`, and the checksum is `6375405d2e34907fc9dad0d11885f4bb87493430c2e2062bef4476457a2c936b`. If either differs, stop and report it.

- [ ] **Step 10: Register the source**

Add to `set(OEHEALTH_SOURCES ...)` in `CMakeLists.txt`, keeping the list sorted:

```cmake
    src/data/engh_huber_2001.cpp
```

- [ ] **Step 11: Build and run the tests**

```bash
cmake --build build-debug -j8
./build-debug/tests/cpp/oehealth_tests --gtest_filter='ReferenceData.*'
ctest --test-dir build-debug --output-on-failure
PYTHONPATH=python python -m pytest tests/python -q
```

Expected: `[  PASSED  ] 6 tests.` from the filtered run; ctest passes every test; pytest reports `12 passed` (7 generator tests and the template's 5).

- [ ] **Step 12: Format and lint**

```bash
clang-format --dry-run --Werror src/data/engh_huber_2001.h src/data/engh_huber_2001.cpp tests/cpp/test_reference_data.cpp
clang-tidy -p build-debug --quiet --extra-arg=-isysroot --extra-arg="$(xcrun --show-sdk-path)" src/data/engh_huber_2001.cpp tests/cpp/test_reference_data.cpp
ruff check .
ruff format --check scripts/generate_engh_huber.py tests/python/test_generate_engh_huber.py
mypy --python-executable "$(command -v python)"
```

Expected: no output from clang-format (the generator writes its tables between `// clang-format off` and `// clang-format on`); no clang-tidy findings; `All checks passed!`; `2 files already formatted`; `Success: no issues found in 5 source files`.

- [ ] **Step 13: Commit**

```bash
git add scripts/generate_engh_huber.py tests/data/standard_geometry_excerpt.cif tests/python/test_generate_engh_huber.py src/data/engh_huber_2001.h src/data/engh_huber_2001.cpp tests/cpp/test_reference_data.cpp pyproject.toml CMakeLists.txt tests/cpp/CMakeLists.txt
git commit -m "Generate Engh & Huber (2001) reference tables from standard_geometry.cif"
```

---
### Task 5: Geometry targets and variant selection

**Risk tier:** standard — the variant rules (CYS, PRO, HIS, peptide links, disulfides) that decide every bond and angle target.

**Depends on:** Task 2, Task 4

**Files:**
- Create: `src/geometry/targets.h`, `src/geometry/targets.cpp`
- Test: `tests/cpp/test_targets.cpp`
- Modify: `CMakeLists.txt` (`OEHEALTH_SOURCES`), `tests/cpp/CMakeLists.txt`

**Interfaces:**
- Consumes: `StructureContext`, `ResidueInfo` (Task 2); `geometry::distance`, `angle_degrees`, `torsion_degrees`, `peptide_successor`, `peptide_predecessor` (Task 2, `src/geometry/measure.h`); `data::residue_bonds()`, `residue_angles()`, `link_bonds()`, `link_angles()`, `data::LinkVariant` (Task 4).
- Produces (private, `src/geometry/targets.h`, namespace `OEHEALTH::geometry`):
  - `struct GeometryTerm { std::vector<const OEChem::OEAtomBase*> atoms; double ideal; double sigma; std::string variant; }` (two atoms for a bond, three for an angle with the vertex second).
  - `std::vector<GeometryTerm> collect_bond_terms(const StructureContext&)` and `std::vector<GeometryTerm> collect_angle_terms(const StructureContext&)`.
  - `std::string residue_variant(const StructureContext&, std::size_t residue_index)`.
  - `bool preceding_peptide_is_cis(const StructureContext&, std::size_t residue_index)`.
  - `std::string_view his_variant_from_hydrogens(bool has_hd1, bool has_he2)`.
  - `std::string lowest_rmsz_his_variant(double rmsz_hise, double rmsz_hisd, double rmsz_hish)`.
  - `double his_variant_rmsz(const StructureContext&, std::size_t residue_index, const std::string& variant)`.
  - `std::string_view link_variant_name(data::LinkVariant)`.

The GLY C–N assumption in the design notes applies to this task: `add_link_bond` keys the GLY row on residue j.

- [ ] **Step 1: Write the failing tests**

`tests/cpp/test_targets.cpp`:

```cpp
#include <gtest/gtest.h>

#include "geometry/measure.h"
#include "geometry/targets.h"
#include "test_helpers.h"

#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace OEHEALTH::test {

namespace {

// Index of a residue number in a crambin context (residues are numbered 1-46 in order).
std::size_t index_of(int residue_number) { return static_cast<std::size_t>(residue_number - 1); }

bool term_has_atoms(const geometry::GeometryTerm& term,
                    const std::vector<const OEChem::OEAtomBase*>& atoms) {
    if (term.atoms.size() != atoms.size()) {
        return false;
    }
    for (std::size_t k = 0; k < atoms.size(); ++k) {
        if (atoms[k] == nullptr || term.atoms[k]->GetIdx() != atoms[k]->GetIdx()) {
            return false;
        }
    }
    return true;
}

const geometry::GeometryTerm* find_term(const std::vector<geometry::GeometryTerm>& terms,
                                        const std::vector<const OEChem::OEAtomBase*>& atoms) {
    for (const auto& term : terms) {
        if (term_has_atoms(term, atoms)) {
            return &term;
        }
    }
    return nullptr;
}

std::size_t count_c_n_cd_terms(const std::vector<geometry::GeometryTerm>& terms) {
    std::size_t count = 0;
    for (const auto& term : terms) {
        if (term.atoms.size() == 3 && trimmed(term.atoms[0]->GetName()) == "C" &&
            trimmed(term.atoms[1]->GetName()) == "N" && trimmed(term.atoms[2]->GetName()) == "CD") {
            ++count;
        }
    }
    return count;
}

// Heavy atoms of HIS from the ideal coordinates in standard_geometry.cif.
void add_ideal_histidine(OEChem::OEMolBase& mol) {
    struct IdealAtom {
        const char* name;
        unsigned int element;
        double x, y, z;
    };
    const std::vector<IdealAtom> atoms = {
        {"N", OEChem::OEElemNo::N, -0.040, -1.210, 0.053},
        {"CA", OEChem::OEElemNo::C, 1.172, -1.709, 0.652},
        {"C", OEChem::OEElemNo::C, 1.083, -3.207, 0.905},
        {"O", OEChem::OEElemNo::O, 0.040, -3.770, 1.222},
        {"CB", OEChem::OEElemNo::C, 1.484, -0.975, 1.962},
        {"CG", OEChem::OEElemNo::C, 2.940, -1.060, 2.353},
        {"ND1", OEChem::OEElemNo::N, 3.380, -2.075, 3.129},
        {"CD2", OEChem::OEElemNo::C, 3.960, -0.251, 2.046},
        {"CE1", OEChem::OEElemNo::C, 4.693, -1.908, 3.317},
        {"NE2", OEChem::OEElemNo::N, 5.058, -0.801, 2.662},
    };
    for (const IdealAtom& atom : atoms) {
        add_atom(mol, atom.element, atom.name, "HIS", 1, atom.x, atom.y, atom.z);
    }
}

} // namespace

TEST(Targets, CrambinCysteinesAreDisulfideBonded) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const StructureContext ctx(mol);
    for (int number : {3, 4, 16, 26, 32, 40}) {
        EXPECT_EQ(geometry::residue_variant(ctx, index_of(number)), "SS") << number;
    }
}

TEST(Targets, CysteineWithoutPartnerIsSH) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    mol.DeleteAtom(find_atom(mol, 40, "SG"));
    const StructureContext ctx(mol);
    EXPECT_EQ(geometry::residue_variant(ctx, index_of(3)), "SH");
    EXPECT_EQ(geometry::residue_variant(ctx, index_of(40)), "SH");
    EXPECT_EQ(geometry::residue_variant(ctx, index_of(4)), "SS");
}

TEST(Targets, CrambinProlinesAreTrans) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const StructureContext ctx(mol);
    for (int number : {5, 19, 22, 36, 41}) {
        EXPECT_FALSE(geometry::preceding_peptide_is_cis(ctx, index_of(number))) << number;
        EXPECT_EQ(geometry::residue_variant(ctx, index_of(number)), "trans") << number;
    }
}

TEST(Targets, ProlineIsCisWhenOmegaIsNearZero) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    // The C40-N41 bond is outside every disulfide ring, so it can be rotated.
    ASSERT_TRUE(OEChem::OESetTorsion(mol, find_atom(mol, 40, "CA"), find_atom(mol, 40, "C"),
                                     find_atom(mol, 41, "N"), find_atom(mol, 41, "CA"), 0.0));
    const StructureContext ctx(mol);
    EXPECT_TRUE(geometry::preceding_peptide_is_cis(ctx, index_of(41)));
    EXPECT_EQ(geometry::residue_variant(ctx, index_of(41)), "cis");

    const auto angles = geometry::collect_angle_terms(ctx);
    const geometry::GeometryTerm* c_n_cd = find_term(
        angles, {find_atom(mol, 40, "C"), find_atom(mol, 41, "N"), find_atom(mol, 41, "CD")});
    ASSERT_NE(c_n_cd, nullptr);
    EXPECT_EQ(c_n_cd->variant, "PRO_CIS");
    EXPECT_DOUBLE_EQ(c_n_cd->ideal, 120.6);
    const geometry::GeometryTerm* c_n_ca = find_term(
        angles, {find_atom(mol, 40, "C"), find_atom(mol, 41, "N"), find_atom(mol, 41, "CA")});
    ASSERT_NE(c_n_ca, nullptr);
    EXPECT_EQ(c_n_ca->variant, "PRO_CIS");
}

TEST(Targets, HistidineVariantFromHydrogens) {
    EXPECT_EQ(geometry::his_variant_from_hydrogens(true, true), "HISH");
    EXPECT_EQ(geometry::his_variant_from_hydrogens(true, false), "HISD");
    EXPECT_EQ(geometry::his_variant_from_hydrogens(false, true), "HISE");
    EXPECT_EQ(geometry::his_variant_from_hydrogens(false, false), "");
}

TEST(Targets, HistidineLowestRmszWithTieOrder) {
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_EQ(geometry::lowest_rmsz_his_variant(1.0, 2.0, 3.0), "HISE");
    EXPECT_EQ(geometry::lowest_rmsz_his_variant(2.0, 1.0, 3.0), "HISD");
    EXPECT_EQ(geometry::lowest_rmsz_his_variant(3.0, 2.0, 1.0), "HISH");
    EXPECT_EQ(geometry::lowest_rmsz_his_variant(1.0, 1.0, 1.0), "HISE");
    EXPECT_EQ(geometry::lowest_rmsz_his_variant(2.0, 1.0, 1.0), "HISD");
    EXPECT_EQ(geometry::lowest_rmsz_his_variant(inf, inf, inf), "HISE");
}

TEST(Targets, HistidineUsesExplicitHydrogens) {
    for (const auto& [hydrogens, expected] :
         std::vector<std::pair<std::vector<std::string>, std::string>>{
             {{"HD1"}, "HISD"}, {{"HE2"}, "HISE"}, {{"HD1", "HE2"}, "HISH"}}) {
        OEChem::OEGraphMol mol;
        add_ideal_histidine(mol);
        for (const std::string& name : hydrogens) {
            add_atom(mol, OEChem::OEElemNo::H, name, "HIS", 1, 9.0, 9.0, 9.0);
        }
        const StructureContext ctx(mol);
        EXPECT_EQ(geometry::residue_variant(ctx, 0), expected);
    }
}

TEST(Targets, HistidineWithoutHydrogensUsesBestFit) {
    OEChem::OEGraphMol mol;
    add_ideal_histidine(mol);
    const StructureContext ctx(mol);
    const std::string expected = geometry::lowest_rmsz_his_variant(
        geometry::his_variant_rmsz(ctx, 0, "HISE"), geometry::his_variant_rmsz(ctx, 0, "HISD"),
        geometry::his_variant_rmsz(ctx, 0, "HISH"));
    EXPECT_EQ(geometry::residue_variant(ctx, 0), expected);
    EXPECT_LT(geometry::his_variant_rmsz(ctx, 0, "HISE"), std::numeric_limits<double>::infinity());

    // The chosen variant's rows, and no other variant's, become terms.
    for (const auto& term : geometry::collect_bond_terms(ctx)) {
        EXPECT_TRUE(term.variant == "." || term.variant == expected) << term.variant;
    }
}

TEST(Targets, CrambinLinkTerms) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const StructureContext ctx(mol);
    const auto bonds = geometry::collect_bond_terms(ctx);
    const auto angles = geometry::collect_angle_terms(ctx);

    std::size_t peptide_bonds = 0;
    std::size_t disulfides = 0;
    for (const auto& term : bonds) {
        if (term.variant == "DISULFIDE") {
            ++disulfides;
        } else if (term.variant == "DEFAULT" || term.variant == "GLY" || term.variant == "PRO") {
            ++peptide_bonds;
        }
    }
    EXPECT_EQ(peptide_bonds, 45U);
    EXPECT_EQ(disulfides, 3U);
    EXPECT_EQ(count_c_n_cd_terms(angles), 5U);

    // C(4)-N(5): residue j is PRO.
    const auto* cys_pro = find_term(bonds, {find_atom(mol, 4, "C"), find_atom(mol, 5, "N")});
    ASSERT_NE(cys_pro, nullptr);
    EXPECT_EQ(cys_pro->variant, "PRO");
    EXPECT_DOUBLE_EQ(cys_pro->ideal, 1.338);
    // C(19)-N(20): residue j is GLY, so the C-N bond uses the GLY row.
    const auto* pro_gly = find_term(bonds, {find_atom(mol, 19, "C"), find_atom(mol, 20, "N")});
    ASSERT_NE(pro_gly, nullptr);
    EXPECT_EQ(pro_gly->variant, "GLY");
    // CA(20)-C(20)-N(21): residue i is GLY.
    const auto* gly_thr = find_term(
        angles, {find_atom(mol, 20, "CA"), find_atom(mol, 20, "C"), find_atom(mol, 21, "N")});
    ASSERT_NE(gly_thr, nullptr);
    EXPECT_EQ(gly_thr->variant, "GLY");
    EXPECT_DOUBLE_EQ(gly_thr->ideal, 116.2);
    // C(20)-N(21)-CA(21): residue j is THR, so DEFAULT.
    const auto* default_angle = find_term(
        angles, {find_atom(mol, 20, "C"), find_atom(mol, 21, "N"), find_atom(mol, 21, "CA")});
    ASSERT_NE(default_angle, nullptr);
    EXPECT_EQ(default_angle->variant, "DEFAULT");
}

TEST(Targets, ProlineWithoutPredecessorHasNoCNCDTerm) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    std::vector<OEChem::OEAtomBase*> residue4;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (OEChem::OEAtomGetResidue(atom).GetResidueNumber() == 4) {
            residue4.push_back(atom);
        }
    }
    for (OEChem::OEAtomBase* atom : residue4) {
        mol.DeleteAtom(atom);
    }
    const StructureContext ctx(mol);
    EXPECT_EQ(count_c_n_cd_terms(geometry::collect_angle_terms(ctx)), 4U);
}

} // namespace OEHEALTH::test
```

- [ ] **Step 2: Register the test and confirm it fails**

Add to `add_executable(oehealth_tests ...)` in `tests/cpp/CMakeLists.txt`, keeping the list sorted:

```cmake
    test_targets.cpp
```

Run: `cmake --build build-debug -j8`
Expected: FAIL to compile with `'geometry/targets.h' file not found`.

- [ ] **Step 3: Write the target selection**

`src/geometry/targets.h`:

```cpp
#ifndef OEHEALTH_SRC_GEOMETRY_TARGETS_H
#define OEHEALTH_SRC_GEOMETRY_TARGETS_H

#include "data/engh_huber_2001.h"
#include "oehealth/context.h"

#include <oechem.h>

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace OEHEALTH::geometry {

/// \brief One bond or angle to assess, with its Engh & Huber (2001) target.
struct GeometryTerm {
    /// Two atoms for a bond; three for an angle, with the vertex second.
    std::vector<const OEChem::OEAtomBase*> atoms;
    double ideal = 0.0; ///< Angstroms or degrees.
    double sigma = 0.0;
    /// The residue variant of the row (".", "SS", "SH", "cis", "trans",
    /// "HISD", "HISE", "HISH") or, for link terms, the link variant name
    /// ("DEFAULT", "GLY", "PRO", "PRO_CIS", "PRO_TRANS", "DISULFIDE").
    std::string variant;
};

/// \brief Every heavy-atom bond term in the structure: intra-residue,
/// peptide link, and disulfide. Terms with a missing atom are left out.
std::vector<GeometryTerm> collect_bond_terms(const StructureContext& ctx);

/// \brief Every heavy-atom angle term in the structure: intra-residue and
/// peptide link. Terms with a missing atom are left out.
std::vector<GeometryTerm> collect_angle_terms(const StructureContext& ctx);

/// \brief The variant used for a residue's intra-residue rows: "SS" or "SH"
/// for CYS, "cis" or "trans" for PRO, "HISD", "HISE", or "HISH" for HIS,
/// and "." for every other residue.
std::string residue_variant(const StructureContext& ctx, std::size_t residue_index);

/// \brief Whether the peptide bond preceding residue j is cis (|omega| < 90).
///
/// False when residue j has no peptide predecessor or a CA atom is missing.
bool preceding_peptide_is_cis(const StructureContext& ctx, std::size_t residue_index);

/// \brief The HIS variant implied by explicit hydrogens: "HISH" for both HD1
/// and HE2, "HISD" for HD1 only, "HISE" for HE2 only, "" for neither.
std::string_view his_variant_from_hydrogens(bool has_hd1, bool has_he2);

/// \brief The HIS variant with the lowest RMSZ; ties go to HISE, then HISD,
/// then HISH.
std::string lowest_rmsz_his_variant(double rmsz_hise, double rmsz_hisd, double rmsz_hish);

/// \brief RMSZ of one HIS residue against one variant's bond and angle rows.
///
/// Only rows with that exact variant are used, and only those whose atoms
/// are present. Returns infinity when no such row can be measured.
double his_variant_rmsz(const StructureContext& ctx, std::size_t residue_index,
                        const std::string& variant);

/// \brief The name of a link variant, e.g. "PRO_CIS".
std::string_view link_variant_name(data::LinkVariant variant);

} // namespace OEHEALTH::geometry

#endif // OEHEALTH_SRC_GEOMETRY_TARGETS_H
```

`src/geometry/targets.cpp`:

```cpp
#include "geometry/targets.h"

#include "geometry/measure.h"

#include <cmath>
#include <cstring>
#include <limits>
#include <optional>
#include <utility>

namespace OEHEALTH::geometry {

namespace {

using data::LinkVariant;

// SG atoms of two CYS residues this close, in angstroms, form a disulfide.
constexpr double DISULFIDE_CUTOFF = 2.5;

bool is_residue(const ResidueInfo& residue, const char* name) { return residue.name == name; }

bool row_applies(const char* row_variant, const std::string& variant) {
    return std::strcmp(row_variant, ".") == 0 || variant == row_variant;
}

bool has_disulfide(const StructureContext& ctx, std::size_t residue_index) {
    const auto& residues = ctx.GetResidues();
    const OEChem::OEAtomBase* sulfur = residues[residue_index].FindAtom("SG");
    if (sulfur == nullptr) {
        return false;
    }
    for (std::size_t other = 0; other < residues.size(); ++other) {
        if (other == residue_index || !is_residue(residues[other], "CYS")) {
            continue;
        }
        const OEChem::OEAtomBase* partner = residues[other].FindAtom("SG");
        if (partner != nullptr &&
            distance(ctx.GetMolecule(), sulfur, partner) <= DISULFIDE_CUTOFF) {
            return true;
        }
    }
    return false;
}

std::string his_variant(const StructureContext& ctx, std::size_t residue_index) {
    const ResidueInfo& residue = ctx.GetResidues()[residue_index];
    const std::string_view from_hydrogens = his_variant_from_hydrogens(
        residue.FindAtom("HD1") != nullptr, residue.FindAtom("HE2") != nullptr);
    if (!from_hydrogens.empty()) {
        return std::string(from_hydrogens);
    }
    return lowest_rmsz_his_variant(his_variant_rmsz(ctx, residue_index, "HISE"),
                                   his_variant_rmsz(ctx, residue_index, "HISD"),
                                   his_variant_rmsz(ctx, residue_index, "HISH"));
}

const data::LinkBond* find_link_bond(const char* atom1, const char* atom2, LinkVariant variant) {
    for (const data::LinkBond& row : data::link_bonds()) {
        if (std::strcmp(row.atom1, atom1) == 0 && std::strcmp(row.atom2, atom2) == 0 &&
            row.link_variant == variant) {
            return &row;
        }
    }
    return nullptr;
}

const data::LinkAngle* find_link_angle(const char* atom1, const char* atom2, const char* atom3,
                                       LinkVariant variant) {
    for (const data::LinkAngle& row : data::link_angles()) {
        if (std::strcmp(row.atom1, atom1) == 0 && std::strcmp(row.atom2, atom2) == 0 &&
            std::strcmp(row.atom3, atom3) == 0 && row.link_variant == variant) {
            return &row;
        }
    }
    return nullptr;
}

LinkVariant pro_link_variant(const StructureContext& ctx, std::size_t j) {
    return preceding_peptide_is_cis(ctx, j) ? LinkVariant::PRO_CIS : LinkVariant::PRO_TRANS;
}

void add_term(std::vector<GeometryTerm>& terms, std::vector<const OEChem::OEAtomBase*> atoms,
              double ideal, double sigma, std::string_view variant) {
    for (const OEChem::OEAtomBase* atom : atoms) {
        if (atom == nullptr) {
            return;
        }
    }
    terms.push_back(GeometryTerm{std::move(atoms), ideal, sigma, std::string(variant)});
}

void add_residue_bonds(const StructureContext& ctx, std::size_t i,
                       std::vector<GeometryTerm>& terms) {
    const ResidueInfo& residue = ctx.GetResidues()[i];
    const std::string variant = residue_variant(ctx, i);
    for (const data::ResidueBond& row : data::residue_bonds()) {
        if (residue.name == row.residue && row_applies(row.variant, variant)) {
            add_term(terms, {residue.FindAtom(row.atom1), residue.FindAtom(row.atom2)}, row.ideal,
                     row.sigma, row.variant);
        }
    }
}

void add_residue_angles(const StructureContext& ctx, std::size_t i,
                        std::vector<GeometryTerm>& terms) {
    const ResidueInfo& residue = ctx.GetResidues()[i];
    const std::string variant = residue_variant(ctx, i);
    for (const data::ResidueAngle& row : data::residue_angles()) {
        if (residue.name == row.residue && row_applies(row.variant, variant)) {
            add_term(terms,
                     {residue.FindAtom(row.atom1), residue.FindAtom(row.atom2),
                      residue.FindAtom(row.atom3)},
                     row.ideal, row.sigma, row.variant);
        }
    }
}

void add_link_bond(const StructureContext& ctx, std::size_t i, std::size_t j,
                   std::vector<GeometryTerm>& terms) {
    const ResidueInfo& first = ctx.GetResidues()[i];
    const ResidueInfo& second = ctx.GetResidues()[j];
    // Assumption from the design: the GLY C-N row is keyed on residue j.
    LinkVariant variant = LinkVariant::DEFAULT;
    if (is_residue(second, "PRO")) {
        variant = LinkVariant::PRO;
    } else if (is_residue(second, "GLY")) {
        variant = LinkVariant::GLY;
    }
    if (const data::LinkBond* row = find_link_bond("C", "N", variant)) {
        add_term(terms, {first.FindAtom("C"), second.FindAtom("N")}, row->ideal, row->sigma,
                 link_variant_name(variant));
    }
}

void add_link_angles(const StructureContext& ctx, std::size_t i, std::size_t j,
                     std::vector<GeometryTerm>& terms) {
    const ResidueInfo& first = ctx.GetResidues()[i];
    const ResidueInfo& second = ctx.GetResidues()[j];
    const bool second_is_pro = is_residue(second, "PRO");

    // CA(i)-C(i)-N(j) and O(i)-C(i)-N(j): the varying atom belongs to residue i.
    LinkVariant before = LinkVariant::DEFAULT;
    if (second_is_pro) {
        before = LinkVariant::PRO;
    } else if (is_residue(first, "GLY")) {
        before = LinkVariant::GLY;
    }
    for (const char* atom : {"CA", "O"}) {
        if (const data::LinkAngle* row = find_link_angle(atom, "C", "N", before)) {
            add_term(terms, {first.FindAtom(atom), first.FindAtom("C"), second.FindAtom("N")},
                     row->ideal, row->sigma, link_variant_name(before));
        }
    }

    // C(i)-N(j)-CA(j), and C(i)-N(j)-CD(j) for PRO: the varying atom belongs to residue j.
    LinkVariant after = LinkVariant::DEFAULT;
    if (second_is_pro) {
        after = pro_link_variant(ctx, j);
    } else if (is_residue(second, "GLY")) {
        after = LinkVariant::GLY;
    }
    if (const data::LinkAngle* row = find_link_angle("C", "N", "CA", after)) {
        add_term(terms, {first.FindAtom("C"), second.FindAtom("N"), second.FindAtom("CA")},
                 row->ideal, row->sigma, link_variant_name(after));
    }
    if (second_is_pro) {
        if (const data::LinkAngle* row = find_link_angle("C", "N", "CD", after)) {
            add_term(terms, {first.FindAtom("C"), second.FindAtom("N"), second.FindAtom("CD")},
                     row->ideal, row->sigma, link_variant_name(after));
        }
    }
}

void add_disulfide_bonds(const StructureContext& ctx, std::vector<GeometryTerm>& terms) {
    const auto& residues = ctx.GetResidues();
    const data::LinkBond* row = find_link_bond("SG", "SG", LinkVariant::DISULFIDE);
    if (row == nullptr) {
        return;
    }
    for (std::size_t a = 0; a < residues.size(); ++a) {
        if (!is_residue(residues[a], "CYS")) {
            continue;
        }
        const OEChem::OEAtomBase* first = residues[a].FindAtom("SG");
        for (std::size_t b = a + 1; b < residues.size() && first != nullptr; ++b) {
            if (!is_residue(residues[b], "CYS")) {
                continue;
            }
            const OEChem::OEAtomBase* second = residues[b].FindAtom("SG");
            if (second != nullptr &&
                distance(ctx.GetMolecule(), first, second) <= DISULFIDE_CUTOFF) {
                add_term(terms, {first, second}, row->ideal, row->sigma,
                         link_variant_name(LinkVariant::DISULFIDE));
            }
        }
    }
}

} // namespace

std::string residue_variant(const StructureContext& ctx, std::size_t residue_index) {
    const ResidueInfo& residue = ctx.GetResidues()[residue_index];
    if (is_residue(residue, "CYS")) {
        return has_disulfide(ctx, residue_index) ? "SS" : "SH";
    }
    if (is_residue(residue, "PRO")) {
        return preceding_peptide_is_cis(ctx, residue_index) ? "cis" : "trans";
    }
    if (is_residue(residue, "HIS")) {
        return his_variant(ctx, residue_index);
    }
    return ".";
}

bool preceding_peptide_is_cis(const StructureContext& ctx, std::size_t residue_index) {
    const auto predecessor = peptide_predecessor(ctx, residue_index);
    if (!predecessor) {
        return false;
    }
    const ResidueInfo& first = ctx.GetResidues()[*predecessor];
    const ResidueInfo& second = ctx.GetResidues()[residue_index];
    const OEChem::OEAtomBase* ca_i = first.FindAtom("CA");
    const OEChem::OEAtomBase* c_i = first.FindAtom("C");
    const OEChem::OEAtomBase* n_j = second.FindAtom("N");
    const OEChem::OEAtomBase* ca_j = second.FindAtom("CA");
    if (ca_i == nullptr || ca_j == nullptr) {
        return false;
    }
    const double omega = torsion_degrees(ctx.GetMolecule(), ca_i, c_i, n_j, ca_j);
    return std::fabs(omega) < 90.0;
}

std::string_view his_variant_from_hydrogens(bool has_hd1, bool has_he2) {
    if (has_hd1 && has_he2) {
        return "HISH";
    }
    if (has_hd1) {
        return "HISD";
    }
    if (has_he2) {
        return "HISE";
    }
    return "";
}

std::string lowest_rmsz_his_variant(double rmsz_hise, double rmsz_hisd, double rmsz_hish) {
    std::string best = "HISE";
    double best_rmsz = rmsz_hise;
    if (rmsz_hisd < best_rmsz) {
        best = "HISD";
        best_rmsz = rmsz_hisd;
    }
    if (rmsz_hish < best_rmsz) {
        best = "HISH";
    }
    return best;
}

double his_variant_rmsz(const StructureContext& ctx, std::size_t residue_index,
                        const std::string& variant) {
    const OEChem::OEMolBase& mol = ctx.GetMolecule();
    const ResidueInfo& residue = ctx.GetResidues()[residue_index];
    double sum_z2 = 0.0;
    std::size_t count = 0;
    for (const data::ResidueBond& row : data::residue_bonds()) {
        if (residue.name != row.residue || variant != row.variant) {
            continue;
        }
        const OEChem::OEAtomBase* a = residue.FindAtom(row.atom1);
        const OEChem::OEAtomBase* b = residue.FindAtom(row.atom2);
        if (a != nullptr && b != nullptr) {
            const double z = (distance(mol, a, b) - row.ideal) / row.sigma;
            sum_z2 += z * z;
            ++count;
        }
    }
    for (const data::ResidueAngle& row : data::residue_angles()) {
        if (residue.name != row.residue || variant != row.variant) {
            continue;
        }
        const OEChem::OEAtomBase* a = residue.FindAtom(row.atom1);
        const OEChem::OEAtomBase* b = residue.FindAtom(row.atom2);
        const OEChem::OEAtomBase* c = residue.FindAtom(row.atom3);
        if (a != nullptr && b != nullptr && c != nullptr) {
            const double z = (angle_degrees(mol, a, b, c) - row.ideal) / row.sigma;
            sum_z2 += z * z;
            ++count;
        }
    }
    if (count == 0) {
        return std::numeric_limits<double>::infinity();
    }
    return std::sqrt(sum_z2 / static_cast<double>(count));
}

std::string_view link_variant_name(data::LinkVariant variant) {
    switch (variant) {
        case LinkVariant::DEFAULT:
            return "DEFAULT";
        case LinkVariant::GLY:
            return "GLY";
        case LinkVariant::PRO:
            return "PRO";
        case LinkVariant::PRO_CIS:
            return "PRO_CIS";
        case LinkVariant::PRO_TRANS:
            return "PRO_TRANS";
        case LinkVariant::PRO_ALL:
            return "PRO_ALL";
        case LinkVariant::DISULFIDE:
            return "DISULFIDE";
    }
    return "UNKNOWN";
}

std::vector<GeometryTerm> collect_bond_terms(const StructureContext& ctx) {
    std::vector<GeometryTerm> terms;
    const auto& residues = ctx.GetResidues();
    for (std::size_t i = 0; i < residues.size(); ++i) {
        if (!residues[i].is_standard_amino_acid) {
            continue;
        }
        add_residue_bonds(ctx, i, terms);
        if (const auto j = peptide_successor(ctx, i)) {
            add_link_bond(ctx, i, *j, terms);
        }
    }
    add_disulfide_bonds(ctx, terms);
    return terms;
}

std::vector<GeometryTerm> collect_angle_terms(const StructureContext& ctx) {
    std::vector<GeometryTerm> terms;
    const auto& residues = ctx.GetResidues();
    for (std::size_t i = 0; i < residues.size(); ++i) {
        if (!residues[i].is_standard_amino_acid) {
            continue;
        }
        add_residue_angles(ctx, i, terms);
        if (const auto j = peptide_successor(ctx, i)) {
            add_link_angles(ctx, i, *j, terms);
        }
    }
    return terms;
}

} // namespace OEHEALTH::geometry
```

- [ ] **Step 4: Register the source**

Add to `set(OEHEALTH_SOURCES ...)` in `CMakeLists.txt`, keeping the list sorted:

```cmake
    src/geometry/targets.cpp
```

- [ ] **Step 5: Build and run the tests**

```bash
cmake --build build-debug -j8
./build-debug/tests/cpp/oehealth_tests --gtest_filter='Targets.*'
ctest --test-dir build-debug --output-on-failure
```

Expected: `[  PASSED  ] 10 tests.` from the filtered run; ctest passes every test.

- [ ] **Step 6: Format and lint**

```bash
clang-format --dry-run --Werror src/geometry/targets.h src/geometry/targets.cpp tests/cpp/test_targets.cpp
clang-tidy -p build-debug --quiet --extra-arg=-isysroot --extra-arg="$(xcrun --show-sdk-path)" src/geometry/targets.cpp tests/cpp/test_targets.cpp
```

Expected: no output from clang-format; no findings from clang-tidy.

- [ ] **Step 7: Commit**

```bash
git add src/geometry/targets.h src/geometry/targets.cpp tests/cpp/test_targets.cpp CMakeLists.txt tests/cpp/CMakeLists.txt
git commit -m "Select Engh & Huber bond and angle targets per residue and peptide link"
```

---

### Task 6: Bond length and bond angle assessments

**Risk tier:** standard — two public assessments and their shared private logic, six new files.

**Depends on:** Task 5

**Files:**
- Create: `include/oehealth/assessments/bond_lengths.h`, `include/oehealth/assessments/bond_angles.h`
- Create: `src/assessments/geometry_outliers.h`, `src/assessments/geometry_outliers.cpp`
- Create: `src/assessments/bond_lengths.cpp`, `src/assessments/bond_angles.cpp`
- Test: `tests/cpp/test_assessments_geometry.cpp`
- Modify: `CMakeLists.txt` (`OEHEALTH_SOURCES`), `tests/cpp/CMakeLists.txt`

**Interfaces:**
- Consumes: `Assessment`, `StructureContext`, result types (Task 2); `make_atom_ref`, `atom_path`, `format_fixed`, `status_from_findings`, `meets_threshold`, `validate_thresholds` (Task 2, `src/common.h`); `geometry::distance`, `geometry::angle_degrees` (Task 2); `geometry::GeometryTerm`, `collect_bond_terms`, `collect_angle_terms` (Task 5).
- Produces (public):
  - `struct BondLengthOptions { double warning_z = 4.0; double error_z = 10.0; }` and `class BondLengthAssessment : public Assessment { explicit BondLengthAssessment(const BondLengthOptions& = BondLengthOptions()); const BondLengthOptions& GetOptions() const; }`, name `"bond_lengths"`, finding code `"bond_length_outlier"`.
  - `struct BondAngleOptions { double warning_z = 4.0; double error_z = 10.0; }` and `class BondAngleAssessment : public Assessment { explicit BondAngleAssessment(const BondAngleOptions& = BondAngleOptions()); const BondAngleOptions& GetOptions() const; }`, name `"bond_angles"`, finding code `"bond_angle_outlier"`.
  - Metrics, in order: `n_assessed`, `n_outliers`, `rmsz`, `max_abs_z`; on NOT_APPLICABLE only `n_assessed` and `n_outliers`, both 0.
- Produces (private, `src/assessments/geometry_outliers.h`): `struct GeometryOutlierSettings { std::string assessment, code, unit; int decimals; double warning_z, error_z; }` and `AssessmentResult assess_geometry_terms(const OEChem::OEMolBase&, const std::vector<geometry::GeometryTerm>&, const GeometryOutlierSettings&)`.

- [ ] **Step 1: Write the failing tests**

`tests/cpp/test_assessments_geometry.cpp`:

```cpp
#include <gtest/gtest.h>

#include "oehealth/assessments/bond_angles.h"
#include "oehealth/assessments/bond_lengths.h"
#include "test_helpers.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace OEHEALTH::test {

namespace {

// Ideal SER CA-CB-OG angle and its sigma in Engh & Huber (2001), degrees.
constexpr double SER_CA_CB_OG = 111.2;
constexpr double SER_CA_CB_OG_SIGMA = 2.7;
// Ideal ALA CA-CB bond length and its sigma, angstroms.
constexpr double ALA_CA_CB = 1.520;
constexpr double ALA_CA_CB_SIGMA = 0.021;
constexpr double PI = 3.14159265358979323846;

std::vector<const Finding*> findings_mentioning(const AssessmentResult& result,
                                                const OEChem::OEAtomBase* atom) {
    std::vector<const Finding*> found;
    for (const Finding& finding : result.findings) {
        if (mentions_atom(finding, atom)) {
            found.push_back(&finding);
        }
    }
    return found;
}

bool has_error_finding(const AssessmentResult& result) {
    for (const Finding& finding : result.findings) {
        if (finding.severity == Severity::ERROR) {
            return true;
        }
    }
    return false;
}

// Move `moving` along the anchor->moving direction to sit `length` from the anchor.
void set_bond_length(OEChem::OEMolBase& mol, const OEChem::OEAtomBase* anchor,
                     OEChem::OEAtomBase* moving, double length) {
    const auto a = coords(mol, anchor);
    const auto m = coords(mol, moving);
    const std::array<double, 3> v{m[0] - a[0], m[1] - a[1], m[2] - a[2]};
    const double norm = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    set_coords(
        mol, moving,
        {a[0] + v[0] / norm * length, a[1] + v[1] / norm * length, a[2] + v[2] / norm * length});
}

// Rotate `moving` about `vertex`, in the plane of the three atoms, so the
// end-vertex-moving angle becomes `degrees` while |vertex-moving| is kept.
void set_bond_angle(OEChem::OEMolBase& mol, const OEChem::OEAtomBase* end,
                    const OEChem::OEAtomBase* vertex, OEChem::OEAtomBase* moving, double degrees) {
    const auto e = coords(mol, end);
    const auto v = coords(mol, vertex);
    const auto m = coords(mol, moving);
    std::array<double, 3> u{e[0] - v[0], e[1] - v[1], e[2] - v[2]};
    const double u_norm = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    for (double& c : u) {
        c /= u_norm;
    }
    const std::array<double, 3> w{m[0] - v[0], m[1] - v[1], m[2] - v[2]};
    const double length = std::sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
    const double along = w[0] * u[0] + w[1] * u[1] + w[2] * u[2];
    std::array<double, 3> p{w[0] - along * u[0], w[1] - along * u[1], w[2] - along * u[2]};
    const double p_norm = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    for (double& c : p) {
        c /= p_norm;
    }
    const double theta = degrees * PI / 180.0;
    set_coords(mol, moving,
               {v[0] + length * (std::cos(theta) * u[0] + std::sin(theta) * p[0]),
                v[1] + length * (std::cos(theta) * u[1] + std::sin(theta) * p[1]),
                v[2] + length * (std::cos(theta) * u[2] + std::sin(theta) * p[2])});
}

// An ALA fragment with only CA and CB, stored in double precision, so the
// CA-CB distance is exactly the requested length.
void make_ala_fragment(OEChem::OEMol& mol, double length) {
    add_atom(mol, OEChem::OEElemNo::C, "CA", "ALA", 1, 0.0, 0.0, 0.0);
    add_atom(mol, OEChem::OEElemNo::C, "CB", "ALA", 1, length, 0.0, 0.0);
}

} // namespace

TEST(BondLengthAssessment, DefaultsAndName) {
    const BondLengthAssessment assessment;
    EXPECT_EQ(assessment.GetName(), "bond_lengths");
    EXPECT_FALSE(assessment.GetDescription().empty());
    EXPECT_DOUBLE_EQ(assessment.GetOptions().warning_z, 4.0);
    EXPECT_DOUBLE_EQ(assessment.GetOptions().error_z, 10.0);
}

TEST(BondAngleAssessment, DefaultsAndName) {
    const BondAngleAssessment assessment;
    EXPECT_EQ(assessment.GetName(), "bond_angles");
    EXPECT_FALSE(assessment.GetDescription().empty());
    EXPECT_DOUBLE_EQ(assessment.GetOptions().warning_z, 4.0);
    EXPECT_DOUBLE_EQ(assessment.GetOptions().error_z, 10.0);
}

TEST(GeometryOptions, RejectInvalidThresholds) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW(BondLengthAssessment(BondLengthOptions{nan, 10.0}), std::invalid_argument);
    EXPECT_THROW(BondLengthAssessment(BondLengthOptions{4.0, inf}), std::invalid_argument);
    EXPECT_THROW(BondLengthAssessment(BondLengthOptions{-1.0, 10.0}), std::invalid_argument);
    EXPECT_THROW(BondLengthAssessment(BondLengthOptions{5.0, 4.0}), std::invalid_argument);
    EXPECT_NO_THROW(BondLengthAssessment(BondLengthOptions{4.0, 4.0}));
    EXPECT_THROW(BondAngleAssessment(BondAngleOptions{nan, 10.0}), std::invalid_argument);
    EXPECT_THROW(BondAngleAssessment(BondAngleOptions{4.0, inf}), std::invalid_argument);
    EXPECT_THROW(BondAngleAssessment(BondAngleOptions{-1.0, 10.0}), std::invalid_argument);
    EXPECT_THROW(BondAngleAssessment(BondAngleOptions{5.0, 4.0}), std::invalid_argument);
    EXPECT_NO_THROW(BondAngleAssessment(BondAngleOptions{4.0, 4.0}));
}

TEST(BondLengthAssessment, CrambinHasNoErrors) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const AssessmentResult result = BondLengthAssessment().Run(mol);
    EXPECT_EQ(result.assessment, "bond_lengths");
    EXPECT_FALSE(has_error_finding(result));
    ASSERT_EQ(result.metrics.size(), 4U);
    EXPECT_EQ(result.metrics[0].name, "n_assessed");
    EXPECT_EQ(result.metrics[1].name, "n_outliers");
    EXPECT_EQ(result.metrics[2].name, "rmsz");
    EXPECT_EQ(result.metrics[3].name, "max_abs_z");
    EXPECT_GT(metric(result, "n_assessed"), 300.0);
    EXPECT_EQ(metric(result, "n_outliers"), static_cast<double>(result.findings.size()));
}

TEST(BondAngleAssessment, CrambinHasNoErrors) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const AssessmentResult result = BondAngleAssessment().Run(mol);
    EXPECT_EQ(result.assessment, "bond_angles");
    EXPECT_FALSE(has_error_finding(result));
    EXPECT_GT(metric(result, "n_assessed"), 400.0);
}

TEST(BondLengthAssessment, StretchedBondIsTheOnlyNewOutlier) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const OEChem::OEAtomBase* ca = find_atom(mol, 27, "CA");
    OEChem::OEAtomBase* cb = find_atom(mol, 27, "CB");
    const AssessmentResult before = BondLengthAssessment().Run(mol);
    ASSERT_TRUE(findings_mentioning(before, cb).empty()) << "pick another ALA";

    // CB is bonded only to CA among heavy atoms, so this changes one bond.
    set_bond_length(mol, ca, cb, ALA_CA_CB + 6.0 * ALA_CA_CB_SIGMA);
    const AssessmentResult after = BondLengthAssessment().Run(mol);
    const auto found = findings_mentioning(after, cb);
    ASSERT_EQ(found.size(), 1U);
    const Finding& finding = *found[0];
    EXPECT_EQ(finding.code, "bond_length_outlier");
    EXPECT_EQ(finding.severity, Severity::WARNING);
    ASSERT_EQ(finding.atoms.size(), 2U);
    EXPECT_EQ(finding.atoms[0].atom_name, "CA");
    EXPECT_EQ(finding.atoms[1].atom_name, "CB");
    EXPECT_NEAR(finding.deviation, 6.0, 1e-3);
    EXPECT_NEAR(finding.expected, ALA_CA_CB, 1e-12);
    EXPECT_EQ(finding.message.rfind("A ALA 27 CA-CB", 0), 0U) << finding.message;
    EXPECT_EQ(metric(after, "n_outliers"), metric(before, "n_outliers") + 1.0);
    EXPECT_EQ(after.status, Status::WARN);
}

TEST(BondLengthAssessment, BondExactlyAtWarningThresholdIsAWarning) {
    OEChem::OEMol mol(OEChem::OEMCMolType::DoubleCartesian);
    make_ala_fragment(mol, ALA_CA_CB + 4.0 * ALA_CA_CB_SIGMA);
    const AssessmentResult result = BondLengthAssessment().Run(mol);
    EXPECT_EQ(metric(result, "n_assessed"), 1.0);
    ASSERT_EQ(result.findings.size(), 1U);
    EXPECT_EQ(result.findings[0].severity, Severity::WARNING);
    EXPECT_NEAR(result.findings[0].deviation, 4.0, 1e-9);
}

TEST(BondLengthAssessment, ValueAtSharedThresholdIsAnError) {
    OEChem::OEMol mol(OEChem::OEMCMolType::DoubleCartesian);
    make_ala_fragment(mol, ALA_CA_CB + 4.0 * ALA_CA_CB_SIGMA);
    const AssessmentResult result = BondLengthAssessment(BondLengthOptions{4.0, 4.0}).Run(mol);
    ASSERT_EQ(result.findings.size(), 1U);
    EXPECT_EQ(result.findings[0].severity, Severity::ERROR);
    EXPECT_EQ(result.status, Status::FAIL);
}

TEST(BondAngleAssessment, DistortedAngleIsTheOnlyNewOutlier) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const OEChem::OEAtomBase* ca = find_atom(mol, 11, "CA");
    const OEChem::OEAtomBase* cb = find_atom(mol, 11, "CB");
    OEChem::OEAtomBase* og = find_atom(mol, 11, "OG");
    const AssessmentResult before = BondAngleAssessment().Run(mol);
    ASSERT_TRUE(findings_mentioning(before, og).empty()) << "pick another SER";

    // OG is bonded only to CB, so the only angle that changes is CA-CB-OG.
    set_bond_angle(mol, ca, cb, og, SER_CA_CB_OG + 6.0 * SER_CA_CB_OG_SIGMA);
    const AssessmentResult after = BondAngleAssessment().Run(mol);
    const auto found = findings_mentioning(after, og);
    ASSERT_EQ(found.size(), 1U);
    const Finding& finding = *found[0];
    EXPECT_EQ(finding.code, "bond_angle_outlier");
    EXPECT_EQ(finding.severity, Severity::WARNING);
    ASSERT_EQ(finding.atoms.size(), 3U);
    EXPECT_EQ(finding.atoms[0].atom_name, "CA");
    EXPECT_EQ(finding.atoms[1].atom_name, "CB");
    EXPECT_EQ(finding.atoms[2].atom_name, "OG");
    EXPECT_NEAR(finding.deviation, 6.0, 1e-3);
    EXPECT_EQ(metric(after, "n_outliers"), metric(before, "n_outliers") + 1.0);
}

TEST(GeometryAssessments, EmptyAndNonProteinMoleculesAreNotApplicable) {
    const OEChem::OEGraphMol empty;
    EXPECT_EQ(BondLengthAssessment().Run(empty).status, Status::NOT_APPLICABLE);
    EXPECT_EQ(BondAngleAssessment().Run(empty).status, Status::NOT_APPLICABLE);

    OEChem::OEGraphMol ethanol;
    OEChem::OESmilesToMol(ethanol, "CCO");
    ethanol.SetDimension(3);
    const AssessmentResult result = BondLengthAssessment().Run(ethanol);
    EXPECT_EQ(result.status, Status::NOT_APPLICABLE);
    EXPECT_FALSE(result.message.empty());
    EXPECT_TRUE(result.findings.empty());
    EXPECT_EQ(metric(result, "n_assessed"), 0.0);
}

} // namespace OEHEALTH::test
```

- [ ] **Step 2: Register the test and confirm it fails**

Add to `add_executable(oehealth_tests ...)` in `tests/cpp/CMakeLists.txt`, keeping the list sorted:

```cmake
    test_assessments_geometry.cpp
```

Run: `cmake --build build-debug -j8`
Expected: FAIL to compile with `'oehealth/assessments/bond_angles.h' file not found`.

- [ ] **Step 3: Write the public headers**

Create the directories if an earlier task has not:

```bash
mkdir -p include/oehealth/assessments src/assessments
```

`include/oehealth/assessments/bond_lengths.h`:

```cpp
#ifndef OEHEALTH_ASSESSMENTS_BOND_LENGTHS_H
#define OEHEALTH_ASSESSMENTS_BOND_LENGTHS_H

#include "oehealth/assessment.h"

namespace OEHEALTH {

/// \brief Thresholds on |Z| for bond length outliers.
struct BondLengthOptions {
    double warning_z = 4.0;
    double error_z = 10.0;
};

/// \brief Compares heavy-atom bond lengths of standard amino acids with
/// Engh & Huber (2001) targets.
///
/// Z = (observed - ideal) / sigma for every intra-residue, peptide link,
/// and disulfide bond. |Z| >= warning_z gives a WARNING finding and
/// |Z| >= error_z an ERROR finding, code "bond_length_outlier".
class BondLengthAssessment : public Assessment {
public:
    /// \throws std::invalid_argument for NaN, infinite, or negative
    ///         thresholds, or error_z < warning_z.
    explicit BondLengthAssessment(const BondLengthOptions& options = BondLengthOptions());

    const BondLengthOptions& GetOptions() const;
    std::string GetName() const override;
    std::string GetDescription() const override;
    AssessmentResult Assess(const StructureContext& ctx) const override;

private:
    BondLengthOptions options_;
};

} // namespace OEHEALTH

#endif // OEHEALTH_ASSESSMENTS_BOND_LENGTHS_H
```

`include/oehealth/assessments/bond_angles.h`:

```cpp
#ifndef OEHEALTH_ASSESSMENTS_BOND_ANGLES_H
#define OEHEALTH_ASSESSMENTS_BOND_ANGLES_H

#include "oehealth/assessment.h"

namespace OEHEALTH {

/// \brief Thresholds on |Z| for bond angle outliers.
struct BondAngleOptions {
    double warning_z = 4.0;
    double error_z = 10.0;
};

/// \brief Compares heavy-atom bond angles of standard amino acids with
/// Engh & Huber (2001) targets.
///
/// Z = (observed - ideal) / sigma, in degrees, for every intra-residue and
/// peptide link angle. |Z| >= warning_z gives a WARNING finding and
/// |Z| >= error_z an ERROR finding, code "bond_angle_outlier".
class BondAngleAssessment : public Assessment {
public:
    /// \throws std::invalid_argument for NaN, infinite, or negative
    ///         thresholds, or error_z < warning_z.
    explicit BondAngleAssessment(const BondAngleOptions& options = BondAngleOptions());

    const BondAngleOptions& GetOptions() const;
    std::string GetName() const override;
    std::string GetDescription() const override;
    AssessmentResult Assess(const StructureContext& ctx) const override;

private:
    BondAngleOptions options_;
};

} // namespace OEHEALTH

#endif // OEHEALTH_ASSESSMENTS_BOND_ANGLES_H
```

- [ ] **Step 4: Write the shared outlier logic**

`src/assessments/geometry_outliers.h`:

```cpp
#ifndef OEHEALTH_SRC_ASSESSMENTS_GEOMETRY_OUTLIERS_H
#define OEHEALTH_SRC_ASSESSMENTS_GEOMETRY_OUTLIERS_H

#include "geometry/targets.h"
#include "oehealth/result.h"

#include <oechem.h>

#include <string>
#include <vector>

namespace OEHEALTH {

/// \brief What differs between the bond length and bond angle assessments.
struct GeometryOutlierSettings {
    std::string assessment; ///< Name of the producing assessment.
    std::string code;       ///< Finding code.
    std::string unit;       ///< Unit used in messages: "A" or "deg".
    int decimals = 3;       ///< Decimals for observed and ideal values in messages.
    double warning_z = 4.0;
    double error_z = 10.0;
};

/// \brief Measure each term, compute Z, and build the result shared by
/// bond_lengths and bond_angles.
///
/// Two-atom terms are measured as distances and three-atom terms as angles.
AssessmentResult assess_geometry_terms(const OEChem::OEMolBase& mol,
                                       const std::vector<geometry::GeometryTerm>& terms,
                                       const GeometryOutlierSettings& settings);

} // namespace OEHEALTH

#endif // OEHEALTH_SRC_ASSESSMENTS_GEOMETRY_OUTLIERS_H
```

`src/assessments/geometry_outliers.cpp`:

```cpp
#include "assessments/geometry_outliers.h"

#include "common.h"
#include "geometry/measure.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace OEHEALTH {

AssessmentResult assess_geometry_terms(const OEChem::OEMolBase& mol,
                                       const std::vector<geometry::GeometryTerm>& terms,
                                       const GeometryOutlierSettings& settings) {
    AssessmentResult result;
    result.assessment = settings.assessment;
    if (terms.empty()) {
        result.status = Status::NOT_APPLICABLE;
        result.message = "no standard amino-acid terms to assess";
        result.metrics = {{"n_assessed", 0.0}, {"n_outliers", 0.0}};
        return result;
    }

    double sum_z2 = 0.0;
    double max_abs_z = 0.0;
    for (const geometry::GeometryTerm& term : terms) {
        const auto& atoms = term.atoms;
        const double observed = atoms.size() == 2
                                    ? geometry::distance(mol, atoms[0], atoms[1])
                                    : geometry::angle_degrees(mol, atoms[0], atoms[1], atoms[2]);
        const double z = (observed - term.ideal) / term.sigma;
        const double abs_z = std::fabs(z);
        sum_z2 += z * z;
        max_abs_z = std::max(max_abs_z, abs_z);

        Severity severity = Severity::INFO;
        if (meets_threshold(abs_z, settings.error_z)) {
            severity = Severity::ERROR;
        } else if (meets_threshold(abs_z, settings.warning_z)) {
            severity = Severity::WARNING;
        } else {
            continue;
        }

        Finding finding;
        finding.code = settings.code;
        finding.severity = severity;
        for (const OEChem::OEAtomBase* atom : atoms) {
            finding.atoms.push_back(make_atom_ref(atom));
        }
        finding.observed = observed;
        finding.expected = term.ideal;
        finding.deviation = z;
        finding.message = atom_path(finding.atoms) + " " +
                          format_fixed(observed, settings.decimals) + " " + settings.unit +
                          " (ideal " + format_fixed(term.ideal, settings.decimals) + ", Z " +
                          format_fixed(z, 1) + ")";
        result.findings.push_back(std::move(finding));
    }

    const auto n_assessed = static_cast<double>(terms.size());
    result.metrics = {
        {"n_assessed", n_assessed},
        {"n_outliers", static_cast<double>(result.findings.size())},
        {"rmsz", std::sqrt(sum_z2 / n_assessed)},
        {"max_abs_z", max_abs_z},
    };
    result.status = status_from_findings(result.findings);
    return result;
}

} // namespace OEHEALTH
```

- [ ] **Step 5: Write the two assessments**

`src/assessments/bond_lengths.cpp`:

```cpp
#include "oehealth/assessments/bond_lengths.h"

#include "assessments/geometry_outliers.h"
#include "common.h"
#include "geometry/targets.h"

namespace OEHEALTH {

BondLengthAssessment::BondLengthAssessment(const BondLengthOptions& options) : options_(options) {
    validate_thresholds("bond length", options_.warning_z, options_.error_z);
}

const BondLengthOptions& BondLengthAssessment::GetOptions() const { return options_; }

std::string BondLengthAssessment::GetName() const { return "bond_lengths"; }

std::string BondLengthAssessment::GetDescription() const {
    return "Compares heavy-atom bond lengths of standard amino acids with Engh & Huber (2001) "
           "targets.";
}

AssessmentResult BondLengthAssessment::Assess(const StructureContext& ctx) const {
    GeometryOutlierSettings settings;
    settings.assessment = GetName();
    settings.code = "bond_length_outlier";
    settings.unit = "A";
    settings.decimals = 3;
    settings.warning_z = options_.warning_z;
    settings.error_z = options_.error_z;
    return assess_geometry_terms(ctx.GetMolecule(), geometry::collect_bond_terms(ctx), settings);
}

} // namespace OEHEALTH
```

`src/assessments/bond_angles.cpp`:

```cpp
#include "oehealth/assessments/bond_angles.h"

#include "assessments/geometry_outliers.h"
#include "common.h"
#include "geometry/targets.h"

namespace OEHEALTH {

BondAngleAssessment::BondAngleAssessment(const BondAngleOptions& options) : options_(options) {
    validate_thresholds("bond angle", options_.warning_z, options_.error_z);
}

const BondAngleOptions& BondAngleAssessment::GetOptions() const { return options_; }

std::string BondAngleAssessment::GetName() const { return "bond_angles"; }

std::string BondAngleAssessment::GetDescription() const {
    return "Compares heavy-atom bond angles of standard amino acids with Engh & Huber (2001) "
           "targets.";
}

AssessmentResult BondAngleAssessment::Assess(const StructureContext& ctx) const {
    GeometryOutlierSettings settings;
    settings.assessment = GetName();
    settings.code = "bond_angle_outlier";
    settings.unit = "deg";
    settings.decimals = 1;
    settings.warning_z = options_.warning_z;
    settings.error_z = options_.error_z;
    return assess_geometry_terms(ctx.GetMolecule(), geometry::collect_angle_terms(ctx), settings);
}

} // namespace OEHEALTH
```

- [ ] **Step 6: Register the sources**

Add to `set(OEHEALTH_SOURCES ...)` in `CMakeLists.txt`, keeping the list sorted:

```cmake
    src/assessments/bond_angles.cpp
    src/assessments/bond_lengths.cpp
    src/assessments/geometry_outliers.cpp
```

- [ ] **Step 7: Build and run the tests**

```bash
cmake --build build-debug -j8
./build-debug/tests/cpp/oehealth_tests --gtest_filter='BondLengthAssessment.*:BondAngleAssessment.*:GeometryOptions.*:GeometryAssessments.*'
ctest --test-dir build-debug --output-on-failure
```

Expected: `[  PASSED  ] 10 tests.` from the filtered run; ctest passes every test.

- [ ] **Step 8: Format and lint**

```bash
clang-format --dry-run --Werror include/oehealth/assessments/bond_lengths.h include/oehealth/assessments/bond_angles.h src/assessments/geometry_outliers.h src/assessments/geometry_outliers.cpp src/assessments/bond_lengths.cpp src/assessments/bond_angles.cpp tests/cpp/test_assessments_geometry.cpp
clang-tidy -p build-debug --quiet --extra-arg=-isysroot --extra-arg="$(xcrun --show-sdk-path)" src/assessments/geometry_outliers.cpp src/assessments/bond_lengths.cpp src/assessments/bond_angles.cpp tests/cpp/test_assessments_geometry.cpp
```

Expected: no output from clang-format; no findings from clang-tidy.

- [ ] **Step 9: Commit**

```bash
git add include/oehealth/assessments/bond_lengths.h include/oehealth/assessments/bond_angles.h src/assessments/geometry_outliers.h src/assessments/geometry_outliers.cpp src/assessments/bond_lengths.cpp src/assessments/bond_angles.cpp tests/cpp/test_assessments_geometry.cpp CMakeLists.txt tests/cpp/CMakeLists.txt
git commit -m "Add bond length and bond angle assessments"
```

---

### Task 7: Clash assessment

**Risk tier:** standard — a public assessment with a neighbor search, pair exclusions, and threshold boundaries, three new files.

**Depends on:** Task 2

**Files:**
- Create: `include/oehealth/assessments/clashes.h`, `src/assessments/clashes.cpp`
- Test: `tests/cpp/test_clashes.cpp`
- Modify: `CMakeLists.txt` (`OEHEALTH_SOURCES`), `tests/cpp/CMakeLists.txt`

**Interfaces:**
- Consumes: `Assessment`, `StructureContext`, result types (Task 2); `make_atom_ref`, `atom_path`, `format_fixed`, `status_from_findings`, `meets_threshold`, `validate_thresholds` (Task 2, `src/common.h`); `geometry::distance` (Task 2); `OEChem::OEGetNearestNbrs`, `OEGetBondiVdWRadius`.
- Produces: `struct ClashOptions { double warning_overlap = 0.4; double error_overlap = 1.0; double polar_contact_distance = 2.9; }` and `class ClashAssessment : public Assessment { explicit ClashAssessment(const ClashOptions& = ClashOptions()); const ClashOptions& GetOptions() const; }`, name `"clashes"`, finding code `"clash"`. Metrics, in order: `n_clashes`, `clashscore`, `max_overlap`, `hydrogens_present`; on NOT_APPLICABLE only `hydrogens_present`.

- [ ] **Step 1: Write the failing tests**

`tests/cpp/test_clashes.cpp`:

```cpp
#include <gtest/gtest.h>

#include "oehealth/assessments/clashes.h"
#include "test_helpers.h"

#include <limits>
#include <stdexcept>

namespace OEHEALTH::test {

namespace {

// Two unbonded atoms on the x axis; coordinates like these are exact in floats.
AssessmentResult assess_pair(unsigned int first, unsigned int second, double separation) {
    OEChem::OEGraphMol mol;
    add_atom(mol, first, "A1", "LIG", 1, 0.0, 0.0, 0.0);
    add_atom(mol, second, "A2", "LIG", 2, separation, 0.0, 0.0);
    return ClashAssessment().Run(mol);
}

} // namespace

TEST(ClashAssessment, DefaultsAndName) {
    const ClashAssessment assessment;
    EXPECT_EQ(assessment.GetName(), "clashes");
    EXPECT_FALSE(assessment.GetDescription().empty());
    EXPECT_DOUBLE_EQ(assessment.GetOptions().warning_overlap, 0.4);
    EXPECT_DOUBLE_EQ(assessment.GetOptions().error_overlap, 1.0);
    EXPECT_DOUBLE_EQ(assessment.GetOptions().polar_contact_distance, 2.9);
}

TEST(ClashAssessment, RejectsInvalidOptions) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW(ClashAssessment(ClashOptions{nan, 1.0, 2.9}), std::invalid_argument);
    EXPECT_THROW(ClashAssessment(ClashOptions{0.4, inf, 2.9}), std::invalid_argument);
    EXPECT_THROW(ClashAssessment(ClashOptions{-0.1, 1.0, 2.9}), std::invalid_argument);
    EXPECT_THROW(ClashAssessment(ClashOptions{1.0, 0.4, 2.9}), std::invalid_argument);
    EXPECT_THROW(ClashAssessment(ClashOptions{0.4, 1.0, 0.0}), std::invalid_argument);
    EXPECT_THROW(ClashAssessment(ClashOptions{0.4, 1.0, -2.9}), std::invalid_argument);
    EXPECT_THROW(ClashAssessment(ClashOptions{0.4, 1.0, nan}), std::invalid_argument);
    EXPECT_NO_THROW(ClashAssessment(ClashOptions{0.4, 0.4, 2.9}));
}

TEST(ClashAssessment, UnbondedCarbonsTwoAngstromsApartAreAnError) {
    const AssessmentResult result = assess_pair(OEChem::OEElemNo::C, OEChem::OEElemNo::C, 2.0);
    ASSERT_EQ(result.findings.size(), 1U);
    const Finding& clash = result.findings[0];
    EXPECT_EQ(clash.code, "clash");
    EXPECT_EQ(clash.severity, Severity::ERROR);
    EXPECT_EQ(clash.atoms.size(), 2U);
    EXPECT_NEAR(clash.observed, 2.0, 1e-9);
    EXPECT_NEAR(clash.expected, 3.4, 1e-9);
    EXPECT_NEAR(clash.deviation, 1.4, 1e-9);
    EXPECT_EQ(result.status, Status::FAIL);
    ASSERT_EQ(result.metrics.size(), 4U);
    EXPECT_EQ(result.metrics[0].name, "n_clashes");
    EXPECT_EQ(result.metrics[1].name, "clashscore");
    EXPECT_EQ(result.metrics[2].name, "max_overlap");
    EXPECT_EQ(result.metrics[3].name, "hydrogens_present");
    EXPECT_EQ(metric(result, "n_clashes"), 1.0);
    EXPECT_DOUBLE_EQ(metric(result, "clashscore"), 500.0);
    EXPECT_NEAR(metric(result, "max_overlap"), 1.4, 1e-9);
    EXPECT_EQ(metric(result, "hydrogens_present"), 0.0);
}

TEST(ClashAssessment, ExactBoundariesAreInclusive) {
    const auto at_warning = assess_pair(OEChem::OEElemNo::C, OEChem::OEElemNo::C, 3.0);
    ASSERT_EQ(at_warning.findings.size(), 1U);
    EXPECT_EQ(at_warning.findings[0].severity, Severity::WARNING);

    const auto at_error = assess_pair(OEChem::OEElemNo::C, OEChem::OEElemNo::C, 2.4);
    ASSERT_EQ(at_error.findings.size(), 1U);
    EXPECT_EQ(at_error.findings[0].severity, Severity::ERROR);

    const auto beyond = assess_pair(OEChem::OEElemNo::C, OEChem::OEElemNo::C, 3.01);
    EXPECT_TRUE(beyond.findings.empty());
    EXPECT_EQ(beyond.status, Status::PASS);
    EXPECT_EQ(metric(beyond, "max_overlap"), 0.0);
}

TEST(ClashAssessment, PolarPairsUsePolarContactDistance) {
    const auto hydrogen_bond = assess_pair(OEChem::OEElemNo::N, OEChem::OEElemNo::O, 2.8);
    EXPECT_TRUE(hydrogen_bond.findings.empty());

    const auto close = assess_pair(OEChem::OEElemNo::N, OEChem::OEElemNo::O, 2.4);
    ASSERT_EQ(close.findings.size(), 1U);
    EXPECT_EQ(close.findings[0].severity, Severity::WARNING);
    EXPECT_NEAR(close.findings[0].expected, 2.9, 1e-9);
    EXPECT_NEAR(close.findings[0].deviation, 0.5, 1e-6);
}

TEST(ClashAssessment, MetalsAreExcluded) {
    const auto result = assess_pair(OEChem::OEElemNo::Zn, OEChem::OEElemNo::N, 2.1);
    EXPECT_TRUE(result.findings.empty());
}

TEST(ClashAssessment, AtomsWithinThreeBondsAreExcluded) {
    OEChem::OEGraphMol mol;
    OEChem::OEAtomBase* c1 = add_atom(mol, OEChem::OEElemNo::C, "C1", "LIG", 1, 0.0, 0.0, 0.0);
    OEChem::OEAtomBase* c2 = add_atom(mol, OEChem::OEElemNo::C, "C2", "LIG", 1, 1.5, 0.0, 0.0);
    OEChem::OEAtomBase* c3 = add_atom(mol, OEChem::OEElemNo::C, "C3", "LIG", 1, 1.5, 1.5, 0.0);
    OEChem::OEAtomBase* c4 = add_atom(mol, OEChem::OEElemNo::C, "C4", "LIG", 1, 0.0, 2.0, 0.0);
    mol.NewBond(c1, c2, 1);
    mol.NewBond(c2, c3, 1);
    mol.NewBond(c3, c4, 1);
    // C1 and C4 are 2.0 A apart, which would be an ERROR clash if not 1-4.
    const AssessmentResult result = ClashAssessment().Run(mol);
    EXPECT_TRUE(result.findings.empty());
}

TEST(ClashAssessment, RecordsWhetherHydrogensArePresent) {
    OEChem::OEGraphMol mol;
    add_atom(mol, OEChem::OEElemNo::C, "C1", "LIG", 1, 0.0, 0.0, 0.0);
    add_atom(mol, OEChem::OEElemNo::H, "H1", "LIG", 1, 10.0, 0.0, 0.0);
    EXPECT_EQ(metric(ClashAssessment().Run(mol), "hydrogens_present"), 1.0);
}

TEST(ClashAssessment, CrambinHasNoErrorClashes) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const AssessmentResult result = ClashAssessment().Run(mol);
    for (const Finding& finding : result.findings) {
        EXPECT_NE(finding.severity, Severity::ERROR) << finding.message;
    }
    EXPECT_NE(result.status, Status::FAIL);
    EXPECT_EQ(metric(result, "hydrogens_present"), 0.0);
}

TEST(ClashAssessment, NothingInScopeIsNotApplicable) {
    const OEChem::OEGraphMol empty;
    const AssessmentResult none = ClashAssessment().Run(empty);
    EXPECT_EQ(none.status, Status::NOT_APPLICABLE);
    EXPECT_FALSE(none.message.empty());
    EXPECT_EQ(metric(none, "hydrogens_present"), 0.0);

    OEChem::OEGraphMol metals;
    add_atom(metals, OEChem::OEElemNo::Zn, "ZN", "ZN", 1, 0.0, 0.0, 0.0);
    EXPECT_EQ(ClashAssessment().Run(metals).status, Status::NOT_APPLICABLE);
}

} // namespace OEHEALTH::test
```

- [ ] **Step 2: Register the test and confirm it fails**

Add to `add_executable(oehealth_tests ...)` in `tests/cpp/CMakeLists.txt`, keeping the list sorted:

```cmake
    test_clashes.cpp
```

Run: `cmake --build build-debug -j8`
Expected: FAIL to compile with `'oehealth/assessments/clashes.h' file not found`.

- [ ] **Step 3: Write the assessment**

Create the directories if an earlier task has not:

```bash
mkdir -p include/oehealth/assessments src/assessments
```

`include/oehealth/assessments/clashes.h`:

```cpp
#ifndef OEHEALTH_ASSESSMENTS_CLASHES_H
#define OEHEALTH_ASSESSMENTS_CLASHES_H

#include "oehealth/assessment.h"

namespace OEHEALTH {

/// \brief Overlap thresholds and the polar contact distance, in angstroms.
struct ClashOptions {
    double warning_overlap = 0.4;
    double error_overlap = 1.0;
    double polar_contact_distance = 2.9; ///< Contact distance for N/O pairs.
};

/// \brief Finds pairs of atoms closer than their contact distance.
///
/// The contact distance is the sum of Bondi radii, or polar_contact_distance
/// when both atoms are N or O. Atoms within three bonds of each other,
/// atoms with atomic number 0, and metals are excluded. The atoms are
/// assessed as given; the result records whether hydrogens were present.
class ClashAssessment : public Assessment {
public:
    /// \throws std::invalid_argument for NaN, infinite, or negative values,
    ///         a zero polar_contact_distance, or error_overlap < warning_overlap.
    explicit ClashAssessment(const ClashOptions& options = ClashOptions());

    const ClashOptions& GetOptions() const;
    std::string GetName() const override;
    std::string GetDescription() const override;
    AssessmentResult Assess(const StructureContext& ctx) const override;

private:
    ClashOptions options_;
};

} // namespace OEHEALTH

#endif // OEHEALTH_ASSESSMENTS_CLASHES_H
```

`src/assessments/clashes.cpp`:

```cpp
#include "oehealth/assessments/clashes.h"

#include "common.h"
#include "geometry/measure.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace OEHEALTH {

namespace {

// Margin added to the neighbor-search cutoff so pairs exactly at the warning
// threshold are candidates whether or not OENearestNbrs treats its cutoff as
// inclusive. The overlap test decides which candidates are findings.
constexpr double SEARCH_MARGIN = 0.01;

bool is_polar(const OEChem::OEAtomBase* atom) {
    const unsigned int element = atom->GetAtomicNum();
    return element == OEChem::OEElemNo::N || element == OEChem::OEElemNo::O;
}

bool within_three_bonds(const OEChem::OEAtomBase* a, const OEChem::OEAtomBase* b) {
    for (OESystem::OEIter<OEChem::OEAtomBase> first = a->GetAtoms(); first; ++first) {
        if (first->GetIdx() == b->GetIdx()) {
            return true;
        }
        for (OESystem::OEIter<OEChem::OEAtomBase> second = first->GetAtoms(); second; ++second) {
            if (second->GetIdx() == a->GetIdx()) {
                continue;
            }
            if (second->GetIdx() == b->GetIdx()) {
                return true;
            }
            for (OESystem::OEIter<OEChem::OEAtomBase> third = second->GetAtoms(); third; ++third) {
                if (third->GetIdx() == b->GetIdx()) {
                    return true;
                }
            }
        }
    }
    return false;
}

struct Clash {
    const OEChem::OEAtomBase* first;
    const OEChem::OEAtomBase* second;
    double distance;
    double contact;
    double overlap;
    Severity severity;
};

} // namespace

ClashAssessment::ClashAssessment(const ClashOptions& options) : options_(options) {
    validate_thresholds("clash overlap", options_.warning_overlap, options_.error_overlap);
    if (!std::isfinite(options_.polar_contact_distance) || options_.polar_contact_distance <= 0.0) {
        throw std::invalid_argument("polar_contact_distance must be finite and positive");
    }
}

const ClashOptions& ClashAssessment::GetOptions() const { return options_; }

std::string ClashAssessment::GetName() const { return "clashes"; }

std::string ClashAssessment::GetDescription() const {
    return "Finds non-bonded atom pairs that overlap by more than the clash threshold.";
}

AssessmentResult ClashAssessment::Assess(const StructureContext& ctx) const {
    const OEChem::OEMolBase& mol = ctx.GetMolecule();
    AssessmentResult result;
    result.assessment = GetName();
    const double hydrogens_present = ctx.HasHydrogens() ? 1.0 : 0.0;

    std::vector<double> radius(mol.GetMaxAtomIdx(), 0.0);
    std::vector<bool> in_scope(mol.GetMaxAtomIdx(), false);
    std::size_t n_in_scope = 0;
    double max_radius = 0.0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        const unsigned int element = atom->GetAtomicNum();
        if (element == 0 || atom->IsMetal()) {
            continue;
        }
        in_scope[atom->GetIdx()] = true;
        radius[atom->GetIdx()] = OEChem::OEGetBondiVdWRadius(element);
        max_radius = std::max(max_radius, radius[atom->GetIdx()]);
        ++n_in_scope;
    }

    if (n_in_scope == 0) {
        result.status = Status::NOT_APPLICABLE;
        result.message = "no atoms in scope for clash detection";
        result.metrics = {{"hydrogens_present", hydrogens_present}};
        return result;
    }

    std::vector<Clash> clashes;
    const double cutoff = std::max(2.0 * max_radius, options_.polar_contact_distance) -
                          options_.warning_overlap + SEARCH_MARGIN;
    if (cutoff > 0.0) {
        for (OESystem::OEIter<OEChem::OENbrs> pair = OEChem::OEGetNearestNbrs(mol, cutoff); pair;
             ++pair) {
            const OEChem::OEAtomBase* a = pair->GetBgn();
            const OEChem::OEAtomBase* b = pair->GetEnd();
            if (!in_scope[a->GetIdx()] || !in_scope[b->GetIdx()] || within_three_bonds(a, b)) {
                continue;
            }
            if (a->GetIdx() > b->GetIdx()) {
                std::swap(a, b);
            }
            const double contact = is_polar(a) && is_polar(b)
                                       ? options_.polar_contact_distance
                                       : radius[a->GetIdx()] + radius[b->GetIdx()];
            const double d = geometry::distance(mol, a, b);
            const double overlap = contact - d;
            if (meets_threshold(overlap, options_.error_overlap)) {
                clashes.push_back({a, b, d, contact, overlap, Severity::ERROR});
            } else if (meets_threshold(overlap, options_.warning_overlap)) {
                clashes.push_back({a, b, d, contact, overlap, Severity::WARNING});
            }
        }
    }

    std::sort(clashes.begin(), clashes.end(), [](const Clash& x, const Clash& y) {
        return std::make_pair(x.first->GetIdx(), x.second->GetIdx()) <
               std::make_pair(y.first->GetIdx(), y.second->GetIdx());
    });

    double max_overlap = 0.0;
    for (const Clash& clash : clashes) {
        max_overlap = std::max(max_overlap, clash.overlap);
        Finding finding;
        finding.code = "clash";
        finding.severity = clash.severity;
        finding.atoms = {make_atom_ref(clash.first), make_atom_ref(clash.second)};
        finding.observed = clash.distance;
        finding.expected = clash.contact;
        finding.deviation = clash.overlap;
        finding.message = atom_path(finding.atoms) + " " + format_fixed(clash.distance, 3) +
                          " A apart (contact " + format_fixed(clash.contact, 3) + ", overlap " +
                          format_fixed(clash.overlap, 3) + ")";
        result.findings.push_back(std::move(finding));
    }

    const auto n_clashes = static_cast<double>(result.findings.size());
    result.metrics = {
        {"n_clashes", n_clashes},
        {"clashscore", n_clashes * 1000.0 / static_cast<double>(n_in_scope)},
        {"max_overlap", max_overlap},
        {"hydrogens_present", hydrogens_present},
    };
    result.status = status_from_findings(result.findings);
    return result;
}

} // namespace OEHEALTH
```

- [ ] **Step 4: Register the source**

Add to `set(OEHEALTH_SOURCES ...)` in `CMakeLists.txt`, keeping the list sorted:

```cmake
    src/assessments/clashes.cpp
```

- [ ] **Step 5: Build and run the tests**

```bash
cmake --build build-debug -j8
./build-debug/tests/cpp/oehealth_tests --gtest_filter='ClashAssessment.*'
ctest --test-dir build-debug --output-on-failure
```

Expected: `[  PASSED  ] 10 tests.` from the filtered run; ctest passes every test.

- [ ] **Step 6: Format and lint**

```bash
clang-format --dry-run --Werror include/oehealth/assessments/clashes.h src/assessments/clashes.cpp tests/cpp/test_clashes.cpp
clang-tidy -p build-debug --quiet --extra-arg=-isysroot --extra-arg="$(xcrun --show-sdk-path)" src/assessments/clashes.cpp tests/cpp/test_clashes.cpp
```

Expected: no output from clang-format; no findings from clang-tidy.

- [ ] **Step 7: Commit**

```bash
git add include/oehealth/assessments/clashes.h src/assessments/clashes.cpp tests/cpp/test_clashes.cpp CMakeLists.txt tests/cpp/CMakeLists.txt
git commit -m "Add the clash assessment"
```

---

### Task 8: Valence assessment

**Risk tier:** standard — a public assessment, three new files plus CMake registration.

**Depends on:** Task 2

**Files:**
- Create: `include/oehealth/assessments/valence.h`, `src/assessments/valence.cpp`
- Test: `tests/cpp/test_valence.cpp`
- Modify: `CMakeLists.txt` (`OEHEALTH_SOURCES`), `tests/cpp/CMakeLists.txt`

**Interfaces:**
- Consumes: `Assessment`, `StructureContext`, result types (Task 2); `make_atom_ref`, `atom_path`, `status_from_findings` (Task 2, `src/common.h`).
- Produces: `class ValenceAssessment : public Assessment { ValenceAssessment(); }`, name `"valence"`, finding code `"invalid_valence"` with details `formal_charge`, `implicit_h_count`, `degree`. Metric: `n_invalid`.

- [ ] **Step 1: Write the failing tests**

`tests/cpp/test_valence.cpp`:

```cpp
#include <gtest/gtest.h>

#include "oehealth/assessments/valence.h"
#include "test_helpers.h"

namespace OEHEALTH::test {

TEST(ValenceAssessment, Name) {
    const ValenceAssessment assessment;
    EXPECT_EQ(assessment.GetName(), "valence");
    EXPECT_FALSE(assessment.GetDescription().empty());
}

TEST(ValenceAssessment, PentavalentCarbonIsReported) {
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(OEChem::OESmilesToMol(mol, "C(C)(C)(C)(C)C"));
    const StructureContext ctx(mol);
    const AssessmentResult result = ValenceAssessment().Assess(ctx);
    ASSERT_EQ(result.findings.size(), 1U);
    const Finding& finding = result.findings[0];
    EXPECT_EQ(finding.code, "invalid_valence");
    EXPECT_EQ(finding.severity, Severity::ERROR);
    ASSERT_EQ(finding.atoms.size(), 1U);
    EXPECT_EQ(finding.observed, 5.0);
    ASSERT_EQ(finding.details.size(), 3U);
    EXPECT_EQ(finding.details[0].name, "formal_charge");
    EXPECT_EQ(finding.details[1].name, "implicit_h_count");
    EXPECT_EQ(finding.details[2].name, "degree");
    EXPECT_EQ(finding.details[2].value, 5.0);
    EXPECT_EQ(result.status, Status::FAIL);
    EXPECT_EQ(metric(result, "n_invalid"), 1.0);
}

TEST(ValenceAssessment, CrambinAsReadIsValid) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const AssessmentResult result = ValenceAssessment().Run(mol);
    EXPECT_TRUE(result.findings.empty());
    EXPECT_EQ(result.status, Status::PASS);
    EXPECT_EQ(metric(result, "n_invalid"), 0.0);
}

TEST(ValenceAssessment, ZeroedImplicitHydrogensAreReported) {
    // Documents the requirement that implicit hydrogen counts be assigned.
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        atom->SetImplicitHCount(0);
    }
    EXPECT_GT(ValenceAssessment().Run(mol).findings.size(), 0U);
}

TEST(ValenceAssessment, EmptyMoleculeIsNotApplicable) {
    const OEChem::OEGraphMol mol;
    const AssessmentResult result = ValenceAssessment().Run(mol);
    EXPECT_EQ(result.status, Status::NOT_APPLICABLE);
    EXPECT_FALSE(result.message.empty());
}

} // namespace OEHEALTH::test
```

- [ ] **Step 2: Register the test and confirm it fails**

Add to `add_executable(oehealth_tests ...)` in `tests/cpp/CMakeLists.txt`, keeping the list sorted:

```cmake
    test_valence.cpp
```

Run: `cmake --build build-debug -j8`
Expected: FAIL to compile with `'oehealth/assessments/valence.h' file not found`.

- [ ] **Step 3: Write the assessment**

Create the directories if an earlier task has not:

```bash
mkdir -p include/oehealth/assessments src/assessments
```

`include/oehealth/assessments/valence.h`:

```cpp
#ifndef OEHEALTH_ASSESSMENTS_VALENCE_H
#define OEHEALTH_ASSESSMENTS_VALENCE_H

#include "oehealth/assessment.h"

namespace OEHEALTH {

/// \brief Reports atoms whose valence OEChem considers invalid.
///
/// Runs OECheckAtomValence on every atom; each failure is an ERROR finding,
/// code "invalid_valence". Implicit hydrogen counts must be assigned, as
/// OpenEye's readers do by default: with zeroed counts, normal atoms are
/// reported as under-valent.
class ValenceAssessment : public Assessment {
public:
    ValenceAssessment();

    std::string GetName() const override;
    std::string GetDescription() const override;
    AssessmentResult Assess(const StructureContext& ctx) const override;
};

} // namespace OEHEALTH

#endif // OEHEALTH_ASSESSMENTS_VALENCE_H
```

`src/assessments/valence.cpp`:

```cpp
#include "oehealth/assessments/valence.h"

#include "common.h"

#include <utility>

namespace OEHEALTH {

ValenceAssessment::ValenceAssessment() = default;

std::string ValenceAssessment::GetName() const { return "valence"; }

std::string ValenceAssessment::GetDescription() const {
    return "Reports atoms whose valence is invalid for their element and formal charge.";
}

AssessmentResult ValenceAssessment::Assess(const StructureContext& ctx) const {
    const OEChem::OEMolBase& mol = ctx.GetMolecule();
    AssessmentResult result;
    result.assessment = GetName();
    if (mol.NumAtoms() == 0) {
        result.status = Status::NOT_APPLICABLE;
        result.message = "the molecule has no atoms";
        result.metrics = {{"n_invalid", 0.0}};
        return result;
    }

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (OEChem::OECheckAtomValence(atom)) {
            continue;
        }
        Finding finding;
        finding.code = "invalid_valence";
        finding.severity = Severity::ERROR;
        finding.atoms = {make_atom_ref(atom)};
        finding.observed = static_cast<double>(atom->GetValence());
        finding.details = {
            {"formal_charge", static_cast<double>(atom->GetFormalCharge())},
            {"implicit_h_count", static_cast<double>(atom->GetImplicitHCount())},
            {"degree", static_cast<double>(atom->GetDegree())},
        };
        finding.message = atom_path(finding.atoms) + " " +
                          OEChem::OEGetAtomicSymbol(atom->GetAtomicNum()) + " with formal charge " +
                          std::to_string(atom->GetFormalCharge()) + " has invalid valence " +
                          std::to_string(atom->GetValence());
        result.findings.push_back(std::move(finding));
    }

    result.metrics = {{"n_invalid", static_cast<double>(result.findings.size())}};
    result.status = status_from_findings(result.findings);
    return result;
}

} // namespace OEHEALTH
```

- [ ] **Step 4: Register the source**

Add to `set(OEHEALTH_SOURCES ...)` in `CMakeLists.txt`, keeping the list sorted:

```cmake
    src/assessments/valence.cpp
```

- [ ] **Step 5: Build and run the tests**

```bash
cmake --build build-debug -j8
./build-debug/tests/cpp/oehealth_tests --gtest_filter='ValenceAssessment.*'
ctest --test-dir build-debug --output-on-failure
```

Expected: `[  PASSED  ] 5 tests.` from the filtered run; ctest passes every test.

- [ ] **Step 6: Format and lint**

```bash
clang-format --dry-run --Werror include/oehealth/assessments/valence.h src/assessments/valence.cpp tests/cpp/test_valence.cpp
clang-tidy -p build-debug --quiet --extra-arg=-isysroot --extra-arg="$(xcrun --show-sdk-path)" src/assessments/valence.cpp tests/cpp/test_valence.cpp
```

Expected: no output from clang-format; no findings from clang-tidy.

- [ ] **Step 7: Commit**

```bash
git add include/oehealth/assessments/valence.h src/assessments/valence.cpp tests/cpp/test_valence.cpp CMakeLists.txt tests/cpp/CMakeLists.txt
git commit -m "Add the valence assessment"
```

---

### Task 9: Ramachandran assessment

**Risk tier:** standard — a public assessment that wraps OEBio's classifier, three new files plus CMake registration.

**Depends on:** Task 2

**Files:**
- Create: `include/oehealth/assessments/ramachandran.h`, `src/assessments/ramachandran.cpp`
- Test: `tests/cpp/test_ramachandran.cpp`
- Modify: `CMakeLists.txt` (`OEHEALTH_SOURCES`), `tests/cpp/CMakeLists.txt`

**Interfaces:**
- Consumes: `Assessment`, `StructureContext`, `ResidueInfo` (Task 2); `make_atom_ref`, `residue_label`, `format_fixed`, `status_from_findings` (Task 2, `src/common.h`); `OEBio::OERamachandranAnalysis` (OEBio, linked in Task 1).
- Produces: `class RamachandranAssessment : public Assessment { RamachandranAssessment(); }`, name `"ramachandran"`, finding code `"ramachandran_outlier"` with details `phi`, `psi` (degrees), `score`, `rama_type`. Metrics, in order: `n_assessed`, `n_favored`, `n_allowed`, `n_outliers`, `pct_favored`, `pct_outliers`; on NOT_APPLICABLE only `n_assessed` = 0.

- [ ] **Step 1: Write the failing tests**

`tests/cpp/test_ramachandran.cpp`:

```cpp
#include <gtest/gtest.h>

#include "oehealth/assessments/ramachandran.h"
#include "test_helpers.h"

namespace OEHEALTH::test {

TEST(RamachandranAssessment, Name) {
    const RamachandranAssessment assessment;
    EXPECT_EQ(assessment.GetName(), "ramachandran");
    EXPECT_FALSE(assessment.GetDescription().empty());
}

TEST(RamachandranAssessment, CrambinMatchesWwpdbReport) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const AssessmentResult result = RamachandranAssessment().Run(mol);
    EXPECT_TRUE(result.findings.empty());
    EXPECT_EQ(result.status, Status::PASS);
    ASSERT_EQ(result.metrics.size(), 6U);
    EXPECT_EQ(result.metrics[0].name, "n_assessed");
    EXPECT_EQ(result.metrics[1].name, "n_favored");
    EXPECT_EQ(result.metrics[2].name, "n_allowed");
    EXPECT_EQ(result.metrics[3].name, "n_outliers");
    EXPECT_EQ(result.metrics[4].name, "pct_favored");
    EXPECT_EQ(result.metrics[5].name, "pct_outliers");
    // The two termini cannot be classified.
    EXPECT_EQ(metric(result, "n_assessed"), 44.0);
    EXPECT_EQ(metric(result, "n_outliers"), 0.0);
    EXPECT_EQ(metric(result, "n_favored") + metric(result, "n_allowed"), 44.0);
    EXPECT_DOUBLE_EQ(metric(result, "pct_favored"), 100.0 * metric(result, "n_favored") / 44.0);
    EXPECT_EQ(metric(result, "pct_outliers"), 0.0);
}

TEST(RamachandranAssessment, DisallowedTorsionsGiveOneOutlier) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    // Residue 44 lies outside every disulfide ring, so its phi and psi can be set.
    // phi = psi = 0 is deep in the disallowed region for a general residue.
    ASSERT_TRUE(OEChem::OESetTorsion(mol, find_atom(mol, 43, "C"), find_atom(mol, 44, "N"),
                                     find_atom(mol, 44, "CA"), find_atom(mol, 44, "C"), 0.0));
    ASSERT_TRUE(OEChem::OESetTorsion(mol, find_atom(mol, 44, "N"), find_atom(mol, 44, "CA"),
                                     find_atom(mol, 44, "C"), find_atom(mol, 45, "N"), 0.0));
    const AssessmentResult result = RamachandranAssessment().Run(mol);
    ASSERT_EQ(result.findings.size(), 1U);
    const Finding& finding = result.findings[0];
    EXPECT_EQ(finding.code, "ramachandran_outlier");
    EXPECT_EQ(finding.severity, Severity::WARNING);
    ASSERT_EQ(finding.atoms.size(), 3U);
    EXPECT_EQ(finding.atoms[0].atom_name, "N");
    EXPECT_EQ(finding.atoms[1].atom_name, "CA");
    EXPECT_EQ(finding.atoms[2].atom_name, "C");
    EXPECT_EQ(finding.atoms[1].residue_number, 44);
    ASSERT_EQ(finding.details.size(), 4U);
    EXPECT_EQ(finding.details[0].name, "phi");
    EXPECT_EQ(finding.details[1].name, "psi");
    EXPECT_EQ(finding.details[2].name, "score");
    EXPECT_EQ(finding.details[3].name, "rama_type");
    EXPECT_NEAR(finding.details[0].value, 0.0, 0.01);
    EXPECT_NEAR(finding.details[1].value, 0.0, 0.01);
    EXPECT_EQ(finding.observed, finding.details[2].value);
    EXPECT_EQ(result.status, Status::WARN);
    EXPECT_EQ(metric(result, "n_outliers"), 1.0);
}

TEST(RamachandranAssessment, NothingToClassifyIsNotApplicable) {
    const OEChem::OEGraphMol empty;
    EXPECT_EQ(RamachandranAssessment().Run(empty).status, Status::NOT_APPLICABLE);

    OEChem::OEGraphMol ethanol;
    OEChem::OESmilesToMol(ethanol, "CCO");
    ethanol.SetDimension(3);
    const AssessmentResult result = RamachandranAssessment().Run(ethanol);
    EXPECT_EQ(result.status, Status::NOT_APPLICABLE);
    EXPECT_FALSE(result.message.empty());
}

} // namespace OEHEALTH::test
```

- [ ] **Step 2: Register the test and confirm it fails**

Add to `add_executable(oehealth_tests ...)` in `tests/cpp/CMakeLists.txt`, keeping the list sorted:

```cmake
    test_ramachandran.cpp
```

Run: `cmake --build build-debug -j8`
Expected: FAIL to compile with `'oehealth/assessments/ramachandran.h' file not found`.

- [ ] **Step 3: Write the assessment**

Create the directories if an earlier task has not:

```bash
mkdir -p include/oehealth/assessments src/assessments
```

`include/oehealth/assessments/ramachandran.h`:

```cpp
#ifndef OEHEALTH_ASSESSMENTS_RAMACHANDRAN_H
#define OEHEALTH_ASSESSMENTS_RAMACHANDRAN_H

#include "oehealth/assessment.h"

namespace OEHEALTH {

/// \brief Classifies standard amino-acid backbone torsions with OEBio's
/// Ramachandran analysis.
///
/// Residues without a CA atom, and residues OEBio cannot classify (termini
/// and chain breaks), are not assessed. Each Outlier is a WARNING finding,
/// code "ramachandran_outlier".
class RamachandranAssessment : public Assessment {
public:
    RamachandranAssessment();

    std::string GetName() const override;
    std::string GetDescription() const override;
    AssessmentResult Assess(const StructureContext& ctx) const override;
};

} // namespace OEHEALTH

#endif // OEHEALTH_ASSESSMENTS_RAMACHANDRAN_H
```

`src/assessments/ramachandran.cpp`:

```cpp
#include "oehealth/assessments/ramachandran.h"

#include "common.h"

#include <oebio.h>

#include <utility>

namespace OEHEALTH {

namespace {

constexpr double RADIANS_TO_DEGREES = 180.0 / 3.14159265358979323846;

} // namespace

RamachandranAssessment::RamachandranAssessment() = default;

std::string RamachandranAssessment::GetName() const { return "ramachandran"; }

std::string RamachandranAssessment::GetDescription() const {
    return "Classifies backbone phi/psi torsions of standard amino acids as favored, allowed, or "
           "outlier.";
}

AssessmentResult RamachandranAssessment::Assess(const StructureContext& ctx) const {
    const OEChem::OEMolBase& mol = ctx.GetMolecule();
    AssessmentResult result;
    result.assessment = GetName();

    std::size_t n_favored = 0;
    std::size_t n_allowed = 0;
    std::size_t n_outliers = 0;
    for (const ResidueInfo& residue : ctx.GetResidues()) {
        if (!residue.is_standard_amino_acid) {
            continue;
        }
        const OEChem::OEAtomBase* alpha = residue.FindAtom("CA");
        if (alpha == nullptr) {
            continue;
        }
        const OEBio::OERamachandranAnalysis analysis(mol, OEChem::OEAtomGetResidue(alpha));
        const unsigned int category = analysis.GetRamaCategory();
        if (category == OEBio::OERamaCategory::Favored) {
            ++n_favored;
        } else if (category == OEBio::OERamaCategory::Allowed) {
            ++n_allowed;
        } else if (category == OEBio::OERamaCategory::Outlier) {
            ++n_outliers;
            const double phi = analysis.GetPhi() * RADIANS_TO_DEGREES;
            const double psi = analysis.GetPsi() * RADIANS_TO_DEGREES;
            Finding finding;
            finding.code = "ramachandran_outlier";
            finding.severity = Severity::WARNING;
            for (const char* name : {"N", "CA", "C"}) {
                if (const OEChem::OEAtomBase* atom = residue.FindAtom(name)) {
                    finding.atoms.push_back(make_atom_ref(atom));
                }
            }
            finding.observed = analysis.GetRamaScore();
            finding.details = {
                {"phi", phi},
                {"psi", psi},
                {"score", analysis.GetRamaScore()},
                {"rama_type", static_cast<double>(analysis.GetRamaType())},
            };
            finding.message = residue_label(make_atom_ref(alpha)) + " phi " + format_fixed(phi, 1) +
                              " psi " + format_fixed(psi, 1) + " (" +
                              OEBio::OEGetRamachandranTypeName(analysis.GetRamaType()) +
                              ") is a Ramachandran outlier";
            result.findings.push_back(std::move(finding));
        }
    }

    const std::size_t n_assessed = n_favored + n_allowed + n_outliers;
    if (n_assessed == 0) {
        result.status = Status::NOT_APPLICABLE;
        result.message = "no standard amino-acid residues could be classified";
        result.metrics = {{"n_assessed", 0.0}};
        return result;
    }
    const auto total = static_cast<double>(n_assessed);
    result.metrics = {
        {"n_assessed", total},
        {"n_favored", static_cast<double>(n_favored)},
        {"n_allowed", static_cast<double>(n_allowed)},
        {"n_outliers", static_cast<double>(n_outliers)},
        {"pct_favored", 100.0 * static_cast<double>(n_favored) / total},
        {"pct_outliers", 100.0 * static_cast<double>(n_outliers) / total},
    };
    result.status = status_from_findings(result.findings);
    return result;
}

} // namespace OEHEALTH
```

- [ ] **Step 4: Register the source**

Add to `set(OEHEALTH_SOURCES ...)` in `CMakeLists.txt`, keeping the list sorted:

```cmake
    src/assessments/ramachandran.cpp
```

- [ ] **Step 5: Build and run the tests**

```bash
cmake --build build-debug -j8
./build-debug/tests/cpp/oehealth_tests --gtest_filter='RamachandranAssessment.*'
ctest --test-dir build-debug --output-on-failure
```

Expected: `[  PASSED  ] 4 tests.` from the filtered run; ctest passes every test.

- [ ] **Step 6: Format and lint**

```bash
clang-format --dry-run --Werror include/oehealth/assessments/ramachandran.h src/assessments/ramachandran.cpp tests/cpp/test_ramachandran.cpp
clang-tidy -p build-debug --quiet --extra-arg=-isysroot --extra-arg="$(xcrun --show-sdk-path)" src/assessments/ramachandran.cpp tests/cpp/test_ramachandran.cpp
```

Expected: no output from clang-format; no findings from clang-tidy.

- [ ] **Step 7: Commit**

```bash
git add include/oehealth/assessments/ramachandran.h src/assessments/ramachandran.cpp tests/cpp/test_ramachandran.cpp CMakeLists.txt tests/cpp/CMakeLists.txt
git commit -m "Add the Ramachandran assessment"
```

---

### Task 10: JSON serialization

**Risk tier:** standard — the public output format, with a private third-party dependency, three new files.

**Depends on:** Task 2, Task 6

**Files:**
- Create: `include/oehealth/json.h`, `src/json.cpp`
- Test: `tests/cpp/test_json.cpp`
- Modify: `CMakeLists.txt` (`OEHEALTH_SOURCES`), `tests/cpp/CMakeLists.txt`

**Interfaces:**
- Consumes: result types (Task 2); `nlohmann_json::nlohmann_json` (Task 1); `OEHEALTH_VERSION_MAJOR`, `_MINOR`, `_PATCH` from `oehealth/oehealth.h` (the template's header still provides them at this point); `BondLengthAssessment` (Task 6, used by the round-trip test only).
- Produces: `constexpr int JSON_SCHEMA_VERSION = 1` and `std::string to_json(const AssessmentResult& result, int indent = 2)`, which throws `std::invalid_argument` when a metric or detail name repeats within one vector.

- [ ] **Step 1: Write the failing tests**

`tests/cpp/test_json.cpp`:

```cpp
#include <gtest/gtest.h>

#include "oehealth/assessments/bond_lengths.h"
#include "oehealth/json.h"
#include "oehealth/oehealth.h"
#include "test_helpers.h"

#include <nlohmann/json.hpp>

#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace OEHEALTH::test {

namespace {

using Json = nlohmann::ordered_json;

AssessmentResult sample_result() {
    AtomRef ca;
    ca.atom_idx = 10;
    ca.atom_name = "CA";
    ca.residue_name = "ALA";
    ca.chain_id = "A";
    ca.residue_number = 5;

    Finding finding;
    finding.code = "bond_length_outlier";
    finding.severity = Severity::WARNING;
    finding.message = "A ALA 5 CA-CB";
    finding.atoms = {ca};
    finding.observed = 1.617;
    finding.expected = 1.520;
    finding.deviation = std::numeric_limits<double>::quiet_NaN();
    finding.details = {{"phi", -60.0}, {"psi", -45.0}};

    AssessmentResult child;
    child.assessment = "child";
    child.status = Status::PASS;

    AssessmentResult result;
    result.assessment = "bond_lengths";
    result.status = Status::WARN;
    result.metrics = {{"n_assessed", 412.0}, {"n_outliers", 1.0}, {"rmsz", 0.83}};
    result.findings = {finding};
    result.children = {child};
    return result;
}

std::vector<std::string> keys(const Json& object) {
    std::vector<std::string> names;
    for (const auto& item : object.items()) {
        names.push_back(item.key());
    }
    return names;
}

} // namespace

TEST(ToJson, RootLayoutAndOrder) {
    const Json root = Json::parse(to_json(sample_result()));
    EXPECT_EQ(keys(root),
              (std::vector<std::string>{"schema_version", "oehealth_version", "assessment",
                                        "status", "message", "metrics", "findings", "children"}));
    EXPECT_EQ(root["schema_version"], 1);
    EXPECT_EQ(root["oehealth_version"], std::to_string(OEHEALTH_VERSION_MAJOR) + "." +
                                            std::to_string(OEHEALTH_VERSION_MINOR) + "." +
                                            std::to_string(OEHEALTH_VERSION_PATCH));
    EXPECT_EQ(root["status"], "WARN");
    EXPECT_EQ(keys(root["metrics"]),
              (std::vector<std::string>{"n_assessed", "n_outliers", "rmsz"}));
    EXPECT_DOUBLE_EQ(root["metrics"]["rmsz"].get<double>(), 0.83);
}

TEST(ToJson, FindingLayout) {
    const Json finding = Json::parse(to_json(sample_result()))["findings"][0];
    EXPECT_EQ(keys(finding),
              (std::vector<std::string>{"code", "severity", "message", "atoms", "observed",
                                        "expected", "deviation", "details"}));
    EXPECT_EQ(finding["severity"], "WARNING");
    EXPECT_EQ(keys(finding["atoms"][0]),
              (std::vector<std::string>{"atom_idx", "atom_name", "residue_name", "residue_number",
                                        "insert_code", "chain_id"}));
    EXPECT_EQ(finding["atoms"][0]["atom_idx"], 10);
    EXPECT_EQ(finding["atoms"][0]["insert_code"], "");
    EXPECT_TRUE(finding["deviation"].is_null());
    EXPECT_EQ(keys(finding["details"]), (std::vector<std::string>{"phi", "psi"}));
}

TEST(ToJson, ChildrenOmitVersionKeys) {
    const Json child = Json::parse(to_json(sample_result()))["children"][0];
    EXPECT_FALSE(child.contains("schema_version"));
    EXPECT_FALSE(child.contains("oehealth_version"));
    EXPECT_EQ(child["assessment"], "child");
    EXPECT_EQ(child["status"], "PASS");
}

TEST(ToJson, DuplicateNamesThrow) {
    AssessmentResult duplicate_metric = sample_result();
    duplicate_metric.metrics.push_back({"rmsz", 1.0});
    EXPECT_THROW(to_json(duplicate_metric), std::invalid_argument);

    AssessmentResult duplicate_detail = sample_result();
    duplicate_detail.findings[0].details.push_back({"phi", 1.0});
    EXPECT_THROW(to_json(duplicate_detail), std::invalid_argument);
}

TEST(ToJson, NegativeIndentIsCompact) {
    const std::string compact = to_json(sample_result(), -1);
    EXPECT_EQ(compact.find('\n'), std::string::npos);
    EXPECT_NE(to_json(sample_result()).find('\n'), std::string::npos);
}

TEST(ToJson, RealResultRoundTrips) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const Json root = Json::parse(to_json(BondLengthAssessment().Run(mol)));
    EXPECT_EQ(root["assessment"], "bond_lengths");
    EXPECT_TRUE(root["metrics"].contains("n_assessed"));
}

} // namespace OEHEALTH::test
```

- [ ] **Step 2: Register the test and confirm it fails**

Add to `add_executable(oehealth_tests ...)` in `tests/cpp/CMakeLists.txt`, keeping the list sorted:

```cmake
    test_json.cpp
```

Run: `cmake --build build-debug -j8`
Expected: FAIL to compile with `'oehealth/json.h' file not found`.

- [ ] **Step 3: Write the serializer**

`include/oehealth/json.h`:

```cpp
#ifndef OEHEALTH_JSON_H
#define OEHEALTH_JSON_H

#include "oehealth/result.h"

#include <string>

namespace OEHEALTH {

/// \brief Version of the JSON layout written by to_json.
constexpr int JSON_SCHEMA_VERSION = 1;

/// \brief Serialize a result, and its children, to JSON.
///
/// The root object carries ``schema_version`` and ``oehealth_version``;
/// children omit them. Metrics and details are objects that keep their
/// order, NaN is written as null, and enums are written by name.
///
/// \param indent Spaces per level; a negative value gives compact output.
/// \throws std::invalid_argument if a metric or detail name repeats within
///         one vector.
std::string to_json(const AssessmentResult& result, int indent = 2);

} // namespace OEHEALTH

#endif // OEHEALTH_JSON_H
```

`src/json.cpp`:

```cpp
#include "oehealth/json.h"

#include "oehealth/oehealth.h"

#include <nlohmann/json.hpp>

#include <cmath>
#include <set>
#include <stdexcept>
#include <string>

namespace OEHEALTH {

namespace {

using Json = nlohmann::ordered_json;

const char* severity_name(Severity severity) {
    switch (severity) {
        case Severity::INFO:
            return "INFO";
        case Severity::WARNING:
            return "WARNING";
        case Severity::ERROR:
            return "ERROR";
    }
    return "UNKNOWN";
}

const char* status_name(Status status) {
    switch (status) {
        case Status::NOT_APPLICABLE:
            return "NOT_APPLICABLE";
        case Status::PASS:
            return "PASS";
        case Status::WARN:
            return "WARN";
        case Status::FAIL:
            return "FAIL";
        case Status::ERROR:
            return "ERROR";
    }
    return "UNKNOWN";
}

Json number(double value) {
    if (std::isnan(value)) {
        return nullptr;
    }
    return value;
}

Json metrics_object(const std::vector<Metric>& metrics, const std::string& where) {
    Json object = Json::object();
    std::set<std::string> seen;
    for (const Metric& metric : metrics) {
        if (!seen.insert(metric.name).second) {
            throw std::invalid_argument("duplicate name '" + metric.name + "' in " + where);
        }
        object[metric.name] = number(metric.value);
    }
    return object;
}

Json atom_object(const AtomRef& atom) {
    Json object = Json::object();
    object["atom_idx"] = atom.atom_idx;
    object["atom_name"] = atom.atom_name;
    object["residue_name"] = atom.residue_name;
    object["residue_number"] = atom.residue_number;
    object["insert_code"] = atom.insert_code;
    object["chain_id"] = atom.chain_id;
    return object;
}

Json finding_object(const Finding& finding) {
    Json object = Json::object();
    object["code"] = finding.code;
    object["severity"] = severity_name(finding.severity);
    object["message"] = finding.message;
    Json atoms = Json::array();
    for (const AtomRef& atom : finding.atoms) {
        atoms.push_back(atom_object(atom));
    }
    object["atoms"] = atoms;
    object["observed"] = number(finding.observed);
    object["expected"] = number(finding.expected);
    object["deviation"] = number(finding.deviation);
    object["details"] =
        metrics_object(finding.details, "details of finding '" + finding.code + "'");
    return object;
}

Json result_object(const AssessmentResult& result, bool is_root) {
    Json object = Json::object();
    if (is_root) {
        object["schema_version"] = JSON_SCHEMA_VERSION;
        object["oehealth_version"] = std::to_string(OEHEALTH_VERSION_MAJOR) + "." +
                                     std::to_string(OEHEALTH_VERSION_MINOR) + "." +
                                     std::to_string(OEHEALTH_VERSION_PATCH);
    }
    object["assessment"] = result.assessment;
    object["status"] = status_name(result.status);
    object["message"] = result.message;
    object["metrics"] = metrics_object(result.metrics, "metrics of '" + result.assessment + "'");
    Json findings = Json::array();
    for (const Finding& finding : result.findings) {
        findings.push_back(finding_object(finding));
    }
    object["findings"] = findings;
    Json children = Json::array();
    for (const AssessmentResult& child : result.children) {
        children.push_back(result_object(child, false));
    }
    object["children"] = children;
    return object;
}

} // namespace

std::string to_json(const AssessmentResult& result, int indent) {
    return result_object(result, true).dump(indent);
}

} // namespace OEHEALTH
```

- [ ] **Step 4: Register the source**

Add to `set(OEHEALTH_SOURCES ...)` in `CMakeLists.txt`, keeping the list sorted:

```cmake
    src/json.cpp
```

- [ ] **Step 5: Build and run the tests**

```bash
cmake --build build-debug -j8
./build-debug/tests/cpp/oehealth_tests --gtest_filter='ToJson.*'
ctest --test-dir build-debug --output-on-failure
```

Expected: `[  PASSED  ] 6 tests.` from the filtered run; ctest passes every test.

- [ ] **Step 6: Format and lint**

```bash
clang-format --dry-run --Werror include/oehealth/json.h src/json.cpp tests/cpp/test_json.cpp
clang-tidy -p build-debug --quiet --extra-arg=-isysroot --extra-arg="$(xcrun --show-sdk-path)" src/json.cpp tests/cpp/test_json.cpp
```

Expected: no output from clang-format; no findings from clang-tidy.

- [ ] **Step 7: Commit**

```bash
git add include/oehealth/json.h src/json.cpp tests/cpp/test_json.cpp CMakeLists.txt tests/cpp/CMakeLists.txt
git commit -m "Serialize assessment results to JSON"
```

---

### Task 11: Protein geometry suite and the umbrella header

**Risk tier:** standard — wires every assessment into one public suite and changes the public umbrella header.

**Depends on:** Task 3, Task 6, Task 7, Task 8, Task 9, Task 10

**Files:**
- Modify: `include/oehealth/suite.h`, `src/suite.cpp`, `include/oehealth/oehealth.h`
- Test: `tests/cpp/test_protein_suite.cpp`
- Modify: `tests/cpp/CMakeLists.txt`

**Interfaces:**
- Consumes: `AssessmentSuite` (Task 3); `BondLengthAssessment`, `BondAngleAssessment` (Task 6); `ClashAssessment` (Task 7); `ValenceAssessment` (Task 8); `RamachandranAssessment` (Task 9); `json.h` (Task 10).
- Produces: `std::shared_ptr<AssessmentSuite> make_protein_geometry_suite()`, a suite named `"protein_geometry"` with children `bond_lengths`, `bond_angles`, `clashes`, `valence`, `ramachandran` in that order; `oehealth/oehealth.h` includes every public header. The header keeps `calculate_molecular_weight` until Task 12, because the SWIG module still wraps it.

`test_protein_suite.cpp` defines `Version.MacrosDefined`, and the template's `test_oehealth.cpp` defines `VersionTest.MacrosDefined`. The names differ, so both build until Task 12 deletes the template test.

- [ ] **Step 1: Write the failing tests**

`tests/cpp/test_protein_suite.cpp`:

```cpp
#include <gtest/gtest.h>

#include "oehealth/oehealth.h"
#include "test_helpers.h"

#include <string>
#include <vector>

namespace OEHEALTH::test {

TEST(Version, MacrosDefined) {
    EXPECT_GE(OEHEALTH_VERSION_MAJOR, 0);
    EXPECT_GE(OEHEALTH_VERSION_MINOR, 0);
    EXPECT_GE(OEHEALTH_VERSION_PATCH, 0);
}

TEST(ProteinGeometrySuite, HasTheFiveAssessmentsInOrder) {
    const auto suite = make_protein_geometry_suite();
    EXPECT_EQ(suite->GetName(), "protein_geometry");
    EXPECT_FALSE(suite->GetDescription().empty());
    std::vector<std::string> names;
    for (const auto& child : suite->GetChildren()) {
        names.push_back(child->GetName());
    }
    EXPECT_EQ(names, (std::vector<std::string>{"bond_lengths", "bond_angles", "clashes", "valence",
                                               "ramachandran"}));
}

TEST(ProteinGeometrySuite, CrambinHasNoErrors) {
    OEChem::OEGraphMol mol;
    read_crambin(mol);
    const AssessmentResult result = make_protein_geometry_suite()->Run(mol);
    ASSERT_EQ(result.children.size(), 5U);
    for (const AssessmentResult& child : result.children) {
        EXPECT_NE(child.status, Status::ERROR) << child.assessment << ": " << child.message;
        EXPECT_NE(child.status, Status::FAIL) << child.assessment;
        for (const Finding& finding : child.findings) {
            EXPECT_NE(finding.severity, Severity::ERROR) << finding.message;
        }
    }
}

TEST(ProteinGeometrySuite, EmptyMoleculeIsNotApplicableEverywhere) {
    const OEChem::OEGraphMol mol;
    const AssessmentResult result = make_protein_geometry_suite()->Run(mol);
    EXPECT_EQ(result.status, Status::NOT_APPLICABLE);
    ASSERT_EQ(result.children.size(), 5U);
    for (const AssessmentResult& child : result.children) {
        EXPECT_EQ(child.status, Status::NOT_APPLICABLE) << child.assessment;
    }
}

} // namespace OEHEALTH::test
```

- [ ] **Step 2: Register the test and confirm it fails**

Add to `add_executable(oehealth_tests ...)` in `tests/cpp/CMakeLists.txt`, keeping the list sorted:

```cmake
    test_protein_suite.cpp
```

Run: `cmake --build build-debug -j8`
Expected: FAIL to compile with `use of undeclared identifier 'make_protein_geometry_suite'` (the template's `oehealth.h` does not include the framework headers yet).

- [ ] **Step 3: Declare and define the suite factory**

```diff
diff --git a/include/oehealth/suite.h b/include/oehealth/suite.h
--- a/include/oehealth/suite.h
+++ b/include/oehealth/suite.h
@@ -36,6 +36,10 @@ private:
     std::vector<std::shared_ptr<Assessment>> children_;
 };
 
+/// \brief The protein geometry suite: bond lengths, bond angles, clashes,
+/// valence, and Ramachandran, each with default options.
+std::shared_ptr<AssessmentSuite> make_protein_geometry_suite();
+
 } // namespace OEHEALTH
 
 #endif // OEHEALTH_SUITE_H
```

```diff
diff --git a/src/suite.cpp b/src/suite.cpp
--- a/src/suite.cpp
+++ b/src/suite.cpp
@@ -1,5 +1,11 @@
 #include "oehealth/suite.h"
 
+#include "oehealth/assessments/bond_angles.h"
+#include "oehealth/assessments/bond_lengths.h"
+#include "oehealth/assessments/clashes.h"
+#include "oehealth/assessments/ramachandran.h"
+#include "oehealth/assessments/valence.h"
+
 #include <exception>
 #include <stdexcept>
 #include <utility>
@@ -52,4 +58,16 @@ AssessmentResult AssessmentSuite::Assess(const StructureContext& ctx) const {
     return result;
 }
 
+std::shared_ptr<AssessmentSuite> make_protein_geometry_suite() {
+    std::vector<std::shared_ptr<Assessment>> children = {
+        std::make_shared<BondLengthAssessment>(),   std::make_shared<BondAngleAssessment>(),
+        std::make_shared<ClashAssessment>(),        std::make_shared<ValenceAssessment>(),
+        std::make_shared<RamachandranAssessment>(),
+    };
+    return std::make_shared<AssessmentSuite>("protein_geometry",
+                                             "Bond lengths, bond angles, clashes, valence, and "
+                                             "Ramachandran checks for protein structures.",
+                                             std::move(children));
+}
+
 } // namespace OEHEALTH
```

- [ ] **Step 4: Include every public header from `oehealth.h`**

Insert the includes above `#include <oechem.h>`. Leave the version macros (vrzn manages them) and the `calculate_molecular_weight` declaration as they are.

```diff
diff --git a/include/oehealth/oehealth.h b/include/oehealth/oehealth.h
--- a/include/oehealth/oehealth.h
+++ b/include/oehealth/oehealth.h
@@ -6,6 +6,17 @@
 #define OEHEALTH_VERSION_MINOR 1
 #define OEHEALTH_VERSION_PATCH 0
 
+#include "oehealth/assessment.h"
+#include "oehealth/assessments/bond_angles.h"
+#include "oehealth/assessments/bond_lengths.h"
+#include "oehealth/assessments/clashes.h"
+#include "oehealth/assessments/ramachandran.h"
+#include "oehealth/assessments/valence.h"
+#include "oehealth/context.h"
+#include "oehealth/json.h"
+#include "oehealth/result.h"
+#include "oehealth/suite.h"
+
 #include <oechem.h>
 
 namespace OEHEALTH {
```

- [ ] **Step 5: Build and run the tests**

```bash
cmake --build build-debug -j8
./build-debug/tests/cpp/oehealth_tests --gtest_filter='Version.*:ProteinGeometrySuite.*'
ctest --test-dir build-debug --output-on-failure
PYTHONPATH=python python -m pytest tests/python -q
```

Expected: `[  PASSED  ] 4 tests.` from the filtered run; `100% tests passed out of 83`; pytest still reports `12 passed`, which shows the SWIG module builds against the new header.

- [ ] **Step 6: Format, lint, and check versions**

```bash
clang-format --dry-run --Werror include/oehealth/suite.h src/suite.cpp include/oehealth/oehealth.h tests/cpp/test_protein_suite.cpp
clang-tidy -p build-debug --quiet --extra-arg=-isysroot --extra-arg="$(xcrun --show-sdk-path)" src/suite.cpp tests/cpp/test_protein_suite.cpp
vrzn get
```

Expected: no output from clang-format; no findings from clang-tidy; vrzn ends with `All version numbers are consistent.`

- [ ] **Step 7: Commit**

```bash
git add include/oehealth/suite.h src/suite.cpp include/oehealth/oehealth.h tests/cpp/test_protein_suite.cpp tests/cpp/CMakeLists.txt
git commit -m "Add the protein geometry suite and include all public headers from oehealth.h"
```

---

### Task 12: Python bindings, template removal, and README

**Risk tier:** high — deletes two tracked files with `git rm` and removes the public `calculate_molecular_weight` API from C++ and Python, alongside the SWIG interface and package changes.

**Depends on:** Task 11

**Files:**
- Modify: `swig/oehealth.i:6-10`, `swig/oehealth.i:355-386`
- Modify: `python/oehealth/__init__.py:559-567`
- Modify: `include/oehealth/oehealth.h`
- Modify: `CMakeLists.txt` (`OEHEALTH_SOURCES`), `tests/cpp/CMakeLists.txt`
- Modify: `README.md`
- Delete: `src/oehealth.cpp`, `tests/cpp/test_oehealth.cpp`
- Test: `tests/python/conftest.py`, `tests/python/test_oehealth.py` (both rewritten)

**Interfaces:**
- Consumes: every public header from Tasks 2-11.
- Produces (Python, `oehealth` package): `Assessment`, `AssessmentResult` (with `to_dict()`), `AssessmentSuite`, `AtomRef`, `Finding`, `Metric`, `BondLengthAssessment`, `BondLengthOptions`, `BondAngleAssessment`, `BondAngleOptions`, `ClashAssessment`, `ClashOptions`, `ValenceAssessment`, `RamachandranAssessment`, `make_protein_geometry_suite`, `to_json`, `worst_status`, `THRESHOLD_TOLERANCE`, `JSON_SCHEMA_VERSION`, and the enum values `Severity_INFO`, `Severity_WARNING`, `Severity_ERROR`, `Status_NOT_APPLICABLE`, `Status_PASS`, `Status_WARN`, `Status_FAIL`, `Status_ERROR`. `std::invalid_argument` becomes `ValueError`; a `None` suite child raises `ValueError` before reaching C++.

How the SWIG changes work: `%feature("naturalvar")` makes struct fields return copies, and the `OEHEALTH_RETURN_COPY` typemaps make `GetOptions()` return a copy, so a value read from a temporary result or assessment stays valid after Python collects the temporary. `test_finding_fields_survive_collection` and the options tests check this.

- [ ] **Step 1: Write the failing Python tests**

`tests/python/conftest.py` (replaces the template's aspirin and ethanol fixtures):

```python
"""Shared fixtures and configuration for oehealth Python tests."""

from pathlib import Path

import pytest

pytest.importorskip("openeye.oechem", reason="OpenEye Toolkits not installed")

DATA_DIR = Path(__file__).resolve().parents[1] / "data"


@pytest.fixture
def crambin():
    """Read PDB entry 1CRN (crambin, 46 residues, no hydrogens)."""
    from openeye import oechem

    mol = oechem.OEGraphMol()
    stream = oechem.oemolistream(str(DATA_DIR / "1crn.pdb"))
    assert oechem.OEReadMolecule(stream, mol)
    return mol


@pytest.fixture
def ethanol_2d():
    """Create ethanol from SMILES, which has atoms but no 3D coordinates."""
    from openeye import oechem

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "CCO")
    return mol
```

`tests/python/test_oehealth.py` (replaces the template's tests; the last two lines of `test_version_available` are the version assertions that vrzn manages, unchanged from the template):

```python
"""Tests for the oehealth Python bindings."""

import gc

import pytest

import oehealth

ASSESSMENTS = [
    (oehealth.BondLengthAssessment, "bond_lengths"),
    (oehealth.BondAngleAssessment, "bond_angles"),
    (oehealth.ClashAssessment, "clashes"),
    (oehealth.ValenceAssessment, "valence"),
    (oehealth.RamachandranAssessment, "ramachandran"),
]


class TestAssessments:
    """Each assessment runs on a native openeye.oechem molecule."""

    @pytest.mark.parametrize(("factory", "name"), ASSESSMENTS)
    def test_run_on_crambin(self, crambin, factory, name):
        result = factory().Run(crambin)
        assert result.assessment == name
        assert result.status in (oehealth.Status_PASS, oehealth.Status_WARN)
        assert len(result.metrics) > 0
        for finding in result.findings:
            assert finding.severity != oehealth.Severity_ERROR

    def test_ramachandran_metrics(self, crambin):
        result = oehealth.RamachandranAssessment().Run(crambin)
        metrics = {metric.name: metric.value for metric in result.metrics}
        assert metrics["n_assessed"] == 44
        assert metrics["n_outliers"] == 0

    def test_finding_fields_survive_collection(self, crambin):
        findings = oehealth.BondLengthAssessment().Run(crambin).findings
        gc.collect()
        assert len(findings) > 0
        atom = findings[0].atoms[0]
        assert atom.chain_id == "A"
        assert atom.residue_number > 0
        assert atom.atom_name != ""

    def test_rejects_2d_molecule(self, ethanol_2d):
        with pytest.raises(ValueError, match="3D"):
            oehealth.BondLengthAssessment().Run(ethanol_2d)

    def test_rejects_non_molecule(self):
        with pytest.raises(TypeError):
            oehealth.BondLengthAssessment().Run("not a molecule")


class TestOptions:
    """Options structs are editable from Python and validated in C++."""

    def test_defaults(self):
        options = oehealth.ClashAssessment().GetOptions()
        assert options.warning_overlap == pytest.approx(0.4)
        assert options.error_overlap == pytest.approx(1.0)
        assert options.polar_contact_distance == pytest.approx(2.9)

    def test_set_fields(self):
        options = oehealth.BondLengthOptions()
        options.warning_z = 3.0
        options.error_z = 6.0
        returned = oehealth.BondLengthAssessment(options).GetOptions()
        gc.collect()
        assert returned.warning_z == pytest.approx(3.0)
        assert returned.error_z == pytest.approx(6.0)

    def test_invalid_options_raise_value_error(self):
        options = oehealth.ClashOptions()
        options.polar_contact_distance = 0.0
        with pytest.raises(ValueError):
            oehealth.ClashAssessment(options)
        reversed_thresholds = oehealth.BondAngleOptions()
        reversed_thresholds.warning_z = 5.0
        reversed_thresholds.error_z = 4.0
        with pytest.raises(ValueError):
            oehealth.BondAngleAssessment(reversed_thresholds)


class TestSuites:
    """Suites compose assessments in C++ and from Python."""

    def test_protein_geometry_suite(self, crambin):
        suite = oehealth.make_protein_geometry_suite()
        assert suite.GetName() == "protein_geometry"
        assert [child.GetName() for child in suite.GetChildren()] == [
            name for _, name in ASSESSMENTS
        ]
        result = suite.Run(crambin)
        assert [child.assessment for child in result.children] == [name for _, name in ASSESSMENTS]
        assert result.status != oehealth.Status_ERROR

    def test_custom_suite_from_python(self, crambin):
        suite = oehealth.AssessmentSuite(
            "custom",
            "Clashes and valence.",
            [oehealth.ClashAssessment(), oehealth.ValenceAssessment()],
        )
        result = suite.Run(crambin)
        assert result.assessment == "custom"
        assert [child.assessment for child in result.children] == ["clashes", "valence"]

    def test_empty_suite_name_raises_value_error(self):
        with pytest.raises(ValueError):
            oehealth.AssessmentSuite("", "", [])

    def test_none_child_raises_value_error(self):
        with pytest.raises(ValueError):
            oehealth.AssessmentSuite("s", "", [oehealth.ClashAssessment(), None])


class TestJson:
    """Results serialize to JSON and to dictionaries."""

    def test_to_dict(self, crambin):
        data = oehealth.BondLengthAssessment().Run(crambin).to_dict()
        assert data["schema_version"] == oehealth.JSON_SCHEMA_VERSION
        assert data["oehealth_version"] == oehealth.__version__
        assert data["assessment"] == "bond_lengths"
        assert list(data["metrics"]) == ["n_assessed", "n_outliers", "rmsz", "max_abs_z"]

    def test_to_json_compact(self, crambin):
        text = oehealth.to_json(oehealth.ValenceAssessment().Run(crambin), -1)
        assert "\n" not in text

    def test_suite_children_omit_version_keys(self, crambin):
        data = oehealth.make_protein_geometry_suite().Run(crambin).to_dict()
        assert len(data["children"]) == 5
        assert "schema_version" not in data["children"][0]


def test_version_available():
    """Verify version info is accessible."""
    assert hasattr(oehealth, "__version__")
    assert hasattr(oehealth, "__version_info__")
    assert oehealth.__version__ == "0.1.0"
    assert oehealth.__version_info__ == (0, 1, 0)
```

- [ ] **Step 2: Run the Python tests to confirm they fail**

Run: `PYTHONPATH=python python -m pytest tests/python/test_oehealth.py -q`
Expected: FAIL during collection with `AttributeError: module 'oehealth' has no attribute 'BondLengthAssessment'`.

- [ ] **Step 3: Wrap the framework in SWIG**

Apply with `git apply`. This adds the `<stdexcept>` include, the vector and shared_ptr typemaps, the `ValueError` catch, and replaces the `calculate_molecular_weight` declaration with the framework API.

```diff
diff --git a/swig/oehealth.i b/swig/oehealth.i
--- a/swig/oehealth.i
+++ b/swig/oehealth.i
@@ -7,6 +7,8 @@
 #include <oechem.h>
 #include <oegrid.h>
 
+#include <stdexcept>
+
 using namespace OEHEALTH;
 %}
 
@@ -355,13 +357,19 @@ OE_CROSS_RUNTIME_REF_TYPEMAPS(OEDocking::OEReceptor, _oehealth_is_oereceptor, "E
 %include "std_string.i"
 %include "stdint.i"
 %include "exception.i"
+%include "std_vector.i"
+%include "std_shared_ptr.i"
 
 // ============================================================================
 // Exception handling
 // ============================================================================
+// std::invalid_argument (bad options, 2D input, null suite child) becomes
+// ValueError; every other C++ exception becomes RuntimeError.
 %exception {
     try {
         $action
+    } catch (const std::invalid_argument& e) {
+        SWIG_exception(SWIG_ValueError, e.what());
     } catch (const std::exception& e) {
         SWIG_exception(SWIG_RuntimeError, e.what());
     } catch (...) {
@@ -379,11 +387,79 @@ OE_CROSS_RUNTIME_REF_TYPEMAPS(OEDocking::OEReceptor, _oehealth_is_oereceptor, "E
 // ============================================================================
 // Wrapped API
 // ============================================================================
+// Assessments are held by shared_ptr so Python and C++ suites can share them.
+%shared_ptr(OEHEALTH::Assessment)
+%shared_ptr(OEHEALTH::AssessmentSuite)
+%shared_ptr(OEHEALTH::BondLengthAssessment)
+%shared_ptr(OEHEALTH::BondAngleAssessment)
+%shared_ptr(OEHEALTH::ClashAssessment)
+%shared_ptr(OEHEALTH::ValenceAssessment)
+%shared_ptr(OEHEALTH::RamachandranAssessment)
+
+// The context borrows a molecule that Python could collect while the context
+// lives, so Python calls Run(mol) instead of Assess(ctx).
+%ignore OEHEALTH::StructureContext;
+%ignore OEHEALTH::ResidueInfo;
+%ignore OEHEALTH::Assessment::Assess;
+
+// Declare the element types before the vector templates, so the templates
+// exist when the result structs' vector members are wrapped.
 namespace OEHEALTH {
+struct AtomRef;
+struct Metric;
+struct Finding;
+struct AssessmentResult;
+class Assessment;
+}
+
+// Member variables are read by value, so a field such as result.metrics stays
+// valid after the result it came from is collected.
+%feature("naturalvar", "1");
+
+// Return options by copy for the same reason: a reference into a temporary
+// assessment would dangle once Python collects it. (GetChildren already
+// returns a tuple, which std_vector.i builds from a copy.)
+%define OEHEALTH_RETURN_COPY(TYPE)
+%typemap(out) const TYPE& {
+    $result = SWIG_NewPointerObj(new TYPE(*$1), $descriptor(TYPE*), SWIG_POINTER_OWN);
+}
+%enddef
+OEHEALTH_RETURN_COPY(OEHEALTH::BondLengthOptions)
+OEHEALTH_RETURN_COPY(OEHEALTH::BondAngleOptions)
+OEHEALTH_RETURN_COPY(OEHEALTH::ClashOptions)
+
+// SWIG cannot convert None to an assessment, so check here and raise the
+// ValueError that C++ raises for a null child.
+%pythonprepend OEHEALTH::AssessmentSuite::AssessmentSuite %{
+    if any(child is None for child in children):
+        raise ValueError("an assessment suite child is None")
+%}
 
-double calculate_molecular_weight(const OEChem::OEMolBase& mol);
+%template(AtomRefVector) std::vector<OEHEALTH::AtomRef>;
+%template(MetricVector) std::vector<OEHEALTH::Metric>;
+%template(FindingVector) std::vector<OEHEALTH::Finding>;
+%template(AssessmentResultVector) std::vector<OEHEALTH::AssessmentResult>;
+%template(AssessmentVector) std::vector<std::shared_ptr<OEHEALTH::Assessment>>;
+
+%include "oehealth/result.h"
+%include "oehealth/assessment.h"
+%include "oehealth/suite.h"
+%include "oehealth/assessments/bond_lengths.h"
+%include "oehealth/assessments/bond_angles.h"
+%include "oehealth/assessments/clashes.h"
+%include "oehealth/assessments/valence.h"
+%include "oehealth/assessments/ramachandran.h"
+%include "oehealth/json.h"
+
+%extend OEHEALTH::AssessmentResult {
+%pythoncode %{
+    def to_dict(self):
+        """Return this result as a dictionary, parsed from ``to_json(self)``."""
+        import json
 
-} // namespace OEHEALTH
+        return json.loads(to_json(self))
+%}
+}
 
 // ============================================================================
 // Module-level Python convenience code
```

- [ ] **Step 4: Re-export the API from the Python package**

```diff
diff --git a/python/oehealth/__init__.py b/python/oehealth/__init__.py
--- a/python/oehealth/__init__.py
+++ b/python/oehealth/__init__.py
@@ -557,11 +557,63 @@ _load_cached_extension_if_needed()
 _check_openeye_version()
 
 from .oehealth import (
-    calculate_molecular_weight,
+    JSON_SCHEMA_VERSION,
+    THRESHOLD_TOLERANCE,
+    Assessment,
+    AssessmentResult,
+    AssessmentSuite,
+    AtomRef,
+    BondAngleAssessment,
+    BondAngleOptions,
+    BondLengthAssessment,
+    BondLengthOptions,
+    ClashAssessment,
+    ClashOptions,
+    Finding,
+    Metric,
+    RamachandranAssessment,
+    Severity_ERROR,
+    Severity_INFO,
+    Severity_WARNING,
+    Status_ERROR,
+    Status_FAIL,
+    Status_NOT_APPLICABLE,
+    Status_PASS,
+    Status_WARN,
+    ValenceAssessment,
+    make_protein_geometry_suite,
+    to_json,
+    worst_status,
 )
 
 __all__ = [
+    "JSON_SCHEMA_VERSION",
+    "THRESHOLD_TOLERANCE",
+    "Assessment",
+    "AssessmentResult",
+    "AssessmentSuite",
+    "AtomRef",
+    "BondAngleAssessment",
+    "BondAngleOptions",
+    "BondLengthAssessment",
+    "BondLengthOptions",
+    "ClashAssessment",
+    "ClashOptions",
+    "Finding",
+    "Metric",
+    "RamachandranAssessment",
+    "Severity_ERROR",
+    "Severity_INFO",
+    "Severity_WARNING",
+    "Status_ERROR",
+    "Status_FAIL",
+    "Status_NOT_APPLICABLE",
+    "Status_PASS",
+    "Status_WARN",
+    "ValenceAssessment",
     "__version__",
     "__version_info__",
-    "calculate_molecular_weight",
+    "make_protein_geometry_suite",
+    "to_json",
+    "worst_status",
 ]
```

- [ ] **Step 5: Remove the template function**

```bash
git rm src/oehealth.cpp tests/cpp/test_oehealth.cpp
```

Replace `include/oehealth/oehealth.h` with this content (the version macros are unchanged):

```cpp
#ifndef OEHEALTH_H
#define OEHEALTH_H

// Version information
#define OEHEALTH_VERSION_MAJOR 0
#define OEHEALTH_VERSION_MINOR 1
#define OEHEALTH_VERSION_PATCH 0

#include "oehealth/assessment.h"
#include "oehealth/assessments/bond_angles.h"
#include "oehealth/assessments/bond_lengths.h"
#include "oehealth/assessments/clashes.h"
#include "oehealth/assessments/ramachandran.h"
#include "oehealth/assessments/valence.h"
#include "oehealth/context.h"
#include "oehealth/json.h"
#include "oehealth/result.h"
#include "oehealth/suite.h"

#endif // OEHEALTH_H
```

Remove `src/oehealth.cpp` from `set(OEHEALTH_SOURCES ...)` in `CMakeLists.txt` and `test_oehealth.cpp` from `add_executable(oehealth_tests ...)` in `tests/cpp/CMakeLists.txt`. The lists then read:

```cmake
set(OEHEALTH_SOURCES
    src/assessment.cpp
    src/assessments/bond_angles.cpp
    src/assessments/bond_lengths.cpp
    src/assessments/clashes.cpp
    src/assessments/geometry_outliers.cpp
    src/assessments/ramachandran.cpp
    src/assessments/valence.cpp
    src/common.cpp
    src/context.cpp
    src/data/engh_huber_2001.cpp
    src/geometry/measure.cpp
    src/geometry/targets.cpp
    src/json.cpp
    src/result.cpp
    src/suite.cpp
)
```

```cmake
add_executable(oehealth_tests
    test_assessments_geometry.cpp
    test_clashes.cpp
    test_context.cpp
    test_json.cpp
    test_protein_suite.cpp
    test_ramachandran.cpp
    test_reference_data.cpp
    test_result.cpp
    test_suite.cpp
    test_targets.cpp
    test_valence.cpp
)
```

- [ ] **Step 6: Build and run every test**

```bash
cmake --build build-debug -j8
ctest --test-dir build-debug --output-on-failure
PYTHONPATH=python python -m pytest tests/python -q
```

Expected: `100% tests passed out of 79`; pytest reports `27 passed`.

- [ ] **Step 7: Update the README**

This replaces the `calculate_molecular_weight` usage, project layout, and "Adding New Functions" sections, and adds the lint commands and the new tools.

````diff
diff --git a/README.md b/README.md
--- a/README.md
+++ b/README.md
@@ -6,6 +6,10 @@ This project provides a C++ library with Python bindings built using SWIG and
 [OpenEye Toolkits](https://www.eyesopen.com/). Molecules created with
 `openeye.oechem` in Python pass natively to C++ without serialization.
 
+The first assessments check protein structures: bond lengths and bond angles
+against Engh & Huber (2001) targets, atom clashes, valence, and Ramachandran
+outliers. Assessments combine into suites, and every result serializes to JSON.
+
 ## Prerequisites
 
 - **OpenEye C++ SDK** -- Headers and libraries (download from
@@ -77,22 +81,57 @@ Python tests:
 pytest tests/python/ -v
 ```
 
-The Python tests verify that molecules created with `openeye.oechem.OEGraphMol`
-pass correctly to the C++ `calculate_molecular_weight` function.
+Without an editable install, put the package on the path instead:
+`PYTHONPATH=python pytest tests/python/ -v`. The Python tests run each
+assessment on crambin (PDB entry 1CRN, `tests/data/1crn.pdb`) through the
+bindings, and test the reference-table generator.
+
+Lint and format checks (clang-tidy needs the build's `compile_commands.json`;
+on macOS also pass `--extra-arg=-isysroot --extra-arg="$(xcrun --show-sdk-path)"`):
+
+```bash
+clang-format --dry-run --Werror include/oehealth/*.h src/*.cpp
+clang-tidy -p build-debug src/*.cpp
+ruff check .
+mypy
+```
 
 ## Usage
 
 ```python
 from openeye import oechem
-from oehealth import calculate_molecular_weight
+import oehealth
 
 mol = oechem.OEGraphMol()
-oechem.OESmilesToMol(mol, "CC(=O)OC1=CC=CC=C1C(=O)O")  # aspirin
+oechem.OEReadMolecule(oechem.oemolistream("tests/data/1crn.pdb"), mol)
 
-mw = calculate_molecular_weight(mol)
-print(f"Molecular weight: {mw:.2f}")  # 180.16
+result = oehealth.make_protein_geometry_suite().Run(mol)
+for child in result.children:
+    print(child.assessment, [(m.name, m.value) for m in child.metrics])
+    for finding in child.findings:
+        print("   ", finding.message)
+
+report = result.to_dict()  # the same content as oehealth.to_json(result)
 ```
 
+The same in C++:
+
+```cpp
+#include "oehealth/oehealth.h"
+
+OEChem::OEGraphMol mol;
+OEChem::oemolistream stream("tests/data/1crn.pdb");
+OEChem::OEReadMolecule(stream, mol);
+const OEHEALTH::AssessmentResult result = OEHEALTH::make_protein_geometry_suite()->Run(mol);
+std::cout << OEHEALTH::to_json(result) << std::endl;
+```
+
+Each assessment can also run on its own, for example
+`oehealth.ClashAssessment().Run(mol)`, and takes an options struct where it has
+thresholds (`BondLengthOptions`, `BondAngleOptions`, `ClashOptions`).
+Molecules are assessed as given: nothing is added or changed. The valence check
+needs implicit hydrogen counts, which OpenEye's readers assign by default.
+
 The `OEGraphMol` object passes directly from Python to C++ via cross-runtime
 SWIG typemaps. No SMILES round-trip or manual pointer handling is needed.
 
@@ -105,9 +144,14 @@ oehealth/
     pyproject.toml                  # Package metadata + [tool.oe-build] config
     vrzn.toml                       # Version locations for vrzn
     include/oehealth/
-        oehealth.h               # Public C++ header
+        oehealth.h                  # Version macros; includes every public header
+        result.h, context.h, assessment.h, suite.h, json.h
+        assessments/                # One header per assessment
     src/
-        oehealth.cpp             # C++ implementation
+        *.cpp                       # Framework implementation
+        assessments/                # Assessment implementations
+        geometry/                   # Measurements and Engh & Huber target selection
+        data/engh_huber_2001.cpp    # Generated reference tables; do not edit
     swig/
         oehealth.i               # SWIG interface with OEMolBase typemaps
         CMakeLists.txt              # SWIG module build rules
@@ -117,72 +161,52 @@ oehealth/
             __init__.py             # Python package (imports, compat layer)
     scripts/
         build_python.py             # Build distributable wheels
+        generate_engh_huber.py      # Regenerate src/data/engh_huber_2001.cpp
     tests/
+        data/                       # 1crn.pdb and a standard_geometry.cif excerpt
         cpp/
             CMakeLists.txt          # C++ test build rules
-            test_oehealth.cpp     # C++ unit tests
+            test_*.cpp              # C++ unit tests
         python/
-            conftest.py             # Pytest fixtures (molecule helpers)
-            test_oehealth.py      # Python tests
+            conftest.py             # Pytest fixtures (crambin)
+            test_*.py               # Python tests
     .github/workflows/
         build-wheels.yml            # CI: multi-platform wheel builds
 ```
 
-## Adding New Functions
-
-The scaffolded `calculate_molecular_weight` function demonstrates the pattern for
-wrapping C++ code that operates on OpenEye molecules. To add your own functions:
+## Adding New Assessments
 
-**1. Declare in the header** (`include/oehealth/oehealth.h`):
+**1. Declare** a class derived from `OEHEALTH::Assessment` in a new header
+under `include/oehealth/assessments/`, and include it from
+`include/oehealth/oehealth.h`. Override `GetName()` (a stable id),
+`GetDescription()`, and `Assess(const StructureContext&)`.
 
-```cpp
-namespace OEHEALTH {
-
-double my_function(const OEChem::OEMolBase& mol);
+**2. Implement** it in `src/assessments/`, and add the file to
+`OEHEALTH_SOURCES` in `CMakeLists.txt`. Build results with the helpers in
+`src/common.h` (`make_atom_ref`, `atom_path`, `status_from_findings`).
 
-} // namespace OEHEALTH
-```
+**3. Expose in SWIG** (`swig/oehealth.i`): add a `%shared_ptr(...)` line for the
+class next to the others, and `%include` its header.
 
-**2. Implement** (`src/oehealth.cpp` or a new `.cpp` file):
+**4. Re-export in Python** (`python/oehealth/__init__.py`): add the class to the
+import list and `__all__`.
 
-```cpp
-namespace OEHEALTH {
-
-double my_function(const OEChem::OEMolBase& mol) {
-    // Your implementation using OpenEye C++ API
-}
-
-} // namespace OEHEALTH
-```
-
-If you add new `.cpp` files, add them to the `OEHEALTH_SOURCES` list in
-`CMakeLists.txt`.
-
-**3. Expose in SWIG** (`swig/oehealth.i`), under the "Wrapped API" section:
+**5. Rebuild and test:**
 
-```swig
-namespace OEHEALTH {
-double my_function(const OEChem::OEMolBase& mol);
-}
+```bash
+cmake --build build-debug
+ctest --test-dir build-debug --output-on-failure
+pytest tests/python/ -v
 ```
 
-Any function accepting `OEMolBase&` or `const OEMolBase&` will automatically use
-the cross-runtime typemaps.
-
-**4. Import in Python** (`python/oehealth/__init__.py`):
-
-```python
-from .oehealth import (
-    calculate_molecular_weight,
-    my_function,
-)
-```
+## Reference Data
 
-**5. Rebuild and test:**
+`src/data/engh_huber_2001.cpp` is generated from `standard_geometry.cif`
+(NDB/wwPDB, distributed with CCP4 as `lib/data/standard_geometry.cif`). To
+regenerate it:
 
 ```bash
-cmake --build build-debug
-pytest tests/python/ -v
+python3 scripts/generate_engh_huber.py /path/to/standard_geometry.cif src/data/engh_huber_2001.cpp
 ```
 
 ## Building Wheels
@@ -270,6 +294,10 @@ This project uses several tools to manage the build, packaging, and versioning:
 | [delocate](https://github.com/matthew-brett/delocate) | Bundles shared libraries into macOS wheels |
 | [auditwheel](https://github.com/pypa/auditwheel) | Bundles shared libraries into Linux wheels |
 | [pytest](https://docs.pytest.org/) | Test framework for the Python test suite |
+| [GoogleTest](https://github.com/google/googletest) | Test framework for the C++ test suite |
+| [nlohmann/json](https://github.com/nlohmann/json) | JSON serialization of results (fetched by CMake) |
+| [clang-format](https://clang.llvm.org/docs/ClangFormat.html), [clang-tidy](https://clang.llvm.org/extra/clang-tidy/) | C++ formatting and linting (`.clang-format`, `.clang-tidy`) |
+| [ruff](https://docs.astral.sh/ruff/), [mypy](https://mypy-lang.org/) | Python linting and type checking (`pyproject.toml`) |
 
 ## Version Management
 
````

Then confirm no stale references remain:

```bash
grep -rn "calculate_molecular_weight" README.md include src swig python/oehealth/__init__.py tests
```

Expected: no output (grep exits 1).

- [ ] **Step 8: Run every lint and version check**

```bash
ruff check .
ruff format --check scripts tests/python
mypy --python-executable "$(command -v python)"
find include src tests/cpp -name '*.h' -o -name '*.cpp' | sort | xargs clang-format --dry-run --Werror
find src tests/cpp -name '*.cpp' | sort | xargs clang-tidy -p build-debug --quiet --extra-arg=-isysroot --extra-arg="$(xcrun --show-sdk-path)"
vrzn get
```

Expected: `All checks passed!`; `4 files already formatted`; `Success: no issues found in 5 source files`; no clang-format output; clang-tidy prints only its progress lines; vrzn ends with `All version numbers are consistent.`

- [ ] **Step 9: Commit**

```bash
git add swig/oehealth.i python/oehealth/__init__.py include/oehealth/oehealth.h CMakeLists.txt tests/cpp/CMakeLists.txt README.md tests/python/conftest.py tests/python/test_oehealth.py
git commit -m "Expose the assessment framework to Python and remove calculate_molecular_weight"
```

`git rm` in Step 5 already staged the two deletions.
