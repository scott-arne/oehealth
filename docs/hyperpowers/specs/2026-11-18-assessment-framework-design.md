# Structure Health Assessment Framework: Design

Date: 2026-11-18
Status: approved design

## Goal

oehealth will hold many ways to assess the health of a molecular structure:
single checks (bond lengths, Ramachandran), named sets of checks (protein
structure health), more involved checks (ligand geometry), and reports. This
first round lays the groundwork and implements five simple checks:

- bond lengths
- bond angles
- clashes
- valence issues
- Ramachandran violations

## Scope

In scope:

- A C++ assessment interface, a shared structure context, a common result
  model, and composable suites.
- The five assessments above, plus one example suite, `protein_geometry`.
- JSON serialization of results in C++.
- SWIG bindings and Python re-exports for all of the above.
- A script that generates the Engh & Huber reference tables, and its
  generated output.
- Tooling: clang-format, clang-tidy, and ruff/mypy configuration.
- Removal of the template's `calculate_molecular_weight` function and its
  tests.

Out of scope for this round:

- Reports of any kind (text, HTML).
- Reading JSON back into result objects.
- Ideal geometry for anything other than the 20 standard amino acids
  (nucleic acids, ligands, modified residues such as MSE, ASX/GLX/UNK).
- Adding hydrogens or otherwise preparing structures.
- Metal coordination geometry.
- Python-defined assessments (SWIG directors).
- An assessment registry or lookup by name.
- Missing-atom and chain-break reporting.

## Decisions

- **C++ first.** Assessments, result types, and suites live in C++ and are
  exposed to Python through SWIG. Rejected: Python-side suites and reports
  over C++ assessments, because C++ callers must be able to run whole suites.
- **Polymorphic assessment objects with a shared context.** An abstract
  `Assessment` class whose `Assess` takes a `StructureContext`, one common
  `AssessmentResult` type, and suites that are themselves assessments.
  Rejected: an enum plus central dispatcher (every new check would edit the
  enum, dispatcher, serializer, and bindings, which scales poorly for a
  growing catalog); a pipeline with finding sinks and injected reference
  providers (lifetimes, dependency ordering, and provider interfaces that
  five simple checks do not need); a separate result type per assessment
  (each needs its own SWIG wrapping and JSON serializer, and suites would
  need type erasure anyway). The context is included now because changing a
  virtual signature later means editing every assessment.
- **The context builds its own residue table instead of holding an
  `OEHierView`.** `OEBio::OEHierView` can only be constructed from a
  non-const `OEMolBase&`, and its iteration methods are non-const, so it
  cannot be built from the borrowed `const OEMolBase&`. Rejected: copying
  the molecule to build a view (doubles memory and complicates atom index
  mapping); `const_cast` (breaks the promise not to modify the input).
- **Ideal geometry covers standard amino acids only, from Engh & Huber
  (2001).** This matches the reference used by MolProbity and earlier wwPDB
  validation. Rejected: adding nucleic acids now (more data curation, no
  design change); MMFF ideal values for all molecules (not crystallographic
  targets, and overlaps the planned ligand-geometry work).
- **Reference tables are generated from `standard_geometry.cif`, not
  transcribed by hand.** The file (NDB/wwPDB, version "Fri Feb 22 17:25:15
  GMT 2013", distributed with CCP4 9 as `lib/data/standard_geometry.cif`)
  states that its amino-acid values come from Engh & Huber (2001), and spot
  checks against the published values match. A committed Python script
  converts it to C++; the generated file is committed so builds and tests
  never need CCP4. Rejected: manual transcription (error-prone and
  unverifiable); the current wwPDB validation targets (they no longer use
  Engh & Huber; for example N–CA is 1.447 ± 0.012 Å there versus 1.459 ±
  0.020 Å in Engh & Huber).
- **Clashes are assessed on the atoms as given.** The input is never
  modified, and the result records whether hydrogens were present.
  Rejected: adding hydrogens to a copy (slower, depends on protonation, and
  adds a failure mode to a simple check); requiring hydrogens (most raw PDB
  input would get no result).
- **Pairs involving metals are excluded from the clash check.** Bondi radii
  make normal metal coordination (Zn–N near 2.1 Å) look like a 0.8 Å
  overlap. Metal coordination belongs in its own future assessment.
- **JSON serialization ships in this round.** nlohmann/json is a private
  dependency used only in `src/json.cpp`. Rejected: deferring serialization
  (the schema would be retrofitted after reports exist); adding a text
  report now (its format would be designed before seeing real results).
- **Suites are immutable and take their children at construction.** Suites
  hold `std::shared_ptr<Assessment>`, which SWIG handles cleanly from
  Python. With no `Add` method, a suite cannot contain itself or form a
  cycle, and every assessment, suites included, is immutable after
  construction, so sharing children between suites and threads is safe.
  Rejected: a mutable `Add` (needs cycle detection and a freeze step);
  `std::unique_ptr` ownership (fragile to transfer from Python).
- **Findings refer to atoms by copied identifiers, not OpenEye pointers.**
  Results then outlive the molecule and serialize without lookups.
- **`StructureContext` is not exposed to Python.** It borrows the molecule,
  and Python could collect the molecule while the context lives. Python
  calls `Run(mol)` instead.
- **The valence check assumes implicit hydrogen counts are assigned.**
  `OECheckAtomValence` flags under-valent atoms as well as over-valent ones.
  Probe on PDB entry 1CRN read with OpenEye's default PDB reader (no
  explicit hydrogens, implicit counts assigned by the reader): 0 failures.
  With implicit counts zeroed: 204 failures. The check reports every
  failure, and the documentation states the requirement. Rejected:
  reporting only over-valent atoms (would hide real under-valence in
  prepared structures).
- **Histidine protonation is inferred from explicit hydrogens, and
  otherwise chosen by best fit.** Without hydrogens the protonation state is
  unknown, and the variants differ by more than one sigma for some bonds.
  Rejected: always using one variant (produces false outliers for the other
  states).
- **Tooling: clang-format, clang-tidy, ruff, and mypy are configured from
  the start.** Rejected: relying on tests alone.

## Architecture

All public code lives in namespace `OEHEALTH`. Names follow the project
conventions: classes and methods in PascalCase, free functions in
lower_snake_case, enum values and constants in CAPITAL_SNAKE_CASE, variables
in lower_snake_case. Doc comments are Doxygen-compatible.

### File layout

```
include/oehealth/
    oehealth.h                  # version macros; includes every public header
    result.h                    # Severity, Status, AtomRef, Metric, Finding, AssessmentResult
    context.h                   # ResidueInfo, StructureContext
    assessment.h                # Assessment
    suite.h                     # AssessmentSuite, make_protein_geometry_suite
    json.h                      # to_json
    assessments/
        bond_lengths.h
        bond_angles.h
        clashes.h
        valence.h
        ramachandran.h
src/
    result.cpp, context.cpp, assessment.cpp, suite.cpp, json.cpp
    assessments/*.cpp
    geometry/
        targets.h, targets.cpp  # private: variant selection and target lookup
    data/
        engh_huber_2001.h       # private; not installed
        engh_huber_2001.cpp     # generated; do not edit by hand
scripts/
    generate_engh_huber.py
tests/
    data/1crn.pdb
    cpp/test_*.cpp
    python/test_*.py
```

Every new source file is added to `OEHEALTH_SOURCES` in `CMakeLists.txt`.

### Build changes

- The library links `OpenEye::OEBio` in addition to `OpenEye::OEChem`
  (`OERamachandranAnalysis` lives in OEBio, which itself depends on OEGrid).
- `swig/CMakeLists.txt` adds `OEBIO_LIBRARY` and `OEGRID_LIBRARY` to
  `EXPECTED_LIB_VARS`, and `pyproject.toml` adds `liboebio` and `liboegrid`
  to `expected-missing-libs`.
- nlohmann/json is fetched with FetchContent from
  `https://github.com/nlohmann/json.git` at tag `v3.11.3`, and linked
  `PRIVATE`.

### Result model (`result.h`)

```cpp
enum class Severity { INFO, WARNING, ERROR };
enum class Status { NOT_APPLICABLE, PASS, WARN, FAIL, ERROR };

struct AtomRef {
    unsigned int atom_idx;
    std::string atom_name;
    std::string residue_name;
    std::string chain_id;
    int residue_number;
    std::string insert_code;    // empty when the residue has none
};

struct Metric {
    std::string name;
    double value;
};

struct Finding {
    std::string code;                 // stable identifier, e.g. "bond_length_outlier"
    Severity severity;
    std::string message;              // human-readable, includes residue and atom names
    std::vector<AtomRef> atoms;
    double observed;                  // NaN when not meaningful
    double expected;                  // NaN when not meaningful
    double deviation;                 // NaN when not meaningful
    std::vector<Metric> details;      // extra numbers, e.g. phi and psi
};

struct AssessmentResult {
    std::string assessment;           // GetName() of the producer
    Status status;
    std::string message;              // set for ERROR and NOT_APPLICABLE
    std::vector<Metric> metrics;      // ordered
    std::vector<Finding> findings;
    std::vector<AssessmentResult> children;   // filled only by suites
};

Status worst_status(Status a, Status b);   // by declaration order
```

Invariants:

- Metric names are unique within `AssessmentResult::metrics` and within each
  `Finding::details`.
- `insert_code` is a string (empty or one character) so SWIG and JSON handle
  it without a `char` special case.

Status rules:

- An assessment that found nothing to assess returns NOT_APPLICABLE with a
  message, no findings, and its metrics where they are still defined (for
  example `hydrogens_present`).
- Otherwise an assessment's status comes from its findings: any ERROR
  finding gives FAIL, otherwise any WARNING finding gives WARN, otherwise
  PASS. INFO findings never affect status.
- ERROR is never derived from findings. Only a suite sets it, for a child
  that threw (see Suites). Outside a suite, exceptions propagate to the
  caller.

### Structure context (`context.h`)

```cpp
struct ResidueInfo {
    std::string name;              // e.g. "ALA"
    std::string chain_id;
    int number;
    std::string insert_code;
    bool is_standard_amino_acid;   // one of the 20 standard residues
    std::vector<const OEChem::OEAtomBase*> atoms;   // in molecule order
    const OEChem::OEAtomBase* FindAtom(const std::string& atom_name) const;  // nullptr if absent
};

class StructureContext {
public:
    explicit StructureContext(const OEChem::OEMolBase& mol);
    const OEChem::OEMolBase& GetMolecule() const;
    const std::vector<ResidueInfo>& GetResidues() const;   // built on first call
    bool HasHydrogens() const;                              // explicit H atoms present
};
```

The context borrows the molecule; the caller keeps the molecule alive and
unchanged for the context's lifetime. Lazily computed members are `mutable`,
so a single context is not thread-safe; the documentation says so.

Residue table: atoms are grouped into residues by OEChem residue identity
(`OEAtomGetResidue` and `OESameResidue`). Residues are ordered by the index
of their first atom, which for PDB input follows file order, so insertion
codes and numbering gaps need no special handling. Atom names are compared
after trimming whitespace. "Standard amino acid" means the residue name is
one of ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR
TRP TYR VAL. Alternate locations: OpenEye's PDB reader keeps one location
per atom by default, and the context assesses whatever atoms are present.

Peptide neighbors: residue j is the peptide successor of residue i when j
is the next entry in the residue table, both are standard amino acids, they
have the same chain id, both C(i) and N(j) are present, and the C(i)–N(j)
distance is under 2.5 Å. At 2.5 Å or more it is treated as a chain break and
no peptide-link terms are assessed.

### Assessment interface (`assessment.h`)

```cpp
class Assessment {
public:
    virtual ~Assessment() = default;
    virtual std::string GetName() const = 0;          // stable id, e.g. "bond_lengths"
    virtual std::string GetDescription() const = 0;
    virtual AssessmentResult Assess(const StructureContext& ctx) const = 0;
    AssessmentResult Run(const OEChem::OEMolBase& mol) const;
};
```

`Run` validates the molecule, builds a context, and calls `Assess`.
Validation order: a molecule with no atoms is passed through (each
assessment then returns NOT_APPLICABLE); a molecule with atoms whose
`GetDimension()` is not 3 throws `std::invalid_argument`. `Run` assesses
the active conformer of a multi-conformer molecule, and its documentation
says so. `Run` has a different name from `Assess` because a subclass that
overrides `Assess` would hide an `Assess(const OEMolBase&)` overload.

Concrete assessments and their options (each in its own header under
`include/oehealth/assessments/`):

```cpp
struct BondLengthOptions { double warning_z = 4.0; double error_z = 10.0; };
struct BondAngleOptions  { double warning_z = 4.0; double error_z = 10.0; };
struct ClashOptions {
    double warning_overlap = 0.4;        // Å
    double error_overlap = 1.0;          // Å
    double polar_contact_distance = 2.9; // Å
};

class BondLengthAssessment : public Assessment {
public:
    explicit BondLengthAssessment(const BondLengthOptions& options = BondLengthOptions());
    const BondLengthOptions& GetOptions() const;
    // GetName() == "bond_lengths"
};
class BondAngleAssessment : public Assessment {
public:
    explicit BondAngleAssessment(const BondAngleOptions& options = BondAngleOptions());
    const BondAngleOptions& GetOptions() const;
    // GetName() == "bond_angles"
};
class ClashAssessment : public Assessment {
public:
    explicit ClashAssessment(const ClashOptions& options = ClashOptions());
    const ClashOptions& GetOptions() const;
    // GetName() == "clashes"
};
class ValenceAssessment : public Assessment {
public:
    ValenceAssessment();
    // GetName() == "valence"
};
class RamachandranAssessment : public Assessment {
public:
    RamachandranAssessment();
    // GetName() == "ramachandran"
};
```

Each class also overrides `GetDescription()` with a one-sentence
description and `Assess()`. Bond lengths and bond angles have separate
options structs, with identical fields for now, so either can gain fields
later without affecting the other.

Constructors validate options and throw `std::invalid_argument` when any
field is NaN, infinite, or negative, when `polar_contact_distance` is zero,
or when an ERROR threshold is less than its WARNING threshold. Equal
thresholds are allowed. When a value meets both thresholds, ERROR takes
precedence.

Threshold comparisons are inclusive and use a tolerance: a value meets a
threshold t when value ≥ t − `THRESHOLD_TOLERANCE`, with
`constexpr double THRESHOLD_TOLERANCE = 1e-6` in `assessment.h`. Without
it, floating-point rounding would decide boundary cases (1.7 + 1.7 − 3.0
evaluates to 0.39999…, just under a 0.4 threshold).

### Suites (`suite.h`)

```cpp
class AssessmentSuite : public Assessment {
public:
    AssessmentSuite(std::string name, std::string description,
                    std::vector<std::shared_ptr<Assessment>> children);
    std::string GetName() const override;
    std::string GetDescription() const override;
    const std::vector<std::shared_ptr<Assessment>>& GetChildren() const;
    AssessmentResult Assess(const StructureContext& ctx) const override;
};

std::shared_ptr<AssessmentSuite> make_protein_geometry_suite();
```

- The constructor throws `std::invalid_argument` if any child is null or if
  the name is empty.
- `Assess` runs each child in order against the same context and stores the
  child results in `children`. The suite's own `metrics` and `findings` are
  empty.
- If a child throws a `std::exception`, the suite records a child result
  with `assessment` set to the child's `GetName()`, `status` ERROR,
  `message` set to the exception's `what()`, and empty metrics, findings,
  and children. It then runs the remaining children.
- Suite status is the `worst_status` of its children; a suite with no
  children is NOT_APPLICABLE. A suite's own `message` is empty unless its
  status is NOT_APPLICABLE.
- Suites nest, since a suite is an `Assessment`.

`make_protein_geometry_suite()` returns a suite named `protein_geometry`
containing, in order: bond lengths, bond angles, clashes, valence,
Ramachandran, each with default options.

## Assessments

The geometry and Ramachandran checks consider only standard amino acids
(as defined above); modified residues and ligands are skipped. Messages
name the chain, residue name, residue number with insertion code, and atom
names.

### Geometry targets (shared by `bond_lengths` and `bond_angles`)

Targets come from the generated tables (see Reference data). Heavy atoms
only: terms involving hydrogens are not assessed. A term with a missing atom
is skipped and not counted.

Intra-residue terms: every bond and angle row for the residue. Rows with a
variant description are selected per residue:

- CYS: "SS" when SG is within 2.5 Å of the SG of another CYS residue,
  otherwise "SH".
- PRO: "cis" when the preceding peptide is cis, otherwise "trans". The
  preceding peptide is cis when PRO is the peptide successor of residue i
  and |ω| < 90°, with ω the CA(i)–C(i)–N–CA torsion. A PRO with no peptide
  predecessor uses "trans".
- HIS: from explicit hydrogens when present: HD1 and HE2 both present gives
  "HISH", only HD1 gives "HISD", only HE2 gives "HISE". With neither
  present, evaluate each of the three variants and use the one with the
  lowest RMSZ over that residue's variant rows (ties go to HISE, then HISD,
  then HISH). The same choice applies to the residue's bonds and angles.

Peptide-link terms, for each residue i with peptide successor j:

- Bond C(i)–N(j); angles CA(i)–C(i)–N(j), O(i)–C(i)–N(j), C(i)–N(j)–CA(j);
  and, when j is PRO, C(i)–N(j)–CD(j).
- Variant selection (by `LinkVariant`, see Reference data): if j is PRO,
  use `PRO` for the C(i)–N(j) bond, CA(i)–C(i)–N(j), and O(i)–C(i)–N(j),
  and `PRO_CIS` or `PRO_TRANS` by ω (as above) for C(i)–N(j)–CA(j) and
  C(i)–N(j)–CD(j); `PRO_ALL` is not used. Otherwise, if the residue that
  owns the varying atom is GLY, use `GLY`: residue i for CA(i)–C(i)–N(j)
  and O(i)–C(i)–N(j), residue j for C(i)–N(j)–CA(j) and for the C(i)–N(j)
  bond. Otherwise use `DEFAULT`.
- C(i)–N(j)–CD(j) is assessed only here, once per PRO that has a peptide
  predecessor; a PRO without one has no C–N–CD term.
- Assumption: the GLY row for the C–N bond is keyed on residue j, validate
  via comparison with the Engh & Huber (2001) table text. The C–N GLY and
  default targets differ by 0.010 Å (about half a sigma), so a wrong choice
  shifts Z without creating outliers on ordinary structures.

Disulfide terms: for each pair of CYS residues whose SG atoms are within
2.5 Å, the SG–SG bond from the link table. This is a bond term only.

### `bond_lengths`: `BondLengthAssessment`

- Options: `warning_z = 4.0`, `error_z = 10.0`.
- Z = (observed − ideal) / sigma for every term above.
- |Z| ≥ `warning_z` gives a WARNING finding and |Z| ≥ `error_z` an ERROR
  finding, code `bond_length_outlier`: two atoms; `observed` and `expected`
  in Å; `deviation` = Z; no details.
- Metrics, in order: `n_assessed`, `n_outliers` (findings of either
  severity), `rmsz` (root mean square of Z over all assessed terms),
  `max_abs_z`.
- NOT_APPLICABLE when no term was assessed.

### `bond_angles`: `BondAngleAssessment`

- Same options, rules, and metrics as `bond_lengths`, applied to angle
  terms, in degrees.
- Code `bond_angle_outlier`: three atoms in angle order (vertex second).

### `clashes`: `ClashAssessment`

- Options: `warning_overlap = 0.4`, `error_overlap = 1.0`,
  `polar_contact_distance = 2.9` (all in Å).
- Scope: every atom in the molecule, including ligands and waters, except
  atoms with atomic number 0 and atoms for which `OEIsMetal` is true.
- Each pair has a contact distance: the sum of the two atoms' radii from
  `OEGetBondiVdWRadius`, except that pairs where both atoms are N or O
  (possible hydrogen bonds) use `polar_contact_distance`. Overlap = contact
  distance − d.
- Excluded pairs: atoms within three bonds of each other (1-2, 1-3, 1-4),
  using the molecule's bonds.
- Candidate pairs come from `OEChem::OENearestNbrs` with cutoff
  max(2 × largest radius among in-scope atoms, `polar_contact_distance`) −
  `warning_overlap` + 0.01 Å. The 0.01 Å margin keeps pairs exactly at the
  warning threshold in the candidate set whether or not `OENearestNbrs`
  treats its cutoff as inclusive; the overlap test below decides which
  candidates are findings.
- Overlap ≥ `warning_overlap` gives a WARNING finding, ≥ `error_overlap` an
  ERROR finding, code `clash`: two atoms; `observed` = d; `expected` =
  contact distance; `deviation` = overlap; no details. With the defaults, a
  polar pair is flagged at 2.5 Å or closer (WARNING) and at 1.9 Å or closer
  (ERROR); an N–O pair at 2.8 Å is not flagged.
- Metrics, in order: `n_clashes` (findings of either severity),
  `clashscore` (`n_clashes` × 1000 / number of in-scope atoms),
  `max_overlap` (largest overlap among findings, 0 when there are none),
  `hydrogens_present` (1 or 0).
- NOT_APPLICABLE when no atom is in scope; `hydrogens_present` is still
  reported.

### `valence`: `ValenceAssessment`

- No options.
- Runs `OEChem::OECheckAtomValence` on every atom. Each failure is an ERROR
  finding, code `invalid_valence`: one atom; `observed` = the atom's
  valence (`GetValence()`, which includes implicit hydrogens); `expected`
  and `deviation` NaN; `details` = `formal_charge`, `implicit_h_count`,
  `degree`. The message gives element, formal charge, and valence.
- Requirement, stated in the documentation: implicit hydrogen counts must
  be assigned, as OpenEye's readers do by default. With zeroed implicit
  counts, normal atoms are reported as under-valent.
- Metrics: `n_invalid`.
- NOT_APPLICABLE when the molecule has no atoms.

### `ramachandran`: `RamachandranAssessment`

- No options; categories come from OEBio.
- For each standard amino-acid residue, constructs
  `OEBio::OERamachandranAnalysis(mol, residue)` using the residue of its CA
  atom. Residues with no CA, or whose category is `OERamaCategory::Unknown`
  (termini and breaks), are not assessed.
- Each Outlier is a WARNING finding, code `ramachandran_outlier`: the
  residue's N, CA, and C atoms (those present); `observed` = the OEBio
  score; `expected` and `deviation` NaN; `details` = `phi` and `psi` in
  degrees (OEBio returns radians), `score`, and `rama_type` (the numeric
  `OERamaType` value). The message includes phi, psi, and the type name
  from `OEGetRamachandranTypeName`.
- Metrics, in order: `n_assessed`, `n_favored`, `n_allowed`, `n_outliers`,
  `pct_favored` and `pct_outliers` (percent of `n_assessed`, 0 to 100).
- NOT_APPLICABLE when no residue was assessed.

## Reference data

`scripts/generate_engh_huber.py` reads `standard_geometry.cif` (path given
on the command line) and writes `src/data/engh_huber_2001.cpp`. It uses only
the Python standard library and runs under the repository's Python tooling
(ruff, mypy). It extracts, for the 20 standard amino acids, heavy-atom rows
from `_chem_comp_bond` and `_chem_comp_angle` (keeping the description
column as the variant: ".", "SS", "SH", "cis", "trans", "HISD", "HISE",
"HISH"), and the peptide (C–N, by link id) and disulfide (SG–SG) rows from
`_chem_link_bond` and `_chem_link_angle`.

Two corrections are applied while extracting:

- Intra-residue rows in which C and N are consecutive atoms (an X–C–N or
  C–N–X pattern) describe the bond to the neighboring residue, not a bond
  inside the residue, and are dropped. In the source file the only such rows
  are PRO C–N–CD ("cis" and "trans"); that angle is assessed as a link term
  from the link table.
- Link ids are normalized to an enum `LinkVariant { DEFAULT, GLY, PRO,
  PRO_CIS, PRO_TRANS, PRO_ALL, DISULFIDE }`. The source spells the default
  row both "NOT GLY PRO" and "NOT PRO GLY"; both map to `DEFAULT`. "CYS"
  maps to `DISULFIDE`. Nucleic-acid link rows are skipped. Any other link id
  on a C–N or SG–SG row makes the script exit with an error naming the row.

The generated file holds `constexpr` arrays of plain structs:

- intra-residue bonds `{residue, atom1, atom2, variant, ideal, sigma}`
- intra-residue angles `{residue, atom1, atom2, atom3, variant, ideal, sigma}`
- link bonds `{atom1, atom2, link_variant, ideal, sigma}`
- link angles `{atom1, atom2, atom3, link_variant, ideal, sigma}`

`src/data/engh_huber_2001.h` declares the structs and accessors. The
generated file's header comment records the Engh & Huber citation, the
source file name and its `_standard_geometry.version` string, and the
command that generated it. `src/geometry/targets.{h,cpp}` implement the
variant selection described under Assessments.

Assumption: redistributing these values in this MIT-licensed repository is
acceptable because they are the published Engh & Huber (2001) numbers;
the project owner validates via review before the first public release.

## JSON (`json.h`)

```cpp
std::string to_json(const AssessmentResult& result, int indent = 2);
```

Output for one result:

```json
{
  "schema_version": 1,
  "oehealth_version": "0.1.0",
  "assessment": "bond_lengths",
  "status": "WARN",
  "message": "",
  "metrics": {"n_assessed": 412, "n_outliers": 1, "rmsz": 0.83, "max_abs_z": 4.6},
  "findings": [
    {
      "code": "bond_length_outlier",
      "severity": "WARNING",
      "message": "A ALA 5 CA-CB 1.617 A (ideal 1.520, Z 4.6)",
      "atoms": [
        {"atom_idx": 10, "atom_name": "CA", "residue_name": "ALA",
         "residue_number": 5, "insert_code": "", "chain_id": "A"},
        {"atom_idx": 13, "atom_name": "CB", "residue_name": "ALA",
         "residue_number": 5, "insert_code": "", "chain_id": "A"}
      ],
      "observed": 1.617,
      "expected": 1.520,
      "deviation": 4.6,
      "details": {}
    }
  ],
  "children": []
}
```

- `schema_version` and `oehealth_version` appear only on the root object;
  children omit them.
- Metrics and details are JSON objects that keep insertion order
  (`nlohmann::ordered_json`). `to_json` throws `std::invalid_argument` if a
  metric or detail name repeats within one vector.
- NaN serializes as `null`.
- Status and severity serialize as their enum names.
- A negative `indent` produces compact single-line output, matching
  nlohmann's `dump` convention.

## Python bindings

`swig/oehealth.i` adds `std_vector.i` and `std_shared_ptr.i`, and declares:

- `%shared_ptr` for `Assessment`, `AssessmentSuite`, and every concrete
  assessment;
- vector templates for `AtomRef`, `Metric`, `Finding`, `AssessmentResult`,
  and `std::shared_ptr<Assessment>` (for the suite constructor);
- the result types, options structs, concrete assessments, the suite,
  `worst_status`, `make_protein_geometry_suite`, and `to_json`;
- `%ignore` for `StructureContext`, `ResidueInfo`, and
  `Assessment::Assess`;
- a `%pythoncode`/`%extend` addition giving `AssessmentResult.to_dict()`,
  which returns `json.loads(to_json(self))`.

The `%exception` block gains a `std::invalid_argument` catch, ahead of the
`std::exception` catch, that raises `ValueError`.
`python/oehealth/__init__.py` re-exports the public names in `__all__` and
drops `calculate_molecular_weight`.

## Error handling

- Nothing to assess (empty molecule, no standard amino acids, no in-scope
  atoms): NOT_APPLICABLE with a message.
- Missing atoms inside a residue: affected terms are skipped and not
  counted.
- Molecule with atoms but no 3D coordinates: `Run` throws
  `std::invalid_argument`; Python sees `ValueError`.
- Invalid options or a null suite child: the constructor throws
  `std::invalid_argument`; Python sees `ValueError`.
- A child throws inside a suite: recorded as an ERROR child as described
  under Suites; the remaining children still run.

## Testing

C++ tests use GoogleTest. `tests/data/1crn.pdb` is PDB entry 1CRN
(crambin, 46 residues, 327 heavy atoms, no hydrogens). The test CMake
passes the data directory to the tests as a compile definition.

Framework:

- `worst_status` follows declaration order.
- Status derivation: no findings gives PASS, a WARNING finding gives WARN,
  an ERROR finding gives FAIL, INFO alone gives PASS.
- Suites: a test-only assessment that throws yields one ERROR child with
  the child's name and the exception message while the other children still
  run; status is the worst child status; an empty suite is NOT_APPLICABLE;
  a null child throws `std::invalid_argument`; `GetName` and
  `GetDescription` return the constructor values.
- `Run` throws `std::invalid_argument` for a 2D molecule with atoms, and an
  empty molecule gives NOT_APPLICABLE from every assessment.
- Options validation: NaN, infinite, negative, and reversed thresholds
  throw, as does a zero `polar_contact_distance`; equal thresholds are
  accepted, and a value at the shared threshold is reported as ERROR;
  default-constructed assessments report the documented defaults through
  `GetOptions()`.
- Context: crambin has 46 residues in file order, all standard; 45 peptide
  successor links; deleting atom N of residue 10 removes the 9→10 link.

Reference data:

- Spot values from Engh & Huber (2001), exact to the table's precision:
  ALA N–CA 1.459 ± 0.020, CA–C 1.525 ± 0.026, C–O 1.229 ± 0.019,
  CA–CB 1.520 ± 0.021 (Å); ALA N–CA–C 111.0 ± 2.7, CA–C–O 120.1 ± 2.1
  (degrees); GLY N–CA–C 113.1 ± 2.5; peptide C–N 1.336 ± 0.023,
  CA–C–N 117.2 ± 2.2, O–C–N 122.7 ± 1.6, C–N–CA 121.7 ± 2.5 ("NOT GLY PRO"
  rows); PRO CA–N–CD 111.7 ± 1.4 (trans) and 111.5 ± 1.4 (cis).
- Coverage: the test lists each standard residue's PDB heavy-atom names
  (for example ALA: N, CA, C, O, CB, OXT) and checks that the residue's bond
  rows form one connected graph over exactly those names, and that every
  angle row's three atoms are among them.
- Variant selection: CYS SS versus SH on crambin's disulfides; PRO cis
  versus trans by ω; HIS by explicit hydrogens and by best fit.
- Link lookup: every peptide term resolves to exactly one row for each
  variant it can take: C–N, CA–C–N, and O–C–N for `DEFAULT`, `GLY`, and
  `PRO`; C–N–CA for `DEFAULT`, `GLY`, `PRO_CIS`, and `PRO_TRANS`; C–N–CD
  for `PRO_CIS` and `PRO_TRANS`.
- The intra-residue PRO angle table has no C–N–CD row. On crambin (PRO at
  residues 5, 19, 22, 36, and 41, none N-terminal) exactly five C–N–CD
  terms are assessed, and a PRO with its predecessor removed gets none.
- The generator script has a pytest test that runs it on a small excerpt of
  `standard_geometry.cif` committed under `tests/data/` and checks the
  emitted rows, the "NOT PRO GLY" alias, the dropped PRO C–N–CD rows, and
  the error on an unknown link id.

Assessments:

- Crambin as deposited produces no ERROR findings from any assessment.
- `bond_lengths`: moving one atom to stretch one bond produces exactly one
  outlier naming the right two atoms with the expected Z.
- `bond_angles`: the same with one distorted angle.
- `clashes`: two unbonded carbons 2.0 Å apart give one ERROR clash; a close
  1-4 pair gives none; an N–O pair at 2.8 Å gives none and at 2.4 Å gives a
  WARNING; a Zn near an N gives none. Exact boundaries: two unbonded
  carbons at exactly 3.0 Å (overlap 0.4) give a WARNING, at exactly 2.4 Å
  (overlap 1.0) an ERROR, and at 3.01 Å nothing.
- `bond_lengths`: a bond placed at exactly ideal + 4σ gives a WARNING.
- `valence`: a five-coordinate carbon gives one `invalid_valence` finding;
  crambin as read gives none; crambin with implicit hydrogen counts zeroed
  gives findings (documents the requirement).
- `ramachandran`: crambin gives 44 assessed residues (termini excluded) and
  no outliers, matching the wwPDB validation report for 1CRN; setting one
  residue's phi and psi into a disallowed region with `OESetTorsion` gives
  one outlier at that residue, with phi and psi in degrees.
- JSON: output parses back; keys and metric order match the schema; NaN
  becomes `null`; duplicate metric names throw; child objects omit the
  version keys.

Python tests (pytest) cover the bindings: `Run` on crambin for each
assessment, `to_dict()`, the protein geometry suite, constructing a custom
suite from Python, setting options fields from Python, `ValueError` for 2D
input and for invalid options, and `TypeError` for a non-molecule.

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
