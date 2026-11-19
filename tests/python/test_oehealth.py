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

    def test_vector_fields_are_owning_copies(self, crambin):
        # SWIG converts std::vector members to tuples of owning element proxies
        # (naturalvar). A regression to references would make these proxies
        # non-owning views into the result, which they would outlive.
        result = oehealth.BondLengthAssessment().Run(crambin)
        assert isinstance(result.findings, tuple)
        assert isinstance(result.metrics, tuple)
        assert len(result.findings) > 0
        assert result.findings[0].this.own()
        assert isinstance(result.findings[0].atoms, tuple)
        assert result.findings[0].atoms[0].this.own()
        assert result.metrics[0].this.own()

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
        assert options.hbond_contact_distance == pytest.approx(1.9)

    def test_set_fields(self):
        options = oehealth.BondLengthOptions()
        options.warning_z = 3.0
        options.error_z = 6.0
        returned = oehealth.BondLengthAssessment(options).GetOptions()
        gc.collect()
        assert returned.warning_z == pytest.approx(3.0)
        assert returned.error_z == pytest.approx(6.0)

    def test_get_options_returns_an_owning_copy(self):
        options = oehealth.ClashAssessment().GetOptions()
        assert options.this.own()
        # Mutating the copy must not change the assessment that produced it.
        options.warning_overlap = 0.9
        assert oehealth.ClashAssessment().GetOptions().warning_overlap == pytest.approx(0.4)

    def test_invalid_options_raise_value_error(self):
        options = oehealth.ClashOptions()
        options.polar_contact_distance = 0.0
        with pytest.raises(ValueError):
            oehealth.ClashAssessment(options)
        bad_hbond = oehealth.ClashOptions()
        bad_hbond.hbond_contact_distance = 0.0
        with pytest.raises(ValueError):
            oehealth.ClashAssessment(bad_hbond)
        reversed_thresholds = oehealth.BondAngleOptions()
        reversed_thresholds.warning_z = 5.0
        reversed_thresholds.error_z = 4.0
        with pytest.raises(ValueError):
            oehealth.BondAngleAssessment(reversed_thresholds)


class TestSuites:
    """Suites compose assessments in C++ and from Python."""

    def test_get_children_returns_a_tuple_of_owning_proxies(self):
        children = oehealth.make_protein_geometry_suite().GetChildren()
        gc.collect()
        assert isinstance(children, tuple)
        assert len(children) == len(ASSESSMENTS)
        assert all(child.this.own() for child in children)
        assert children[0].GetName() == ASSESSMENTS[0][1]

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


STATUS_ORDER = [
    oehealth.Status_NOT_APPLICABLE,
    oehealth.Status_PASS,
    oehealth.Status_WARN,
    oehealth.Status_FAIL,
    oehealth.Status_ERROR,
]


class TestWorstStatus:
    """worst_status returns the later of two statuses in declaration order."""

    @pytest.mark.parametrize("first", range(5))
    @pytest.mark.parametrize("second", range(5))
    def test_returns_the_later_status(self, first, second):
        a, b = STATUS_ORDER[first], STATUS_ORDER[second]
        assert oehealth.worst_status(a, b) == STATUS_ORDER[max(first, second)]

    def test_error_dominates_and_not_applicable_yields(self):
        assert oehealth.worst_status(oehealth.Status_ERROR, oehealth.Status_PASS) == (
            oehealth.Status_ERROR
        )
        assert oehealth.worst_status(oehealth.Status_NOT_APPLICABLE, oehealth.Status_PASS) == (
            oehealth.Status_PASS
        )
