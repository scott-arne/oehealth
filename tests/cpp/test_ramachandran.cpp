#include <gtest/gtest.h>

#include "oehealth/assessments/ramachandran.h"
#include "test_helpers.h"

#include <oebio.h>

#include <array>

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
    // Count the categories directly with OEBio's (mol, residue) constructor, so a swap of
    // Favored and Allowed in the assessment cannot pass.
    double favored = 0.0;
    double allowed = 0.0;
    double outliers = 0.0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (trimmed(atom->GetName()) != "CA") {
            continue;
        }
        const OEBio::OERamachandranAnalysis analysis(mol, OEChem::OEAtomGetResidue(atom));
        const unsigned int category = analysis.GetRamaCategory();
        if (category == OEBio::OERamaCategory::Favored) {
            favored += 1.0;
        } else if (category == OEBio::OERamaCategory::Allowed) {
            allowed += 1.0;
        } else if (category == OEBio::OERamaCategory::Outlier) {
            outliers += 1.0;
        }
    }
    EXPECT_EQ(metric(result, "n_favored"), favored);
    EXPECT_EQ(metric(result, "n_allowed"), allowed);
    EXPECT_EQ(metric(result, "n_outliers"), outliers);
    EXPECT_DOUBLE_EQ(metric(result, "pct_favored"),
                     100.0 * favored / (favored + allowed + outliers));
    EXPECT_EQ(metric(result, "pct_outliers"), 0.0);
}

TEST(RamachandranAssessment, ResidueWithAlternateCaIsClassifiedFromItsFirstCa) {
    OEChem::OEGraphMol plain;
    read_crambin(plain);
    const AssessmentResult baseline = RamachandranAssessment().Run(plain);

    OEChem::OEGraphMol mol;
    read_crambin(mol);
    OEChem::OEAtomBase* ca = find_atom(mol, 20, "CA");
    ASSERT_NE(ca, nullptr);
    // This builds the alternate in memory; OEBio's (mol, residue) constructor still classifies
    // the residue here, so the test pins the new behavior and not the old Unknown result.
    // The new atom comes after the original in index order, so the original stays the first CA.
    OEChem::OEResidue alternate = OEChem::OEAtomGetResidue(ca);
    OEChem::OEResidue original = alternate;
    original.SetAlternateLocation('A');
    alternate.SetAlternateLocation('B');
    OEChem::OEAtomSetResidue(ca, original);
    const std::array<double, 3> xyz = coords(mol, ca);
    OEChem::OEAtomBase* second =
        add_atom(mol, OEChem::OEElemNo::C, "CA", "", 20, xyz[0] + 0.3, xyz[1], xyz[2]);
    OEChem::OEAtomSetResidue(second, alternate);
    mol.NewBond(second, find_atom(mol, 20, "N"), 1);
    mol.NewBond(second, find_atom(mol, 20, "C"), 1);

    const AssessmentResult result = RamachandranAssessment().Run(mol);
    EXPECT_EQ(metric(result, "n_assessed"), metric(baseline, "n_assessed"));
    EXPECT_EQ(metric(result, "n_favored"), metric(baseline, "n_favored"));
    EXPECT_EQ(metric(result, "n_allowed"), metric(baseline, "n_allowed"));
    EXPECT_EQ(metric(result, "n_outliers"), metric(baseline, "n_outliers"));
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
