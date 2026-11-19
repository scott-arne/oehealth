#include <gtest/gtest.h>

#include "oehealth/assessments/ramachandran.h"
#include "test_helpers.h"

#include <oebio.h>

#include <array>
#include <cstdio>
#include <fstream>
#include <string>

namespace OEHEALTH::test {

namespace {

/// Return crambin PDB text in which the CA of one residue has A and B alternate locations.
std::string crambin_with_alternate_ca(int residue_number) {
    std::ifstream file(std::string(OEHEALTH_TEST_DATA_DIR) + "/1crn.pdb");
    std::string text;
    std::string line;
    while (std::getline(file, line)) {
        const bool is_target = line.compare(0, 6, "ATOM  ") == 0 && line.size() >= 54 &&
                               trimmed(line.substr(12, 4)) == "CA" &&
                               std::stoi(line.substr(22, 4)) == residue_number;
        if (!is_target) {
            text.append(line).append("\n");
            continue;
        }
        // Columns 31-38 hold x; the B conformer is shifted by 0.3 A along it.
        std::string a_line = line;
        a_line[16] = 'A';
        std::string b_line = a_line;
        b_line[16] = 'B';
        std::array<char, 9> x{};
        std::snprintf(x.data(), x.size(), "%8.3f", std::stod(line.substr(30, 8)) + 0.3);
        b_line.replace(30, 8, x.data());
        text.append(a_line).append("\n").append(b_line).append("\n");
    }
    return text;
}

} // namespace

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

TEST(RamachandranAssessment, ResidueWithAlternateCaIsStillAssessed) {
    OEChem::OEGraphMol plain;
    read_crambin(plain);
    const AssessmentResult baseline = RamachandranAssessment().Run(plain);

    // The alternate must come from PDB text read with the ALTLOC flavor: only then does OEBio's
    // (mol, residue) constructor return Unknown for the residue, which is the case the atom
    // constructor in the assessment is meant to handle.
    OEChem::OEGraphMol mol;
    OEChem::oemolistream ifs;
    ifs.SetFlavor(OEChem::OEFormat::PDB, OEChem::OEIFlavor::Generic::Default |
                                             OEChem::OEIFlavor::PDB::Default |
                                             OEChem::OEIFlavor::PDB::ALTLOC);
    ifs.SetFormat(OEChem::OEFormat::PDB);
    ASSERT_TRUE(ifs.openstring(crambin_with_alternate_ca(20)));
    ASSERT_TRUE(OEChem::OEReadMolecule(ifs, mol));

    OEChem::OEAtomBase* ca = find_atom(mol, 20, "CA");
    ASSERT_NE(ca, nullptr);
    EXPECT_EQ(OEChem::OEAtomGetResidue(ca).GetAlternateLocation(), 'A');
    EXPECT_EQ(OEBio::OERamachandranAnalysis(mol, OEChem::OEAtomGetResidue(ca)).GetRamaCategory(),
              OEBio::OERamaCategory::Unknown);
    const unsigned int first_ca_category = OEBio::OERamachandranAnalysis(ca).GetRamaCategory();
    ASSERT_NE(first_ca_category, OEBio::OERamaCategory::Unknown);

    const AssessmentResult result = RamachandranAssessment().Run(mol);
    EXPECT_EQ(metric(result, "n_assessed"), 44.0);
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
