import unittest
import pandas as pd
import os

# Import functions
from glycopeptide_sequence_finder_cmd import (
    cleave_sequence,
    find_glycopeptides,
    calculate_peptide_mass,
    predict_hydrophobicity,
    calculate_pI,
    compute_mz,
    process_glycopeptides,
    process_fasta,
    write_csv
)

class TestGlycopeptideSequenceFinder(unittest.TestCase):
    """Unit tests for glycopeptide_sequence_finder_cmd.py."""

    def test_cleave_sequence(self):
        """Test the cleave_sequence function."""
        sequence = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"
        protease = "trypsin"
        expected_peptides = ["MK", "WVTFISLLFLFSSAYSR", "GVFR", "R", "DTHK", "SEIAHR", "FK", "DLGE"]
        result = cleave_sequence(sequence, protease)
        self.assertEqual(result, expected_peptides)

    def test_find_glycopeptides(self):
        """Test the find_glycopeptides function."""
        peptides = ["MKWVTFISLLFLFSSAYSR", "GVFR", "RDTHK", "SEIAHR", "FKDLGE"]
        full_sequence = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"
        glycosylation_type = "N"
        expected_glycopeptides = []
        result = find_glycopeptides(peptides, full_sequence, glycosylation_type)
        self.assertEqual(result, expected_glycopeptides)

    def test_calculate_peptide_mass(self):
        """Test the calculate_peptide_mass function."""
        sequence = "MKWVTFISLLFLFSSAYSR"
        expected_mass = 2295.21  # Approximate mass
        result = calculate_peptide_mass(sequence)
        self.assertAlmostEqual(result, expected_mass, places=2)

    def test_predict_hydrophobicity(self):
        """Test the predict_hydrophobicity function."""
        peptide_sequence = "MKWVTFISLLFLFSSAYSR"
        expected_hydrophobicity = 0.93  # Approximate value
        result = predict_hydrophobicity(peptide_sequence)
        self.assertAlmostEqual(result, expected_hydrophobicity, places=2)

    def test_calculate_pI(self):
        """Test the calculate_pI function."""
        peptide_sequence = "MKWVTFISLLFLFSSAYSR"
        expected_pI = 9.99  # Approximate value
        result = calculate_pI(peptide_sequence)
        self.assertAlmostEqual(result, expected_pI, places=2)

    def test_compute_mz(self):
        """Test the compute_mz function."""
        mass = 2376.27
        charge = 2
        expected_mz = 1189.14  # Approximate value
        result = compute_mz(mass, charge)
        self.assertAlmostEqual(result, expected_mz, places=2)

    def test_process_glycopeptides(self):
        """Test the process_glycopeptides function."""
        peptide_data = {
            'ProteinID': ['P12345'],
            'Site': [1],
            'Peptide': ['MKWVTFISLLFLFSSAYSR'],
            'Start': [1],
            'End': [20],
            'Length': [20],
            'Sequon': ['N/A'],
            'PredictedMass': [2376.27],
            'Hydrophobicity': [1.18],
            'pI': [9.88],
            'Protease': ['trypsin'],
            'GlycosylationType': ['N'],
            'MissedCleavages': [0],
            'Species': ['Homo sapiens'],
            'TaxonID': [9606],
            'GeneName': ['ALB'],
            'ProteinEvidence': [1],
            'SequenceVersion': [1]
        }
        glycan_data = {
            'glytoucan_ac': ['G80920RR'],
            'byonic': ['HexNAc(2)Hex(9) % 1864.634157'],
            'composition': ['HexNAc(2)Hex(9)'],
            'mass': [1864.634157]
        }
        peptides = pd.DataFrame(peptide_data)
        glycans = pd.DataFrame(glycan_data)
        max_charge = 5

        peptides.to_csv('test_peptides.csv', index=False)
        glycans.to_csv('test_glycans.csv', index=False)

        result = process_glycopeptides('test_peptides.csv', glycans, max_charge)
        self.assertEqual(len(result), 1)
        self.assertIn('z2', result.columns)

        os.remove('test_peptides.csv')
        os.remove('test_glycans.csv')

    def test_process_fasta(self):
        """Test the process_fasta_function."""
        fasta_content = """>sp|P12345|ALBU_HUMAN Serum albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=1
MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGE
"""
        with open('test.fasta', 'w') as f:
            f.write(fasta_content)

        protease = "trypsin"
        missed_cleavages = 0
        glycosylation_type = "N"
        result = process_fasta('test.fasta', protease, missed_cleavages, glycosylation_type)
        self.assertEqual(len(result), 0)

        os.remove('test.fasta')

    def test_write_csv(self):
        """Test write_csv function."""
        data = [
            {
                "ProteinID": "P12345",
                "Site": 1,
                "Peptide": "MKWVTFISLLFLFSSAYSR",
                "Start": 1,
                "End": 20,
                "Length": 20,
                "Sequon": "N/A",
                "PredictedMass": 2376.27,
                "Hydrophobicity": 1.18,
                "pI": 9.88,
                "Protease": "trypsin",
                "GlycosylationType": "N",
                "MissedCleavages": 0,
                "Species": "Homo sapiens",
                "TaxonID": 9606,
                "GeneName": "ALB",
                "ProteinEvidence": 1,
                "SequenceVersion": 1
            }
        ]
        output_file = 'test_output.csv'
        write_csv(output_file, data)
        self.assertTrue(os.path.exists(output_file))

        with open(output_file, 'r') as f:
            content = f.read()
            self.assertIn("ProteinID", content)
            self.assertIn("P12345", content)

        os.remove(output_file)

if __name__ == '__main__':
    unittest.main(verbosity=2)
