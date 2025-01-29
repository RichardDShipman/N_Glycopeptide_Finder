"""
n_glycopeptide_finder_cmd.py

This script processes a FASTA file to identify N-glycopeptides based on protease cleavage rules, 
predicts their masses, calculates their hydrophobicity, and calculates their pI. The results are written to a CSV file.

Functions:
    cleave_sequence(sequence, protease, missed_cleavages):
        Cleaves a sequence based on protease rules.
    
    find_n_glycopeptides(peptides, full_sequence):
        Identifies peptides containing N-glycosylation sequons and maps the sites to the full protein sequence.
    
    calculate_peptide_mass(sequence):
        Calculates the mass of a peptide.
    
    predict_hydrophobicity(peptide_sequence):
        Predicts the hydrophobicity of a peptide sequence using the Kyte-Doolittle scale.
    
    calculate_pI(peptide_sequence):
        Calculates the isoelectric point (pI) of a peptide sequence.
    
    process_fasta(file, protease, missed_cleavages):
        Processes the input FASTA file and extracts glycopeptides.
    
    write_csv(output_file, data):
        Writes results to a CSV file.
    
    main():
        Main function to parse arguments and execute the glycopeptide finding process.

Usage:
    python n_glycopeptide_finder_cmd.py -i <input_fasta_file> -o <output_csv_file> -p <protease> -c <missed_cleavages>

Author:
    Richard Shipman
"""
import argparse
import csv
import re
from Bio import SeqIO
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from pyteomics.mass import calculate_mass

# Define protease cleavage rules
proteases = {
    "trypsin": ("[KR]", "P"),  # Cleaves after K or R unless followed by P
    "chymotrypsin": ("[FWY]", "P"),  # Cleaves after F, W, or Y unless followed by P
    "glu-c": ("E", None),  # Cleaves after E
    "lys-c": ("K", None),  # Cleaves after K
    "arg-c": ("R", None),  # Cleaves after R
    "pepsin": ("[FLWY]", None),  # Cleaves after F, W, or Y
    "asp-n": ("D", None),  # Cleaves **before** Asp (D)
    "proteinase-k": ("[AFILVWY]", None),  # Cleaves after A, F, I, L, V, W, Y
}

def cleave_sequence(sequence, protease, missed_cleavages=0):
    """Cleaves a sequence based on protease rules."""
    cleavage_pattern, exclusion = proteases[protease.lower()]

    # Handle Asp-N separately because it cleaves **before** D
    if protease.lower() == "asp-n":
        regex = rf"(?={cleavage_pattern})"  # Cleaves before D
    else:
        regex = rf"(?<={cleavage_pattern})(?!{exclusion})"  # Cleaves after other residues

    # Perform cleavage
    fragments = re.split(regex, sequence)

    # Generate peptides including missed cleavages
    peptides = []
    for i in range(len(fragments)):
        for j in range(i + 1, min(i + 2 + missed_cleavages, len(fragments) + 1)):
            peptides.append("".join(fragments[i:j]))
    
    return peptides

def find_n_glycopeptides(peptides, full_sequence):
    """Identifies peptides containing N-glycosylation sequons and maps the sites to the full protein sequence."""
    n_glyco_sequon = re.compile("N[^P][ST]")
    glycopeptides = []
    for pep in peptides:
        for match in n_glyco_sequon.finditer(pep):
            # Map the position in the peptide to the full protein sequence
            start_in_protein = full_sequence.find(pep) + match.start() + 1  # Adjust to 1-based indexing
            glycopeptides.append((pep, start_in_protein))
    return glycopeptides

def calculate_peptide_mass(sequence):
    """Calculates the mass of a peptide, ignoring sequences with unknown residues."""
    invalid_residues = {"X", "B", "Z", "J", "U", "O"}  # Common ambiguous residues
    if any(aa in invalid_residues for aa in sequence):
        return "Unknown"  # Or return None if you prefer
    return calculate_mass(sequence=sequence)

def predict_hydrophobicity(peptide_sequence):
    """Predicts the hydrophobicity of a peptide sequence using the Kyte-Doolittle scale."""
    hydrophobicity_values = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }

    total_hydrophobicity = sum(hydrophobicity_values.get(aa, 0) for aa in peptide_sequence)
    average_hydrophobicity = round(total_hydrophobicity / len(peptide_sequence), 2)

    return average_hydrophobicity

def calculate_pI(peptide_sequence):
    """Calculates the isoelectric point (pI) of a peptide sequence."""
    pI_calculator = IsoelectricPoint(peptide_sequence)
    return round(pI_calculator.pi(), 2)

def process_fasta(file, protease, missed_cleavages):
    """Processes the input FASTA file and extracts glycopeptides."""
    results = []
    for record in SeqIO.parse(file, "fasta"):
        protein_id = record.id
        sequence = str(record.seq)
        peptides = cleave_sequence(sequence, protease, missed_cleavages)
        n_glycopeptides = find_n_glycopeptides(peptides, sequence)
        for peptide, site in n_glycopeptides:
            mass = calculate_peptide_mass(peptide)
            hydrophobicity = predict_hydrophobicity(peptide)
            pI = calculate_pI(peptide)
            start_pos = sequence.find(peptide) + 1  # 1-based indexing
            end_pos = start_pos + len(peptide)
            results.append({
                "ProteinID": protein_id,
                "Site": site,
                "Peptide": peptide,
                "Start": start_pos,
                "End": end_pos,
                "Length": len(peptide),
                "NSequon": sequence[site - 1:site + 2],
                "PredictedMass": mass,
                "Hydrophobicity": hydrophobicity,
                "pI": pI,
            })
    return results

def write_csv(output_file, data):
    """Writes results to a CSV file."""
    with open(output_file, mode="w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=["ProteinID", "Site", "Peptide", "Start", "End", "Length", "NSequon", "PredictedMass", "Hydrophobicity", "pI"])
        writer.writeheader()
        writer.writerows(data)

def main():
    parser = argparse.ArgumentParser(description="Glycopeptide Finder")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", help="Output CSV file prefix")
    parser.add_argument("-p", "--protease", default="trypsin", help="Protease to use for cleavage ('all' for all proteases)")
    parser.add_argument("-c", "--missed_cleavages", type=int, default=0, help="Number of missed cleavages allowed")

    args = parser.parse_args()

    input_file = args.input
    base_filename = input_file.rsplit(".", 1)[0]
    missed_cleavages = args.missed_cleavages

    # If "all" is selected, process all proteases
    if args.protease.lower() == "all":
        selected_proteases = list(proteases.keys())  # All available proteases
    elif args.protease.lower() in proteases:
        selected_proteases = [args.protease.lower()]
    else:
        print(f"Protease {args.protease} is not supported. Supported proteases: {', '.join(proteases.keys())}")
        return

    for protease in selected_proteases:
        output_file = args.output or f"{base_filename}_predicted_{protease}_glycopeptides.csv"
        results = process_fasta(input_file, protease, missed_cleavages)
        write_csv(output_file, results)
        print(f"Results for {protease} written to {output_file}")

if __name__ == "__main__":
    main()