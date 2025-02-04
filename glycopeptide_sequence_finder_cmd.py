"""
glycopeptide_finder_cmd.py

This script processes a FASTA file to identify glycopeptides based on protease cleavage rules and glycosylation sequons, 
predicts their masses, calculates their hydrophobicity, and calculates their pI. The results are written to a CSV file.

Functions:
    cleave_sequence(sequence, protease, missed_cleavages):
        Cleaves a sequence based on protease rules.
    
    find_glycopeptides(peptides, full_sequence, glycosylation_type):
        Identifies peptides containing glycosylation sequons and maps the sites to the full protein sequence.
    
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
    python glycopeptide_finder_cmd.py -i <input_fasta_file> -o <output_csv_file> -p <protease> -g <glycosylation_type> -c <missed_cleavages> -l <log_file> -v

Author:
    Richard Shipman -- 2025
"""
import argparse
import csv
import re
from Bio import SeqIO
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from pyteomics.mass import calculate_mass
import os
import logging

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

# Define glycosylation rules
glycosylation = {
    "N": ("N[^P][STC]"),  # N-glycosylation sequon
    "O": ("[ST]"),  # O-glycosylation sequon (experimental! Creates large number of O-glycopeptides.)
    "C": ("W")  # C-glycosylation sequon - tryoptophan (W) only.
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

def find_glycopeptides(peptides, full_sequence, glycosylation_type):
    """Identifies peptides containing sequons and maps the sites to the full protein sequence."""
    glyco_sequon = re.compile(glycosylation[glycosylation_type])
    glycopeptides = []
    for pep in peptides:
        for match in glyco_sequon.finditer(pep):
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
    # Calculate the average hydrophobicity of the peptide
    total_hydrophobicity = sum(hydrophobicity_values.get(aa, 0) for aa in peptide_sequence)
    average_hydrophobicity = round(total_hydrophobicity / len(peptide_sequence), 2)

    return average_hydrophobicity

def calculate_pI(peptide_sequence):
    """Calculates the isoelectric point (pI) of a peptide sequence."""
    pI_calculator = IsoelectricPoint(peptide_sequence)
    return round(pI_calculator.pi(), 2)

def setup_logging(log_file):
    """Sets up logging to a file."""
    logging.basicConfig(filename=log_file, level=logging.INFO, 
                        format='%(asctime)s - %(levelname)s - %(message)s')

def process_fasta(file, protease, missed_cleavages, glycosylation_type):
    """Processes the input FASTA file and extracts glycopeptides."""
    results = []

    # Regular expression to capture OS and OX from the header
    os_ox_pattern = r"OS=([^\s]+(?: [^\s]+)*)\s+OX=(\d+)\s+GN=([^\s]+)\s+PE=(\d+)\s+SV=(\d+)"

    for record in SeqIO.parse(file, "fasta"):
        protein_id = record.id
        sequence = str(record.seq)

        header = record.description  # Get the full description line
        
        # Use regex to find OS, OX, GN, PE, and SV
        match = re.search(os_ox_pattern, header)
        if match:
            os_value = match.group(1)  # Organism Species
            ox_value = match.group(2)  # Organism Taxonomy ID
            gn_value = match.group(3)  # Gene Name
            pe_value = match.group(4)  # Protein Evidence
            sv_value = match.group(5)  # Sequence Version
        else:
            os_value = "Unknown"
            ox_value = "Unknown"
            gn_value = "Unknown"
            pe_value = "Unknown"
            sv_value = "Unknown"

        logging.info(f"Processing {protein_id} with {len(sequence)} amino acids.")        
        peptides = cleave_sequence(sequence, protease, missed_cleavages)
        logging.info(f"Found {len(peptides)} peptides after {protease} cleavage. The peptides were: {peptides}")
        x_glycopeptides = find_glycopeptides(peptides, sequence, glycosylation_type)
        logging.info(f"Found {len(x_glycopeptides)} {glycosylation_type}-glycopeptides. The glycopeptides were: {x_glycopeptides}")
        for peptide, site in x_glycopeptides:
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
                "Sequon": sequence[site - 1:site + 3], # Extract the sequon amino acid sequence + 1 flanking residue
                "PredictedMass": mass,
                "Hydrophobicity": hydrophobicity,
                "pI": pI,
                "Protease": protease,
                "GlycosylationType": glycosylation_type,
                "MissedCleavages": missed_cleavages,
                "Species": os_value,
                "TaxonID": ox_value,
                "GeneName": gn_value,
                "ProteinEvidence": pe_value,
                "SequenceVersion": sv_value
            })

    return results

# Add this function to write the results to a CSV file
def write_csv(output_file, data):
    """Writes results to a CSV file."""
    with open(output_file, mode="w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=["ProteinID", "Site", "Peptide", "Start", "End", "Length", "Sequon", "PredictedMass", "Hydrophobicity", "pI", "Protease", "GlycosylationType", "MissedCleavages", "Species", "TaxonID", "GeneName", "ProteinEvidence", "SequenceVersion"])
        writer.writeheader()
        writer.writerows(data)

# Add this to the main function to set up logging
def main():
    parser = argparse.ArgumentParser(description="Glycopeptide Finder")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file. Can be found in the test_proteomes folder.")
    parser.add_argument("-g", "--glycosylation", default="N", help="Glycosylation type (N, O, or C). Default is N. Large file sizes may result from selecting O or C.")
    parser.add_argument("-o", "--output", help="Output CSV file prefix. Default output directory for files is 'digested_glycopeptide_library'.")
    parser.add_argument("-p", "--protease", default="trypsin", help="Protease to use for cleavage ('all' for all proteases). Default is trypsin. Proteases: trypsin, chymotrypsin, glu-c, lys-c, arg-c, pepsin, asp-n, proteinase-k.")
    parser.add_argument("-c", "--missed_cleavages", type=int, default=0, help="Number of missed cleavages allowed. Default is 0.")
    parser.add_argument("-l", "--log", help="Provide log file name. (suggestion: -l log.txt)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose output.")

    # Parse arguments
    args = parser.parse_args()

    # Set up logging if log file is provided
    if args.log:
        setup_logging(args.log)

    # Process the input file
    input_file = args.input
    base_filename = input_file.rsplit(".", 1)[0]
    missed_cleavages = args.missed_cleavages
    glycosylation_type = args.glycosylation

    # If "all" is selected, process all proteases
    if args.protease.lower() == "all":
        selected_proteases = list(proteases.keys())  # All available proteases
    elif args.protease.lower() in proteases:
        selected_proteases = [args.protease.lower()] # Single protease
    else:
        if args.log: # Log the error
            logging.error(f"Protease {args.protease} is not supported. Supported proteases: {', '.join(proteases.keys())}")
        return

    # Ensure the output directory exists
    output_dir = "digested_glycopeptide_library"
    os.makedirs(output_dir, exist_ok=True)

    # Extract the base filename without any directory path
    base_filename = os.path.basename(base_filename)

    for protease in selected_proteases:
        output_file = args.output or f"{output_dir}/{base_filename}_predicted_{protease}_{glycosylation_type}-glycopeptides.csv"
        results = process_fasta(input_file, protease, missed_cleavages, glycosylation_type)
        write_csv(output_file, results)
        if args.log:
            logging.info(f"Results for {protease} written to {output_file}. Processing complete.")
        if args.verbose:
            print(f"Results for {protease} written to {output_file}. Processing complete.")


# main()
if __name__ == "__main__":
    main()