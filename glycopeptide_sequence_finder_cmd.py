"""
glycopeptide_sequence_finder_cmd.py

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
    
    process_fasta(file, protease, missed_cleavages, glycosylation_type):
        Processes the input FASTA file and extracts glycopeptides.
    
    write_csv(output_file, data):
        Writes results to a CSV file.
    
    compute_mz(mass, charge):
        Compute m/z value for a given mass and charge state.

    process_glycopeptides(peptide_file, glycan_file, max_charge):
        Generate glycopeptides and compute m/z values.
    
    setup_logging(log_file):
        Sets up logging to a file.
    
    main():
        Main function to parse arguments and execute the glycopeptide finding process.

Usage:
    python glycopeptide_sequence_finder_cmd.py -i <input_fasta_file> -o <output_csv_file> -p <protease> -g <glycosylation_type> -c <missed_cleavages> -l <log_file> -v -y <glycan_file> -z <max_charge>

Author:
    Richard Shipman -- 2025
"""
import argparse
import csv
import re
from Bio import SeqIO
import os
import logging
import pandas as pd

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

# Define glycosylation sequon rules
glycosylation = {
    "N": ("N[^P][STC]"),  # N-glycosylation sequon - N-X-S/T/C (X is any amino acid except P).
    "O": ("[ST]"),  # O-glycosylation sequon - S/T  (experimental! Creates large number of O-glycopeptides.)
    "C": ("W..[WCF]")  # C-glycosylation sequon - W-X-X-W, W-X-X-C, or W-X-X-F (X is any amino acid).
}

# Define amino acid mass values for calculating peptide mass
amino_acid_masses = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694, 'C': 103.00919,
    'Q': 128.05858, 'E': 129.04259, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
    'L': 113.08406, 'K': 128.09496, 'M': 131.04049, 'F': 147.06841, 'P': 97.05276,
    'S': 87.03203, 'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
}

# Kyte-Doolittle hydrophobicity values for amino acids
hydrophobicity_values = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# Define amino acid pKa values for calculating pI (placeholder)
pKa_values = {
    'C': 8.18, 'D': 3.65, 'E': 4.25, 'H': 6.00, 'K': 10.53,
    'R': 12.48, 'Y': 10.07, 'N': 3.22, 'Q': 3.22, 'S': 3.70,
    'T': 3.70, 'W': 10.07
}

# Define default N-glycan mass library as a DataFrame (N-Glycans) Using most common N-glycans as default, HexNAc(2)Hex(8) with no stereochemistry
default_n_glycan_library = pd.DataFrame([
    #{"glytoucan_ac": "G59324HL", "byonic": "HexNAc(2)Hex(12) % 2350.792627", "composition": "HexNAc(2)Hex(12)", "mass": 2350.792627, "shorthand_glycan": "N2H12"}, # N2H12
    #{"glytoucan_ac": "G58087IP", "byonic": "HexNAc(2)Hex(11) % 2188.739804", "composition": "HexNAc(2)Hex(11)", "mass": 2188.739804, "shorthand_glycan": "N2H11"}, # N2H11
    #{"glytoucan_ac": "G83460ZZ", "byonic": "HexNAc(2)Hex(10) % 2026.686980", "composition": "HexNAc(2)Hex(10)", "mass": 2026.686980, "shorthand_glycan": "N2H10"}, # N2H10
    #{"glytoucan_ac": "G80920RR", "byonic": "HexNAc(2)Hex(9) % 1864.634157", "composition": "HexNAc(2)Hex(9)", "mass": 1864.634157, "shorthand_glycan": "N2H9"}, # N2H9
    {"glytoucan_ac": "G62765YT", "byonic": "HexNAc(2)Hex(8) % 1702.581333", "composition": "HexNAc(2)Hex(8)", "mass": 1702.581333, "shorthand_glycan": "N2H8"}, # N2H8
    #{"glytoucan_ac": "G31852PQ", "byonic": "HexNAc(2)Hex(7) % 1540.528510", "composition": "HexNAc(2)Hex(7)", "mass": 1540.528510, "shorthand_glycan": "N2H7"}, # N2H7
    #{"glytoucan_ac": "G41247ZX", "byonic": "HexNAc(2)Hex(6) % 1378.475686", "composition": "HexNAc(2)Hex(6)", "mass": 1378.475686, "shorthand_glycan": "N2H6"}, # N2H6
])

# Define default O-glycan mass library as a DataFrame (O-GalNAc Glycans) Using most common O-glycans as default, HexNac(1) with no stereochemistry
default_o_glycan_library = pd.DataFrame([
    {"glytoucan_ac": "G14843DJ", "byonic": "HexNAc(1) % 221.089937305", "composition": "HexNAc(1)", "mass": 221.089937305, "shorthand_glycan": "N1"}, # N1
])

# Define default C-glycan mass library as a DataFrame (C-Mannosylation) Using most common C-glycans as default, Hex(1) with no stereochemistry
default_c_glycan_library = pd.DataFrame([
    {"glytoucan_ac": "G81399MY", "byonic": "Hex(1) % 180.0633882", "composition": "Hex(1)", "mass": 180.0633882, "shorthand_glycan": "H1"}, # H1
])

# Define monosaccharide mass library with multiple properties (later use)
monosaccharide_library = {
    "Hex": {"mass": 162.0528, "formula": "C6H10O5", "symbol": "H"},
    "HexA": {"mass": 176.0321, "formula": "C6H8O6", "symbol": "Ha"},
    "HexNAc": {"mass": 203.0794, "formula": "C8H13NO5", "symbol": "N"},
    "Fuc": {"mass": 146.0579, "formula": "C6H12O5", "symbol": "F"},
    "dHex": {"mass": 146.0579, "formula": "C6H12O5", "symbol": "dH"},
    "NeuAc": {"mass": 291.0954, "formula": "C11H17NO8", "symbol": "S"},
    "NeuGc": {"mass": 307.0903, "formula": "C11H17NO9", "symbol": "G"},
    "Xyl": {"mass": 150.0423, "formula": "C5H10O5", "symbol": "X"},
    "Pent": {"mass": 132.0423, "formula": "C5H10O5", "symbol": "P"},
    "Sulpho": {"mass": 79.9568, "formula": "SO3", "symbol": "Su"},
    "Phospho": {"mass": 79.9663, "formula": "PO3", "symbol": "Ph"},
    "Methyl": {"mass": 14.0157, "formula": "CH3", "symbol": "Me"},
    "Acetyl": {"mass": 42.0106, "formula": "C2H3O", "symbol": "Ac"},
    "Deoxy": {"mass": -18.0106, "formula": "H2O", "symbol": "d"},
    "Amine": {"mass": 1.0078, "formula": "H", "symbol": "NH2"}
}

# Define hydrophobicity values for glycan library (experimental) 
# (between -1 and 1) in relation to the peptide backbone and rank spread of glycopeptides with the same backdone.
# compute_glycan_hydrophobicity.py script will be used to calculate the hydrophobicity values.
# averaged aggregated from human and mouse datasets
default_glycan_hydrophobicity = {
    "N2H8": -0.0062899,
    "N2H7": -0.005723031,
    "N2H6": -0.004347122,
    "N2H9": -0.001982055,
    "N2H10": -0.001554309,
    "N2H11": -2.53305E-05,
    "N2H12": 5.36847E-06
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

def find_glycopeptides(peptides_df, glycosylation_type):
    """Identifies peptides containing sequons and maps the sites to the full protein sequence."""
    
    # Process pandas Dataframe
    peptides = [pep for sublist in peptides_df["Peptides"].tolist() for pep in sublist]
    full_sequence = peptides_df["Sequence"].tolist()[0]
    glycosylation_type = glycosylation_type
    
    glyco_sequon = re.compile(glycosylation[glycosylation_type])
    glycopeptides = []
    for pep in peptides:
        for match in glyco_sequon.finditer(pep):
            # Map the position in the peptide to the full protein sequence
            start_in_protein = full_sequence.find(pep) + match.start() + 1  # Adjust to 1-based indexing
            glycopeptides.append((pep, start_in_protein))
    return glycopeptides

def calculate_peptide_mass(sequence):
    """Calculates the mass of a peptide using predefined amino acid masses."""
    
    # Common ambiguous residues
    invalid_residues = {"X", "B", "Z", "J", "U", "O"}  
    if any(aa in invalid_residues for aa in sequence):
        return "Unknown"  # Or return unknown if you prefer
    
    # water
    water  = 18.010565

    # Calculate the mass using the amino_acid_masses dictionary
    mass = sum(amino_acid_masses.get(aa, 0) for aa in sequence) + water
    return mass

def predict_hydrophobicity(peptide_sequence):
    """Predicts the hydrophobicity of a peptide sequence using the Kyte-Doolittle scale."""
    
    if len(peptide_sequence) == 0:
        return 0.0

    # Calculate the average hydrophobicity of the peptide
    total_hydrophobicity = sum(hydrophobicity_values.get(aa, 0) for aa in peptide_sequence)
    average_hydrophobicity = round(total_hydrophobicity / len(peptide_sequence), 5)

    return average_hydrophobicity

def calculate_pI(peptide_sequence):
    """
    Calculate the isoelectric point (pI) of a peptide sequence.
    
    Parameters:
        peptide (str): Peptide sequence (1-letter amino acid codes)
    
    Returns:
        float: Estimated isoelectric point (pI)
    """
    # peptide sequence to string
    peptide_sequence = str(peptide_sequence)

    # Terminal group pKa values (assuming N-term = 9.6, C-term = 2.3)
    N_term_pKa = 9.6
    C_term_pKa = 2.3

    # Count occurrences of ionizable residues
    residue_counts = {aa: peptide_sequence.count(aa) for aa in pKa_values}
    
    # Function to calculate charge at a given pH
    def calculate_net_charge(pH):
        charge = 0.0

        # N-terminal charge
        charge += 1 / (1 + 10**(pH - N_term_pKa))

        # C-terminal charge
        charge -= 1 / (1 + 10**(C_term_pKa - pH))

        # Side chain charges
        for aa, count in residue_counts.items():
            if count > 0:
                pKa = pKa_values[aa]
                if aa in ['D', 'E', 'Y', 'C']:  # Acidic side chains
                    charge -= count / (1 + 10**(pKa - pH))
                elif aa in ['H', 'K', 'R']:  # Basic side chains
                    charge += count / (1 + 10**(pH - pKa))
        
        return charge

    # Use bisection method to find the pH where net charge is closest to zero
    low, high = 0.0, 14.0
    while high - low > 0.01:  # Precision threshold
        mid = (low + high) / 2
        net_charge = calculate_net_charge(mid)
        if net_charge > 0:
            low = mid
        else:
            high = mid

    return round((low + high) / 2, 2)

def compute_mz(mass, charge):
    """Compute m/z value for a given mass and charge state."""
    proton_mass = 1.007276
    return (mass + (charge * proton_mass)) / charge

# experimental glycopeptide hydrophobicity calculation
def compute_hf_experimental(peptide_hydrophobicity, glycan, hf_weight=10, rt_scale=60):
    """Compute HF_experimental and scaled retention time (rt_HF_experimental) values."""
    glycan_hydrophobicity = default_glycan_hydrophobicity.get(glycan, None)
    
    # Compute HF_experimental and scaled retention time (rt_HF_experimental) values
    if glycan_hydrophobicity is not None:
        hf_experimental = peptide_hydrophobicity + glycan_hydrophobicity * hf_weight
        rt_hf_experimental = (glycan_hydrophobicity * hf_weight + 1) * (rt_scale / 2)
    else:
        hf_experimental = ""
        rt_hf_experimental = ""
    
    return round(hf_experimental, 5), rt_hf_experimental

def process_glycopeptides(peptide_file, glycans, max_charge):
    """Generate glycopeptides and compute m/z values."""
    # Load peptide and glycan data
    peptides = pd.read_csv(peptide_file, low_memory=False)
    
    # Convert relevant columns to float
    peptides = peptides[pd.to_numeric(peptides['PredictedMass'], errors='coerce').notnull()]
    peptides['PredictedMass'] = peptides['PredictedMass'].astype(float)
    glycans['mass'] = glycans['mass'].astype(float)
    
    # Initialize list to store results
    results = []
    
    # Generate glycopeptides, compute m/z values, compute HF_experimental
    for _, pep in peptides.iterrows():
        for _, gly in glycans.iterrows():

            # Compute glycopeptide mass
            glycopeptide_mass = pep['PredictedMass'] + gly['mass']

            # Compute m/z values for charge states from 2 to max_charge
            mz_values = {f'z{z}': compute_mz(glycopeptide_mass, z) for z in range(2, max_charge + 1)}

            # Create a dictionary with the results
            result = {
                'ProteinID': pep['ProteinID'],
                'Site': pep['Site'],
                'GlyToucan_AC': gly['glytoucan_ac'],
                'Composition': gly['composition'],
                #'converted_glycan': gly['converted_glycan'],
                'Peptide': pep['Peptide'],
                'Start': pep['Start'],
                'End': pep['End'],
                'Length': pep['Length'],
                'Sequon': pep['Sequon'],
                'GlycopeptideMass': glycopeptide_mass,
                'PeptideMass': pep['PredictedMass'],
                'GlycanMass': gly['mass'],
                'Hydrophobicity': pep['Hydrophobicity'],
                'pI': pep['pI'],
                'Protease': pep['Protease'],
                'GlycosylationType': pep['GlycosylationType'],
                'MissedCleavages': pep['MissedCleavages'],
                #'Species': pep['Species'],
                #'TaxonID': pep['TaxonID'],
                'GeneName': pep['GeneName'],
                #'ProteinEvidence': pep['ProteinEvidence'],
                #'SequenceVersion': pep['SequenceVersion'],
                **mz_values, # z charge states values
            }

            results.append(result)
    
    return pd.DataFrame(results)

def setup_logging(log_file):
    """Sets up logging to a file."""
    logging.basicConfig(filename=log_file, level=logging.INFO, 
                        format='%(asctime)s - %(levelname)s - %(message)s')

def process_fasta(file, protease, missed_cleavages, glycosylation_type):
    """Processes the input FASTA file and extracts glycopeptides."""
    
    # Initialize list to store results
    results = []

    # Regular expression to capture OS and OX from the header
    os_ox_pattern = r"OS=([^\s]+(?: [^\s]+)*)\s+OX=(\d+)\s+GN=([^\s]+)\s+PE=(\d+)\s+SV=(\d+)"

    # Process each record in the FASTA file
    for record in SeqIO.parse(file, "fasta"):
        protein_id = record.id
        sequence = str(record.seq)

        # Get the full description line
        header = record.description  
        
        # Use regex to find OS, OX, GN, PE, and SV
        match = re.search(os_ox_pattern, header)
        if match:
            os_value = match.group(1)  # Organism Species
            ox_value = match.group(2)  # Organism Taxonomy ID
            gn_value = match.group(3)  # Gene Name
            pe_value = match.group(4)  # Protein Evidence
            sv_value = match.group(5)  # Sequence Version
        else:
            # blank values if not found
            os_value = ""
            ox_value = ""
            gn_value = ""
            pe_value = ""
            sv_value = ""

        # Log the digestion processing of the protein with protease and missed cleavages
        logging.info(f"Processing {protein_id} with {len(sequence)} amino acids.")        
        peptides = cleave_sequence(sequence, protease, missed_cleavages)
        logging.info(f"Found {len(peptides)} peptides after {protease} cleavage. The peptides were: {peptides}")
        
        # Write protein, peptide and data above list of string results to a pandas DataFrame
        results.append({
            "ProteinID": protein_id,
            "Peptides": peptides,
            "Sequence": sequence,
            "Protease": protease,
            "MissedCleavages": missed_cleavages,
            "GlycosylationType": glycosylation_type,
            "Species": os_value,
            "TaxonID": ox_value,
            "GeneName": gn_value,
            "ProteinEvidence": pe_value,
            "SequenceVersion": sv_value
        })

    return pd.DataFrame(results)

# Function to Calculate glycopeptide ion series m/z values
def calculate_n_glycopeptide_ions(peptide, glycan_composition, glycan_frag_order=None, charge=1):
    """
    Calculate theoretical m/z values for b, y, Y, B, oxonium ions, and sugar chunks of a glycopeptide.
    
    Parameters:
      peptide (str): The peptide sequence (e.g., "NTSK").
      glycan_composition (str): A string representing the glycan composition, e.g., "HexNAc(2)Hex(8)".
      glycan_frag_order (list, optional): A list specifying the order of sugar losses.
                                          If provided, it is used to simulate sequential 
                                          loss for Y and B ion series and to generate sugar chunks.
      charge (int): The charge state (default is 1).
      
    Returns:
      dict: A dictionary with keys 'b', 'y', 'Y', 'B', 'oxonium', and 'chunks' containing
            the calculated m/z values.
    """
    # Constants (monoisotopic masses in Da)
    proton = 1.007276
    water  = 18.010565

    # Compute the neutral mass of the peptide (include water for the full peptide mass)
    peptide_mass = sum(amino_acid_masses[aa] for aa in peptide) + water

    # Parse glycan composition string into a dictionary
    glycan_dict = {}
    for part in glycan_composition.split(')'):
        if part:
            sugar, count = part.split('(')
            glycan_dict[sugar] = int(count)

    # --- Calculate b ions ---
    b_ions = []
    cumulative = 0.0
    # b ions: cumulative from the N-terminus (excluding the full peptide)
    for i in range(len(peptide) - 1):
        cumulative += amino_acid_masses[peptide[i]]
        b_val = cumulative + proton  # b ion mass = sum + proton
        b_ions.append(round(b_val / charge, 4))
    
    # --- Calculate y ions ---
    y_ions = []
    cumulative = 0.0
    # y ions: cumulative from the C-terminus (excluding the full peptide)
    for i in range(len(peptide) - 1, 0, -1):
        cumulative += amino_acid_masses[peptide[i]]
        y_val = cumulative + water + proton  # y ion mass = sum + water + proton
        y_ions.insert(0, round(y_val / charge, 4))
    
    # --- Glycan Calculations ---
    glycan_total_mass = sum(monosaccharide_library[sugar]['mass'] * count for sugar, count in glycan_dict.items())
    intact_glycopeptide_mass = peptide_mass + glycan_total_mass

    # --- Calculate Y ions (peptide + residual glycan fragment) ---
    Y_ions = {}
    if glycan_frag_order:
        current_mass = intact_glycopeptide_mass
        Y_ions['Y0'] = round((current_mass + proton) / charge, 4)
        for i, sugar in enumerate(glycan_frag_order, start=1):
            current_mass -= monosaccharide_library[sugar]['mass']
            Y_ions[f'Y{i}'] = round((current_mass + proton) / charge, 4)
    else:
        # Default approach for many N-glycopeptides:
        Y_ions['Y0'] = round((peptide_mass + proton) / charge, 4)
        if glycan_dict.get('HexNAc', 0) > 0:
            Y_ions['Y1'] = round((peptide_mass + monosaccharide_library['HexNAc']['mass'] + proton) / charge, 4)
        if glycan_dict.get('HexNAc', 0) > 1:
            Y_ions['Y2'] = round((peptide_mass + 2 * monosaccharide_library['HexNAc']['mass'] + proton) / charge, 4)
    
    # # --- Calculate B ions (glycan fragment ions) ---
    # B_ions = {}
    # if glycan_frag_order:
    #     cumulative = 0.0
    #     for i, sugar in enumerate(glycan_frag_order, start=1):
    #         cumulative += monosaccharide_library[sugar]['mass']
    #         B_ions[f'B{i}'] = round((cumulative + proton) / charge, 4)
    # else:
    #     # Report common oxonium ions as B ions if no fragmentation order is given.
    #     if glycan_dict.get('HexNAc', 0) > 0:
    #         B_ions['B_HexNAc'] = round((monosaccharide_library['HexNAc']['mass'] + proton) / charge, 4)
    #     if glycan_dict.get('Hex', 0) > 0:
    #         B_ions['B_Hex'] = round((monosaccharide_library['Hex']['mass'] + proton) / charge, 4)
    
    # # --- Calculate Oxonium ions ---
    # # These are common low–m/z ions originating from individual sugar units.
    # oxonium_ions = {}
    # for sugar, count in glycan_dict.items():
    #     if count > 0:
    #         oxonium_ions[f'ox_{sugar}'] = round(monosaccharide_library[sugar]['mass'] + proton, 4)
    
    # --- Calculate Sugar Chunks (contiguous groups of sugars) ---
    # This simulates the possibility that multiple sugars break off together.
    # sugar_chunks = {}
    # if glycan_frag_order:
    #     n = len(glycan_frag_order)
    #     # Enumerate all contiguous subsequences (chunks)
    #     for i in range(n):
    #         for j in range(i, n):
    #             chunk_label = f"chunk_{i+1}-{j+1}"
    #             chunk_mass = sum(monosaccharide_library[sugar]['mass'] for sugar in glycan_frag_order[i:j+1])
    #             # Adding a proton to simulate a charged fragment
    #             sugar_chunks[chunk_label] = round((chunk_mass + proton) / charge, 4)

    return {
        'b': b_ions,
        'y': y_ions,
        'Y': Y_ions,
        #'B': B_ions,
        #'oxonium': oxonium_ions,
        #'chunks': sugar_chunks
    }

# Add this function to write the results to a CSV file
def write_csv(output_file, data):
    """Writes results to a CSV file."""
    with open(output_file, mode="w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=["ProteinID", "Site", "Peptide", "Start", "End", "Length", "Sequon", "PredictedMass", "Hydrophobicity", "pI", "Protease", "GlycosylationType", "MissedCleavages", "Species", "TaxonID", "GeneName", "ProteinEvidence", "SequenceVersion"])
        writer.writeheader()
        writer.writerows(data)

# Add this to the main function to set up logging
def main():
    
    # SETUP ARGUMENT PARSER
    
    parser = argparse.ArgumentParser(description="Glycopeptide Finder")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file. Can be found in the test_proteomes folder.")
    parser.add_argument("-g", "--glycosylation", default="N", help="Glycosylation type (N, O, or C). Default is N. Large file sizes may result from selecting O or C.")
    parser.add_argument("-o", "--output", help="Output CSV file prefix. Default output directory for files is 'digested_glycopeptide_library'.")
    parser.add_argument("-p", "--protease", default="trypsin", help="Protease to use for cleavage ('all' for all proteases). Default is trypsin. Proteases: trypsin, chymotrypsin, glu-c, lys-c, arg-c, pepsin, asp-n, proteinase-k.")
    parser.add_argument("-c", "--missed_cleavages", type=int, default=0, help="Number of missed cleavages allowed. Default is 0.")
    parser.add_argument("-y", "--glycan", default=None, help="Path to glycan file (CSV). Default is 'default_glycan_library.csv'.")
    parser.add_argument("-l", "--log", help="Provide log file name. (suggestion: -l log.txt)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose output.")
    parser.add_argument("-z", "--charge", type=int, default=5, help="Maximum charge state (default: 5).")

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
    charge_state = args.charge
    glycan_library = args.glycan

    # Set default glycan library based on glycosylation type only if no glycan library is provided
    if glycan_library is None:
        if glycosylation_type == "N":
            default_glycan_library = default_n_glycan_library
        elif glycosylation_type == "O":
            default_glycan_library = default_o_glycan_library
        elif glycosylation_type == "C":
            default_glycan_library = default_c_glycan_library
        else:
            if args.log:
                logging.error(f"Glycosylation type {glycosylation_type} is not supported. Supported glycosylation types: N, O, C.")
            return
    else:
        default_glycan_library = None

    # Load glycan library
    if glycan_library is None:
        glycans = default_glycan_library
    else:
        glycans = pd.read_csv(glycan_library)

    # Ensure the output directory exists (digested_glycopeptide_library)
    output_dir = "digested_glycopeptide_library"
    os.makedirs(output_dir, exist_ok=True)

    # Ensure the output directory exists (digested_peptide_library)
    peptide_output_dir = "digested_peptide_library"
    os.makedirs(peptide_output_dir, exist_ok=True)

    # Extract the base filename without any directory path
    base_filename = os.path.basename(base_filename)

    # Process the input file with the selected protease(s)
    all_results = []

    # If "all" is selected, process all proteases
    if args.protease.lower() == "all":
        selected_proteases = list(proteases.keys())  # All available proteases
    elif args.protease.lower() in proteases:
        selected_proteases = [args.protease.lower()] # Single protease
    else:
        if args.log: # Log the error
            logging.error(f"Protease {args.protease} is not supported. Supported proteases: {', '.join(proteases.keys())}")
        return

    # WORKFLOW STARTS HERE

    # Process the input file with the selected protease(s)
    for protease in selected_proteases:
        
        # Log the start of the process
        print(f"Processing {input_file} with protease {protease} and {missed_cleavages} missed cleavages...")
        if args.log:
            logging.info(f"Processing {input_file} with protease {protease} and {missed_cleavages} missed cleavages...")
        
        # Process the input file and return pandas DataFrame
        results_df = process_fasta(input_file, protease, missed_cleavages, glycosylation_type)
        #print(results_df)

        # Flatten results_df to extract peptides
        digest_peptide_library = pd.DataFrame([pep for sublist in results_df["Peptides"].tolist() for pep in sublist])
        
        # Set name of first column to "Peptide"
        digest_peptide_library.columns = ["Peptide"]

        # Add other columns to the DataFrame
        digest_peptide_library["ProteinID"] = results_df["ProteinID"].tolist()[0]
        #digest_peptide_library["Protease"] = protease
        #digest_peptide_library["MissedCleavages"] = missed_cleavages
        #digest_peptide_library["Species"] = results_df["Species"].tolist()[0]
        #digest_peptide_library["TaxonID"] = results_df["TaxonID"].tolist()[0]
        digest_peptide_library["GeneName"] = results_df["GeneName"].tolist()[0]
        #digest_peptide_library["ProteinEvidence"] = results_df["ProteinEvidence"].tolist()[0]
        #digest_peptide_library["SequenceVersion"] = results_df["SequenceVersion"].tolist()[0]

        # Remove duplicate peptide entries in digest_peptide_library
        digest_peptide_library = digest_peptide_library.drop_duplicates(subset=["Peptide"])

        # Remove empty peptide entries in digest_peptide_library
        digest_peptide_library = digest_peptide_library[digest_peptide_library["Peptide"].str.strip().astype(bool)]

        # Compute peptide mass, hydrophobicity, pI, and m/z values per row
        digest_peptide_library["PredictedMass"] = digest_peptide_library["Peptide"].apply(calculate_peptide_mass)
        digest_peptide_library["Hydrophobicity"] = digest_peptide_library["Peptide"].apply(predict_hydrophobicity)
        digest_peptide_library["pI"] = digest_peptide_library["Peptide"].apply(calculate_pI)
        digest_peptide_library["PredictedMass"] = pd.to_numeric(digest_peptide_library["PredictedMass"], errors='coerce')
        digest_peptide_library = digest_peptide_library.dropna(subset=["PredictedMass"])
        digest_peptide_library["z1"] = digest_peptide_library["PredictedMass"].apply(lambda mass: compute_mz(mass, 1)) # charge state = 1

        # Write results_df to a CSV file in peptide_library folder
        digest_peptide_library.to_csv(f"{peptide_output_dir}/{base_filename}_{protease}_digested_mc{missed_cleavages}_peptides.csv", index=False)

        # Extract the protein sequence
        all_glycopeptides = []
        for _, row in results_df.iterrows():
            protein_id = row["ProteinID"]
            sequence = row["Sequence"]

            # Find glycopeptides
            x_glycopeptides = find_glycopeptides(pd.DataFrame([row]), glycosylation_type)
            
            # Process each glycopeptide
            for peptide, site in x_glycopeptides:
                mass = calculate_peptide_mass(peptide)
                hydrophobicity = predict_hydrophobicity(peptide)
                pI = calculate_pI(peptide)
                start_pos = sequence.find(peptide) + 1  # 1-based indexing
                end_pos = start_pos + len(peptide) - 1
                all_glycopeptides.append({
                    "ProteinID": protein_id,
                    "Site": site,
                    "Peptide": peptide,
                    "Start": start_pos,
                    "End": end_pos,
                    "Length": len(peptide),
                    "Sequon": sequence[site - 1:site + 2], # Extract the sequon amino acid sequence + 1 flanking residue
                    "PredictedMass": mass,
                    "Hydrophobicity": hydrophobicity,
                    "pI": pI,
                    #"Protease": protease,
                    #"GlycosylationType": glycosylation_type,
                    #"MissedCleavages": missed_cleavages,
                    #"Species": row["Species"], # Uncomment to include 
                    #"TaxonID": row["TaxonID"],
                    "GeneName": row["GeneName"],
                    #"ProteinEvidence": row["ProteinEvidence"],
                    #"SequenceVersion": row["SequenceVersion"]
                })

        # Concatenate the results to the main DataFrame
        results_df = pd.concat([results_df, pd.DataFrame(all_glycopeptides)], ignore_index=True)

        # Log the completion of the process
        all_results.extend(results_df.to_dict(orient='records'))
        if args.verbose:
            print(f"Processed {len(results_df)} peptides for protease {protease}.")
        
        # Write the results to a new CSV file
        output_file = args.output or f"{output_dir}/{base_filename}_{protease}_digested_mc{missed_cleavages}_z{charge_state}_{glycosylation_type}-glycopeptides.csv"
        
        # Filter out unwanted fields before writing to CSV
        filtered_results = [{k: v for k, v in result.items() if k in ["ProteinID", "Site", "Peptide", "Start", "End", "Length", "Sequon", "PredictedMass", "Hydrophobicity", "pI", "Protease", "GlycosylationType", "MissedCleavages", "Species", "TaxonID", "GeneName", "ProteinEvidence", "SequenceVersion"]} for result in all_results]
        write_csv(output_file, filtered_results)
        #print(f"Results written to {output_file}. Processing complete.")

        # Log the completion of the process
        if args.log:
            logging.info(f"Results written to {output_file}. Processing complete.")
        if args.verbose:
            print(f"Results written to {output_file}. Processing complete.")

        # Set the output file for glycopeptides
        glycopeptide_output_file = output_file
        #print(f"Generating glycopeptides and computing m/z values...")

        # Log the start of the glycopeptide processing
        if args.log:
            logging.info(f"Generating glycopeptides and computing m/z values with range of +2 to +{args.charge} charge states using {args.glycan} glycan library.")
        if args.verbose:
            print(f"Generating glycopeptides and computing m/z values with range of +2 to +{args.charge} charge states using {args.glycan} glycan library.")

        # Process glycopeptides and compute m/z values
        glycopeptide_results = process_glycopeptides(output_file, glycans, charge_state)

        # Compute Ion Series
        glycopeptide_results["IonSeries"] = glycopeptide_results.apply(lambda row: calculate_n_glycopeptide_ions(row["Peptide"], row["Composition"]), axis=1)

        # Write the results to a new CSV file
        glycopeptide_results.to_csv(glycopeptide_output_file, index=False)
        #print(glycopeptide_results)
        #print(f"Glycopeptide results written to {glycopeptide_output_file}. Processing complete.")

        # Log the completion of the glycopeptide processing
        if args.log:
            logging.info(f"Glycopeptide results written to {glycopeptide_output_file}. Processing complete.")
        if args.verbose:
            print(f"Glycopeptide results written to {glycopeptide_output_file}. Processing complete.")

# main function
if __name__ == "__main__":
    main()
