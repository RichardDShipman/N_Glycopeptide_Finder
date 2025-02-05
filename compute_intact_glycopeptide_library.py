import argparse
import pandas as pd
import os

def compute_mz(mass, charge):
    """Compute m/z value for a given mass and charge state."""
    proton_mass = 1.007276
    return (mass + (charge * proton_mass)) / charge

def process_glycopeptides(peptide_file, glycan_file, max_charge):
    """Generate glycopeptides and compute m/z values."""
    # Load peptide and glycan data
    peptides = pd.read_csv(peptide_file, low_memory=False)
    glycans = pd.read_csv(glycan_file)
    
    # Convert relevant columns to float
    peptides = peptides[pd.to_numeric(peptides['PredictedMass'], errors='coerce').notnull()]
    peptides['PredictedMass'] = peptides['PredictedMass'].astype(float)
    glycans['mass'] = glycans['mass'].astype(float)
    
    results = []
    
    for _, pep in peptides.iterrows():
        for _, gly in glycans.iterrows():
            glycopeptide_mass = pep['PredictedMass'] + gly['mass']
            mz_values = {f'z{z}': compute_mz(glycopeptide_mass, z) for z in range(2, max_charge + 1)}
            
            result = {
                'ProteinID': pep['ProteinID'],
                'Site': pep['Site'],
                'glytoucan_ac': gly['glytoucan_ac'],
                'composition': gly['composition'],
                'Peptide': pep['Peptide'],
                'Start': pep['Start'],
                'End': pep['End'],
                'Length': pep['Length'],
                'Sequon': pep['Sequon'],
                'GlycopeptideMass': glycopeptide_mass,
                'PeptideMass': pep['PredictedMass'],
                'glycan_mass': gly['mass'],
                'Hydrophobicity': pep['Hydrophobicity'],
                'pI': pep['pI'],
                'Protease': pep['Protease'],
                'GlycosylationType': pep['GlycosylationType'],
                'MissedCleavages': pep['MissedCleavages'],
                'Species': pep['Species'],
                'TaxonID': pep['TaxonID'],
                'GeneName': pep['GeneName'],
                'ProteinEvidence': pep['ProteinEvidence'],
                'SequenceVersion': pep['SequenceVersion'],
                **mz_values
            }
            results.append(result)
    
    return pd.DataFrame(results)

def main():
    '''Main function to handle command-line arguments and call the process_glycopeptides function.'''
    parser = argparse.ArgumentParser(description='Compute glycopeptide mass and m/z values.')
    parser.add_argument('-i', '--input', required=True, help='Path to peptide file (CSV).')
    parser.add_argument('-g', '--glycan', required=True, help='Path to glycan file (CSV).')
    parser.add_argument('-o', '--output', help='Output file path.')
    parser.add_argument('-z', '--charge', type=int, default=6, help='Maximum charge state (default: 6).')
    
    args = parser.parse_args()
    
    input_filename = os.path.splitext(os.path.basename(args.input))[0]
    output_dir = 'predicted_intact_glycopeptide_library'
    os.makedirs(output_dir, exist_ok=True)
    
    if args.output:
        output_filename = args.output
    else:
        if input_filename.endswith('glycopeptides'):
            input_filename = input_filename[:-13]  # Remove 'glycopeptides'
        output_filename = os.path.join(output_dir, f'{input_filename}glycopeptide_library.csv')
    
    output_df = process_glycopeptides(args.input, args.glycan, args.charge)
    output_df.to_csv(output_filename, index=False)
    print(f'Results saved to {output_filename}')

if __name__ == '__main__':
    main()
