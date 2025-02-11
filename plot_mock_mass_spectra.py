import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import os
import argparse
import ast  # To safely parse the IonSeries string into a dictionary

def read_input_data(csv_file):
  """"""
  return pd.read_csv(csv_file)

def create_output_directory(output_dir):
  """Create output directory function"""
  if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def define_colors(df):
    """Define a color mapping function for ion types in mass spectrometry data."""
    
    # Predefined color mapping for common ion types
    color_mapping = {
        "b": "blue",      # HCD
        "y": "green",     # HCD
        "c": "gold",      # ETD
        "z": "hotpink",   # ETD
        "Y": "red",       # Large fragment
        "B": "darkorange",# Large fragment
        "oxonium": "purple" # Small diagnostic ions
    }
    
    # Convert IonSeries from string to dictionary if necessary
    df['IonSeries'] = df['IonSeries'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

    # Identify unique ion types present in the dataset
    unique_ions = {ion for ion_series in df['IonSeries'] for ion in ion_series.keys()}

    # Additional colors for unknown ion types
    additional_colors = [
        "cyan", "magenta", "lime", "brown", "gray", "teal"
    ]

    # Assign colors to any new ion types
    for ion in unique_ions:
        if ion not in color_mapping:
            color_mapping[ion] = additional_colors.pop(0) if additional_colors else "black"

    return color_mapping

def compute_ion_number(row):
    """Compute Ion numbers for each ion series and return in a formatted string"""
    ion_series = row['IonSeries']
    ion_numbers = []

    # Process b ions
    if 'b' in ion_series:
        b_ions = [f"b{i+1}" for i in range(len(ion_series['b']))]
        ion_numbers.extend(b_ions)
    
    # Process y ions (order reversed)
    if 'y' in ion_series:
        y_ions = [f"y{i+1}" for i in range(len(ion_series['y']))]
        ion_numbers.extend(y_ions[::-1])  # Reverse the y ions
    
    # Process c ions
    if 'c' in ion_series:
        c_ions = [f"c{i+1}" for i in range(len(ion_series['c']))]
        ion_numbers.extend(c_ions)
    
    # Process z ions (order reversed)
    if 'z' in ion_series:
        z_ions = [f"z{i+1}" for i in range(len(ion_series['z']))]
        ion_numbers.extend(z_ions[::-1])
    
    # Process Y ions
    if 'Y' in ion_series:
        Y_ions = [f"Y{i}" for i in range(len(ion_series['Y']))]
        ion_numbers.extend(Y_ions)
    
    # Process B ions
    if 'B' in ion_series:
        B_ions = [f"B{i}" for i in range(len(ion_series['B']))]
        ion_numbers.extend(B_ions)
    
    # Process Oxonium ions
    if 'oxonium' in ion_series:
        oxonium_ions = [f"ox{i}" for i in range(len(ion_series['oxonium']))]
        ion_numbers.extend(oxonium_ions)

    return ', '.join(ion_numbers)

def create_mass_spectrum_plot(df, color_mapping, title):
    """Plot mass spectrum with ion numbers as labels and variable ion line length, distinguishing ETD and HCD ions."""
    fig, ax = plt.subplots(figsize=(12, 8))

    # Define custom line lengths for each ion type
    ion_length_map = {
        'b': 10,  # b ions (HCD)
        'y': 20,  # y ions (HCD)
        'c': 30,  # c ions (ETD)
        'z': 40,   # z ions (ETD)
        'Y': 60,  # Y ions
        '2Y':65, # Y ions 2+ charge state
        'B': 50,  # B ions
        'oxonium': 70,  # Oxonium ions
    }

    # Define line styles to differentiate ion types
    ion_style_map = {
        'b': '-',  # Solid line for HCD
        'y': '-',  # Solid line for HCD
        'c': '--', # Dashed line for ETD
        'z': '--', # Dashed line for ETD
        'Y': ':',  # Dotted line
        '2Y': ':', # Plus for 2Y+ plus charge state
        'B': '-.', # Dash-dot line
        'oxonium': ':' # Dotted line
    }

    # Store legend handles
    legend_handles = {}

    # Iterate over the DataFrame
    for _, row in df.iterrows():
        ion_series = row['IonSeries']

        for ion_type, values in ion_series.items():
            if isinstance(values, list):  # For b, y, c, z ions
                total_length = len(values)

                for i, mz in enumerate(values, start=0):
                    if ion_type in ['b', 'c']:  # Normal numbering
                        ion_number = i + 1
                    elif ion_type in ['y', 'z']:  # Reverse numbering
                        ion_number = total_length - i

                    ion_label = f"{ion_type}{ion_number}"
                    line_length = ion_length_map.get(ion_type, 50)
                    color = color_mapping.get(ion_type, "black")
                    linestyle = ion_style_map.get(ion_type, "-")  # Default to solid

                    ax.vlines(mz, 0, line_length, color=color, lw=2, linestyle=linestyle)

                    # Add text label with white outline
                    text = ax.text(mz, line_length, ion_label, rotation=90, verticalalignment='bottom', fontsize=8, color='black')
                    text.set_path_effects([path_effects.Stroke(linewidth=3, foreground="white"), path_effects.Normal()])

                    # Add to legend if not already present
                    if ion_type not in legend_handles:
                        label = f"{ion_type} ({'HCD' if ion_type in ['b', 'y'] else 'ETD'})"
                        legend_handles[ion_type] = plt.Line2D([0], [0], color=color, lw=2, linestyle=linestyle, label=label)

            elif isinstance(values, dict):  # For Y, B, oxonium ions
                for ion_name, mz in values.items():
                    line_length = ion_length_map.get(ion_type, 50)
                    color = color_mapping.get(ion_type, "black")
                    linestyle = ion_style_map.get(ion_type, "-")

                    ax.vlines(mz, 0, line_length, color=color, lw=2, linestyle=linestyle)

                    # Add text label with white outline
                    text = ax.text(mz, line_length, ion_name, rotation=90, verticalalignment='bottom', fontsize=8, color='black')
                    text.set_path_effects([path_effects.Stroke(linewidth=3, foreground="white"), path_effects.Normal()])

                    # Add to legend if not already present
                    if ion_type not in legend_handles:
                        legend_handles[ion_type] = plt.Line2D([0], [0], color=color, lw=2, linestyle=linestyle, label=ion_type)

    # Add the legend
    ax.legend(handles=legend_handles.values(), title="Ion Types", loc="upper right", fontsize=8)

    ax.set_xlabel("Calculated m/z values (Da)")
    ax.set_ylabel("Intensity (%) (Mock values)")
    ax.set_title(title)
    ax.set_ylim(0, 100)  # Adjust based on the plot
    plt.tight_layout()
    return fig, ax

def save_plot(fig, output_file):
  """Save plot to output location"""
  fig.savefig(output_file)
  plt.close(fig)
  print(f"Plot saved to {output_file}")

def main():
  """Main function"""
  parser = argparse.ArgumentParser(description="Plot mock mass spectra from glycopeptide ion series CSV file.")
  parser.add_argument('-i', '--input', required=True, help="Input CSV file containing glycopeptide ion series info.")
  parser.add_argument('-o', '--output', default="mock_mass_spectra", help="Output directory to save the plot.")
  args = parser.parse_args()

  # Arguement parser
  csv_file = args.input
  output_dir = args.output

  # Read csv to dataframe from input -i flag
  df = read_input_data(csv_file)

  # Review file for input data
  required_columns = ['IonSeries', 'ProteinID', 'Peptide', 'Composition', 'GlyToucan_AC']
  print("Columns in CSV file:", df.columns.tolist())
  for col in required_columns:
    if col not in df.columns:
      raise KeyError(f"Missing required column: {col}")

  # Add ion_number column based on IonSeries
  # Ensure IonSeries is parsed as dictionary
  df['IonSeries'] = df['IonSeries'].apply(ast.literal_eval)
  df['ion_number'] = df.apply(compute_ion_number, axis=1)

  #print(df)

  # Creat output directory
  create_output_directory(output_dir)
  
  # Define column mapping
  color_mapping = define_colors(df)

  # Iterate through each row in the dataframe
  for index, row in df.iterrows():
    # Extract metadata for plot title and filename
    protein = row['ProteinID']
    peptide_sequence = row['Peptide']
    glytoucan_ac = row['GlyToucan_AC']
    glycan_composition = row['Composition']
    glycopeptide_mass = row ['GlycopeptideMass']
    glycan_mass = row['GlycanMass']
    peptide_mass = row['PeptideMass']
    site = pd.to_numeric(row['Site'])
    if site.is_integer():
      site = int(site)
    
    # Create descriptive title with metadata
    title = f"Glycopeptide Sequence Finder: Mock Glycopeptide Mass Spectrum (HCD: b, y ions / ETD: c, z ions)\nProtein Identifer: {protein}, Glycopeptide Mass: {glycopeptide_mass}\nSite: {site}, Peptide: {peptide_sequence}, Peptide Mass: {peptide_mass}\nGlycan Composition: {glycan_composition}, GlyToucan: {glytoucan_ac}, Glycan Mass: {glycan_mass}"
    
    # Create unique output filename for each row using metadata
    output_file = os.path.join(output_dir, f"{protein}_{site}_{peptide_sequence}_{glytoucan_ac}_mock_mass_spectrum.png").replace('|', '_') # safe file name

    # Create single row dataframe for plotting
    row_df = pd.DataFrame([row])

    # If peptide (peptide_sequence) length is greater than 50 stop for loop
    if len(peptide_sequence) > 50:
       continue
    
    # Generate plot for this specific row
    fig, _ = create_mass_spectrum_plot(row_df, color_mapping, title)
    
    # Save plot to unique filename
    save_plot(fig, output_file)

if __name__ == "__main__":
  main()