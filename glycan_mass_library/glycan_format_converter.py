import pandas as pd
import argparse
import re
import os

# Mapping definitions
SHORTHAND_TO_LONGHAND = {
    'H': 'Hex',
    'N': 'HexNAc',
    'F': 'Fuc',
    'A': 'NeuAc',
    'G': 'NeuGc',
    'X': 'Pent',
    'S': 'Sulpho',
    'Su': 'Sulpho',
    'P': 'Phospho',
    'dF': 'dHex'  # Added dHex for F
}

LONGHAND_TO_SHORTHAND = {v: k for k, v in SHORTHAND_TO_LONGHAND.items()}

def convert_shorthand_to_longhand(shorthand):
    """Converts shorthand glycan notation (e.g., H5N4A1) to longhand (e.g., HexNAc(4)Hex(5)NeuAc(1))."""
    matches = re.findall(r'([A-Za-z]+)(\d*)', shorthand)
    longhand_parts = []
    for code, count in matches:
        if code in SHORTHAND_TO_LONGHAND:
            longhand_parts.append(f"{SHORTHAND_TO_LONGHAND[code]}({count or '1'})")
    return ''.join(longhand_parts)

def convert_longhand_to_shorthand(longhand):
    """Converts longhand glycan notation (e.g., HexNAc(4)Hex(5)NeuAc(1)) to shorthand (e.g., H5N4A1)."""
    matches = re.findall(r'(HexNAc|Hex|Fuc|NeuAc|NeuGc|Pent|Sulpho|Phospho|dHex)\((\d+)\)', longhand)
    shorthand = ''.join(f"{LONGHAND_TO_SHORTHAND[name]}{count}" for name, count in matches)
    return shorthand
    
def detect_format(glycan):
    if re.match(r"^H\d+N\d+(A\d+|F\d+|G\d+|X\d+|A\d+)*$", glycan):  # Adjust regex for NeuGc, Xyl
        return "shorthand"
    elif re.match(r'^(HexNAc|Hex|Fuc|NeuAc|NeuGc|Pent|Sulpho|Phospho|dHex)\(\d+\)', glycan):
        return 'longhand'
    else:
        raise ValueError("Unknown glycan format.")
    raise ValueError("Unknown glycan format.")

def process_csv(input_file, output_file, glycan_column):
    try:
        print(f"Reading file: {input_file}")
        df = pd.read_csv(input_file, dtype=str)  # Force all columns to be read as strings
        print(df.head())

        # Detect format of the first glycan
        first_glycan = df[glycan_column].iloc[0]
        print(f"First glycan detected: {first_glycan}")

        format_type = detect_format(first_glycan)
        print(f"Detected format: {format_type}")

        # Add new columns for converted glycans
        df['shorthand_glycan'] = df[glycan_column].apply(
            convert_longhand_to_shorthand if format_type == "longhand" else lambda x: x
        )
        df['longhand_glycan'] = df[glycan_column].apply(
            convert_shorthand_to_longhand if format_type == "shorthand" else lambda x: x
        )

        # Write the updated dataframe to the output file
        df.to_csv(output_file, index=False)
        print(f"File written to: {output_file}")

    except Exception as e:
        print("Error:", e)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert glycan formats between shorthand and longhand notation.")
    parser.add_argument("-i", "--input_file", required=True, help="Path to input CSV file.")
    parser.add_argument("-o", "--output_file", help="Path to output CSV file. Default is {input_filename}_converted.csv in the same directory as the input file.")
    parser.add_argument("-g", "--glycan_column", default="Glycans", help="Name of the glycan column in the CSV file (default: Glycans).")
    parser.add_argument("-d", "--direction", choices=["to_longhand", "to_shorthand"], help="Conversion direction: to_longhand (H5N4 -> HexNAc(4)Hex(5)) or to_shorthand (HexNAc(4)Hex(5) -> H5N4). If not provided, the format will be detected automatically.")
    args = parser.parse_args()

    if not args.output_file:
        input_filename = os.path.splitext(os.path.basename(args.input_file))[0]
        input_dir = os.path.dirname(args.input_file)
        args.output_file = os.path.join(input_dir, f"{input_filename}_converted.csv")
    
    process_csv(args.input_file, args.output_file, args.glycan_column)
    print(f"File written to: {args.output_file}")
