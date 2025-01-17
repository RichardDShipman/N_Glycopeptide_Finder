# N-Glycopeptide Finder

17JAN2025 -- Richard Shipman

## Overview

Glycopeptide Finder is a Python script designed to process protein sequences from a FASTA file, cleave them using user-specified proteases, identify N-linked glycopeptides, and calculate their properties, such as mass and glycosylation sites. The output is written to a CSV file, enabling easy integration into downstream analyses.

## Features

1. **Protease-Specific Cleavage:**
    - Supports several commonly used proteases, including:
        - **Trypsin:** Cleaves after K or R, except if followed by P.
        - **Chymotrypsin:** Cleaves after F, W, or Y, except if followed by P.
        - **Glu-C:** Cleaves after E.
        - **Lys-C:** Cleaves after K.
        - **Arg-C:** Cleaves after R.
        - **Asp-N:** Cleaves before D.
        - **Pepsin:** Cleaves after F, L, W, or Y.
2. **Missed Cleavages:**
    - Allows specifying the number of missed cleavages to simulate incomplete digestion.

## Requirements

- Python 3.7 or later
- Libraries:
    - argparse
    - csv
    - re
    - biopython
    - pyteomics

## Installation

Install the required Python libraries using pip:

```sh
pip install biopython pyteomics
```

## Usage

Run the script from the command line with the following arguments:

```sh
python glycopeptide_finder.py -i <input_fasta> [-o <output_csv>] [-p <protease>] [-c <missed_cleavages>]
```

### Arguments

- `-i`, `--input` (required): Path to the input FASTA file.
- `-o`, `--output` (optional): Path to the output CSV file. If omitted, a default name is generated.
- `-p`, `--protease` (optional): Protease to use for cleavage. Default is trypsin.
- `-c`, `--missed_cleavages` (optional): Number of missed cleavages allowed. Default is 0.

### Example

```sh
python glycopeptide_finder.py -i example.fasta -p chymotrypsin -c 1
```

The output file will be named:

`example_predicted_chymotrypsin_glycopeptides.csv`

### Example CSV Content

| ProteinID | Peptide | Site | Length | NSequon | PredictedMass |
|-----------|---------|------|--------|---------|---------------|
| P12345    | NTTKPN  | 10   | 6      | NTT     | 748.35        |
| P12345    | WFNST   | 25   | 5      | NST     | 643.28        |

## Protease Rules

The following proteases are supported:

| Protease     | Cleavage Rule                        |
|--------------|--------------------------------------|
| Trypsin      | After K or R, not P                  |
| Chymotrypsin | After F, W, or Y, not P              |
| Glu-C        | After E                              |
| Lys-C        | After K                              |
| Arg-C        | After R                              |
| Asp-N        | Before D (and sometimes B)           |
| Pepsin       | After F, L, W, or Y                  |

## Notes

1. The script assumes well-formatted FASTA input files.
2. Only N-linked glycosylation sequons are detected (no O-linked or other modifications).

## License

This script is released under the MIT License. Feel free to use, modify, and distribute it as needed.

## Acknowledgments

- BioPython for handling FASTA files.

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., … others. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423.

- Pyteomics for accurate peptide mass calculations.

Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6

Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717