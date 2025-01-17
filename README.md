# N_Glycopeptide_Sequence_Finder

17JAN2025 -- Richard Shipman

## Overview

Glycopeptide Sequence Finder is a Python script designed to process protein sequences from a FASTA file, digest/cleave them using user-specified proteases, identify N-linked glycopeptides, and calculate their properties, such as mass, hydrophobicity, and glycosylation sites. The output is written to a CSV file, enabling easy integration into downstream analyses.

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
3. **Peptide Property Calculation:**
    - Calculates peptide mass, hydrophobicity, and isoelectric point (pI).

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
python n_glycopeptide_finder_cmd.py -i <input_fasta> [-o <output_csv>] [-p <protease>] [-c <missed_cleavages>]
```

### Arguments

- `-i`, `--input` (required): Path to the input FASTA file.
- `-o`, `--output` (optional): Path to the output CSV file. If omitted, a default name is generated.
- `-p`, `--protease` (optional): Protease to use for cleavage. Default is trypsin.
- `-c`, `--missed_cleavages` (optional): Number of missed cleavages allowed. Default is 0.

### Example

```sh
python n_glycopeptide_finder_cmd.py -i example.fasta -p chymotrypsin -c 1
```

The output file will be named:

`example_predicted_chymotrypsin_glycopeptides.csv`

### Example CSV Content

```csv
ProteinID,Peptide,Site,Start,End,Length,NSequon,PredictedMass,Hydrophobicity,pI
sp|A0A0B7P3V8|YP41B_YEAST,NVIDDNISAR,16,11,21,10,NIS,1115.55710069347,-0.43,4.21
sp|A0A0B7P3V8|YP41B_YEAST,TNDTVR,28,27,33,6,NDT,704.3453170220799,-1.45,5.5
sp|A0A0B7P3V8|YP41B_YEAST,EGLGESLDIMNTNTTDIFR,211,199,218,19,NTT,2124.99974879832,-0.42,4.05
sp|A0A0B7P3V8|YP41B_YEAST,ELRPDSTNFSK,368,361,372,11,NFS,1292.63607929016,-1.47,6.17
sp|A0A0B7P3V8|YP41B_YEAST,LVIIDTGSGVNITNDK,420,410,426,16,NIT,1657.88866514437,0.3,4.21
```

## Protease Rules

The following proteases are supported:

| Protease     | Cleavage Rule                        |
|--------------|--------------------------------------|
| Trypsin      | After K or R, not P                  |
| Chymotrypsin | After F, W, or Y, not P              |
| Glu-C        | After E                              |
| Lys-C        | After K                              |
| Arg-C        | After R                              |
| Pepsin       | After F, L, W, or Y                  |

## Notes

1. The script assumes well-formatted FASTA input files.
2. Only N-linked glycosylation sequons are detected (no O-linked or other modifications).

## License

This script is released under the MIT License. 

## Acknowledgments

- BioPython for handling FASTA files.

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., … others. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423.

- Pyteomics for accurate peptide mass calculations.

Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6

Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717