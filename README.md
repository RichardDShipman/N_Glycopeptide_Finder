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
        - **Proteinase K:** Cleaves after A, F, I, L, V, W, or Y.
        - **all:** Runs all proteases above.
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
ProteinID,Site,Peptide,Start,End,Length,NSequon,PredictedMass,Hydrophobicity,pI
sp|A0A087X1C5|CP2D7_HUMAN,416,GTTLITNLSSVLK,410,423,13,NLS,1345.78168089469,0.66,8.75
sp|A0A0K2S4Q6|CD3CH_HUMAN,100,SDQVIITDHPGDLTFTVTLENLTADDAGK,80,109,29,NLT,3085.50915394004,-0.23,4.05
sp|A0A1B0GTW7|CIROP_HUMAN,333,ENCSTR,332,338,6,NCS,708.2860878938,-1.75,6.09
sp|A0A1B0GTW7|CIROP_HUMAN,425,DSGWYQVNHSAAEELLWGQGSGPEFGLVTTCGTGSSDFFCTGSGLGCHYLHLDK,418,472,54,NHS,5720.534757730061,-0.22,4.56
sp|A0A1B0GTW7|CIROP_HUMAN,491,MYKPLANGSECWK,485,498,13,NGS,1525.70575615645,-0.75,7.93
sp|A0A1B0GTW7|CIROP_HUMAN,524,CFFANLTSQLLPGDKPR,520,537,17,NLT,1905.9771033092898,-0.16,8.22
sp|A0A1B0GTW7|CIROP_HUMAN,713,KPLEVYHGGANFTTQPSK,703,721,18,NFT,1973.00067520484,-0.91,8.51
sp|A0AV02|S12A8_HUMAN,221,LQLLLLFLLAVSTLDFVVGSFTHLDPEHGFIGYSPELLQNNTLPDYSPGESFFTVFGVFFPAATGVMAGFNMGGDLR,182,259,77,NNT,8353.194115527771,0.56,4.14
sp|A0AV02|S12A8_HUMAN,561,SEGTQPEGTYGEQLVPELCNQSESSGEDFFLK,542,574,32,NQS,3504.5514893278296,-0.92,4.05
sp|A0AV02|S12A8_HUMAN,645,ASPGLHLGSASNFSFFR,634,651,17,NFS,1793.88491717661,0.16,9.8
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