# N_Glycopeptide_Sequence_Finder

17JAN2025 -- Richard Shipman

## Overview

Glycopeptide Sequence Finder is a Python script designed to process protein sequences from a FASTA file, digest/cleave them using user-specified proteases, identify N-linked glycopeptides using N-sequon "N[^P][STC]" (NX[STC], X NOT P), and calculate their properties, such as mass, hydrophobicity, and glycosylation sites. The output is written to a CSV file, enabling easy integration into downstream analyses.

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
        - **All:** Runs all proteases above.
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
sp|A0A087X1C5|CP2D7_HUMAN,416,GTTLITNLSSVLK,410,423,13,NLSS,1345.78168089469,0.66,8.75
sp|A0A0K2S4Q6|CD3CH_HUMAN,100,SDQVIITDHPGDLTFTVTLENLTADDAGK,80,109,29,NLTA,3085.50915394004,-0.23,4.05
sp|A0A1B0GTW7|CIROP_HUMAN,333,ENCSTR,332,338,6,NCST,708.2860878938,-1.75,6.09
sp|A0A1B0GTW7|CIROP_HUMAN,425,DSGWYQVNHSAAEELLWGQGSGPEFGLVTTCGTGSSDFFCTGSGLGCHYLHLDK,418,472,54,NHSA,5720.534757730061,-0.22,4.56
sp|A0A1B0GTW7|CIROP_HUMAN,491,MYKPLANGSECWK,485,498,13,NGSE,1525.70575615645,-0.75,7.93
sp|A0A1B0GTW7|CIROP_HUMAN,524,CFFANLTSQLLPGDKPR,520,537,17,NLTS,1905.9771033092898,-0.16,8.22
sp|A0A1B0GTW7|CIROP_HUMAN,713,KPLEVYHGGANFTTQPSK,703,721,18,NFTT,1973.00067520484,-0.91,8.51
sp|A0AV02|S12A8_HUMAN,221,LQLLLLFLLAVSTLDFVVGSFTHLDPEHGFIGYSPELLQNNTLPDYSPGESFFTVFGVFFPAATGVMAGFNMGGDLR,182,259,77,NNTL,8353.194115527771,0.56,4.14
sp|A0AV02|S12A8_HUMAN,561,SEGTQPEGTYGEQLVPELCNQSESSGEDFFLK,542,574,32,NQSE,3504.5514893278296,-0.92,4.05
sp|A0AV02|S12A8_HUMAN,645,ASPGLHLGSASNFSFFR,634,651,17,NFSF,1793.88491717661,0.16,9.8
sp|A0AV96|RBM47_HUMAN,302,EDAVHAMNNLNGTELEGSCLEVTLAKPVDK,292,322,30,NGTE,3196.53802862351,-0.32,4.35
sp|A0AVF1|IFT56_HUMAN,82,ALEEYENATK,76,86,10,NATK,1166.5455329502602,-1.25,4.25
sp|A0AVF1|IFT56_HUMAN,248,SLMDNASSSFEFAK,244,258,14,NASS,1532.6817091532398,-0.19,4.37
sp|A0AVF1|IFT56_HUMAN,408,AATGNTSEGEEAFLLIQSEK,404,424,20,NTSE,2094.01169338101,-0.42,4.09
```

## Protease Rules

The following proteases are supported:

| Protease      | Cleavage Rule                        |
|---------------|--------------------------------------|
| Trypsin       | After K or R, not P                  |
| Chymotrypsin  | After F, W, or Y, not P              |
| Glu-C         | After E                              |
| Lys-C         | After K                              |
| Arg-C         | After R                              |
| Asp-N         | Before D                             |
| Pepsin        | After F, L, W, or Y                  |
| Proteinase K  | After A, F, I, L, V, W, or Y         |
| All           | Runs all proteases above             |

## Notes

1. The script assumes well-formatted FASTA input files.
2. Only N-linked glycosylation sequons are detected (no O-linked or other modifications).

## License

This script is released under the MIT License. 

## Acknowledgments

- The Hitchhiker’s Guide to Glycoproteomics

Oliveira, Tiago, Morten Thaysen-Andersen, Nicolle Packer, and Daniel Kolarich. “The Hitchhiker’s Guide to Glycoproteomics.” Biochemical Society Transactions 49 (July 20, 2021). https://doi.org/10.1042/BST20200879.

- BioPython for handling FASTA files.

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., … others. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423.

- Pyteomics for accurate peptide mass calculations.

Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6

Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717
