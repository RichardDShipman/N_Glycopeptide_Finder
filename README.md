# Glycopeptide_Sequence_Finder

17JAN2025 -- Richard Shipman

## Overview

Glycopeptide Sequence Finder is a Python script that processes protein sequences from a FASTA file to find amino acid sequences which may or may not contain the post-translational modification glycosylation, the attachment of glycans (polysaccharides) to protein sequences. It uses user-specified proteases to digest and cleave protein sequences into amino acid sequences. The script then identifies N-linked glycopeptides using glycosylation sequon (motifs) like the N-sequon “N[^P][STC]” (NX[STC], where X is not P). It calculates the properties of these glycopeptides, including mass, hydrophobicity, and glycosylation sites. Additionally, the script gathers information from the inputted FASTA file to create a predicted digested glycopeptide (peptide sequence backbone) library. The output is written to a CSV file, making it easy to integrate into downstream analyses.

Here is a table of contents for your README file:

# Table of Contents
1. [Overview](#overview)
2. [Features](#features)
    - [Protease-Specific Cleavage](#protease-specific-cleavage)
    - [Missed Cleavages](#missed-cleavages)
    - [Glycosylation Type](#glycosylation-type)
    - [Peptide Property Calculation](#peptide-property-calculation)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Usage](#usage)
    - [Command-Line Arguments](#arguments)
    - [Example Usage](#example)
    - [Example CSV Output](#example-csv-content)
6. [Protease Rules](#protease-rules)
7. [Glycosylation Type Rules](#glycosylation-type-rules)

### Additional Features

8. [Create Glycan Mass Library](#create-glycan-mass-library)
    - [Input File Format](#input-file-format)
    - [Output File Format](#output-file-format)
9. [Glycan Hydrophobicity Ranking](#glycan-hydrophobicity-ranking)
10. [Glycopeptide Hydrophobicity Calculation](#glycopeptide-hydrophobicity-calculation)
11. [Batch Processing Scripts](#batch-processing-scripts)
    - [Batch Run for FASTA Processing](#batch-run)
    - [Batch Run for Glycopeptide Library Calculation](#batch-run)
12. [Dockerfile](#dockerfile)
13. [Merging CSV Files](#merging-csv-files)

### Reference Materials

14. [License](#license)
15. [Acknowledgments](#acknowledgments)
16. [Appendix](#appendix)
    - [Log File Details](#log-file)
    - [Test Proteomes List](#test-proteomes)
    - [Glycan Mass Library](#glycan-mass-library)

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
3. **Glycosylation Type:** 
    - Select from N-linked (N), O-linked (O), or C-linked (C) glycopeptides. (Adjust or add sequon)
        - **N-linked:** N-sequon “N[^P][STC]"
        - **O-linked:** O-sequon “[ST]" 
        - **C-linked:** C-sequon "W..[WCF]"
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
pip install argparse biopython pyteomics pandas
```

## Usage

Run the script from the command line with the following arguments:

```sh
python glycopeptide_finder_cmd.py -i <input_fasta> [-o <output_csv>] [-p <protease>] [-g <glycosylation>] [-c <missed_cleavages>] [-l <log.txt>] [-v]
```

### Arguments

- `-i`, `--input` (required): Path to the input FASTA file.
- `-o`, `--output` (optional): Path to the output CSV file. If omitted, a default name is generated.
- `-p`, `--protease` (optional): Protease to use for cleavage. Default is trypsin.
- `-g`, `--glycosylation` (optional): Glycosylation sequon to find in peptides. Default is N-linked. (N, O, C) Warning when using O or C, experimental.
- `-c`, `--missed_cleavages` (optional): Number of missed cleavages allowed. Default is 0.
- `-l log.txt`, `--log log.txt` (optional): Path to the log file. If omitted, logging is disabled.
- `-v`, `--verbose` (optional): Enable verbose output. Default is False.
- `-y`, `--glycan`: Path to the glycan file (CSV format) (Default, 4 glycans stored in file). 
- `-z`, `--charge`: (Optional) Maximum charge state to compute (default: 5).

### Example

```sh
python glycopeptide_finder_cmd.py -i test_proteomes/human_uniprotkb_proteome_UP000005640_AND_revi_2025_01_17.fasta -p trypsin -g N -c 0 -z 3
```

The output file will be named:

`example_predicted_trypsin_glycopeptides.csv`

### Example CSV Content

```CSV
ProteinID,Site,GlyToucan_AC,Composition,shorthand_glycan,Peptide,Start,End,Length,Sequon,GlycopeptideMass,PeptideMass,GlycanMass,Hydrophobicity,pI,Protease,GlycosylationType,MissedCleavages,Species,TaxonID,GeneName,ProteinEvidence,SequenceVersion,z2,HF_experimental
sp|A0A087X1C5|CP2D7_HUMAN,416,G80920RR,HexNAc(2)Hex(9),N2H9,GTTLITNLSSVLK,410,423,13,NLSS,3210.41583789469,1345.78168089469,1864.634157,0.66154,8.75,trypsin,N,0,Homo sapiens,9606.0,CYP2D7,5.0,1.0,1606.215194947345,0.64172
sp|A0A087X1C5|CP2D7_HUMAN,416,G62765YT,HexNAc(2)Hex(8),N2H8,GTTLITNLSSVLK,410,423,13,NLSS,3048.36301389469,1345.78168089469,1702.581333,0.66154,8.75,trypsin,N,0,Homo sapiens,9606.0,CYP2D7,5.0,1.0,1525.188782947345,0.59864
sp|A0A087X1C5|CP2D7_HUMAN,416,G31852PQ,HexNAc(2)Hex(7),N2H7,GTTLITNLSSVLK,410,423,13,NLSS,2886.31019089469,1345.78168089469,1540.52851,0.66154,8.75,trypsin,N,0,Homo sapiens,9606.0,CYP2D7,5.0,1.0,1444.162371447345,0.60431
sp|A0A087X1C5|CP2D7_HUMAN,416,G41247ZX,HexNAc(2)Hex(6),N2H6,GTTLITNLSSVLK,410,423,13,NLSS,2724.2573668946898,1345.78168089469,1378.475686,0.66154,8.75,trypsin,N,0,Homo sapiens,9606.0,CYP2D7,5.0,1.0,1363.135959447345,0.61807
sp|A0A0K2S4Q6|CD3CH_HUMAN,100,G80920RR,HexNAc(2)Hex(9),N2H9,SDQVIITDHPGDLTFTVTLENLTADDAGK,80,109,29,NLTA,4950.1433109400405,3085.50915394004,1864.634157,-0.22759,4.05,trypsin,N,0,Homo sapiens,9606.0,CD300H,1.0,1.0,2476.07893147002,-0.24741
sp|A0A0K2S4Q6|CD3CH_HUMAN,100,G62765YT,HexNAc(2)Hex(8),N2H8,SDQVIITDHPGDLTFTVTLENLTADDAGK,80,109,29,NLTA,4788.09048694004,3085.50915394004,1702.581333,-0.22759,4.05,trypsin,N,0,Homo sapiens,9606.0,CD300H,1.0,1.0,2395.05251947002,-0.29049
sp|A0A0K2S4Q6|CD3CH_HUMAN,100,G31852PQ,HexNAc(2)Hex(7),N2H7,SDQVIITDHPGDLTFTVTLENLTADDAGK,80,109,29,NLTA,4626.03766394004,3085.50915394004,1540.52851,-0.22759,4.05,trypsin,N,0,Homo sapiens,9606.0,CD300H,1.0,1.0,2314.02610797002,-0.28482
sp|A0A0K2S4Q6|CD3CH_HUMAN,100,G41247ZX,HexNAc(2)Hex(6),N2H6,SDQVIITDHPGDLTFTVTLENLTADDAGK,80,109,29,NLTA,4463.98483994004,3085.50915394004,1378.475686,-0.22759,4.05,trypsin,N,0,Homo sapiens,9606.0,CD300H,1.0,1.0,2232.9996959700197,-0.27106
sp|A0A1B0GTW7|CIROP_HUMAN,333,G80920RR,HexNAc(2)Hex(9),N2H9,ENCSTR,332,338,6,NCST,2572.9202448938,708.2860878938,1864.634157,-1.75,6.09,trypsin,N,0,Homo sapiens,9606.0,CIROP,1.0,1.0,1287.4673984469,-1.76982
sp|A0A1B0GTW7|CIROP_HUMAN,333,G62765YT,HexNAc(2)Hex(8),N2H8,ENCSTR,332,338,6,NCST,2410.8674208938,708.2860878938,1702.581333,-1.75,6.09,trypsin,N,0,Homo sapiens,9606.0,CIROP,1.0,1.0,1206.4409864469,-1.8129
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

## Glycosylation Type Rules

The following glycosylation types sequons (motifs) are supported:

| Glycosylation Type | Sequon Pattern |
|--------------------|----------------|
| N-linked           | N[^P][STC]     |
| O-linked           | [ST]           |
| C-linked           | W..[WCF]       |

## Glycan Library

The default glycan mass library is defined as a DataFrame containing a set of glycans with their respective compositions and masses. This library is used to calculate the properties of glycopeptides. Alter if you wish to change the glycan mass library in the script

```python
import pandas as pd

default_glycan_library = pd.DataFrame([
    {"glytoucan_ac": "G80920RR", "byonic": "HexNAc(2)Hex(9) % 1864.634157", "composition": "HexNAc(2)Hex(9)", "mass": 1864.634157}, # N2H9
    {"glytoucan_ac": "G62765YT", "byonic": "HexNAc(2)Hex(8) % 1702.581333", "composition": "HexNAc(2)Hex(8)", "mass": 1702.581333}, # N2H8
    {"glytoucan_ac": "G31852PQ", "byonic": "HexNAc(2)Hex(7) % 1540.528510", "composition": "HexNAc(2)Hex(7)", "mass": 1540.528510}, # N2H7
    {"glytoucan_ac": "G41247ZX", "byonic": "HexNAc(2)Hex(6) % 1378.475686", "composition": "HexNAc(2)Hex(6)", "mass": 1378.475686}, # N2H6
])
```

This DataFrame includes the following columns:
- `glytoucan_ac`: The glycosylation context identifier.
- `byonic`: The peptide sequence and mass in Byonic format.
- `composition`: The glycan composition.
- `mass`: The mass of the glycan.

The default glycan mass library can be expanded or customized as needed for specific analyses.

### Example glycan mass library data

This data is stored in the `glycan_mass_library` directory.

```csv
glytoucan_ac,byonic,composition,mass
G62765YT,HexNAc(2)Hex(8) % 1702.581333,HexNAc(2)Hex(8),1702.581333
G31852PQ,HexNAc(2)Hex(7) % 1540.528510,HexNAc(2)Hex(7),1540.528510
G41247ZX,HexNAc(2)Hex(6) % 1378.475686,HexNAc(2)Hex(6),1378.475686
```

## Create Glycan Mass Library 

create_glycan_library.py

This script processes glycan data from a CSV file and splits it into columns based on specific formatting rules.

## Overview
The script takes an input CSV file with glyc肽 data formatted in two columns:
1. **glytoucan_ac**: A string representing the glycosylation context.
2. **byonic** / **byonic_sequence**: A composite column containing two pieces of information separated by a '%' character:
    - The peptide sequence (composition).
    - A numerical value representing mass.

The script splits the `byonic` column into its constituent parts, creating three new columns in the output file:
1. **glytoucan_ac**: The glycosylation context.
2. **composition**: The peptide sequence.
3. **mass**: The numerical mass value.

## Key Features

### Input File Format
The input CSV file must have exactly two columns per row, with the second column formatted as `<peptide> % <mass>`.

Example:
```
"glytoucan_ac","sequence_byonic","name_source"
"G00002CF","Hex(2)NeuGc(2) % 956.29687423","G93218EI"
"G00009BX","HexNAc(2)Hex(2)dHex(1) % 894.33286578","G93579XB"
"G00012MO","HexNAc(1)Hex(3) % 707.248407805","G08590QR"
"G00012RZ","HexNAc(5)Hex(6)dHex(3)NeuAc(1) % 2734.99351197","G04784US"
"G00013MO","Hex(4) % 666.2218587","G90306QO"
"G00014MO","HexNAc(2)Hex(2) % 748.27495691","G53434XO"
"G00015MO","HexNAc(1)Hex(3) % 707.248407805","G08590QR"
"G00016MO","Hex(2) % 342.1162117","G90627TW"
"G00017IP","HexNAc(4)Hex(3)Sulpho(1) % 1396.44334021","G67486RJ"
"G00024MO","Hex(3) % 504.1690352","G39365VM"

```

### Output File Format
The output file will have three columns:
- The first column (`glytoucan_ac`) remains unchanged.
- The second column contains the peptide sequence from `byonic`.
- The third column contains the numerical mass value.

Example:
```
glytoucan_ac,byonic,name_source,composition,mass
G00002CF,Hex(2)NeuGc(2) % 956.29687423,G93218EI,Hex(2)NeuGc(2),956.29687423
G00009BX,HexNAc(2)Hex(2)dHex(1) % 894.33286578,G93579XB,HexNAc(2)Hex(2)dHex(1),894.33286578
G00012MO,HexNAc(1)Hex(3) % 707.248407805,G08590QR,HexNAc(1)Hex(3),707.248407805
G00012RZ,HexNAc(5)Hex(6)dHex(3)NeuAc(1) % 2734.99351197,G04784US,HexNAc(5)Hex(6)dHex(3)NeuAc(1),2734.99351197
G00013MO,Hex(4) % 666.2218587,G90306QO,Hex(4),666.2218587
G00014MO,HexNAc(2)Hex(2) % 748.27495691,G53434XO,HexNAc(2)Hex(2),748.27495691
G00015MO,HexNAc(1)Hex(3) % 707.248407805,G08590QR,HexNAc(1)Hex(3),707.248407805
G00016MO,Hex(2) % 342.1162117,G90627TW,Hex(2),342.1162117
G00017IP,HexNAc(4)Hex(3)Sulpho(1) % 1396.44334021,G67486RJ,HexNAc(4)Hex(3)Sulpho(1),1396.44334021
G00024MO,Hex(3) % 504.1690352,G39365VM,Hex(3),504.1690352
```

## Glycan Hydrophobicity Ranking

This script ranks glycans based on their adjusted hydrophobicity factor (HF). The analysis removes the peptide effect, computes a weighted HF score considering glycan frequency, and normalizes the adjusted HF values using Z-scores. The glycans are then ranked by their weighted adjusted HF. Note, more glycopeptide data is needed for this to hold any value, work in progress.

### Functions:

- `compute_adjusted_hf(df)`: Adjusts HF by removing the peptide effect.
- `compute_weighted_adjusted_hf(df)`: Normalizes the adjusted HF and computes a weighted HF based on glycan frequency.
- `count_glycan_frequency(df)`: Counts the frequency of each glycan in the dataset.
- `rank_glycans_by_hydrophobicity(input_csv, output_csv)`: Main function that processes the input CSV, computes the adjusted HF, and outputs the glycan ranking by hydrophobicity.

```sh
python script.py -i input_file.csv -o output_file.csv
```

Arguments:
- `-i`, `--input`: Input CSV file path (required).
- `-o`, `--output`: Output CSV file path (optional, defaults to `<input_file>_glycan_hydrophobicity_index.csv`).

The output CSV will contain glycans ranked by their weighted adjusted HF score.

## Glycopeptide Hydrophobicity Calculation

This script calculates hydrophobicity scores for glycopeptides based on peptide hydrophobicity and glycan composition. It takes as input a glycopeptide data file and a glycan hydrophobicity data file, then outputs the results to a CSV file. Note, more glycopeptide data is needed for this to hold any value, work in progress.

Usage

```sh
python script.py -i <glycopeptide_file> -gh <glycan_hf_file> -o <output_file>
```

Arguments:
- `-i`: Input file containing glycopeptide data (CSV).
- `-gh`: Input file containing glycan hydrophobicity data (CSV, optional if default is used).
- `-o`: Output file name (optional, defaults to `<input_file_name>_HF.csv`).

Workflow
	1.	Load glycopeptide and glycan hydrophobicity data.
	2.	Calculate hydrophobicity scores for each glycopeptide.
	3.	Save the results to a CSV file.

# Batch Processing Scripts

Shell scripts for batch processing.

## Batch Run for FASTA Processing

batch_glycopeptide_sequence_finder.sh

To process multiple FASTA files in parallel using all proteases, run the following command:

```sh
./batch_glycopeptide_sequence_finder.sh
```

Parameters can be adjusted in the shell script.

### Parameters

- `ls test_proteomes/*.fasta`: Lists all FASTA files in the `test_proteomes` directory.
- `xargs -I {} -P 4`: Executes the command in parallel with up to 4 processes. The `{}` is a placeholder for each file name.
- `python glycopeptide_finder_cmd.py`: The script to run for each FASTA file.
- `-i "{}"`: Specifies the input FASTA file, where `{}` is replaced by each file name.
- `-p all`: Uses all proteases for cleavage.
- `-g N`: Searches for N-linked glycosylation sequons.
- `-c 0`: Allows 0 missed cleavages.
- `-v`: Enables verbose output.

This command allows you to efficiently process multiple FASTA files in parallel, reducing the overall processing time.

## Batch Run for Glycopeptide Library Calculation

batch_compute_intact_glycopeptide_library.sh

To process multiple digested glycopeptide libraries in parallel, run the following command:

```sh
./batch_compute_intact_glycopeptide_library.sh
```

This command will process all CSV files in the `digested_glycopeptide_library` directory using the specified glycan library and compute intact glycopeptide masses for charge states up to 6 in parallel with up to 4 processes.

### Parameters

- `ls test_proteomes/*.fasta`: Lists all FASTA files in the `test_proteomes` directory.
- `xargs -I {} -P 4`: Executes the command in parallel with up to 4 processes. The `{}` is a placeholder for each file name.
- `python n_glycopeptide_sequence_finder_cmd.py`: The script to run for each FASTA file.
- `-i "{}"`: Specifies the input FASTA file, where `{}` is replaced by each file name.
- `-p trypsin`: Uses trypsin as the protease for cleavage.
- `-g N`: Searches for N-linked glycosylation sequons.
- `-c 0`: Allows 0 missed cleavages.
- `-v`: Enables verbose output.

This command allows you to efficiently process multiple FASTA files in parallel, reducing the overall processing time.

## Merging CSV Files

The script includes a function to merge all CSV files from a specified directory into a single CSV file. This can be useful for consolidating the results of multiple digestions into one file for easier analysis.

```sh
python merge_digested_glycopeptide_library.py
```

## Dockerfile

- Docker Setup for Glycopeptide Sequence Finder

This section explains how to build and run the Docker container for the Glycopeptide Sequence Finder.

1. Build the Docker Image

To create the Docker image, run the following command in the directory containing your Dockerfile and requirements.txt:

```sh
docker build -t gsf .
```

This will:
	•	Use the official Python 3.10-slim image.
	•	Set /app as the working directory.
	•	Install dependencies from requirements.txt.
	•	Copy the glycopeptide_sequence_finder_cmd.py script into the container.
	•	Set the entrypoint so that the script can be executed with arguments.

2. Run the Docker Container

To execute the script with test data, use:

```sh
docker run --rm \
    -v "$(pwd)/test_proteomes:/app/test_proteomes" \
    -v "$(pwd)/output:/app/digested_glycopeptide_library" \
    gsf \
    -i test_proteomes/apple_uniprotkb_proteome_UP000290289_AND_revi_2025_02_04.fasta \
    -g N \
    -o digested_glycopeptide_library/test.csv \
    -p chymotrypsin \
    -c 0 \
    -v
```

Explanation of Flags:
	•	--rm → Removes the container after execution.
	•	-v "$(pwd)/test_proteomes:/app/test_proteomes" → Mounts the input FASTA files.
	•	-v "$(pwd)/output:/app/digested_glycopeptide_library" → Mounts the output directory.
	•	gsf → Runs the built image.
	•	-i → Specifies the input FASTA file.
	•	-g → Sets the glycosylation type (default: N).
	•	-o → Defines the output file.
	•	-p → Specifies the protease (e.g., chymotrypsin).
	•	-c → Defines the missed cleavages.
	•	-v → Enables verbose mode.

3. Access the Output

The output files will be saved in the mounted directory on your local machine:

```sh
ls output/digested_glycopeptide_library/
```

Your results should be inside output/digested_glycopeptide_library/test.csv.

## License

This script is released under the MIT License. 

## Acknowledgments

- The Hitchhiker’s Guide to Glycoproteomics.

Oliveira, Tiago, Morten Thaysen-Andersen, Nicolle Packer, and Daniel Kolarich. “The Hitchhiker’s Guide to Glycoproteomics.” Biochemical Society Transactions 49 (July 20, 2021). https://doi.org/10.1042/BST20200879.

- In Silico Platform for Prediction of N-, O- and C-Glycosites in Eukaryotic Protein Sequences.

Chauhan, Jagat Singh, Alka Rao, and Gajendra P. S. Raghava. “In Silico Platform for Prediction of N-, O- and C-Glycosites in Eukaryotic Protein Sequences.” PLoS ONE 8, no. 6 (June 28, 2013): e67008. https://doi.org/10.1371/journal.pone.0067008.

- Large-Scale Identification of N-Linked Intact Glycopeptides in Human Serum Using HILIC Enrichment and Spectral Library Search.

Shu, Qingbo, Mengjie Li, Lian Shu, Zhiwu An, Jifeng Wang, Hao Lv, Ming Yang, et al. “Large-Scale Identification of N-Linked Intact Glycopeptides in Human Serum Using HILIC Enrichment and Spectral Library Search.” Molecular & Cellular Proteomics : MCP 19, no. 4 (April 2020): 672–89. https://doi.org/10.1074/mcp.RA119.001791.

- Assessing the Hydrophobicity of Glycopeptides Using Reversed-Phase Liquid Chromatography and Tandem Mass Spectrometry.

Wang, Junyao, Aiying Yu, Byeong Gwan Cho, and Yehia Mechref. “Assessing the Hydrophobicity of Glycopeptides Using Reversed-Phase Liquid Chromatography and Tandem Mass Spectrometry.” Journal of Chromatography. A 1706 (September 13, 2023): 464237. https://doi.org/10.1016/j.chroma.2023.464237.

- Molecular Basis of C-Mannosylation – a Structural Perspective.

Crine, Samuel L., and K. Ravi Acharya. “Molecular Basis of C-Mannosylation – a Structural Perspective.” The FEBS Journal 289, no. 24 (2022): 7670–87. https://doi.org/10.1111/febs.16265.

- BioPython for handling FASTA files.

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., … others. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423.

- Pyteomics for accurate peptide mass calculations.

Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6

Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717

- GlyGen: Computational and Informatics Resources for Glycoscience.

York WS, Mazumder R, Ranzinger R, Edwards N, Kahsay R, Aoki-Kinoshita KF, Campbell MP, Cummings RD, Feizi T, Martin M, Natale DA, Packer NH, Woods RJ, Agarwal G, Arpinar S, Bhat S, Blake J, Castro LJG, Fochtman B, Gildersleeve J, Goldman R, Holmes X, Jain V, Kulkarni S, Mahadik R, Mehta A, Mousavi R, Nakarakommula S, Navelkar R, Pattabiraman N, Pierce MJ, Ross K, Vasudev P, Vora J, Williamson T, Zhang W. GlyGen: Computational and Informatics Resources for Glycoscience. Glycobiology. 2020 Jan 28;30(2):72-73. doi: 10.1093/glycob/cwz080. PMID: 31616925; PMCID: PMC7335483.

# Appendix

Additional information, logging runs, and references.

## Log File

The log file provides detailed information about the processing steps, including the number of peptides found after cleavage and the identified N-glycopeptides.

```bash
python glycopeptide_finder_cmd.py -i test_proteomes/human_uniprotkb_proteome_UP000005640_AND_revi_2025_01_17.fasta -p trypsin -c 0 -l log.txt
```

The script generates a log file that records the processing details of each protein sequence. Logging to a text file can be activated with the `-l log.txt` flag. Below are some example log entries:

```log
2025-02-04 00:57:21,942 - INFO - Processing sp|A0A087X1C5|CP2D7_HUMAN with 515 amino acids.
2025-02-04 00:57:21,942 - INFO - Found 50 peptides after trypsin cleavage. The peptides were: ['MGLEALVPLAMIVAIFLLLVDLMHR', 'HQR', 'WAAR', 'YPPGPLPLPGLGNLLHVDFQNTPYCFDQLR', 'R', 'R', 'FGDVFSLQLAWTPVVVLNGLAAVR', 'EAMVTR', 'GEDTADRPPAPIYQVLGFGPR', 'SQGVILSR', 'YGPAWR', 'EQR', 'R', 'FSVSTLR', 'NLGLGK', 'K', 'SLEQWVTEEAACLCAAFADQAGRPFRPNGLLDK', 'AVSNVIASLTCGR', 'R', 'FEYDDPR', 'FLR', 'LLDLAQEGLK', 'EESGFLR', 'EVLNAVPVLPHIPALAGK', 'VLR', 'FQK', 'AFLTQLDELLTEHR', 'MTWDPAQPPR', 'DLTEAFLAK', 'K', 'EK', 'AK', 'GSPESSFNDENLR', 'IVVGNLFLAGMVTTSTTLAWGLLLMILHLDVQR', 'GR', 'R', 'VSPGCPIVGTHVCPVR', 'VQQEIDDVIGQVR', 'RPEMGDQAHMPCTTAVIHEVQHFGDIVPLGVTHMTSR', 'DIEVQGFR', 'IPK', 'GTTLITNLSSVLK', 'DEAVWK', 'KPFR', 'FHPEHFLDAQGHFVKPEAFLPFSAGR', 'R', 'ACLGEPLAR', 'MELFLFFTSLLQHFSFSVAAGQPRPSHSR', 'VVSFLVTPSPYELCAVPR', '']
2025-02-04 00:57:21,943 - INFO - Found 1 N-glycopeptides. The glycopeptides were: [('GTTLITNLSSVLK', 416)]
2025-02-04 00:57:21,943 - INFO - Processing sp|A0A0B4J2F0|PIOS1_HUMAN with 54 amino acids.
2025-02-04 00:57:21,945 - INFO - Found 9 peptides after trypsin cleavage. The peptides were: ['MFR', 'R', 'LTFAQLLFATVLGIAGGVYIFQPVFEQYAK', 'DQK', 'ELK', 'EK', 'MQLVQESEEK', 'K', 'S']
2025-02-04 00:57:21,945 - INFO - Found 0 N-glycopeptides. The glycopeptides were: []
2025-02-04 00:57:21,945 - INFO - Processing sp|A0A0C5B5G6|MOTSC_HUMAN with 16 amino acids.
2025-02-04 00:57:21,945 - INFO - Found 5 peptides after trypsin cleavage. The peptides were: ['MR', 'WQEMGYIFYPR', 'K', 'LR', '']
2025-02-04 00:57:21,945 - INFO - Found 0 N-glycopeptides. The glycopeptides were: []
2025-02-04 00:57:21,945 - INFO - Processing sp|A0A0K2S4Q6|CD3CH_HUMAN with 201 amino acids.
2025-02-04 00:57:21,945 - INFO - Found 13 peptides after trypsin cleavage. The peptides were: ['MTQR', 'AGAAMLPSALLLLCVPGCLTVSGPSTVMGAVGESLSVQCR', 'YEEK', 'YK', 'TFNK', 'YWCR', 'QPCLPIWHEMVETGGSEGVVR', 'SDQVIITDHPGDLTFTVTLENLTADDAGK', 'YR', 'CGIATILQEDGLSGFLPDPFFQVQVLVSSASSTENSVK', 'TPASPTRPSQCQGSLPSSTCFLLLPLLK', 'VPLLLSILGAILWVNRPWR', 'TPWTES']
2025-02-04 00:57:21,945 - INFO - Found 1 N-glycopeptides. The glycopeptides were: [('SDQVIITDHPGDLTFTVTLENLTADDAGK', 100)]
```

### Notes

1. The script assumes well-formatted FASTA input files.
2. Only N-linked glycosylation sequons are detected (no O-linked or other modifications).
3. FASTA protein files contain new lines and or return carrages. When returning to the FASTA, remember this when searching for peptide in original sequence. 

```
for file in ./test_proteomes/*; do
  filename=$(basename "$file")
  part_before_underscore="${filename%%_*}"
  echo "$part_before_underscore"
done
```

## Test Proteomes

Test proteome files from UniProt are available in the `test_proteomes` folder. Below is a list of species gathered. Only reviewed proteins were downloaded, and not every sequence available for a species is included. 

| Common Name   | Scientific Name                                                                 | Taxon ID                |
|---------------|---------------------------------------------------------------------------------|-------------------------|
| SARS-CoV      | SARS-CoV (Severe Acute Respiratory Syndrome Coronavirus)                        | 694009                  |
| Apple         | Malus domestica                                                                 | 3750                    |
| Arabidopsis   | Arabidopsis thaliana                                                            | 3702                    |
| Bat           | Myotis lucifugus                                                                | 59463                   |
| Cat           | Felis catus                                                                     | 9685                    |
| C. elegans    | Caenorhabditis elegans                                                          | 6239                    |
| Chicken       | Gallus gallus                                                                   | 9031                    |
| Chimpanzee    | Pan troglodytes                                                                 | 9598                    |
| C. jejuni     | Campylobacter jejuni                                                            | 1951                    |
| Cow           | Bos taurus                                                                      | 9913                    |
| Dictyostelium | Dictyostelium discoideum                                                        | 44689                   |
| Dog           | Canis lupus familiaris                                                          | 9615                    |
| Donkey        | Equus asinus                                                                    | 9796                    |
| Duck          | Cairina moschata                                                                | 8855                    |
| Elephant      | Loxodonta africana (African Elephant)                                           | 9785                    |
| Fruit Fly     | Drosophila melanogaster                                                         | 7227                    |
| Goat          | Capra hircus                                                                    | 9925                    |
| Guinea Pig    | Cavia porcellus                                                                 | 10141                   |
| Honeybee      | Apis mellifera                                                                  | 7460                    |
| Horse         | Equus caballus                                                                  | 9796                    |
| Human         | Homo sapiens                                                                    | 9606                    |
| Mouse         | Mus musculus                                                                    | 10090                   |
| Orangutan     | Pongo abelii                                                                    | 9601                    |
| Pig           | Sus scrofa domesticus                                                           | 9823                    |
| Rat           | Rattus norvegicus                                                               | 10116                   |
| Rice          | Oryza sativa subsp. japonica                                                    | 39947                   |
| Sheep         | Ovis aries                                                                      | 9940                    |
| Sorghum       | Sorghum bicolor                                                                 | 4558                    |
| Squirrel      | Ictidomys tridecemlineatus                                                      | 43179                   |
| Yeast         | Saccharomyces cerevisiae (strain ATCC 204508 / S288c)                           | 559292                  |
| Zebrafish     | Danio rerio                                                                     | 7955                    |

## Glycan Mass Library

Glycan mass libraries were gathered from GlyGen. The follow data was processed and used in this tool.

It is advisable to create a targeted glycan library along with a list of peptides to compute an intact glycan library. The size of the library can grow rapidly, so it is important to manage it effectively.

| File Name          | Glycan Count |
|--------------------|--------------|
| glycan_database     | 44686        |
| glycan_type_n_linked_byonic.csv | 369    |
| test_glycan_library.csv | 3       |
