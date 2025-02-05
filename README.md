# Glycopeptide_Sequence_Finder

17JAN2025 -- Richard Shipman

## Overview

Glycopeptide Sequence Finder is a Python script that processes protein sequences from a FASTA file to find amino acid sequences which may or may not contain the post-translational modification glycosylation, the attachment of glycans (polysaccharides) to protein sequences. It uses user-specified proteases to digest and cleave protein sequences into amino acid sequences. The script then identifies N-linked glycopeptides using glycosylation sequon (motifs) like the N-sequon “N[^P][STC]” (NX[STC], where X is not P). It calculates the properties of these glycopeptides, including mass, hydrophobicity, and glycosylation sites. Additionally, the script gathers information from the inputted FASTA file to create a predicted digested glycopeptide (peptide sequence backbone) library. The output is written to a CSV file, making it easy to integrate into downstream analyses.

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

### Example

```sh
python glycopeptide_finder_cmd.py -i test_proteomes/human_uniprotkb_proteome_UP000005640_AND_revi_2025_01_17.fasta -p trypsin -g N -c 0
```

The output file will be named:

`example_predicted_trypsin_glycopeptides.csv`

### Example CSV Content

```CSV
ProteinID,Site,Peptide,Start,End,Length,Sequon,PredictedMass,Hydrophobicity,pI,Protease,GlycosylationType,MissedCleavages,Species,TaxonID,GeneName,ProteinEvidence,SequenceVersion
sp|A0A087X1C5|CP2D7_HUMAN,416,GTTLITNLSSVLK,410,423,13,NLSS,1345.78168089469,0.66,8.75,trypsin,N,0,Homo sapiens,9606,CYP2D7,5,1
sp|A0A0K2S4Q6|CD3CH_HUMAN,100,SDQVIITDHPGDLTFTVTLENLTADDAGK,80,109,29,NLTA,3085.50915394004,-0.23,4.05,trypsin,N,0,Homo sapiens,9606,CD300H,1,1
sp|A0A1B0GTW7|CIROP_HUMAN,333,ENCSTR,332,338,6,NCST,708.2860878938,-1.75,6.09,trypsin,N,0,Homo sapiens,9606,CIROP,1,1
sp|A0A1B0GTW7|CIROP_HUMAN,425,DSGWYQVNHSAAEELLWGQGSGPEFGLVTTCGTGSSDFFCTGSGLGCHYLHLDK,418,472,54,NHSA,5720.534757730061,-0.22,4.56,trypsin,N,0,Homo sapiens,9606,CIROP,1,1
sp|A0A1B0GTW7|CIROP_HUMAN,491,MYKPLANGSECWK,485,498,13,NGSE,1525.70575615645,-0.75,7.93,trypsin,N,0,Homo sapiens,9606,CIROP,1,1
sp|A0A1B0GTW7|CIROP_HUMAN,524,CFFANLTSQLLPGDKPR,520,537,17,NLTS,1905.9771033092898,-0.16,8.22,trypsin,N,0,Homo sapiens,9606,CIROP,1,1
sp|A0A1B0GTW7|CIROP_HUMAN,713,KPLEVYHGGANFTTQPSK,703,721,18,NFTT,1973.00067520484,-0.91,8.51,trypsin,N,0,Homo sapiens,9606,CIROP,1,1
sp|A0AV02|S12A8_HUMAN,221,LQLLLLFLLAVSTLDFVVGSFTHLDPEHGFIGYSPELLQNNTLPDYSPGESFFTVFGVFFPAATGVMAGFNMGGDLR,182,259,77,NNTL,8353.194115527771,0.56,4.14,trypsin,N,0,Homo sapiens,9606,SLC12A8,1,4
sp|A0AV02|S12A8_HUMAN,561,SEGTQPEGTYGEQLVPELCNQSESSGEDFFLK,542,574,32,NQSE,3504.5514893278296,-0.92,4.05,trypsin,N,0,Homo sapiens,9606,SLC12A8,1,4
sp|A0AV02|S12A8_HUMAN,645,ASPGLHLGSASNFSFFR,634,651,17,NFSF,1793.88491717661,0.16,9.8,trypsin,N,0,Homo sapiens,9606,SLC12A8,1,4
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

## Create Glycan Mass Library (create_glycan_library.py)

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

## Calculate Predicted Intact (Peptide+Glycan mass) Glycopeptide Library (compute_intact_glycopeptide_library.py)

This script computes the intact glycopeptide mass by combining a peptide library with a glycan library. It calculates m/z values for charge states ranging from 2 to a user-defined maximum (default: 6).

- Reads peptide sequences and glycan compositions from CSV files.
- Computes glycopeptide masses by adding peptide and glycan masses.
- Calculates m/z values for multiple charge states.
- Outputs results as a CSV file.

```sh
python script.py -i digested_glycopeptide_library/human_uniprotkb_proteome_UP000005640_AND_revi_2025_01_17_predicted_trypsin_mc_0_N-glycopeptides.csv -g lycan_mass_library/glycan_type_n_linked_byonic_glycan_library.csv -z 6
```

**Note:** It is recommended to use small, targeted glycan databases. Larger glycan databases can significantly increase the size and complexity of the resulting glycopeptide combinations, which may lead to longer processing times and higher computational resource requirements. By focusing on specific glycans of interest, you can streamline the analysis and improve the efficiency of the glycopeptide identification process.

#### Arguments

- `-i`, `--input`: Path to the peptide file (CSV format).
- `-g`, `--glycan`: Path to the glycan file (CSV format).
- `-o`, `--output`: (Optional) Output file path. Defaults to `/predicted_intact_glycopeptide_library/<input_filename>_glycopeptide_library.csv`.
- `-z`, `--charge`: (Optional) Maximum charge state to compute (default: 6).

The script follows these steps:

1. Load peptide and glycan CSV files.
2. Iterate over each peptide and glycan to compute glycopeptide mass.
3. Calculate m/z values for charge states.
4. Write results to a CSV file.

### Example

#### Input Peptide File (peptides.csv):

```csv
ProteinID,Site,Peptide,Start,End,Length,Sequon,PredictedMass,Hydrophobicity,pI,Protease,GlycosylationType,MissedCleavages,Species,TaxonID,GeneName,ProteinEvidence,SequenceVersion
sp|A0A498KFL4|ATGSA_MALDO,163,LNHTMCDAAGLLLFLTAIAEMAR,162,185,23,NHTM,2474.248,0.94,5.32,trypsin,N,0,Malus domestica,3750,AAT1GSA,1,1
```

#### Input Glycan File (glycans.csv):

```csv
glytoucan_ac,byonic,composition,mass
G62765YT,HexNAc(2)Hex(8) % 1702.581333,HexNAc(2)Hex(8),1702.581333
G31852PQ,HexNAc(2)Hex(7) % 1540.528510,HexNAc(2)Hex(7),1540.528510
G41247ZX,HexNAc(2)Hex(6) % 1378.475686,HexNAc(2)Hex(6),1378.475686
```

#### Output File (output.csv):

```csv
ProteinID,Site,glytoucan_ac,composition,glycan_mass,Peptide,Start,End,Length,Sequon,PeptideMass,GlycopeptideMass,Hydrophobicity,pI,Protease,GlycosylationType,MissedCleavages,Species,TaxonID,GeneName,ProteinEvidence,SequenceVersion,z2,z3,z4,z5,z6
sp|A0A087X1C5|CP2D7_HUMAN,416,G62765YT,HexNAc(2)Hex(8),1702.581333,GTTLITNLSSVLK,410,423,13,NLSS,1345.78168089469,3048.36301389469,0.66,8.75,trypsin,N,0,Homo sapiens,9606,CYP2D7,5,1,1525.188782947345,1017.1282806315634,763.0980294736726,610.679878778938,509.06777831578165
sp|A0A087X1C5|CP2D7_HUMAN,416,G31852PQ,HexNAc(2)Hex(7),1540.52851,GTTLITNLSSVLK,410,423,13,NLSS,1345.78168089469,2886.31019089469,0.66,8.75,trypsin,N,0,Homo sapiens,9606,CYP2D7,5,1,1444.162371447345,963.1106729648967,722.5848237236726,578.269314178938,482.0589744824483
sp|A0A087X1C5|CP2D7_HUMAN,416,G41247ZX,HexNAc(2)Hex(6),1378.475686,GTTLITNLSSVLK,410,423,13,NLSS,1345.78168089469,2724.2573668946898,0.66,8.75,trypsin,N,0,Homo sapiens,9606,CYP2D7,5,1,1363.135959447345,909.0930649648966,682.0716177236725,545.858749378938,455.05017048244827
sp|A0A0K2S4Q6|CD3CH_HUMAN,100,G62765YT,HexNAc(2)Hex(8),1702.581333,SDQVIITDHPGDLTFTVTLENLTADDAGK,80,109,29,NLTA,3085.50915394004,4788.09048694004,-0.23,4.05,trypsin,N,0,Homo sapiens,9606,CD300H,1,1,2395.05251947002,1597.0374383133467,1198.02989773501,958.6253733880079,799.0223571566734
sp|A0A0K2S4Q6|CD3CH_HUMAN,100,G31852PQ,HexNAc(2)Hex(7),1540.52851,SDQVIITDHPGDLTFTVTLENLTADDAGK,80,109,29,NLTA,3085.50915394004,4626.03766394004,-0.23,4.05,trypsin,N,0,Homo sapiens,9606,CD300H,1,1,2314.02610797002,1543.01983064668,1157.51669198501,926.214808788008,772.01355332334
sp|A0A0K2S4Q6|CD3CH_HUMAN,100,G41247ZX,HexNAc(2)Hex(6),1378.475686,SDQVIITDHPGDLTFTVTLENLTADDAGK,80,109,29,NLTA,3085.50915394004,4463.98483994004,-0.23,4.05,trypsin,N,0,Homo sapiens,9606,CD300H,1,1,2232.9996959700197,1489.00222264668,1117.00348598501,893.8042439880079,745.00474932334
sp|A0A1B0GTW7|CIROP_HUMAN,333,G62765YT,HexNAc(2)Hex(8),1702.581333,ENCSTR,332,338,6,NCST,708.2860878938,2410.8674208938,-1.75,6.09,trypsin,N,0,Homo sapiens,9606,CIROP,1,1,1206.4409864469,804.6297496312667,603.7241312234501,483.18076017876,402.8185128156333
sp|A0A1B0GTW7|CIROP_HUMAN,333,G31852PQ,HexNAc(2)Hex(7),1540.52851,ENCSTR,332,338,6,NCST,708.2860878938,2248.8145978938,-1.75,6.09,trypsin,N,0,Homo sapiens,9606,CIROP,1,1,1125.4145749469,750.6121419646,563.2109254734501,450.77019557876,375.8097089823
sp|A0A1B0GTW7|CIROP_HUMAN,333,G41247ZX,HexNAc(2)Hex(6),1378.475686,ENCSTR,332,338,6,NCST,708.2860878938,2086.7617738937997,-1.75,6.09,trypsin,N,0,Homo sapiens,9606,CIROP,1,1,1044.3881629469,696.5945339645999,522.69771947345,418.35963077875994,348.8009049822999
sp|A0A1B0GTW7|CIROP_HUMAN,425,G62765YT,HexNAc(2)Hex(8),1702.581333,DSGWYQVNHSAAEELLWGQGSGPEFGLVTTCGTGSSDFFCTGSGLGCHYLHLDK,418,472,54,NHSA,5720.534757730061,7423.116090730061,-0.22,4.56,trypsin,N,0,Homo sapiens,9606,CIROP,1,1,3712.5653213650303,2475.3793062433538,1856.7862986825153,1485.630494146012,1238.1932911216768
```

## Batch Run (batch_glycopeptide_sequence_finder.sh)

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

## Batch Run (batch_compute_intact_glycopeptide_library.sh)

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

## License

This script is released under the MIT License. 

## Acknowledgments

- The Hitchhiker’s Guide to Glycoproteomics.

Oliveira, Tiago, Morten Thaysen-Andersen, Nicolle Packer, and Daniel Kolarich. “The Hitchhiker’s Guide to Glycoproteomics.” Biochemical Society Transactions 49 (July 20, 2021). https://doi.org/10.1042/BST20200879.

- In Silico Platform for Prediction of N-, O- and C-Glycosites in Eukaryotic Protein Sequences.

Chauhan, Jagat Singh, Alka Rao, and Gajendra P. S. Raghava. “In Silico Platform for Prediction of N-, O- and C-Glycosites in Eukaryotic Protein Sequences.” PLoS ONE 8, no. 6 (June 28, 2013): e67008. https://doi.org/10.1371/journal.pone.0067008.

- Large-Scale Identification of N-Linked Intact Glycopeptides in Human Serum Using HILIC Enrichment and Spectral Library Search.

Shu, Qingbo, Mengjie Li, Lian Shu, Zhiwu An, Jifeng Wang, Hao Lv, Ming Yang, et al. “Large-Scale Identification of N-Linked Intact Glycopeptides in Human Serum Using HILIC Enrichment and Spectral Library Search.” Molecular & Cellular Proteomics : MCP 19, no. 4 (April 2020): 672–89. https://doi.org/10.1074/mcp.RA119.001791.

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