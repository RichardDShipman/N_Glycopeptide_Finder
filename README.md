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
pip install argparse biopython pyteomics
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

## Batch Run

To process multiple FASTA files in parallel using all proteases, run the following command:

```sh
./batch_glycopeptide_sequence_finder.sh
```

Parameters can be adjusted in the shell script.

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

- The Hitchhiker’s Guide to Glycoproteomics

Oliveira, Tiago, Morten Thaysen-Andersen, Nicolle Packer, and Daniel Kolarich. “The Hitchhiker’s Guide to Glycoproteomics.” Biochemical Society Transactions 49 (July 20, 2021). https://doi.org/10.1042/BST20200879.

- In Silico Platform for Prediction of N-, O- and C-Glycosites in Eukaryotic Protein Sequences.

Chauhan, Jagat Singh, Alka Rao, and Gajendra P. S. Raghava. “In Silico Platform for Prediction of N-, O- and C-Glycosites in Eukaryotic Protein Sequences.” PLoS ONE 8, no. 6 (June 28, 2013): e67008. https://doi.org/10.1371/journal.pone.0067008.

- BioPython for handling FASTA files.

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., … others. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423.

- Pyteomics for accurate peptide mass calculations.

Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6

Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717

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
