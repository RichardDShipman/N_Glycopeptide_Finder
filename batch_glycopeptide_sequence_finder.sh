#!/bin/bash
start_time=$(date +%s)

# define parameters
# cores
cores=4

# Input file directory (-i)
input_dir="test_proteomes"

# Protease (-p)
protease="trypsin"  # Cleaves after K or R unless followed by P
#protease="chymotrypsin"  # Cleaves after F, W, or Y unless followed by P
#protease="pepsin"  # Cleaves after F, W, or Y
#protease="glu-c"  # Cleaves after E
#protease="lys-c"  # Cleaves after K
#protease="arg-c"  # Cleaves after R
#protease="asp-n"  # Cleaves before Asp (D)
#protease="proteinase-k"  # Cleaves after A, F, I, L, V, W, Y
#protease="all" # all proteases

# missed cleavages and max peptide length filter (-c, -m)
missed_cleavages=0
max_peptide_length=25

# Glycosylation type (-g)
glycosylation_type="N"
#glycosylation_type="O"
#glycosylation_type="C"

# Charge states (-z)
charge_state=2

# welcome message
ascii_glycopeptide1="
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  🔬 F I N D I N G   G L Y C O P E P T I D E   S E Q U E N C E S ! ! ! 🔬
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                 GLYCO-
      H-N-C-C-O--PEPTIDE--N-C-C-O-H-N-C-C-O-H
          |      SEQUENCE   |         |
          R      FINDER     R         R
         /                  \\          \\
        N-Glycan             O-Glycan   C-Glycan
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  📊 P R E D I C T I N G   G L Y C O P E P T I D E   P R O P E R T I E S ! ! ! 📊
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Welcome message
echo "Welcome to the Batch Run of the Glycopeptide Sequence Finder!"
echo "$ascii_glycopeptide1"
echo "Starting glycopeptide sequence finder..."
echo "Please wait while the glycopeptide sequence finder is running..."

# Run the glycopeptide sequence finder script
time ls ${input_dir}/*.fasta | xargs -I {} -P ${cores} python glycopeptide_sequence_finder_cmd.py \
    -i "{}" -p ${protease} -g ${glycosylation_type} -c ${missed_cleavages} -z ${charge_state} -m ${max_peptide_length} -v 

echo "Digested and tasted the glycoproteome. Yummy! 🍽️"

# calculate elapsed time 
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
elapsed_time_mins=$(echo "scale=2; $elapsed_time / 60" | bc)
echo "Total elapsed time: $elapsed_time_mins minutes."

# Count the number of peptide sequences, row entries, in all csv files in the digested_peptide_library directory, ignoring headers
peptide_count=$(tail -n +2 digested_peptide_library/*.csv | wc -l)

# Count the number of glycopeptide sequences, row entries, in all csv files in the digested_glycopeptide_library directory, ignoring headers
glycopeptide_count=$(tail -n +2 digested_glycopeptide_library/*.csv | wc -l)

ascii_glycopeptide2="
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    🔬 F O U N D   G L Y C O P E P T I D E   S E Q U E N C E S ! ! ! 🔬
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━                                   
    I D E N T I F I E D: $glycopeptide_count $protease-digested $glycosylation_type-glycopeptide sequences.                                      
    I D E N T I F I E D: $peptide_count $protease-digested peptide sequences.
        ╭────────────────────────────────────────────╮
                        GLYCO-
            H-N-C-C-O--PEPTIDE--N-C-C-O-H-N-C-C-O-H
                    |  SEQUENCE   |         |
                    R   FINDER    R         R
                   /               \\         \\
                  N-Glycan          O-Glycan   C-Glycan
        ╰────────────────────────────────────────────╯
                    ╲
                    💬 \"I crave more glycopeptides! Bring me more proteomes!\"
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    📊 P R E D I C T E D   G L Y C O P E P T I D E   P R O P E R T I E S ! ! ! 📊
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# outro message
echo "$ascii_glycopeptide2"
sleep 2

# print the number of peptide sequences found
echo "The number of $protease digested peptide sequences found is: $peptide_count."
echo "Wow! That is a LOT of peptides, buddy! You really know how to slice and dice! 😎🍴"
# print the number of glycopeptide sequences found
echo "The number of $glycosylation_type-glycopeptide sequences found is: $glycopeptide_count. Great job, Glycoscientist! You've done it again!"
echo "Wow! That is a TON of glycopeptides! You are officially the Glycoprotein Master. 🏆"
echo "You are a glycopeptide finding MACHINE! Keep up the epic work, my friend! 💪🚀"
echo "Glycopeptide sequence finder has finished.  🎉"
echo "But... Wait! Please feed me more proteome data! I'm hungry for more glycopeptides! 🍽️ Bring on the next batch!"
sleep 1
echo "Results written to summary_batch_run.txt."

# Create a summary report
summary_report="summary_batch_run.txt"

# Write the summary to the file
{
    echo "$ascii_glycopeptide1"
    echo "==================== Glycopeptide Sequence Finder Summary ===================="
    echo "After some finding and searching for glycopeptides in input protein/proteome FASTA..."
    echo "Total Peptide Sequences Found: $peptide_count"
    echo "Peptides generated using $protease digestion with $missed_cleavages missed cleavages."
    echo "Total Glycopeptide Sequences Found: $glycopeptide_count"
    echo "---------------------------------------------------------------------------"
    echo "Batch Run Summary:"
    echo "Total Elapsed Time: $elapsed_time_mins minutes"
    echo "Protease Used: $protease"
    echo "Missed Cleavages: $missed_cleavages"
    echo "Glycosylation Type: $glycosylation_type-Glycosylation"
    echo "Number of Peptide Sequences: $peptide_count"
    echo "Number of $glycosylation_type-Glycopeptide Sequences: $glycopeptide_count"
    echo "---------------------------------------------------------------------------"
    echo "Thank you for using the Glycopeptide Sequence Finder!" 
    echo "Please come back when you're ready to discover more glycopeptides!"
    echo "---------------------------------------------------------------------------"
} > "$summary_report"

# Notify user of the report
echo "A summary of the batch run has been saved to $summary_report."


