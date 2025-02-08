#!/bin/bash

start_time=$(date +%s)

# define parameters

# cores
cores=4

# Input file directory
input_dir="test_proteomes"

# Protease
protease="trypsin"  # Cleaves after K or R unless followed by P
#protease="chymotrypsin"  # Cleaves after F, W, or Y unless followed by P
#protease="pepsin"  # Cleaves after F, W, or Y
#protease="glu-c"  # Cleaves after E
#protease="lys-c"  # Cleaves after K
#protease="arg-c"  # Cleaves after R
#protease="asp-n"  # Cleaves before Asp (D)
#protease="proteinase-k"  # Cleaves after A, F, I, L, V, W, Y
#protease="all" # all proteases

# missed cleavages
missed_cleavages=0

# Glycosylation type
glycosylation_type="N"
#glycosylation_type="O"
#glycosylation_type="C"

# Charge states
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
echo "LETS GET READY TO FIND SOME GLYCOPEPTIDE SEQUENCES!!!"
sleep 1
echo "Wasting time... loading peptides... loading glycans... loading glycopeptides... testing user's patience..."
sleep 1
echo "Glycopeptide sequence finder is running... maybe... I think... I hope... I'm not sure... I'm just a computer... "
echo "I don't know what I'm doing... I'm just finding... Must find glycopeptides. Must find glycopeptides. Must find glycopeptides."
sleep 1
echo "Stay tuned! We are preparing the ultimate glycopeptide discovery mission!"
sleep 1
echo "Running super advanced algorithms... in the background... somewhere... probably..."
echo "💀 ERROR: WHOA, BRO! Glycosylation type 'X' just crashed the party! ❌"
echo "What are you even doing? Only 'N', 'O', or 'C' are allowed in this twisted code dimension! 🔥"
echo "Check yourself before you wreck yourself—fix it ASAP, or this thing's gonna get REAL weird. 🤖💥"
echo "💀 DEBUG: LOOKING FOR GLYCOSYLATION TYPES... N, O, C ONLY... X? BRUHHHHHH!"
sleep 1
echo "Don't panic, I'm sure we'll find those glycopeptides. Or maybe not. I'm just a computer after all."
sleep 3
echo "LET'S GO ALREADY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
sleep 1
echo "Ope! Here I got glycopeptiding again! 🤣"
sleep 3
echo "WARNING: CPU TEMP THRESHOLD EXCEEDED! 🔥"
echo ""
echo "Traceback (most recent call last):"
echo "  File \"<system>\", line 42, in monitor_cpu"
echo "    check_temperature()"
echo "  File \"<system>\", line 17, in check_temperature"
echo "    raise OverheatedCPUError(\"CPU temperature is too high! 🚨\")"
echo "OverheatedCPUError: CPU Temp: 108°C | Safe Range: 30°C - 80°C 🌡️"
echo ""
echo "💀 ACTION REQUIRED: Cool down your system or prepare for thermal throttling! 🔥"
sleep 1
echo "Rebooting finder. Beeep booop. *sounds*"

# Run the glycopeptide sequence finder script
time ls ${input_dir}/*.fasta | xargs -I {} -P ${cores} python glycopeptide_sequence_finder_cmd.py \
    -i "{}" -p ${protease} -g ${glycosylation_type} -c ${missed_cleavages} -z ${charge_state} -v 

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

# Mess with user
echo "🔥 SYSTEM OVERHEAT DETECTED 🔥"
echo ""
echo "Traceback (most recent call last):"
echo "  File \"<system>\", line 1337, in cool_down"
echo "    prevent_meltdown()"
echo "  File \"<system>\", line 404, in prevent_meltdown"
echo "    raise ThermalRunawayError(\"🔥 CPU TEMP CRITICAL 🔥\")"
echo "ThermalRunawayError: System is now entering emergency cooldown mode... ☠️"
echo ""
echo "[ERROR] CPU TEMP: 105°C 🌡️ (Max Safe: 80°C)"
echo "[SYSTEM] Engaging failsafe protocol... ⏳"
echo "[WARNING] All unsaved data is about to vanish into the void... 🫠"
echo "[REBOOTING] Hold on to your circuits... 💀"
echo ""
echo "🔄 Rebooting in 3... 2... 1... 🚀"

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
                   /               \\        \\
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
sleep 1
echo "Wow! That is a LOT of peptides, buddy! You really know how to slice and dice! 😎🍴"
sleep 1
# print the number of glycopeptide sequences found
echo "The number of $glycosylation_type-glycopeptide sequences found is: $glycopeptide_count. Great job, Glycoscientist! You've done it again!"
echo "Wow! That is a TON of glycopeptides! You are officially the Glycoprotein Master. 🏆"
sleep 1
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
    echo "You really know how to slice and dice those peptides! 😎🍴"
    echo "Total Glycopeptide Sequences Found: $glycopeptide_count"
    echo "You've discovered an impressive number of glycopeptides! Keep it up, Glycoscientist! 🏆"
    echo "---------------------------------------------------------------------------"
    echo "After performing an extensive computational search for glycopeptides in the input protein/proteome FASTA dataset, a total of $peptide_count peptide sequences were identified in the digested_peptide_library."
    echo "These peptides were generated through $protease digestion, with 0 to $missed_cleavages missed cleavages observed during the process."
    echo "Additionally, glycopeptide sequences were detected, yielding a total of $glycopeptide_count glycopeptides in the digested_glycopeptide_library."
    echo "This analysis highlights the successful identification of glycopeptides, which serves as a critical step in the ongoing study of glycoproteomics and the characterization of glycosylation patterns across various proteins."
    echo "The findings provide valuable insight into the glycosylation landscape of the proteome by calculating the mass spectrometry properties of peptides and glycopeptides, laying the groundwork for more advanced analyses in future research."
    echo "---------------------------------------------------------------------------"
    echo "Batch Run Summary:"
    echo "Total Elapsed Time: $elapsed_time_mins minutes"
    echo "Protease Used: $protease"
    echo "Missed Cleavages: $missed_cleavages"
    echo "Glycosylation Type: $glycosylation_type-Glycosylation"
    echo "Number of Peptide Sequences: $peptide_count"
    echo "Number of $glycosylation_type-Glycopeptide Sequences: $glycopeptide_count"
    echo "---------------------------------------------------------------------------"
    echo "Thank you for using the Glycopeptide Sequence Finder! Please come back when you're ready to discover more glycopeptides!"
    echo "$ascii_glycopeptide2"
    echo "---------------------------------------------------------------------------"
} > "$summary_report"

# Notify user of the report
echo "A summary of the batch run has been saved to $summary_report."


