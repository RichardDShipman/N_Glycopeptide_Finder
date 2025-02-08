#!/bin/bash

start_time=$(date +%s)

echo "Starting glycopeptide sequence finder..."

time ls test_proteomes/*.fasta | xargs -I {} -P 4 python glycopeptide_sequence_finder_cmd.py -i "{}" -p trypsin -g N -c 0 -z 2 -v 
#time ls test_proteomes/*.fasta | xargs -I {} -P 4 python glycopeptide_sequence_finder_cmd.py -i "{}" -p trypsin -g O -c 0 -z 2 -v 
time ls test_proteomes/*.fasta | xargs -I {} -P 4 python glycopeptide_sequence_finder_cmd.py -i "{}" -p trypsin -g C -c 0 -z 2 -v 

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

elapsed_time_mins=$(echo "scale=2; $elapsed_time / 60" | bc)
echo "Total elapsed time: $elapsed_time_mins minutes"
