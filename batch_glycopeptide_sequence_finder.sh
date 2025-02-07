ls test_proteomes/*.fasta | xargs -I {} -P 4 python glycopeptide_sequence_finder_cmd.py -i "{}" -p trypsin -g N -c 0 -z 2 -v 
