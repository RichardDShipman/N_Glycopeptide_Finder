ls test_proteomes/*.fasta | xargs -I {} -P 4 python n_glycopeptide_sequence_finder_cmd.py -i "{}" -p all -c 0
