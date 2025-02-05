ls predicted_intact_glycopeptide_library/*.csv | xargs -I {} -P 4 python compute_glycopeptide_hydrophobicity.py -i "{}" 
