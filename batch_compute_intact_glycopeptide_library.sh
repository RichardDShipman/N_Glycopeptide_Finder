ls digested_glycopeptide_library/*.csv | xargs -I {} -P 4 python compute_intact_glycopeptide_library.py -i "{}" -g glycan_mass_library/test_glycan_library.csv -z 6
