#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob, os
import sys

if __name__ == "__main__":
    data = sys.argv[1]
    results = sys.argv[2]

# Read in genome list
all_genomes = pd.read_csv(f"{data}/00_file_list.txt", sep = "\t", header = None, names = ['Genome'])

# Read in functions
functs = pd.read_csv(f"{data}/PA1430_functions_incomplete.csv", sep = ",")

# Merge to account for all genomes
merge = all_genomes.merge(functs, how = "left", on = "Genome")

# Replace "No function" with 1
merge['Mut_Status'] = merge['Mut_Status'].str.replace('No function', '1')

# Replace "Functional" with 0
merge['Mut_Status'] = merge['Mut_Status'].str.replace('Functional', '0')

# Fill all not accounted for with NA
merge_filled = merge.fillna('NA')

# Save
merge_filled.to_csv(f"{results}/PA_1430_phenotypes.tsv", sep = "\t", index = None)