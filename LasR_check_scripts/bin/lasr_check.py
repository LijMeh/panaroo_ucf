#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob, os
import sys

if __name__ == "__main__":
    data = sys.argv[1]
    panaroo = sys.argv[2]
    results = sys.argv[3]
    variant = sys.argv[4]

# Read in gene presence absence table
table = pd.read_csv(f"{panaroo}/gene_presence_absence.csv", sep = ",")

# Only keep the row with the variant
filtered_df = table[table['Gene'].str.contains(f'{variant}')]

# Drop metadata columns
filtered_df_sub = filtered_df.drop(["Gene", "Non-unique Gene name", "Annotation"], axis=1)

# Transform
filtered_df_sub_T = filtered_df_sub.T.reset_index().rename(columns={"index": "Genome"})

# Rename
filtered_df_sub_T.columns = ['Genome', 'LasR_panaroo']

# Drop NA
filtered_df_sub_T_na = filtered_df_sub_T.dropna(subset=["LasR_panaroo"])

# Read in function data
lasr_funct = pd.read_csv(f"{data}/PA1430_all_functions_incomplete.csv", sep = ",")

# Merge panaroo output with functions data
merge = filtered_df_sub_T_na.merge(lasr_funct, how = "left", on = "Genome")

# Save
merge.to_csv(f"{results}/LasR_functions_panaroo.csv", index = None)
