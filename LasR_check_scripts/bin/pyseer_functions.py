#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob, os
import sys

if __name__ == "__main__":
    path = sys.argv[1]
    panaroo = sys.argv[2]
    results = sys.argv[3]
    cutoff = sys.argv[4]
    high_cutoff = sys.argv[5]

# Read in pyseer table
pyseer_table = pd.read_csv(f"{results}/Pyseer/LasR_pyseer_gwas.txt", sep = "\t")

# Get counts of presence and absence of features
pyseer_table['k-count'] = pyseer_table['k-samples'].fillna('').apply(lambda x: len(x.split(',')) if x else 0)
pyseer_table['nk-count'] = pyseer_table['nk-samples'].fillna('').apply(lambda x: len(x.split(',')) if x else 0)

# Drop the original columns
pyseer_table = pyseer_table.drop(columns=['k-samples', 'nk-samples'])

# Read in gene presence absence table from panaroo
functions_table = pd.read_csv(f"{panaroo}/gene_presence_absence_roary.csv", sep = ",")

# Only keep functions
functions_table_sub = functions_table[['Gene', 'Non-unique Gene name', 'Annotation']]
functions_table_sub.rename(columns = {'Gene':'variant'}, inplace = True)

# Identify the functions per variant
merged = pyseer_table.merge(functions_table_sub, how = "left", on = "variant")

# Sort by lrt-pvalue
merged_sort = merged.sort_values(by='lrt-pvalue', ascending=True, inplace=False)

# Only keep values greater than cutoff
merged_sort_filtered = merged_sort[merged_sort['lrt-pvalue'] <= float(cutoff)]

# Save all
merged_sort_filtered.to_csv(f"{results}/Pyseer/LasR_pyseer_gwas_functions.csv", index = None)

# Only keep values greater than high cutoff for further analysis
merged_sort_filtered_high = merged_sort_filtered[merged_sort_filtered['lrt-pvalue'] <= float(high_cutoff)]

# Save subset
merged_sort_filtered_high.to_csv(f"{results}/Pyseer/LasR_pyseer_gwas_functions_high.csv", index = None)

# Only keep variants of interest
merged_sort_filtered_high_sub = merged_sort_filtered_high[["variant"]]

# Save variants
merged_sort_filtered_high_sub.to_csv(f"{results}/Pyseer/LasR_pyseer_gwas_genes_of_interest.txt", sep = "\t", index = None, header = None)

# Identify variant of interest CHECK MANUALLY
value = merged_sort_filtered_high_sub.iloc[0, 0]

# Save to file
with open(f"{results}/LasR_variant.txt", "w") as f:
    f.write(str(value) + "\n")