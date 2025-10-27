#!/usr/bin/env python3

import pandas as pd
import numpy as np
import glob, os
import sys

if __name__ == "__main__":
    panaroo = sys.argv[1]
    results = sys.argv[2]
    pyseer = sys.argv[3]

# Open panaroo initial cd-hit output
with open(f"{panaroo}/combined_protein_cdhit_out.txt.clstr", "r") as f:
    lines = f.read().strip().split("\n")

# Split into clusters based on lines starting with ">Cluster"
clusters = []
current_cluster = []
current_name = None

for line in lines:
    if line.startswith(">Cluster"):
        # Save previous cluster if it exists
        if current_cluster:
            df = pd.DataFrame(current_cluster, columns=["raw"])
            df["Cluster"] = current_name
            clusters.append(df)
            current_cluster = []
        current_name = line.strip().split()[1]  # e.g., "0", "1", "2"
    else:
        current_cluster.append([line.strip()])

# Add the last cluster
if current_cluster:
    df = pd.DataFrame(current_cluster, columns=["raw"])
    df["Cluster"] = current_name
    clusters.append(df)

# Combine all clusters into a single dataframe (optional)
combined_df = pd.concat(clusters, ignore_index=True)

# Split to get the correct reference
ident = combined_df['raw'].str.split(" ", expand = True)

ident.columns = ['Length', 'clustering_id', 'Reference', 'Threshold']

# Keep only reference proteins
refs = ident[ident['Reference'].str.contains(r"\*")]

# Reformat reference
refs['clustering_id'] = refs['clustering_id'].str.replace('>','')
refs['clustering_id'] = refs['clustering_id'].str.replace('...','')

# Subset for clustering id
refs_sub = refs[['clustering_id']]

# refs_sub.to_csv("refs_sub_no_refind_changes.csv", index = None)

# refs_sub = pd.read_csv("refs_sub_no_refind_changes.csv", sep = ",")

# Read in gene data
gene_data = pd.read_csv(f"{panaroo}/gene_data.csv", sep = ",")

# Remove the sequences
gene_data_sub = gene_data.drop(["prot_sequence", "dna_sequence"], axis=1)

# Merge with reference clusters to assign function to references
merge = refs_sub.merge(gene_data_sub, how = "left", on = "clustering_id")

# Read in gene presence absence
pa_df = pd.read_csv(f"{panaroo}/gene_presence_absence.csv", sep = ",")

# Obtain the gff_files and annotation ids
col = merge["gff_file"].iloc[0]
val = merge["annotation_id"].iloc[0]

# Map gff columns
col_mappings = {col: dict(zip(pa_df[col], pa_df["Gene"])) for col in pa_df.columns if col != "Gene"}

# Lookup using the mapping
merge["variant"] = merge.apply(lambda row: col_mappings[row["gff_file"]].get(row["annotation_id"], None), axis=1)

# Save
merge.to_csv(f"{results}/Column_mappings_no_refind_changes.csv", index = None)

# Read in pyseer gene results
pyseer_genes = pd.read_csv(f"{pyseer}/LasR_pyseer_gwas_functions_high.csv")

# Merge results with functions
pyseer_merge = pyseer_genes.merge(merge, how = "left", on = "variant")

#Save
pyseer_merge.to_csv(f"{results}/pyseer_merge_no_refind_changes.csv", index = None)