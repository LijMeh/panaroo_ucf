#!/bin/bash
#SBATCH --job-name=pyseer_and_lasR_check
#SBATCH --account=gts-sbrown365-paid
#SBATCH -N1 --ntasks-per-node=24
#SBATCH --mem-per-cpu=4G
#SBATCH --time=30:00:00
#SBATCH -q inferno
#SBATCH -o Report-%x.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=iirby3@gatech.edu

# Set paths and parameters
file_dir="/storage/home/hcoda1/2/iirby3/scratch/plasmid_genomes/dp_comparison/dp_ident/all_host_seq/P_aeruginosa/Mutant_calling/panaroo_ucf"
data_dir="${file_dir}/data"
panaroo_dir="${file_dir}/panaroo_100_results_no_collapse_families"
lrt_pvalue_cutoff="0.05"
lrt_pvalue_cutoff_high="0.00001"

export PATH="${file_dir}/scripts/bin:$PATH"

mkdir -p ${file_dir}/Results
results_dir="${file_dir}/Results"

module load mamba/1.4.9

# Create presence absence input
pres_abs_input.py "$data_dir" "$results_dir"

mamba activate pyseer

mkdir -p ${results_dir}/Pyseer

# Create distance matrix
python /storage/home/hcoda1/2/iirby3/scratch/pyseer/scripts/phylogeny_distance.py \
    --lmm ${data_dir}/iqtree_new_1500_5000_max_pruned.contree \
    > ${results_dir}/Pyseer/phylogeny_distances.tsv

# Run pyseer
python /storage/home/hcoda1/2/iirby3/scratch/pyseer/pyseer-runner.py \
    --lmm \
    --phenotypes ${results_dir}/PA_1430_phenotypes.tsv \
    --pres ${panaroo_dir}/gene_presence_absence.Rtab \
    --similarity ${results_dir}/Pyseer/phylogeny_distances.tsv \
    --print-samples \
    > ${results_dir}/Pyseer/LasR_pyseer_gwas.txt

# Reformat pyseer output
pyseer_functions.py "$file_dir" "$panaroo_dir" "$results_dir" "$lrt_pvalue_cutoff" "$lrt_pvalue_cutoff_high"

# Get LasR variant of interest
value=$(cat ${results_dir}/LasR_variant.txt)

# Check the LasR variants called by panaroo
lasr_check.py "$data_dir" "$panaroo_dir" "$results_dir" "$value"

mkdir -p ${results_dir}/function_mapping

function_mapping_dir="${results_dir}/function_mapping"

# Assign reference genes to pyseer genes of interest
reformat_cd_hit.py "$panaroo_dir" "$function_mapping_dir" "${results_dir}/Pyseer"