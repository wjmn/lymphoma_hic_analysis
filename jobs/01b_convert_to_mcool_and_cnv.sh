#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --partition=cpu_short,fn_short,gpu4_short,gpu8_short
#SBATCH --mem=16GB
#SBATCH --mail-type=ALL

sample=$1

hic=/gpfs/home/$USER/hic/data/hic/${sample}_inter_30.hic
input=/gpfs/home/$USER/hic/data/mcool/${sample}.mcool

# Convert
# Load conda environment with HiCExplorer installed
source ~/labspace/hic/scripts/conda_env.sh
hicConvertFormat -m $hic --inputFormat hic --outputFormat cool -o $input

# Prefix and balance (weight)
module load condaenvs/new/cooler
python3 ~/labspace/hic/scripts/add_prefix_to_mcool.py $input
for res in 5000 10000 25000 50000; do
    cooler balance ${input}::/resolutions/${res} -f
done
module unload condaenvs/new/cooler

# CNV correction (sweight)
module load condaenvs/gpu/neoloop
for res in 5000 10000 25000 50000; do
    profile=/gpfs/home/$USER/hic/data/cnv_profile/${res}/${sample}_CNV_profile_${res}.bedgraph
    output=/gpfs/home/$USER/hic/data/cnv_segment/${res}/${sample}_CNV_segment_${res}.bedgraph
    calculate-cnv -H ${input}::/resolutions/${res} -g hg38 -e Arima --output ${profile}
    segment-cnv --cnv-file ${profile} --binsize ${res} --ploidy 2 --output ${output}
    correct-cnv -H ${input}::/resolutions/${res} --cnv-file ${output} -f
done

# Chromosomal arm-level CNV calculation
res=500000
profile=/gpfs/data/$USER/hic/data/cnv_profile/${res}/${sample}_CNV_profile_${res}.bedgraph
output=/gpfs/data/$USER/hic/data/cnv_segment/${res}/${sample}_CNV_segment_${res}.bedgraph
calculate-cnv -H ${input}::/resolutions/${res} -g hg38 -e Arima --output ${profile}
segment-cnv --cnv-file ${profile} --binsize ${res} --ploidy 2 --output ${output}
