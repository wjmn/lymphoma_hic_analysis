#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --partition=cpu_short,fn_short,cpu_dev,gpu4_dev,gpu4_short
#SBATCH --mem=16GB
#SBATCH --mail-type=ALL

# This job uses NeoLoopFinder to call neoloops. This step occurs after
# breakpoints have been manually reviewed.

module load condaenvs/gpu/neoloop

sample=$1
prob=90
neotad_res=25000

input=/gpfs/home/$USER/labspace/hic/data/mcool/${sample}.mcool
breakpoints=/gpfs/home/$USER/labspace/hic/data/curated_breakpoints/${sample}_curated_breakpoints.tsv
assembly_prefix=/gpfs/home/$USER/labspace/hic/data/assemblies/${sample}
assembly_out=/gpfs/home/$USER/labspace/hic/data/assemblies/${sample}.assemblies.txt
neoloop_out=/gpfs/home/$USER/labspace/hic/data/neoloops/${prob}/${sample}.neoloops.bedpe
neotad_out=/gpfs/home/$USER/labspace/hic/data/neotads/${neotad_res}/${sample}_${neotad_res}.neotad.txt

assemble-complexSVs -O ${assembly_prefix} -B ${breakpoints} \
    --balance-type CNV --nproc 8 \
    -H ${input}::/resolutions/25000 ${input}::/resolutions/10000 ${input}::/resolutions/5000

neoloop-caller -O ${neoloop_out} --assembly ${assembly_out} \
    --balance-type CNV --prob 0.${prob} --nproc 8 \
    -H ${input}::/resolutions/25000 ${input}::/resolutions/10000 ${input}::/resolutions/5000

neotad-caller -O ${neotad_out} --assembly ${assembly_out} \
    --balance-type CNV --nproc 8 \
    -H ${input}::/resolutions/${neotad_res}
