#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --partition=fn_short,cpu_short,gpu4_short,gpu8_short,gpu4_dev,cpu_dev
#SBATCH --mem=8GB

# This job uses the eigenvector utility from Juicer Tools to call compartments
# from a processed Hi-C file.

module load java/1.8

# Resolution for calling compartments
RES=100000

# Sample file prefix
PREFIX=$1

JUICER_PATH=/gpfs/scratch/$USER/juicer/juicer-1.6/SLURM/scripts

PARENT=/gpfs/home/$USER/labspace/hic/data
HIC=$PARENT/hic/${PREFIX}_inter_30.hic

OUTDIR=/gpfs/home/$USER/labspace/hic/data/compartments_juicer_raw/$RES/$PREFIX/
mkdir -p $OUTDIR

for CHR in {1..22}; do
    $JUICER_PATH/juicer_tools -p eigenvector KR $HIC $CHR BP $RES $OUTDIR/$PREFIX.$RES.$CHR.txt
done

# Process compartments to follow convention positive score = active chromatin state
BIGWIG_H3K27AC=$2 # H3K27ac data to compare against to flip compartments
python3 /gpfs/home/$USER/labspace/hic/scripts/process_chr_compartments.py \
    $RES \
    $BIGWIG_H3K27AC \
    "/gpfs/home/$USER/labspace/hic/temp/compartments_juicer_raw/$RES/$PREFIX" \
    $PREFIX \
    "/gpfs/home/$USER/labspace/hic/data/compartments/$RES/$PREFIX.$RES.bedgraph"
