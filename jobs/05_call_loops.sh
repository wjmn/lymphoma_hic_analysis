#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=fn_short,cpu_short,cpu_dev,gpu4_dev,gpu4_short
#SBATCH --time=01:00:00
#SBATCH --mem=16GB
#SBATCH --mail-type=ALL

# This job calls loops from the Hi-C file using HiCExplorer.

PREFIX=$1

HIC_DIR=~/labspace/hic/data/hic

COOL_DIR=~/labspace/hic/data/cool
OUT_DIR=~/labspace/hic/data/hicexplorer_loops
mkdir -p $COOL_DIR
mkdir -p $OUT_DIR

# Load conda environment with HiCExplorer installed
source ~/labspace/hic/scripts/conda_env.sh

for RES in 10000; do
    mkdir -p $COOL_DIR/$RES
    mkdir -p $OUT_DIR/$RES

    hicConvertFormat \
        --inputFormat hic \
        --outputFormat cool \
        -r $RES \
        -m "${HIC_DIR}/${PREFIX}_inter_30.hic" \
        -o "${COOL_DIR}/${RES}/${PREFIX}.cool"

    hicDetectLoops \
        --matrix "${COOL_DIR}/${RES}/${PREFIX}_${RES}.cool" \
        --outFileName "${OUT_DIR}/${RES}/${PREFIX}_hicexplorer_loops_${RES}.bedpe" \
        --maxLoopDistance 5000000 \
        --peakInteractionsThreshold 10 \
        --windowSize 10 \
        --peakWidth 6 \
        --pValue 1 \
        --pValuePreselection 1

done
