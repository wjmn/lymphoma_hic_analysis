#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_dev,fn_short,cpu_short,gpu4_short,gpu8_short,gpu4_dev
#SBATCH --time=01:00:00
#SBATCH --mem=8GB
#SBATCH --mail-type=ALL

# This job is a wrapper over a script to merge hic_breakfinder and EagleC SV
# calls.

PREFIX=$1

HIC=~/labspace/hic/data/hic/${PREFIX}_inter_30.hic
QC=~/labspace/hic/data/qc_deep/${PREFIX}_v1.3_Arima_QC_deep.txt
BREAKFINDER=~/labspace/hic/data/hic_breakfinder/${PREFIX}.breaks.bedpe
EAGLEC=~/labspace/hic/data/eaglec/${PREFIX}.CNN_SVs.5K_combined.txt

OUT_DIR=~/labspace/hic/transfer/merged_breakpoints
OUTPUT=$OUT_DIR/${PREFIX}.merged.breaks.bedpe

mkdir -p $OUT_DIR

source ~/labspace/hic/scripts/conda_env.sh

python3 ~/labspace/hic/scripts/merge_breakfinder_eaglec.py $BREAKFINDER $EAGLEC $OUTPUT
