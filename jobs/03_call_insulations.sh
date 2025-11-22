#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --partition=fn_short,cpu_short,gpu4_short,gpu8_short,gpu4_dev,cpu_dev
#SBATCH --mem=8GB

# This job is a wrapper over using cooltools to calculate insulation scores (at a set resolution/window)
# from a processed Hi-C file.


PREFIX=$1
RESOLUTION=10000
WINDOW=100000

OUT_DIR=~/labspace/hic/data/insulations/res_${RESOLUTION}_window_${WINDOW}/
mkdir -p $OUT_DIR

source ~/labspace/hic/scripts/conda_env.sh

python3 ~/labspace/hic/scripts/get_insulations.py $PREFIX $RESOLUTION $WINDOW