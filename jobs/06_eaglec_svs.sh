#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=gpu4_medium,gpu8_medium
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:00:00
#SBATCH --mem=16G
#SBATCH --mail-type=ALL

# This job uses EagleC to call SVs on .mcool files.

PREFIX=$1

MCOOL_DIR=/gpfs/home/$USER/labspace/hic/data/mcool
MCOOL_PATH=$MCOOL_DIR/$PREFIX.mcool
OUT_DIR=/gpfs/home/$USER/labspace/hic/temp/eaglec/$PREFIX
mkdir -p $OUT_DIR
OUT_PREFIX=$OUT_DIR/$PREFIX

module load cuda/11.8
module load condaenvs/new/EagleC

predictSV --hic-5k $MCOOL_PATH::/resolutions/5000 \
	--hic-10k $MCOOL_PATH::/resolutions/10000 \
	--hic-50k $MCOOL_PATH::/resolutions/50000 \
	-O $OUT_PREFIX -g hg38 --balance-type Raw --output-format NeoLoopFinder \
	--prob-cutoff-5k 0.8 --prob-cutoff-10k 0.8 --prob-cutoff-50k 0.99999

cp $OUT_DIR/${PREFIX}.CNN_SVs.5K_combined.txt ~/labspace/hic/data/eaglec/

module unload cuda/11.8
module unload condaenvs/new/EagleC

