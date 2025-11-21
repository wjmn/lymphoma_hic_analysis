#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=5-00:00:00
#SBATCH --partition=cpu_medium,fn_medium
#SBATCH --mem=64GB
#SBATCH --mail-type=ALL


# This job uses hitad from TADLib to call TADs

# Source conda environment with TADLib installed
source ~/labspace/hic/scripts/conda_env.sh

PARENT=/gpfs/home/$USER/labspace/hic/data
RES=25000

PREFIX=$1

HIC_FILE=$PARENT/hic/${PREFIX}_inter_30.hic
MCOOL_ABSOLUTE=$PARENT/mcool/$PREFIX.mcool

# Create metadata
METADATA=/gpfs/home/$USER/labspace/hic/temp/tad_metadata/$PREFIX.$RES.metadata.txt
python3 ~/labspace/hic/scripts/make_tad_metadata.py $MCOOL_ABSOLUTE::/resolutions/$RES $METADATA $RES

# Run HiTAD
OUTPUT=$PARENT/tads/$RES/$PREFIX.txt
hitad -O $OUTPUT -d $METADATA -p 24 -W weight
