#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=25
#SBATCH --partition=fn_medium,cpu_medium
#SBATCH --time=4-00:00:00
#SBATCH --mem=88GB

# This job processes raw FASTQ files into Hi-C files using the Arima SV
# Pipeline. Candidate SVs are also called with hic_breakfinder through the
# pipeline.

PREFIX=$1

# Setup up analysis folder
RAW_DIR=/gpfs/scratch/$USER/lymphoma
BOUND_DIR=/gpfs/scratch/$USER/lymphoma_bound/$PREFIX
FASTQ_DIR=$BOUND_DIR/fastq
F1=$FASTQ_DIR/${PREFIX}_R1.fastq.gz
F2=$FASTQ_DIR/${PREFIX}_R2.fastq.gz

mkdir -p $FASTQ_DIR
mkdir -p $BOUND_DIR/output

cp $RAW_DIR/${PREFIX}_R1.fastq.gz $F1
cp $RAW_DIR/${PREFIX}_R2.fastq.gz $F2

LABSPACE=/gpfs/home/$USER/labspace/hic/data

# Run Arima SV Pipeline
module load singularity/3.7.1

SINGULARITY_DIR=/gpfs/home/$USER/labspace/hic/singularity

singularity exec -B $BOUND_DIR:/FFPE/$USER \
    $SINGULARITY_DIR/Arima-SV-Pipeline-singularity-v1.3.sif \
    bash /FFPE/Arima-SV-Pipeline-v1.3.sh \
    -I /FFPE/$USER/fastq/${PREFIX}_R1.fastq.gz,/FFPE/$USER/fastq/${PREFIX}_R2.fastq.gz \
    -o /FFPE/$USER/output/ \
    -p $PREFIX \
    -t 24 \
    -W 1 \
    -B 1 \
    -J 1 \
    -H 0 \
     -a /root/anaconda3/bin/bowtie2 \
    -b /usr/local/bin/ \
    -w /FFPE/HiCUP-0.8.0 \
    -j /FFPE/juicer-1.6/ \
    -r /FFPE/Arima_files/reference/hg38/hg38.fa \
    -s /FFPE/Arima_files/Juicer/hg38.chrom.sizes \
    -c /FFPE/Arima_files/Juicer/hg38_GATC_GANTC.txt \
    -x /FFPE/Arima_files/reference/hg38/hg38 \
    -d /FFPE/Arima_files/HiCUP/Digest_hg38_Arima.txt \
    -e /FFPE/Arima_files/hic_breakfinder/intra_expect_100kb.hg38.txt \
    -E /FFPE/Arima_files/hic_breakfinder/inter_expect_1Mb.hg38.txt
