#!/usr/bin/env python3
#
# Usage:
#     python3 process_chr_compartments.py resolution bw_file chr_dir sample_id output_path

# Juicer Tools calculates compartments for each chromosome separately and
# arbitrarily assigns a sign.  This script compare compartment scores for a
# chromosome against H3K27ac reference data and flips the sign of compartment
# scores for the whole chromosome to follow the convention of positive
# compartments = active chromatin states.

import sys
import pandas as pd
import pyBigWig
import numpy as np
import scipy.stats

# https://github.com/wjmn/hicdash
from hicdash.constants import *
from hicdash.utilities import *

resolution = int(sys.argv[1])
bw_file = sys.argv[2]
chr_dir = sys.argv[3]
sample_id = sys.argv[4]
output_path = sys.argv[5]

# Load data files - mask segdups and blacklist
SEGDUPS_FILE="~/labspace/hic/annotations/segdups/GRCh38_SegmentalDuplications_GRCh38_segdups.bed"
BLACKLIST_FILE="~/labspace/hic/annotations/blacklist/hg38-blacklist.v2.bed"

segdups = pd.read_csv(SEGDUPS_FILE, sep="\s+", skiprows=1, names=["chr", "start", "end"])
blacklist = pd.read_csv(BLACKLIST_FILE, sep="\t", names=["chr", "start", "end", "reason"])
all_blacklist = pd.concat([segdups, blacklist.iloc[:, :3]], ignore_index=True)

bw = pyBigWig.open(bw_file)

def make_genomic_bins(chr, bin_size=100000):
    start = 0
    end = CHROM_SIZES[chr]
    bins = []
    while start < end:
        bins.append((chr, start, start+bin_size-1))
        start += bin_size
    return bins

results = []
all_bins = []
for chr in CHROMS:

    # Ignore chrX and chrY
    if chr == "chrX" or chr == "chrY":
        continue

    eigenvalues = pd.read_csv(f"{chr_dir}/{sample_id}.{resolution}.{chr_unprefix(chr)}.txt", sep="\s+", header=None).values.flatten()
    bins = make_genomic_bins(chr, bin_size=resolution)
    bigwig_values = np.array(bw.stats(chr, 0, CHROM_SIZES[chr], type="mean", nBins=len(bins)), dtype=np.float64)
    assert (len(eigenvalues) == len(bins))
    assert (len(bigwig_values) == len(bins))

    # Make an array of 1s corresponding to regions blacklisted
    chr_blacklist = all_blacklist[all_blacklist["chr"] == chr]
    blacklist_mask = np.zeros(len(bins), dtype=np.int8)
    # For every region in chr_blacklist, set the corresponding range of bins to 1
    for i, row in chr_blacklist.iterrows():
        start = row["start"] // resolution
        end = row["end"] // resolution
        blacklist_mask[start:end] = 1

    nan_bigwig = np.isnan(bigwig_values)
    nan_eigen = np.isnan(eigenvalues)
    mask = ~(nan_bigwig | nan_eigen) & (blacklist_mask == 0)
    corr_noflip = scipy.stats.pearsonr(eigenvalues[mask], bigwig_values[mask])[0]
    corr_flip = scipy.stats.pearsonr(eigenvalues[mask], -bigwig_values[mask])[0]
    if corr_flip > corr_noflip:
        print(f"Flipped {chr} for corr from {corr_noflip} to {corr_flip}")
        eigenvalues = -eigenvalues
    results.append(eigenvalues)
    all_bins.extend(bins)

results = np.concatenate(results)
eigenvalues = np.nan_to_num(results)

with open(output_path, "w") as f:
    for (chr, start, end), value in zip(all_bins, eigenvalues):
        f.write(f"{chr}\t{start}\t{end}\t{value}\n")
