#!/usr/bin/env python3
#
# Merge hic_breakfinder and eaglec calls.  Essentially a union of the two
# callsets, with the exception that breakfinder calls are skipped if there is an
# eaglec call with both anchors within 1MB of the breakfinder call.
#
# Note: designed to work with the default output format from Arima-SV-Pipeline's
# hic_breakfinder output and EagleC combined output (*.CNN_SVs.5K_combined.txt), i.e.:
#   - Skips first row of hic_breakfinder call file (usually the header beginning with "#")
#   - Does NOT skip the first row of the EagleC output file (as there is no header by default)
#   - Expects hic_breakfinder chromosomes to be prefixed with "chr"
#   - Expects EagleC chromosomes NOT to be prefixed with "chr"
# Modify these assumptions as needed below.
#
# IMPORTANT NOTE:
#    - Outputs merged calls in hic_breakfinder's .breaks.bedpe format
#    - Gives each EagleC call a resolution of 5kb by default (this may not always be the best, but is done for simplicity)
#    - Gives every call a -logP of 0 (i.e. blank value to indicate these are merged calls; does NOT use the hic_breakfinder -logP value)
#    - Sets "start" and "end" values for EagleC calls to the same value (i.e. the "pos1" and "pos2" values). This may
#      cause issues with software that expects the start and end values to be different for proper .bedpe format.
#      Trivial to change (e.g. set end to start +/- resolution to set it to the width of a bin)
# Modify these as needed below.
#
# Dependencies:
#     - python 3.10 or higher
#     - pandas
#
# Usage:
#     python3 merge_breakfinder_eaglec.py breakfinder_filepath eaglec_filepath output_filepath

import sys
import pandas as pd

breakfinder_file = sys.argv[1]
eaglec_file = sys.argv[2]
output = sys.argv[3]

breakfinder_columns = ["chr1", "x1", "x2", "chr2", "y1", "y2", "strand1", "strand2", "resolution", "-logP"]
eaglec_columns = ["unp_chr1", "unp_chr2", "strands", "pos1", "pos2", "category"]

breakfinder = pd.read_csv(breakfinder_file, sep="\s+", skiprows=1, names=breakfinder_columns)
eaglec = pd.read_csv(eaglec_file, sep="\s+", names=eaglec_columns)

def chr_prefix(chr: str) -> str:
    if not isinstance(chr, str):
        chr = str(chr)
    if chr.startswith("chr"):
        return chr
    else:
        return f"chr{chr}"

threshold = 1000000

merged = []

# Go through EagleC calls first and include them all
for _, breakpoint in eaglec.iterrows():
    chr1 = chr_prefix(breakpoint["unp_chr1"])
    x1 = breakpoint["pos1"]
    x2 = breakpoint["pos1"]
    chr2 = chr_prefix(breakpoint["unp_chr2"])
    y1 = breakpoint["pos2"]
    y2 = breakpoint["pos2"]
    strand1, strand2 = breakpoint["strands"]
    # NOTE : uses a default resolution of 5000 for EagleC calls (and a dummy -logP of 0)
    resolution = "5kb"
    negLogP = 0
    # Note: exclude intra-chromosomal breakpoints less than 1MB apart and with -+ strandness (very tricky to distinguish from TADs)
    if chr1 == chr2 and abs( ((x1 + x2) / 2) - ((y1 + y2) / 2)) <= 1000000 and strand1 == "-" and strand2 == "+":
        continue
    merged.append((chr1, x1, x2, chr2, y1, y2, strand1, strand2, resolution, negLogP))

 # Go through all breakfinder calls next, and include them if they're not already in the EagleC calls (+/- threshold default 1Mb)
for _, breakpoint in breakfinder.iterrows():
    chrA = breakpoint["chr1"]
    strandA = breakpoint["strand1"]
    posA = breakpoint["x2"] if strandA == "+" else breakpoint["x1"]
    chrB = breakpoint["chr2"]
    strandB = breakpoint["strand2"]
    posB = breakpoint["y2"] if strandB == "+" else breakpoint["y1"]
    found_nearby_eaglec = False
    # Check if there is an eagle c breakpoint within threshold (default 1Mb)
    for _, eaglec_breakpoint in eaglec.iterrows():
        eaglec_chrA = chr_prefix(eaglec_breakpoint["unp_chr1"])
        eaglec_chrB = chr_prefix(eaglec_breakpoint["unp_chr2"])
        eaglec_posA = eaglec_breakpoint["pos1"]
        eaglec_posB = eaglec_breakpoint["pos2"]
        if chrA == eaglec_chrA and chrB == eaglec_chrB:
            if abs(posA - eaglec_posA) < threshold and abs(posB - eaglec_posB) < threshold:
                found_nearby_eaglec = True
                break
        elif chrA == eaglec_chrB and chrB == eaglec_chrA:
            if abs(posA - eaglec_posB) < threshold and abs(posB - eaglec_posA) < threshold:
                found_nearby_eaglec = True
                break
    if not found_nearby_eaglec:
        x1 = breakpoint["x1"]
        x2 = breakpoint["x2"]
        y1 = breakpoint["y1"]
        y2 = breakpoint["y2"]
        # NOTE: uses a dummy value of 0 for -logP REGARDLESS of hic_breakfinder -logP value
        merged.append((chrA, x1, x2, chrB, y1, y2, strandA, strandB, breakpoint["resolution"], 0))

with open(output, "w") as f:
    f.write("#chr1\tx1\tx2\tchr2\ty1\ty2\tstrand1\tstrand2\tresolution\t-logP\n")
    for (chr1, x1, x2, chr2, y1, y2, strand1, strand2, resolution, negLogP) in merged:
        f.write(f"{chr1}\t{x1}\t{x2}\t{chr2}\t{y1}\t{y2}\t{strand1}\t{strand2}\t{resolution}\t{negLogP}\n")
