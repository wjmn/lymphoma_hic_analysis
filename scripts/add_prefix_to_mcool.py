#!/usr/bin/env python3
# Adds chr prefix to chromosomes in mcool file

import cooler, sys

in_mcool = sys.argv[1]

for resolution in [5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000, 2500000]:
    clr = cooler.Cooler(f"{in_mcool}::/resolutions/{resolution}")
    chromlist = clr.chromnames
    if chromlist[0].startswith("chr") or chromlist[1].startswith("chr"):
        break
    else:
        chrom_map = {c:'chr'+c for c in chromlist}
        cooler.rename_chroms(clr, chrom_map)
