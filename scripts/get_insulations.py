#!/usr/bin/env python
# Calculates insulation scores (args sample_id, resolution, window)

import sys
import cooler
import cooltools

sample_id = sys.argv[1]
resolution = int(sys.argv[2])
window = int(sys.argv[3])

def read_cooler(sample_id: str, resolution: int):
    return cooler.Cooler(f"../data/mcool/{sample_id}.mcool::/resolutions/{resolution}")

cool_file = read_cooler(sample_id, resolution)
insulation_table = cooltools.insulation(cool_file, [window], verbose=False, nproc=8)
insulation_table.to_csv(f"../data/insulations/res_{resolution}_window_{window}/{sample_id}_insulation_r{resolution}_w{window}.csv")