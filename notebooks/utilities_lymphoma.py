"""Initial utilities for Hi-C notebook analyses - later deprecated"""

################################################################################
# IMPORTS
################################################################################

from initial_utilities.definitions import *
from initial_utilities.constants import *
from initial_utilities.utilities import *
from initial_utilities.plotters import *
from initial_utilities.readers import *
from initial_utilities.generator import *

import pandas as pd
import numpy as np

import os

from collections import Counter

#from scipy.spatial.distance import jaccard as jaccard_distance
#from sklearn.decomposition import PCA

import pyBigWig

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Arc
#from matplotlib_venn import venn2

import cooler
#from neoloop.visualize.core import *

#from pyryotype import GENOME, plot_ideogram

from scipy import stats
import scipy
import scipy.cluster.hierarchy as sch 

#import seaborn as sns
from dataclasses import field

#import cooltools

################################################################################
# TYPES
################################################################################

@dataclass
class Position: 
    """A single genomic position."""
    chr: str
    pos: int

class FeatureCall: 

    def __init__(self, chr: str, start: int, end: int, pos: int):
        self.chr = chr
        self.start = start
        self.end = end
        self.pos = pos

    def to_position(self) -> Position:
        return Position(self.chr, self.pos)

    def to_region(self) -> Region:
        return Region(self.chr, self.start, self.end)

    def make_window(self, window_size=40*25000) -> tuple[Region, Region]:
        position = self.to_position()
        region = Region(position.chr, position.pos - window_size//2, position.pos + window_size//2)
        return (region, region)
    
    def __str__(self):
        return f"{self.chr}:{self.start}-{self.end} ({self.chr}:{self.pos})"
    
    def __repr__(self):
        return self.__str__()

class DisjointFeatureCall:

    def __init__(self, chrA: str, startA: int, endA: int, posA: int, chrB: str, startB: int, endB: int, posB: int):
        self.chrA = chrA
        self.startA = startA
        self.endA = endA
        self.posA = posA
        self.chrB = chrB
        self.startB = startB
        self.endB = endB
        self.posB = posB
    
    def to_positions(self) -> tuple[Position, Position]:
        return (Position(self.chrA, self.posA), Position(self.chrB, self.posB))
    
    def to_regions(self) -> tuple[Region, Region]:
        return (Region(self.chrA, self.startA, self.endA), Region(self.chrB, self.startB, self.endB))
    
    def make_window(self, window_size=40*25000) -> tuple[Region, Region]:
        regionA = Region(self.chrA, self.posA - window_size//2, self.posA + window_size//2)
        regionB = Region(self.chrB, self.posB - window_size//2, self.posB + window_size//2)
        return (regionA, regionB)

    def __str__(self):
        return f"{self.chrA}:{self.startA}-{self.endA} ({self.chrA}:{self.posA}) - {self.chrB}:{self.startB}-{self.endB} ({self.chrB}:{self.posB})"

    def __repr__(self):
        return self.__str__()
    



class Loop(DisjointFeatureCall):
    def __init__(self, chrA: str, startA: int, endA: int, posA: int, chrB: str, startB: int, endB: int, posB: int):
        super().__init__(chrA, startA, endA, posA, chrB, startB, endB, posB)


def read_loops(loops_file: str) -> list[Loop]:
    #loops_file = f"../data/loops/{id}/{id}_merged_loops.bedpe"
    if not os.path.isfile(loops_file):
        return None
    try:
        df = pd.read_csv(loops_file, sep="\t", header=None, skiprows=0, names=["chr1", "start1", "end1", "chr2", "start2", "end2", "pvalue"],
                        dtype={"chr1": str, "start1": int, "end1": int, "chr2": str, "start2": int, "end2": int, "pvalue": float})
        # Filter by p value < 0.05
        df = df[df.pvalue < 0.05]
    except pd.errors.EmptyDataError:
        return []
    loops = []
    for _, row in df.iterrows():
        chrA, startA, endA, chrB, startB, endB = row[:6]
        posA = (startA + endA) // 2
        posB = (startB + endB) // 2
        loops.append(Loop(chr_prefix(chrA), startA, endA, posA, chr_prefix(chrB), startB, endB, posB))
    return loops


class EagleCall(DisjointFeatureCall):

    def __init__(self, breakpointA: Breakpoint, breakpointB: Breakpoint, pairing: Pairing, category: VariantCategory):
        super().__init__(
            breakpointA.chr, breakpointA.start, breakpointA.end, breakpointA.pos,
            breakpointB.chr, breakpointB.start, breakpointB.end, breakpointB.pos,
        )
        self.strandA = breakpointA.strand
        self.strandB = breakpointB.strand
        self.pairing = pairing
        self.category = category


class CuratedCall(DisjointFeatureCall):

    def __init__(self, breakpointA: Breakpoint, breakpointB: Breakpoint, pairing: Pairing, category: VariantCategory):
        super().__init__(
            breakpointA.chr, breakpointA.start, breakpointA.end, breakpointA.pos,
            breakpointB.chr, breakpointB.start, breakpointB.end, breakpointB.pos,
        )
        self.strandA = breakpointA.strand
        self.strandB = breakpointB.strand
        self.pairing = pairing
        self.category = category


def breakfinder_to_curated(sv: BreakfinderCall):
    return CuratedCall(sv.breakpointA, sv.breakpointB, sv.pairing, sv.category)

def read_curated_calls(curated_file: str, threshold=5000000) -> list[CuratedCall] | None:
    #curated_file = f"../data/curated_breakpoints/{id}_curated_breakpoints.tsv"
    if not os.path.isfile(curated_file):
        return None
    try:
        df = pd.read_csv(curated_file, sep="\t")
    except pd.errors.EmptyDataError:
        return []
    calls = []
    if curated_file.endswith(".breaks.bedpe"):
        return list(map(breakfinder_to_curated, read_breakfinder_data(curated_file)))
    for _, (chrA, chrB, strands, posA, posB, category) in df.iterrows():
        strandA, strandB = strands
        breakpointA = Breakpoint(chr_prefix(chrA), posA, posA, posA, to_strand(strandA))
        breakpointB = Breakpoint(chr_prefix(chrB), posB, posB, posB, to_strand(strandB))
        if category == "duplication":
            category = VariantCategory.DUPLICATION
        elif category == "deletion":
            category = VariantCategory.DELETION
        elif category == "inversion":
            category = VariantCategory.INVERSION_OR_OTHER
        elif category == "translocation":
            category = VariantCategory.TRANSLOCATION
        pairing = Pairing.INTRA if chrA == chrB else Pairing.INTER

        # NOTE: IMPORTANT
        #Preprocess by removing all SVs with anchors within threshold of each other
        if chrA == chrB and abs(posA - posB) <= threshold:
            continue

        calls.append(CuratedCall(breakpointA, breakpointB, pairing, category))
    return calls

def read_cooler(sample_id: str, resolution: int):
    return cooler.Cooler(f"../data/mcool/{sample_id}.mcool::/resolutions/{resolution}")

def read_eagle_calls(eaglec_file: str) -> list[EagleCall] | None:
    #eaglec_file = f"../data/eaglec/{id}.CNN_SVs.5K_combined.txt"
    if not os.path.isfile(eaglec_file):
        return None
    try:
        df = pd.read_csv(eaglec_file, sep="\t")
    except pd.errors.EmptyDataError:
        return []
    calls = []
    for _, (chrA, chrB, strands, posA, posB, category) in df.iterrows():
        strandA, strandB = strands
        breakpointA = Breakpoint(chr_prefix(chrA), chr_prefix(posA), posA, posA, to_strand(strandA))
        breakpointB = Breakpoint(chr_prefix(chrB), chr_prefix(posB), posB, posB, to_strand(strandB))
        if category == "duplication":
            category = VariantCategory.DUPLICATION
        elif category == "deletion":
            category = VariantCategory.DELETION
        elif category == "inversion":
            category = VariantCategory.INVERSION_OR_OTHER
        elif category == "translocation":
            category = VariantCategory.TRANSLOCATION
        pairing = Pairing.INTRA if chrA == chrB else Pairing.INTER
        calls.append(EagleCall(breakpointA, breakpointB, pairing, category))


    return calls


class Neoloop(DisjointFeatureCall):
    
        def __init__(self, chrA: str, startA: int, endA: int, posA: int, chrB: str, startB: int, endB: int, posB: int, resolution: int, assembly: str):
            super().__init__(chrA, startA, endA, posA, chrB, startB, endB, posB)
            self.resolution = resolution
            self.assembly = assembly

def read_neoloops(neoloop_file: str, prob=90) -> list[Neoloop] | None:
    #neoloop_file = f"../data/neoloops/{prob}/{id}.neoloops.bedpe"
    if not os.path.isfile(neoloop_file):
        return None
    try:
        df = pd.read_csv(neoloop_file, sep='\s+')
    except pd.errors.EmptyDataError:
         return []
    neoloops = []
    for _, (chrA, startA, endA, chrB, startB, endB, assembly) in df.iterrows():
        posA = (startA + endA) // 2
        posB = (startB + endB) // 2
        resolution = endA - startA
        assert resolution == endB - startB
        neoloops.append(Neoloop(chrA, startA, endA, posA, chrB, startB, endB, posB, resolution, assembly))
    return neoloops



class Neotad(DisjointFeatureCall):
    
    def __init__(self, chrA: str, startA: int, endA: int, posA: int, chrB: str, startB: int, endB: int, posB: int, assembly: str):
        super().__init__(chrA, startA, endA, posA, chrB, startB, endB, posB)
        self.assembly = assembly

def read_neotads(neotad_file: str, res=25000) -> list[Neotad] | None:
    #neotad_file = f"../data/neotads/{res}/{id}_{res}.neotad.txt"
    if not os.path.isfile(neotad_file):
        return None
    try:
        df = pd.read_csv(neotad_file, sep='\s+')
    except pd.errors.EmptyDataError:
        return []
    neotads = []
    for _, (chrA, startA, endA, chrB, startB, endB, assembly) in df.iterrows():
        posA = (startA + endA) // 2
        posB = (startB + endB) // 2
        neotads.append(Neotad(chrA, startA, endA, posA, chrB, startB, endB, posB, assembly))
    return neotads

def read_compartments(compartment_file: str, resolution=100000) -> NDArray | None:
    # if not os.path.isfile(compartment_file):
    #compartment_file = f"../data/compartments/{resolution}/{id}.{resolution}.bedgraph"
    try:
        df = pd.read_csv(compartment_file, sep='\s+', header=None)#, names=["chr", "start", "end", "score"])
    except pd.errors.EmptyDataError:
        return []
    except FileNotFoundError:
        return None
    compartments = np.array(df.iloc[:, 3])#["score"])
    return compartments


example_compartments_100kb = pd.read_csv("../data/compartments/100000/12_MGH.100000.bedgraph", sep="\t", header=None, names=["chr", "start", "end", "score", "bin"])
COMPARTMENT_REGIONS_100kb = example_compartments_100kb[["chr", "start", "end"]]

def read_hicrep_result(file_1, file_2):
    name_1 = f"../data/hicrep/100000/{file_1}_vs_{file_2}_scc.txt"
    name_2 = f"../data/hicrep/100000/{file_2}_vs_{file_1}_scc.txt"
    if os.path.exists(name_1):
        with open(name_1, "r") as f:
            lines = f.readlines()
    elif os.path.exists(name_2):
        with open(name_2, "r") as f:
            lines = f.readlines()
    else:
        raise FileNotFoundError(f"Could not find {name_1} or {name_2}")
    # Skip first two lines
    lines = lines[2:]
    assert len(lines) == 22
    return [float(line.strip()) for line in lines]

class Tad(FeatureCall):

    def __init__(self, chr: str, start: int, end: int, pos: int, hierarchy: int):
        super().__init__(chr, start, end, pos)
        self.hierarchy = hierarchy

def read_tads(tad_file: str, res=25000) -> list[Tad] | None: 
    #tad_file = f"../data/tads/{res}/{id}.txt"
    if not os.path.isfile(tad_file):
        return None
    try:
        df = pd.read_csv(tad_file, sep='\s+')
    except pd.errors.EmptyDataError:
        return []
    tads = []
    for _, (chr, start, end, hierarchy) in df.iterrows():
        pos = (start + end) // 2
        tads.append(Tad(chr_prefix(chr), start, end, pos, hierarchy))
    return tads

def tad_to_region(tad: Tad) -> Region:
    return Region(tad.chr, tad.start, tad.end)


def read_assemblies(assembly_file: str) -> list[str] | None:
    #assembly_file = f"../data/assemblies/{id}.assemblies.txt"
    if not os.path.isfile(assembly_file):
        return None
    with open(assembly_file, "r") as f:
        assemblies = f.read().splitlines()
    return assemblies


# Next, sort by unique_valid_pairs from qc
def get_uvp(sample_id):
    return pd.read_csv(f"../data/qc_deep/{sample_id}_v1.3_Arima_QC_deep.txt", sep='\s+')["Unique_valid_pairs"].values[0]

def get_hic_path(sample_id: str): return f"../data/hic/{sample_id}_inter_30.hic"
def get_qc_path(sample_id: str): return  f"../data/qc_deep/{sample_id}_v1.3_Arima_QC_deep.txt"
def get_breakfinder_path(sample_id: str): return f"../data/hic_breakfinder/{sample_id}.breaks.bedpe"


@dataclass
class DataProfile:
    hic: str | None
    qc: str | None
    breakfinder: str | None
    mcool: str | None
    cnv_profile: str | None
    cnv_segment: str | None
    eaglec: str | None
    curated: str | None 
    compartments: str | None 
    tads: str | None 
    loops_hiccups: str | None
    loops_hicexplorer: str | None
    neoloops: str | None

def exists_or_none(path: str):
    return path if os.path.exists(path) else None

def make_profile(sample_id: str, comp_res=100000):
    return DataProfile(
        hic=exists_or_none(f"../data/hic/{sample_id}_inter_30.hic"),
        qc=exists_or_none(f"../data/qc_deep/{sample_id}_v1.3_Arima_QC_deep.txt"),
        breakfinder=exists_or_none(f"../data/hic_breakfinder/{sample_id}.breaks.bedpe"),
        mcool=exists_or_none(f"../data/mcool/{sample_id}.mcool"),
        cnv_profile=exists_or_none(f"../data/cnv_profile/25000/{sample_id}_CNV_profile_25000.bedgraph"),
        cnv_segment=exists_or_none(f"../data/cnv_segment/25000/{sample_id}_CNV_segment_25000.bedgraph"),
        eaglec=exists_or_none(f"../data/eaglec/{sample_id}.CNN_SVs.5K_combined.txt"),
        curated=exists_or_none(f"../data/curated_breakpoints/{sample_id}_curated_breakpoints.tsv") or exists_or_none(f"../data/curated_bedpe/{sample_id}.breaks.bedpe"),
        compartments=exists_or_none(f"../data/compartments/{comp_res}/{sample_id}.{comp_res}.bedgraph"),
        tads=exists_or_none(f"../data/tads/25000/{sample_id}.bed"),
        loops_hiccups=exists_or_none(f"../data/loops/{sample_id}/{sample_id}_merged_loops.bedpe") or exists_or_none(f"../data/loops/{sample_id}/merged_loops.bedpe"),
        loops_hicexplorer=exists_or_none(f"../data/hicexplorer_loops/10000/{sample_id}_hicexplorer_loops_10000.bedpe") or exists_or_none(f"../data/loops/{sample_id}/merged_loops.bedpe"),
        neoloops=exists_or_none(f"../data/neoloops/90/{sample_id}.neoloops.bedpe"),
    )

################################################################################
# OVERRIDES

# THIS FUNCTION HAS BEEN MODIFIED from hicdash: arrowheads not shown at end, but instead shown above
def plot_gene_track(
    chr: str,
    start: int,
    end: int,
    gene_filter: list[str] | None = None,
    hide_axes=True,
    vertical=False,
    ax=None,
    fontsize=8,
    max_rows=6,
    min_rows=3,
    protein_coding_only=True,
    crosshairs: bool=False, 
    plotted_crosshairs: list[tuple[str, int, str, float]]=[],
    centered_names=False,
    arrowhead_length=None,
    arrowhead_length_proportion=None,
    arrowhead_width=None,
    arrow_length=None,
    all_same_line=False,
    show_arrows=True,
) -> plt.Axes:
    """Plot a gene track (based on GENE_ANNOTATIONS) for a given range.

    Note: if too many genes are in the given range, then only a subset of genes will be plotted.
    (otherwise the track will be too crowded to read).
    """

    center = (start + end) // 2

    # Requires unprefixed chromosome
    # Get genes in gene regions (only protein-coding genes or IG genes)
    # If there's a gene filter, check through all genes
    # If not, then choose only protein-coding genes (if specified as True)
    candidates = GENE_ANNOTATIONS.genes_at_locus(
        contig=chr_unprefix(chr), position=start, end=end
    )
    if gene_filter:
        gene_filter = set([gene.upper() for gene in gene_filter])
        genes = list(filter(lambda g: g.gene_name in gene_filter, candidates))
    else:
        if protein_coding_only:
            genes = list(
                filter(
                    lambda g: g.biotype == "protein_coding" and g.gene_name != "",
                    candidates,
                )
            )
        else:
            genes = list(filter(lambda g: g.gene_name != "", candidates))

        # Narrow down the genes if there are too many to plot, keeping only a number of genes around the center
        # Keep genes whose midpoints are closest to the center
        genes.sort(key=lambda g: abs(center - (g.start + g.end) / 2))
        genes = genes[:max_rows]

    # Now sort genes by their start position
    genes.sort(key=lambda gene: gene.start)

    # Plot genes
    if ax is None:
        ax = plt.gca()

    # Prepare axes
    if vertical:
        ax.set_ylim(start, end)
        ax.invert_yaxis()
        ax.invert_xaxis()
    else:
        ax.set_xlim(start, end)

    if hide_axes:
        ax.yaxis.set_visible(False)
        ax.xaxis.set_visible(False)
        ax.spines[["top", "right", "left", "bottom"]].set_visible(False)

    # Calculate values for plotting

    # Keep track of how many genes are plotted
    plot_counter = 0

    # Use 5 percent as a temporary padding value
    pct_5 = 5 * (end - start) / 100
    if arrowhead_length is None:
        if arrowhead_length_proportion is None:
            arrowhead_length = (end - start) / 75
        else:
            arrowhead_length = (end - start) * arrowhead_length_proportion
    plot_line_width = 0.4
    if arrowhead_width is None:
        arrowhead_width = 0.25
    if arrow_length is None:
        arrow_length = (end - start) / 50

    # Get genes which intersect the center directly
    direct_genes = GENE_ANNOTATIONS.genes_at_locus(
        contig=chr_unprefix(chr), position=center
    )
    direct_genes_set = set([gene.gene_name for gene in direct_genes])

    for gene in genes:

        exon_width_start = plot_counter - (plot_line_width / 2)
        exon_width = plot_line_width
        # Plot base arrow indicating strandness
        arr_start, full_length, arr_length = (
            (gene.start, gene.end - gene.start, arrow_length)
            if gene.strand == "+"
            else (gene.end, gene.start - gene.end, - arrow_length)
        )
        safe_start = max(start, gene.start) if gene.strand == "+" else (min(end, gene.end))
        ec, fc, ls = ("gray", "white", "solid") if (safe_start != arr_start) else ("black", "black", "solid")
        arr_axis_line = plot_counter+plot_line_width/2 + arrowhead_width/2
        if vertical:
            ax.vlines(plot_counter, gene.start, gene.end, colors="black")
        else:
            ax.hlines(plot_counter, gene.start, gene.end, colors="black")
        if show_arrows:
            if vertical:
                ax.hlines(arr_start, plot_counter, plot_counter+plot_line_width/2, colors="black", ls=":", lw=0.5)
                ax.hlines(arr_start, plot_counter+plot_line_width/2, arr_axis_line, colors="black", lw=0.5)
                ax.arrow(
                    arr_axis_line,
                    safe_start,
                    0,
                    arr_length,
                    head_width=arrowhead_width,
                    head_length=arrowhead_length,
                    length_includes_head=False,
                    lw=0.5,
                    ls=ls,
                    ec=ec,
                    fc=fc,
                )
            else:
                ax.vlines(arr_start, plot_counter, plot_counter+plot_line_width/2, colors="black", ls=":", lw=0.5)
                ax.vlines(arr_start, plot_counter+plot_line_width/2, arr_axis_line, colors="black", lw=0.5)
                ax.arrow(
                    safe_start,
                    arr_axis_line,
                    arr_length,
                    0,
                    head_width=arrowhead_width,
                    head_length=arrowhead_length,
                    length_includes_head=False,
                    lw=0.5,
                    ls=ls,
                    ec=ec,
                    fc=fc,
                )

        # Plot each exon as a rectangle on the gene line
        for exon in gene.exons:
            exon_length = exon.end - exon.start
            if vertical:
                ax.add_patch(
                    Rectangle(
                        (exon_width_start, exon.start),
                        exon_width,
                        exon_length,
                        edgecolor="black",
                        facecolor="black",
                    )
                )
            else:
                ax.add_patch(
                    Rectangle(
                        (exon.start, exon_width_start),
                        exon_length,
                        exon_width,
                        edgecolor="black",
                        facecolor="black",
                    )
                )

        # Get text position - and avoid putting the text out of axes
        # Bolden the text if passes directly through the center of the range (assuming a breakpoint is in the center)
        # fontweight = "bold" if gene.gene_name in direct_genes_set else "normal"
        fontweight = "normal"
        color = "black"
        if vertical:
            if centered_names:
                text_position = (gene.end + gene.start) / 2
                ha = "right"
                va = "center"
                if text_position + pct_5 > ax.get_ylim()[0]:
                    text_position = ax.get_ylim()[0]
                    va = "bottom"
                elif  text_position - pct_5 < ax.get_ylim()[1]:
                    text_position = ax.get_ylim()[1]
                    va = "top"
                ax.text(
                    plot_counter+plot_line_width/2,
                    text_position,
                    gene.gene_name,
                    ha=ha,
                    va=va,
                    fontsize=fontsize,
                    rotation=90,
                    color=color,
                    fontweight=fontweight,
                )
            else:
                row_position = plot_counter
                va = "bottom"
                text_position = gene.start - (pct_5 * 0.35)
                # if gene.strand == "-":
                    # text_position -= arrowhead_length
                if text_position - pct_5 < ax.get_ylim()[1]:
                    text_position = gene.end + (pct_5 * 0.25)
                    va = "top"
                    # if gene.strand == "+":
                        # text_position += arrowhead_length
                    if text_position + pct_5 > ax.get_ylim()[0]:
                        text_position = ax.get_ylim()[0]
                        va = "bottom"
                        row_position = plot_counter + plot_line_width / 2
                ax.text(
                    row_position,
                    text_position,
                    gene.gene_name,
                    ha="center",
                    va=va,
                    fontsize=fontsize,
                    rotation=90,
                    color=color,
                    fontweight=fontweight,
                ).set_clip_on(True)
        else:
            if centered_names:
                text_position = (gene.end + gene.start) / 2
                ha = "center"
                va = "bottom"
                if text_position - pct_5 < ax.get_xlim()[0]:
                    text_position = ax.get_xlim()[0]
                    ha = "left"
                elif text_position + pct_5 > ax.get_xlim()[1]:
                    text_position = ax.get_xlim()[1]
                    ha = "right"
                ax.text(
                    text_position,
                    plot_counter+plot_line_width/2,
                    gene.gene_name,
                    ha=ha,
                    va=va,
                    fontsize=fontsize,
                    color=color,
                    fontweight=fontweight,
                )
            else:
                ha = "right"
                row_position = plot_counter
                text_position = gene.start - (pct_5 * 0.25)
                # if gene.strand == "-":
                    # text_position -= arrowhead_length
                if text_position - pct_5 < ax.get_xlim()[0]:
                    text_position = gene.end + (pct_5 * 0.25)
                    ha = "left"
                    # if gene.strand == "+":
                        # text_position += arrowhead_length
                    if text_position + pct_5 > ax.get_xlim()[1]:
                        text_position = ax.get_xlim()[1]
                        ha = "right"
                        row_position = plot_counter + plot_line_width / 2
                ax.text(
                    text_position,
                    row_position,
                    gene.gene_name,
                    ha=ha,
                    va="center",
                    fontsize=fontsize,
                    color=color,
                    fontweight=fontweight,
                ).set_clip_on(True)

        # Increment plot counter
        if not all_same_line:
            plot_counter += 1

    # Extend the axis limits if centered_gene_names is used
    if centered_names and plot_counter >= min_rows:
        if vertical:
            xmax, xmin = ax.get_xlim()
            ax.set_xlim(xmax + plot_line_width / 2, xmin)
        else:
            ymin, ymax = ax.get_ylim()
            ax.set_ylim(ymin, ymax + plot_line_width / 2)

    # If only a few genes were plotted, then add a bit of space padding to the plot
    if plot_counter < min_rows and not all_same_line:
        if vertical:
            ax.set_xlim(min_rows + plot_line_width + arrowhead_width/2, 0 - plot_line_width)
        else:
            ax.set_ylim(0 - plot_line_width, min_rows + plot_line_width + arrowhead_width/2)

    # Add crosshairs if specified
    if crosshairs:
        for chr, pos, col, alpha in plotted_crosshairs:
            if chr == chr:
                if vertical:
                    ax.axhline(
                        pos,
                        color=col,
                        linestyle=(0, (1, 5)),
                        linewidth=1,
                        alpha=alpha,
                    )
                else:
                    ax.axvline(
                        pos,
                        color=col,
                        linestyle=(0, (1, 5)),
                        linewidth=1,
                        alpha=alpha,
                    )

    return ax

def plot_bigwig_track(
    bw_handle: object,
    chr: str,
    start: int,
    end: int,
    num_bins=1000,
    hide_axes=True,
    vertical=False,
    ax=None,
    fontsize=8,
    label: str = "",
    color="blue",
    crosshairs: bool=False, 
    ymax=None,
    plotted_crosshairs: list[tuple[str, int, str, float]]=[],
) -> plt.Axes:

    # Check if chromosomes are prefixed or unprefixed
    # MODIFIED HERE -----------------------------------------------------------------------------------------------------
    if "chr8" in bw_handle.chroms().keys():
    # END MODIFIED HERE -----------------------------------------------------------------------------------------------------
        chr_prefixed = chr
    else:
        chr_prefixed = chr_prefix(chr)

    # Get the data from the bigwig file
    # FIrst ensure bounds are safe 
    safe_start = max(0, start)
    safe_end = min(CHROM_SIZES[chr_prefixed], end)
    safe_nbins = (safe_end - safe_start) // ((end - start) // num_bins)
    data = bw_handle.stats(chr, safe_start, safe_end, type="mean", nBins=safe_nbins, numpy=True)
    # Pad with extra zeros if was out of bounds initially
    bin_width = (end - start) / num_bins
    if start < 0:
        data = np.pad(data, (int(-start // bin_width), 0), "constant")
    if end > CHROM_SIZES[chr_prefixed]:
        data = np.pad(data, (0, int((end - CHROM_SIZES[chr]) // bin_width)+1), "constant")

    if ymax is None:
        ymax = np.sqrt(bw_handle.header()["maxVal"])

    positions = np.linspace(start, end, num_bins)

    # Plot the data depending on if horizontal or vertical axis
    if ax is None:
        ax = plt.gca()

    if vertical:
        # Fill from right to data peak on left
        ax.fill_betweenx(positions, np.zeros(num_bins), data, color=color)
        ax.set_ylim(start, end)
        ax.invert_yaxis()
        ax.set_xlim(0, ymax)
        ax.invert_xaxis()
    else:
        ax.fill_between(positions, np.zeros(num_bins), data, color=color)
        ax.set_xlim(start, end)
        ax.set_ylim(0, ymax)

    if hide_axes:
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.spines[["top", "right", "left", "bottom"]].set_visible(False)

    # Add crosshairs if specified
    if crosshairs:
        for chr, pos, col, alpha in plotted_crosshairs:
            if chr == chr:
                if vertical:
                    ax.axhline(
                        pos,
                        color=col,
                        linestyle=(0, (1, 5)),
                        linewidth=1,
                        alpha=alpha,
                    )
                else:
                    ax.axvline(
                        pos,
                        color=col,
                        linestyle=(0, (1, 5)),
                        linewidth=1,
                        alpha=alpha,
                    )


    # Add label
    if vertical:
        ax.text(
            0,
            1,
            label,
            ha="left",
            va="top",
            transform=ax.transAxes,
            color="gray",
            fontdict={"fontsize": fontsize},
            rotation=90,
        )
    else:
        ax.text(
            0,
            1,
            label,
            ha="left",
            va="top",
            transform=ax.transAxes,
            color="gray",
            fontdict={"fontsize": fontsize},
        )


    
################################################################################
# HELPERS
################################################################################
def overlap(regionA: Region, regionB: Region) -> bool:
    """Determine if two regions overlap."""
    if regionA.chr != regionB.chr:
        return False
    if regionA.start <= regionB.start <= regionA.end or regionB.start <= regionA.start <= regionB.end \
        or regionA.start <= regionB.end <= regionA.end or regionB.start <= regionA.end <= regionB.end:
        return True   

def get_mcool_region_data(mcool_filepath, regionX: Region, regionY: Region, resolution: int, balance=None):
    """Balance options: None, sweight, weight"""
    clr = cooler.Cooler(f"{mcool_filepath}::/resolutions/{resolution}")
    fetchX = f"{regionX.chr}:{regionX.start}-{regionX.end}"
    fetchY = f"{regionY.chr}:{regionY.start}-{regionY.end}"
    data = clr.matrix(balance=balance).fetch(fetchY, fetchX)
    return data
################################################################################
# CONVERSIONS
################################################################################

def curated_calls_to_breakfinder_calls(calls: list[CuratedCall]) -> list[BreakfinderCall]:
    breakfinder_calls = []
    for call in calls:
        breakfinder_calls.append(BreakfinderCall(
            breakpointA=Breakpoint(call.chrA, call.posA, call.posA, call.posA, call.strandA),
            breakpointB=Breakpoint(call.chrB, call.posB, call.posB, call.posB, call.strandB),
            resolution=5000,
            neg_log_pval=0,
            pairing=call.pairing,
            category=call.category
        ))
    return breakfinder_calls
    
################################################################################
# NEARBY GENES
################################################################################

def show_nearby_genes(chr, pos, width=5000, protein_coding_only=False):
    genes = [ gene for gene in GENE_ANNOTATIONS.genes_at_locus(contig=chr_unprefix(chr), position=pos-width, end=pos+width) if gene.gene_name != "" ]

    if protein_coding_only:
        genes = [ gene for gene in genes if gene.biotype == "protein_coding" ]

    genes.sort(key=lambda x: min(abs(x.start - pos), abs(x.end - pos)))
    return genes



def show_nearby_genes_stranded(
    call: CuratedCall | BedpeLine, width=300000, buffer=25000, protein_coding=True
):
    """Gets a list of genes around a breakpoint (based on strandness).

    Looks in direction of strandness for WIDTH bp and in opposite direction of strandness for BUFFER bp.
    """

    # Unpack breakfinder call
    chrA, posA, strandA = call.chrA, call.posA, call.strandA
    chrB, posB, strandB = call.chrB, call.posB, call.strandB

    if strandA == Strand.POS:
        startA = posA - width
        endA = posA + buffer
    else:
        startA = posA - buffer
        endA = posA + width

    if strandB == Strand.POS:
        startB = posB - width
        endB = posB + buffer
    else:
        startB = posB - buffer
        endB = posB + width

    # Get start and end ranges based on strandness
    non_empty = lambda x: x.gene_name != ""
    genesA = list(filter(
        non_empty,
        GENE_ANNOTATIONS.genes_at_locus(
            contig=chr_unprefix(chrA), position=startA, end=endA
        ),
    ))
    genesB = list(filter(
        non_empty,
        GENE_ANNOTATIONS.genes_at_locus(
            contig=chr_unprefix(chrB), position=startB, end=endB
        ),
    ))

    if protein_coding:
        genesA = [gene for gene in genesA if gene.biotype == "protein_coding"]
        genesB = [gene for gene in genesB if gene.biotype == "protein_coding"]

    # Sort genes by proximity to breakpoint
    genesA.sort(key=lambda x: min(abs(x.end - posA), abs(x.start - posA)))
    genesB.sort(key=lambda x: min(abs(x.end - posB), abs(x.start - posB)))

    return ([g.gene_name for g in genesA], [g.gene_name for g in genesB])

def show_nearby_genes_stranded_unrolled(chrA, posA, strandA, chrB, posB, strandB, width=300000, buffer=25000, protein_coding=True):
    """Gets a list of genes around a breakpoint (based on strandness).

    Looks in direction of strandness for WIDTH bp and in opposite direction of strandness for BUFFER bp.
    """

    if strandA == Strand.POS:
        startA = posA - width
        endA = posA + buffer
    else:
        startA = posA - buffer
        endA = posA + width

    if strandB == Strand.POS:
        startB = posB - width
        endB = posB + buffer
    else:
        startB = posB - buffer
        endB = posB + width

    # Get start and end ranges based on strandness
    non_empty = lambda x: x.gene_name != ""
    genesA = list(filter(
        non_empty,
        GENE_ANNOTATIONS.genes_at_locus(
            contig=chr_unprefix(chrA), position=startA, end=endA
        ),
    ))
    genesB = list(filter(
        non_empty,
        GENE_ANNOTATIONS.genes_at_locus(
            contig=chr_unprefix(chrB), position=startB, end=endB
        ),
    ))

    if protein_coding:
        genesA = [gene for gene in genesA if gene.biotype == "protein_coding"]
        genesB = [gene for gene in genesB if gene.biotype == "protein_coding"]

    # Sort genes by proximity to breakpoint
    genesA.sort(key=lambda x: min(abs(x.end - posA), abs(x.start - posA)))
    genesB.sort(key=lambda x: min(abs(x.end - posB), abs(x.start - posB)))

    return ([g.gene_name for g in genesA], [g.gene_name for g in genesB])

################################################################################
# PLOTTING TRIANGLES
################################################################################


@dataclass
class Segment:
    chr: str 
    start: int 
    end: int 
    is_minus_strand: bool = False

@dataclass
class BinAlignedSegment:
    segment: Segment
    resolution: int

@dataclass 
class BinAlignedSegmentPair: 
    segmentX: Segment
    segmentY: Segment
    resolution: int

@dataclass
class Assembly:
    segments: list[BinAlignedSegment]

class PlotRegionType(Enum):
    DATA_SEGMENT = 0 
    GAP = 1

@dataclass
class PlotRegion:
    bin_range: tuple[int, int]
    plot_region_type: PlotRegionType
    matched_segment: BinAlignedSegment | None = None

@dataclass
class AssembledHic: 
    data: NDArray
    plot_regions: list[PlotRegion]

class TrackType(Enum):
    TRACK_GENES = 0 
    TRACK_BIGWIG = 1
    TRACK_BED = 2
    TRACK_BED_DIRECTIONAL = 3
    TRACK_CHR_ARROWS = 4
    TRACK_VIRTUAL_4C = 5
    TRACK_BLANK = 6
    
@dataclass 
class TrackSpec:
    track_type: TrackType
    track_kwargs: dict = field(default_factory=dict)
    ylabel: str = ""
    rel_height: float = 1.0

def blank_axis(ax: plt.Axes):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[['top', 'right', 'left', 'bottom']].set_visible(False)


CHROM_COLORS = {'chr1': '#B15928',
            'chr2': '#6A3D9A',
            'chr3': '#CAB2D6',
            'chr4': '#FDBF6F',
            'chr5': '#A6CEE3',
            'chr6': '#B2DF8A',
            'chr7': '#1B9E77',
            'chr8': '#A6CEE3',
            'chr9': '#33A02C',
            'chr10':'#E6AB02',
            'chr11':'#E58606',
            'chr12':'#5D69B1',
            'chr13':'#FF7F00',
            'chr14':'#52BCA3',
            'chr15':'#99C945',
            'chr16':'#CC61B0',
            'chr17':'#88CCEE',
            'chr18':'#1D6996',
            'chr19':'#117733',
            'chr20':'#661100',
            'chr21':'#882255',
            'chr22':'#1F78B4',
            'chrX':'#666666',
            'chrY':'#5F4690'
            }

def bin_align_value(value: int, resolution: int) -> int:
    return value // resolution * resolution

def pair_segments(segmentX: BinAlignedSegment, segmentY: BinAlignedSegment) -> BinAlignedSegmentPair:
    assert segmentX.resolution == segmentY.resolution
    return BinAlignedSegmentPair(
        segmentX = segmentX.segment,
        segmentY = segmentY.segment,
        resolution = segmentX.resolution
    )


def bin_align_segment(segment: Segment, resolution: int) -> BinAlignedSegment:
    return BinAlignedSegment(
        segment = Segment(
            chr = segment.chr,
            start = bin_align_value(segment.start, resolution),
            end = bin_align_value(segment.end, resolution),
            is_minus_strand = segment.is_minus_strand
        ),
        resolution = resolution
    )

def switch_segments(segment_pair: BinAlignedSegmentPair) -> BinAlignedSegmentPair:
    return BinAlignedSegmentPair(
        segmentX = segment_pair.segmentY,
        segmentY = segment_pair.segmentX,
        resolution = segment_pair.resolution
    )
    
def get_mcool_for_segment_pair(
        mcool_path,
        segment_pair: BinAlignedSegmentPair, # coords must be bin-aligned to resolution
        balance = None,
        measure: str = "observed"
    ) -> NDArray:
    """
    Args:
        - hic: HiC handle
        - segment_pair: pair of bin-aligned segments with same resolution
        - norm (default "NONE"): also "VC" or "VC_SQRT" or "KR"
        - measure (default "observed"): also "oe" or "expected"

    Gives raw Hi-C matrix output (not normalized)

    Fails if any coordinates are outside range (0, chrom_size).
    """
    # Get matrix zoom data
    x = segment_pair.segmentX
    y = segment_pair.segmentY
    res = segment_pair.resolution
    regionX = Region(x.chr, x.start, x.end)
    regionY = Region(y.chr, y.start, y.end)
    data = np.nan_to_num(get_mcool_region_data(mcool_path, regionX, regionY, res, balance))
    return data

def assemble_mcool(
    mcool_path,
    assembly: Assembly,
    balance = None,
    gap_size: int = 0,
) -> NDArray:
    """
    Args:
        - hic: HiC handle
        - assembly: list of bin-aligned segments with same resolution
        - norm (default "NONE"): also "VC" or "VC_SQRT" or "KR"
        - measure (default "observed"): also "oe" or "expected"
        - gap_size (default 0): number of bins of nan to insert between segments
    """

    # Double check all resolutions match 
    resolution = None 
    for segment in assembly.segments:
        if resolution is None:
            resolution = segment.resolution
        assert segment.resolution == resolution

    # Get matrix data for each pair of segments; assemble row by row 
    rows = [] 
    for y, segmentY in enumerate(assembly.segments):
        row = []
        row_divisions = []
        col_pointer = 0 
        for x, segmentX in enumerate(assembly.segments):
            pair_data = get_mcool_for_segment_pair(mcool_path, pair_segments(segmentX, segmentY), balance=balance)
            if segmentX.segment.is_minus_strand:
                pair_data = np.flip(pair_data, axis=1)
            if segmentY.segment.is_minus_strand:
                pair_data = np.flip(pair_data, axis=0)
            row.append(pair_data)
            row_divisions.append((segmentX, (col_pointer, col_pointer + pair_data.shape[1])))
            col_pointer += pair_data.shape[1]
            if x < len(assembly.segments)-1 and gap_size > 0:
                gap_matrix = np.empty((pair_data.shape[0], gap_size))
                gap_matrix[:] = np.nan
                row.append(gap_matrix)
                row_divisions.append((None, (col_pointer, col_pointer + gap_size)))
                col_pointer += gap_size
        stacked = np.hstack(row)
        rows.append(stacked)
        if y < len(assembly.segments)-1 and gap_size > 0:
            gap_matrix = np.empty((gap_size, stacked.shape[1]))
            gap_matrix[:] = np.nan
            rows.append(gap_matrix)
    assembled = np.vstack(rows)

    plot_regions = [ 
        PlotRegion(
            bin_range = (start, end),
            plot_region_type = PlotRegionType.GAP if seg is None else PlotRegionType.DATA_SEGMENT,
            matched_segment = seg
        )
        for seg, (start, end) in row_divisions
    ]

    return AssembledHic(
        data = assembled,
        plot_regions = plot_regions
    )


def get_hic_for_segment_pair(
        hic: HiCFile,
        segment_pair: BinAlignedSegmentPair, # coords must be bin-aligned to resolution
        norm: str = "NONE",
        measure: str = "observed"
    ) -> NDArray:
    """
    Args:
        - hic: HiC handle
        - segment_pair: pair of bin-aligned segments with same resolution
        - norm (default "NONE"): also "VC" or "VC_SQRT" or "KR"
        - measure (default "observed"): also "oe" or "expected"

    Gives raw Hi-C matrix output (not normalized)

    Fails if any coordinates are outside range (0, chrom_size).
    """
    
    # Ensure matrix pair chromosomes are ordered
    transpose = False 
    retrieval_pair = segment_pair
    if CHROM_INDICES[segment_pair.segmentY.chr] > CHROM_INDICES[segment_pair.segmentX.chr]:
        retrieval_pair = switch_segments(segment_pair)
        transpose = True

    # Get matrix zoom data
    mzd = hic.getMatrixZoomData(
        chr_unprefix(retrieval_pair.segmentY.chr),
        chr_unprefix(retrieval_pair.segmentX.chr),
        measure,
        norm, 
        "BP",
        retrieval_pair.resolution,
    )

    # Get matrix data
    data = mzd.getRecordsAsMatrix(
        retrieval_pair.segmentY.start,
        retrieval_pair.segmentY.end-1,
        retrieval_pair.segmentX.start,
        retrieval_pair.segmentX.end-1,
    )

    if transpose:
        data = data.T
    
    return data

def assemble_hic(
    hic: HiCFile,
    assembly: Assembly,
    norm: str = "NONE",
    measure: str = "observed",
    gap_size: int = 0,
) -> NDArray:
    """
    Args:
        - hic: HiC handle
        - assembly: list of bin-aligned segments with same resolution
        - norm (default "NONE"): also "VC" or "VC_SQRT" or "KR"
        - measure (default "observed"): also "oe" or "expected"
        - gap_size (default 0): number of bins of nan to insert between segments
    """

    # Double check all resolutions match 
    resolution = None 
    for segment in assembly.segments:
        if resolution is None:
            resolution = segment.resolution
        assert segment.resolution == resolution

    # Get matrix data for each pair of segments; assemble row by row 
    rows = [] 
    for y, segmentY in enumerate(assembly.segments):
        row = []
        row_divisions = []
        col_pointer = 0 
        for x, segmentX in enumerate(assembly.segments):
            pair_data = get_hic_for_segment_pair(hic, pair_segments(segmentX, segmentY), norm, measure)
            if segmentX.segment.is_minus_strand:
                pair_data = np.flip(pair_data, axis=1)
            if segmentY.segment.is_minus_strand:
                pair_data = np.flip(pair_data, axis=0)
            row.append(pair_data)
            row_divisions.append((segmentX, (col_pointer, col_pointer + pair_data.shape[1])))
            col_pointer += pair_data.shape[1]
            if x < len(assembly.segments)-1 and gap_size > 0:
                gap_matrix = np.empty((pair_data.shape[0], gap_size))
                gap_matrix[:] = np.nan
                row.append(gap_matrix)
                row_divisions.append((None, (col_pointer, col_pointer + gap_size)))
                col_pointer += gap_size
        stacked = np.hstack(row)
        rows.append(stacked)
        if y < len(assembly.segments)-1 and gap_size > 0:
            gap_matrix = np.empty((gap_size, stacked.shape[1]))
            gap_matrix[:] = np.nan
            rows.append(gap_matrix)
    assembled = np.vstack(rows)

    plot_regions = [ 
        PlotRegion(
            bin_range = (start, end),
            plot_region_type = PlotRegionType.GAP if seg is None else PlotRegionType.DATA_SEGMENT,
            matched_segment = seg
        )
        for seg, (start, end) in row_divisions
    ]

    return AssembledHic(
        data = assembled,
        plot_regions = plot_regions
    )



def plot_assembled_triangle(assembled: AssembledHic, ax: plt.Axes, vmax: float | None =None, aspect="equal", mask=None, mask_value=np.nan, rasterized=False, cmap=REDMAP, plot_points: list[DisjointFeatureCall]=[], plot_point_kwargs={"ec": "blue"},):

    data = assembled.data

    resolution = assembled.plot_regions[0].matched_segment.resolution
    if mask is not None:
        chrA = mask.chrA 
        chrB = mask.chrB
        posA = bin_align_value(mask.posA, resolution)
        posB = bin_align_value(mask.posB, resolution)

        has_gaps = len(list(p for plot_region in assembled.plot_regions if plot_region.plot_region_type == PlotRegionType.GAP)) > 0
        if has_gaps:
            iterator = zip(assembled.plot_regions[::2], assembled.plot_regions[2::2])
        else:
            iterator = zip(assembled.plot_regions, assembled.plot_regions[1:])
        for regA, regB in iterator:
            if regA.matched_segment.segment.chr == chrB and regB.matched_segment.segment.chr == chrA:
                chrA, chrB = chrB, chrA
                posA, posB = posB, posA
            if regA.matched_segment.segment.chr == chrA and regB.matched_segment.segment.chr == chrB:

                if not regA.matched_segment.segment.is_minus_strand:
                    bin_diff = (regA.matched_segment.segment.end - posA ) // resolution
                    regEnd = regA.bin_range[1]
                    if bin_diff > 0: 
                        data[:, regEnd-bin_diff:regEnd] = mask_value
                        data[regEnd-bin_diff:regEnd, :] = mask_value
                else:
                    bin_diff = (posA - regA.matched_segment.segment.start) // resolution
                    regStart = regA.bin_range[0]
                    if bin_diff > 0: 
                        data[:, regStart:regStart+bin_diff] = mask_value
                        data[regStart:regStart+bin_diff, :] = mask_value

                if not regB.matched_segment.segment.is_minus_strand:
                    bin_diff = (posB - regB.matched_segment.segment.start) // resolution
                    regStart = regB.bin_range[0]
                    if bin_diff > 0:
                        data[regStart:regStart+bin_diff, :] = mask_value
                        data[:, regStart:regStart+bin_diff] = mask_value
                else:
                    bin_diff = (regB.matched_segment.segment.end - posB) // resolution
                    regEnd = regB.bin_range[1]
                    if bin_diff > 0:
                        data[:, regEnd-bin_diff:regEnd] = mask_value
                        data[regEnd-bin_diff:regEnd, :] = mask_value

    def align_plot_point(p: DisjointFeatureCall):
        mat_pos = []
        cum_bins = 0 
        for plot_region in assembled.plot_regions:
            bin_start, bin_end = plot_region.bin_range
            if plot_region.matched_segment:
                s = plot_region.matched_segment.segment
                if s.chr == p.chrA:
                    new_pos = (p.posA - s.start) / (s.end + resolution - s.start)
                    new_pos = 1 - new_pos if s.is_minus_strand else new_pos 
                    new_pos = new_pos * (bin_end - bin_start) + cum_bins
                    mat_pos.append(new_pos)
                elif s.chr == p.chrB:
                    new_pos = (p.posB - s.start) / (s.end + resolution - s.start)
                    new_pos = 1 - new_pos if s.is_minus_strand else new_pos 
                    new_pos = new_pos * (bin_end - bin_start) + cum_bins
                    mat_pos.append(new_pos)
            cum_bins +=  bin_end - bin_start
        assert(len(mat_pos) == 2)
        return data.shape[0] - np.min(mat_pos), np.max(mat_pos)
                    
    aligned_plot_points = []
    if len(plot_points) > 0:
        for point in plot_points: 
            py, px = align_plot_point(point)
            aligned_plot_points.append((px, py))

    dim = data.shape[0]
    col, row = np.meshgrid(np.arange(dim), np.arange(dim))
    xy_offsets = np.stack([col.flatten(), row.flatten()]).T
    
    # Invert the y axis
    xy_offsets[:,1] = dim - xy_offsets[:,1] 
    
    # Calculate offsets with 45 degree rotation
    cos45 = np.cos(np.pi/4)
    affine = matplotlib.transforms.Affine2D().translate(0, -data.shape[0]).rotate_around(0, 0, np.pi/4).scale(cos45, cos45)
    xy_offsets = affine.transform(xy_offsets)
    
    col = xy_offsets[:, 0].reshape((dim, dim))
    row = xy_offsets[:, 1].reshape((dim, dim))

    vmax = np.nanmax(data) if vmax is None else vmax
    # cmap = cmap.copy()
    # cmap.set_under("#eee")
    ax.pcolormesh(col, row, data, cmap=cmap, vmax=vmax, vmin=0, rasterized=rasterized)

    # Add plot points
    if len(aligned_plot_points) > 0:
        # print(aligned_plot_points)
        aligned_plot_points = affine.transform(np.array(aligned_plot_points))
        # print(aligned_plot_points)
        for (x, y) in aligned_plot_points:
            ellipse = matplotlib.patches.Ellipse((x, y), dim/20, dim/20, fc="none", **plot_point_kwargs)
            ax.add_patch(ellipse)
        

    # Set axis limits
    ax.set_xlim(0-0.5, data.shape[0]-1+0.5)
    ax.set_ylim(0, (data.shape[1])/2)
    ax.set_aspect(aspect)

    # Add anchored size bar in top right corner of plot
    all_resolutions = [seg.matched_segment.resolution for seg in assembled.plot_regions if seg.plot_region_type == PlotRegionType.DATA_SEGMENT]
    assert len(set(all_resolutions)) == 1
    resolution = assembled.plot_regions[0].matched_segment.resolution
    scalebar = AnchoredSizeBar(
        ax.transData,
        0.5, int_to_resolution(resolution), 'upper right', 
        pad=0.1,
        color='black',
        frameon=False,
        size_vertical=0.5,
    )
    ax.add_artist(scalebar)
    blank_axis(ax)


def plot_region_plot_gene_track(plot_region: PlotRegion, ax: plt.Axes, **kwargs):
    if plot_region.plot_region_type == PlotRegionType.GAP:
        return 
    segment = plot_region.matched_segment
    chr = segment.chr 
    start = segment.start 
    end = segment.end
    is_minus_strand = segment.is_minus_strand

    plot_gene_track(chr, start, end, ax=ax, **kwargs)

    # Invert x axis if minus strand
    if is_minus_strand:
        ax.invert_xaxis()

def plot_cnv_track(
        sample_id: str,
        chr: str,
        ax=None,
        vertical=False,
        cnv_lim: tuple[float,float]=None,
        linewidth=3,
        dot_alpha=0.5,
        dot_size=0.5,
        plot_scatter=False,
        plot_segments=True,
        cnv_res=25000,
        locus_lim: tuple[int, int]=None,
        show_zero_line=False,
):

    # Read
    columns = ["chr", "start", "end", "value"]
    df_profile = pd.read_csv(f"../data/cnv_profile/{cnv_res}/{sample_id}_CNV_profile_{cnv_res}.bedgraph", sep="\s+", names=columns)
    df_segment = pd.read_csv(f"../data/cnv_segment/{cnv_res}/{sample_id}_CNV_segment_{cnv_res}.bedgraph", sep="\s+", names=columns)

    subset_profile = df_profile[(df_profile.chr == chr)]
    subset_segment = df_segment[(df_segment.chr == chr)]

    # Plot scatter plot of subset_profile
    positions = subset_profile[subset_profile.value>0].start
    cnv_values = np.log2(subset_profile[subset_profile.value>0].value)

    cmin = cnv_values.min()
    cmax = cnv_values.max()
    cabs = max(abs(cmin), abs(cmax))
    cmap = LinearSegmentedColormap.from_list('rg',["red", "gray", "green"], N=256)
    cvalues = [ cmap((v+2) / 4) for v in cnv_values]

    if ax is None:
        ax = plt.gca()

    if show_zero_line:
        ax.axhline(0, ls=":", color="gray")

    if plot_scatter:
        if vertical:
            ax.scatter(cnv_values, positions, alpha=dot_alpha, c=cvalues, marker=".", s=dot_size, rasterized=True)
        else:
            ax.scatter(positions, cnv_values, alpha=dot_alpha, c=cvalues, marker=".", s=dot_size, rasterized=True)

    if plot_segments:
        minlen = 1
        res = cnv_res
        for _, (_, start, end, value) in subset_segment.iterrows():
            si = start // res
            ei = end // res
            if ei - si >= minlen:
                tmp = subset_profile.value[si:ei]
                mask = tmp == 0
                zero_ratio = mask.sum() / mask.size
                if zero_ratio > 0.80:
                    continue
                seg_cnv_value = np.log2(np.median(tmp[tmp!=0]))
                color=cmap((seg_cnv_value + 2)/4)
                if vertical:
                    ax.vlines(seg_cnv_value, start, end, color=color, linewidth=linewidth)
                else:
                    ax.hlines(seg_cnv_value, start, end, color=color, linewidth=linewidth)

    if vertical:
        if locus_lim is not None:
            ax.set_ylim(*locus_lim)
        else:
            ax.set_ylim(0, CHROM_SIZES[chr])
        ax.set_yticks([])
        ax.set_xlabel("$\log_2$(CN)")
        if cnv_lim is not None:
            ax.set_xlim(*cnv_lim)
    else:
        if locus_lim is not None:
            ax.set_xlim(*locus_lim)
        else:
            ax.set_xlim(0, CHROM_SIZES[chr])
        ax.set_xticks([])
        ax.set_ylabel("$\log_2$(CN)")
        if cnv_lim is not None:
            ax.set_ylim(*cnv_lim)

    if vertical:
        ax.invert_yaxis()
        ax.invert_xaxis()

    return ax

def plot_track(ax: plt.Axes, plot_region: PlotRegion, track: TrackSpec, gap_color: str = REDMAP.get_bad(), is_first: bool = False, is_last: bool = False):
    blank_axis(ax)
    if plot_region.plot_region_type == PlotRegionType.GAP:
        if track.track_type == TrackType.TRACK_CHR_ARROWS:
            ax.plot([0, 1], [0.5, 0.5], color="black", linestyle=":")
            ax.set_ylim(0, 0.5)
            ax.set_xlim(0, 1)
        else:
            ax.set_facecolor(gap_color)
        return 

    segment = plot_region.matched_segment.segment
    resolution = plot_region.matched_segment.resolution
    chr = segment.chr
    start = segment.start
    end = segment.end
    is_minus_strand = segment.is_minus_strand

    match track.track_type: 
        case TrackType.TRACK_GENES:
            plot_gene_track(chr, start, end, ax=ax, **track.track_kwargs)
            if ("centered_names" not in track.track_kwargs.keys() or not track.track_kwargs["centered_names"]) and is_minus_strand:
                # Change ha of each text object: ha=left, then change to right and vice versa
                for text in ax.texts:
                    if text.get_ha() == "left":
                        text.set_ha("right")
                    elif text.get_ha() == "right":
                        text.set_ha("left")
        case TrackType.TRACK_BIGWIG:
            if not "bw_handle" in track.track_kwargs.keys():
                raise ValueError("Missing bigwig handle in track_kwargs")
            kwargs = track.track_kwargs.copy()
            bw_handle = kwargs["bw_handle"]
            del kwargs["bw_handle"]
            if "ymax" in kwargs:
                ymax = kwargs["ymax"]
                del kwargs["ymax"]
            else:
                ymax = None
            if "show_ylim" in kwargs and ymax is not None:
                show_ylim = kwargs["show_ylim"]
                del kwargs["show_ylim"]
            else:
                show_ylim = False
            plot_bigwig_track(bw_handle, chr, start, end, ax=ax, **kwargs)
            if "vertical" not in kwargs:
                ax.set_ylim((0, ymax))
                if show_ylim and is_first and ymax is not None:
                    ax.yaxis.set_visible(True)
                    ax.set_yticks([0, ymax], [0, ymax], fontsize=8)
                    ax.get_yticklabels()[0].set_va("bottom")
                    ax.get_yticklabels()[-1].set_va("top")
                    ax.spines[["left"]].set_visible(True)
            else:
                ax.set_xlim((ymax, 0))
                if show_ylim and is_first and ymax is not None:
                    ax.xaxis.set_visible(True)
                    ax.set_xticks([ymax, 0], [ymax, 0], fontsize=8)
                    ax.get_xticklabels()[0].set_ha("left")
                    ax.get_xticklabels()[-1].set_ha("right")
                    ax.spines[["bottom"]].set_visible(True)
                    ax.xaxis.tick_bottom()
        case TrackType.TRACK_BED:
            raise NotImplementedError("BED track plotting not yet implemented")
        case TrackType.TRACK_BED_DIRECTIONAL:
            necessary_kwargs = ["bed_file"]
            for kwarg in necessary_kwargs:
                if not kwarg in track.track_kwargs.keys():
                    raise ValueError(f"Missing {kwarg} in track_kwargs")
            kwargs = track.track_kwargs.copy()
            bed_file = kwargs["bed_file"]
            del kwargs["bed_file"]
            color = "seagreen" if not "color" in kwargs.keys() else kwargs["color"]
            motifs = pd.read_csv(bed_file, sep="\s+", header=None, names=["chr", "start", "end", "strand"])
            motifs = motifs[motifs["chr"] == chr]
            motifs = motifs[(motifs["start"] >= start) & (motifs["start"] <= end) | (motifs["end"] >= start) & (motifs["end"] <= end)]
            style = "scatter" if not "style" in kwargs.keys() else kwargs["style"]
            s = 20 if not "s" in kwargs.keys() else kwargs["s"]
            if style == "scatter":
                pos = motifs[motifs["strand"] == "+"].start
                neg = motifs[motifs["strand"] == "-"].end
                marker_pos = ">"
                marker_neg = "<"
                if is_minus_strand:
                    marker_pos, marker_neg = marker_neg, marker_pos
                ax.scatter(pos, [0.5] * len(pos), color=color, s=s, marker=marker_pos)
                ax.scatter(neg, [0.5] * len(neg), color=color, s=s, marker=marker_neg)
            elif style == "interval":
                head_length = ((end - start) / 25) if not "head_length" in kwargs.keys() else kwargs["head_length"]
                for i, row in motifs.iterrows():
                    rectangle = matplotlib.patches.Rectangle((row["start"], 0), row["end"] - row["start"], 1, linewidth=1, color=color)
                    ax.add_patch(rectangle)
                    # Add either a left or right arrow depending on strand 
                    if row["strand"] == "+":
                        ax.arrow(row["start"], 0.5, 0.5, 0, ec=color, fc="none", head_length=head_length, head_width=0.5)
                    elif row["strand"] == "-":
                        ax.arrow(row["end"], 0.5, -0.5, 0, ec=color, fc="none", head_length=head_length, head_width=0.5)
            else:
                raise ValueError(f"Unknown style {style}")
            ax.set_xlim(start, end)
            ax.set_ylim(0, 1)
        case TrackType.TRACK_CHR_ARROWS:
            arrowhead_style = "left"
            head_length = ((end - start) / 25)
            label_chr = False if not "label_chr" in track.track_kwargs.keys() else track.track_kwargs["label_chr"]
            adjust_with_sv = None if not "adjust_with_sv" in track.track_kwargs.keys() else track.track_kwargs["adjust_with_sv"]

            arrow_end = end
            arrow_start = start
            dotted_start = None
            dotted_end = None
            if adjust_with_sv is not None:
                # Assume the SV pos is closest to the end to truncate
                if adjust_with_sv.chrA == chr:
                    pos = adjust_with_sv.posA
                elif adjust_with_sv.chrB == chr:
                    pos = adjust_with_sv.posB
                else:
                    raise ValueError("SV does not match track chromosome")
                pos = bin_align_value(pos, resolution)
                if is_first:
                    if not is_minus_strand:
                        if pos < end:
                            arrow_end = pos
                            dotted_start = pos 
                            dotted_end = end
                    elif is_minus_strand:
                        if pos > start:
                            arrow_start = pos
                            dotted_start = start
                            dotted_end = pos
                elif is_last:
                    if not is_minus_strand:
                        if pos > start:
                            arrow_start = pos
                            dotted_start = start
                            dotted_end = pos
                    elif is_minus_strand:
                        if pos < end:
                            arrow_end = pos
                            dotted_start = pos
                            dotted_end = end
            ax.arrow(arrow_start, 0.5, arrow_end-arrow_start, 0, ec="none", fc=CHROM_COLORS[chr], width=0.5, head_length=head_length, head_width=1, shape=arrowhead_style, length_includes_head=True)
            if dotted_start is not None and dotted_end is not None:
                ax.plot([dotted_start, dotted_end], [0.5, 0.5], color="black", linestyle=":")
            ax.set_xlim(start, end)
            if label_chr:
                ax.set_ylim(-1, 0.5)
                ax.text(start, -1, f"{start:,}", ha="left", va="bottom")
                ax.text(end, -1, f"{end:,}", ha="right", va="bottom")
                ax.text((start + end) / 2, -1, f"{chr}", ha="center", va="bottom")
                # If is minus strand, flip all ha of text
                if is_minus_strand:
                    for text in ax.texts:
                        if text.get_ha() == "left":
                            text.set_ha("right")
                        elif text.get_ha() == "right":
                            text.set_ha("left")
            else:
                ax.set_ylim(0, 0.5)
        case TrackType.TRACK_VIRTUAL_4C:
            # Get hic data 
            necessary_kwargs = ["hic_file", "anchor_chr", "anchor_pos"]
            for kwarg in necessary_kwargs:
                if not kwarg in track.track_kwargs.keys():
                    raise ValueError(f"Missing {kwarg} in track_kwargs")
            kwargs = track.track_kwargs.copy()
            hic_file = kwargs["hic_file"]

            anchor_chr = track.track_kwargs["anchor_chr"]
            anchor_pos = track.track_kwargs["anchor_pos"]
            anchor_bin_radius = 0 if not "anchor_bin_radius" in track.track_kwargs.keys() else track.track_kwargs["anchor_bin_radius"]
            anchor_pos_bin = bin_align_value(anchor_pos, resolution)
            anchor_start = anchor_pos_bin - (anchor_bin_radius * resolution)
            anchor_end = anchor_pos_bin + resolution + (anchor_bin_radius * resolution)

            norm = "NONE" if not "norm" in kwargs.keys() else kwargs["norm"]
            measure = "observed" if not "measure" in kwargs.keys() else kwargs["measure"]

            cmap = "Purples" if not "cmap" in kwargs.keys() else kwargs["cmap"]

            adjust_with_sv = None if not "adjust_with_sv" in kwargs.keys() else kwargs["adjust_with_sv"]

            data = get_hic_for_segment_pair(
                hic_file,
                BinAlignedSegmentPair(
                    segmentX = segment,
                    segmentY = Segment(anchor_chr, anchor_start, anchor_end),
                    resolution = resolution
                ),
                norm=norm,
                measure=measure,
            ).mean(axis=0)

            if adjust_with_sv is not None:
                # Assume the SV pos is closest to the end to truncate
                if adjust_with_sv.chrA == chr:
                    pos = adjust_with_sv.posA
                elif adjust_with_sv.chrB == chr:
                    pos = adjust_with_sv.posB
                else:
                    raise ValueError("SV does not match track chromosome")
                pos_bin = bin_align_value(pos, resolution)
                if is_first:
                    if not is_minus_strand:
                        bin_diff = (end - pos_bin) // resolution
                        if bin_diff > 0:
                            data[-bin_diff:] = np.nan
                    elif is_minus_strand:
                        bin_diff = (pos_bin - start) // resolution
                        if bin_diff > 0:
                            data[:bin_diff] = np.nan
                elif is_last:
                    if not is_minus_strand:
                        bin_diff = (pos_bin - start) // resolution
                        if bin_diff > 0:
                            data[:bin_diff] = np.nan
                    elif is_minus_strand:
                        bin_diff = (end - pos_bin) // resolution
                        if bin_diff > 0:
                            data[-bin_diff:] = np.nan
            if isinstance(cmap, str):
                cmap = plt.get_cmap(cmap)
            cmap.set_bad(REDMAP.get_bad())
            ax.matshow(data.reshape(1, -1), cmap=cmap, aspect="auto")
            blank_axis(ax)
        case TrackType.TRACK_BLANK:
            blank_axis(ax)

    if is_minus_strand:
        ax.invert_xaxis()

    ax.yaxis.set_visible(True)
    ax.set_ylabel(track.ylabel, rotation=0, ha="right", va="center")


    return 

def plot_triangle_and_tracks(sample_id: str , assembled: AssembledHic, tracks: list[TrackSpec], figsize=None, vmax=None, size_base=2, plot_points=[], plot_point_kwargs={"ec": "blue"}):
    width_ratios = [abs(r.bin_range[1] - r.bin_range[0]) for r in assembled.plot_regions]
    height_ratios_tracks = [t.rel_height for t in tracks]
    if figsize is None:
        size_base = size_base
        figheight = 2.5 + sum(height_ratios_tracks)
        figwidth = 5
        figsize = (size_base*figwidth, size_base*figheight)
    height_ratios = [2.5] + height_ratios_tracks
    fig = plt.figure(figsize=figsize)
    nrows = 1 + len(tracks)
    ncols = len(assembled.plot_regions)
    gs = fig.add_gridspec(nrows, ncols, wspace=0, hspace=0, height_ratios=height_ratios, width_ratios=width_ratios)

    ax_triangle = fig.add_subplot(gs[0, :])
    if vmax is None:
        vmax = np.nanmax(assembled.data) / 4
    plot_assembled_triangle(assembled, ax_triangle, vmax=vmax, rasterized=True, aspect="auto", plot_points=plot_points, plot_point_kwargs=plot_point_kwargs)
    ax_triangle.set_ylabel(sample_id, rotation=0, ha="right", va="center")

    for row, track in enumerate(tracks):
        for col, plot_region in enumerate(assembled.plot_regions):
            ax = fig.add_subplot(gs[row+1, col])
            is_first = col == 0
            plot_track(ax, plot_region, track, is_first=is_first)
            if col > 0:
                ax.set_ylabel("")
            
    return fig

DEFAULT_TRACKS = [
    TrackSpec(
        TrackType.TRACK_GENES,
        ylabel="Genes",
        track_kwargs={
            "max_rows": 10
        },
        rel_height=1,
    ),
    TrackSpec(
        TrackType.TRACK_CHR_ARROWS,
        track_kwargs={
            "label_chr": True
        },
        ylabel="Region",
        rel_height=0.2,
    ),
]

def make_assembly_for_sv(sv: BreakfinderCall | CuratedCall, res=10000, width=20*25000) -> Assembly:
    if isinstance(sv, BreakfinderCall):
        chrA = sv.breakpointA.chr
        posA = sv.breakpointA.pos
        strandA = sv.breakpointA.strand
        chrB = sv.breakpointB.chr
        posB = sv.breakpointB.pos
        strandB = sv.breakpointB.strand
    elif isinstance(sv, CuratedCall) or isinstance(sv, EagleCall):
        chrA = sv.chrA
        posA = sv.posA
        strandA = sv.strandA
        chrB = sv.chrB
        posB = sv.posB
        strandB = sv.strandB

    if strandA == Strand.POS: 
        startA = posA - width 
        endA = posA
        is_minus_strand_A = False
    else:
        startA = posA
        endA = posA + width
        is_minus_strand_A = True
    if strandB == Strand.NEG: 
        startB = posB
        endB = posB + width
        is_minus_strand_B = False
    else:
        startB = posB - width
        endB = posB
        is_minus_strand_B = True

    segmentA = bin_align_segment(Segment(chrA, startA, endA, is_minus_strand=is_minus_strand_A), res)
    segmentB = bin_align_segment(Segment(chrB, startB, endB, is_minus_strand=is_minus_strand_B), res)

    assembly = Assembly([
        segmentA,
        segmentB,
    ])

    return assembly



def to_arima_pipeline_sample(sample_id: str, profile: DataProfile) -> ArimaPipelineSample:
    if profile.curated is not None:
        return read_sample(sample_id, profile.hic, profile.qc, profile.curated)
    else:
        return read_sample(sample_id, profile.hic, profile.qc, profile.breakfinder)

def sv_involves_gene(sv: CuratedCall, gene: str, width=1000000):
    if isinstance(sv, BreakfinderCall):
        sv = breakfinder_to_curated(sv)
    nearbyA, nearbyB = show_nearby_genes_stranded(sv, width=width, protein_coding=False)
    if gene in nearbyA or gene in nearbyB:
        return True
    return False