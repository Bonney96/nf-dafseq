#!/usr/bin/env python3
import pysam
import pandas as pd
import os
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import code
import numpy as np
from scipy.sparse import coo_matrix
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as colors
import matplotlib.cm as cm
import argparse


def parse_region(region):
    """
    Parses a genomic region in the format chr:start-end.
    Returns a tuple of chromosome, start, and end positions.
    """
    chrom, positions = region.split(":")
    start, end = map(int, positions.split("-"))
    return chrom, start, end


def strand_metrics(read, cutoff=0.9, min_deamination_count=50):
    """
    Determines the strand acted upon by DddA based on the proportion of C->T & G->A,
    and counts deamination for each two base pair context (AC, CC, GC, TC).
    For G->A deamination strands, the complement is recorded.

    Args:
        read (pysam.AlignedSegment): A read from a BAM file.
        cutoff (float): Proportion threshold to determine strand type.

    Returns a tuple containing:
        - strand: "CT", "GA", "chimera", "undetermined", or "none"
        - doublets: Dictionary with counts of deamination doublets
        - mutation_count: Number of mutations in the read
        - deamination_pos: List of positions in the read where deamination occurred
    """

    seq = read.query_sequence
    pair = read.get_aligned_pairs(matches_only=False, with_seq=True)

    if not read.has_tag('st'):
        strand='none'
    else: strand=read.get_tag('st')

    mutation_count = -1
    
    # Count and record position of deaminations for CT and GA strands
    if strand in ["CT", "GA"]:

        deamination_pos = []  # deaminated positions in read coordinates
        deamination_ref = []

        for i, pos in enumerate(pair):
            if pos[0] is None or pos[1] is None:  # indel, ignore
                continue

            ref_base = pos[2].upper()

            if ref_base != strand[0]:  # ignore non C/G bases
                continue

            strand_base = seq[pos[0]]
            if strand == 'CT' and strand_base=='Y':
                deamination_pos.append(pos[0])
                deamination_ref.append(pos[1])
            elif strand == 'GA' and strand_base=='R':
                deamination_pos.append(pos[0])
                deamination_ref.append(pos[1])
    else:
        deamination_pos = None
        deamination_ref = None

    return strand, mutation_count, deamination_pos, deamination_ref


def strand_metrics_table(
        psfile, chrom, start, end, output_prefix, max_distance=100, max_error_rate=0.015, cutoff=0.9, min_deamination_count=50, include_readnames=None
):
    """
    Iterates through reads in a specified region of a BAM file,
    calculates strand metrics, and returns a DataFrame with the results.

    Args:
        psfile (pysam.AlignmentFile): The BAM file to read.
        chrom (str): Chromosome name.
        start (int): Start position of the region.
        end (int): End position of the region.
        chimera_cutoff (float): Proportion threshold to determine strand type.
        min_deamination_count (int): Minimum number of deaminations to designate a strand.
        include_readnames (list, optional): List of read names to include in the analysis.
            If None, all reads in the region are included.
    """

    read_collector = []

    # check that this region has reads
    if psfile.count(chrom, start, end) == 0:
        return pd.DataFrame()

    for read in psfile.fetch(chrom, start, end):
        if include_readnames is not None and read.query_name not in include_readnames:
            continue
        if read.is_secondary or read.is_supplementary:
            continue

        if read.reference_start > start + max_distance or  read.reference_end < end - max_distance :
            continue

        strand, mutation_count, deam_pos, ref_pos = strand_metrics(
            read, cutoff=cutoff, min_deamination_count=min_deamination_count
        )

        duplicate = read.get_tag("du") if read.has_tag("du") else "None"

        read_data = {
            "read_name": read.query_name,
            "chrom": chrom,
            "reg_start": start,
            "reg_end": end,
            "strand_start": read.reference_start,
            "strand_end": read.reference_end,
            "length": len(read.query_sequence),
            "strand": strand,
            "duplicate": duplicate,
            "map_qual": read.mapq,
            "mutation_count": mutation_count,
            "deamination_positions": ",".join(map(str, deam_pos)) if deam_pos is not None else "",
            "ref_positions": ",".join(map(str, ref_pos)) if ref_pos is not None else ""
        }

        read_collector.append(read_data)

    reads_table = pd.DataFrame(read_collector)

    # No reads passed the filters for this region (e.g. a pooled amplicon with no
    # full-length spanning reads). Skip the region instead of raising KeyError.
    if reads_table.empty or 'strand' not in reads_table.columns:
        return pd.DataFrame()

    reads_table= reads_table[(reads_table['strand']=="CT") | (reads_table['strand']=="GA")].copy().reset_index()

    # No top/bottom (CT/GA) reads in this region: nothing to cluster, skip it.
    if reads_table.empty:
        return pd.DataFrame()
    positions=[]

    for ii in range(reads_table.shape[0]):
        positions.extend(reads_table['ref_positions'][ii].split(','))

    positions=list(set(positions))
    positions=[int(kk) for kk in positions]
    positions.sort()
    pos_pos={}
    for ii in range(len(positions)):
        pos_pos[positions[ii]]=ii
    sp=[]
    for ii in range(reads_table.shape[0]):
        xx=reads_table['ref_positions'][ii].split(',')
        strand=1
        if reads_table['strand'][ii]=="GA":
            strand=-1
        for xx1 in xx:
            sp.append([ii, pos_pos[int(xx1)], strand])

    sp=np.array(sp)
    sparse_matrix = coo_matrix((sp[:,2], (sp[:,0], sp[:,1])), shape=(max(sp[:,0])+1, max(sp[:,1])+1))
    dense_array = sparse_matrix.toarray()

    linked = linkage(dense_array, method='complete', metric='cityblock') 
    df=pd.DataFrame(dense_array)
    df['read_name']=reads_table['read_name']
    df['cluster']=fcluster(linked, max_error_rate*dense_array.shape[1], criterion='distance')
    df.sort_values(by=['cluster'], inplace=True)
    
    clusters1=df[['cluster', 'read_name']].copy()
    clusters1['gp']=clusters1['cluster']%30
    clusters1['fracs']=clusters1['gp']/np.max(clusters1['gp'])
    norm = colors.Normalize(np.min(clusters1['fracs']), np.max(clusters1['fracs']))
    colors2 = cm.gist_ncar(norm(clusters1['fracs']))
    sns.clustermap(df.iloc[:, :(df.shape[1]-2)], row_cluster=False, col_cluster=False,   figsize=(40,100), row_colors=colors2)

    plt.savefig(output_prefix+'.'+chrom+'_'+str(start)+'_'+str(end)+'.clusters.png')
    plt.close()

    clusters=df[['read_name', 'cluster']].copy()
    clusters=pd.merge(clusters, reads_table[['read_name', 'map_qual']], on='read_name', how='left')
    clusters.sort_values(by=['cluster', 'map_qual'], inplace=True)
    leads=clusters.groupby('cluster')[[ 'read_name']].first()
    leads.rename(columns={'read_name':'lead'}, inplace=True) 
    clusters=pd.merge(clusters, leads, on='cluster', how='left')
    return clusters

def make_dedup_bam(input_bam, output_bam, regions, clusters):
    # Reads that were clustered, and their cluster lead. Empty when no region
    # produced clusters.
    if clusters.empty or 'read_name' not in clusters.columns:
        return
    lead_lookup = dict(zip(clusters['read_name'], clusters['lead']))

    # Fetch every region (not just the last one) and write each alignment once.
    # The (name, flag, pos) key dedups records seen via overlapping amplicons.
    written = set()
    for region in regions:
        chrom, start, end = parse_region(region)
        for read in input_bam.fetch(chrom, start, end):
            if read.query_name not in lead_lookup:
                continue
            key = (read.query_name, read.flag, read.reference_start)
            if key in written:
                continue
            written.add(key)
            lead = lead_lookup[read.query_name]
            if read.query_name != lead:
                read.set_tag("du", lead)
            output_bam.write(read)


## Main script execution
#bam_path = snakemake.input.data
#regions = snakemake.params.regions
#targeting_metrics = snakemake.input.targeting_data
#chimera_cutoff = snakemake.params.chimera_cutoff
#min_deamination_count = snakemake.params.min_deamination_count
#read_metrics = snakemake.output.read_metrics
#summary_metrics = snakemake.output.summary_metrics
#threads = snakemake.threads



parser = argparse.ArgumentParser(description='mark dups for dafseq')
parser.add_argument('-i','--inbam', help='input bam', required=True)
parser.add_argument('-o','--outprefix', help='output prefix', required=True)
parser.add_argument('-b','--regions', help='regions bedfile', required=True)
parser.add_argument('-m','--max_error_rate', help='max error rate', default=0.015, required=False)
parser.add_argument('-d', '--max_distance', help='max distance of read start/end inside region start/end, to filter short reads', default=100, required=False)
args = parser.parse_args()


bam_file = args.inbam
output_prefix = args.outprefix
output_bam = output_prefix + '.mkdup.bam'
max_error_rate=args.max_error_rate
max_distance=args.max_distance

regions=[]
with open (args.regions, 'r') as f:
  for line in f:
      ll=line.strip().split('\t')
      regions.append(ll[0]+':'+ll[1]+'-'+ll[2])


chimera_cutoff = 0.9
min_deamination_count = 50
clusters_list=[]


input_bam = pysam.AlignmentFile(bam_file, "rb")
output_bam = pysam.AlignmentFile(output_bam, "wb", template=input_bam)
tables = []

for region in regions:

    chrom, start, end = parse_region(region)
    clusters = strand_metrics_table(
        input_bam, chrom, start, end, output_prefix, max_distance, max_error_rate, cutoff=chimera_cutoff, min_deamination_count=min_deamination_count
    )
    clusters_list.append(clusters)
    
clusters=pd.concat(clusters_list, ignore_index=True) if clusters_list else pd.DataFrame()

make_dedup_bam(input_bam, output_bam, regions, clusters)

clusters.to_csv(output_prefix+'.clusters.txt', sep="\t", index=False)

