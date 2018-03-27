
import os as path
import sys
import operator
import configparser
import pysam
import src.umi_network_collapse as network

def umi_count(contig, region_start, region_end, bam_file):
    """Returns tally of UMIs in given region."""

    umi_counts = {}
    
    with pysam.AlignmentFile(bam_file, "rb") as bam_reader:

        for read in bam_reader.fetch(contig, region_start, region_end):

            umi = read.query_name.split(':')[-1]

            if umi in umi_counts:
                umi_counts[umi] += 1
            else:
                umi_counts[umi] = 1

    return umi_counts


def group_position(contig, region_start, region_end, bam_file, umi_groups, pos_threshold):
    """Further groups umi_groups by position."""

    umi_positions = {umi_group: {'r1': {}, 'r2': {}} for umi_group in umi_groups}
    
    with pysam.AlignmentFile(bam_file, "rb") as bam_reader:

        for read in bam_reader.fetch(contig, region_start, region_end):

            umi = read.query_name.split(':')[-1]
            flag = read.flag 
            is_r1 = (flag & 64) >> 6 ## 1 if first in pair, else 0
            read_type = 'r1' if is_r1 else 'r2'
            pos = read.reference_start

            for umi_group in umi_groups:
                if umi in umi_group:
                    if pos in umi_positions[umi_group][read_type]:
                        umi_positions[umi_group][read_type][pos] += 1
                    else:
                        umi_positions[umi_group][read_type][pos] = 1

    families = {}
    for umi_group in umi_positions:
        for read_type in ('r1', 'r2'):
            while umi_positions[umi_group][read_type]:

                ## Add UMI group + most common position to families
                most_frequent_pos = max(umi_positions[umi_group][read_type].items(), key=operator.itemgetter(1))[0]
                families[(umi_group, most_frequent_pos)] = 1

                for pos in umi_positions[umi_group][read_type].copy():
                    if abs(most_frequent_pos - pos) <= pos_threshold:
                        families[(umi_group, most_frequent_pos)] += 1
                        del umi_positions[umi_group][read_type][pos]

    return families


def get_umi_families(contig, region_start, region_end, bam_file, config):
    """
    Returns a list of (umi_group, pos) pairs representing
    error-corrected UMI families.
    """ 

    pos_threshold = int(config['SETTINGS']['umi_family_pos_threshold']) if config else 10

    counts = umi_count(contig, region_start, region_end, bam_file)
    umis = counts.keys()

    clusterer = network.UMIClusterer(cluster_method="directional")
    umi_groups = clusterer(umis, counts, config)

    families = group_position(contig, region_start, region_end, bam_file, umi_groups, pos_threshold)

    return families


