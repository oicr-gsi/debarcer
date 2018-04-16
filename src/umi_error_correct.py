
import os as path
import sys
import operator
import configparser
import pysam
import functools
import src.umi_network_collapse as network

class UMIGroup:
    """Contains position and count info for all UMIs representing one group."""

    def __init__(self, key):
        self.families = {}
        self.key = key

    def addNew(self, pos):
        self.families[pos] = 1

    def add(self, pos):
        self.families[pos] += 1

    @functools.lru_cache(maxsize=4, typed=False)
    def getClosest(self, pos, pos_threshold):
        dist = lambda x, y: abs(x - y)
        closest = [fampos for fampos in self.families if dist(pos, fampos) <= pos_threshold]
        return closest[0] if closest else None

    def key(self):
        return key


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
    """Splits umi_groups into families (umi + position pairs)."""

    umi_table = {}

    ## Initialize the table with references to UMIGroup objs
    for group in umi_groups:
        new_group = UMIGroup(key=group[0])
        for umi in group:
            umi_table[umi] = new_group

    with pysam.AlignmentFile(bam_file, "rb") as bam_reader:

        for read in bam_reader.fetch(contig, region_start, region_end):

            umi = read.query_name.split(':')[-1]
            pos = read.reference_start

            if umi in umi_table:

                umi_group = umi_table[umi]
                families = umi_group.families

                if not families:
                    umi_group.addNew(pos)

                elif pos in families:
                    umi_group.add(pos)

                else:
                    closest = umi_table[umi].getClosest(pos, pos_threshold)
                    if closest:
                        umi_group.add(closest)
                    else:
                        umi_group.addNew(pos)

    return umi_table


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

    umi_table = group_position(contig, region_start, region_end, bam_file, umi_groups, pos_threshold)

    return umi_table


