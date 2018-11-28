
import os as path
import sys
import operator
import configparser
import pysam
import functools
import src.umi_network_collapse as network
import csv
from utils.generate_plots import umi_plot
import matplotlib.pyplot as plt
import pandas as pd
import fnmatch
import itertools

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
        for fampos in self.families:
            if abs(pos - fampos) <= pos_threshold:
                return fampos

        return None

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

    print("Counting UMIs...")
    counts = umi_count(contig, region_start, region_end, bam_file)
    umis = counts.keys()

    print("Clustering UMIs...")
    clusterer = network.UMIClusterer(cluster_method="directional")  
    
    umi_groups = clusterer(umis, counts, config)

    #Writing datafile
    #umi_datafile(umi_groups, contig, region_start, region_end, filename)

    print("Grouping UMIs by position...")
    umi_table = group_position(contig, region_start, region_end, bam_file, umi_groups, pos_threshold)
    
    return umi_table, umi_groups, temp_str




def umi_datafile(umi_groups):

    #Each 'row' in umi_groups is a tuple, where the first element is the parent umi and all following elements in the tuple are the child umis
    elements_per_row = [len(row) for row in umi_groups]
    total_parent_umi_count = len(elements_per_row)
    total_child_umi_count = sum(elements_per_row)-total_parent_umi_count

    parent_umis = [row[0] for row in umi_groups] 
    children_per_parent = [count-1 for count in elements_per_row]

    #Below dictionary has key=parent umi sequence, and value=num of child umis corresponding to that parent umi
    umis_and_numofchildren = {}
    for i in range(len(parent_umis)):
        umis_and_numofchildren[(parent_umis[i])] = children_per_parent[i]
   
    #var is a list of (key,value) pairs where key=number of children (0,1,2,3...) & value=no. of parents having that number of children
    var = [(k, len(list(v))) for k, v in itertools.groupby(sorted(umis_and_numofchildren.values()))]

    num_of_children = ""
    freq_of_parent_umis = ""
    for i in range(len(var)):
        if i == 0:
            num_of_children = str(var[i][0])
            freq_of_parent_umis = str(var[i][1])
        else:
            num_of_children = num_of_children+","+str(var[i][0])
            freq_of_parent_umis = freq_of_parent_umis+","+str(var[i][1])

    return total_parent_umi_count, total_child_umi_count, num_of_children, freq_of_parent_umis



    """
    consdp = str(len(umi_groups))
    elements_per_row = [len(row) for row in umi_groups]
    sum_child_count = 0
    regions_child_count = "" #Comma-separated string containing number of child umis per parent umi

    for child_count in elements_per_row:
        child_count = child_count-1

        if regions_child_count == "":
            regions_child_count = str(child_count)
        else:
            regions_child_count = regions_child_count+","+str(child_count)

        sum_child_count = sum_child_count+child_count

    childdp = str(sum_child_count)

    csvrow = {'CHR' : contig, 'START' : str(region_start), 'END' : str(region_end), 'PTU_NUM' : consdp, 'CTU_NUM' : childdp}
    headers = ['CHR', 'START', 'END', 'PTU_NUM', 'CTU_NUM']

    #df = pd.DataFrame(csvrow, index=[idx])

    csv.register_dialect('myDialect', delimiter='\t', quoting=csv.QUOTE_NONE)

    dirpath = "/u/iwarikoo/Debarcer2/d_output/haloplex_9538005/umifiles/"
    num_umifiles = len(fnmatch.filter(path.listdir(dirpath), '*.umis'))
    print("nu. of umifiles: "+str(num_umifiles))

    file = open(filename, "a")
    with file:
        writer = csv.DictWriter(file, dialect='myDialect', fieldnames=headers)
        writer.writerow(csvrow)

    #df = pd.read_csv(filename)
    #var = df['START']
    #print(str(var[0]))
    #print(str(var[1]))

    start = []
    parent_umis = []
    child_umis = []
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        row_count = sum(1 for row in reader)

        for line in f.readlines():
            array = line.split('\t')
            start.append(array[1])
            parent_umis.append(array[3])
            child_umis.append(array[4])
        f.close()

    print("no. of rows in datafile.csv: " + str(row_count))
    #umi_plot(filename, start, parent_umis)

    temp_str = "hello"
    """

