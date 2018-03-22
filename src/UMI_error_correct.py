
import os as path
import sys
import operator
import argparse
import configparser
import itertools
from handle_args import handle_arg
from _edit_distance import edit_distance

class UMIClusterer:
    """
    A functor that clusters a dictionary of UMIs and their counts.
    The primary return value is a list of representative UMIs.
    
    Adapted from https://github.com/CGATOxford/UMI-tools/blob/master/umi_tools/network.py
    """
    
    def __init__(self, cluster_method="directional"):
        """Select the required class methods for the cluster_method."""
        
        self.max_umis_per_position = 0
        self.total_umis_per_position = 0
        self.positions = 0
        
        if cluster_method == "directional":
            self.get_adj_list = self._get_adj_list_directional
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_groups = self._group_adjacency
        

    def __call__(self, umis, counts, threshold):
        """Counts is a dictionary that maps UMIs to their counts."""
        
        umis = list(umis)
        
        self.positions += 1
        
        number_of_umis = len(umis)
        self.total_umis_per_position += number_of_umis
        
        if number_of_umis > self.max_umis_per_position:
            self.max_umis_per_position = number_of_umis
            
        len_umis = [len(x) for x in umis]
        
        assert max(len_umis) == min(len_umis), (
            "Not all UMIs are the same length(!): {} - {}"
            .format(min(len_umis), max(len_umis)))
        
        adj_list = self.get_adj_list(umis, counts, threshold)
        clusters = self.get_connected_components(umis, adj_list, counts)
        final_umis = [list(x) for x in self.get_groups(clusters, adj_list, counts)]
        
        return final_umis
        
if __name__=='__main__':
    
    parser = argparse.ArgumentParser()
    #parser.add_argument('-bam', '--bam_file',    help='Path to your BAM file.')

