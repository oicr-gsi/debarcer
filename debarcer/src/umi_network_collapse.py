
import os as path
import sys
import operator
import configparser
import itertools
import collections

from src.edit_distance import edit_distance
 
"""
Adapted from https://github.com/CGATOxford/UMI-tools/blob/master/umi_tools/network.py
(Retrieved March 21, 2018)
"""

sys.setrecursionlimit(10000)

def breadth_first_search(node, adj_list):
    searched = set()
    queue = set()
    queue.update((node,))
    searched.update((node,))

    while len(queue) > 0:
        node = queue.pop()
        for next_node in adj_list[node]:
            if next_node not in searched:
                queue.update((next_node,))
                searched.update((next_node,))

    return searched


def recursive_search(node, adj_list):
    children = adj_list[node]

    children = [x for x in children if x not in recursive_search.component]
    for child in children:
        recursive_search.component.update((child,))
        recursive_search.component.update(
            recursive_search(child, adj_list))
    return recursive_search.component


def breadth_first_search_recursive(node, adj_list):
    try:
        recursive_search.component = set((node,))
        return recursive_search(node, adj_list)

    except RecursionError as error:
        print('Recursion Error: %s' % error)
        return breadth_first_search(node, adj_list)


def remove_umis(adj_list, cluster, nodes):
    """
    Removes the specified nodes from the cluster 
    and returns the remaining nodes.
    """

    # list incomprehension: for x in nodes: for node in adj_list[x]: yield node
    nodes_to_remove = set([node
                           for x in nodes
                           for node in adj_list[x]] + nodes)

    return cluster - nodes_to_remove


def get_substr_slices(umi_length, idx_size):
    """
    Create slices to split a UMI into approximately equal size substrings
    Returns a list of tuples that can be passed to slice function.
    """
    cs, r = divmod(umi_length, idx_size)
    sub_sizes = [cs + 1] * r + [cs] * (idx_size - r)
    offset = 0
    slices = []
    for s in sub_sizes:
        slices.append((offset, offset + s))
        offset += s
    return slices


def build_substr_idx(umis, umi_length, min_edit):
    """
    Build a dictionary of nearest neighbours using substrings, can be used
    to reduce the number of pairwise comparisons.
    """
    substr_idx = collections.defaultdict(
        lambda: collections.defaultdict(set))
    slices = get_substr_slices(umi_length, min_edit + 1)
    for idx in slices:
        for u in umis:
            u_sub = u[slice(*idx)]
            substr_idx[idx][u_sub].add(u)
    return substr_idx


def iter_nearest_neighbours(umis, substr_idx):
    """
    Use substring dict to get (approximately) all the nearest neighbours to
    each in a set of UMIs.
    """
    for u in umis:
        neighbours = set()
        for idx, substr_map in substr_idx.items():
            u_sub = u[slice(*idx)]
            neighbours = neighbours.union(substr_map[u_sub])
        neighbours.remove(u)
        for nbr in neighbours:
            yield u, nbr




class UMIClusterer:
    """
    A functor that clusters a dictionary of UMIs and their counts.
    The primary return value is a list of representative UMIs.    
    """

    def _get_adj_list_directional(self, umis, counts, threshold):
        """
        Identify all UMIs within the hamming distance threshold, and where
        the counts of the first UMI is > (2 * second UMI count) - 1.
		The hamming distance between two strings (of equal length) is equal to 
		the number of positions at which the string characters are different.
        """

        adj_list = {umi: [] for umi in umis}

        if len(umis) > 25:
            umi_length = len(umis[0])
            substr_idx = build_substr_idx(umis, umi_length, threshold)
            iter_umi_pairs = iter_nearest_neighbours(umis, substr_idx)
        else:
            iter_umi_pairs = itertools.combinations(umis, 2)
        for umi1, umi2 in iter_umi_pairs:
            if edit_distance(umi1, umi2) <= threshold:
                if counts[umi1] >= (counts[umi2]*2)-1:
                    adj_list[umi1].append(umi2)
                if counts[umi2] >= (counts[umi1]*2)-1:
                    adj_list[umi2].append(umi1)

        return adj_list
    

    def _get_connected_components_adjacency(self, umis, graph, counts):
        """Find the connected UMIs within an adjacency dictionary."""

        found = set()
        components = list()

        for node in sorted(graph, key=lambda x: counts[x], reverse=True):
            if node not in found:
                # component = self.search(node, graph)
                component = breadth_first_search(node, graph)
                found.update(component)
                components.append(component)       
        return components



    def _group_directional(self, clusters, adj_list, counts):
        """Return groups for directional method."""

	#'groups' is a 2D list structure, which contains the parent umi at the first position of each row, and all child nodes
	#(in order of highest-occuring to lowest-occuring) as subsequent elements  

        observed = set()
        groups = []


        for cluster in clusters:
            if len(cluster) == 1:
                groups.append(list(cluster))
                observed.update(cluster)
            else:
                cluster = sorted(cluster, key=lambda x: counts[x],
                                 reverse=True)
                # need to remove any node which has already been observed
                temp_cluster = []
                for node in cluster:
                    if node not in observed:
                        temp_cluster.append(node)
                        observed.add(node)
                groups.append(temp_cluster)

        return groups


    def __init__(self, cluster_method="directional"):
        """Select the required class methods for the cluster_method."""
        
        self.max_umis_per_position = 0
        self.total_umis_per_position = 0
        self.positions = 0
        
        if cluster_method == "directional":
            self.get_adj_list = self._get_adj_list_directional #Represent graph data structure by creating an adjacency list 
            self.get_connected_components = self._get_connected_components_adjacency
            self.get_groups = self._group_directional #Create 2D list structure to represent graph's groups 


    def __call__(self, umis, counts, config):
        """Counts is a dictionary that maps UMIs to their counts."""

        
        threshold = int(config['SETTINGS']['umi_edit_distance_threshold']) if config else 1

        umis = list(umis)
        
        self.positions += 1
        
        number_of_umis = len(umis)
        self.total_umis_per_position += number_of_umis
        
        if number_of_umis > self.max_umis_per_position:
            self.max_umis_per_position = number_of_umis
            
        len_umis = [len(x) for x in umis]
        
        assert len_umis, "No UMIs present(!)"

        assert max(len_umis) == min(len_umis), (
            "Not all UMIs are the same length(!): {} - {}"
            .format(min(len_umis), max(len_umis)))
        
        adj_list = self.get_adj_list(umis, counts, threshold)
        clusters = self.get_connected_components(umis, adj_list, counts)
        final_umis = [tuple(x) for x in self.get_groups(clusters, adj_list, counts)]

        return final_umis
