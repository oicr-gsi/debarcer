import pysam
import src.umi_network_collapse as network

#from utils.generate_plots import umi_plot


#class UMIGroup:
#    """Contains position and count info for all UMIs representing one group."""
#    
#    def __init__(self, key):
#        self.families = {}
#        self.key = key
#    
#    def addNew(self, pos):
#        self.families[pos] = 1
#    
#    def add(self, pos):
#        self.families[pos] += 1
#    
#    @functools.lru_cache(maxsize=4, typed=False)
#    def getClosest(self, pos, pos_threshold):
#        for fampos in self.families:
#            if abs(pos - fampos) <= pos_threshold:
#                return fampos
#        return None
#    
#    def key(self):
#        return key
#    
#    
#    ######## added for viewing
#    def __str__(self):
#        return "{0}".format([(i, self.families[i]) for i in self.families])





def umi_count(contig, region_start, region_end, bam_file):
    '''
    (str, int, int, file) -> dict
    
    :param contig: Chromosome, eg chrN
    :param region_start: Start index of the region, 1-based inclusive
    :param region_end: End index of the region, 1-based inclusive
    :bam_file: Bam file with umi in read names
    
    Returns a dictionary of umi tally for each umi in a given region
    start and end are 1-based inclusive and get converted to 0-based by pysam
    '''
    
    umi_counts = {}
    
    with pysam.AlignmentFile(bam_file, "rb") as bam_reader:
        for read in bam_reader.fetch(contig, region_start, region_end):
            # extract the umi sequence from read
            # umi <- list of umi sequences
            umis = read.query_name.split(':')[-1].split(';')
            for i in umis:
                if i in umi_counts:
                    umi_counts[i] += 1
                else:
                    umi_counts[i] = 1
    return umi_counts


#def group_position(contig, region_start, region_end, bam_file, umi_groups, pos_threshold):
#    """
#    (str, int, int, file, list, int) -> dict
#    
#    :param contig: Chromosome, eg ChrN
#    :param region_start: Start index of the region, 1-based inclusive
#    :param region_end: End index of the region, 1-based inclusive
#    :bam_file: Bam file with umi in read names
#    :umi_groups: List with groups of umi sequences separated by given hamming distance
#    :pos_threshold: 
#    
#    
#    Splits umi_groups into families (umi + position pairs)."""
#
#    umi_table = {}
#
#    ## Initialize the table with references to UMIGroup objs
#    for group in umi_groups:
#        new_group = UMIGroup(key=group[0])
#        for umi in group:
#            umi_table[umi] = new_group
#
#    with pysam.AlignmentFile(bam_file, "rb") as bam_reader:
#        for read in bam_reader.fetch(contig, region_start, region_end):
#            # umi <- list of umi sequences
#            umis = read.query_name.split(':')[-1].split(';')
#            # get the start position 0-based
#            pos = read.reference_start
#            # for each umi sequence
#            for umi in umis:
#                if umi in umi_table:
#                    # get the umi group 
#                    umi_group = umi_table[umi]
#                    # form families
#                    families = umi_group.families
#                    if not families:
#                        umi_group.addNew(pos)
#                    elif pos in families:
#                        umi_group.add(pos)
#                    else:
#                        closest = umi_table[umi].getClosest(pos, pos_threshold)
#                        if closest:
#                            umi_group.add(closest)
#                        else:
#                            umi_group.addNew(pos)
#
#    return umi_table





def map_umi_to_parent(umi_groups):
    '''
    (list) -> dict
    
    :umi_groups: List with groups of umi sequences separated by given hamming distance
    
    Return a dictionary of umi sequences: parent umi sequence value, pairs
    '''
    
    # create a dict {umi : parent}
    D = {}
    for group in umi_groups:
        # parent is always the 1st umi in the group tuple
        parent = group[0]
        for umi in group:
            D[umi] = parent
    return D


def cluster_umis(umis, counts, dist_threshold):
    '''
    (list, dict, int) -> list
    
    :param umis: List of umis sequences, keys of dictionary counts
    :param counts" Dictionary of umi:count value, pairs
    :param dist_threshold: hamming distance to connect 2 umis
    
    Return a list of tuples, each containing a group of umi sequences distant 
    from each other by dist_threshold mutations
    '''

    # initialize clusterer
    clusterer = network.UMIClusterer(cluster_method="directional")  
    # cluster all umis by dist_threshold
    umi_groups = clusterer(umis, counts, dist_threshold)
    return umi_groups



def extract_umi_from_read(contig, region_start, region_end, bam_file, umi_groups):
    """
    (str, int, int, file, list) -> dict
    
    :param contig: Chromosome, eg ChrN
    :param region_start: Start index of the region, 1-based inclusive
    :param region_end: End index of the region, 1-based inclusive
    :bam_file: Bam file with umi in read names
    :umi_groups: List with groups of umi sequences separated by given hamming distance
        
    Return a dictionary of dictionaries with parent umi and all children umis
    with their count at each position
    """
    
    # set a dict {parent: {umi: {pos: count}}}
    D = {}

    # map all umi to their parent
    parent_umi = map_umi_to_parent(umi_groups)

    # loop over aligned reads
    with pysam.AlignmentFile(bam_file, "rb") as bam_reader:
        for read in bam_reader.fetch(contig, region_start, region_end):
            # umi <- list of umi sequences
            umis = read.query_name.split(':')[-1].split(';')
            # get the start position 0-based
            pos = read.reference_start
            # for each umi sequence
            for umi in umis:
                # get the parent umi
                parent = parent_umi[umi]
                # initialize inner dict if parent not in umi_families
                if parent not in D:
                    D[parent] = {}
                # check if umi is recorded for that group
                if umi not in D[parent]:
                    D[parent][umi] = {}
                # check if position is recorded
                if pos in D[parent][umi]:
                    D[parent][umi][pos] += 1
                else:
                    D[parent][umi][pos] = 1
    return D


def most_frequent(L):
    '''
    (list) -> tuple
    
    :param L: A list of umi (pos, count) sorted by position
    
    Return the (position, count) of the most abundant umi (ties broken arbitrarily)
    Pre-condition: L is sorted by the first value in tuples
    '''
    
    pos, frequent = L[0][0], L[0][1]
    for i in L:
        if i[1] >= frequent:
            pos, frequent = i[0], i[1]
    return (pos, frequent)
    

def find_group_families(umi_families, pos_threshold, ignore_others):
    '''
    (dict, int, bool) -> dict
    
    :param umi_families: A dictionary holding parent umi, their children and their count at positions where they are found
    :param pos_threshold: Window size to group indivual umis into families within groups 
    :param ignore_others: Ignore families distant from the most abundant family
    
    Return a dictionary with umi key and dictionary with parent and umi family count at each position as value     
    '''
    
    # create a dict {umi: {'parent': parent, 'positions': {pos: count}}} 
    C = {}
    
    # identify the most frequent position for given umi parent
    for parent in umi_families:
        # make a list of (pos, count) from all children umis
        L = []
        for umi in umi_families[parent]:
            for pos in umi_families[parent][umi]:
                L.append((pos, umi_families[parent][umi][pos]))
        # sort list of positions, counts by position
        L.sort()
        # create a dict {pos: count} recording umi counts within family merging
        # positions onto the position of the most frequent umi within the pos_threshold
        D = {}
        # stop when all positions have been recorded
        while len(L) != 0:
            # identify the position with the most abundant umi
            most_abundant = most_frequent(L)
            assert most_abundant[0] not in D
            D[most_abundant[0]] = most_abundant[1]
            # remove most abundant from list
            L.remove(most_abundant)
            # identify all positions within pos_threshold
            # make a list of positions to remove
            to_remove = []
            for i in L:
                if abs(i[0] - most_abundant[0]) <= pos_threshold:
                    # add count to most_abundant positon
                    D[most_abundant[0]] += i[1]
                    to_remove.append(i)
            # remove posisions from L
            for i in to_remove:
                L.remove(i)
        # remove families distant from the most abundant family if option is actioned
        if ignore_others == True:
            if len(D) > 1:
                # make a list of (position, count) sorted by positions
                L = [(pos, D[pos]) for pos in D]
                L.sort()
                keep_pos = most_frequent(L)[0]
                # make a list of positions to remove
                to_remove = [i for i in D if i != keep_pos]
                for i in to_remove:
                    del D[i]
        # add position and counts to individual umi within group   
        for umi in umi_families[parent]:
            assert umi not in C
            C[umi] = {'parent': parent, 'positions': D}
    return C


def get_umi_families(contig, region_start, region_end, bam_file, pos_threshold, dist_threshold, ignore_others):
    """
    
    (str, int, int, file, int, int) -> list
    
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 1-based inclusive
    :param region_end: End index of the region of interest. 1-based inclusive
    :param bam_file: Path to the bam file
    :param pos_threshold: Position threshold to group umis together 
    :param dist_threshold: The hamming distance threshold to connect parent and child umis     
    :param ignore_others: Ignore families distant from the most abundant family    
        
    Returns a list of tuples (umi_group, pos) pairs representing error-corrected UMI families.
    start and end are 1-based inclusive and get converted to 0-based by pysam
    """ 
    
    print("Counting UMIs...")
    # count umi sequences -> dict {umi_seq: count}
    counts = umi_count(contig, region_start, region_end, bam_file)
    umis = counts.keys()
    
    print("Clustering UMIs...")
    # cluster umis -> list of umi groups connected by dist_threshold eg. [('AAA', 'ATA', 'AAG'), ...] 
    umi_groups = cluster_umis(umis, counts, dist_threshold)
        
    print("Grouping UMIs by position...")
    # record all umis per group and position
    umi_positions = extract_umi_from_read(contig, region_start, region_end, bam_file, umi_groups)
    
    # get the positions and counts of umi families within each group (position is from the most abundant family)  
    umi_families = find_group_families(umi_positions, pos_threshold, ignore_others)
    
    return umi_families, umi_groups




#def umi_datafile(umi_groups):
#    '''
#    
#    
#    
#    
#    '''
#    
#    
#    
#
#    #Each 'row' in umi_groups is a tuple, where the first element is the parent umi and all following elements in the tuple are the child umis
#    elements_per_row = [len(row) for row in umi_groups]
#    total_parent_umi_count = len(elements_per_row)
#    total_child_umi_count = sum(elements_per_row)-total_parent_umi_count
#
#    parent_umis = [row[0] for row in umi_groups] 
#    children_per_parent = [count-1 for count in elements_per_row]
#
#    #Below dictionary has key=parent umi sequence, and value=num of child umis corresponding to that parent umi
#    umis_and_numofchildren = {}
#    for i in range(len(parent_umis)):
#        umis_and_numofchildren[(parent_umis[i])] = children_per_parent[i]
#   
#    #var is a list of (num,freq) pairs where num=number of children (0,1,2,3...) & freq=no. of parents having that number of children
#    var = [(k, len(list(v))) for k, v in itertools.groupby(sorted(umis_and_numofchildren.values()))]
#
#    num_of_children = ""
#    freq_of_parent_umis = ""
#    for i in range(len(var)):
#        if i == 0:
#            num_of_children = str(var[i][0])
#            freq_of_parent_umis = str(var[i][1])
#        else:
#            num_of_children = num_of_children+","+str(var[i][0])
#            freq_of_parent_umis = freq_of_parent_umis+","+str(var[i][1])
#
#
#    return total_parent_umi_count, total_child_umi_count, num_of_children, freq_of_parent_umis


def umi_datafile(umi_groups):
    '''
    (list) -> list
    
    :umi_groups: List with groups of umi sequences separated by given hamming distance
    
    Return a list of strings representing total number of parent umis, total of children umis,
    a comma-separated sorted list of number of children per parent, and a corresponding
    comma-separated list of counts of parents having a given number of children
    '''
    
    # count the number of children for each parent
    children = {}
    for group in umi_groups:
        parent = group[0]
        assert parent not in children
        children[parent] = len(group[1:])
    
    # count the total number of parents
    total_parent = len(children)
    # count the total number of children
    total_children = sum(children.values())
    
    # count the number of parents having a given number of children
    freq_parent = {}
    # loop over parent sequences
    for i in children:
        # get the number of children for that parent
        j = children[i]
        if j in freq_parent:
            # update parent count
            freq_parent[j] += 1
        else:
            # initialize counter
            freq_parent[j] = 1
    
    # make a sorted list of number of children per parent
    num_children = [i for i in freq_parent]
    num_children.sort()
    # make a parallel list of number of parent having x number of children
    num_parent = [freq_parent[i] for i in num_children]
    
    num_children = ','.join(list(map(lambda x: str(x), num_children)))
    num_parent = ','.join(list(map(lambda x: str(x), num_parent)))
    
    return [str(total_parent), str(total_children), num_children, num_parent]
    