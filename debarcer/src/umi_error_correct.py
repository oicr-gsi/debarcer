import pysam
import src.umi_network_collapse as network

def umi_count(contig, region_start, region_end, bam_file):
    '''
    (str, int, int, file) -> dict
    
    :param contig: Chromosome, eg chrN
    :param region_start: Start index of the region, 0-based half-opened
    :param region_end: End index of the region, 0-based half opened
    :bam_file: Bam file with umi in read names
    
    Returns a dictionary of umi tally for each umi in a given region
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



def is_overlapping(read_start, read_end, region_start, region_end):
    '''
    (int, int, int, int) -> bool
    
    :param read_start: read start position, 0-based half open
    :param read_end: read end position. Points to one past the last aligned residue
    :param region_start: Start index of the region, 0-based half opened
    :param region_end: End index of the region, 0-based half opened
    
    Return True if the read is contained within the genomic interval specified by contig,
    region_start and region_end
    '''
    
    # read_pos is None if unmapped read or missing cigar
    try:
        if region_start <= read_start < region_end and read_end <= region_end:
            overlap = True
        else:
            overlap = False
    except:
        overlap = False
    return overlap


def extract_umi_from_read(contig, region_start, region_end, bam_file, umi_groups, truncate):
    """
    (str, int, int, file, list, bool) -> dict
    
    :param contig: Chromosome, eg ChrN
    :param region_start: Start index of the region, 0-based half opened
    :param region_end: End index of the region, 0-based half opened
    :param bam_file: Bam file with umi in read names
    :param umi_groups: List with groups of umi sequences separated by given hamming distance
    :param truncate: Skip reads overlapping with the genomic interval if True    
        
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
            pos = int(read.reference_start)
            end = read.reference_end
            
            
            if end == None:
                print(pos, read, end, read.query_alignment_length, read.query_length, read.reference_length, read.cigarstring, read.cigartuples)
                raise ValueError('end is None')
            else:
                end = int(end)
            
            
            
            
            
            # skip reads overlapping with region if truncate is True
            if truncate == True and is_overlapping(pos, end, region_start, region_end) == True:
                continue
            else:
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
    

def find_group_families(contig, umi_families, pos_threshold, ignore_others):
    '''
    (str, dict, int, bool) -> dict
    
    :param contig: Chromosome, eg chrN
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
        # create a dict {'chr:pos': count} recording umi counts within family merging
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
            d = {}
            for i in D:
                d[contig + ':' + str(i)] = D[i]
            C[umi] = {'parent': parent, 'positions': d}
    return C


def get_umi_families(contig, region_start, region_end, bam_file, pos_threshold, dist_threshold, ignore_others, truncate):
    """
    
    (str, int, int, file, int, int, bool) -> tuple
    
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param bam_file: Path to the bam file
    :param pos_threshold: Position threshold to group umis together 
    :param dist_threshold: The hamming distance threshold to connect parent and child umis     
    :param ignore_others: Ignore families distant from the most abundant family
    :param truncate: Skip reads overlapping with the genomic interval if True    
        
    Returns a tuple of dictionaries with umi information before and after grouping
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
    umi_positions = extract_umi_from_read(contig, region_start, region_end, bam_file, umi_groups, truncate)
    
    # get the positions and counts of umi families within each group (position is from the most abundant family)  
    umi_families = find_group_families(contig, umi_positions, pos_threshold, ignore_others)
    
    # add contig to positions in umi_positions dict
    D = {}
    for parent in umi_positions:
        for umi in umi_positions[parent]:
            for pos in umi_positions[parent][umi]:
                if parent not in D:
                    D[parent] = {}
                if umi not in D[parent]:
                    D[parent][umi] = {}
                newpos = contig + ':' + str(pos)
                D[parent][umi][newpos] = umi_positions[parent][umi][pos]
       
    return umi_families, umi_groups, D


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
    