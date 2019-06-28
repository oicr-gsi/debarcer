import os
import sys
import pysam
import configparser
import argparse
import operator
import functools



def get_ref_seq(contig, region_start, region_end, reference):
    '''
    (str, int, int, str) -> str
    
    :param contig: Chromosome, eg chrN
    :param region_start: Start index of the region, 0-based half-opened
    :param region_end: End index of the region, 0-based half opened
    :param reference: Path to the reference genome
    
    Returns the sequence of the reference genome on contig between region_start and region_end included
    '''
    
    with pysam.FastaFile(reference) as reader:
        ref_seq = reader.fetch(contig, region_start, region_end).upper()
    return ref_seq



def find_closest(pos, L):
    '''
    (int, list) -> tuple
    
    :param pos: Position of interest along chromosome 0-based
    :param L: List of (positions, counts) 
    
    Returns a tuple (i, k) corresponding to the closest position i from pos
    and the highest count k if multiple counts exist for the smallest distance
    between pos and i
    '''
    
    # make a dict {distance: count}
    # with dist being the distance between pos and each positions
    D = {}
    for i in L:
        dist = abs(pos - i[0])
        if dist in D:
            D[dist].append(i[1])
        else:
            D[dist] = [i[1]]
    # sort all counts from smallest to highest
    for i in D:
        D[i].sort()
    # make a sorted list of distances from smallest to highest
    distances = [i for i in D]
    distances.sort()
    # get the (distance, count) for the smallest distance from pos
    # retrieve the highest count if multiple counts recorded per distance
    return (distances[0], D[distances[0]][-1])


def get_consensus_seq(umi_families, fam_size, ref_seq, contig, region_start, region_end, bam_file, pos_threshold, max_depth=1000000, truncate=True, ignore_orphans=True):
    '''
    
    (dict, int, str, str, int, int, str, int, int, bool, bool) -> (dict, dict)
    
    
    :param umi_families: Information about each umi: parent umi and positions, counts of each family within a given group
    :param fam_size: Minimum umi family size 
    :param ref_seq: Sequence of the reference corresponding to the given region
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param bam_file: Path to the bam file
    :param pos_threshold: Window size to group indivual umis into families within groups 
    :param max_depth: Maximum read depth. Default is 1000000 reads
    :param truncate: Consider only pileup columns within interval defined by region start and end. Default is True
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair). Default is True
    
    Returns a tuple with a dictionary representing consensus info at each base position in the given region
    and a dictionary to keep track of family size for each position
    '''


    # create a dict to store umi family size at each position {position: {parent: {distance: count}}}
    FamSize = {}

    # create a dict to store consensus seq info
    consensus_seq = {}
    
    with pysam.AlignmentFile(bam_file, "rb") as reader:
        # loop over pileup columns
        for pileupcolumn in reader.pileup(contig, region_start, region_end, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans):
            # get column position. by default consider only positions within region
            # however, number of reads in families consider reads overlapping with region
            # not only contained within region
            pos = int(pileupcolumn.reference_pos)  
            assert pos != region_end
            # loop over reads in pileup column
            for read in pileupcolumn.pileups:
                # get read information
                read_data = read.alignment
                read_name, start_pos = read_data.query_name, int(read_data.reference_start)
                # get all recorded umis
                umis = read_name.split(":")[-1].split(';')
                for umi in umis:
                    # check that umi is recorded
                    if umi in umi_families:
                        # find closest family from umi
                        # make a list of (positions counts)
                        L = [(int(i), umi_families[umi]['positions'][i]) for i in umi_families[umi]['positions']]
                        closest, count = find_closest(start_pos, L)

                        # check if closest family is within the position threshold
                        if closest <= pos_threshold:
                            # found a umi family. check if family count is greater than family threshold
                            if count >= fam_size:
                                # get the parent sequence                                
                                parent = umi_families[umi]['parent']
                                                                
                                # keep track of family size used to derive alt
                                # for a given position , parent and read position
                                if pos not in FamSize:
                                    FamSize[pos] = {}
                                if parent not in FamSize[pos]:
                                    FamSize[pos][parent] = {}
                                FamSize[pos][parent][closest] = count
                                
                                # use family key to count allele. collapsing is done within families. not per position
                                family_key = parent + str(closest)
                                
                                ref_pos = pos - region_start
                    
                                # read.indel is indel length of next position 
                                # 0 --> not indel; > 0 --> insertion; < 0 --> deletion
                                                
                                # get reference and alternative bases  
                                if not read.is_del and read.indel == 0:
                                    ref_base = ref_seq[ref_pos]
                                    alt_base = read_data.query_sequence[read.query_position]
                                elif read.indel > 0:
                                    # Next position is an insert (current base is ref)
                                    ref_base = ref_seq[ref_pos]
                                    alt_base = read_data.query_sequence[read.query_position:read.query_position + abs(read.indel)+1]
                                elif read.indel < 0:
                                    # Next position is a deletion (current base + next bases are ref)
                                    ref_base = ref_seq[ref_pos:ref_pos + abs(read.indel) + 1]
                                    alt_base = read_data.query_sequence[read.query_position]
                            
                                # query position is None if is_del or is_refskip is set
                                if not read.is_del and not read.is_refskip:
                                    # add base info
                                    allele = (ref_base, alt_base)
                                    # count the number of reads supporting this allele
                                    if pos not in consensus_seq:
                                        consensus_seq[pos] = {}
                                    if family_key not in consensus_seq[pos]:
                                        consensus_seq[pos][family_key] = {}
                                    if allele in consensus_seq[pos][family_key]:
                                        consensus_seq[pos][family_key][allele] += 1
                                    else:
                                        consensus_seq[pos][family_key][allele] = 1
    return consensus_seq, FamSize


def get_uncollapsed_seq(ref_seq, contig, region_start, region_end, bam_file, max_depth=1000000, truncate=True, ignore_orphans=True):
    '''
    (str, str, int, int, str, str, int, bool, bool) -> dict
    
    :param ref_seq: Sequence of the reference corresponding to the given region
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param bam_file: Path to the bam file
    :param max_depth: Maximum read depth. Default is 1000000 reads
    :param truncate: Consider only pileup columns within interval defined by region start and end. Default is True
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair). Default is True
    
    Returns a nested dictionary representing counts of each base at each base position.
    '''

    uncollapsed_seq = {}

    with pysam.AlignmentFile(bam_file, "rb") as reader:
        # loop over pileup columns 
        for pileupcolumn in reader.pileup(contig, region_start, region_end, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans):
            # get column position. by default consider only positions within region
            # however, number of reads in families consider reads overlapping with region
            # not only contained within region
            pos = int(pileupcolumn.reference_pos) 
            assert pos != region_end  
            # loop over reads in pileup column
            for read in pileupcolumn.pileups:
                # read.indel is indel length of next position 
                # 0 --> not indel; > 0 --> insertion; < 0 --> deletion
                
                # get reference and alternative bases
                if not read.is_del and read.indel == 0:
                    ref_base = ref_seq[pos - region_start]
                    alt_base = read.alignment.query_sequence[read.query_position]
                elif read.indel > 0:
                    # Next position is an insert (current base is ref)
                    ref_base = ref_seq[pos - region_start]
                    alt_base = read.alignment.query_sequence[read.query_position:read.query_position + abs(read.indel) + 1]
                elif read.indel < 0:
                    # Next position is a deletion (current base + next bases are ref)
                    ref_base = ref_seq[read.query_position:read.query_position + abs(read.indel) + 1]
                    alt_base = read.alignment.query_sequence[read.query_position]
                
                # query position is None if is_del or is_refskip is set
                if not read.is_del and not read.is_refskip:
                    # add base info
                    allele = (ref_base, alt_base)
                    # count the number of reads supporting this allele
                    if pos not in uncollapsed_seq:
                        uncollapsed_seq[pos] = {}
                    if allele not in uncollapsed_seq[pos]:
                        uncollapsed_seq[pos][allele] = 1
                    else:
                        uncollapsed_seq[pos][allele] += 1
    return uncollapsed_seq


#class ConsDataRow:
#    """Holds consensus data for one base position."""
#    ref_info  = {}
#    cons_info = {}
#    stats     = {}
#    
#    def __init__(self, ref_info, cons_info, stats):
#        self.ref_info = ref_info
#        self.cons_info = cons_info
#        self.stats = stats
#
#    def impute_allele_depths(self):
#        """Returns allele depths (dict of str tuples -> ints)."""
#        return self.cons_info
#    
#    @functools.lru_cache(maxsize=2, typed=False)
#    def impute_allele_freqs(self, threshold):
#        """Returns allele frequencies (dict of str tuples -> floats)."""
#        
#        is_ref = lambda allele: allele[0] is allele[1]
#        
#        freqs = {}
#        for allele in self.cons_info:
#            freq = (self.cons_info[allele] / sum(self.cons_info.values())) * 100
#            
#            if not is_ref(allele) and freq > threshold:
#                freqs[allele] = freq
#        
#        return freqs
#    
#    def get_alleles(self, threshold):
#        """Returns alt alleles with their associated refs (list of str tuples)."""
#        
#        freqs   = self.impute_allele_freqs(threshold)
#        alleles = []
#        
#        for allele in self.cons_info:
#            if allele in freqs:
#                alleles.append(allele)
#        
#        return alleles
#
#    def get_ref_info(self):
#        return self.ref_info
#    
#    def get_cons_info(self):
#        return self.cons_info
#
#    def get_stats(self):
#        return self.stats


def get_fam_size(FamSize, position):
    '''
    (dict, int) -> tuple
    
    :param FamSize: A dictionary with family size for each position, and parent umi
    :param position: A given position in a genomic region. 0-based
    
    Returns a tuple with the minimum and mean family size at position
    '''
    
    # make a list of family sizes
    L = []
    # loop over parent at given position
    for i in FamSize[position]:
         # loop over distance between read and umi family        
         for j in FamSize[position][i]:
             L.append(FamSize[position][i][j])
    min_fam, mean_fam = min(L), sum(L) / len(L)
    return (min_fam, mean_fam)
    
    
def generate_consensus(umi_families, fam_size, ref_seq, contig, region_start, region_end, bam_file, pos_threshold, percent_threshold, count_threshold, max_depth=1000000, truncate=True, ignore_orphans=True):
    '''
    (dict, int, str, str, int, int, str, int, int, int, int, bool, bool) -> dict
    
    :param umi_families: Information about each umi: parent umi and positions,
                         counts of each family within a given group
                         positions are 0-based half opened
    :param fam_size: Minimum umi family size
    :param ref_seq: Sequence of the reference corresponding to the given region
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param bam_file: Path to the bam file
    :param pos_threshold: Window size to group indivual umis into families within groups
    :param percent_threshold: Percent consensus threshold for alternative base within family 
    :param count_threshold: Count consensus threshold for alternative base within family
    :param max_depth: Maximum read depth. Default is 1000000 reads
    :param truncate: Consider only pileup columns within interval defined by region start and end. Default is True
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair). Default is True
    
    Generates consensus data for a given family size and genomic region
    '''

    # get consensus info for each base position and umi group in the given region {pos: {fam_key: {(ref, alt):count}}}
    # get family size at each position 
    consensus_seq, FamSize = get_consensus_seq(umi_families, fam_size, ref_seq, contig, region_start, region_end, bam_file, pos_threshold, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans)

    # create a dict to store consensus info
    cons_data = {}

    # loop over positions in region
    for base_pos in range(region_start, region_end):
        # extract ref base
        ref_base = ref_seq[base_pos-region_start]
        # check if base pos has been recorded 
        if base_pos in consensus_seq:
            # record raw depth and consensus info at position
            consensuses = {}
            raw_depth = 0
            
            # compute minimum and mean family size
            min_fam, mean_fam = get_fam_size(FamSize, base_pos) 
                        
            for family in consensus_seq[base_pos]:
                # get the allele with highest count       
                cons_allele = max(consensus_seq[base_pos][family].items(), key = operator.itemgetter(1))[0]
                # compute allele frequency within umi family
                cons_denom = sum(consensus_seq[base_pos][family].values())
                cons_percent = (consensus_seq[base_pos][family][cons_allele]/cons_denom) * 100
                # compute raw depth                
                raw_depth += cons_denom
                # check if allele frequencyand allele count > thresholds
                if cons_percent >= percent_threshold and consensus_seq[base_pos][family][cons_allele] >= count_threshold:
                    # count allele
                    if cons_allele in consensuses:
                        consensuses[cons_allele] += 1
                    else:
                        consensuses[cons_allele] = 1
            # compute consensus depth across all alleles
            cons_depth = sum(consensuses.values())
            
            # compute ref frequency
            if (ref_base, ref_base) in consensuses:
                ref_freq = (consensuses[(ref_base, ref_base)] / cons_depth) * 100
            else:
                ref_freq = 0
            
            # record ref, consensus and stats info
            ref_info = {"contig": contig, "base_pos": base_pos, "ref_base": ref_base}
            cons_info = consensuses
            stats = {"rawdp": raw_depth, "consdp": cons_depth, "min_fam": min_fam, "mean_fam": mean_fam, "ref_freq": ref_freq}
                    
            cons_data[base_pos] = {'ref_info': ref_info, 'cons_info': cons_info, 'stats': stats}
    return cons_data



def generate_uncollapsed(ref_seq, contig, region_start, region_end, bam_file, max_depth=1000000, truncate=True, ignore_orphans=True):
    '''
    (str, str, int, int, str, int, bool, bool) -> dict
    
    :param ref_seq: Sequence of the reference corresponding to the given region
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param bam_file: Path to the bam file
    :param max_depth: Maximum read depth. Default is 1000000 reads
    :param truncate: Consider only pileup columns within interval defined by region start and end. Default is True
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair). Default is True
    
    Generates uncollapsed consensus data for the genomic region
    '''
    
    # get uncolapased seq info {pos: {(ref, atl): count}}
    uncollapsed_seq = get_uncollapsed_seq(ref_seq, contig, region_start, region_end, bam_file, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans)
    
    # create a dict to store consensus info
    cons_data = {}
    
    # loop over positions in genomic region
    for base_pos in range(region_start, region_end):
        # extract ref base    
        ref_base = ref_seq[base_pos-region_start]
        # check if base pos has been recorded     
        if base_pos in uncollapsed_seq:
            # compute depth at position        
            depth = sum(uncollapsed_seq[base_pos].values())
            # compute ref frequency
            if (ref_base, ref_base) in uncollapsed_seq[base_pos]:
                ref_freq = (uncollapsed_seq[base_pos][(ref_base, ref_base)] / depth) * 100
            else:
                ref_freq = 0

            # record ref, consensus and stats info 
            ref_info = {"contig": contig, "base_pos": base_pos, "ref_base": ref_base}
            cons_info = uncollapsed_seq[base_pos]
            stats = {"rawdp": depth, "consdp": depth, "min_fam": 0, "mean_fam": 0, "ref_freq": ref_freq}
            
            cons_data[base_pos] = {'ref_info': ref_info, 'cons_info': cons_info, 'stats': stats}
    return cons_data


def raw_table_output(cons_data, ref_seq, contig, region_start, region_end, output_path, ref_threshold, all_threshold):
    '''
    (dict, str, str, int, int, str, num, num) -> None
    
    :param ref_seq: Sequence of the reference corresponding to the given region
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param output_path: Output directory
    
    :param ref_threshold: ??
        
    :param all_threshold: ??
    
    Writes a long-form consensus file for every event detected in the collapsed data
    '''
    
    # get the path to the output file
    OutputFile = os.path.join(output_path, '{}:{}-{}.cons'.format(contig, region_start, region_end))
    newfile = open(OutputFile, 'w')

    Header = ['CHROM', 'POS', 'REF', 'A', 'C', 'G', 'T', 'I', 'D', 'N', 'RAWDP', 'CONSDP', 'FAM', 'REF_FREQ', 'MEAN_FAM']
    newfile.write('\t'.join(Header) + '\n')

    # loop over positions         
    for base_pos in range(region_start, region_end):
        # get the reference
        ref_base = ref_seq[base_pos-region_start]
        # loop over fam size
        for f_size in cons_data:
            # check if position recorded for given fam size
            if base_pos in cons_data[f_size]:
                # get ref, stats and consensus info
                ref = cons_data[f_size][base_pos]['ref_info']
                
                assert ref_base == ref['ref_base']
                                
                cons = cons_data[f_size][base_pos]['cons_info']
                stats = cons_data[f_size][base_pos]['stats']
                # count each bases
                counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'I': 0, 'D': 0, 'N': 0}
                for allele in cons:
                    # ref > 1 => deletion
                    if len(allele[0]) > 1:
                        counts['D'] += cons[allele]
                    # allele > 1 => insertion
                    elif len(allele[1]) > 1:
                        counts['I'] += cons[allele]
                    else:
                        counts[allele[1]] += cons[allele]
                # write line to file
                line = [contig, base_pos + 1, ref_base, counts['A'], counts['C'],
                        counts['G'], counts['T'], counts['I'], counts['D'],
                        counts['N'], stats['rawdp'], stats['consdp'], f_size,
                        stats['ref_freq'], stats['mean_fam']]
                newfile.write('\t'.join(list(map(lambda x: str(x), line))) + '\n')
                
    
    # close file after writing
    newfile.close()
            

#                row = cons_data[f_size][base_pos]
#                ref = row.get_ref_info()
#                cons = row.get_cons_info()
#                stats = row.get_stats()
#
#                if stats['ref_freq'] <= ref_threshold:
#
#                    alleles = row.get_alleles(all_threshold)
#                    ref_bases = set([allele[0] for allele in alleles])
#                    ref_allele = (ref_seq[base_pos - region_start], ref_seq[base_pos - region_start])
#                    depths = row.impute_allele_depths()
#                    ref_depth = depths[ref_allele] if ref_allele in depths else 0
#                    alt_freqs = row.impute_allele_freqs(all_threshold)
#
#                    info = "RDP={};CDP={};MIF={};MNF={:.1f}".format(
#                        stats['rawdp'], stats['consdp'], stats['min_fam'], stats['mean_fam'])
#                    fmt_string = "AD:AL:AF" # Allele depth, alt allele depth, reference frequency
#
#                    for ref_base in ref_bases:
#
#                        #Handle error where ref_base is assigned to multiple ref_alleles
#                        if len(ref_base) > 1:
#                            ref_base = ref_base[0]
#                        if (base_pos < 170837514) and (base_pos > 170837510):
#                            print("--BASE_POS: "+str(base_pos)+" ref_base: "+ref_base+"--")
#
#                        snips = []
#                        for allele in alleles:
#                            if allele[0] == ref_base:
#                                snips.append(allele)
#
#                        alt_string = ','.join( [allele[1] for allele in snips] )
#                        depth_string = ','.join( [str(depths[allele]) for allele in snips] )
#                        freq_string = ','.join( ["{:.2f}".format(alt_freqs[allele]) for allele in snips] )
#                        smp_string = "{}:{}:{}".format(ref_depth, depth_string, freq_string)
#                        filt = "PASS" if any( [depths[alt] > 10 for alt in snips] ) else "a10"
#                        #writer.write("# {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(contig, base_pos, ".", ref_base, alt_string, "0", filt, info, fmt_string, smp_string))


def generate_consensus_output(reference, contig, region_start, region_end, bam_file, umi_families, output_path, fam_size, pos_threshold, percent_threshold, count_threshold, ref_threshold, all_threshold, max_depth=1000000, truncate=True, ignore_orphans=True):
    """(Main) generates consensus output file."""

    # get minimum umi family sizes
    family_sizes = list(map(lambda x: int(x.strip()), fam_size.split(',')))

    ## Get reference sequence for the region 
    print("Getting reference sequence...")
    ref_seq = get_ref_seq(contig, region_start, region_end, reference)

    ## Get consensus data for each f_size + uncollapsed data
    print("Building consensus data...")
    cons_data = {}
    
    for f_size in family_sizes:
        # check if 0 is passed as fam_size argument
        if f_size == 0:
            # compute consensus for uncollapsed data
            cons_data[f_size] = generate_uncollapsed(ref_seq, contig, region_start, region_end, bam_file, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans)
        else:
            cons_data[f_size] = generate_consensus(umi_families, f_size, ref_seq, contig, region_start, region_end, bam_file, pos_threshold, percent_threshold, count_threshold, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans)
    
    # compute consensus for uncollapsed data if not in fam_size argument
    if 0 not in family_sizes:
        cons_data[0] = generate_uncollapsed(ref_seq, contig, region_start, region_end, bam_file, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans)

    ## Output
    print("Writing output...")
    raw_table_output(cons_data, ref_seq, contig, region_start, region_end, output_path, ref_threshold, all_threshold)
    

if __name__=="__main__":
    
    ## Argument + config parsing and error handling
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam_file', help='Path to your BAM file.')
    parser.add_argument('-r', '--region', help='Region to analyze (string of the form chrX:posA-posB).')
    parser.add_argument('-o', '--output_path', help='Path to write output files to.')
    parser.add_argument('-c', '--config', help='Path to your config file.')
    parser.add_argument('-t', '--tally', help='Path to your tally (output of UMI_count.py).')

    args = parser.parse_args()

    if args.config:
        config = configparser.ConfigParser()
        config.read(config_file)
    else:
        config = None
    
    region = args.region
    if any(x not in region for x in ["chr", ":", "-"]):
        raise ValueError('Incorrect region string (should look like chr1:1200000-1250000).')
        sys.exit(1)

    contig = region.split(":")[0]
    region_start = int(region.split(":")[1].split("-")[0])
    region_end = int(region.split(":")[1].split("-")[1])

    bam_file = handle_arg(args.bam_file, config['PATHS']['bam_file'] if config else None, 
                    'No BAM file provided in args or config.')
    output_path = handle_arg(args.output_path, config['PATHS']['output_path'] if config else None, 
                    'No output path provided in args or config.')
    tally_file = handle_arg(args.tally, output_path + '/' + region + '.tally' if config else None, 
                    'No tally file provided.')

    ## Output
    generate_consensus_output(contig, region_start, region_end, bam_file, tally_file, output_path, config)
