import os
import pysam
import operator
import yaml
from src.utilities import get_umi_from_name

def find_closest(pos, L):
    '''
    (int, list) -> tuple
    
    :param pos: Position of interest along chromosome 0-based
    :param L: List of (positions, counts) 
    
    Returns a tuple (i, k) corresponding to the closest position i from pos
    and the highest count k if multiple counts exist for the smallest distance
    between pos and i
    '''
    
    # make a dict {distance: [(count, position)]}
    # with dist being the distance between pos and each positions
    D = {}
    for i in L:
        dist = abs(pos - i[0])
        if dist in D:
            D[dist].append((i[1], i[0]))
        else:
            D[dist] = [(i[1], i[0])]
    # sort all counts from smallest to highest
    for i in D:
        # sort on count
        D[i].sort(key=lambda x: x[0])
    # make a sorted list of distances from smallest to highest
    distances = sorted(D.keys())
    # get the (distance, count, position) for the smallest distance from pos
    # retrieve the highest count if multiple counts recorded per distance
    smallest_dist = distances[0]
    return (smallest_dist, D[smallest_dist][-1][0], D[smallest_dist][-1][1])


def get_consensus_seq(umi_families, fam_size, contig, region_start, region_end, bam_file, pos_threshold, max_depth, truncate, ignore_orphans, stepper):
    '''
    
    (dict, int, str, int, int, str, int, int, bool, bool, str) -> (dict, dict)
    
    
    :param umi_families: Information about each umi: parent umi and positions, counts of each family within a given group
    :param fam_size: Minimum umi family size 
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param bam_file: Path to the bam file
    :param pos_threshold: Window size to group indivual umis into families within groups 
    :param max_depth: Maximum read depth
    :param truncate: Consider only pileup columns within interval defined by region start and end if True
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair) if True
    :param stepper: Controls how the iterator advances. Accepeted values:
                    'all': skip reads with following flags: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
                    'nofilter': uses every single read turning off any filtering
        
    Returns a tuple with a dictionary representing consensus info at each base position in the given region
    and a dictionary to keep track of family size for each position
    '''

    # create a dict to store umi family size at each position {position: {parent: {distance: count}}}
    FamSize = {}

    # create a dict to store consensus seq info
    # {pos: {'ref_base': ref_base, 'families': {famkey: {allele: count}}}}
    consensus_seq = {}
    
    with pysam.AlignmentFile(bam_file, "rb") as reader:
        # loop over pileup columns
        for pileupcolumn in reader.pileup(contig, region_start, region_end, truncate=truncate, ignore_orphans=ignore_orphans, max_depth=max_depth, stepper=stepper):
            # get column position
            pos = int(pileupcolumn.reference_pos)  

            # restict pileup columns to genomic region
            if region_start <= pos < region_end:
                # loop over reads in pileup column
                for read in pileupcolumn.pileups:
                    # pileupread obj: represention of an aligned read
                    # get read information as AlignedSegment
                    read_data = read.alignment
                    
                    # skip unmapped, secondary and supplementary reads/alignments
                    if read_data.is_unmapped == False and read_data.is_secondary == False and read_data.is_supplementary == False:
                        read_name, start_pos = read_data.query_name, int(read_data.reference_start)
                        # extract umi. expecting a single umi
                        umi = get_umi_from_name(read_name)
                        # check that umi is recorded
                        if umi in umi_families:
                            # find closest family from umi
                            # make a list of (positions counts)
                            L = [(int(i.split(':')[1]), umi_families[umi]['positions'][i]) for i in umi_families[umi]['positions']]
                            closest, count, position_closest = find_closest(start_pos, L)
                        
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
                                    family_key = parent + str(position_closest)
                                
                                    # skip positions with deletions or ref not defined
                                    # events are captured at the position before they occur
                                    if not read.is_del and not read.is_refskip:
                                        # read.indel looks ahead to see if indel at next position(s) 
                                        # 0 --> not indel; > 0 --> insertion; < 0 --> deletion
    
                                        # get aligned read, ref pos and ref base 
                                        pairs = read_data.get_aligned_pairs(with_seq=True)
                                        # read.indel looks ahead to see if indel at next position(s)  
                                        if read.indel == 0:
                                            # no indel, record ref and alt 
                                            # get index of pileupcolumn pos in aligned pairs
                                            j = [i[1] for i in pairs].index(pos)
                                            # record base on ref and read
                                            ref_base = pairs[j][-1].upper()
                                            alt_base = read_data.query_sequence[read.query_position].upper()
                                        elif read.indel > 0:
                                            # next position is an insertion
                                            # get index of pileupcolumn pos in aligned pairs
                                            j = [i[1] for i in pairs].index(pos)
                                            # record base on ref and insertion on read
                                            ref_base = pairs[j][-1].upper()
                                            alt_base = read_data.query_sequence[read.query_position:read.query_position + abs(read.indel) + 1].upper()
                                        elif read.indel < 0:
                                            # next position is deletion
                                            # get index of pileupcolumn pos in aligned pairs
                                            j = [i[1] for i in pairs].index(pos)
                                            # record base on ref at pos + ref bases deleted on read and base on read
                                            ref_base = ''.join([i[-1] for i in pairs[j: j +  abs(read.indel) + 1]]).upper()
                                            alt_base = read_data.query_sequence[read.query_position]
                                                                                
                                        # add base info
                                        allele = (ref_base, alt_base)
                                        # keep track of ref_base at pos
                                        if pos not in consensus_seq:
                                            consensus_seq[pos] = {}
                                        if 'ref_base' not in consensus_seq[pos]:
                                            consensus_seq[pos]['ref_base'] = ref_base
                                        # count the number of reads supporting this allele
                                        if 'families' not in consensus_seq[pos]:
                                            consensus_seq[pos]['families'] = {}
                                        if family_key not in consensus_seq[pos]['families']:
                                            consensus_seq[pos]['families'][family_key] = {}
                                        if allele in consensus_seq[pos]['families'][family_key]:
                                            consensus_seq[pos]['families'][family_key][allele] += 1
                                        else:
                                            consensus_seq[pos]['families'][family_key][allele] = 1
    return consensus_seq, FamSize


def get_uncollapsed_seq(contig, region_start, region_end, bam_file, max_depth, truncate, ignore_orphans, stepper):
    '''
    (str, int, int, str, str, int, bool, bool, str) -> (dict, float)
    
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param bam_file: Path to the bam file
    :param max_depth: Maximum read depth
    :param truncate: Consider only pileup columns within interval defined by region start and end if True
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair) if True
    :param stepper: Controls how the iterator advances. Accepeted values:
                    'all': skip reads with following flags: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
                    'nofilter': uses every single read turning off any filtering
    
    Returns a tuple with a nested dictionary representing counts of each base at each base position,
    and the average read depth for the given region
    '''

    # create a dict {pos: {'ref_base': reference base}, {'alleles': {allele: count}}}
    uncollapsed_seq = {}

    # make a list to store number of reads 
    covArray = []
    
    with pysam.AlignmentFile(bam_file, "rb") as reader:
        # loop over pileup columns 
        for pileupcolumn in reader.pileup(contig, region_start, region_end, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans, stepper=stepper):
            # get column position
            pos = int(pileupcolumn.reference_pos) 

            # record number of non-filtered reads in pileup
            read_count = 0

            # restict pileup columns to genomic region
            if region_start <= pos < region_end:
                # loop over reads in pileup column
                for read in pileupcolumn.pileups:
                    # pileupread obj: represention of an aligned read
                    # read.indel is indel length of next position 
                    # 0 --> not indel; > 0 --> insertion; < 0 --> deletion
                    
                    # get reference and alternative bases
                    
                    # skip unmapped, secondary and supplementary reads/alignments
                    # get the AlignedSegment obj
                    read_data = read.alignment
                    if read_data.is_unmapped == False and read_data.is_secondary == False and read_data.is_supplementary == False:
                        # update read counter
                        read_count += 1
                        
                        # skip positions with deletions or ref not defined
                        # events are captured at the position before they occur
                        if not read.is_del and not read.is_refskip:
                            # get aligned read, ref pos and ref base 
                            pairs = read_data.get_aligned_pairs(with_seq=True)
    
                            # read.indel looks ahead to see if indel at next position(s)
                            if read.indel == 0:
                                # no indel, record ref and alt 
                                # get index of pileupcolumn pos in aligned pairs
                                j = [i[1] for i in pairs].index(pos)
                                # record base on ref and read
                                ref_base = pairs[j][-1].upper()
                                alt_base = read_data.query_sequence[read.query_position].upper()
                            elif read.indel > 0:
                                # next position is an insertion
                                # get index of pileupcolumn pos in aligned pairs
                                j = [i[1] for i in pairs].index(pos)
                                # record base on ref and insertion on read
                                ref_base = pairs[j][-1].upper()
                                alt_base = read_data.query_sequence[read.query_position:read.query_position + abs(read.indel) + 1].upper()
                            elif read.indel < 0:
                                # next position is deletion
                                # get index of pileupcolumn pos in aligned pairs
                                j = [i[1] for i in pairs].index(pos)
                                # record base on ref at pos + ref bases deleted on read and base on read
                                ref_base = ''.join([i[-1] for i in pairs[j: j +  abs(read.indel) + 1]]).upper()
                                alt_base = read.alignment.query_sequence[read.query_position]
                
                            # add base info
                            allele = (ref_base, alt_base)
                            # record ref base
                            if pos not in uncollapsed_seq:
                                uncollapsed_seq[pos] = {}
                            if 'ref_base' not in uncollapsed_seq[pos]:
                                uncollapsed_seq[pos]['ref_base'] = ref_base
                            # count the number of reads supporting this allele
                            if 'alleles' not in uncollapsed_seq[pos]:
                                uncollapsed_seq[pos]['alleles'] = {}
                            if allele not in uncollapsed_seq[pos]['alleles']:
                                uncollapsed_seq[pos]['alleles'][allele] = 1
                            else:
                                uncollapsed_seq[pos]['alleles'][allele] += 1
                covArray.append(read_count)                          
    # compute coverage
    try:
        coverage = sum(covArray) / len(covArray)
    except:
        coverage = 0
    
    return uncollapsed_seq, round(coverage, 2)


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
    
    
def generate_consensus(umi_families, fam_size, contig, region_start, region_end, bam_file, pos_threshold, consensus_threshold, count_threshold, max_depth, truncate, ignore_orphans, stepper):
    '''
    (dict, int, str, int, int, str, int, float, int, int, bool, bool, str) -> dict
        
    :param umi_families: Information about each umi: parent umi and positions,
                         counts of each family within a given group
                         positions are 0-based half opened
    :param fam_size: Minimum umi family size
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param bam_file: Path to the bam file
    :param pos_threshold: Window size to group indivual umis into families within groups
    :param consensus_threshold: Majority rule consensus threshold for alternative base within family 
    :param count_threshold: Count consensus threshold for alternative base within family
    :param max_depth: Maximum read depth
    :param truncate: Consider only pileup columns within interval defined by region start and end if True
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair) if True
    :param stepper: Controls how the iterator advances. Accepeted values:
                    'all': skip reads with following flags: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
                    'nofilter': uses every single read turning off any filtering    

    Generates consensus data for a given family size and genomic region
    '''

    # get consensus info for each base position and umi group in the given region
    # {pos: {'ref_base': ref_base, 'families': {famkey: {allele: count}}}}
    # get family size at each position 
    consensus_seq, FamSize = get_consensus_seq(umi_families, fam_size, contig, region_start, region_end, bam_file, pos_threshold, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans, stepper=stepper)

    # create a dict to store consensus info
    cons_data = {}

    # loop over positions in region. positions already recorded in consensus_seq
    for pos in consensus_seq:
        # extract ref base
        ref_base = consensus_seq[pos]['ref_base']
        # record raw depth and consensus info at position
        consensuses = {}
        raw_depth = 0
            
        # compute minimum and mean family size
        min_fam, mean_fam = get_fam_size(FamSize, pos) 
                        
        for family in consensus_seq[pos]['families']:
            # get the allele with highest count       
            cons_allele = max(consensus_seq[pos]['families'][family].items(), key = operator.itemgetter(1))[0]
            # compute allele frequency within umi family
            cons_denom = sum(consensus_seq[pos]['families'][family].values())
            cons_freq = (consensus_seq[pos]['families'][family][cons_allele]/cons_denom) * 100
            # compute raw depth                
            raw_depth += cons_denom
            # check if allele frequencyand allele count > thresholds
            if cons_freq >= consensus_threshold and consensus_seq[pos]['families'][family][cons_allele] >= count_threshold:
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
        ref_info = {"contig": contig, "base_pos": pos, "ref_base": ref_base}
        cons_info = consensuses
        stats = {"rawdp": raw_depth, "consdp": cons_depth, "min_fam": min_fam, "mean_fam": mean_fam, "ref_freq": ref_freq}
                  
        cons_data[pos] = {'ref_info': ref_info, 'cons_info': cons_info, 'stats': stats}
    
    return cons_data



def generate_uncollapsed(contig, region_start, region_end, bam_file, max_depth, truncate, ignore_orphans, stepper):
    '''
    (str, int, int, str, int, bool, bool, str) -> (dict, float)
    
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param bam_file: Path to the bam file
    :param max_depth: Maximum read depth
    :param truncate: Consider only pileup columns within interval defined by region start and end if True
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair) if True
    :param stepper: Controls how the iterator advances. Accepeted values:
                    'all': skip reads with following flags: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
                    'nofilter': uses every single read turning off any filtering
        
    Returns a dictionary with uncollapsed consensus data for the genomic region,
    and the average read depth per position for the region
    '''
    
    # get uncolapased seq info {pos: {'ref_base': reference base}, {'alleles': {allele: count}}}
    uncollapsed_seq, coverage = get_uncollapsed_seq(contig, region_start, region_end, bam_file, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans, stepper=stepper)
    
    # create a dict to store consensus info
    cons_data = {}
    
    # loop over positions in region. positions already recorded in uncollapsed_seq
    for pos in uncollapsed_seq:
        # get ref base    
        ref_base = uncollapsed_seq[pos]['ref_base']
        # compute depth at position        
        depth = sum(uncollapsed_seq[pos]['alleles'].values())
        # compute ref frequency
        if (ref_base, ref_base) in uncollapsed_seq[pos]['alleles']:
            ref_freq = (uncollapsed_seq[pos]['alleles'][(ref_base, ref_base)] / depth) * 100
        else:
            ref_freq = 0
        # record ref, consensus and stats info 
        ref_info = {"contig": contig, "base_pos": pos, "ref_base": ref_base}
        cons_info = uncollapsed_seq[pos]['alleles']
        stats = {"rawdp": depth, "consdp": depth, "min_fam": 0, "mean_fam": 0, "ref_freq": ref_freq}
        
        cons_data[pos] = {'ref_info': ref_info, 'cons_info': cons_info, 'stats': stats}
    
    return cons_data, coverage


def raw_table_output(cons_data, contig, region_start, region_end, outdir):
    '''
    (dict, str, int, int, str) -> None
    
    :param cons_data: Consensus data for each family size and position
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param outdir: Output directory
    
    Writes a long-form consensus file for every event detected in the collapsed data
    '''
    
    # get the path to the output file
    OutputFile = os.path.join(outdir, '{0}:{1}-{2}.cons'.format(contig, region_start + 1, region_end))
    newfile = open(OutputFile, 'w')

    Header = ['CHROM', 'POS', 'REF', 'A', 'C', 'G', 'T', 'N', 'I_(ref,ins)', 'I_counts', 'D_(ref,del)', 'D_counts', 'RAWDP', 'CONSDP', 'FAM', 'REF_FREQ', 'MEAN_FAM']
    newfile.write('\t'.join(Header) + '\n')

    # make a sorted list of positions
    positions = []
    for i in cons_data:
        positions.extend(list(cons_data[i].keys()))
    positions = sorted(list(map(lambda x: int(x), list(set(positions)))))

    # loop over positions
    for pos in positions:
        # loop over fam size
        for f_size in sorted(cons_data.keys()):
            # check that position is recorded
            if pos in cons_data[f_size]:
                # get ref, stats and consensus info
                ref = cons_data[f_size][pos]['ref_info']
                # get the reference
                ref_base = ref['ref_base']
                               
                cons = cons_data[f_size][pos]['cons_info']
                stats = cons_data[f_size][pos]['stats']
            
                # count each allele, initiate with single nucleotides
                counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
                for allele in cons:
                    # count single nucleotides
                    if len(allele[0]) == 1 and len(allele[1]) == 1:
                        # snv or no change
                        counts[allele[1]] += cons[allele]
                    else:
                        # indel, record allele and its count
                        if allele in counts:
                            counts[allele] += 1
                        else:
                            counts[allele] = 1
                # make lists of indels and indel counts 
                D, I = [], []
                for allele in counts:
                    if len(allele) == 2:
                        # indel
                        if len(allele[1]) > 1:
                            # insertions
                            I.append([counts[allele], allele])
                        elif len(allele[0]) > 1:
                            # deletions
                            D.append([counts[allele], allele])
                D.sort()
                I.sort()
                        
                line = [contig, pos + 1, ref_base, counts['A'], counts['C'],
                        counts['G'], counts['T'], counts['N'],
                        ';'.join([str(i[1]) for i in I]), ';'.join([str(i[0]) for i in I]),
                        ';'.join([str(i[1]) for i in D]), ';'.join([str(i[0]) for i in D]),
                        stats['rawdp'], stats['consdp'], f_size,
                        stats['ref_freq'], stats['mean_fam']]
                newfile.write('\t'.join(list(map(lambda x: str(x), line))) + '\n')
    # close file after writing
    newfile.close()
            


def generate_consensus_output(contig, region_start, region_end, bam_file, umi_families, outdir, fam_size, pos_threshold, consensus_threshold, count_threshold, max_depth, truncate, ignore_orphans, stepper):
    '''
    (str, str, int, int, str, dict, str, str, int, num, int, num, num, int, bool, bool, str) -> None
    
    
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param bam_file: Path to the bam file
    :param umi_families: Information about each umi: parent umi and positions,
                         counts of each family within a given group
                         positions are 0-based half opened
    :param outdir: Output directory
    :param fam_size: A comma-separated list of minimum umi family size
    :param pos_threshold: Window size to group indivual umis into families within groups
    :param consensus_threshold: Majority rule consensus threshold for alternative base within family 
    :param count_threshold: Count consensus threshold for alternative base within family
    :param max_depth: Maximum read depth
    :param truncate: Consider only pileup columns within interval defined by region start and end if True
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair) if True
    :param stepper: Controls how the iterator advances. Accepeted values:
                    'all': skip reads with following flags: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
                    'nofilter': uses every single read turning off any filtering
    
    Generates consensus output file and yaml file with coverage for each region
    '''
    
    # get minimum umi family sizes
    family_sizes = list(map(lambda x: int(x.strip()), fam_size.split(',')))

    # get consensus data for each f_size + uncollapsed data
    print("Building consensus data...")
    cons_data = {}
    for f_size in family_sizes:
        # check if 0 is passed as fam_size argument
        if f_size == 0:
            # compute consensus for uncollapsed data, and get coverage
            cons_data[f_size], coverage = generate_uncollapsed(contig, region_start, region_end, bam_file, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans, stepper=stepper)
        else:
            cons_data[f_size] = generate_consensus(umi_families, f_size, contig, region_start, region_end, bam_file, pos_threshold, consensus_threshold, count_threshold, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans, stepper=stepper)
    # compute consensus for uncollapsed data if not in fam_size argument, and get coverage
    if 0 not in family_sizes:
        cons_data[0], coverage = generate_uncollapsed(contig, region_start, region_end, bam_file, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans, stepper=stepper)

    # write output consensus file
    print("Writing output...")
    raw_table_output(cons_data, contig, region_start, region_end, outdir)
    # save coverage to a yaml in outdir/Stats    
    StatsDir = os.path.join(os.path.dirname(outdir), 'Stats')
    if os.path.isdir(StatsDir) == False:
        os.mkdir(StatsDir)
    covdata = {contig + ':' + str(region_start+1) + '-' + str(region_end): coverage}
        
    with open(os.path.join(StatsDir, 'CoverageStats.yml'), 'a') as newfile:
        yaml.dump(covdata, newfile, default_flow_style=False)
    
    