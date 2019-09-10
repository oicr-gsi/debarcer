import os
import pysam
import configparser
import argparse
import operator
import json
from src.utilities import CheckRegionFormat, GetOutputDir, GetInputFiles, GetThresholds


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


def get_consensus_seq(umi_families, fam_size, ref_seq, contig, region_start, region_end, bam_file, pos_threshold, max_depth, truncate, ignore_orphans):
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
    :param max_depth: Maximum read depth
    :param truncate: Consider only pileup columns within interval defined by region start and end if True
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair) if True
    
    Returns a tuple with a dictionary representing consensus info at each base position in the given region
    and a dictionary to keep track of family size for each position
    '''


    # create a dict to store umi family size at each position {position: {parent: {distance: count}}}
    FamSize = {}

    # create a dict to store consensus seq info
    consensus_seq = {}
    
    with pysam.AlignmentFile(bam_file, "rb") as reader:
        # loop over pileup columns
        for pileupcolumn in reader.pileup(contig, region_start, region_end, truncate=truncate, ignore_orphans=ignore_orphans, max_depth=max_depth):
            # get column position. by default consider only positions within region
            # however, number of reads in families consider reads overlapping with region
            # not only contained within region
            pos = int(pileupcolumn.reference_pos)  
            
            print('pos', pos)
            
            
            
            assert pos != region_end
            # loop over reads in pileup column
            for read in pileupcolumn.pileups:
                # get read information
                read_data = read.alignment
                read_name, start_pos = read_data.query_name, int(read_data.reference_start)
                # get all recorded umis
                umis = read_name.split(":")[-1].split(';')
                
                
                print(umis)
                print(read)
                print(read_name, start_pos)
                
                
                for umi in umis:
                    
                    print(umi, umi in umi_families, len(umi_families))
                    
                    
                    # check that umi is recorded
                    if umi in umi_families:
                        
                        
                        print(umi, umi in umi_families)
                        
                        
                        # find closest family from umi
                        # make a list of (positions counts)
                        L = [(int(i.split(':')[1]), umi_families[umi]['positions'][i]) for i in umi_families[umi]['positions']]
                        closest, count = find_closest(start_pos, L)



                        print(closest, pos_threshold, closest <= pos_threshold)
                        print(count, fam_size, count >= fam_size)


                        # check if closest family is within the position threshold
                        if closest <= pos_threshold:
                            # found a umi family. check if family count is greater than family threshold
                            if count >= fam_size:
                                
                                
                                
                                print('umi', umi, 'closest', closest, 'pos_threshold', pos_threshold, closest <= pos_threshold,   'closest < threshold', 'count', count, 'fam_size', fam_size, count >= fam_size)
                                
                                
                                
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


def get_uncollapsed_seq(ref_seq, contig, region_start, region_end, bam_file, max_depth, truncate, ignore_orphans):
    '''
    (str, str, int, int, str, str, int, bool, bool) -> dict
    
    :param ref_seq: Sequence of the reference corresponding to the given region
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param bam_file: Path to the bam file
    :param max_depth: Maximum read depth
    :param truncate: Consider only pileup columns within interval defined by region start and end if True
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair) if True
    
    Returns a nested dictionary representing counts of each base at each base position.
    '''

    uncollapsed_seq = {}

    with pysam.AlignmentFile(bam_file, "rb") as reader:
        # loop over pileup columns 
        for pileupcolumn in reader.pileup(contig, region_start, region_end, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans):
           
            
            
            print('uncollapsed', pileupcolumn)
            
            
            
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
    
    
def generate_consensus(umi_families, fam_size, ref_seq, contig, region_start, region_end, bam_file, pos_threshold, percent_threshold, count_threshold, max_depth, truncate, ignore_orphans):
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
    :param max_depth: Maximum read depth
    :param truncate: Consider only pileup columns within interval defined by region start and end if True
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair) if True
    
    Generates consensus data for a given family size and genomic region
    '''

    # get consensus info for each base position and umi group in the given region {pos: {fam_key: {(ref, alt):count}}}
    # get family size at each position 
    consensus_seq, FamSize = get_consensus_seq(umi_families, fam_size, ref_seq, contig, region_start, region_end, bam_file, pos_threshold, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans)




    print('cons seq', len(consensus_seq))
    print('fam size', len(FamSize))




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
    
    
    print('cons data per base', len(cons_data))
    
    
    return cons_data



def generate_uncollapsed(ref_seq, contig, region_start, region_end, bam_file, max_depth, truncate, ignore_orphans):
    '''
    (str, str, int, int, str, int, bool, bool) -> dict
    
    :param ref_seq: Sequence of the reference corresponding to the given region
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param bam_file: Path to the bam file
    :param max_depth: Maximum read depth
    :param truncate: Consider only pileup columns within interval defined by region start and end if True
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair) if True
    
    Generates uncollapsed consensus data for the genomic region
    '''
    
    # get uncolapased seq info {pos: {(ref, atl): count}}
    uncollapsed_seq = get_uncollapsed_seq(ref_seq, contig, region_start, region_end, bam_file, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans)
    
    
    
    print(ref_seq, contig, region_start, region_end, bam_file, max_depth, truncate, ignore_orphans)
    
    print('uncollapsed')
    print(uncollapsed_seq)
    
    
    
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


def raw_table_output(cons_data, ref_seq, contig, region_start, region_end, outdir, ref_threshold, all_threshold):
    '''
    (dict, str, str, int, int, str, num, num) -> None
    
    :param ref_seq: Sequence of the reference corresponding to the given region
    :param contig: Chromosome name, eg. chrN
    :param region_start: Start index of the region of interest. 0-based half opened
    :param region_end: End index of the region of interest. 0-based half opened
    :param outdir: Output directory
    
    :param ref_threshold: ??
        
    :param all_threshold: ??
    
    Writes a long-form consensus file for every event detected in the collapsed data
    '''
    
    # get the path to the output file
    OutputFile = os.path.join(outdir, '{}:{}-{}.cons'.format(contig, region_start + 1, region_end))
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


def generate_consensus_output(reference, contig, region_start, region_end, bam_file, umi_families, outdir, fam_size, pos_threshold, percent_threshold, count_threshold, ref_threshold, all_threshold, max_depth, truncate, ignore_orphans):
    '''
    (str, str, int, int, str, dict, str, str, int, num, int, num, num, int, bool, bool) -> None
    
    
    :param reference: Path to the reference genome
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
    :param percent_threshold: Percent consensus threshold for alternative base within family 
    :param count_threshold: Count consensus threshold for alternative base within family
    :param ref_threshold: ??
    :param all_threshold: ??
    :param max_depth: Maximum read depth
    :param truncate: Consider only pileup columns within interval defined by region start and end if True
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair) if True
    
    Generates consensus output file
    '''
    
    # get minimum umi family sizes
    family_sizes = list(map(lambda x: int(x.strip()), fam_size.split(',')))


    print(family_sizes)



    # get reference sequence for the region 
    print("Getting reference sequence...")
    ref_seq = get_ref_seq(contig, region_start, region_end, reference)

    # get consensus data for each f_size + uncollapsed data
    print("Building consensus data...")
    cons_data = {}
    for f_size in family_sizes:
        
        print(f_size, type(f_size))
        
        
        
        # check if 0 is passed as fam_size argument
        if f_size == 0:
            # compute consensus for uncollapsed data
            cons_data[f_size] = generate_uncollapsed(ref_seq, contig, region_start, region_end, bam_file, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans)
        else:
            cons_data[f_size] = generate_consensus(umi_families, f_size, ref_seq, contig, region_start, region_end, bam_file, pos_threshold, percent_threshold, count_threshold, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans)
    # compute consensus for uncollapsed data if not in fam_size argument
    if 0 not in family_sizes:
        cons_data[0] = generate_uncollapsed(ref_seq, contig, region_start, region_end, bam_file, max_depth=max_depth, truncate=truncate, ignore_orphans=ignore_orphans)



    print('cons data', len(cons_data))
    for i in cons_data:
        print(i, len(cons_data[i]))



    # write output consensus file
    print("Writing output...")
    raw_table_output(cons_data, ref_seq, contig, region_start, region_end, outdir, ref_threshold, all_threshold)
    

if __name__=="__main__":
    
    # Base collapse command
    parser = argparse.ArgumentParser(description='generate_consensus.py', help='Base collapsing from given UMI families file')
    parser.add_argument('-c', '--Config', dest='config', help='Path to the config file')
    parser.add_argument('-o', '--Outdir', dest='outdir', help='Directory where consensus files are written')
    parser.add_argument('-b', '--Bamfile', dest='bamfile', help='Path to the BAM file')
    parser.add_argument('-rf', '--Reference', dest='reference', help='Path to the refeence genome')
    parser.add_argument('-r', '--Region', dest='region', help='Region coordinates to search for UMIs. chrN:posA-posB. posA and posB are 1-based included', required=True)
    parser.add_argument('-u', '--Umi', dest='umifile', help='Path to the .umis file', required=True)
    parser.add_argument('-f', '--Famsize', dest='famsize', help='Comma-separated list of minimum umi family size to collapase on')
    parser.add_argument('-ct', '--CountThreshold', dest='countthreshold', help='Base count threshold in pileup column')
    parser.add_argument('-pt', '--PercentThreshold', dest='percentthreshold', help='Base percent threshold in pileup column')
    parser.add_argument('-rt', '--RefThreshold', dest='refthreshold', help='Reference threshold')
    parser.add_argument('-at', '--AlleleThreshold', dest='allthreshold', help='Allele threshold')
    parser.add_argument('-p', '--Position', dest='postthreshold', help='Umi position threshold for grouping umis together')
    parser.add_argument('-m', '--MaxDepth', dest='maxdepth', default=1000000, help='Maximum read depth. Default is 1000000')
    parser.add_argument('-t', '--Truncate', dest='truncate', action='store_false',
                       help='If truncate is True and a region is given,\
                       only pileup columns in the exact region specificied are returned. Default is True')
    parser.add_argument('-i', '--IgnoreOrphans', dest='ignoreorphans', action='store_false',
                       help='Ignore orphans (paired reads that are not in a proper pair). Default is True')
    
    args = parser.parse_args()
    
    # get output directory from the config or command. set to current dir if not provided
    outdir = GetOutputDir(args.config, args.outdir)
    # create outputdir if doesn't exist
    if os.path.isdir(outdir) == False:
        if os.path.isfile(outdir) == True:
            raise ValueError('ERR: Output directory cannot be a file')
        else:
            os.makedirs(outdir)
    
    # get input bam from config or command
    bam_file = GetInputFiles(args.config, args.bamfile, 'bam_file')

    # check that region is properly formatted
    region = args.region
    CheckRegionFormat(region)
    # get chromosome 
    contig = region.split(":")[0]
    # get 1-based inclusive region coordinates
    region_start, region_end = int(region.split(":")[1].split("-")[0]), int(region.split(":")[1].split("-")[1])
    # convert coordinates to 0-based half opened coordinates
    region_start = region_start -1
    
    # load json with count of umi families per position and umi group
    try:
        infile = open(args.umifile)
        umi_families = json.load(infile)
        infile.close()
    except:
        raise ValueError("ERR: Unable to load .umi json file")
        
    # get percent threshold 
    percent_threshold = GetThresholds(args.config, 'percent_consensus_threshold', args.percentthreshold)
    # get count threshold
    count_threshold = GetThresholds(args.config, 'count_consensus_threshold', args.countthreshold)
    # get reference threshold
    ref_threshold = GetThresholds(args.config, 'percent_ref_threshold', args.refthreshold)
    # get allele threshold
    all_threshold = GetThresholds(args.config, 'percent_allele_threshold', args.allthreshold)
    # get umi position threshold 
    pos_threshold = GetThresholds(args.config, 'umi_family_pos_threshold', args.postthreshold)
    
    # get reference
    reference = GetInputFiles(args.config, args.reference, 'reference_file')
    
    # get comma-separated list of minimum family sizes 
    try:
        config = configparser.ConfigParser()
        config.read(args.config)
        fam_size = config['SETTINGS']['min_family_sizes']
    except:
        # check if provided in command
        fam_size = args.famsize
    finally:
        # check if fam_size is defined
        if fam_size in [None, '']:
            raise ValueError('ERR: Missing minimum family sizes')
    
    # write consensus output file
    generate_consensus_output(reference, contig, region_start, region_end, bam_file, umi_families, outdir, fam_size, pos_threshold, percent_threshold, count_threshold, ref_threshold, all_threshold, max_depth=args.maxdepth, truncate=args.truncate, ignore_orphans=args.ignoreorphans)
 