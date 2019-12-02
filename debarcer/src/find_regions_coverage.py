# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 09:42:42 2019

@author: rjovelin
"""

import pysam


def GetContigs(bamfile):
    '''
    (str)- > list
    
    :param bamfile: Path to the bam file
    
    Returns a list of contig names in the bam header
        
    Precondition: bamfile is coordinate-sorted and has 'SQ' fields
    '''
    
    infile = pysam.AlignmentFile(bamfile)
    # convert header object to dict
    header = dict(infile.header)
    # make a list of contigs
    chromo = [i['SN'] for i in header['SQ']]
    infile.close()
    return chromo 


def FindRegionsCoverage(bamfile, contig, min_cov, region_size, max_depth, ignore_orphans, stepper):
    '''
    (str, str, int, int, int, bool, str) -> dict
    
    :param bamfile: Path to the bam file
    :param contig: Chromosome name, eg. chrN
    :param min_cov: Minimum read depth for all positions in genomic interval
    :param region_size: Minimum length of the genomic interval    
    :param max_depth: Maximum read depth
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair) if True
    :param stepper: Controls how the iterator advances. Accepeted values:
                    'all': skip reads with following flags: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
                    'nofilter': uses every single read turning off any filtering    
    
    Find all genomic intervals on contig of minimum length region_size for which
    all positions have read depth equals to min_cov or greater
        
    Precondition: bamfile is coordinate-sorted and has 'SQ' fields    
    '''
    
    # open file for reading
    infile = pysam.AlignmentFile(bamfile, 'rb')
    # store region coordinates in a dict {region_number: [pos, end]
    D = {}
    # initiate list with region coordinates
    L = []
    # initiate region number    
    n = 1
        
    for pileupcolumn in infile.pileup(contig, max_depth=max_depth, ignore_orphans=ignore_orphans, stepper=stepper):
        # get column position
        pos = int(pileupcolumn.reference_pos)
        # compare number of reads at pileup position with minimum read depth required
        if pileupcolumn.nsegments >= min_cov:
            # read depth is greater than requirement, check if region has been recorded
            if len(L) == 0:
                # new region, record start position, set end to 0 
                L = [pos, 0]
            else:
                # update end position
                L[-1] = pos
        else:
            # read depth requirement not met, record current region, and start new region
            D[n] = L
            L = []
            n += 1
            
    # remove empty lists
    to_remove = [i for i in D if len(D[i]) == 0]
    for i in to_remove:
        del D[i]
    # remove regions with length < region_size
    if len(D) != 0:
        to_remove = [i for i in D if len(range(D[i][0], D[i][-1]+1)) < region_size]
        for i in to_remove:
            del D[i]
    infile.close()
    return D
    


def WriteTargetsBed(bamfile, outputfile, min_cov, region_size, max_depth, ignore_orphans, stepper):
    '''
    (str, str, int, int, int, bool, str) -> None
    
    :param bamfile: Path to the bam file
    :param outputfile: Path to the output bed file
    :param contig: Chromosome name, eg. chrN
    :param min_cov: Minimum read depth for all positions in genomic interval
    :param region_size: Minimum length of the genomic interval    
    :param max_depth: Maximum read depth
    :param ignore_orphans: Ignore orphan reads (paired reads not in proper pair) if True
    :param stepper: Controls how the iterator advances. Accepeted values:
                    'all': skip reads with following flags: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
                    'nofilter': uses every single read turning off any filtering    
        
    Write a bed file (1-based coordinates) with all genomic intervals of minimum
    length region_size for which all positions have read depth equals to min_cov or greater
        
    Precondition: bamfile is coordinate-sorted and has 'SQ' fields    
    '''
    
    Regions = {}
    # make a list of chromosomes
    chromos = GetContigs(bamfile)
    # loop over chromosomes
    for contig in chromos:
        
        print(contig)
        
        # find genomic intervals on given chromosome
        D = FindRegionsCoverage(bamfile, contig, min_cov, region_size, max_depth, ignore_orphans, stepper)
        # collect intervals for chromo
        if len(D) != 0:
            Regions[contig] = D
    
    
    print(len(Regions))
    for i in Regions:
        print(i, len(Regions[i]))
    
    
    
    # write header to bed file
    newfile = open(outputfile, 'w')
    # loop over all chromosomes in bam header
    for contig in chromos:
        # check that genomic intervals are found on chromo
        if contig in Regions:
            for interval in sorted(Regions[contig]):
                # adjust interval positions to be 1-based
                start, end = Regions[contig][interval][0] + 1, Regions[contig][interval][1] + 1
                line = [contig, str(start), str(end)]
                newfile.write('\t'.join(line) + '\n')
    newfile.close()
    
