# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 14:25:48 2019

@author: RJovelin
"""

import time
from src.version import __version__


def GetConsData(consfile):
    '''
    (str) -> dict
    
    :param consfile: Path to the consensus file (merged or not)
           
    Returns a dictionary with consensus file info organized by chromo, and family
    size for each position
    '''

    # create a dict with consensus info per contig, region and umi family size
    # {contig :{fam: {pos: info}}}
    data = {}

    infile = open(consfile)
    header = infile.readline().rstrip().split('\t')
    
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            contig = line[header.index('CHROM')]
            if contig not in data:
                data[contig] = {}
            # get the position
            pos = int(line[header.index('POS')])
            # get umi fam size              
            famsize = int(line[header.index('FAM')])
            if famsize not in data[contig]:
                data[contig][famsize] = {}
            # collect info 
            assert pos not in data[contig][famsize]
            data[contig][famsize][pos] = line
    
    infile.close()
    return data
        
    

def WriteVCF(consfile, outputfile, reference, ref_threshold, alt_threshold, filter_threshold, famsize):
    '''
    (str, str, str, float, float, int, int) -> None
    
    :param consfile: Path to the consensus file (merged or not)
    :param outputfile: Path to the output VCF file
    :param reference" Path to the reference genome 
    :param ref_threshold: Maximum reference frequency (in %) to consider alternative variants
                          (ie. position with ref freq <= ref_threshold is considered variable)
    :param alt_threshold: Minimum allele frequency (in %) to consider an alternative allele at a variable position 
                          (ie. allele freq >= alt_threshold and ref freq <= ref_threshold --> record alternative allele)
    :param filter_threshold: minimum number of reads to pass alternative variants 
                             (ie. filter = PASS if variant depth >= alt_threshold)
    :param famsize: Minimum umi family size 
        
    Write a VCF from the consensus file for a given umi family size.
    Allow multiple records per position for SNVs and indels
    '''
    
    # parse consensus file -> consensus info for all recorded umi fam size
    consdata = GetConsData(consfile)

    # get the header of the consfile
    infile = open(consfile)
    header = infile.readline().rstrip().split('\t')
    infile.close()

    # get debarcer version
    version = __version__

    # create a list with VCF records for SNVs and indels, sorted by chromosome and position
    Records = []

    # loop over sorted contigs and sorted positions in cons data for given famsize 
    
    # make a sorted list of contigs
    Chromosomes = [i.replace('chr', '') for i in consdata.keys()]
    # make a list of non numerical contigs
    others = []
    for i in range(len(Chromosomes)):
        if Chromosomes[i].isnumeric() == False:
            others.append(Chromosomes[i])
    # sort other
    others.sort()
    # make list of numerical chromos
    if len(others) != 0:
        for i in others:
            Chromosomes.remove(i)
    Chromosomes = sorted(list(map(lambda x: int(x), Chromosomes)))
    Chromosomes = list(map(lambda x: 'chr' + str(x), Chromosomes))
    others = list(map(lambda x: 'chr' + str(x), others))
    # add back non-numerical contigs
    Chromosomes.extend(others)
    
    # make a sorted list of positions
    positions = []
    for i in consdata:
        for j in consdata[i]:
            positions.extend(list(consdata[i][j].keys()))
    positions = sorted(list(map(lambda x: int(x), list(set(positions)))))
       
    
    for contig in Chromosomes:
        # check membership of famsize
        if famsize in consdata[contig]:
            for pos in positions:
                # check pos membership for merged consensus files
                if pos in consdata[contig][famsize]:
                    L = consdata[contig][famsize][pos]
                    # get reference frequency
                    ref_freq = float(L[header.index('REF_FREQ')]) 
                    # create VCF record if ref freq low enough to consider variant at position 
                    if ref_freq <= ref_threshold:
                        # get consensus and raw depth       
                        consdepth = int(L[header.index('CONSDP')])
                        rawdepth = int(L[header.index('RAWDP')])
                        # get minimum and mean family size
                        minfam = int(L[header.index('FAM')])
                        meanfam = float(L[header.index('MEAN_FAM')])
                        # set up info
                        info = 'RDP={0};CDP={1};MIF={2};MNF={3};AD={4};AL={5};AF={6}'
        
                        # get the reference allele
                        ref = L[header.index('REF')]
                        # get the list of single nucleotides
                        alleles = [header[header.index(i)] for i in 'ACGTN']
                        
                        # make lists of deletions and counts
                        deletions = L[header.index('D_(ref,del)')].split(';')
                        delcounts = L[header.index('D_counts')].split(';')
                        
                        # make lists of insertions and counts
                        insertions = L[header.index('I_(ref,ins)')].split(';')
                        inscounts = L[header.index('I_counts')].split(';')
                        
                        # get the read depth for each allele and indels
                        depth = {i:int(L[header.index(i)]) for i in alleles}
                        if deletions != ['']:
                            for i in range(len(deletions)):                      
                                depth[deletions[i]] = int(delcounts[i])
                        if insertions != ['']:
                            for i in range(len(insertions)):
                                depth[insertions[i]] = int(inscounts[i])
                                               
                        # compute frequencies for each allele and indel
                        
                        
                        
                        if sum(depth.values()) == 0:
                            print(consfile)
                            print(contig)
                            print(pos)
                            print(depth)
                            print(deletions)
                            print(delcounts)
                            print(insertions)
                            print(inscounts)
                                
                                
                                
                                
                        
                        freq = {i: (depth[i]/sum(depth.values())) * 100 for i in depth}
                
                        # record snvs and indels on different lines
                                    
                        # make a list of alternative alleles with frequency >= alt_threshold
                        alt_alleles = [i for i in freq if i in alleles and i != ref and freq[i] >= alt_threshold]
                        # make a list of read depth for alternative alleles passing alt_threshold
                        alt_depth = [str(depth[i]) for i in alt_alleles]
                        # make a list of frequencies for alternative alelles passing alt_threshold 
                        alt_freq = [str(round(freq[i], 4)) for i in alt_alleles]
                        # record info
                        alt_info = info.format(rawdepth, consdepth, minfam, round(meanfam, 2), depth[ref], ','.join(alt_depth), ','.join(alt_freq))
                        
                        # make list of deletions with frequency >= alt_threshold
                        del_alleles = [i for i in freq if i in deletions and freq[i] >= alt_threshold]
                        # make a list of read depth for deletions passing alt_threshold
                        del_depth = [str(depth[i]) for i in del_alleles]
                        # make a list of frequencies for deletions passing alt_threshold
                        del_freq = [str(round(freq[i], 4)) for i in del_alleles]
                        
                        # make list of insertions with frequency >= alt_threshold
                        ins_alleles = [i for i in freq if i in insertions and freq[i] >= alt_threshold]
                        # make a list of read depth for insertions passing alt_threshold
                        ins_depth = [str(depth[i]) for i in ins_alleles] 
                        # make a list of frequencies for insertions passing alt_threshold
                        ins_freq = [str(round(freq[i], 4)) for i in ins_alleles]
                        
                        # check that alernative alleles are recorded
                        if len(alt_alleles) != 0:
                            # get the filter value based on min_read_depth
                            if True in [depth[i] >= filter_threshold for i in alt_alleles]:
                                filt = 'PASS' 
                            else:
                                filt = 'a{0}'.format(filter_threshold)
            
                            Records.append('\t'.join([contig, str(pos), '.', ref, ','.join(alt_alleles), '0', filt, alt_info]) + '\n')
                        # check that deletions are recorded
                        if len(del_alleles) != 0:
                            # record deletions seperately on distinct lines
                            for i in range(len(del_alleles)):
                                # get the filter value based on min_read_depth
                                if depth[del_alleles[i]] >= filter_threshold == True:
                                    filt == 'PASS'
                                else:
                                    filt = 'a{0}'.format(filter_threshold)
                                # record info
                                del_info = info.format(rawdepth, consdepth, minfam, round(meanfam, 2), depth[ref], del_depth[i], del_freq[i])
                                # extract ref allele and alt allele
                                k = list(map(lambda x: x.strip(), del_alleles[i].replace("'", '').replace('(', '').replace(')', '').split(',')))
                                Records.append('\t'.join([contig, str(pos), '.', k[0], k[1], '0', filt, del_info]) + '\n')
                        # check that insertions are recorded
                        if len(ins_alleles) != 0:
                            # record insertions seperately on distinct lines
                            for i in range(len(ins_alleles)):
                                # get the filter value based on min_read_depth
                                if depth[ins_alleles[i]] >= filter_threshold == True:
                                    filt == 'PASS'
                                else:
                                    filt = 'a{0}'.format(filter_threshold)
                                # record info
                                ins_info = info.format(rawdepth, consdepth, minfam, round(meanfam, 2), depth[ref], ins_depth[i], ins_freq[i])
                                # extract ref allele and alt allele
                                k = list(map(lambda x: x.strip(), ins_alleles[i].replace("'", '').replace('(', '').replace(')', '').split(',')))
                                Records.append('\t'.join([contig, str(pos), '.', k[0], k[1], '0', filt, ins_info]) + '\n')
    
    
    # write VCF only if positions are recorded
    if len(Records) != 0:
        # open file for writing
        newfile = open(outputfile, 'w')

        # write VCF header 
        newfile.write('##fileformat=VCFv4.1\n')
        newfile.write('##fileDate={0}\n'.format(time.strftime('%Y%m%d', time.localtime())))
        newfile.write('##reference={0}\n'.format(reference))
        newfile.write('##source=Debarcer v. {0}\n'.format(version))
            
        # write info/filter/format metadata
        newfile.write('##INFO=<ID=RDP,Number=1,Type=Integer,Description=\"Raw Depth\">\n')
        newfile.write('##INFO=<ID=CDP,Number=1,Type=Integer,Description=\"Consensus Depth\">\n')
        newfile.write('##INFO=<ID=MIF,Number=1,Type=Integer,Description=\"Minimum Family Size\">\n')
        newfile.write('##INFO=<ID=MNF,Number=1,Type=Float,Description=\"Mean Family Size\">\n')
        newfile.write('##INFO=<ID=AD,Number=1,Type=Integer,Description=\"Reference allele Depth\">\n')
        newfile.write('##INFO=<ID=AL,Number=A,Type=Integer,Description=\"Alternate Allele Depth\">\n')
        newfile.write('##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate Allele Frequency\">\n')
        newfile.write('##FILTER=<ID=a{0},Description=\"Alternate allele depth below {0}\">\n'.format(filter_threshold))
        
        # write data header 
        newfile.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']) + '\n')
        # write data
        newfile.write('\n'.join(Records))
        newfile.close()        

