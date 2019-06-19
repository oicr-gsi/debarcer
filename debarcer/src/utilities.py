# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 10:05:10 2019

@author: rjovelin
"""



def CheckRegionFormat(region):
    '''
    (str) -> None
    
    :param region: A string with expected format chrN:posA-posB. posA and posB are 1-based inclusive
    
    Checks that region is properly formatted and raise ValueError if not
    '''
    
    if any(i not in region for i in ["chr", ":", "-"]):
        raise ValueError('ERR: Incorrect region string (should look like chr1:1200000-1250000)')
    
    # get chromosome and check format 
    contig = region.split(":")[0]
    chromos = [str(i) for i in range(23)] + ['X', 'Y']
    if contig[:len('chr')] != 'chr' and contig[len('chr'):] not in chromos:
        raise ValueError('ERR: Incorrect chromosome name (should look like chr1:1200000-1250000)')
    region_start, region_end = region.split(":")[1].split("-")[0], region.split(":")[1].split("-")[1]
    if region_start.isnumeric() == False or region_end.isnumeric() == False:
        raise ValueError('ERR: Incorrect start and end coordinates (should look like chr1:1200000-1250000)')
    