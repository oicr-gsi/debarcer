# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 10:05:10 2019

@author: rjovelin
"""

import configparser
import os
import itertools
import numpy as np

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


def GetInputFiles(configfile, file, filetype):
    '''
    (str, str, str) -> str
    
    :param configfile: Path to config file
    :param file: Path to input file from the command
    "param filetype: Type of input file expected in the config. Accepted values:
                     'prep_file', 'bam_file', 'reference_file'
    '''
    
    # get file from config in priority
    try:
        config = configparser.ConfigParser()
        config.read(configfile)
        inputfile = config['PATHS'][filetype]
    except:
        # get inputfile from command
        inputfile = file 
    finally:
        # check that input is a valid file
        if inputfile in ['', None]:
            raise ValueError('ERR: Invalid path to input {0}'.format(filetype))
        elif os.path.isfile(inputfile) == False:
            raise ValueError('ERR: Invalid path to input {0}'.format(filetype))
    return inputfile

def GetOutputDir(configfile, outputdir):
    '''
    (str, str) -> str
    
    :param configfile: Path to config file
    :param outputdir: Path to output directory from the command
    
    Return the path to output directory or raise a ValueError if path
    not provided or path doesn't exist
    '''
    
    # get outdir from config in priority
    try:
        config = configparser.ConfigParser()
        config.read(configfile)
        outdir = config['PATHS']['outdir']
    except:
        # get outdir from the command
        outdir = outputdir
    finally:
        if outdir in ['', None]:
            # set output directy as current directory
            outdir = os.getcwd()
    return outdir



def GetThresholds(configfile, parameter, threshold):
    '''
    (str, str) -> float
    
    :param configfile: Path to config file
    :param parameter: Parameter name in config file. Accepted values:
                      'umi_family_pos_threshold'
                      'umi_edit_distance_threshold'
                      'percent_consensus_threshold'
                      'count_consensus_threshold'
    :param threshold: Setting threshold passed from command    
    
    Return a setting threshold
    '''
    
    # get parameter from the config file in priority
    if parameter in ['umi_family_pos_threshold', 'umi_edit_distance_threshold',
                     'percent_consensus_threshold', 'count_consensus_threshold']:
        Level = 'SETTINGS'
    elif parameter in ['percent_ref_threshold', 'percent_allele_threshold']:
        Level = 'REPORT'
        
    
    if parameter in ['umi_family_pos_threshold', 'umi_edit_distance_threshold']:
        try:
            config = configparser.ConfigParser()
            config.read(configfile)
            ThresholdVal = int(config[Level][parameter])
        except:
            # check if threshold provided in command
            try:
                ThresholdVal = int(threshold)
            except:
                # raise error and exit
                raise ValueError('ERR: Missing setting threshold')
        finally:
            # check that threshold is float
            if type(ThresholdVal) != int:
                raise ValueError('ERR: Setting threshold should be integer')
    elif parameter in ['percent_consensus_threshold', 'count_consensus_threshold', 'percent_ref_threshold', 'percent_allele_threshold']:
        try:
            config = configparser.ConfigParser()
            config.read(configfile)
            ThresholdVal = float(config[Level][parameter])
        except:
            # check if threshold provided in command
            try:
               ThresholdVal = float(threshold)
            except:
                # raise error and exit
                raise ValueError('ERR: Missing setting threshold')
        finally:
            # check that threshold is float
            if type(ThresholdVal) != float:
                raise ValueError('ERR: Setting threshold should be float')
    return ThresholdVal


def GetFamSize(configfile, famsize):
    '''
    (str, str) -> str
    
    :param configfile: Path to config file
    :param famsize: Comma-separated list of family size
      
    Returns a comma-separated flist of minimum family size for grouping and collapsing
    '''
        
    # get comma-separated list of minimum family sizes 
    try:
        config = configparser.ConfigParser()
        config.read(configfile)
        fam_size = config['SETTINGS']['min_family_sizes']
    except:
        # check if provided in command
        fam_size = famsize
    finally:
        # check if fam_size is defined
        if fam_size in [None, '']:
            raise ValueError('ERR: Missing minimum family sizes')
    return fam_size


def edit_distance(a, b):
    '''
    (str, str) -> int
    
    :param a": String of length L
    :param b": String of length L
    
    Returns the Hamming edit distance between a and b
    Pre-condition: a and b have identical length
    '''
    
    return sum(letter_a != letter_b for letter_a, letter_b in zip(a, b))



def get_consecutive_items(L):
    '''
    (list) -> generator
    
    :param L: List
    
    Return a generator with 1st and last item of consecutive items in L
    '''
    
    # remove duplicate and sort L
    L = sorted(set(L))
    for key, group in itertools.groupby(enumerate(L), lambda t: t[1] - t[0]):
        group = list(group)
        # return a generator with 1st and last item of consecutive items in L
        yield group[0][1], group[-1][1]


### use this function to extract regions from a merged consensus file
#def ExtractRegions(ConsFile):
#    '''
#    (str) -> dict
#    
#    :param ConsFile: Path to consensus file (merged or not)
#    
#    Returns all intervals in ConsFile in the form a dictionary with interval range
#    (1-based coordinates) and a list of chromosomes with these coordinates 
#    '''
#    
#    # extract all positions per chromosome
#    D = {}
#    
#    # open file for reading, skip header
#    infile = open(ConsFile)
#    Header = infile.readline().rstrip().split()
#    for line in infile:
#        # consider lines with valid chromosome format
#        if line.startswith('chr'):
#            line = line.rstrip().split()
#            chromo, pos = line[Header.index('CHROM')], int(line[Header.index('POS')])
#            if chromo in D:
#                D[chromo].append(pos)
#            else:
#                D[chromo] = [pos]
#    infile.close()
#    
#    # find intervals for each chromo {range(start, end) : [chromos]}
#    I = {}
#    for chromo in D:
#        # get the consecutive positions for each chromosome
#        intervals = list(get_consecutive_items(D[chromo]))
#        # convert positions to ranges 
#        for i in intervals:
#            # adjust last position, range is not inclusive
#            positions = range(i[0], i[1] + 1)
#            if positions in I:
#                I[positions].append(chromo)
#            else:
#                I[positions] = [chromo]
#    return I    
        

#def FormatRegion(ConsFile):
#    '''
#    (str) -> str
#    
#    :param ConsFile: Path to consensus file (not merged)
#    
#    Returns a region in the format 'chr:start-end'  
#    
#    Pre-condition: The consensus file is not merged and coordinates are from a single region
#    '''
#    
#    
#    # extract interval from file
#    I = ExtractRegions(ConsFile)
#    # get start and end positions
#    start, end = list(I.keys())[0][0], list(I.keys())[0][-1]
#    # get chromosome
#    chromo = I[list(I.keys())[0]][0]
#    # adjust end to get region
#    region = chromo + ':' + str(start) + '-' + str(end + 1)
#    return region




def FormatRegion(File):
    '''
    (str) -> str
    
    :param ConsFile: Path to file (not merged)
    
    Returns a region in the format 'chr:start-end'  
    Pre-condition: File is named *chr*:*-*.* 
    '''
    
    region = os.path.basename(File)
    if 'chr' in region and '.' in region:
        region = region[region.index('chr'):region.index('.')]
        if '_' in region:
            region = region.replace('_', ':')
    else:
        region = 'NA'
    return region


def GroupQCWriter(umi_positions, Outputfile):
    '''
    (dict, str) ->  None
    
    :param umi_positions: A dictionary with individual umi information (before grouping) output from get_umi_families function 
    :param Outputfile: Name of the output file
    
    Generates a table with information summarized from the umi data resulting from Grouping
    '''
    
    # umi_positions is a dict in the form {parent: {umi: {pos: count}}}
    
    # write a summary file of UMI relationships
    newfile = open(Outputfile, 'w')
    newfile.write('\t'.join(['UMI', 'Type', 'Count', 'Parent', 'Distance', 'Position']) + '\n')
    for parent in umi_positions:
        for umi in umi_positions[parent]:
            # compute hamming distance between i and umi
            distance = edit_distance(parent, umi)
            if distance == 0:
                # umi is parent
                umi_type = 'parent'
            else:
                # umi type is children
                umi_type = 'children'
            # compute median of umi occurence
            count = np.median([umi_positions[parent][umi][pos] for pos in umi_positions[parent][umi]]) 
            # position is the median of positions. remove chromosome from positions
            position = np.median(list(map(lambda x: int(x), list(map(lambda x: x.split(':')[1], list(umi_positions[parent][umi].keys()))))))
            # write umi info to file       
            newfile.write('\t'.join([umi, umi_type, str(count), parent, str(distance), str(position)]) + '\n')   
    newfile.close()    


def CreateDirTree(directory):
    '''
    (str) -> None
    
    :param directory: directory where d=subdirectories should be created
    
    Create a directory structure under directory
    '''
    
    # check that input directory is valid
    if os.path.isdir(directory) == False:
        raise ValueError('not a valid directory {0}'.format(directory))
    
    # create subdirectoy structure
    L = ['Umifiles', 'Datafiles', 'Stats', 'Consfiles', 'Qsubs', 'Figures', 'VCFfiles']
    Tree = []
    for i in L:
        j = os.path.join(directory, i)
        if os.path.isdir(j) == False:
            os.makedirs(j)
        Tree.append(j)
    # unpack variables
    UmiDir, DataDir, StatsDir, ConsDir, QsubDir, FigDir, VCFDir = Tree 
    # create log directory    
    LogDir = os.path.join(QsubDir, 'Logs')
    if os.path.isdir(LogDir) == False:
        os.makedirs(LogDir)
         
       
def CheckFileContent(File):
    '''
    (str) -> bool
    
    :param File: Path to File. First line is header
    
    Return True if content of the file (excluding header) is not the empty string
    '''
    
    infile = open(File)
    # skip header
    header = infile.readline()
    content = infile.read().rstrip()
    infile.close()
    if content != '':
        return True
    else:
        return False
    

def DropEmptyFiles(L):
    '''
    (list) -> list
    
    :param L: List of files
    
    return a modified list in which empty files (files with header but missing data) are removed
    '''
    
    # remove empty files
    to_remove = []
    for i in L:
        if CheckFileContent(i) == False:
            to_remove.append(i)
            print('ERR: Missing data in {0}'.format(i))
    for i in to_remove:
        L.remove(i)
    return L

def CheckFilePath(L):
    '''
    (list) -> None
    
    :param L: List of files
    
    Raise a ValueError if any file in L doesn't have a valid path
    '''
    
    # check that paths to files are valid
    for i in L:
        if os.path.isfile(i) == False:
            raise ValueError('ERR: Invalid file path {0}'.format(i))

