# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 10:05:10 2019

@author: rjovelin
"""

import configparser
import os


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
    "param filetype: Type of input file expected in the config.
                     Accepted values: prep_file or bam_file
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
            raise ValueError('ERR: Invalid path to input {0}'.format(filetype.replace('_', ' ')))
        elif os.path.isfile(inputfile) == False:
            raise ValueError('ERR: Invalid path to input {0}'.format(filetype.replace('_', ' ')))
    
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