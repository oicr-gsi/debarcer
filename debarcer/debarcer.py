import argparse
import sys
import os
import json
import mistune
import collections 


from debarcer.preprocess_fastqs import reheader_fastqs, check_library_prep
from debarcer.umi_error_correct import get_umi_families, umi_datafile
from debarcer.generate_consensus import generate_consensus_output
from debarcer.generate_vcf import WriteVCF
from debarcer.run_analyses import submit_jobs
from debarcer.utilities import CheckRegionFormat, GetOutputDir, GetInputFiles, GetThresholds, GetFamSize, \
 FormatRegion, GroupQCWriter, CreateDirTree, DropEmptyFiles, CheckFilePath, ConvertArgToBool, GetCurrentTime, \
 get_read_count, MergeDataFiles, MergeUmiFiles, MergeConsensusFiles
from debarcer.generate_plots import PlotMeanFamSize, PlotNonRefFreqData, PlotConsDepth,\
 PlotParentsToChildrenCounts, PlotParentFreq, PlotNetworkDegree, PlotUMiFrequency,\
 GetUmiCountFromPreprocessing, PlotFamSizeReadDepth, PlotReadDepth, GetIndividualUmiInfo,\
 PlotIncorrectReads, PlotDataPerRegion
from debarcer.generate_report import WriteReport   
from debarcer.find_regions_coverage import WriteTargetsBed

import matplotlib.pyplot as plt
import numpy as np



"""
debarcer.py - main interface for Debarcer
=========================================

Purpose
-------
Debarcer (De-Barcoding and Error Correction) is a package 
for working with next-gen sequencing data that contains
molecular barcodes.

Authors: Alexandra Bodak, Lawrence Heisler, Richard Jovelin, Isha Warikoo 
Copyright (c) 2018 GSI, Ontario Institute for Cancer Research
"""



def preprocess_reads(outdir, read1, read2, read3, prepname, prepfile, config, prefix):
    '''
    (list) -> None
    Preprocesses fastq files by removing UMIs from reads, appending them to 
    the read names and writing new fastqs and also write QC reports in a Stats directory
	 
    Parameters
    ----------
    - outdir (str): Directory where new fastqs are written. From command or config
    - read1 (str): Path to first FASTQ file
    - read2 (str or None): Path to second FASTQ file
    - read3 (str or None): Path to third FASTQ file
    - prepname (str): Name of the library preparation
    - prepfile (str): Path to the library preparation ini file. From command or config
    - config (str or None): Path to the config file. 
    - prefix (str): Prefix for naming umi-reheradered fastqs 
    '''
        
    # get output directory from the config or command. set to current dir if not provided
    outdir = GetOutputDir(config, outdir)
    # create outputdir
    if os.path.isfile(outdir) == True:
        raise ValueError('ERR: Output directory cannot be a file')
    else:
        os.makedirs(outdir, exist_ok=True)
    
    # get input prep file from config or command
    prepfile = GetInputFiles(config, prepfile, 'prep_file')
    
    # check library prep format
    check_library_prep(prepname, prepfile)
        
    # check that files are valid
    for i in [read1, read2, read3]:
        # check that argument (or file) is provided/exists 
        if i:
            # check if provided file is file
            if os.path.isfile(i) == False:
                # raise error and exit
                raise ValueError('ERR: read file is not a valid file')
                sys.exit(1)
    
    # reheader fastqs and add umi in new fastqs header
    # get a json with read counts and list of umi sequences
    D, UmiSequences = reheader_fastqs(read1, outdir, prepname, prepfile, r2=read2, r3=read3, prefix=prefix)
	 
    # create subdirectoy structure
    CreateDirTree(outdir)
    
    # write json with correct/incorrect reads
    with open(os.path.join(outdir, 'Stats/Processing_Read_Info.json'), 'w') as newfile:
        json.dump(D, newfile, indent=4)
       
    # count umi sequences
    umi_counts = dict(collections.Counter(UmiSequences))
    # write table with umi sequences counts
    Outpufile = os.path.join(outdir, 'Stats/Umi_counts.txt')
    newfile = open(Outpufile, 'w')
    newfile.write('\t'.join(['Umi', 'Count']) + '\n')
    for i in umi_counts:
        newfile.write('\t'.join([i, str(umi_counts[i])]) + '\n')
    newfile.close()



def group_umis(outdir, region, bamfile, config, distthreshold, postthreshold, ignore, truncate, separator, readcount):
    '''
    (str, str, str, str, int, int, bool, bool, str, int) -> None
    
    Groups by hamming distance and form families based on physical distances within groups
    
    Parameters
    ----------
    - outdir (str): Output directory where subdirectories are created
    - region (str): A string with region coordinates chrN:posA-posB. posA and posB are 1-based included
    - bamfile (str): Path to the bam file
    - config (str): Path to your config file
    - distthreshold (int): Hamming distance threshold for connecting parent-children umis
    - postthreshold (int): Distance threshold in bp for defining families within groups
    - ignore (bool): Keep the most abundant family and ignore families at other positions within each group if True. Default is False
    - truncate (bool): Skip reads overlapping with the genomic interval if True. Default is False
    - separator (str): String separating the UMI from the remaining of the read name
    - readcount (int): Minimum number of reads in region required for grouping. Default is 0  
    '''
    
    # get output directory from the config or command. set to current dir if not provided
    outdir = GetOutputDir(config, outdir)
    # create outputdir
    if os.path.isfile(outdir) == True:
        raise ValueError('ERR: Output directory cannot be a file')
    else:
        os.makedirs(outdir, exist_ok=True)
    
    # create subdirectoy structure
    CreateDirTree(outdir)

    # get input bam from config or command
    bam_file = GetInputFiles(config, bamfile, 'bam_file')
    
    # check that region is properly formatted
    CheckRegionFormat(bam_file, region)
    # get chromosome
    contig = region.split(":")[0]
    # get 1-based inclusive region coordinates   
    region_start, region_end = int(region.split(":")[1].split("-")[0]), int(region.split(":")[1].split("-")[1])
    # convert coordinates to 0-based hal opened coordinates
    region_start = region_start -1
    
    # get umi position and distance thresholds 
    pos_threshold = GetThresholds(config, 'umi_family_pos_threshold', postthreshold)
    dist_threshold = GetThresholds(config, 'umi_edit_distance_threshold', distthreshold)
          
    print(GetCurrentTime() + "Grouping UMIs...")
    
    # count reads in genomic region
    read_count = get_read_count(bam_file, contig, region_start, region_end)
    
    print(GetCurrentTime() + "{0} reads in {1}".format(read_count, region))
        
    # check that read count is greater than threshold
    if read_count > readcount:
        # Generate UMI families within groups using the position of the most frequent umi as reference for each family
        # keep the most abundant family within group and ignore others if args.ignore is True
        umi_families, umi_groups, umi_positions, mapped_reads = get_umi_families(contig, region_start, region_end, bam_file, pos_threshold, dist_threshold, ignore, truncate, separator)
    
        # get the number of parent umis, number of children and number of parent given a number of children
        filename= os.path.join(outdir, 'Datafiles/datafile_{}.csv'.format(region))
        header = ['CHR', 'START', 'END', 'PTU', 'CTU', 'CHILD_NUMS', 'FREQ_PARENTS']
        # use 1-based inclusive coordinates in datafile output
        info = [contig, str(region_start + 1), str(region_end)] + umi_datafile(umi_groups)
        with open(filename, 'w') as newfile:
            newfile.write('\t'.join(header) + '\n')
            newfile.write('\t'.join(info) + '\n')
    
        # save umi families as a json. positions in the json are 0-based half opened
        umi_file = os.path.join(outdir, 'Umifiles/{0}.json'.format(region))
        with open(umi_file, 'w') as newfile:
            json.dump(umi_families, newfile, sort_keys = True, indent=4)
    
        # write a summary file of UMI relationships and count for each individual umi before grouping
        Outputfile = os.path.join(outdir, 'Stats/UMI_relationships_{0}.txt'.format(region))
        GroupQCWriter(umi_positions, Outputfile)
    
        # save information about individual UMIs as a json in Stats directory
        umi_file = os.path.join(outdir, 'Stats/Umis_{0}_before_grouping.json'.format(region))
        with open(umi_file, 'w') as newfile:
            json.dump(umi_positions, newfile, sort_keys = True, indent=4)
    
        # save counts of unmapped and mapped reads as a json in Stats directory
        read_file = os.path.join(outdir, 'Stats/Mapped_read_counts_{0}.json'.format(region))
        with open(read_file, 'w') as newfile:
            json.dump(mapped_reads, newfile, sort_keys = True, indent=4)
        
        print(GetCurrentTime() + 'UMI grouping complete. CSV files written to {0}.'.format(os.path.join(outdir, 'Datafiles')))
        print(GetCurrentTime() + 'UMI grouping complete. UMI files written to {0}.'.format(os.path.join(outdir, 'Umifiles')))
        print(GetCurrentTime() + 'UMI grouping complete. QC files written to {0}.'.format(os.path.join(outdir, 'Stats')))
    
    else:
        print(GetCurrentTime() + 'Not enough reads in region {0}. Found {1} reads and more than {2} reads are required for grouping'.format(region, read_count, readcount))

def collapse(config, outdir, bamfile, reference, region, umifile, famsize,
             countthreshold, percentthreshold, postthreshold, maxdepth, truncate,
             ignoreorphans, stepper, separator, base_quality_score):
    '''
    (str, str, str, str, str, str, str, str, str, str, int, bool, bool, bool, str, int) -> None
    
    Base collapses from given BAM and umi family file
    
    Parameters
    ----------
    - config (str): Path to the config file
    - outdir (str): Output directory where subdirectories are created
    - bamfile (str): Path to the BAM file
    - reference (str): Path to the reference genome
    - region (str): Region coordinates (1-based included) to search for UMIs (eg. chrN:posA-posB)
    - umifile (str): Path to the .umis file
    - famsize (str): Comma-separated list of minimum umi family size to collapase on
    - countthreshold (str): Base count threshold in pileup column
    - percentthreshold (str): Base percent threshold in pileup column
    - postthreshold (str): Umi position threshold for grouping umis together
    - maxdepth (int): Maximum read depth. Default is 1000000
    - truncate (bool): Only consider pileup columns in given region. Default is False
    - ignoreorphans (bool): Ignore orphans (paired reads that are not in a proper pair). Default is True
    - stepper (bool): Controls how the iterator advances. Accepeted values:
                      'all': skip reads with following flags: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
                      'nofilter': uses every single read turning off any filtering
    - separator (str): String separating the UMI from the remaining of the read name
    - base_quality_score (int): Base quality score threshold. No offset of 33 needs to be subtracted
    '''
    
    # get output directory from the config or command. set to current dir if not provided
    outdir = GetOutputDir(config, outdir)
    # create outputdir
    if os.path.isfile(outdir) == True:
        raise ValueError('ERR: Output directory cannot be a file')
    else:
        os.makedirs(outdir, exist_ok=True)
    
    # create subdirectoy structure
    CreateDirTree(outdir)

    # get input bam from config or command
    bam_file = GetInputFiles(config, bamfile, 'bam_file')

    # check that region is properly formatted
    CheckRegionFormat(bam_file, region)
    # get chromosome 
    contig = region.split(":")[0]
    # get 1-based inclusive region coordinates
    region_start, region_end = int(region.split(":")[1].split("-")[0]), int(region.split(":")[1].split("-")[1])
    # convert coordinates to 0-based half opened coordinates
    region_start = region_start -1
    
    # load json with count of umi families per position and umi group
    try:
        infile = open(umifile)
        umi_families = json.load(infile)
        infile.close()
    except:
        raise ValueError("ERR: Unable to load umi json file")
        
    print(GetCurrentTime() + 'Generating consensus...')

    # get percent threshold 
    consensus_threshold = GetThresholds(config, 'percent_consensus_threshold', percentthreshold)
    # get count threshold
    count_threshold = GetThresholds(config, 'count_consensus_threshold', countthreshold)
    # get umi position threshold 
    pos_threshold = GetThresholds(config, 'umi_family_pos_threshold', postthreshold)
    
    # get comma-separated list of minimum family sizes 
    fam_size = GetFamSize(config, famsize)
    
    # write consensus output file
    ConsDir = os.path.join(outdir, 'Consfiles')
    generate_consensus_output(contig, region_start, region_end, bam_file, umi_families, ConsDir, fam_size, pos_threshold, consensus_threshold, count_threshold, separator, base_quality_score, max_depth=maxdepth, truncate=truncate, ignore_orphans=ignoreorphans, stepper=stepper)
 
    print(GetCurrentTime() + 'Consensus generated. Consensus file written to {0}.'.format(ConsDir))


def VCF_converter(config, outdir, reference, refthreshold, altthreshold, filterthreshold, famsize, consFiles):
    '''
    (str, str, str, float, float, int, int, list) --> None

    Converts consensus files into VCF format for a given family size

    Parameters
    ----------
    - config (str): Path to the config file
    - outdir (str): Output directory where subdirectories are created
    - reference (str): Path to the reference genome 
    - refthreshold (float): Maximum reference frequency (in %) to consider alternative variants
                          (ie. position with ref freq <= ref_threshold is considered variable)
    - altthreshold (float): Minimum allele frequency (in %) to consider an alternative allele at a variable position 
                          (ie. allele freq >= alt_threshold and ref freq <= ref_threshold --> record alternative allele)
    - filterthreshold (int): Minimum number of reads to pass alternative variants 
                             (ie. filter = PASS if variant depth >= alt_threshold)
    - famsize (int): Minimum UMI family size
    - consFiles (list): List of consensus files generated during UMI collapse 
    '''

    # get output directory from the config or command. set to current dir if not provided
    outdir = GetOutputDir(config, outdir)
    # create outputdir if doesn't exist
    if os.path.isfile(outdir) == True:
        raise ValueError('ERR: Output directory cannot be a file')
    else:
        os.makedirs(outdir, exist_ok=True)
    
    # make a list of consensus files
    # remove empty files in place and print a warning
    DropEmptyFiles(consFiles)
        
    # create vcf dir if doesn't exist already
    VCFDir = os.path.join(outdir, 'VCFfiles')
    os.makedirs(VCFDir, exist_ok=True)
       
    # get reference threshold
    ref_threshold = GetThresholds(config, 'percent_ref_threshold', refthreshold)
    # get allele threshold
    alt_threshold = GetThresholds(config, 'percent_alt_threshold', altthreshold)
    # get filter threshold
    filter_threshold = GetThresholds(config, 'filter_threshold', filterthreshold)
    
    print(GetCurrentTime() + 'Generating VCFs...')

    # loop over consensus files
    for filename in consFiles:
        # write a VCF per consensus file for fiven family size
        if os.path.basename(filename)[:-5] == '.cons':
            outputfile = os.path.join(VCFDir, os.path.basename(filename)[:-5] + '_famsize_{0}.vcf'.format(famsize))
        else:
            outputfile = os.path.join(VCFDir, os.path.basename(filename) + '_famsize_{0}.vcf'.format(famsize))
        WriteVCF(filename, outputfile, reference, ref_threshold, alt_threshold, filter_threshold, famsize)

    print(GetCurrentTime() + 'VCFs generated. VCF files written to {0}'.format(VCFDir))


def merge_files(outdir, files, datatype):
    '''
    (str, list, str) -> None
    
    Merge all files of given datatype in a single file
        
    Parameters
    ----------
    - outdir (str): Output directory where subdirectories are created
    - files (list): List of files to merge. Must be of the same datatatype
    - datatype (str): Type of files to be merged.
                     Valid options are 'datafiles', 'consensusfiles', 'umifiles'
    '''
    
    # create outputdir if doesn't exist
    os.makedirs(outdir, exist_ok=True)
    # create vcf dir if doesn't exist already
    subdirs = []
    for i in ['VCFfiles', 'Consfiles', 'Datafiles', 'Consfiles']:
        os.makedirs(os.path.join(outdir, i), exist_ok=True)    
        subdirs.append(os.path.join(outdir, i))
    VCFDir, ConsDir, DataDir, UmiDir = subdirs 
    
    # check which files need to be merged
    if datatype == 'datafiles':
        content = MergeDataFiles(files)
        filename = os.path.join(DataDir, 'Merged_DataFile.csv')
    elif datatype == 'consensusfiles':
        content = MergeConsensusFiles(files)
        filename = os.path.join(ConsDir, 'Merged_ConsensusFile.cons')
    elif datatype == 'umifiles':
        content = MergeUmiFiles(files)
        filename = os.path.join(UmiDir, 'Merged_UmiFile.json')
        
    newfile = open(filename, 'w')
    if datatype in ['datafiles', 'consensusfiles']:
        newfile.write(content)
    elif datatype == 'umifiles':
        json.dump(content, newfile, sort_keys = True, indent=4)
    newfile.close()
        

def run_scripts(outdir, config, bamfile, reference, famsize, bedfile, countthreshold,
                consensusthreshold, postthreshold, distthreshold, refthreshold,
                altthreshold, filterthreshold, percentthreshold, maxdepth, truncate, ignoreorphans,
                ignore, stepper, merge, plot, report, call, mincov, minratio, minumis,
                minchildren, extension, sample, mem, mypython, mydebarcer, project,
                separator, base_quality_score, readcount):
    '''
    (str, str, str, str, str, str, int, float, int, int, float, float, int, int, float,
    bool, bool, bool, str, bool, bool, bool, bool, int, float, int, int, str, str | None,
    str, str, str, str, str, int, int) -> None
       
    Submits jobs to run Umi Grouping, Collapsing and Plotting and Reporting if activated
    
    Parameters
    ----------
    - outdir (str): Output directory where subdirectories are created
    - config (str): Path to the config file
    - bamfile (str): Path to the BAM file
    - reference (str): Path to the refeence genome
    - famsize (str): Comma-separated list of minimum umi family size to collapase on
    - bedfile (str): Path to the bed file
    - countthreshold (int): Base count threshold in pileup column
    - consensus_threshold (float): Majority rule consensus threshold in pileup column
    - postthreshold (int): Umi position threshold for grouping umis together
    - distthreshold (int): Hamming distance threshold for connecting parent-children umis
    - refthreshold (float): Maximum reference frequency (in %) to consider alternative variants
                          (ie. position with ref freq <= ref_threshold is considered variable)
    - altthreshold (float): Minimum allele frequency (in %) to consider an alternative allele at a variable position 
                          (ie. allele freq >= alt_threshold and ref freq <= ref_threshold --> record alternative allele)
    - filterthreshold (int): Minimum number of reads to pass alternative variants 
                             (ie. filter = PASS if variant depth >= alt_threshold)
    - percentthreshold (float): Base percent threshold in pileup column
    - maxdepth (int): Maximum read depth. Default is 1000000
    - truncate (bool): If truncate is True and a region is given, only pileup columns
                     in the exact region specificied are returned. Default is False
    - ignoreorphans (bool): Ignore orphans (paired reads that are not in a proper pair). Default is True'
    - ignore (bool): Keep the most abundant family and ignore families at other positions within each group. Default is False
    - stepper (str): Controls how the iterator advances. Accepeted values:
                    'all': skip reads with following flags: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
                    'nofilter': uses every single read turning off any filtering
    - merge (bool): Merge data, json and consensus files respectively into a 1 single file. Default is True
    - plot (bool): Generate figure plots if True
    - report (bool): Generate analysis report if True
    - call (bool): Convert consensus files to VCF if True
    - mincov (int): Minimum read depth to label regions
    - minratio (float): Minimum ratio to label regions    
    - minumis (int): Minimum number of umis to label regions
    - minchildren (int): Minimum number of umi children to label regions
    - extension (str): Figure file extension
    - sample (str or None): Sample name to appear in report. If empty str, outdir basename is used
    - mem (str): Requested memory for submiiting jobs to SGE. Default is 10g
    - mypython (str): Path to python. Default is: /.mounts/labs/PDE/Modules/sw/python/Python-3.6.4/bin/python3.6
    - mydebarcer (str): Path to the file debarcer.py. Default is /.mounts/labs/PDE/Modules/sw/python/Python-3.6.4/lib/python3.6/site-packages/debarcer/debarcer.py
    - project (str): Project name to submit jobs on univa
    - separator (str): String separating the UMI from the remaining of the read name
    - base_quality_score (int): Base quality score threshold. No offset of 33 needs to be subtracted
    - readcount (int): Minimum number of reads in region required for grouping. Default is 0  
    '''
    
    # get bam file from config or command
    bamfile = GetInputFiles(config, bamfile, 'bam_file')
    # get reference
    reference = GetInputFiles(config, reference, 'reference_file')
       
    # get output directory from the config or command. set to current dir if not provided
    outdir = GetOutputDir(config, outdir)
    # create outputdir
    if os.path.isfile(outdir) == True:
        raise ValueError('ERR: Output directory cannot be a file')
    else:
        os.makedirs(outdir, exist_ok=True)
    # create subdirectoy structure
    CreateDirTree(outdir)

    # get comma-separated list of minimum family size
    famsize = GetFamSize(config, famsize)
    
    # get thresholds from command or config
    count_threshold = GetThresholds(config, 'count_consensus_threshold', countthreshold)
    consensus_threshold = GetThresholds(config, 'percent_consensus_threshold', percentthreshold)
    dist_threshold = GetThresholds(config, 'umi_edit_distance_threshold', distthreshold)
    post_threshold = GetThresholds(config, 'umi_family_pos_threshold', postthreshold)
    ref_threshold = GetThresholds(config, 'percent_ref_threshold', refthreshold)
    alt_threshold = GetThresholds(config, 'percent_alt_threshold', altthreshold)
    filter_threshold = GetThresholds(config, 'filter_threshold', filterthreshold)
    
    # create shell scripts and run qsubs to Group and Collapse umis 
    submit_jobs(bamfile, outdir, reference, famsize, bedfile, count_threshold,
                consensus_threshold, dist_threshold, post_threshold, ref_threshold,
                alt_threshold, filter_threshold, maxdepth, truncate, ignoreorphans,
                ignore, stepper, merge, plot, report, call, mincov, minratio, minumis,
                minchildren, extension, sample, mydebarcer, mypython, mem, project, separator,
                base_quality_score, readcount)
  
    
def generate_figures(directory, config, extension, report, sample, min_cov, min_ratio, min_umis, min_children, ref_threshold):
    '''
    (str, str, str, bool, str | None, int, float, int, int, int) -> None
    
    Generate plots in Figures directory
    
    Parameters
    ----------
    - directory (str): Directory with subdirectories ConsFiles and Datafiles 
    - config (str): Path to the config file
    - extension (str): Figure format. Accepted values: png, pdf, jpeg, tiff
    - report (bool): Boolean, generate a report if True
    - sample (str or None): Optional parameter, sample name to appear in report
    - min_cov (int): Minimum read depth to label regions
    - min_ratio (float): Minimum ratio to label regions    
    - min_umis (int): Minimum number of umis to label regions
    - min_children (int): Minimum number of umi children to label regions
    - ref_threshold (int): Cut Y axis at 100 - ref_threshold
    '''
    
    # get the reference threshold to consider variable positions
    # cut Y axis at non-ref-freq
    ref_threshold = GetThresholds(config, 'percent_ref_threshold', ref_threshold)
    non_ref_freq = 100 - ref_threshold
        
    # get subdirectories
    L = ['Consfiles', 'Umifiles', 'Stats', 'Datafiles']
    T = [os.path.join(directory, i) for i in L]
    for i in T:
        if os.path.isdir(i) == False:
            raise ValueError('ERR: Expecting directory {0}'.format(i))
      
    # unpack directories
    ConsDir, UmiDir, StatsDir, DataDir =  T 
        
    # create directory to save figures if it doesn't exist
    FigDir = os.path.join(directory, 'Figures')
    os.makedirs(FigDir, exist_ok=True)
        
    # make a list of consensus files that are not empty
    ConsFiles = [os.path.join(ConsDir, i) for i in os.listdir(ConsDir) if i.startswith('chr') and i[-5:] == '.cons']
    # remove empty files in place and print a warning
    DropEmptyFiles(ConsFiles)
    # make a list of datafiles with umis that are not empty
    DataFiles = [os.path.join(DataDir, i) for i in os.listdir(DataDir) if i.startswith('datafile') and 'chr' in i and i[-4:] == '.csv']
    # remove empty files in place and print a warning
    DropEmptyFiles(DataFiles)
    # make a list of umi files
    UmiFiles = [os.path.join(UmiDir, i) for i in os.listdir(UmiDir) if i.startswith('chr') and i[-5:] == '.json']
    # make a list of files with individual umi information
    UmiInfoFiles = [os.path.join(StatsDir, i) for i in os.listdir(StatsDir) if i.startswith('Umis_') and i[-5:] == '.json' and 'before_grouping' in i]
    # make a list of json files with mapping reads information
    MappingInfo = [os.path.join(StatsDir, i) for i in os.listdir(StatsDir) if i.startswith('Mapped') and i[-5:] == '.json']

    # get the file with coverage stats
    CovStats = os.path.join(StatsDir, 'CoverageStats.yml')
    
    # check that paths to files are valid. raise ValueError if file invalid
    CheckFilePath(ConsFiles + DataFiles + UmiFiles + UmiInfoFiles)
    
    # make a list of colors. each color is used for plotting data for a given family size
    Colors = ['black', '#4B0082', '#7B68EE', '#c3baf7', '#8A2BE2',
              '#b54dff', '#BA55D3', '#ce85e0', '#DDA0DD', '#f8ecf8',
              '#066024', '#0ba83f', '#0ff05a', '#57f48c', '#9ff9bd',
              '#b35900', '#e67300', '#ff8c1a', '#ffa64d', '#ffbf80',
              '#b30000', '#e60000', '#ff1a1a', '#ff4d4d', '#ff8080']
    
    # plot proportions of correct and incorrect reads seen during pre-processing
    Inputfile = os.path.join(StatsDir, 'Processing_Read_Info.json')
    CheckFilePath([Inputfile])
    Outputfile = os.path.join(FigDir, 'Proportion_correct_reads.' + extension)
    PlotIncorrectReads(Inputfile, Outputfile, 'preprocessing', 6, 6)
        
    # plot UMI occurence resulting from pre-processing
    Inputfile = os.path.join(StatsDir, 'Umi_counts.txt')
    CheckFilePath([Inputfile])
    # get umi occurence
    umi_occurence = GetUmiCountFromPreprocessing(Inputfile)
    Outputfile = os.path.join(FigDir, 'UMI_occurence_preprocessing.' + extension)
    PlotUMiFrequency(umi_occurence, Outputfile, 'UMI distribution after pre-processing', False, 8, 5)
    
    # plot coverage
    PlotDataPerRegion(CovStats, DataFiles, outputfile=os.path.join(FigDir, 'Coverage_Umi_Count'), mincov=min_cov, datatype='coverage')

    # plot graphs for each consensus file
    for filename in ConsFiles:
        # plot mean family size for each consensus file/region
        region = FormatRegion(filename).replace(':', '-')
        Outputfile = os.path.join(FigDir, 'MeanFamilySize_{0}.{1}'.format(region, extension))
        PlotMeanFamSize(filename, Colors[1:], Outputfile, 9, 6)
            
        # plot non-reference frequency
        Outputfile = os.path.join(FigDir, 'NonRefFreq_{0}.{1}'.format(region, extension))
        PlotNonRefFreqData(filename, Colors, Outputfile, 8, 10, ylabel='Non-reference allele frequency')
    
        # plot non-reference frequency limiting Y axis to 20% for visualization of low-frequency variants 
        Outputfile = os.path.join(FigDir, 'NonRefFreq_low_freq_{0}.{1}'.format(region, extension))
        PlotNonRefFreqData(filename, Colors, Outputfile, 8, 10, YLimit=non_ref_freq, title='Y axis cut at {0}%'.format(non_ref_freq), legend='legend')
        
        # plot raw and consensus depth
        Outputfile = os.path.join(FigDir, 'RawConsensusDepth_{0}.{1}'.format(region, extension))    
        PlotConsDepth(filename, Colors, Outputfile, 9, 6)
       
    # plot network and network degree for each umi file/region
    for filename in UmiFiles:
        # get region from file name
        region = os.path.basename(filename)
        region = region[:-5]
        region = '-'.join(list(map(lambda x: x.strip(), region.split(':'))))
        
        # plot network and degree
        Outputfile = os.path.join(FigDir, 'UMI_network_degree_{0}.{1}'.format(region, extension))        
        PlotNetworkDegree(filename, Outputfile, 9, 6)
        
        # plot marginal distributions of UMI family size and read depth
        Outputfile = os.path.join(FigDir, 'UMI_size_depth_marginal_distribution_{0}.{1}'.format(region, extension))
        PlotFamSizeReadDepth(filename, Outputfile)
        
        # plot distribution of read depth for each umi families
        Outputfile = os.path.join(FigDir, 'Read_depth_per_umi_family_{0}.{1}'.format(region, extension))
        PlotReadDepth(filename, Outputfile, 10, 6)

    # plot umi frequency for individual umis before grouping
    for filename in UmiInfoFiles:
        region = os.path.basename(filename)
        region = region[region.index('chr'): region.index('_before')].replace(':', '-')
        # get parent+children and parent only counts
        all_umis, parent_umis = GetIndividualUmiInfo(filename)
        Outputfile = os.path.join(FigDir, 'UMI_freq_distribution_{0}.{1}'.format(region, extension)) 
        PlotUMiFrequency([all_umis, parent_umis], Outputfile, 'UMI distribution before grouping', True, 9, 6)
    
    # plot proportion of mapped/unmapped reads
    for filename in MappingInfo:
        region = os.path.basename(filename)
        region = region[region.rindex('_')+1:-5].replace(':', '-')
        Outputfile = os.path.join(FigDir, 'Proportion_unmapped_reads_{0}.{1}'.format(region, extension))
        PlotIncorrectReads(filename, Outputfile, 'mapping', 5, 5)
        
    # plot children to parent umi count ratio
    PlotDataPerRegion(CovStats, DataFiles, outputfile=os.path.join(FigDir, 'Child_Parent_Umis_Ratio'), minval=min_ratio, datatype='ratio')
    
    # plot total umi counts
    PlotDataPerRegion(CovStats, DataFiles, outputfile=os.path.join(FigDir, 'Total_Umis'), minval=min_umis, datatype='umis')

    # plot children umi counts
    PlotDataPerRegion(CovStats, DataFiles, outputfile=os.path.join(FigDir, 'Children_Umis'), minval=min_children, datatype='children')

    # plot children vs parent umis for each interval
    PlotParentsToChildrenCounts(DataFiles, os.path.join(FigDir, 'PTU_vs_CTU.' + extension), 9, 6)

    # plot parent frequencies vs children UMI counts
    PlotParentFreq(DataFiles, os.path.join(FigDir, 'Children_vs_ParentFreq.' + extension), 7, 4)
    
    # check if reporting
    if report == True:
        if extension != 'pdf':
            # create subdirectory
            ReportDir = os.path.join(directory, 'Report')
            os.makedirs(ReportDir, exist_ok=True)
            report = os.path.join(ReportDir, 'debarcer_report.html')
            WriteReport(directory, extension, report, min_cov, min_ratio, min_umis, min_children, renderer=mistune.Markdown(), sample=sample)
            

def report(directory, extension, min_cov, min_ratio, min_umis, min_children, sample=None):
    '''
    (str,  str, int, float, int, int, str | None) -> None

    Write an html report of debarcer analysis for a given sample

    Parameters
    ----------
    - directory (str): Directory with subfolders including Figures 
    - extension (str): Extension of the figure files
    - min_cov (int): Minimum coverage value. Values below are plotted in red
    - min_ratio (float) Minimum children to parent umi ratio. Values below are plotted in red
    - min_umis (int): Minimum umi count. Values below are plotted in red
    - min_children (int): Minimum children umi count. Values below are plotted in red
    - sample (str): Sample name 
    '''
    
    # create subdirectory
    ReportDir = os.path.join(directory, 'Report')
    os.makedirs(ReportDir, exist_ok=True)
    report = os.path.join(ReportDir, 'debarcer_report.html')
    
    # check that subdirectories with required files exist
    T = [os.path.join(directory, i) for i in ['Stats', 'Datafiles', 'Figures']]
    for i in T:
        if os.path.isdir(i) == False:
            raise ValueError('ERR: Expecting directory {0}'.format(i))
    WriteReport(directory, extension, report, min_cov, min_ratio, min_umis, min_children, renderer=mistune.Markdown(), sample=sample)
    

def generate_bed(bamfile, outputfile, contig, min_cov, region_size, max_depth, ignore_orphans, stepper):
    '''
    (str, str, str, int, int, int, bool, str) -> None
    
    Write a bed file (1-based coordinates) with all genomic intervals of minimum
    length region_size for which all positions have read depth equals to min_cov or greater
        
    Precondition: bamfile is coordinate-sorted and has 'SQ' fields    
    
    Parameters
    ----------
    - bamfile (str): Path to the bam file
    - outputfile (str): Path to the output bed file
    - contig (str): Chromosome name, eg. chrN
    - min_cov (int): Minimum read depth for all positions in genomic interval
    - region_size (int): Minimum length of the genomic interval    
    - max_depth (int): Maximum read depth
    - ignore_orphans (bool): Ignore orphan reads (paired reads not in proper pair) if True
    - stepper (str): Controls how the iterator advances. Accepeted values:
                    'all': skip reads with following flags: BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP
                    'nofilter': uses every single read turning off any filtering    
    '''
    
    WriteTargetsBed(bamfile, outputfile, min_cov, region_size, max_depth, ignore_orphans, stepper)
    

def main():
    
    ## Argument + config parsing and error handling
    parser = argparse.ArgumentParser(prog='debarcer.py', description="A package for De-Barcoding\
                                     and Error Correction of sequencing data containing molecular barcodes")
    subparsers = parser.add_subparsers(help='sub-command help', dest='subparser_name')
    
    ## Preprocess command
    p_parser = subparsers.add_parser('preprocess', help="Preprocess mode for processing fastq files")
    p_parser.add_argument('-o', '--OutDir', dest='outdir', help='Output directory. Available from command or config')
    p_parser.add_argument('-r1', '--Read1', dest='read1', help='Path to first FASTQ file.', required=True)
    p_parser.add_argument('-r2', '--Read2', dest='read2', help='Path to second FASTQ file, if applicable')
    p_parser.add_argument('-r3', '--Read3', dest='read3', help='Path to third FASTQ file, if applicable')
    p_parser.add_argument('-p', '--Prepname', dest='prepname', choices=['HALOPLEX', 'SURESELECT', 'EPIC-DS', 'SIMSENSEQ-PE', 'SIMSENSEQ-SE'], 
                          help='Name of library prep to  use (defined in library_prep_types.ini)', required=True)
    p_parser.add_argument('-pf', '--Prepfile', dest='prepfile', help='Path to your library_prep_types.ini file')
    p_parser.add_argument('-c', '--Config', dest='config', help='Path to your config file')
    p_parser.add_argument('-px', '--Prefix', dest= 'prefix', help='Prefix for naming umi-reheradered fastqs. Use Prefix from Read1 if not provided') 
        
    ## Bed command
    b_parser = subparsers.add_parser('bed', help='Generate a bed file by scanning input bam for regions of coverage')
    b_parser.add_argument('-b', '--Bamfile', dest='bamfile', help='Path to the BAM file', required=True)
    b_parser.add_argument('-bd', '--Bedfile', dest='bed', help='Path to the output bed file', required=True)
    b_parser.add_argument('-mv', '--MinCov', dest='mincov', type=int, help='Minimum read depth value at all positions in genomic interval', required=True)
    b_parser.add_argument('-r', '--RegionSize', dest='regionsize', type=int, help='Minimum length of the genomic interval (in bp)', required=True)
    b_parser.add_argument('-m', '--MaxDepth', dest='maxdepth', default=1000000, type=int, help='Maximum read depth. Default is 1000000')
    b_parser.add_argument('-io', '--IgnoreOrphans', dest='ignoreorphans', action='store_true', help='Ignore orphans (paired reads that are not in a proper pair). Default is False, becomes True if used')
    b_parser.add_argument('-stp', '--Stepper', dest='stepper', choices=['all', 'nofilter'], default='nofilter',
                          help='Filter or include reads in the pileup. Options all: skip reads with BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP flags,\
                          nofilter: uses every single read turning off any filtering')
        
    ## UMI group command
    g_parser = subparsers.add_parser('group', help="Groups UMIs into families.")
    g_parser.add_argument('-o', '--Outdir', dest='outdir', help='Output directory where subdirectories are created')
    g_parser.add_argument('-r', '--Region', dest='region', help='Region coordinates to search for UMIs. chrN:posA-posB. posA and posB are 1-based included', required=True)
    g_parser.add_argument('-b', '--Bamfile', dest='bamfile', help='Path to the BAM file')
    g_parser.add_argument('-c', '--Config', dest='config', help='Path to the config file')
    g_parser.add_argument('-d', '--Distance', dest='distthreshold', type=int, help='Hamming distance threshold for connecting parent-children umis')
    g_parser.add_argument('-p', '--Position', dest='postthreshold', type=int, help='Umi position threshold for grouping umis together')
    g_parser.add_argument('-i', '--Ignore', dest='ignore', choices=[True, False], type=ConvertArgToBool, default=False, help='Keep the most abundant family and ignore families at other positions within each group. Default is False')
    g_parser.add_argument('-t', '--Truncate', dest='truncate', choices=[True, False], default=False, type=ConvertArgToBool, help='Discard reads overlapping with the genomic region if True. Default is False')
    g_parser.add_argument('-s', '--Separator', dest='separator', default=':', help = 'String separating the UMI from the remaining of the read name')
    g_parser.add_argument('-rc', '--ReadCount', dest='readcount', default=0, type=int, help = 'Minimum number of reads in region required for grouping. Default is 0')
        
    ## Base collapse command
    c_parser = subparsers.add_parser('collapse', help="Base collapsing from given UMI families file.")
    c_parser.add_argument('-c', '--Config', dest='config', help='Path to the config file')
    c_parser.add_argument('-o', '--Outdir', dest='outdir', help='Output directory where subdirectories are created')
    c_parser.add_argument('-b', '--Bamfile', dest='bamfile', help='Path to the BAM file')
    c_parser.add_argument('-rf', '--Reference', dest='reference', help='Path to the refeence genome')
    c_parser.add_argument('-r', '--Region', dest='region', help='Region coordinates to search for UMIs. chrN:posA-posB. posA and posB are 1-based included', required=True)
    c_parser.add_argument('-u', '--Umi', dest='umifile', help='Path to the .umis file', required=True)
    c_parser.add_argument('-f', '--Famsize', dest='famsize', help='Comma-separated list of minimum umi family size to collapase on')
    c_parser.add_argument('-ct', '--CountThreshold', dest='countthreshold', type=int, help='Base count threshold in pileup column')
    c_parser.add_argument('-pt', '--PercentThreshold', dest='percentthreshold', type=float, help='Majority rule consensus threshold in pileup column')
    c_parser.add_argument('-p', '--Position', dest='postthreshold', type=int, help='Umi position threshold for grouping umis together')
    c_parser.add_argument('-m', '--MaxDepth', dest='maxdepth', default=1000000, type=int, help='Maximum read depth. Default is 1000000')
    c_parser.add_argument('-t', '--Truncate', dest='truncate', choices=[True, False], default=False, type=ConvertArgToBool, help='If truncate is True and a region is given,\
                          only pileup columns in the exact region specificied are returned. Default is False')
    c_parser.add_argument('-i', '--IgnoreOrphans', dest='ignoreorphans', choices=[True, False], default=False, type=ConvertArgToBool, help='Ignore orphans (paired reads that are not in a proper pair). Default is False')
    c_parser.add_argument('-stp', '--Stepper', dest='stepper', choices=['all', 'nofilter'], default='nofilter',
                          help='Filter or include reads in the pileup. Options all: skip reads with BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP flags,\
                          nofilter: uses every single read turning off any filtering')
    c_parser.add_argument('-s', '--Separator', dest='separator', default=':', help = 'String separating the UMI from the remaining of the read name')
    c_parser.add_argument('-bq', '--Quality', dest='base_quality_score', type=int, default=25, help = 'Base quality score threshold. Bases with quality scores below the threshold are not used in the consensus. Default is 25')
    
    ## Variant call command - requires cons file (can only run after collapse)
    v_parser = subparsers.add_parser('call', help="Convert consensus file into VCF format.")
    v_parser.add_argument('-o', '--Outdir', dest='outdir', help='Output directory where subdirectories are created')
    v_parser.add_argument('-c', '--Config', dest='config', help='Path to the config file')
    v_parser.add_argument('-rf', '--Reference', dest='reference', help='Path to the refeence genome')
    v_parser.add_argument('-rt', '--RefThreshold', dest='refthreshold', default=95, type=float, 
                          help='Maximum reference frequency to consider (in percent) alternative variants\
                          (ie. position with ref freq <= ref_threshold is considered variable)')
    v_parser.add_argument('-at', '--AlternativeThreshold', dest='altthreshold', type=float, default=2,
                          help='Minimum allele frequency (in percent) to consider an alternative allele at a variable position\
                          (ie. allele freq >= alt_threshold and ref freq <= ref_threshold: alternative allele)')
    v_parser.add_argument('-ft', '--FilterThreshold', dest='filterthreshold', type=int, default=10,
                          help='Minimum number of reads to pass alternative variants\
                          (ie. filter = PASS if variant depth >= alt_threshold)')
    v_parser.add_argument('-f', '--Famsize', dest='famsize', type=int, help='Minimum UMI family size', required=True)
    v_parser.add_argument('-cf', '--Consfiles', dest='consfiles', nargs = '*', help='List of consensus files', required=True)
    
    ## Run scripts command 
    r_parser = subparsers.add_parser('run', help="Generate scripts for umi grouping, collapsing and VCF formatting for target regions specified by the BED file.")
    r_parser.add_argument('-o', '--Outdir', dest='outdir', help='Output directory where subdirectories are created')
    r_parser.add_argument('-c', '--Config', dest='config', help='Path to the config file')
    r_parser.add_argument('-b', '--Bamfile', dest='bamfile', help='Path to the BAM file')
    r_parser.add_argument('-rf', '--Reference', dest='reference', help='Path to the refeence genome')
    r_parser.add_argument('-f', '--Famsize', dest='famsize', help='Comma-separated list of minimum umi family size to collapase on')
    r_parser.add_argument('-bd', '--Bedfile', dest='bedfile', help='Path to the bed file', required=True)
    r_parser.add_argument('-ct', '--CountThreshold', dest='countthreshold', type=int, help='Base count threshold in pileup column')
    r_parser.add_argument('-pt', '--PercentThreshold', dest='percentthreshold', type=float, help='Base percent threshold in pileup column')
    r_parser.add_argument('-p', '--Position', dest='postthreshold', type=int, help='Umi position threshold for grouping umis together')
    r_parser.add_argument('-d', '--Distance', dest='distthreshold', type=int, help='Hamming distance threshold for connecting parent-children umis')
    r_parser.add_argument('-rt', '--RefThreshold', dest='refthreshold', default=95, type=float, help='A position is considered variable of reference frequency is <= ref_threshold')
    r_parser.add_argument('-at', '--AlternativeThreshold', dest='altthreshold', default=2, type=float, help='Variable position is labeled PASS if allele frequency >= alt_threshold')
    r_parser.add_argument('-ft', '--FilterThreshold', dest='filterthreshold', default=10, type=int, help='Minimum number of reads to pass alternative variants')
    r_parser.add_argument('-m', '--MaxDepth', dest='maxdepth', default=1000000, type=int, help='Maximum read depth. Default is 1000000')
    r_parser.add_argument('-t', '--Truncate', dest='truncate', action='store_true', help='Only pileup columns in the exact region specificied are returned. Default is False, becomes True is used')
    r_parser.add_argument('-io', '--IgnoreOrphans', dest='ignoreorphans', action='store_true', help='Ignore orphans (paired reads that are not in a proper pair). Default is False, becomes True if used')
    r_parser.add_argument('-i', '--Ignore', dest='ignore', action='store_true', help='Keep the most abundant family and ignore families at other positions within each group. Default is False, becomes True if used')
    r_parser.add_argument('-mg', '--Merge', dest='merge', action='store_false', help='Merge data, json and consensus files respectively into a 1 single file. Default is True, becomes False if used')
    r_parser.add_argument('-pl', '--Plot', dest='plot',  action='store_false', help='Generate figure plots. Default is True, becomes False if used')
    r_parser.add_argument('-rp', '--Report', dest='report', action='store_false', help='Generate report. Default is True, becomes False if used')
    r_parser.add_argument('-cl', '--Call', dest='call', action='store_false', help='Convert consensus files to VCF format. Default is True, becomes False if used')
    r_parser.add_argument('-ex', '--Extension', dest='extension', choices=['png', 'jpeg', 'pdf'], default='png', help='Figure format. Does not generate a report if pdf, even with -r True. Default is png')
    r_parser.add_argument('-sp', '--Sample', dest='sample', help='Sample name to appear to report. Optional, use Output directory basename if not provided')
    r_parser.add_argument('-pr', '--Project', dest='project', default='gsi', help='Project for submitting jobs on Univa')
    r_parser.add_argument('-mm', '--Memory', dest='mem', default=20, type=int, help='Requested memory for submitting jobs to SGE. Default is 20g')
    r_parser.add_argument('-py', '--MyPython', dest='mypython', default='/.mounts/labs/PDE/Modules/sw/python/Python-3.6.4/bin/python3.6',
                          help='Path to python. Default is /.mounts/labs/PDE/Modules/sw/python/Python-3.6.4/bin/python3.6')
    r_parser.add_argument('-db', '--MyDebarcer', dest='mydebarcer', default='/.mounts/labs/PDE/Modules/sw/python/Python-3.6.4/lib/python3.6/site-packages/debarcer/debarcer.py',
                          help='Path to the file debarcer.py. Default is /.mounts/labs/PDE/Modules/sw/python/Python-3.6.4/lib/python3.6/site-packages/debarcer/debarcer.py')
    r_parser.add_argument('-mv', '--MinCov', dest='mincov', type=int, default=1000, help='Minimum coverage value. Values below are plotted in red')
    r_parser.add_argument('-mr', '--MinRatio', dest='minratio', type=float, default=0.1, help='Minimum children to parent umi ratio. Values below are plotted in red')
    r_parser.add_argument('-mu', '--MinUmis', dest='minumis', type=int, default=1000, help='Minimum umi count. Values below are plotted in red')
    r_parser.add_argument('-mc', '--MinChildren', dest='minchildren', type=int, default=500, help='Minimum children umi count. Values below are plotted in red')
    r_parser.add_argument('-stp', '--Stepper', dest='stepper', choices=['all', 'nofilter'], default='nofilter',
                          help='Filter or include reads in the pileup. Options all: skip reads with BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP flags,\
                          nofilter: uses every single read turning off any filtering')
    r_parser.add_argument('-s', '--Separator', dest='separator', default=':', help = 'String separating the UMI from the remaining of the read name')
    r_parser.add_argument('-bq', '--Quality', dest='base_quality_score', type=int, default=25, help = 'Base quality score threshold. Bases with quality scores below the threshold are not used in the consensus. Default is 25')
    r_parser.add_argument('-rc', '--ReadCount', dest='readcount', default=0, type=int, help = 'Minimum number of reads in region required for grouping. Default is 0')
        
    ## Merge files command 
    m_parser = subparsers.add_parser('merge', help="Merge files from each region into a single file")
    m_parser.add_argument('-o', '--Outdir', dest='outdir', help='Output directory where subdirectories are created')
    m_parser.add_argument('-f', '--Files', dest='files', nargs='*', help='List of files to be merged', required=True)
    m_parser.add_argument('-dt', '--DataType', dest='datatype', choices=['datafiles', 'consensusfiles', 'umifiles'], help='Type of files to be merged', required=True)
        
    ## Generate graphs	
    plot_parser = subparsers.add_parser('plot', help="Generate graphs for umi and cons data files", add_help=True)
    plot_parser.add_argument('-c', '--Config', dest='config', help='Path to the config file')
    plot_parser.add_argument('-d', '--Directory', dest='directory', help='Directory with subdirectories ConsFiles and Datafiles', required=True)
    plot_parser.add_argument('-e', '--Extension', dest='extension', default='png', choices=['pdf', 'png', 'jpeg'], help='Figure format. Does not generate a report if pdf, even with -r True. Default is png')
    plot_parser.add_argument('-s', '--Sample', dest='sample', help='Sample name to apear in the report is reporting flag activated. Optional')
    plot_parser.add_argument('-r', '--Report', dest='report', choices=[False, True], type=ConvertArgToBool, default=True, help='Generate a report if activated. Default is True')
    plot_parser.add_argument('-mv', '--MinCov', dest='mincov', type=int, default=1000, help='Minimum coverage value. Values below are plotted in red')
    plot_parser.add_argument('-mr', '--MinRatio', dest='minratio', type=float, default=0.1, help='Minimum children to parent umi ratio. Values below are plotted in red')
    plot_parser.add_argument('-mu', '--MinUmis', dest='minumis', type=int, default=1000, help='Minimum umi count. Values below are plotted in red')
    plot_parser.add_argument('-mc', '--MinChildren', dest='minchildren', type=int, default=500, help='Minimum children umi count. Values below are plotted in red')
    plot_parser.add_argument('-rt', '--RefThreshold', dest='refthreshold', default=95, type=float, help='Cut Y axis at non-ref frequency, the minimum frequency to consider a position variable')
        
    ## Generate report
    report_parser = subparsers.add_parser('report', help="Generate report", add_help=True)
    report_parser.add_argument('-d', '--Directory', dest='directory', help='Directory with subdirectories including Figures', required=True)
    report_parser.add_argument('-e', '--Extension', dest='extension', default='png', choices=['pdf', 'png', 'jpeg'], help='Figure format. Does not generate a report if pdf, even with -r True. Default is png')
    report_parser.add_argument('-s', '--Sample', dest='sample', help='Sample name. Optional. Directory basename is sample name if not provided')
    report_parser.add_argument('-mv', '--MinCov', dest='mincov', type=int, default=1000, help='Minimum coverage value. Values below are plotted in red')
    report_parser.add_argument('-mr', '--MinRatio', dest='minratio', type=float, default=0.1, help='Minimum children to parent umi ratio. Values below are plotted in red')
    report_parser.add_argument('-mu', '--MinUmis', dest='minumis', type=int, default=1000, help='Minimum umi count. Values below are plotted in red')
    report_parser.add_argument('-mc', '--MinChildren', dest='minchildren', type=int, default=500, help='Minimum children umi count. Values below are plotted in red')
    
        
    args = parser.parse_args()
    
    
    if args.subparser_name == 'preprocess':
        try:
            preprocess_reads(args.outdir, args.read1, args.read2, args.read3, args.prepname, args.prepfile, args.config, args.prefix)
        except AttributeError as e:
            print('AttributeError: {0}\n'.format(e))
            print(parser.format_help())
    elif args.subparser_name == 'bed':
        try:
            generate_bed(args.bamfile, args.bed, args.contig, args.mincov, args.regionsize, args.maxdepth, args.ignoreorphans, args.stepper)
        except AttributeError as e:
            print('AttributeError: {0}\n'.format(e))
            print(parser.format_help())
    elif args.subparser_name == 'group':
        try:
            group_umis(args.outdir, args.region, args.bamfile, args.config, args.distthreshold, args.postthreshold, args.ignore, args.truncate, args.separator, args.readcount)
        except AttributeError as e:
            print('AttributeError: {0}\n'.format(e))
            print(parser.format_help())
    elif args.subparser_name == 'collapse':
        try:
            collapse(args.config, args.outdir, args.bamfile, args.reference, args.region,
                     args.umifile, args.famsize, args.countthreshold, args.percentthreshold,
                     args.postthreshold, args.maxdepth, args.truncate, args.ignoreorphans, 
                     args.stepper, args.separator, args.base_quality_score)
        except AttributeError as e:
            print('AttributeError: {0}\n'.format(e))
            print(parser.format_help())
    elif args.subparser_name == 'call':
        try:
            VCF_converter(args.config, args.outdir, args.reference, args.refthreshold, args.altthreshold, args.filterthreshold, args.famsize, args.consfiles)
        except AttributeError as e:
            print('AttributeError: {0}\n'.format(e))
            print(parser.format_help())
    elif args.subparser_name == 'run':
        try:
            run_scripts(args.outdir, args.config, args.bamfile, args.reference,
                        args.famsize, args.bedfile, args.countthreshold, args.consensusthreshold,
                        args.postthreshold, args.distthreshold, args.refthreshold,
                        args.altthreshold, args.filterthreshold, args.percentthreshold, args.maxdepth,
                        args.truncate, args.ignoreorphans, args.ignore, args.stepper,
                        args.merge, args.plot, args.report, args.call, args.mincov,
                        args.minratio, args.minumis, args.minchildren, args.extension,
                        args.sample, args.mem, args.mypython, args.mydebarcer, args.project,
                        args.separator, args.base_quality_score, args.readcount)
        except AttributeError as e:
            print('AttributeError: {0}\n'.format(e))
            print(parser.format_help())
    elif args.subparser_name == 'merge':
        try:
            merge_files(args.outdir, args.files, args.datatype)
        except AttributeError as e:
            print('AttributeError: {0}\n'.format(e))
            print(parser.format_help())
    elif args.subparser_name == 'plot':
        try:
            generate_figures(args.directory, args.config, args.extension, args.report, args.sample, args.mincov, args.minratio, args.minumis, args.minchildren, args.refthreshold)
        except AttributeError as e:
            print('AttributeError: {0}\n'.format(e))
            print(parser.format_help())
    elif args.subparser_name == 'report':
        try:
            report(args.directory, args.extension, args.mincov, args.minratio, args.minumis, args.minchildren, sample=args.sample)
        except AttributeError as e:
            print('AttributeError: {0}\n'.format(e))
            print(parser.format_help())
    elif args.subparser_name is None:
        print(parser.format_help())
    


