import argparse
import configparser
import sys
import os
import datetime
import json
import time
from src.preprocess_fastqs import reheader_fastqs
from src.umi_error_correct import get_umi_families, umi_datafile
from src.generate_consensus import generate_consensus_output
from src.generate_vcf import generate_vcf_output
from src.generate_vcf import create_consensus_output, get_vcf_output, check_consfile
from src.get_run_data import merge_umi_datafiles, concat_cons, modify_cons, submit_jobs, check_job_status, find_pos
from src.create_plots import umi_plot, cons_plot, check_file
from src.utilities import CheckRegionFormat, GetOutputDir, GetInputFiles, GetThresholds
    
    

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

def timestamp():
    """Returns the current time in a nice format for log files."""
    return "[{}] ".format(str(datetime.datetime.now()).split('.')[0])


def preprocess_reads(args):
    '''
    (list) -> None
    Preprocesses fastq files by removing UMIs from reads, appending them to 
    the read names and writing new fastqs 
	 
    :param outdir: Directory where new fastqs are written. From command or config
    :param read1: Path to first FASTQ file
    :param read2: Path to second FASTQ file
    :param read3: Path to third FASTQ file
    :param prepname: Name of the library preparation
    :param prepfile: Path to the library preparation ini file. From command or config
    :param config: Path to the config file. 
    :param prefix: Prefix for naming umi-reheradered fastqs 
    '''
        
    # get output directory from the config or command. set to current dir if not provided
    outdir = GetOutputDir(args.config, args.outdir)
    # create outputdir if doesn't exist
    if os.path.isdir(outdir) == False:
        if os.path.isfile(outdir) == True:
            raise ValueError('ERR: Output directory cannot be a file')
        else:
            os.makedirs(outdir)
    
    # get input prep file from config or command
    prepfile = GetInputFiles(args.config, args.prepfile, 'prep_file')
    
    # check that files are valid
    for i in [args.read1, args.read2, args.read3]:
        # check that argument (or file) is provided/exists 
        if i:
            # check if provided file is file
            if os.path.isfile(i) == False:
                # raise error and exit
                raise ValueError('ERR: prepfile is not a valid file')
                sys.exit(1)
    
    # reheader fastqs and add umi in new fastqs header
    reheader_fastqs(args.read1, outdir, args.prepname, prepfile, r2=args.read2, r3=args.read3, prefix=args.prefix)
	 
    
    

def group_umis(args):
    '''
    (list) -> None
       
    :param outdir: Directory where .umis, and datafiles are written
    :param region: A string with region coordinates chrN:posA-posB. posA and posB are 1-based included
    :param bamfile: Path to the bam file
    :param config: Path to your config file
    :param distthreshold: Hamming distance threshold for connecting parent-children umis
    :param postthreshold: Distance threshold in bp for defining families within groups
    :param ignore: Keep the most abundant family and ignore families at other positions within each group if True. Default is False
    
    Groups by hamming distance and form families based on physical distances within groups
    '''
    
    # get output directory from the config or command. set to current dir if not provided
    outdir = GetOutputDir(args.config, args.outdir)
    # create outputdir if doesn't exist
    if os.path.isdir(outdir) == False:
        if os.path.isfile(outdir) == True:
            raise ValueError('ERR: Output directory cannot be a file')
        else:
            os.makedirs(outdir)
    
    # create subdirectoy structure
    UmiDir = os.path.join(outdir, 'Umifiles')
    DataDir = os.path.join(outdir, 'Datafiles')
    for i in [UmiDir, DataDir]:
        if os.path.isdir(i) == False:
            os.mkdir(i)
    
    # get input bam from config or command
    bam_file = GetInputFiles(args.config, args.bamfile, 'bam_file')
    
    # check that region is properly formatted
    region = args.region
    CheckRegionFormat(region)
    # get chromosome
    contig = region.split(":")[0]
    # get 1-based inclusive region coordinates   
    region_start, region_end = int(region.split(":")[1].split("-")[0]), int(region.split(":")[1].split("-")[1])
    # convert coordinates to 0-based hal opened coordinates
    region_start = region_start -1
    
    # get umi position and distance thresholds 
    pos_threshold = GetThresholds(args.config, 'umi_family_pos_threshold', args.postthreshold)
    dist_threshold = GetThresholds(args.config, 'umi_edit_distance_threshold', args.distthreshold)
                
    print(timestamp() + "Grouping UMIs...")
    
    # Generate UMI families within groups using the position of the most frequent umi as reference for each family
    # keep the most abundant family within group and ignore others if args.ignore is True
    umi_families, umi_groups = get_umi_families(contig, region_start, region_end, bam_file, pos_threshold, dist_threshold, args.ignore)
        
    # get the number of parent umis, number of children and number of parent given a number of children
    filename= os.path.join(DataDir, 'datafile_{}.csv'.format(region))
    header = ['CHR', 'START', 'END', 'PTU', 'CTU', 'CHILD_NUMS', 'FREQ_PARENTS']
    # use 1-based inclusive coordinates in datafile output
    info = [contig, str(region_start + 1), str(region_end)] + umi_datafile(umi_groups)
    with open(filename, 'w') as newfile:
        newfile.write('\t'.join(header) + '\n')
        newfile.write('\t'.join(info) + '\n')
    
    # save umi families as a json. positions in the json are 0-based half opened
    umi_file = os.path.join(UmiDir, '{}.json'.format(region))
    with open(umi_file, 'w') as newfile:
        json.dump(umi_families, newfile, sort_keys = True, indent=4)
        
    print(timestamp() + "UMI grouping complete. Output written to {}.".format(outdir))



def collapse(args):
    '''
    (list) -> None
    
    :param config: Path to the config file
    :param outdir: Directory where consensus files are written
    :param bamfile: Path to the BAM file
    :param reference: Path to the reference genome
    :param region: Region coordinates (1-based included) to search for UMIs (eg. chrN:posA-posB)
    :param umifile: Path to the .umis file
    :param famsize: Comma-separated list of minimum umi family size to collapase on
    :param countthreshold: Base count threshold in pileup column
    :param percentthreshold: Base percent threshold in pileup column
    :param refthreshold: Reference threshold
    :param allthreshold: Allele threshold
    :param postthreshold: Umi position threshold for grouping umis together
    :param maxdepth: Maximum read depth. Default is 1000000
    :param truncate: Only consider pileup columns in given region. Default is True\
    :param ignoreorphans: Ignore orphans (paired reads that are not in a proper pair). Default is True
       
    Base collapses from given BAM and umi family file
    '''
    
    # get output directory from the config or command. set to current dir if not provided
    outdir = GetOutputDir(args.config, args.outdir)
    # create outputdir if doesn't exist
    if os.path.isdir(outdir) == False:
        if os.path.isfile(outdir) == True:
            raise ValueError('ERR: Output directory cannot be a file')
        else:
            os.makedirs(outdir)
    
    # create subdirectoy structure
    ConsDir = os.path.join(outdir, 'Consfiles')
    if os.path.isdir(ConsDir) == False:
        os.mkdir(ConsDir)
    
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
        raise ValueError("ERR: Unable to load umi json file")
        
    print(timestamp() + "Generating consensus...")

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
    generate_consensus_output(reference, contig, region_start, region_end, bam_file, umi_families, ConsDir, fam_size, pos_threshold, percent_threshold, count_threshold, ref_threshold, all_threshold, max_depth=args.maxdepth, truncate=args.truncate, ignore_orphans=args.ignoreorphans)
 
    print(timestamp() + "Consensus generated. Consensus file written to {}.".format(ConsDir))


def VCF_converter(args):
    '''
    
    
    
    
    
    
    Converts consensus files tino VCF format
    '''

    # get output directory from the config or command. set to current dir if not provided
    outdir = GetOutputDir(args.config, args.outdir)
    # create outputdir if doesn't exist
    if os.path.isdir(outdir) == False:
        if os.path.isfile(outdir) == True:
            raise ValueError('ERR: Output directory cannot be a file')
        else:
            os.makedirs(outdir)
    
    # create vcf dir
    VCFDir = os.path.join(outdir, 'VCFfiles')
    if os.path.isdir(VCFDir) == False:
        os.mkdir(VCFDir)
    
    # check that region is properly formatted
    region = args.region
    CheckRegionFormat(region)
    # get chromosome 
    contig = region.split(":")[0]
    # get 1-based inclusive region coordinates
    region_start, region_end = int(region.split(":")[1].split("-")[0]), int(region.split(":")[1].split("-")[1])
    # convert coordinates to 0-based half opened coordinates
#    region_start = region_start -1
        
    
    
#    cons_is_merged = check_consfile(cons_file)
#
#    if cons_is_merged:
#        region_start = region.split("_")[0]
#        region_end = region.split("_")[1]		
#    else:
#        if any(x not in region for x in ["chr", ":", "-"]):
#            raise ValueError('Incorrect region string (should look like chr1:1200000-1250000).')
#            sys.exit(1)
#
#        contig = region.split(":")[0]
#        region_start = int(region.split(":")[1].split("-")[0])
#        region_end = int(region.split(":")[1].split("-")[1])

    print(timestamp() + "Generating VCFs...")

    get_vcf_output(args.consfile, region_start, region_end, outdir, args.config)
    
    print(timestamp() + "VCFs generated. VCF files written to {}.".format(outdir))

def run_scripts(args):
    '''
    (list) -> 
    
    
    
    
    
    '''


    # get bam file from config or command
    bamfile = GetInputFiles(args.config, args.bamfile, 'bam_file')
    
    # get output directory from the config or command. set to current dir if not provided
    outdir = GetOutputDir(args.config, args.outdir)
    # create outputdir if doesn't exist
    if os.path.isdir(outdir) == False:
        if os.path.isfile(outdir) == True:
            raise ValueError('ERR: Output directory cannot be a file')
        else:
            os.makedirs(outdir)
    
    # create subdirectoy structure
    UmiDir = os.path.join(outdir, 'Umifiles')
    ConsDir = os.path.join(outdir, 'Consfiles')
    VCFDir = os.path.join(outdir, 'VCFfiles')
    DataDir = os.path.join(outdir, 'Datafiles')
    QsubDir = os.path.join(outdir, 'Qsubs')
    for i in [UmiDir, ConsDir, VCFDir, DataDir, QsubDir]:
        if os.path.isdir(i) == False:
            os.mkdir(i)
    LogDir = os.path.join(QsubDir, 'Logs')
    if os.path.isdir(LogDir) == False:
        os.mkdir(LogDir)
    
    bedfile = args.bed_file
    #dir = args.output_path
    id = str(args.run_id)
    output_dir = dir+id+"/"
    
    debarcer_path = os.getcwd()+"/"
    #Read bedfile
    with open(bedfile) as textFile:
        lines = [line.split() for line in textFile]
    index = find_pos(lines)

    #Create and run scripts for all subprocesses
    submit_jobs(bamfile, bedfile, output_dir, config_path, index, debarcer_path)
    
    
    
    
    
    
    
    
    
    
    
    
    
    #Check UMI job status before merging files
    print("Checking UMI job status...")
    umi_job_flag = False
    while umi_job_flag == False:
        umi_job_flag = check_job_status(output_dir, flag='umi', file='temp_umi_jobs.txt')
    
    print("Merging UMI datafiles...")
    merge_umi_datafiles(output_dir, id)
    
    print("Checking CONS job status...")
    cons_job_flag = False
    while cons_job_flag == False:
        cons_job_flag = check_job_status(output_dir, flag='cons', file='temp_cons_jobs.txt')
    
    print("Merging cons...")
    concat_cons(output_dir, config, id)
    print("Finished. Output written to: "+output_dir)

def generate_plots(args):
    if not umi_flag and not cons_flag:
        print("ERR: Specify umi or cons flag, and provide the appropriate data file")		

    if umi_flag:
        if umi_datafile:
            umi_plot_type = args.umi_flag
            umi_file = args.umi_datafile
            if check_file(umi_file, 'csv'):
                if umi_plot_type == 'rs':
                    umi_plot(output_path, file_name, umi_plot_type)
                elif umi_plot == 'all or umi_plot == None':
                    #Default: Generate all umi plots
                    umi_plot_type = 'all'
                    umi_plot(output_path, file_name, umi_plot_type)
                else:
                    print("ERR: Incorrect argument passed to the umi flag option")

            print("ERR: Incorrect or non-existing file specified. Expecting file with extension '.vcf'")
        else:
            print("ERR: Missing CSV umi data file")

    if cons_flag:
        if cons_file:
            cons_plot_type = args.cons_flag
            cons_file = args.cons_file
            if check_file(cons_file, 'cons'):
                if cons_plot_type == 'all' or cons_plot_type == None:
                    cons_plot(output_path, file_name, cons_plot_type)
                else:
                    print("ERR: Incorrect argument passed to the cons flag option")
            else:
                print("ERR: Incorrect or non-existing file specified. Expecting file with extension '.cons'")
        else:
            print("ERR: Missing CONS data file")

if __name__ == '__main__':
        
    ## Argument + config parsing and error handling
    parser = argparse.ArgumentParser(prog='debarcer.py', description="A package for De-Barcoding\
                                     and Error Correction of sequencing data containing molecular barcodes")
    subparsers = parser.add_subparsers()
       		
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
    p_parser.set_defaults(func=preprocess_reads)
    
    ## UMI group command
    g_parser = subparsers.add_parser('group', help="Groups UMIs into families.")
    g_parser.add_argument('-o', '--Outdir', dest='outdir', help='Directory where .umis and datafiles are written')
    g_parser.add_argument('-r', '--Region', dest='region', help='Region coordinates to search for UMIs. chrN:posA-posB. posA and posB are 1-based included', required=True)
    g_parser.add_argument('-b', '--Bamfile', dest='bamfile', help='Path to the BAM file')
    g_parser.add_argument('-c', '--Config', dest='config', help='Path to the config file')
    g_parser.add_argument('-d', '--Distance', dest='distthreshold', help='Hamming distance threshold for connecting parent-children umis')
    g_parser.add_argument('-p', '--Position', dest='postthreshold', help='Umi position threshold for grouping umis together')
    g_parser.add_argument('-i', '--Ignore', dest='ignore', action='store_true', help='Keep the most abundant family and ignore families at other positions within each group. Default is False')
    g_parser.set_defaults(func=group_umis)
    
    ## Base collapse command
    c_parser = subparsers.add_parser('collapse', help="Base collapsing from given UMI families file.")
    c_parser.add_argument('-c', '--Config', dest='config', help='Path to the config file')
    c_parser.add_argument('-o', '--Outdir', dest='outdir', help='Directory where consensus files are written')
    c_parser.add_argument('-b', '--Bamfile', dest='bamfile', help='Path to the BAM file')
    c_parser.add_argument('-rf', '--Reference', dest='reference', help='Path to the refeence genome')
    c_parser.add_argument('-r', '--Region', dest='region', help='Region coordinates to search for UMIs. chrN:posA-posB. posA and posB are 1-based included', required=True)
    c_parser.add_argument('-u', '--Umi', dest='umifile', help='Path to the .umis file', required=True)
    c_parser.add_argument('-f', '--Famsize', dest='famsize', help='Comma-separated list of minimum umi family size to collapase on')
    c_parser.add_argument('-ct', '--CountThreshold', dest='countthreshold', help='Base count threshold in pileup column')
    c_parser.add_argument('-pt', '--PercentThreshold', dest='percentthreshold', help='Base percent threshold in pileup column')
    c_parser.add_argument('-rt', '--RefThreshold', dest='refthreshold', help='Reference threshold')
    c_parser.add_argument('-at', '--AlleleThreshold', dest='allthreshold', help='Allele threshold')
    c_parser.add_argument('-p', '--Position', dest='postthreshold', help='Umi position threshold for grouping umis together')
    c_parser.add_argument('-m', '--MaxDepth', dest='maxdepth', default=1000000, help='Maximum read depth. Default is 1000000')
    c_parser.add_argument('-t', '--Truncate', dest='truncate', action='store_false',
                          help='If truncate is True and a region is given,\
                          only pileup columns in the exact region specificied are returned. Default is True')
    c_parser.add_argument('-i', '--IgnoreOrphans', dest='ignoreorphans', action='store_false',
                          help='Ignore orphans (paired reads that are not in a proper pair). Default is True')
    c_parser.set_defaults(func=collapse)

    ## Variant call command - requires cons file (can only run after collapse)
    v_parser = subparsers.add_parser('call', help="Convert consensus file into VCF format.")
    v_parser.add_argument('-o', '--output_path', help='Path to writer output files to.', required=True)
    v_parser.add_argument('-r', '--region', help='Region to analyze (string of the form chrX:posA-posB).', required=True)
    v_parser.add_argument('-cf', '--cons_file', help='Path to your cons file.', required=True)
    v_parser.add_argument('-f', '--f_sizes', help='Comma-separated list of family sizes to make VCF files for.', required=True)
    v_parser.add_argument('-c', '--config', help='Path to your config file.')
    v_parser.set_defaults(func=VCF_converter)
    
    ##Run scripts command - requires bed file, and generates scripts for umi grouping, collapse and call functions
    s_parser = subparsers.add_parser('run', help="Generate scripts for umi grouping, collapse and call functions for target regions specified by the BED file.")
    s_parser.add_argument('-o', '--output_path', help='Path to write output files to.', required=True)
    s_parser.add_argument('-be', '--bed_file', help='Path to your BED fle.', required=True)
    s_parser.add_argument('-b', '--bam_file', help='Path to your BAM file.', required=True)
    s_parser.add_argument('-c', '--config', help='Path to your config file.')
    s_parser.add_argument('-id', '--run_id', help='Run id.', required=True)
    s_parser.set_defaults(func=run_scripts)
    
    ##Generate graphs	
    g_parser = subparsers.add_parser('plot', help="Generate graphs for umi and cons data files.")
    g_parser.add_argument('-o', '--output_path', help='Path to write output files to.', required=True)
    g_parser.add_argument('-c', '--cons_file', help='Path to your CONS fle.')
    g_parser.add_argument('-u', '--umi_datafile', help='Path to your umi CSV data file.')
    g_parser.add_argument('-uf', '--umi_flag', help="Pass 'all' or 'rs' as a parameter, to generate all or region-specific umi plots, respectively.")
    g_parser.add_argument('-cf', '--cons_flag', help="Pass 'all' as a parameter, to generate all cons plots.")
    
    args = parser.parse_args()
    try:
        args.func(args)
    except AttributeError as e:
        print(parser.format_help())
    
