import argparse
import configparser
import pickle
import sys
import os
import datetime
import json
import csv
import time
from src.preprocess_fastqs import reheader_fastqs
from src.umi_error_correct import get_umi_families, umi_count, umi_datafile
from src.generate_consensus import generate_consensus_output
from src.generate_vcf import generate_vcf_output
from src.generate_vcf import create_consensus_output, get_vcf_output, check_consfile
from src.get_run_data import merge_umi_datafiles, concat_cons, modify_cons, submit_jobs, check_job_status, find_pos
from src.create_plots import umi_plot, cons_plot, check_file


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
        
    # get prepfile and outdir from config in priority
    try:
        config = configparser.ConfigParser()
        config.read(args.config)
        prepfile = config['PATHS']['prep_file']         
        outdir = config['PATHS']['output_dir']
    except:
        # check if prepfile and outdir are provided in the command
        try:
            prepfile, outdir = args.prepfile, args.outdir
        except:
            # raise error and exit if no prepfile and outdir are provided
            raise ValueError('ERR: Missing prepfile and/or output directory')
            sys.exit(1)
    finally:
        # check that prepfile is a valid file
        if os.path.isfile(prepfile) == False:
            raise ValueError('ERR: Invalid path to prepfile')
        
    # code below is executed only if prepfile and outdir are provided  
    # create outputdir if doesn't exist
    if os.path.isdir(outdir) == False:
        if os.path.isfile(outdir) == True:
            raise ValueError('ERR: Output directory cannot be a file')
        else:
            os.makedirs(outdir)
    
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
    reheader_fastqs(args.read1, outdir, args.prepname, prepfile, r2=args.read2, r3=args.read3, prefix=prefix)
	 
    
    

def group_umis(args):
    '''
    (list) -> None
       
    :param outdir: Directory where .umis, and datafiles are written
    :param region: A string with region coordinates chrN:posA-posB. posA and posB are 1-based included
    :param bamfile: Path to the bam file
    :param config: Path to your config file
    :param distance: Hamming distance threshold for connecting parent-children umis
    :param position: Distance threshold in bp for defining families within groups
    :param ignore: Keep the most abundant family and ignore families at other positions within each group if True. Default is False
    
    Groups by hamming distance and form families based on physical distances within groups
    '''
    
    # get bam and outdir from config in priority
    try:
        config = configparser.ConfigParser()
        config.read(args.config)
        bam_file = config['PATHS']['bam_file']
        outdir = config['PATHS']['outdir']
    except:
        # check if bam file and outdir are provided in the command
        try:
            bam_file, outdir = args.bamfile, args.outdir
        except:
            # raise error and exit
            raise ValueError('ERR: Missing input bam and/or output directory')
    finally:
        # check that bam is a valid file
        if os.path.isfile(bam_file) == False:
            raise ValueError('ERR: Invalid path to input bam file')
        
    # create outputdir if doesn't exist
    if os.path.isdir(outdir) == False:
        if os.path.isfile(outdir) == True:
            raise ValueError('ERR: Output directory cannot be a file')
        else:
            os.makedirs(outdir)
    
    # check that region is properly formatted
    region = args.region
    if any(i not in region for i in ["chr", ":", "-"]):
        raise ValueError('ERR: Incorrect region string (should look like chr1:1200000-1250000)')
    # get chromosome and check format 
    contig = region.split(":")[0]
    chromos = [str(i) for i in range(23)] + ['X', 'Y']
    if contig[:len('chr')] != 'chr' and contig[len('chr'):] not in chromos:
        raise ValueError('ERR: Incorrect chromosome name (should look like chr1:1200000-1250000)')
    # get region coordinates. use 1-based inclusive. this will be converted to 0-based by pysam   
    region_start, region_end = region.split(":")[1].split("-")[0], region.split(":")[1].split("-")[1]
    if region_start.isnumeric() == False or region_end.isnumeric() == False:
        raise ValueError('ERR: Incorrect start and end coordinates (should look like chr1:1200000-1250000)')
    region_start, region_end = int(region_start), int(region_end)
        
    
    # get umi position and distance thresholds from config
    try:
        config = configparser.ConfigParser()
        config.read(args.config)
        pos_threshold = int(config['SETTINGS']['umi_family_pos_threshold'])
        dist_threshold = int(config['SETTINGS']['umi_edit_distance_threshold'])
    except:
        # check pos and dist threshold provided in the command
        try:
            pos_threshold, dist_threshold = int(args.postthreshold), int(args.distthreshold)
        except:
            # raise error and exit
            raise ValueError('ERR: Missing umi position and/or distance thresholds')
    finally:
        # check that threshold are integers
        if type(pos_threshold) != int and type(dist_threshold) != int:
            raise ValueError('ERR: Umi position and distance thresholds should be integers')
        
    print(timestamp() + "Grouping UMIs...")
    
    # Generate UMI families within groups using the position of the most frequent umi as reference for each family
    # keep the most abundant family within group and ignore others if args.ignore is True
    umi_families, umi_groups = get_umi_families(contig, region_start, region_end, bam_file, pos_threshold, dist_threshold, args.ignore)
        
    # get the number of parent umis, number of children and number of parent given a number of children
    filename="{}/datafile_{}.csv".format(outdir,region)
    header = ['CHR', 'START', 'END', 'PTU', 'CTU', 'CHILD_NUMS', 'FREQ_PARENTS']
    info = [contig, str(region_start), str(region_end)] + umi_datafile(umi_groups)
    with open(filename, 'w') as newfile:
        newfile.write('\t'.join(header) + '\n')
        newfile.write('\t'.join(info) + '\n')
    
    umi_file = "{}/{}.umis".format(outdir, region)
    with open(umi_file, 'w') as newfile:
        json.dump(umi_families, newfile, sort_keys = True, indent=4)
        
    print(timestamp() + "UMI grouping complete. Output written to {}.".format(outdir))



def collapse(args):
	"""Base collapses from given BAM and umi family file."""

	if args.config:
		config = configparser.ConfigParser()
		config.read(args.config)
		config_validation(conf_paths = dict(config.items('PATHS'))) ##Check whether PATHS in config file exist
	else:
		config = None

	region = args.region
	

	if any(item not in region for item in ["chr", ":", "-"]):
		raise ValueError('ERR: Incorrect region string (should look like chr1:1200000-1250000).')
		sys.exit(1)

	contig = region.split(":")[0]
	region_start = int(region.split(":")[1].split("-")[0])
	region_end = int(region.split(":")[1].split("-")[1])

	bam_file = handle_arg(args.bam_file, config['PATHS']['bam_file'] if config else None, 
					'ERR: No BAM file provided in args or config.')
	output_path = handle_arg(args.output_path, config['PATHS']['output_path'] if config else None, 
					'ERR: No output path provided in args or config.')

	arg_exists(sys.argv) ##Check whether args directories/files exist

	if args.umi_file:
		umi_file = args.umi_file
	elif config:
		umi_file = config['PATHS']['umi_file'] if 'umi_file' in config['PATHS'] else None

	if umi_file:
		try:
			umi_table = pickle.load(open(umi_file, "rb"))
		except IOError:
			print("ERR: Unable to load .umis file.", file=sys.stderr)
			sys.exit(1)
	else:
		umi_table = None

	print(timestamp() + "Generating consensus...")

	generate_consensus_output(
		contig=contig,
		region_start=region_start,
		region_end=region_end,
		bam_file=bam_file,
		umi_table=umi_table,
		output_path=output_path,
		config=config)

	print(timestamp() + "Consensus generated. Consensus file written to {}.".format(output_path))




def call_variants(args):
	"""Generates VCF files from given cons file."""

	if args.config:
		config = configparser.ConfigParser()
		config.read(args.config)
		config_validation(conf_paths = dict(config.items('PATHS'))) ##Check whether PATHS in config file exist
	else:
		config = None
	

	cons_file = args.cons_file
	f_sizes = args.f_sizes.split(',')

	output_path = handle_arg(args.output_path, config['PATHS']['output_path'] if config else None,
					'No output path provided in args or config.')
	
	region = args.region
	arg_exists(sys.argv) ##Check whether args directories/files exist

	cons_is_merged = check_consfile(cons_file)

	if cons_is_merged:
		region_start = region.split("_")[0]
		region_end = region.split("_")[1]		

	else:
		if any(x not in region for x in ["chr", ":", "-"]):
			raise ValueError('Incorrect region string (should look like chr1:1200000-1250000).')
			sys.exit(1)

		contig = region.split(":")[0]
		region_start = int(region.split(":")[1].split("-")[0])
		region_end = int(region.split(":")[1].split("-")[1])

	print(timestamp() + "Generating VCFs...")

	get_vcf_output(cons_file=cons_file, region_start=region_start, region_end=region_end, output_path=output_path, config=config)
		
	print(timestamp() + "VCFs generated. VCF files written to {}.".format(output_path))





def run_scripts(args):

	bamfile = args.bam_file
	bedfile = args.bed_file
	dir = args.output_path
	id = str(args.run_id)
	output_dir = dir+id+"/"

	if args.config:
		config_path = args.config
		config = configparser.ConfigParser()
		config.read(args.config)
		config_validation(conf_paths = dict(config.items('PATHS'))) ##Check whether PATHS in config file exist
	else:
		config_path = None
		config = None


	arg_exists(sys.argv) ##Check whether args directories/files exist
	
	#Make directories 
	if not os.path.exists(output_dir+"umifiles"):
		os.makedirs(output_dir+"umifiles")
	if not os.path.exists(output_dir+"consfiles"):
		os.makedirs(output_dir+"consfiles")
	if not os.path.exists(output_dir+"vcffiles"):
		os.makedirs(output_dir+"vcffiles")


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
       		
    ## Preprocess command - requires unprocessed fastq file(s)
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
    
    ## UMI group command - requires BAM file
    g_parser = subparsers.add_parser('group', help="Groups and error-corrects UMIs into families.")
    g_parser.add_argument('-o', '--Outdir', dest='outdir', help='Directory where .umis and datafiles are written')
    g_parser.add_argument('-r', '--Region', dest='region', help='Region coordinates to search for UMIs. chrN:posA-posB. posA and posB are 1-based included', required=True)
    g_parser.add_argument('-b', '--Bamfile', dest='bamfile', help='Path to the BAM file')
    g_parser.add_argument('-c', '--Config', dest='config', help='Path to the config file')
    g_parser.add_argument('-d', '--Distance', dest='distthreshold', help='Hamming distance threshold for connecting parent-children umis')
    g_parser.add_argument('-p', '--Position', dest='postthreshold', help='Umi position threshold for grouping umis together')
    g_parser.add_argument('-i', '--Ignore', dest='ignore', action='store_true', help='Keep the most abundant family and ignore families at other positions within each group. Default is False')
    g_parser.set_defaults(func=group_umis)
    
    ## Base collapse command - requires BAM file, UMI family file optional
    c_parser = subparsers.add_parser('collapse', help="Base collapsing from given UMI families file.")
    c_parser.add_argument('-o', '--output_path', help='Path to write output files to.')
    c_parser.add_argument('-r', '--region', help='Region to analyze (string of the form chrX:posA-posB).', required=True)
    c_parser.add_argument('-b', '--bam_file', help='Path to your BAM file.')
    c_parser.add_argument('-u', '--umi_file', help='Path to your .umis file.')
    c_parser.add_argument('-c', '--config', help='Path to your config file.')
    c_parser.add_argument('-f', '--f_sizes', help='Comma-separated list of family sizes to collapse on.') ##implement
    c_parser.set_defaults(func=collapse)
    
    ## Variant call command - requires cons file (can only run after collapse)
    v_parser = subparsers.add_parser('call', help="Variant calling from analyzed BAM file.")
    v_parser.add_argument('-o', '--output_path', help='Path to writer output files to.', required=True)
    v_parser.add_argument('-r', '--region', help='Region to analyze (string of the form chrX:posA-posB).', required=True)
    v_parser.add_argument('-cf', '--cons_file', help='Path to your cons file.', required=True)
    v_parser.add_argument('-f', '--f_sizes', help='Comma-separated list of family sizes to make VCF files for.', required=True)
    v_parser.add_argument('-c', '--config', help='Path to your config file.')
    v_parser.set_defaults(func=call_variants)
    
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
    
