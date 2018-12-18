import os.path
import argparse
import configparser
import pickle
import sys
import os
import datetime
import json
import csv
import time
from src.handle_args import handle_arg
from src.handle_args import arg_exists
from src.handle_args import config_validation
from src.preprocess_fastqs import reheader_fastqs
from src.umi_error_correct import get_umi_families, umi_count, umi_datafile
from src.generate_consensus import generate_consensus_output
from src.generate_vcf import generate_vcf_output
from src.generate_vcf import create_consensus_output
from src.generate_vcf import get_vcf_output2, get_vcf_output, check_consfile
from src.generate_merge import merge_umi_datafiles2, concat_cons, modify_cons, temp_cons, submit_jobs, check_merge_flag
from src.create_plots import umi_plot, cons_plot, check_file


"""
debarcer.py - main interface for Debarcer
=========================================

Purpose
-------
Debarcer (De-Barcoding and Error Correction) is a package
for working with next-gen sequencing data that contains
molecular barcodes.

Authors: Theodore Bodak, Isha Warikoo
Copyright (c) 2018 GSI, Ontario Institute for Cancer Research
"""

def timestamp():
	"""Returns the current time in a nice format for log files."""
	return "[{}] ".format(str(datetime.datetime.now()).split('.')[0])


def preprocess_reads(args):
	"""
	Preprocesses fastq files by removing UMIs from reads and appending
	them to the read names.
	"""

	if args.config:
		config = configparser.ConfigParser()
		config.read(args.config)
		#config.read('../config/demo_config.ini')
		config_validation(conf_paths = dict(config.items('PATHS'))) ##Check whether PATHS in config file exist
	else:
		config = None

	prepfile = handle_arg(args.prepfile, config['PATHS']['prep_file'] if config else None,
					'ERR: No prepfile provided in args or config.')
	output_path = handle_arg(args.output_path, config['PATHS']['output_path'] if config else None,
					'ERR: No output path provided in args or config.')

	arg_exists(sys.argv) ##Check whether args directories/files exist	

	
	reheader_fastqs(
		r1_file=args.read1, 
		r2_file=args.read2, 
		r3_file=args.read3, 
		output_path=output_path, 
		prepname=args.prepname,
		prepfile=prepfile)
	


def group_umis(args):
	"""Groups and error-corrects UMIs into families."""

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

	print(timestamp() + "Grouping UMIs...")

	## Generate an error-corrected list of UMI families
	umi_families, umi_groups = get_umi_families(
		contig=contig,
		region_start=region_start,
		region_end=region_end,
		bam_file=bam_file,
		config=config)

	total_parent_umi_count, total_child_umi_count, num_of_children, freq_of_parent_umis = umi_datafile(umi_groups)

	filename="{}/datafile_{}.csv".format(output_path,region)
	headers = ['CHR', 'START', 'END', 'PTU', 'CTU', 'CHILD_NUMS', 'FREQ_PARENTS']
	csv.register_dialect('myDialect', delimiter='\t', quoting=csv.QUOTE_NONE)
	csvrow = {'CHR' : contig, 'START' : str(region_start), 'END' : str(region_end), 'PTU' : str(total_parent_umi_count), 'CTU' : str(total_child_umi_count), 'CHILD_NUMS': num_of_children, 'FREQ_PARENTS' : freq_of_parent_umis}
	info = [contig, region_start, region_end, total_parent_umi_count, total_child_umi_count, num_of_children, freq_of_parent_umis]

	file = open(filename, "w")
	writer = csv.DictWriter(file, dialect='myDialect', fieldnames=headers)
	writer.writeheader()
	writer.writerow(csvrow)
	

	umi_file = "{}/{}.umis".format(output_path, region)
	pickle.dump(umi_families, open(umi_file, "wb"))
	
	print(timestamp() + "UMI grouping complete. Output written to {}.".format(output_path))




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

	region = args.region
	if any(x not in region for x in ["chr", ":", "-"]):
		raise ValueError('Incorrect region string (should look like chr1:1200000-1250000).')
		sys.exit(1)

	#contig = region.split(":")[0]
	#region_start = int(region.split(":")[1].split("-")[0])
	#region_end = int(region.split(":")[1].split("-")[1])

	contig = region.split("_")[0].split(":")[0]
	region_start = region.split("_")[0]
	region_end = region.split("_")[1]

	output_path = handle_arg(args.output_path, config['PATHS']['output_path'] if config else None, 
					'No output path provided in args or config.')

	arg_exists(sys.argv) ##Check whether args directories/files exist
	print(timestamp() + "Generating VCFs...")

	#generate_vcf_output(cons_file=cons_file, f_sizes=f_sizes, contig=contig, region_start=region_start, region_end=region_end, output_path=output_path, config=config)
	#create_consensus_output(cons_file=cons_file, f_sizes=f_sizes, contig=contig, region_start=region_start, region_end=region_end, bam_file=bam_file, output_path=output_path, config=config)
  
	get_vcf_output(cons_file=cons_file, contig=contig, region_start=region_start, region_end=region_end, output_path=output_path, config=config)

	print(timestamp() + "VCFs generated. VCF files written to {}.".format(output_path))




def call_variants2(args):
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

	get_vcf_output2(cons_file=cons_file, region_start=region_start, region_end=region_end, output_path=output_path, config=config)
		
	#call the merged version of the function
	#get_vcf_output2(cons_file=cons_file, region_start=region_start, region_end=region_end, output_path=output_path, config=config, is_cons_merged=is_cons_merged)

	print(timestamp() + "VCFs generated. VCF files written to {}.".format(output_path))





def call_variants_modified(args):
	"""Generates VCF files from given cons file."""

	if args.config:
		config = configparser.ConfigParser()
		config.read(args.config)
		config_validation(conf_paths = dict(config.items('PATHS'))) ##Check whether PATHS in config file exist
	else:
		config = None

	cons_file = args.cons_file

	if args.f_sizes:
		f_sizes = args.f_sizes.split(',')
	else:
		f_sizes = [n for n in config['SETTINGS']['min_family_sizes'].split(',')] if config else [1, 2, 5]

	#If the cons file is not the merged file
	if check_consfile(cons_file):
		
		region = args.region
		if any(x not in region for x in ["chr", ":", "-"]):
			raise ValueError('Incorrect region string (should look like chr1:1200000-1250000).')
			sys.exit(1)

		contig = region.split(":")[0]
		region_start = int(region.split(":")[1].split("-")[0])
		region_end = int(region.split(":")[1].split("-")[1])

		output_path = handle_arg(args.output_path, config['PATHS']['output_path'] if config else None, 
						'No output path provided in args or config.')


		arg_exists(sys.argv) ##Check whether args directories/files exist

		print(timestamp() + "Generating VCFs...")

		#generate_vcf_output(cons_file=cons_file, f_sizes=f_sizes, contig=contig, region_start=region_start, region_end=region_end, output_path=output_path, config=config)

		#create_consensus_output(cons_file=cons_file, f_sizes=f_sizes, contig=contig, region_start=region_start, region_end=region_end, bam_file=bam_file, output_path=output_path, config=config)

		get_vcf_output(cons_file=cons_file, contig=contig, region_start=region_start, region_end=region_end, output_path=output_path, config=config)

		print(timestamp() + "VCFs generated. VCF files written to {}.".format(output_path))


	else:
		merged_file = cons_file

		headers = ['INTVL', 'CHROM', 'POS', 'REF', 'A', 'C', 'G',' T', 'I', 'D', 'N', 'RAWDP', 'CONSDP', 'FAM', 'REF_FREQ', 'MEAN_FAM']
		regions = []
		region_rows = []
		f = open(merged_file, "r")
		reader = csv.DictReader(f, delimiter='\t', fieldnames=headers)
		next(reader)
		row_counter = 1
		last_row = None
		index = 0		

		#region_rows is a list of tuples: [(start_row_no, end_row_no)...] listing start and end row numbers of the merged file, of data specific to unique regions 
		for row in reader:
			if row['INTVL'] not in regions:
				regions.append(row['INTVL'])
				if last_row is not None:
					region_rows.append(tuple((first_row,last_row)))
				first_row = row_counter
			else:
				last_row = row_counter
			row_counter+=1

		for region in regions:

			cons_file = temp_cons(output_dir, region, row[index], merged_file)

			contig = region.split(":")[0]
			region_start = int(region.split(":")[1].split("-")[0])
			region_end = int(region.split(":")[1].split("-")[1])
		
			print(consfile)

			#arg_exists(sys.argv) ##Check whether args directories/files exist

			#print(timestamp() + "Generating VCFs...")

			#get_vcf_output(cons_file=cons_file, contig=contig, region_start=region_start, region_end=region_end, output_path=output_path, config=config)

			#print(timestamp() + "VCFs generated. VCF files written to {}.".format(output_path))

			index+=1


def find_pos(lines):
	for i in range(len(lines)):
		if (i < len(lines) and 'chr' in lines[i][0]):
			first_pos = i
			return first_pos



def run_scripts(args):

	bamfile = args.bam_file
	bedfile = args.bed_file
	output_dir = args.output_path

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
	print("Making directories")
	if not os.path.exists(output_dir+"umifiles"):
		os.makedirs(output_dir+"umifiles")
	if not os.path.exists(output_dir+"consfiles"):
		os.makedirs(output_dir+"consfiles")
	if not os.path.exists(output_dir+"vcffiles"):
		os.makedirs(output_dir+"vcffiles")
	if not os.path.exists(output_dir+"logs"):
		os.makedirs(output_dir+"logs")


	debarcer_path = os.getcwd()+"/"
	#print(debarcer_path)

	"""
	#Read bedfile
	print("Reading bed files")
	with open(bedfile) as textFile:
		lines = [line.split() for line in textFile]
	index = find_pos(lines)


	#Create scripts for all subprocesses
	print("Creating scripts for all subprocesses")
	submit_jobs(bamfile, bedfile, output_dir, config_path, index, debarcer_path)

	"""

	#Check UMI job status before merging files
	print("Checking UMI job status...")
	umi_job_flag = False
	while umi_job_flag == False:
		umi_job_flag = check_umi_status(output_dir)
		print(umi_job_flag)

	print("Merging UMI datafiles...")
	merge_umi_datafiles2(output_dir)

	print("Checking CONS job status...")
	cons_job_flag = False
	while cons_job_flag == False:
		cons_job_flag = check_cons_status(output_dir)
		print(cons_job_flag)


	print("Merging cons...")
	concat_cons(output_dir, config)
	print("Finished. Output written to: "+output_dir)




	"""
	#Not needed
	if merge == 'umi':
		merge_umi_datafiles2(output_dir)
	elif merge == 'cons':
		print("Merging cons...")
		#merge_cons_files2(output_dir)
		#modify_cons("/u/iwarikoo/Debarcer2/d_output/haloplex_9538005/consfiles/chr13:28592509-28592706.cons", output_dir)
		#concat_cons(output_dir, config_path)
		concat_cons(output_dir, config)
		print("Finished")

	print("---QRLOG Jobs---")
	os.system("qstat -r | grep 'QRLOG*' > jobs.txt")
	print("---UMI Jobs---")
	os.system("qstat -r | grep 'UMI*' >> jobs.txt")
	print("---CONS Jobs---")
	os.system("qstat -r | grep 'CONS*' >> jobs.txt")

	job_flag = os.stat("./jobs.txt").st_size == 0
	print(job_flag)
	"""

def check_umi_status(output_dir):
	time.sleep(3)
	os.system("qstat -r | grep 'UMI*' > "+output_dir+"temp_umi_jobs.txt")
	job_flag = os.stat(output_dir+"temp_umi_jobs.txt").st_size == 0

	return job_flag

def check_cons_status(output_dir):
	time.sleep(3)
	os.system("qstat -r | grep 'CONS*' > "+output_dir+"temp_cons_jobs.txt")
	job_flag = os.stat(output_dir+"temp_cons_jobs.txt").st_size == 0

	return job_flag


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
	parser = argparse.ArgumentParser(description="A package for De-Barcoding and Error Correction" \
	 											 " of sequencing data containing molecular barcodes.")
	
	subparsers = parser.add_subparsers()

		
	## Preprocess command - requires unprocessed fastq file(s)
	p_parser = subparsers.add_parser('preprocess', help="Preprocess mode for processing fastq files.")
	p_parser.add_argument('-o', '--output_path', help='Path to write updated fastq files to.', required=True)
	p_parser.add_argument('-r1', '--read1', help='Path to first FASTQ file.', required=True)
	p_parser.add_argument('-r2', '--read2', help='Path to second FASTQ file, if applicable.')
	p_parser.add_argument('-r3', '--read3', help='Path to third FASTQ file, if applicable.')
	p_parser.add_argument('-p', '--prepname', help='Name of library prep to  use (defined in library_prep_types.ini).', required=True)
	p_parser.add_argument('-pf', '--prepfile', help='Path to your library_prep_types.ini file.')
	p_parser.add_argument('-c', '--config', help='Path to your config file.')
	p_parser.set_defaults(func=preprocess_reads)

	## UMI group command - requires BAM file
	g_parser = subparsers.add_parser('group', help="Groups and error-corrects UMIs into families.")
	g_parser.add_argument('-o', '--output_path', help='Path to write output files to.')
	g_parser.add_argument('-r', '--region', help='Region to find UMIs in (string of the form chrX:posA-posB).', required=True)
	g_parser.add_argument('-b', '--bam_file', help='Path to your BAM file.')
	g_parser.add_argument('-c', '--config', help='Path to your config file.')
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
	v_parser.set_defaults(func=call_variants2)

	##Run scripts command - requires bed file, and generates scripts for umi grouping, collapse and call functions
	s_parser = subparsers.add_parser('run', help="Generate scripts for umi grouping, collapse and call functions for target regions specified by the BED file.")
	s_parser.add_argument('-o', '--output_path', help='Path to write output files to.', required=True)
	s_parser.add_argument('-be', '--bed_file', help='Path to your BED fle.', required=True)
	s_parser.add_argument('-b', '--bam_file', help='Path to your BAM file.', required=True)
	s_parser.add_argument('-c', '--config', help='Path to your config file.')
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
		print(e)
