
import argparse
import configparser
import pickle
import sys
import os
import datetime
import json
from src.handle_args import handle_arg
from src.handle_args import arg_exists
from src.handle_args import config_validation
from src.preprocess_fastqs import reheader_fastqs
from src.umi_error_correct import get_umi_families, umi_count
from src.generate_consensus import generate_consensus_output
from src.generate_vcf import generate_vcf_output
from src.generate_vcf import create_consensus_output
from src.generate_vcf import get_vcf_output


"""
debarcer.py - main interface for Debarcer
=========================================

Purpose
-------
Debarcer (De-Barcoding and Error Correction) is a package
for working with next-gen sequencing data that contains
molecular barcodes.

Author: Theodore Bodak
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
	umi_families = get_umi_families(
		contig=contig,
		region_start=region_start,
		region_end=region_end,
		bam_file=bam_file,
		config=config)

	umi_file = "{}/{}.umis".format(output_path, region)
	pickle.dump(umi_families, open(umi_file, "wb"))
	
	#os.system("python3.6 -m pickletools "+umi_file+" -o /u/iwarikoo/Debarcer2/d_output/haloplex_9538005/umifiles/pickletest_"+region+".txt") #Create text file
	umi_counts = umi_count(contig, region_start, region_end, bam_file)
	
	with open("/u/iwarikoo/Debarcer2/d_output/haloplex_9538005/umifiles/datafile_"+region+".txt","w") as file:
		file.write(json.dumps(umi_counts))
	#txt_writer = open("/u/iwarikoo/Debarcer2/d_output/haloplex_9538005/umifiles/datafile_"+region+".txt","w") #Pull specific data from file
	#txt_writer.write(umi_counts) 

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

	contig = region.split(":")[0]
	region_start = int(region.split(":")[1].split("-")[0])
	region_end = int(region.split(":")[1].split("-")[1])

	output_path = handle_arg(args.output_path, config['PATHS']['output_path'] if config else None, 
					'No output path provided in args or config.')
	bam_file = handle_arg(args.bam_file, config['PATHS']['bam_file'] if config else None,
                                        'ERR: No BAM file provided in args or config.')

	arg_exists(sys.argv) ##Check whether args directories/files exist

	print(timestamp() + "Generating VCFs...")

	#generate_vcf_output(cons_file=cons_file, f_sizes=f_sizes, contig=contig, region_start=region_start, region_end=region_end, output_path=output_path, config=config)
	
	#create_consensus_output(cons_file=cons_file, f_sizes=f_sizes, contig=contig, region_start=region_start, region_end=region_end, bam_file=bam_file, output_path=output_path, config=config)
  
	get_vcf_output(cons_file=cons_file, contig=contig, region_start=region_start, region_end=region_end, output_path=output_path, config=config)

	print(timestamp() + "VCFs generated. VCF files written to {}.".format(output_path))



def find_pos(lines):
	for i in range(len(lines)):
		if (i < len(lines) and 'chr' in lines[i][0]):
			first_pos = i
			return first_pos



def generate_scripts(args):
	
	bamfile = args.bed_file
	bedfile = args.bed_file
	output_dir = args.output_path
	deb_path = args.deb_path

	#Make directories 
	if not os.path.exists(output_dir+"umifiles"):
		os.makedirs(output_dir+"umifiles")
	if not os.path.exists(output_dir+"consfiles"):
		os.makedirs(output_dir+"consfiles")
	if not os.path.exists(output_dir+"vcffiles"):
		os.makedirs(output_dir+"vcffiles")


	#Read bedfile
	with open(bedfile) as textFile:
		lines = [line.split() for line in textFile]
	index = find_pos(lines)


	#Create scripts for all subprocesses
	for i in range(index,len(lines)):
		chromosome = lines[i][0]
		pos1 = lines[i][1]
		pos2 = lines[i][2]

		print(pos1)


		#Create umi scripts
		f = open(output_dir+"umifiles/umigrp_"+chromosome+"_"+pos1+".sh","w")
		os.system("chmod +x "+output_dir+"umifiles/umigrp_"+chromosome+"_"+pos1+".sh")
		f.write("module load /.mounts/labs/PDE/Modules/modulefiles/python-gsi/3.6.4")
		f.write("\npython3.6 "+deb_path+"debarcer.py group -o "+output_dir+"umifiles/ -r "+chromosome+"\:"+pos1+"-"+pos2+" -b "+bamfile+" -c ./config/demo_config.ini")

		#Create cons scripts
		f = open(output_dir+"consfiles/cons_"+chromosome+"_"+pos1+".sh","w")
		os.system("chmod +x "+output_dir+"consfiles/cons_"+chromosome+"_"+pos1+".sh")
		f.write("module load /.mounts/labs/PDE/Modules/modulefiles/python-gsi/3.6.4")
		f.write("\npython3.6 "+deb_path+"debarcer.py collapse -o "+output_dir+"consfiles/ -r "+chromosome+"\:"+pos1+"-"+pos2+" -b "+bamfile+" -u "+output_dir+"umifiles/"+chromosome+"\:"+pos1+"-"+pos2+".umis -c ./config/demo_config.ini")

		#Create vcf scripts
		f = open(output_dir+"vcffiles/call_"+chromosome+"_"+pos1+".sh","w")
		os.system("chmod +x "+output_dir+"vcffiles/call_"+chromosome+"_"+pos1+".sh")
		f.write("module load /.mounts/labs/PDE/Modules/modulefiles/python-gsi/3.6.4")
		f.write("\npython3.6 "+deb_path+"debarcer.py call -o "+output_dir+"vcffiles/ -r "+chromosome+"\:"+pos1+"-"+pos2+" -cf "+output_dir+"consfiles/"+chromosome+"\:"+pos1+"-"+pos2+".cons -f 1,2,5 -c ./config/demo_config.ini")




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
	v_parser.add_argument('-b', '--bam_file', help='Path to your BAM file.')
	v_parser.set_defaults(func=call_variants)

	##Generate scripts command - requires bed file, and generates scripts for umi grouping, collapse and call functions
	s_parser = subparsers.add_parser('generate', help="Generate scripts for umi grouping, collapse and call functions for target regions specified by the BED file.")
	s_parser.add_argument('-o', '--output_path', help='Path to write output files to.', required=True)
	s_parser.add_argument('-be', '--bed_file', help='Path to your BED fle.', required=True)
	s_parser.add_argument('-b', '--bam_file', help='Path to your BAM file.', required=True)
	s_parser.add_argument('-d', '--deb_path', help='Path to the debarcer.py script.', required=True)
	s_parser.set_defaults(func=generate_scripts)

	args = parser.parse_args()
	
	try:
		args.func(args)
	except AttributeError as e:
		print(parser.format_help())
		print(e)
