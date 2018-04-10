
import argparse
import configparser
import pickle
import sys
from src.handle_args import handle_arg
from src.umi_error_correct import get_umi_families
from src.generate_consensus import generate_consensus_output
from src.preprocess_fastqs import reheader_fastqs

"""
debarcer.py - main wrapper for Debarcer
=======================================

Purpose
-------
Debarcer (De-Barcoding and Error Correction) is a package
for working with next-gen sequencing data that contains
molecular barcodes.

Author: Theodore Bodak
Copyright (c) 2018 GSI, Ontario Institute for Cancer Research
"""

def preprocess_reads(args):
	"""
	Preprocesses fastq files by removing UMIs from reads and appending
	them to the read names.
	"""

	if args.config:
		config = configparser.ConfigParser()
		config.read(args.config)
	else:
		config = None

	prepfile = handle_arg(args.prepfile, config['PATHS']['prep_file'] if config else None,
					'ERR: No prepfile provided in args or config.')
	output_path = handle_arg(args.output_path, config['PATHS']['output_path'] if config else None,
					'ERR: No output path provided in args or config.')

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

	## Generate an error-corrected list of UMI families
	umi_families = get_umi_families(
		contig=contig,
		region_start=region_start,
		region_end=region_end,
		bam_file=bam_file,
		config=config)

	pickle.dump(umi_families, open("{}/{}.umis".format(output_path, region), "wb"))


def call_variants(args):
	"""Calls variants from given BAM and umi family files."""

	if args.config:
		config = configparser.ConfigParser()
		config.read(args.config)
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

	if config:
		umi_file = config['PATHS']['umi_file'] if 'umi_file' in config['PATHS'] else None

	if umi_file in args:
		umi_file = args.umi_file 

	print(umi_file)

	if umi_file:
		try:
			umi_table = pickle.load(open(umi_file, "rb"))
		except IOError:
			print("ERR: Unable to load .umis file.", file=sys.stderr)
			sys.exit(1)
	else:
		umi_table = None

	generate_consensus_output(
		contig=contig,
		region_start=region_start,
		region_end=region_end,
		bam_file=bam_file,
		umi_table=umi_table,
		output_path=output_path,
		config=config)


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

	## Variant call command - requires BAM file, UMI family file
	c_parser = subparsers.add_parser('call', help="Variant calling from given UMI families.")
	c_parser.add_argument('-o', '--output_path', help='Path to write output files to.')
	c_parser.add_argument('-r', '--region', help='Region to analyze (string of the form chrX:posA-posB).', required=True)
	c_parser.add_argument('-b', '--bam_file', help='Path to your BAM file.')
	c_parser.add_argument('-u', '--umi_file', help='Path to your .umis file.')
	c_parser.add_argument('-c', '--config', help='Path to your config file.')
	c_parser.set_defaults(func=call_variants)

	args = parser.parse_args()
	args.func(args)
