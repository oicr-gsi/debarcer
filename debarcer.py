
"""
debarcer.py - main wrapper for Debarcer
=======================================

Purpose
-------
Debarcer (De-Barcoding and Error Correction) is a package
for working with next-gen sequencing data that contains
molecular barcodes.

Usage
-----
--tally        Run a UMI tally (UMI_count.py). Necessary for later steps.
--bam_file     Path to your BAM file.
--bed_file     Path to your BED file.
--region       Region to analyze (string of the form chrX:posA-posB).
--output_path  Path to write output files to.
--config       Path to your config file.

Author: Theodore Bodak
Copyright (c) 2018 GSI, Ontario Institute for Cancer Research
"""

import argparse
import configparser
from src.handle_args import handle_arg
from src.UMI_count import generate_tally_output
from src.generate_consensus import generate_consensus_output

if __name__ == '__main__':

	## Argument + config parsing and error handling
	parser = argparse.ArgumentParser()
	parser.add_argument('-t',  '--tally',       help='Run a UMI tally (UMI_count.py).', action='store_true')
	parser.add_argument('-b',  '--bam_file',    help='Path to your BAM file.')
	parser.add_argument('-be', '--bed_file',    help='Path to your BED file.')
	parser.add_argument('-r',  '--region',      help='Region to analyze (string of the form chrX:posA-posB).')
	parser.add_argument('-o',  '--output_path', help='Path to write output files to.')
	parser.add_argument('-c',  '--config',      help='Path to your config file.')

	args = parser.parse_args()

	if args.config:
   		config = configparser.ConfigParser()
   		config.read(args.config)
	else:
  		config = None

	region = args.region
	if any(item not in region for item in ["chr", ":", "-"]):
	    raise ValueError('Incorrect region string (should look like chr1:1200000-1250000).')
	    sys.exit(1)

	contig       = region.split(":")[0]
	region_start = int(region.split(":")[1].split("-")[0])
	region_end   = int(region.split(":")[1].split("-")[1])

	bam_file    = handle_arg(args.bam_file, config['PATHS']['bam_file'] if config else None, 
					'No BAM file provided in args or config.')
	bed_file    = handle_arg(args.bed_file, config['PATHS']['bed_file'] if config else None, 
					'No BED file providid in args or config.')
	output_path = handle_arg(args.output_path, config['PATHS']['output_path'] if config else None, 
					'No output path provided in args or config.')

	## Generate tally file if one does not exist
	if args.tally:
		generate_tally_output(
			contig=contig,
			region_start=region_start,
			region_end=region_end,
			bed_file=bed_file,
			bam_file=bam_file,
			output_path=output_path)

	tally_file = "{}/{}:{}-{}.tally".format(output_path, contig, region_start, region_end)

	## Generate consensus output files (.cons and .vcf)
	generate_consensus_output(
		contig=contig,
		region_start=region_start,
		region_end=region_end,
		bam_file=bam_file,
		tally_file=tally_file,
		output_path=output_path,
		config=config)

	## TODO plots and additional output/stats
	## ...





	




