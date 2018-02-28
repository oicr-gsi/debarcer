
import sys
import argparse
import configparser

from report import coverage
from report import variants


def handle_arg(var, alt, config, error):
    """Argument parsing and error handling"""
    
    if var is None:
 
        if config is not None:
            var = alt
        
        if var is None:
            raise ValueError(error)
            sys.exit(1)
        
    return var


parser = argparse.ArgumentParser()
parser.add_argument('-t',  '--tally',       help='Path to your tally file (output of UMI_count.py).')
parser.add_argument('-o',  '--output_path', help='Path to write output files to.')
parser.add_argument('-b',  '--bam_file',    help='Path to your BAM file.')
parser.add_argument('-r',  '--region',      help='Region to tally (string of the form chrX:posA-posB).')
parser.add_argument('-c',  '--config',      help='Path to your config file.')

args = parser.parse_args()
config_file = args.config

if config_file:
    config = configparser.ConfigParser()
    config.read(config_file)
else:
    config = None
    
region = args.region
if any(x not in region for x in ["chr", ":", "-"]):
    raise ValueError('Incorrect region string (should look like chr1:1200000-1250000).')
    sys.exit(1)

contig       = region.split(":")[0]
region_start = region.split(":")[1].split("-")[0]
region_end   = region.split(":")[1].split("-")[1]

output_path = handle_arg(args.output_path, config['PATHS']['output_path'], config, 'No output path provided in args or config.')
tally_file  = handle_arg(args.tally, output_path + '/' + region + '.tally', config, 'No tally file provided.')
bam_file    = handle_arg(args.bam_file, config['PATHS']['bam_file'], config, 'No BAM file provided in args or config.')


## -- Output --
contents = []

## Variants
if config['REPORT']['run_variants'] == 'TRUE':
    contents.append(variants.get_variants(tally_file, config_file))
    contents.append(["\n"])

## Coverage
if config['REPORT']['run_coverage'] == 'TRUE':
    contents.append(coverage.get_coverage(contig, region_start, region_end, bam_file, config_file))
    contents.append(["\n"])

## --   Report plots   --
## TODO TODO TODO TODO TODO


## Write output
with open("{}/{}:{}-{}.report".format(output_path, contig, region_start, region_end), "w") as writer:

    for content in contents:
        for line in content:
            writer.write(line)