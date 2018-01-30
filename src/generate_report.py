
import sys
import configparser

from report import variants
from report import coverage

region_file  = sys.argv[1]
output_path  = sys.argv[2]
bam_file     = sys.argv[3]
contig       = sys.argv[4]
region_start = sys.argv[5]
region_end   = sys.argv[6]
config_file  = sys.argv[7]

config = configparser.ConfigParser()
config.read(config_file)

contents = []

## Variants
if config['REPORT']['run_variants'] == 'TRUE':
    contents.append(variants.get_variants(region_file, config_file))
    contents.append(["\n"])

## Coverage
if config['REPORT']['run_coverage'] == 'TRUE':
    contents.append(coverage.get_coverage(contig, region_start, region_end, bam_file, config_file))
    contents.append(["\n"])

## --   Report plots   --
## TODO TODO TODO TODO TODO


## Write output
with open(output_path + "/report.txt", "w") as writer:

    for content in contents:
        for line in content:
            writer.write(line)