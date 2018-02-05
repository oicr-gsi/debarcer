
import sys
import configparser

from report import coverage
from report import variants

region_file  = sys.argv[1]
output_path  = sys.argv[2]
bam_file     = sys.argv[3]
region       = sys.argv[4]
config_file  = sys.argv[5]

config = configparser.ConfigParser()
config.read(config_file)

contents = []

contig       = region.split(":")[0]
region_start = region.split(":")[1].split("-")[0]
region_end   = region.split(":")[1].split("-")[1]

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
with open("{}/{}:{}-{}.report".format(output_path, contig, region_start, region_end), "w") as writer:

    for content in contents:
        for line in content:
            writer.write(line)