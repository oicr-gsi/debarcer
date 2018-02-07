
import sys
import pysam
import configparser
import argparse
import importlib.util
import operator


## Argument parsing and error handling
def handle_arg(var, alt, config, error):
    
    if var is None:
 
        if config is not None:
            var = alt
        
        if var is None:
            raise ValueError(error)
            sys.exit(1)
        
    return var


parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bam_file',    help='Path to your BAM file.')
parser.add_argument('-r', '--region',      help='Region to analyze (string of the form chrX:posA-posB).')
parser.add_argument('-o', '--output_path', help='Path to write output files to.')
parser.add_argument('-c', '--config',      help='Path to your config file.')
parser.add_argument('-t', '--tally',       help='Path to your tally (output of UMI_count.py).')

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
region_start = int(region.split(":")[1].split("-")[0])
region_end   = int(region.split(":")[1].split("-")[1])

bam_file    = handle_arg(args.bam_file, config['PATHS']['bam_file'], config, 'No BAM file provided in args or config.')
output_path = handle_arg(args.output_path, config['PATHS']['output_path'], config, 'No output path provided in args or config.')
tally_file  = handle_arg(args.tally, output_path + '/' + region + '.tally', config, 'No tally file provided.')

routine_file = config['PATHS']['routine_file']
spec         = importlib.util.spec_from_file_location("routine", routine_file)
routine      = importlib.util.module_from_spec(spec)
spec.loader.exec_module(routine)

with pysam.FastaFile(config['PATHS']['reference_file']) as reader:
    ref_seq = reader.fetch(contig, region_start, region_end).upper()


## Lists of UMI+Pos pairs with count >= cutoff
families = {}
cutoff   = int(config['SETTINGS']['min_family_size']) if config else 1

with open(tally_file, "r") as reader:
    lines = reader.readlines()
    
for line in lines[1:]:

    umi, pos, count, *rest = line.split("\t")

    if(int(count) >= cutoff):
        families[umi + pos] = 0
        

## Keys: each base position in the region
## Values: tables of A,T,C,G (etc) counts from each UMI+Pos family
consensus_seq = routine.routine(families, contig, region_start, region_end, bam_file, config_file)

## Build output
output_file = output_path + "/{}:{}-{}.cons".format(contig, region_start, region_end)
writer = open(output_file, "w")

##TEST
with open(output_file + ".testdict", "w") as testwrite:
    for family in consensus_seq[36164450]:
        testwrite.write("{}\n".format(family))
        for base in consensus_seq[36164450][family]:
            testwrite.write("{} = {}, ".format(base, consensus_seq[36164450][family][base]))

for base_pos in range(region_start, region_end):

    if base_pos in consensus_seq:

        consensuses = {} ## each family's consensus base (incl. indels, "N"s)

        for family in consensus_seq[base_pos]:
        
            cons_base    = max(consensus_seq[base_pos][family].items(), key = operator.itemgetter(1))[0]
            cons_denom   = sum(consensus_seq[base_pos][family].values())
            cons_percent = (consensus_seq[base_pos][family][cons_base]/cons_denom) * 100
    
            percent_threshold = float(config['SETTINGS']['percent_consensus_threshold']) if config else 70.0
            count_threshold   = int(config['SETTINGS']['count_consensus_threshold']) if config else 3

            if cons_percent >= percent_threshold and consensus_seq[base_pos][family][cons_base] >= count_threshold:
            
                if cons_base in consensuses:
                    consensuses[cons_base] += 1

                else:
                    consensuses[cons_base] = 1

            else:

                if "U" in consensuses: ## "U" (unclear) - placeholder for programmatically ambiguous bases, TODO revisit
                    consensuses["U"] += 1

                else:
                    consensuses["U"] = 1
        
        ref_base = ref_seq[base_pos-region_start]

        writer.write("{} (ref. {}): ".format(base_pos + 1, ref_base))
        
        for base, count in sorted(consensuses.items(), key = operator.itemgetter(1), reverse=True):
            writer.write("{}={} ".format(base, str(count)))

        writer.write("\n")

    else:

        writer.write("{}: Missing\n".format(base_pos))

writer.close()

