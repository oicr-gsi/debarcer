
import sys
import pysam
import configparser
import argparse
import operator
from get_consensus_seq import get_consensus_seq

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


with pysam.FastaFile(config['PATHS']['reference_file']) as reader:
    ref_seq = reader.fetch(contig, region_start, region_end).upper()


## Lists of UMI+Pos pairs with count >= f_size
families = {}
f_sizes  = [int(n) for n in config['SETTINGS']['min_family_sizes'].split(',')] if config else [1, 2, 5, 10]

with open(tally_file, "r") as reader:
    lines = reader.readlines()
    
for line in lines[1:]:

    umi, pos, count, *rest = line.split("\t")
    
    for f_size in f_sizes:
        
        if(int(count) >= f_size):
            
            if f_size not in families:
                families[f_size] = {}
                
            families[f_size][umi + pos] = 0
        

def generate_consensus(families, f_size, ref_seq, contig, region_start, region_end, bam_file, config_file):

    ## Keys: each base position in the region
    ## Values: tables of A,T,C,G (etc) counts from each UMI+Pos family
    consensus_seq = get_consensus_seq(families, contig, region_start, region_end, bam_file, config_file)

    ## Build output
    output_file = output_path + "/{}:{}-{}.fsize{}.cons".format(contig, region_start, region_end, f_size)

    with open(output_file, "w") as cons_writer:
        for base_pos in range(region_start, region_end):

            ref_base = ref_seq[base_pos-region_start]

            if base_pos in consensus_seq:

                consensuses = {}
                min_fam     = max([sum(consensus_seq[base_pos][fam].values()) for fam in consensus_seq[base_pos]]) ## start with the largest
            
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
                    
                        if sum(consensus_seq[base_pos][family].values()) < min_fam:
                            min_fam = sum(consensus_seq[base_pos][family].values())
                    
                cons_depth = sum(consensuses.values())
                n_fam      = len(consensus_seq[base_pos])
                ref_freq   = (consensuses[ref_base] / sum(consensuses.values())) * 100 if ref_base in consensuses else 0
            
                get_cons = lambda base: consensuses[base] if base in consensuses else 0

                cons_writer.write("{}\t{}\t{}\t".format(contig, base_pos, ref_base))
                cons_writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(get_cons('A'), get_cons('C'), get_cons('G'), get_cons('T'), get_cons('I'), get_cons('D'), get_cons('N')))
                cons_writer.write("{}\t{}\t{}\t{}".format(cons_depth, n_fam, min_fam, ref_freq))
                cons_writer.write("\n")

            else:
                if config['REPORT']['keep_missing_pos'] == 'TRUE':
                    cons_writer.write("{}\t{}\t{}\tMissing\n".format(contig, base_pos, ref_base))



generate_consensus(0, 0, ref_seq, contig, region_start, region_end, bam_file, config_file)

for f_size in f_sizes:
    try:
        generate_consensus(families[f_size], f_size, ref_seq, contig, region_start, region_end, bam_file, config_file)
        
    except:
        print("f_size " + str(f_size) + " not present!", file=sys.stderr)
        pass

    
lines = []
for f_size in f_sizes:
    try:
        lines.extend(open("{}/{}:{}-{}.fsize{}.cons".format(output_path, contig, region_start, region_end, f_size), "r").readlines())
    except FileNotFoundError:
        pass
        
with open("{}/{}:{}-{}.cons".format(output_path, contig, region_start, region_end), "w") as writer:
    
    for base_pos in range(region_start, region_end):
        
        for line in lines:
            if str(base_pos) in line and "Missing" not in line:
                writer.write(line)
            
    
    