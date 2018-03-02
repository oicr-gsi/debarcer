
import sys
import pysam
import configparser
import argparse
import operator
from get_consensus_seq import get_consensus_seq, get_uncollapsed_seq

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
    """Generates a consensus output file (.cons) for the given family size and region."""

    ## Keys: each base position in the region
    ## Values: tables of A,T,C,G (etc) counts from each UMI+Pos family
    consensus_seq = get_consensus_seq(families, contig, region_start, region_end, bam_file, config_file)

    ## Build output
    output_file = "{}/{}:{}-{}.fsize{}.cons".format(output_path, contig, region_start, region_end, f_size)

    with open(output_file, "w") as cons_writer:
        for base_pos in range(region_start, region_end):

            ref_base = ref_seq[base_pos-region_start]

            if base_pos in consensus_seq:

                consensuses = {}
                cons_depth  = 0
                #min_fam     = max([sum(consensus_seq[base_pos][fam].values()) for fam in consensus_seq[base_pos]]) ## start with the largest
            
                for family in consensus_seq[base_pos]:
        
                    cons_base    = max(consensus_seq[base_pos][family].items(), key = operator.itemgetter(1))[0]
                    cons_denom   = sum(consensus_seq[base_pos][family].values())
                    cons_percent = (consensus_seq[base_pos][family][cons_base]/cons_denom) * 100
    
                    percent_threshold = float(config['SETTINGS']['percent_consensus_threshold']) if config else 70.0
                    count_threshold   = int(config['SETTINGS']['count_consensus_threshold']) if config else 3
    
                    if cons_percent >= percent_threshold and consensus_seq[base_pos][family][cons_base] >= count_threshold:
                        
                        ## Add the consensus base reads to consensus depth
                        cons_depth += consensus_seq[base_pos][family][cons_base] 
                        
                        if cons_base in consensuses:
                            consensuses[cons_base] += 1

                        else:
                            consensuses[cons_base] = 1
                        
                        ## min_fam logic, may not be worth reporting
                        #if sum(consensus_seq[base_pos][family].values()) < min_fam:
                        #    min_fam = sum(consensus_seq[base_pos][family].values())
                    
                n_fam      = len(consensus_seq[base_pos])
                ref_freq   = (consensuses[ref_base] / n_fam) * 100 if ref_base in consensuses else 0
            
                get_cons = lambda base: consensuses[base] if base in consensuses else 0

                cons_writer.write("{}\t{}\t{}\t".format(contig, base_pos, ref_base))
                cons_writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(get_cons('A'), get_cons('C'), get_cons('G'), get_cons('T'), get_cons('I'), get_cons('D'), get_cons('N')))
                cons_writer.write("{}\t{}\t{}\t{}\n".format(cons_depth, n_fam, f_size, ref_freq))

            else:
                if config['REPORT']['keep_missing_pos'] == 'TRUE':
                    cons_writer.write("{}\t{}\t{}\tMissing\n".format(contig, base_pos, ref_base))


def generate_uncollapsed(ref_seq, contig, region_start, region_end, bam_file, config_file):
    """Generates an uncollapsed consensus output file (.fsize0.cons) for the given family size and region."""
    
    ## Keys: each base position in the region
    ## Values: tables of A,T,C,G (etc) counts from each UMI+Pos family
    uncollapsed_seq = get_uncollapsed_seq(contig, region_start, region_end, bam_file, config_file)
    
    ## Build output
    output_file = "{}/{}:{}-{}.fsize0.cons".format(output_path, contig, region_start, region_end)
    
    with open(output_file, "w") as seq_writer:
        for base_pos in range(region_start, region_end):
            
            ref_base = ref_seq[base_pos-region_start]
            
            if base_pos in uncollapsed_seq:
                
                get_cons = lambda base: uncollapsed_seq[base_pos][base] if base in uncollapsed_seq[base_pos] else 0
                depth    = sum(uncollapsed_seq[base_pos].values())
                ref_freq = (uncollapsed_seq[base_pos][ref_base] / depth) * 100 if ref_base in uncollapsed_seq[base_pos] else 0

                seq_writer.write("{}\t{}\t{}\t".format(contig, base_pos, ref_base))
                seq_writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(get_cons('A'), get_cons('C'), get_cons('G'), get_cons('T'), get_cons('I'), get_cons('D'), get_cons('N')))
                seq_writer.write("{}\t{}\t{}\t{}\n".format(depth, "N/A", 0, ref_freq))
    
            else:
                if config['REPORT']['keep_missing_pos'] == 'TRUE':
                    seq_writer.write("{}\t{}\t{}\tMissing\n".format(contig, base_pos, ref_base))
    

generate_uncollapsed(ref_seq, contig, region_start, region_end, bam_file, config_file)


for f_size in f_sizes:
    try:
        generate_consensus(families[f_size], f_size, ref_seq, contig, region_start, region_end, bam_file, config_file)
        
    except:
        print("f_size " + str(f_size) + " not present!", file=sys.stderr)
        pass
    
cons_files = []
cons_files.append( [open("{}/{}:{}-{}.fsize0.cons".format(output_path, contig, region_start, region_end), "r").readlines(), 0] )

for f_size in f_sizes:
    try:
        cons_files.append( [open("{}/{}:{}-{}.fsize{}.cons".format(output_path, contig, region_start, region_end, f_size), "r").readlines(), 0] )
    except FileNotFoundError:
        pass
        
with open("{}/{}:{}-{}.cons".format(output_path, contig, region_start, region_end), "w") as writer:
    
    for base_pos in range(region_start, region_end):
        
        for idx, cons_file in enumerate(cons_files):
            
            counter = cons_file[1]
            if counter < len(cons_file[0]):
                
                line = cons_file[0][counter]
            
                if str(base_pos) in line:
                
                    if "Missing" not in line or config['REPORT']['keep_missing_pos'] == 'TRUE':
                        writer.write(line)
                    
                    cons_files[idx][1] += 1 ## increment counter
                
            
            
            
    
    