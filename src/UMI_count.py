
import os as path
import sys
import operator
import argparse
import configparser
import pysam 


## Argument parsing and error handling
def handle_arg(var, alt, config, error):
    
    if var is None:
 
        var = alt
        
        if var is None:
            raise ValueError(error)
            sys.exit(1)
        
    return var
            

parser = argparse.ArgumentParser()
parser.add_argument('-bam', '--bam_file',    help='Path to your BAM file.')
parser.add_argument('-bed', '--bed_file',    help='Path to your BED file.')
parser.add_argument('-r',   '--region',      help='Region to tally (string of the form chrX:posA-posB).')
parser.add_argument('-o',   '--output_path', help='Path to write output files to.')
parser.add_argument('-c',   '--config',      help='Path to your config file.')

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

bam_file    = handle_arg(args.bam_file, config['PATHS']['bam_file'] if config else None, config, 'No BAM file provided in args or config.')
bed_file    = handle_arg(args.bed_file, config['PATHS']['bed_file'] if config else None, config, 'No BED file provided in args or config.')
output_path = handle_arg(args.output_path, config['PATHS']['output_path'] if config else None, config, 'No output path provided in args or config.')
    

## Parses BED file into a CSV of positions
def get_bed_regions(bed_file, output_path):

    regions = []

    with open(bed_file, "r") as bed_reader:

        for line in bed_reader.readlines()[2:]:

            contig, start, end, *rest = line.split("\t")
            regions.append([contig, start, end])


    with open("{}/{}:{}-{}.regions".format(output_path, contig, region_start, region_end), "w") as output:

        for region in regions:
            output.write("{},{},{}\n".format(region[0], region[1], region[2]))
            
        return output.name
            

## Updates umi_table to include isBedRegion column
def update_bed(umi_table, amps):
    
    for u_id in umi_table:
        
        posA, posB = umi_table[u_id]['posn'].split('-')
        posA = int(posA)
        posB = int(posB)
        
        for amp in amps:
        
            bed_posA = amp[0]
            bed_posB = amp[1]
            
            if posA >= bed_posA and posB <= bed_posB:
                umi_table[u_id]['isBedRegion'] = True
                
            ## Assumes coordinate-sorted BED file, revisit this?
            if bed_posB >= posA:
                break


## Returns tally of UMIs in given region (includes mates)
def UMI_count(contig, start, end):

    umi_table = {}
    amps      = []
    
    bed_regions = get_bed_regions(bed_file, output_path)
    
    with open(bed_regions, "r") as bed_reader:
        for line in bed_reader.readlines():
            amps.append([int(line.split(",")[1]), int(line.split(",")[2])])
    
    with pysam.AlignmentFile(bam_file, "rb") as bam_reader: 
    
        for read in bam_reader.fetch(contig, start, end):

            umi  = read.query_name.split(':')[-1]
            posA = read.reference_start
            posB = read.reference_end
            posn = str(posA) + '-' + str(posB)
            
            u_id = umi + posn
            
            if posA and posB: ## Make sure no null values sneak in
                
                if u_id in umi_table:
                    umi_table[u_id]['count'] += 1
                else:
                    umi_table[u_id] = {'umi': umi, 'posn': posn, 'count': 1, 'isBedRegion': False}
    
    
    ## Build output table
    update_bed(umi_table, amps)
    posn_table = {}
    
    for u_id in umi_table:
        
        umi = umi_table[u_id]
        pos = umi_table[u_id]['posn']
        
        if pos in posn_table:
            posn_table[pos]['tcount'] += 1
            posn_table[pos]['treads'] += umi['count']
            posn_table[pos]['umi_csv'].append(umi['umi'])
            posn_table[pos]['count_csv'].append(str(umi['count']))
        else:
            posn_table[pos] = {'tcount': 1, 'treads': umi['count'], 'umi_csv': [umi['umi']], 'count_csv': [str(umi['count'])]}
            
    return umi_table, posn_table


## Write output (if script is called standalone)
umi_table, posn_table = UMI_count(contig, int(region_start), int(region_end))

with open("{}/{}:{}-{}.tally".format(output_path, contig, region_start, region_end), "w") as out:
    #out.write("UMI       \tPosn           \tCount\tisBedRegion?\n")

    for u_id in sorted(umi_table, key=lambda x: (umi_table[x]['count']), reverse=True):
        out.write("{}\t{}\t{}\t{}\t\n".format(umi_table[u_id]['umi'], umi_table[u_id]['posn'], 
                                              umi_table[u_id]['count'], umi_table[u_id]['isBedRegion']))

with open("{}/{}:{}-{}.posns".format(output_path, contig, region_start, region_end), "w") as out:
    #out.write("Posn | UMI Count | Read Count | UMIs (CSV) | Counts per UMI (CSV)\n")

    for posn in posn_table:
        out.write("{}\t{}\t{}\t{}\t{}\t\n".format(posn, posn_table[posn]['tcount'], posn_table[posn]['treads'], 
                                                  ",".join(posn_table[posn]['umi_csv']), ",".join(posn_table[posn]['count_csv'])))

