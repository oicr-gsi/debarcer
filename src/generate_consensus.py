
import sys
import pysam
import configparser
import importlib.util
import operator

bam_file    = sys.argv[1]
region_file = sys.argv[2]
config_file = sys.argv[3]
region      = sys.argv[4]

config = configparser.ConfigParser()
config.read(config_file)

contig       = region.split(":")[0]
region_start = int(region.split(":")[1].split("-")[0])
region_end   = int(region.split(":")[1].split("-")[1])

routine_file = config['PATHS']['routine_file']
spec         = importlib.util.spec_from_file_location("routine", routine_file)
routine      = importlib.util.module_from_spec(spec)
spec.loader.exec_module(routine)

def get_reference():

    with pysam.FastaFile(config['PATHS']['reference_file']) as reader:
        return(reader.fetch(contig, region_start, region_end).upper())

ref_seq = get_reference()


## Lists of UMI+Pos pairs with count >= cutoff
families = {}
cutoff   = int(config['SETTINGS']['min_family_size'])

with open(region_file, "r") as reader:
    lines = reader.readlines()
    
for line in lines[1:]:

    umi, pos, count, *rest = line.split("\t")

    if(int(count) >= cutoff):
        families[umi + pos] = 0
        

## Keys: each base position in the region
## Values: tables of A,T,C,G (etc) counts from each UMI+Pos family
consensus_seq = routine.routine(families, contig, region_start, region_end, bam_file, config_file)

## Build output
output_file = config['PATHS']['output_folder'] + "/{}:{}-{}.cons".format(contig, region_start, region_end)
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
    
            percent_threshold = float(config['SETTINGS']['percent_consensus_threshold'])
            count_threshold   = int(config['SETTINGS']['count_consensus_threshold'])

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

