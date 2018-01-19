
import sys
import pysam
import operator
from collections import OrderedDict

bam_file    = sys.argv[1]
region_file = sys.argv[2]

name = region_file.split("output_")[1]
name = name.strip(".txt")

contig, region_start, region_end = name.split("-")
region_start = int(region_start)
region_end   = int(region_end)

## Keys: each base position in the region
## Values: tables of A,T,C,G (etc) counts from each UMI+Pos family
region_sequence = {}


## Lists of UMI+Pos pairs with count >= cutoff
families = set()
cutoff   = 15

with open(region_file, "r") as reader:
    lines = reader.readlines()
    
for line in lines[1:]:

    umi, pos, count, rest = line.split("\t")

    if(int(count) >= cutoff):
        families.add(umi + pos)


def add_base(pos, family, base):

    added = False

    while not added:
                            
        if base_pos in region_sequence:

            if family_key in region_sequence[base_pos]:

                if base in region_sequence[base_pos][family_key]:
                    region_sequence[base_pos][family_key][base] += 1

                else:
                    region_sequence[base_pos][family_key][base] = 1

                added = True

            else:
                region_sequence[base_pos][family_key] = {}

        else:
            region_sequence[base_pos] = {}


# Build the region_sequence table
with pysam.AlignmentFile(bam_file, "rb") as reader:

    for read in reader.fetch(contig, region_start, region_end):
        
        umi   = read.query_name.split("HaloplexHS-")[1][:10]
        start = read.reference_start
        end   = start + read.infer_read_length() if read.infer_read_length() else start
        pos   = "{}-{}".format(str(start), str(end))
        cigar = read.cigartuples if read.cigartuples else False

        sequence = read.query_sequence

        family_key = umi + pos
        
        if family_key in families and cigar:

            ## TODO: Need to incorporate CIGAR logic
            ## Notes:
            ## Sequence should usually be 150 long, DOES NOT REMOVE soft clips from start and end
            
            curr_pos = start ## we "start" at start
            seq_pos  = 0     ## position within the sequence

            #stored = "" ## TEST

            for tuple in cigar:

                op     = tuple[0]
                op_len = tuple[1]

                if op in [0, 3]: ## Match or CREF_SKIP ("N")

                    for base_pos in range(curr_pos, curr_pos + op_len):

                        base  = sequence[seq_pos]
                        
                        add_base(base_pos, family_key, base)
                        #stored += base ## TEST

                        seq_pos += 1

                    curr_pos += op_len

                elif op == 1: ## Insertion

                    base = sequence[seq_pos : seq_pos + op_len]
                    
                    add_base(base_pos, family_key, base)
                    #stored += base ## TEST

                    curr_pos += 1 # Insertion only counts as one "base"
                    seq_pos  += op_len

                elif op == 2: ## Deletion

                    base = "-"

                    for base_pos in range(curr_pos, curr_pos + op_len):

                        #print("Deletion\n") ## TEST

                        add_base(base_pos, family_key, base)
                        #stored += base ## TEST

                    curr_pos += op_len

                elif op == 4: ## Soft clip
                    curr_pos += op_len
                    seq_pos  += op_len

                elif op == 5: ## Hard clip
                    ## do nothing I guess?
                    curr_pos = curr_pos

                else:
                    curr_pos += op_len
                    seq_pos  += op_len

            #memory = ""
            #for base_pos in range(start, end):
            #    for letter in region_sequence[base_pos][family_key]:
            #        memory += letter
            #        memory += str(region_sequence[base_pos][family_key][letter])
            #    memory += "~"

            #print("Seq: {}\nStored: {}\nMemory: {}\n".format(sequence, stored, memory))

            ## END TODO
 
# Generate consensus sequence and statistics
## /TODO/

#TEST
print(region_sequence[36164450])


read_cutoff = 5 ## min. number of consensus reads required in a family 

writer = open(region_file + ".cons", "w")

for base_pos in range(region_start, region_end):

    if base_pos in region_sequence:

        consensuses = {} ## each family's consensus base (incl. indels, "N"s)

        for family in region_sequence[base_pos]:
        
            base = max(region_sequence[base_pos][family].items(), key = operator.itemgetter(1))[0]
            
            if region_sequence[base_pos][family][base] >= read_cutoff:
            
                if base in consensuses:
                    consensuses[base] += 1

                else:
                    consensuses[base] = 1

        writer.write("{}: ".format(base_pos + 1))
        
        for base, count in sorted(consensuses.items(), key = operator.itemgetter(1), reverse=True):
            writer.write("{}={} ".format(base, str(count)))

        writer.write("\n")

    else:

        writer.write("{}: Missing\n".format(base_pos))

writer.close()

