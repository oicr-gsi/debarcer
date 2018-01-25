
import sys
import pysam
import operator
from collections import OrderedDict

bam_file    = sys.argv[1]
ref_file    = sys.argv[2]
region_file = sys.argv[3]

name = region_file.split("output_")[1]
name = name.strip(".txt")

contig, region_start, region_end = name.split("-")
region_start = int(region_start)
region_end   = int(region_end)

## Keys: each base position in the region
## Values: tables of A,T,C,G (etc) counts from each UMI+Pos family
region_sequence = {}

## Lists of UMI+Pos pairs with count >= cutoff
families = {}
cutoff   = 3

with open(region_file, "r") as reader:
    lines = reader.readlines()
    
for line in lines[1:]:

    umi, pos, count, *rest = line.split("\t")

    if(int(count) >= cutoff):
        families[umi + pos] = 0


def get_reference():

    with pysam.FastaFile(ref_file) as reader:
        return(reader.fetch(contig, region_start, region_end).upper())


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

    read_threshold = 50 ## max number of reads to consider in a family before finding consensus

    for read in reader.fetch(contig, region_start, region_end):
        
        umi   = read.query_name.split("HaloplexHS-")[1][:10]
        start = read.reference_start
        end   = start + len(str(read).split("\t")[9])
        pos   = "{}-{}".format(str(start), str(end))
        cigar = read.cigartuples if read.cigartuples else False

        sequence = read.query_alignment_sequence

        family_key = umi + pos
        
        if cigar and family_key in families and families[family_key] <= read_threshold:

            families[family_key] += 1

            ## Add read bases to region_sequence

            reg_pos = start ## position within the region
            seq_pos = 0     ## position within the sequence
            
            overwrite_next_match = False
            overwrite_length     = 0
            
            for tuple in cigar:

                op     = tuple[0]
                op_len = tuple[1]

                if overwrite_next_match:
                    op_len -= overwrite_length

                overwrite_next_match = False
                
                if op == 0: ## Match

                    for base_pos in range(reg_pos, reg_pos + op_len):

                        base  = sequence[seq_pos]
                        add_base(base_pos, family_key, base)
                        
                        seq_pos += 1

                    reg_pos += op_len

                elif op == 1: ## Insertion
                    
                    ## Do something with this, for the future
                    base = sequence[seq_pos : seq_pos + op_len + 1]
                    add_base(reg_pos, family_key, "I")
                    reg_pos += 1 # Insertion only counts as one "base"
                    seq_pos  += op_len + 1

                    overwrite_next_match = True ## Since we have insert, next match is shortened
                    overwrite_length     = op_len

                elif op == 2: ## Deletion

                    for base_pos in range(reg_pos, reg_pos + op_len):
                        add_base(base_pos, family_key, "-")

                    reg_pos += op_len

                elif op == 3: ## CREF_SKIP ("N")

                    for base_pos in range(reg_pos, reg_pos + op_len):
                        add_base(base_pos, family_key, "N")
                    
                    reg_pos += op_len 

                elif op == 4: ## Soft clip do nothing
                    reg_pos = reg_pos

                elif op == 5: ## Hard clip
                    ## do nothing I guess?
                    reg_pos = curr_pos

                else:
                    raise ValueError("Unknown operation in CIGAR tuple.")


ref_seq = get_reference()

writer = open(region_file + ".cons", "w")

for base_pos in range(region_start, region_end):

    if base_pos in region_sequence:

        consensuses = {} ## each family's consensus base (incl. indels, "N"s)

        for family in region_sequence[base_pos]:
        
            cons_base    = max(region_sequence[base_pos][family].items(), key = operator.itemgetter(1))[0]
            cons_denom   = sum(region_sequence[base_pos][family].values())
            cons_percent = (region_sequence[base_pos][family][cons_base]/cons_denom) * 100
    
            percent_threshold = 70.0 ## How "confident" a family must be in the consensus
            count_threshold   = 3    ## How large the consensus must be

            if cons_percent >= percent_threshold and region_sequence[base_pos][family][cons_base] >= count_threshold:
            
                if cons_base in consensuses:
                    consensuses[cons_base] += 1

                else:
                    consensuses[cons_base] = 1

            else:

                if "U" in consensuses: ## "U" (unclear) - placeholder for programmatically ambiguous bases, revisit
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

