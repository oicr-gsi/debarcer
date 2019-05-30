
import pysam
import configparser


def add_base(mode, seq, pos, family, allele):
    """(Helper) Adds a base to a seq dictionary."""

    added = False
    while not added and mode == "uncollapsed":
        if pos in seq:
            if allele in seq[pos]:
                seq[pos][allele] += 1

            else:
                seq[pos][allele] = 1

            added = True

        else:
            seq[pos] = {}

    while not added and mode == "consensus":
        if pos in seq:
            if family in seq[pos]:
                if allele in seq[pos][family]:
                    seq[pos][family][allele] += 1

                else:
                    seq[pos][family][allele] = 1

                added = True

            else:
                seq[pos][family] = {}

        else:
            seq[pos] = {}
    


def get_consensus_seq(umi_table, f_size, ref_seq, contig, region_start, region_end, bam_file, config):
    """Returns consensus info for each family at each base position in the given region."""
    consensus_seq = {}

    pos_threshold = int(config['SETTINGS']['umi_family_pos_threshold']) if config else 10

    # Need to reconsider this for pileup, significant performance improvements possible
    ## read_threshold = int(config['SETTINGS']['max_family_size']) if config else 50

    region_start-=1

    with pysam.AlignmentFile(bam_file, "rb") as reader:

        for pileupcolumn in reader.pileup(contig, region_start, region_end, max_depth=1000000, ignore_orphans=False):

            #pos = pileupcolumn.reference_pos-1 
            pos = pileupcolumn.reference_pos
            
            if pos >= region_start and pos < region_end:

                for read in pileupcolumn.pileups:

                    read_data = read.alignment
                    read_name = read_data.query_name
                    start_pos = read_data.reference_start
                    umi = read_name.split(":")[-1]

                    if umi in umi_table:

                        ## Find the most proximate UMI family
                        umi_group = umi_table[umi]
                        closest_fam = umi_group.getClosest(start_pos, pos_threshold)                      
                        count = umi_group.families[closest_fam]

                        if count >= f_size:

                            family_key = umi_group.key + str(closest_fam)

                            ref_pos = pos - region_start
                    
                            if not read.is_del and not read.indel:
                                ref_base = ref_seq[ref_pos]
                                #alt_base = read_data.query_sequence[read.query_position-1]
                                alt_base = read_data.query_sequence[read.query_position]

                            # Next position is an insert (current base is ref)
                            elif read.indel > 0:
                                
                                ref_base = ref_seq[ref_pos]
                                alt_base = read_data.query_sequence[
                                    read.query_position:read.query_position + abs(read.indel)+1]

                            # Next position is a deletion (current base + next bases are ref)
                            elif read.indel < 0:
                                
                                ref_base = ref_seq[ref_pos:ref_pos + abs(read.indel) + 1]
                                alt_base = read_data.query_sequence[read.query_position]

                            if not read.is_del:
                                #add_base(mode="consensus", seq=consensus_seq, pos=pos+1,family=family_key, allele=(ref_base, alt_base))
                                add_base(mode="consensus", seq=consensus_seq, pos=pos+1,family=family_key, allele=(ref_base, alt_base))


    return consensus_seq


def get_uncollapsed_seq(ref_seq, contig, region_start, region_end, bam_file, config):
    """
    Returns a nested dictionary representing counts of each base at each base pos'n.
     - Keys: each base position in the region
     - Values: tables of A,T,C,G (etc) counts
    """
    #input(region_start)
    uncollapsed_seq = {}

    region_start-=1
    #input(region_start)
    

    with pysam.AlignmentFile(bam_file, "rb") as reader:

        ### this goes across the interval and grabs pileupcolumns at each position
        ### each pileupcolumn has a position, and a set of reads associted with it, and characteristics for each read as determined by the cigar string
        for pileupcolumn in reader.pileup(contig, region_start, region_end, max_depth=1000000, ignore_orphans=False):
            
            #pos = pileupcolumn.reference_pos - 1
            ### this is the poistion in the interval
            pos = pileupcolumn.reference_pos
            regionpos=pos-region_start
            
            
            if pos >= region_start and pos < region_end:
                #print("pileup pos is " + str(pos) + " region_start is " + str(region_start))
                #input(pos)
                #input(pileupcolumn)
                ### get each red in the pileupp colum at this interval
                readN=0
                for read in pileupcolumn.pileups:
                    readN += 1
                    #print("read is\n")
                    #input(read)
                    if not read.is_del and not read.indel:
                        
                        ### is this alway a match or a mismatch?? can it be a softclipped region?
                        ref_base = ref_seq[pos - region_start]
                        alt_base = read.alignment.query_sequence[read.query_position]
                        
                    # Next position is an insert (current base is ref)
                    elif read.indel > 0:
                        #print("read.indel>0 " + str(read.indel) + " readN=" + str(readN))
                        #print("position: chromosome:" + str(pos) + " region:" + str(regionpos) + " read:" + str(read.query_position))
                        #print("refseq: " + ref_seq)
                        #print("read: " + read.alignment.query_sequence)
                        #print("readname: " + read.alignment.query_name)
                        #print("cigarstring: " + read.alignment.cigarstring)
                        
                        
                        #ref_base = ref_seq[pos - region_start]
                        ref_pos=pos-region_start+1
                        ref_base = ref_seq[ref_pos]
                        alt_base = read.alignment.query_sequence[read.query_position + 1 :read.query_position + abs(read.indel) + 2]
                        #input(ref_base + " " + alt_base)
                    # Next position is a deletion (current base + next bases are ref)
                    elif read.indel < 0:
                        ### a deletion at this position
                        #print("read.indel<0 " + str(read.indel) + " readN=" + str(readN))
                        #print("position: chromosome:" + str(pos) + " region:" + str(regionpos) + " read:" + str(read.query_position))
                        #print("refseq: " + ref_seq)
                        #print("read: " + read.alignment.query_sequence)
                        #print("readname: " + read.alignment.query_name)
                        #print("cigarstring: " + read.alignment.cigarstring)
                        
                        
                        ### THIS IS WRONG, A was using read.query_position on the ref_seq, they are not the same
                        #ref_base = ref_seq[read.query_position+1:(read.query_position + abs(read.indel)) + 5]
                        ### THIS IS CORRECT
                        ref_pos=pos-region_start+1
                        ref_base = ref_seq[ref_pos:ref_pos +abs(read.indel)+1]
                        
                        #ref_base = ref_seq[read.query_position:(read.query_position + abs(read.indel) + 1)]
                        #alt_base = read.alignment.query_sequence[read.query_position+1]
                        alt_base = read.alignment.query_sequence[read.query_position+1 :read.query_position+abs(read.indel)+1]
                        #input(ref_base + " " + alt_base)
                    
                    
                    
                    #add_base(mode="uncollapsed", seq=uncollapsed_seq, pos=pos+1, family=None, allele=(ref_base, alt_base))
                    ## my correction, below
                    
                    add_base(mode="uncollapsed", seq=uncollapsed_seq, pos=pos+1, family=None, allele=(ref_base, alt_base))




    return uncollapsed_seq
