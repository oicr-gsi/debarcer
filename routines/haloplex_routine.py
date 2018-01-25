
import pysam
import configparser

def add_base(seq, pos, family, base):

    added = False

    while not added:
                            
        if pos in seq:

            if family in seq[pos]:

                if base in seq[pos][family]:
                    seq[pos][family][base] += 1

                else:
                    seq[pos][family][base] = 1

                added = True

            else:
                seq[pos][family] = {}

        else:
            seq[pos] = {}

def routine(families, contig, region_start, region_end, bam_file, config_file):
    
    consensus_seq = {}
    
    config = configparser.ConfigParser()
    config.read(config_file)
    
    with pysam.AlignmentFile(bam_file, "rb") as reader:
        
        read_threshold = int(config['SETTINGS']['max_family_size'])
        
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
                            add_base(consensus_seq, base_pos, family_key, base)
                                
                            seq_pos += 1
                                
                        reg_pos += op_len
                                
                    elif op == 1: ## Insertion
                                    
                        ## Do something with this, for the future
                        base = sequence[seq_pos : seq_pos + op_len + 1]
                        add_base(consensus_seq, reg_pos, family_key, "I")
                        reg_pos += 1 # Insertion only counts as one "base"
                        seq_pos  += op_len + 1
                                    
                        overwrite_next_match = True ## Since we have insert, next match is shortened
                        overwrite_length     = op_len
                                    
                    elif op == 2: ## Deletion
                                        
                        for base_pos in range(reg_pos, reg_pos + op_len):
                            add_base(consensus_seq, base_pos, family_key, "-")
                                    
                        reg_pos += op_len
                                            
                    elif op == 3: ## CREF_SKIP ("N")
                                                
                        for base_pos in range(reg_pos, reg_pos + op_len):
                            add_base(consensus_seq, base_pos, family_key, "N")
                                                    
                        reg_pos += op_len 
                                                    
                    elif op == 4: ## Soft clip do nothing (?)
                        reg_pos = reg_pos
                                                    
                    elif op == 5: ## Hard clip do nothing (?)
                        reg_pos = reg_pos
                                                            
                    else:
                        raise ValueError("Unknown operation in CIGAR tuple.")

    return consensus_seq
