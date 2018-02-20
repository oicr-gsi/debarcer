
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


def get_consensus_seq(families, contig, region_start, region_end, bam_file, config_file):
    
    consensus_seq = {}
    
    config = configparser.ConfigParser()
    config.read(config_file)
    
    ## read_threshold = int(config['SETTINGS']['max_family_size']) if config else 50 ## Need to reconsider this
    
    with pysam.AlignmentFile(bam_file, "rb") as reader:
        
        for pileupcolumn in reader.pileup(contig, region_start, region_end, max_depth=1000000):
            
            pos = pileupcolumn.reference_pos
            
            for read in pileupcolumn.pileups:
                
                read_data = read.alignment
                read_name = read_data.query_name
                start     = read_data.reference_start
                end       = read_data.reference_end
                
                umi        = read_name.split(":")[-1]
                family_key = umi + str(start) + '-' + str(end)
                                
                if families == 0 or family_key in families:
                    
                    if not read.is_del:
                        base = read_data.query_sequence[read.query_position]
                        add_base(consensus_seq, pos, family_key, base)
                        
                    if not read.indel == 0: ## *next* position(s) are an indel
                        for i in range(pos + 1, pos + abs(read.indel) + 1):
                            add_base(consensus_seq, i, family_key, 'I' if read.indel > 0 else 'D')                                                

    return consensus_seq
