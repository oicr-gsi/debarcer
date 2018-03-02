
import pysam
import configparser

def add_base(mode, seq, pos, family, base):
    """(Helper) Adds a base to a seq dictionary."""

    added = False
    
    while not added and mode == "uncollapsed":
        if pos in seq:
            if base in seq[pos]:
                seq[pos][base] += 1
            
            else:
                seq[pos][base] = 1
                
            added = True
        
        else:
            seq[pos] = {}
    
    while not added and mode == "consensus":
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
    """
    Returns a nested dictionary representing counts of each base in each family at each base pos'n.
     - Keys: each base position in the region
     - Values: tables of A,T,C,G (etc) counts from each UMI+Pos family
    """
    
    consensus_seq = {}
    
    config = configparser.ConfigParser()
    config.read(config_file)
    
    ## Need to reconsider this
    ## read_threshold = int(config['SETTINGS']['max_family_size']) if config else 50
    
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
        
                if family_key in families:
                    
                    if not read.is_del:
                        base = read_data.query_sequence[read.query_position]
                        add_base("consensus", consensus_seq, pos, family_key, base)
                    else:
                        add_base("consensus", consensus_seq, pos, None, "D")
                        
                    ## Next position is an insert
                    if read.indel > 0: 
                        add_base("consensus", consensus_seq, pos + 1, family_key, 'I')
                        ## TODO code goes here for identifying the sequence of these inserts
                        
    return consensus_seq


def get_uncollapsed_seq(contig, region_start, region_end, bam_file, config_file):
    """
    Returns a nested dictionary representing counts of each base at each base pos'n.
     - Keys: each base position in the region
     - Values: tables of A,T,C,G (etc) counts
    """
    
    uncollapsed_seq = {}
    
    config = configparser.ConfigParser()
    config.read(config_file)
    
    with pysam.AlignmentFile(bam_file, "rb") as reader:
        
        for pileupcolumn in reader.pileup(contig, region_start, region_end, max_depth=1000000):
            
            pos = pileupcolumn.reference_pos
            
            for read in pileupcolumn.pileups:
                
                if not read.is_del:
                    base = read.alignment.query_sequence[read.query_position]
                    add_base("uncollapsed", uncollapsed_seq, pos, None, base)
                else:
                    add_base("uncollapsed", uncollapsed_seq, pos, None, "D")
                    
                ## Next position is an insert
                if read.indel > 0:
                    add_base("uncollapsed", uncollapsed_seq, pos + 1, None, 'I')
                    ## TODO code goes here for identifying the sequence of these inserts
                    
                        
    return uncollapsed_seq
