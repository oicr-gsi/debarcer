
import sys
import pysam
import configparser

def get_coverage(contig, region_start, region_end, bam_file, config_file):

    config = configparser.ConfigParser()
    config.read(config_file)
    
    coverages = []
    coverages.append("-- Coverage --\n")
    
    with pysam.AlignmentFile(bam_file, "rb") as reader:
    
        for pileupcolumn in reader.pileup(contig, int(region_start), int(region_end)):
            coverages.append("Coverage at base {}: {}\n".format(pileupcolumn.pos + 1, pileupcolumn.n))
        
    return coverages

    
    
    
