
import sys
import configparser

def get_variants(tally_file, config_file):
    
    config = configparser.ConfigParser()
    config.read(config_file)

    with open(tally_file, "r") as reader:
        lines = reader.readlines()

    variants = []
    variants.append("-- Variants --\n")
    
    for line in lines:
    
        values   = line.split("\t")
        base_pos = values[1]

        if values[3] != "Missing\n":

            ref_base = values[2]
    
            bases    = {}
            base_tbl = ['A', 'C', 'G', 'T', 'I', 'D', 'N']
            total    = 0
    
            for base_count, idx in enumerate(values[3:10]):
                bases[idx] = int(count)
                total     += int(base_count)
    
            percent_ref =  (bases[ref_base] / total) * 100 if ref_base in bases else 0
            percent_ref_threshold = float(config['REPORT']['percent_ref_threshold'])
                                          
            if percent_ref < percent_ref_threshold:
                variants.append(line)

    return variants
    
    
    
