
import sys
import configparser

region_file = sys.argv[1]
output_path = sys.argv[2]
config_file = sys.argv[3]

config = configparser.ConfigParser()
config.read(config_file)

with open(region_file, "r") as reader:
    lines = reader.readlines()

variants = []
    
for line in lines:
    
    values   = line.split(" ")
    base_pos = values[0].strip(":")

    if values[1] != "Missing\n":

        ref_base = values[2].strip("):")
    
        bases = {}
        total = 0
    
        for base_count in values[3:-1]:

            base, count = base_count.split("=")
            bases[base] = int(count)
            total += int(count)
    
        percent_ref =  (bases[ref_base] / total) * 100
        percent_ref_threshold = float(config['REPORT']['percent_ref_threshold']
    
        if percent_ref < percent_ref_threshold:
            variants.append(line)

with open(output_path + "/report.txt", "w") as writer:

    writer.write("Variants:\n")

    for variant in variants:
        writer.write(variant)
    
    
    
