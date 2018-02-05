import os as path
import sys

bed_file    = sys.argv[1]
output_path = sys.argv[2]
region      = sys.argv[3]

contig       = region.split(":")[0]
region_start = region.split(":")[1].split("-")[0]
region_end   = region.split(":")[1].split("-")[1]

def get_bed_regions(bed_file, output_path):

    regions = []

    with open(sys.argv[1], "r") as bed_reader:

        for line in bed_reader.readlines()[2:]:

            contig, start, end, *rest = line.split("\t")
            regions.append([contig, start, end])


    with open("{}/{}:{}-{}.regions".format(output_path, contig, region_start, region_end), "w") as output:

        for region in regions:
            output.write("{},{},{}\n".format(region[0], region[1], region[2]))


get_bed_regions(bed_file, output_path)
