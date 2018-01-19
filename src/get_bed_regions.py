import os as path
import sys

regions = []

with open(sys.argv[1], "r") as bed_reader:

    for line in bed_reader.readlines()[2:]:

        contig, start, end, *rest = line.split("\t")
        regions.append([contig, start, end])


with open("bed_regions.csv", "w") as output:

    for region in regions:

        output.write("{},{},{}\n".format(region[0], region[1], region[2]))


