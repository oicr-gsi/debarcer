import os as path
import sys

regions = []
contigs = [] 
starts  = []
ends    = []

with open(sys.argv[1], "r") as bed_reader:

    for line in bed_reader.readlines()[2:]:

        contig, start, end, *rest = line.split("\t")
        contigs.append(contig)
        starts.append(int(start))
        ends.append(int(end))


curr_start  = starts[0]
curr_end    = ends[0]
curr_contig = contigs[0]

for curr_amp, start in enumerate(starts):

    if(curr_end + 20 < start or curr_contig != contigs[curr_amp]):
        regions.append([curr_contig, curr_start - 20, curr_end + 20])
        curr_start = start
        curr_contig = contigs[curr_amp]

    curr_end = ends[curr_amp]

with open("regions.csv", "w") as output:

    for region in regions:
        output.write("{},{},{}\n".format(region[0], region[1], region[2]))


