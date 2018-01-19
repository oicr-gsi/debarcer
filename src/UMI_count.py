import os as path
import sys
import operator
import pysam 

output = sys.argv[2]
contig = sys.argv[3]
start  = int(sys.argv[4])
end    = int(sys.argv[5])
bed    = sys.argv[6]

with open(bed, "r") as bed_reader:
    regions = bed_reader.readlines()

## Returns TRUE if given region matches
## a region defined in the BED file.
def isBedRegion(posA, posB):
    
    for region in regions:
    
        values = region.split(",")

        if(values[0] == contig and int(values[1]) == posA and int(values[2]) == posB):
            return(True)

    return(False)


## Returns tally of UMIs in given region (includes mates)
def UMI_count(start, end):

    umis = {}
    
    bam_reader = pysam.AlignmentFile(sys.argv[1], "rb") 
    
    for read in bam_reader.fetch(contig, start, end):

        umi  = str(read).split("HaloplexHS-")[1][:10]
        posA = read.reference_start
        rlen = read.infer_read_length()
        
        # TODO weird edge case here, re-examine later
        posB = posA + rlen if rlen else posA

        current_umi = "{},{}-{},{}".format(umi, str(posA), str(posB), isBedRegion(posA, posB))
            
        if current_umi not in umis:
            umis[current_umi] = 1
        else:
            umis[current_umi] += 1
    
    bam_reader.close()

    return(umis)


result = UMI_count(start, end)

with open(output + "/output_{}-{}-{}.txt".format(contig, start, end), "w") as out:
    out.write("UMI       \tPosn                \tCount\tisBedRegion?\n")

    for umi, count in sorted(result.items(), key=operator.itemgetter(1), reverse=True):
        out.write("{}\t{}\t{}\t{}\n".format(umi.split(",")[0], umi.split(",")[1], count, umi.split(",")[2]))

