import os as path
import sys
import operator
import pysam 

output_path = sys.argv[1]
bam_file    = sys.argv[2]
bed_regions = sys.argv[3]
region      = sys.argv[4]

contig       = region.split(":")[0]
region_start = region.split(":")[1].split("-")[0]
region_end   = region.split(":")[1].split("-")[1]

## Returns tally of UMIs in given region (includes mates)
def UMI_count(contig, start, end):

    umis = {}
    
    bam_reader = pysam.AlignmentFile(bam_file, "rb") 
    
    for read in bam_reader.fetch(contig, start, end):

        umi  = str(read).split("HaloplexHS-")[1][:10]
        posA = read.reference_start
        rlen = len(str(read).split("\t")[9])
        
        posB = posA + rlen

        current_umi = "{},{}-{}".format(umi, str(posA), str(posB))
            
        if current_umi not in umis:
            umis[current_umi] = 1
        else:
            umis[current_umi] += 1
    
    bam_reader.close()

    return(umis)


result = UMI_count(contig, int(region_start), int(region_end))

with open("{}/{}:{}-{}.tally".format(output_path, contig, region_start, region_end), "w") as out:
    out.write("UMI       \tPosn                \tCount\n")

    for umi, count in sorted(result.items(), key=operator.itemgetter(1), reverse=True):
        out.write("{}\t{}\t{}\n".format(umi.split(",")[0], umi.split(",")[1], count))

