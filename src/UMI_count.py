import os as path
import sys
import operator
import pysam 

output = sys.argv[2]
contig = sys.argv[3]
start  = int(sys.argv[4])
end    = int(sys.argv[5])


## Returns tally of UMIs in given region (includes mates)
def UMI_count(start, end):

    umis = {}
    
    bam_reader = pysam.AlignmentFile(sys.argv[1], "rb") 
    
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


result = UMI_count(start, end)

with open(output + "/output_{}-{}-{}.txt".format(contig, start, end), "w") as out:
    out.write("UMI       \tPosn                \tCount\n")

    for umi, count in sorted(result.items(), key=operator.itemgetter(1), reverse=True):
        out.write("{}\t{}\t{}\n".format(umi.split(",")[0], umi.split(",")[1], count))

