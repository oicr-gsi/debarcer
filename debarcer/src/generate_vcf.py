
import sys
import pysam
import configparser
import argparse
import operator
import functools
import csv
from src.generate_consensus import ConsDataRow
from src.get_ref_seq import get_ref_seq
from src.umi_error_correct import UMIGroup
from src.get_consensus_seq import get_consensus_seq, get_uncollapsed_seq
from src.handle_args import handle_arg
from src.get_ref_seq import get_ref_seq
from src.generate_consensus import temp_umi_table
from src.generate_consensus import generate_uncollapsed
from src.generate_consensus import generate_consensus
from src.generate_consensus import memoize
#from src.generate_merge import check_consfile


def parse_raw_table(cons_file, f_sizes):
    """Parses a .cons file generated by generate_consensus into VCF entries."""
    
    rows = {f_size: [] for f_size in f_sizes}

    with open(cons_file, "r") as reader:

        f_size = 0

        for line in reader:

            if line.startswith('#'):
                f_size = line.split('\t')[11]

            elif f_size in f_sizes:
                rows[f_size].append(line)

    return rows


def write_rows(cons_data, f_size, contig, region_start, region_end, output_path, config):

    with open("{}/{}:{}-{}.{}.vcf".format(output_path, contig, region_start, region_end, f_size), "w") as writer:

        writer.write("##fileformat=VCFv4.2\n")
        writer.write("##reference={}\n".format(config['PATHS']['reference_file']))
        writer.write("##source=Debarcer2\n")
        writer.write("##f_size={}\n".format(f_size))
        
        ## INFO/FILTER/FORMAT metadata
        writer.write("##INFO=<ID=RDP,Number=1,Type=Integer,Description=\"Raw Depth\">\n")
        writer.write("##INFO=<ID=CDP,Number=1,Type=Integer,Description=\"Consensus Depth\">\n")
        writer.write("##INFO=<ID=MIF,Number=1,Type=Integer,Description=\"Minimum Family Size\">\n")
        writer.write("##INFO=<ID=MNF,Number=1,Type=Float,Description=\"Mean Family Size\">\n")
        writer.write("##FILTER=<ID=a10,Description=\"Alt allele depth below 10\">\n")
        writer.write("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n")
        writer.write("##FORMAT=<ID=AL,Number=R,Type=Integer,Description=\"Alternate Allele Depth\">\n")
        writer.write("##FORMAT=<ID=AF,Number=R,Type=Float,Description=\"Alternate Allele Frequency\">\n")
  
        writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        
        for row in cons_data[f_size]:
            writer.write(row)


def vcf_output(cons_file, f_size, contig, region_start, region_end, output_path, config):
    """Writes a .vcf consensus file."""
        
    with open("{}/{}:{}-{}.fsize{}.vcf".format(output_path, contig, region_start, region_end, f_size), "w") as writer:
        
        writer.write("##fileformat=VCFv4.2\n")
        writer.write("##reference={}\n".format(config['PATHS']['reference_file']))
        writer.write("##source=Debarcer2\n")
        writer.write("##f_size={}\n".format(f_size))
        
        ## INFO/FILTER/FORMAT metadata
        writer.write("##INFO=<ID=RDP,Number=1,Type=Integer,Description=\"Raw Depth\">\n")
        writer.write("##INFO=<ID=CDP,Number=1,Type=Integer,Description=\"Consensus Depth\">\n")
        writer.write("##INFO=<ID=MIF,Number=1,Type=Integer,Description=\"Minimum Family Size\">\n")
        writer.write("##INFO=<ID=MNF,Number=1,Type=Float,Description=\"Mean Family Size\">\n")
        writer.write("##FILTER=<ID=a10,Description=\"Alt allele depth below 10\">\n")
        writer.write("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n")
        writer.write("##FORMAT=<ID=AL,Number=R,Type=Integer,Description=\"Alternate Allele Depth\">\n")
        writer.write("##FORMAT=<ID=AF,Number=R,Type=Float,Description=\"Alternate Allele Frequency\">\n")
        
        ref_threshold = float(config['REPORT']['percent_ref_threshold']) if config else 95.0
        all_threshold = float(config['REPORT']['percent_allele_threshold']) if config else 2.0
        
        writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        for base_pos in range(region_start, region_end):
            if base_pos in cons_data[f_size]:
                
                row = cons_data[f_size][base_pos]
                ref = row.get_ref_info()
                cons = row.get_cons_info()
                stats = row.get_stats()

                if stats['ref_freq'] <= ref_threshold:
                        
                    alleles = row.get_alleles(all_threshold)
                    ref_bases = set([allele[0] for allele in alleles])
                    ref_allele = (ref_seq[base_pos - region_start], ref_seq[base_pos - region_start])
                    depths = row.impute_allele_depths()
                    ref_depth = depths[ref_allele] if ref_allele in depths else 0
                    alt_freqs = row.impute_allele_freqs(all_threshold)

                    info = "RDP={};CDP={};MIF={};MNF={:.1f}".format(
                        stats['rawdp'], stats['consdp'], stats['min_fam'], stats['mean_fam'])
                    fmt_string = "AD:AL:AF" # Allele depth, alt allele depth, reference frequency

                    for ref_base in ref_bases:
                        snips = []
                        for allele in alleles:
                            if allele[0] == ref_base:
                                snips.append(allele)

                        alt_string = ','.join( [allele[1] for allele in snips] )
                        depth_string = ','.join( [str(depths[allele]) for allele in snips] )
                        freq_string = ','.join( ["{:.2f}".format(alt_freqs[allele]) for allele in snips] )
                        smp_string = "{}:{}:{}".format(ref_depth, depth_string, freq_string)
                        filt = "PASS" if any( [depths[alt] > 10 for alt in snips] ) else "a10"

                        writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            contig, base_pos, ".", ref_base, alt_string, "0", filt, info, fmt_string, smp_string))
    


def create_vcf_output(cons_file, f_size, cons_data, ref_seq, contig, region_start, region_end, output_path, config):
    #base_positions, ref_bases = {}, {}
    ref_threshold = float(config['REPORT']['percent_ref_threshold']) if config else 95.0
    all_threshold = float(config['REPORT']['percent_allele_threshold']) if config else 2.0

	
    #Read cons file
    with open(cons_file, "r") as reader:
        print("Reading cons file")

        #f_size = 0
        i = region_start
        for line in reader:
            ref_freq = line.split('\t')[13]

            if i <= region_end:

                #contig[i] = line.split('\t')[1] #### CHANGE: Added the following lines
                base_positions[i] = line.split('\t')[2]
                ref_bases[i] = line.split('\t')[3]
                #rawdps[i] = line.split('\t')[10]
                #consdps[i] = line.split('\t')[11]
                #f_size[i] = line.split('\t')[12] ##### CHANGE: from [11] to [12]
 
                i=+1



    with open("{}/{}:{}-{}.{}.vcf".format(output_path, contig, region_start, region_end, f_size), "w") as writer:
        
        writer.write("##fileformat=VCFv4.2\n")
        writer.write("##reference={}\n".format(config['PATHS']['reference_file']))
        writer.write("##source=Debarcer2\n")
        writer.write("##f_size={}\n".format(f_size))
        
        ## INFO/FILTER/FORMAT metadata
        writer.write("##INFO=<ID=RDP,Number=1,Type=Integer,Description=\"Raw Depth\">\n")
        writer.write("##INFO=<ID=CDP,Number=1,Type=Integer,Description=\"Consensus Depth\">\n")
        writer.write("##INFO=<ID=MIF,Number=1,Type=Integer,Description=\"Minimum Family Size\">\n")
        writer.write("##INFO=<ID=MNF,Number=1,Type=Float,Description=\"Mean Family Size\">\n")
        writer.write("##FILTER=<ID=a10,Description=\"Alt allele depth below 10\">\n")
        writer.write("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n")
        writer.write("##FORMAT=<ID=AL,Number=R,Type=Integer,Description=\"Alternate Allele Depth\">\n")
        writer.write("##FORMAT=<ID=AF,Number=R,Type=Float,Description=\"Alternate Allele Frequency\">\n")
        
        writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        for base_pos in range(region_start, region_end):
            #if any( [base_pos in cons_data[f_size] for f_size in cons_data]):
                #for f_size in cons_data:
            if base_pos in cons_data[f_size]:

                row = cons_data[f_size][base_pos]
                ref = row.get_ref_info()
                cons = row.get_cons_info()
                stats = row.get_stats()


                if stats['ref_freq'] <= ref_threshold:

                    alleles = row.get_alleles(all_threshold)
                    depths = row.impute_allele_depths()
                    ref_bases = set([allele[0] for allele in alleles])
                    ref_allele = (ref_seq[base_pos - region_start], ref_seq[base_pos - region_start])
                    depths = row.impute_allele_depths()
                    ref_depth = depths[ref_allele] if ref_allele in depths else 0
                    alt_freqs = row.impute_allele_freqs(all_threshold)


                    info = "RDP={};CDP={};MIF={};MNF={:.1f}".format(
                        stats['rawdp'], stats['consdp'], stats['min_fam'], stats['mean_fam'])
                    fmt_string = "AD:AL:AF" # Allele depth, alt allele depth, reference frequency

                    for ref_base in ref_bases:
                        snips = []
                        for allele in alleles:
                            if allele[0] == ref_base:
                                snips.append(allele)


                        alt_string = ','.join( [allele[1] for allele in snips] )
                        depth_string = ','.join( [str(depths[allele]) for allele in snips] )
                        freq_string = ','.join( ["{:.2f}".format(alt_freqs[allele]) for allele in snips] )
                        smp_string = "{}:{}:{}".format(ref_depth, depth_string, freq_string)
                        filt = "PASS" if any( [depths[alt] > 10 for alt in snips] ) else "a10"

                        writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            contig, base_pos, ".", ref_base, alt_string, "0", filt, info, fmt_string, smp_string))



        


def create_consensus_output(cons_file, f_sizes, contig, region_start, region_end, bam_file, output_path, config):
    
    print("Building temporary UMI table...")
    umi_table = temp_umi_table(contig, region_start, region_end, bam_file, config)
        
    ## Lists of umi families with count >= f_size
    f_sizes = [int(n) for n in config['SETTINGS']['min_family_sizes'].split(',')] if config else [1, 2, 5]

    ## Get reference sequence for the region 
    print("Getting reference sequence...")
    ref_seq = get_ref_seq(contig, region_start, region_end, config)

    ## Get consensus data for each f_size + uncollapsed data
    print("Building consensus data...")
    cons_data = {}
    
    
    cons_data[0] = generate_uncollapsed(ref_seq, contig, region_start, region_end, bam_file, config)

    for f_size in f_sizes:
        cons_data[f_size] = generate_consensus(umi_table, f_size, ref_seq, contig, region_start, region_end, bam_file, config)
        if f_size in cons_data:
            ## Output
            print("Writing output...")
            create_vcf_output(cons_file, f_size, cons_data, ref_seq, contig, region_start, region_end, output_path, config)





def generate_vcf_output(cons_file, f_sizes, contig, region_start, region_end, output_path, config):
    """(Main) generates VCF output file(s)."""

    ## Get reference sequence for the region 
    ref_seq = get_ref_seq(contig, region_start, region_end, config)

    ## Parse table to extract VCF events
    cons_data = parse_raw_table(cons_file, f_sizes)

    ## Generate vcf files
    for f_size in f_sizes:
        if f_size in cons_data:
            write_rows(cons_data, f_size, contig, region_start, region_end, output_path, config)




def check_consfile(cons_file):
    f = open(cons_file, "r")
    #reader = csv.reader(f, delimiter='\t')
    line = f.readline()
    if 'INTVL' in line:
        return True
    else:
        return False




def write_vcf(config, contigs, f_sizes, base_positions, region_start, region_end, ref_base, alt_string, filt, info, fmt_string, smp_string, output_path, cons_is_merged):

    ref_threshold = 100.0
    #fam_idx = 0

    if not cons_is_merged:
        for i in range(len(contigs)):
            for size in f_sizes:
                if contigs[i][size] != None:
                    contig = contigs[i][size]
                    break

    #Create vcf file for a given base position for specified min fam size
    for f_size in f_sizes:

        if cons_is_merged:
            writer = open("{}/{}_{}.{}.vcf".format(output_path, region_start, region_end, str(f_size)), "w")
        else:
            writer = open("{}/{}:{}-{}.{}.vcf".format(output_path, contig, region_start, region_end, str(f_size)), "w")

        writer.write("##fileformat=VCFv4.2\n")
        writer.write("##reference={}\n".format(config['PATHS']['reference_file']))
        writer.write("##source=Debarcer2\n")
        writer.write("##f_size={}\n".format(f_size))

        ## INFO/FILTER/FORMAT metadata
        writer.write("##INFO=<ID=RDP,Number=1,Type=Integer,Description=\"Raw Depth\">\n")
        writer.write("##INFO=<ID=CDP,Number=1,Type=Integer,Description=\"Consensus Depth\">\n")
        writer.write("##INFO=<ID=MIF,Number=1,Type=Integer,Description=\"Minimum Family Size\">\n")
        writer.write("##INFO=<ID=MNF,Number=1,Type=Float,Description=\"Mean Family Size\">\n")
        writer.write("##FILTER=<ID=a10,Description=\"Alt allele depth below 10\">\n")
        writer.write("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth (ref allele depth, alt allele depth(s))\">\n")
        writer.write("##FORMAT=<ID=AF,Number=R,Type=Float,Description=\"Allele Frequency (ref allele frequency, alt allele freq(s))\">\n")
        writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        for index in range(len(contigs)):
            if contigs[index][f_size] != None:
                writer.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                       contigs[index][f_size], base_positions[index][f_size], ".", ref_base[index][f_size], alt_string[index][f_size], "0", filt[index][f_size], info[index][f_size], fmt_string, smp_string[index][f_size]))
        #fam_idx += 1




def get_data(f_sizes, ref_base, base_A, base_C, base_G, base_T, filt, consdps, alt_string_info, allele_data, sample_data, rawdps, mean_fams, min_fam):
    for index in range(len(ref_base)):
        for f_size in f_sizes:
            if ref_base[index][f_size] != None:
                AD_string = ""
                alt_string = ""
                freq_string = ""

                alleles = {}
                ref_allele = ref_base[index][f_size]
               
                if len(ref_allele) > 1:
                    print("----REF ALLELE: "+ref_allele)               
 
                #Temporarily assigning keys (base) to value (allele depth)
                alleles['A'] = int(base_A[index][f_size])
                alleles['C'] = int(base_C[index][f_size])
                alleles['G'] = int(base_G[index][f_size])
                alleles['T'] = int(base_T[index][f_size])

                #Sort dict 'alleles' from most-frequent allele to least-frequent allele
                temp_tuple = sorted(alleles.items(), key=operator.itemgetter(1))
                sorted_alleles = [temp_tuple[3][0], temp_tuple[2][0], temp_tuple[1][0], temp_tuple[0][0]]
                sorted_alleles_depth = [int(temp_tuple[3][1]), int(temp_tuple[2][1]), int(temp_tuple[1][1]), int(temp_tuple[0][1])]

                #getting index of the reference allele, and ref depth
                #print(str(index)+" "+str(f_size)+" "+ref_allele+" INDEX: "+str(sorted_alleles.index(ref_allele))) ###################

                ref_index = int(sorted_alleles.index(ref_allele))
                ref_allele_depth = sorted_alleles_depth[ref_index]
                AD_string = str(ref_allele_depth) #Sample string: ref allele depth, alt allele depth(s)
                ref_AF = ((ref_allele_depth)/(int(consdps[index][f_size])))*100
                freq_string = "{:.2f}".format(ref_AF)

                filt[index][f_size] = "a10"

                for count in range(4):
                    #Check that the allele is not the ref allele and that it has a depth greater than 0, and append to the appropriate string
                    if count != ref_index and sorted_alleles_depth[count] != 0:

                        #AD_string consists of allele depths as: ref_allele_depth, alt_alleles_depth(s)
                        #alt_string consists of all alternate alleles at each base position, ordered from most-frequenct to least-frequent
                        #freq_string consists of the frequency that each alternate allele occurs at a given base position

                        AD_string = AD_string+","+str(sorted_alleles_depth[count])

                        if alt_string != "":
                            alt_string = alt_string+","+sorted_alleles[count]
                        else:
                            alt_string = sorted_alleles[count]

                        AF_value = ((sorted_alleles_depth[count])/(int(consdps[index][f_size])))*100

                        freq_string = freq_string+","+"{:.2f}".format(AF_value)

                        if sorted_alleles_depth[count] >= 10:
                            filt[index][f_size] = "PASS"

                alt_string_info[index][f_size] = alt_string
                allele_data[index][f_size] = "AD="+AD_string+":AF="+freq_string

                sample_data[index][f_size] = "RDP="+rawdps[index][f_size]+":CDP="+consdps[index][f_size]+":MNF="+mean_fams[index][f_size]+":"+allele_data[index][f_size]
                if f_size != None:
                    min_fam[index][f_size] = "MIF="+str(min(i for i in f_sizes if int(i) > 0))
                else: 
                    min_fam[index][f_size] = "MIF=ERROR"





def get_data2(f_sizes, ref_base, base_A, base_C, base_G, base_T, filt, consdps, alt_string_info, allele_data, sample_data, rawdps, mean_fams, min_fam, intvls):
    for index in range(len(ref_base)):      
        for f_size in f_sizes:
            if ref_base[index][f_size] != None:
                AD_string = ""
                alt_string = ""
                freq_string = ""

                alleles = {}
                ref_allele = ref_base[index][f_size]
		
                #Temporarily assigning keys (base) to value (allele depth)
                alleles['A'] = int(base_A[index][f_size])
                alleles['C'] = int(base_C[index][f_size])
                alleles['G'] = int(base_G[index][f_size])
                alleles['T'] = int(base_T[index][f_size])

                #Sort dict 'alleles' from most-frequent allele to least-frequent allele
                temp_tuple = sorted(alleles.items(), key=operator.itemgetter(1))                          
                sorted_alleles = [temp_tuple[3][0], temp_tuple[2][0], temp_tuple[1][0], temp_tuple[0][0]]
                sorted_alleles_depth = [int(temp_tuple[3][1]), int(temp_tuple[2][1]), int(temp_tuple[1][1]), int(temp_tuple[0][1])]    
             
                #getting index of the reference allele, and ref depth
                #print(str(index)+" "+str(f_size)+" "+ref_allele+" INDEX: "+str(sorted_alleles.index(ref_allele))) ###################

                ref_index = int(sorted_alleles.index(ref_allele))
                ref_allele_depth = sorted_alleles_depth[ref_index]
                AD_string = str(ref_allele_depth) #Sample string: ref allele depth, alt allele depth(s)
                ref_AF = ((ref_allele_depth)/(int(consdps[index][f_size])))*100
                freq_string = "{:.2f}".format(ref_AF)

                filt[index][f_size] = "a10"

                for count in range(4):

                    if sorted_alleles_depth[count] >= 10:
                        filt[index][f_size] = "PASS"

                    #Check that the allele is not the ref allele and that it has a depth greater than 0, and append to the appropriate string
                    if count != ref_index and sorted_alleles_depth[count] != 0:

                        #AD_string consists of allele depths as: ref_allele_depth, alt_alleles_depth(s)
                        #alt_string consists of all alternate alleles at each base position, ordered from most-frequenct to least-frequent
                        #freq_string consists of the frequency that each alternate allele occurs at a given base position
 
                        AD_string = AD_string+","+str(sorted_alleles_depth[count])
                    
                        if alt_string != "":
                            alt_string = alt_string+","+sorted_alleles[count] 
                        else:
                            alt_string = sorted_alleles[count]
                   
                        AF_value = ((sorted_alleles_depth[count])/(int(consdps[index][f_size])))*100

                        freq_string = freq_string+","+"{:.2f}".format(AF_value) 
 
                if alt_string != "":
                    alt_string_info[index][f_size] = alt_string
                else:
                    alt_string_info[index][f_size] = "N/A"
                allele_data[index][f_size] = "AD="+AD_string+":AF="+freq_string

                sample_data[index][f_size] = "RDP="+rawdps[index][f_size]+":CDP="+consdps[index][f_size]+":MNF="+mean_fams[index][f_size]+":"+allele_data[index][f_size]
                if (f_size != None) and (intvls == None):
                    min_fam[index][f_size] = "MIF="+str(min(i for i in f_sizes if int(i) > 0))
                elif (f_size != None) and (intvls != None):
                    intvl_str = "INTVL="+intvls[index][f_size]
                    min_fam_str = "MIF="+str(min(i for i in f_sizes if int(i) > 0))
                    min_fam[index][f_size] = min_fam_str+":"+intvl_str




def get_vcf_output(cons_file, contig, region_start, region_end, output_path, config):

    #f_sizes = [int(n) for n in config['SETTINGS']['min_family_sizes'].split(',')] if config else [1, 2, 5]
    f_sizes = [1,2,5]
    max_size = max(f_sizes)+1

    fams, contigs, ref_base, base_positions, base_A, base_C, base_G, base_T, rawdps, consdps, mean_fams = ([] for i in range(11))
    temp_fams, temp_cont, temp_ref, temp_pos, temp_baseA, temp_baseC, temp_baseG, temp_baseT, temp_rawdps, temp_consdps, temp_mf = ([None for x in range(max_size)] for i in range(11))

    fmt_string = "RDP:CDP:MNF:AD:AF"
    #ref_threshold = float(config['REPORT']['percent_ref_threshold']) if config else 95.0
    ref_threshold = 100.0

    #Read data from each line of the cons file into the appropriate 2-D list, and ignore/skip cons file if it does not exist
    try:
        with open(cons_file, "r") as file:
            first_line = file.readline().strip() #Get the first line in the cons file 
            headers = first_line.split('\t')
            reader = csv.DictReader(file, delimiter='\t', fieldnames=headers)
            next(reader)

            for line in reader:
                #Create vcf record for base position only if ref_freq <= ref_threshold when reads are uncollapsed
                fam_size = int(line['FAM'])
                ref_freq = float(line['REF_FREQ'])
                if ref_freq <= ref_threshold:
                    if fam_size in f_sizes:

                        mfam = float(line['MEAN_FAM'])
                        mfam = "{:.2f}".format(mfam)

                        #temp_fams[fam_size]=line['FAM']; temp_cont[fam_size]=line['CHROM']; temp_ref[fam_size]=line['REF']; temp_pos[fam_size]=line['POS']; temp_baseA[fam_size]=line["'A's"]; temp_baseC[fam_size]=line["'C's"]; temp_baseG[fam_size]=line["'G's"]; temp_baseT[fam_size]=line["'T's"]; temp_rawdps[fam_size]=line['RAWDP']; temp_consdps[fam_size]=line['CONSDP']; temp_mf[fam_size]=str(mfam)

                        temp_fams[fam_size]=line['FAM']; temp_cont[fam_size]=line['CHROM']; temp_ref[fam_size]=line['REF']; temp_pos[fam_size]=line['POS']; temp_baseA[fam_size]=line['A']; temp_baseC[fam_size]=line['C']; temp_baseG[fam_size]=line['G']; temp_baseT[fam_size]=line['T']; temp_rawdps[fam_size]=line['RAWDP']; temp_consdps[fam_size]=line['CONSDP']; temp_mf[fam_size]=str(mfam)

                        if fam_size == max_size-1:
                            contigs.append(temp_cont); base_positions.append(temp_pos); ref_base.append(temp_ref); base_A.append(temp_baseA); base_C.append(temp_baseC); base_G.append(temp_baseG); base_T.append(temp_baseT); rawdps.append(temp_rawdps); consdps.append(temp_consdps); fams.append(temp_fams); mean_fams.append(temp_mf)
                            temp_fams, temp_cont, temp_ref, temp_pos, temp_baseA, temp_baseC, temp_baseG, temp_baseT, temp_rawdps, temp_consdps, temp_mf = ([None for x in range(max_size)] for i in range(11))

    except FileNotFoundError:
        pass


    filt, alt_string_info, allele_data, sample_data, min_fam = ([[0 for x in range(max_size)] for y in range(len(contigs))] for i in range(5))

    get_data(f_sizes, ref_base, base_A, base_C, base_G, base_T, filt, consdps, alt_string_info, allele_data, sample_data, rawdps, mean_fams, min_fam)

    write_vcf(config, contigs, f_sizes, base_positions, region_start, region_end, ref_base, alt_string_info, filt, min_fam, fmt_string, sample_data, output_path, cons_is_merged=False)



#----------------------------MOD------------------------------[


def get_vcf_output2(cons_file, region_start, region_end, output_path, config):

    cons_is_merged = check_consfile(cons_file)
    print("---vcf cons_is_merged: "+str(cons_is_merged))

    #f_sizes = [int(n) for n in config['SETTINGS']['min_family_sizes'].split(',')] if config else [1, 2, 5]
    f_sizes=[1,2,5]
    max_size = max(f_sizes)+1

    fams, contigs, ref_base, base_positions, base_A, base_C, base_G, base_T, rawdps, consdps, mean_fams = ([] for i in range(11))
    temp_fams, temp_cont, temp_ref, temp_pos, temp_baseA, temp_baseC, temp_baseG, temp_baseT, temp_rawdps, temp_consdps, temp_mf = ([None]*max_size for i in range(11))

    if cons_is_merged:
        intvls = []	
        temp_intvls = [None]*max_size
    else:
        intvls = None

    fmt_string = "RDP:CDP:MNF:AD:AF"
    #ref_threshold = float(config['REPORT']['percent_ref_threshold']) if config else 95.0
    ref_threshold = 100.0

    #Read data from each line of the cons file into the appropriate 2-D list, and ignore/skip cons file if it does not exist
    try:
        with open(cons_file, "r") as file:
            first_line = file.readline().strip() #Get the first line in the cons file
            headers = first_line.split('\t')
            reader = csv.DictReader(file, delimiter='\t', fieldnames=headers)
            next(reader)

            for line in reader:
                #Create vcf record for base position only if ref_freq <= ref_threshold when reads are uncollapsed
                fam_size = int(line['FAM'])
                ref_freq = float(line['REF_FREQ'])

                if ref_freq <= ref_threshold:
                    if fam_size in f_sizes:

                        mfam = float(line['MEAN_FAM'])
                        mfam = "{:.2f}".format(mfam)

                        temp_fams[fam_size]=line['FAM']; temp_cont[fam_size]=line['CHROM']; temp_ref[fam_size]=line['REF']; temp_pos[fam_size]=line['POS']; temp_baseA[fam_size]=line['A']; temp_baseC[fam_size]=line['C']; temp_baseG[fam_size]=line['G']; temp_baseT[fam_size]=line['T']; temp_rawdps[fam_size]=line['RAWDP']; temp_consdps[fam_size]=line['CONSDP']; temp_mf[fam_size]=str(mfam)


                        #temp_fams[fam_size]=line['FAM']; temp_cont[fam_size]=line['CHROM']; temp_ref[fam_size]=line['REF']; temp_pos[fam_size]=line['POS']; temp_baseA[fam_size]=line["'A's"]; temp_baseC[fam_size]=line["'C's"]; temp_baseG[fam_size]=line["'G's"]; temp_baseT[fam_size]=line["'T's"]; temp_rawdps[fam_size]=line['RAWDP']; temp_consdps[fam_size]=line['CONSDP']; temp_mf[fam_size]=str(mfam)


                        if cons_is_merged:
                            temp_intvls[fam_size]=line['INTVL']

                        if fam_size == max(f_sizes):
                            if cons_is_merged == False:
                                contigs.append(temp_cont); base_positions.append(temp_pos); ref_base.append(temp_ref); base_A.append(temp_baseA); base_C.append(temp_baseC); base_G.append(temp_baseG); base_T.append(temp_baseT); rawdps.append(temp_rawdps); consdps.append(temp_consdps); fams.append(temp_fams); mean_fams.append(temp_mf)

                                temp_fams, temp_cont, temp_ref, temp_pos, temp_baseA, temp_baseC, temp_baseG, temp_baseT, temp_rawdps, temp_consdps, temp_mf = ([None for x in range(max_size)] for i in range(11))


                            elif cons_is_merged == True:
                                contigs.append(temp_cont); base_positions.append(temp_pos); ref_base.append(temp_ref); base_A.append(temp_baseA); base_C.append(temp_baseC); base_G.append(temp_baseG); base_T.append(temp_baseT); rawdps.append(temp_rawdps); consdps.append(temp_consdps); fams.append(temp_fams); mean_fams.append(temp_mf); intvls.append(temp_intvls)
                                temp_fams, temp_cont, temp_ref, temp_pos, temp_baseA, temp_baseC, temp_baseG, temp_baseT, temp_rawdps, temp_consdps, temp_mf, temp_intvls = ([None for x in range(max_size)] for i in range(12))


    except FileNotFoundError:
        pass


    filt, alt_string_info, allele_data, sample_data, min_fam = ([[0 for x in range(max_size)] for y in range(len(contigs))] for i in range(5))

    print("---Calling get_data---")
    get_data2(f_sizes, ref_base, base_A, base_C, base_G, base_T, filt, consdps, alt_string_info, allele_data, sample_data, rawdps, mean_fams, min_fam, intvls)

    print("---Calling write_vcf---")
    write_vcf(config, contigs, f_sizes, base_positions, region_start, region_end, ref_base, alt_string_info, filt, min_fam, fmt_string, sample_data, output_path, cons_is_merged)


#---------------------------MOD--------------------------------]




def get_vcf_output_original(cons_file, contig, region_start, region_end, output_path, config):

    f_sizes = [str(n) for n in config['SETTINGS']['min_family_sizes'].split(',')] if config else ['1', '2', '5']

    fams, contigs, ref_base, base_positions, base_A, base_C, base_G, base_T, rawdps, consdps, mean_fams = ([] for i in range(11))
    fmt_string = "RDP:CDP:MNF:AD:AF"

    ref_threshold = float(config['REPORT']['percent_ref_threshold']) if config else 95.0

    #Read data from each line of the cons file into the appropriate 2-D list, and ignore/skip cons file if it does not exist 
    try:
        with open(cons_file, "r") as reader:
            next(reader)
            index = 0
            for line in reader:
                fam_size = int(line.split('\t')[12])
                if fam_size == 0:
                    #Create vcf record for base position only if ref_freq <= ref_threshold when reads are uncollapsed
                    ref_freq = float(line.split('\t')[13])
                    continue
                if ref_freq <= ref_threshold: 
                    if fam_size == int(f_sizes[0]):
                        index += 1
                        contigs.append([]); base_positions.append([]); ref_base.append([]); base_A.append([]); base_C.append([]); base_G.append([]); base_T.append([]); rawdps.append([]); consdps.append([]); fams.append([]); mean_fams.append([])
                        
                    contigs[index-1].append(line.split('\t')[0])
                    base_positions[index-1].append(int(line.split('\t')[1]))
                    ref_base[index-1].append(line.split('\t')[2])
                    base_A[index-1].append(line.split('\t')[3])
                    base_C[index-1].append(line.split('\t')[4])
                    base_G[index-1].append(line.split('\t')[5])
                    base_T[index-1].append(line.split('\t')[6])
                    rawdps[index-1].append(line.split('\t')[10])
                    consdps[index-1].append(line.split('\t')[11]) 
                    fams[index-1].append(line.split('\t')[12])
                    mfam = float(line.split('\t')[14])
                    mfam = "{:.2f}".format(mfam)
                    mean_fams[index-1].append(str(mfam))


    except FileNotFoundError:
        pass


    filt, alt_string_info, allele_data, sample_data, min_fam = ([[0 for x in range(len(f_sizes))] for y in range(len(contigs))] for i in range(5))    

    get_data(f_sizes, ref_base, base_A, base_C, base_G, base_T, filt, consdps, alt_string_info, allele_data, sample_data, rawdps, mean_fams, min_fam)

    write_vcf(config, contigs, f_sizes, base_positions, region_start, region_end, ref_base, alt_string_info, filt, min_fam, fmt_string, sample_data, output_path)




if __name__=="__main__":
    ## Argument + config parsing and error handling
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cons_file', help='Path to your cons_file.')
    parser.add_argument('-f', '--f_sizes', help='Comma-separated list of f_sizes to make VCF files for.')
    parser.add_argument('-r', '--region', help='Region to analyze (string of the form chrX:posA-posB).')
    parser.add_argument('-o', '--output_path', help='Path to write output files to.')
    parser.add_argument('-c', '--config', help='Path to your config file.')

    args = parser.parse_args()

    if args.config:
        config = configparser.ConfigParser()
        config.read(args.config)
    else:
        config = None
    
    cons_file = args.cons_file

    f_sizes = args.f_sizes.split(',')

    region = args.region
    if any(x not in region for x in ["chr", ":", "-"]):
        raise ValueError('Incorrect region string (should look like chr1:1200000-1250000).')
    sys.exit(1)

    contig = region.split(":")[0]
    region_start = int(region.split(":")[1].split("-")[0])
    region_end = int(region.split(":")[1].split("-")[1])

    output_path = handle_arg(args.output_path, config['PATHS']['output_path'] if config else None, 
                    'No output path provided in args or config.')

    bam_file = handle_arg(args.bam_file, config['PATHS']['bam_file'] if config else None, 
                    'No BAM file provided in args or config.')

    ## Output
    #Original
    #generate_vcf_output(cons_file, f_sizes, contig, region_start, region_end, output_path, config)

    get_vcf_output(cons_file, contig, region_start, region_end, output_path, config)
