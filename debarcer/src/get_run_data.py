import sys
import pysam
import configparser
import argparse
import operator
import functools
import os
import csv
import matplotlib.pyplot as plt
import pandas as pd
import fnmatch
import itertools
import glob
import numpy as np
import time
import json
from src.generate_vcf import get_vcf_output
import subprocess
from src.utilities import CheckRegionFormat


def ExtractRegions(bedfile):
    '''
    (str) -> list
    
    :param bedfile: Path to bedfile. Positions are 1-based included
    
    Returns a list of regions extracted from the bed file
    Precondition: The 1st 3 columns must be chrN start end.
                  Other columns are ignored
                  Chromsome name must starts with 'chr'.
    '''
    
    # create a list to store regions
    L = []
    
    # extract regions from bed file
    infile = open(bedfile)
    for line in infile:
        if line.startswith('chr'):
            line = line.rstrip().split()
            chromo, start, end = line[0], line[1], line[2]
            region = chromo + ':' + start + '-' + end
            L.append(region)
    infile.close()
    
    return L



def check_umi_status(output_dir):
    time.sleep(3)
    os.system("qstat -r | grep 'UMI*' > "+output_dir+"temp_umi_jobs.txt")
    job_flag = os.stat(output_dir+"temp_umi_jobs.txt").st_size == 0
    
    return job_flag


def check_job_status(output_dir, flag, file):
    time.sleep(3)
    
    if flag == 'umi':
        os.system("qstat -r | grep 'UMI*' > "+output_dir+file)
    elif flag == 'cons':
        os.system("qstat -r | grep 'CONS*' > "+output_dir+file)
    
    job_flag = os.stat(output_dir+file).st_size == 0
    
    return job_flag


def check_cons_status(output_dir):
    time.sleep(3)
    os.system("qstat -r | grep 'CONS*' > "+output_dir+"temp_cons_jobs.txt")
    job_flag = os.stat(output_dir+"temp_cons_jobs.txt").st_size == 0
    
    return job_flag


def MergeUmiFiles(UmiDir):
    '''
    (str) -> None

    :param UmiDir: Directory Umifiles containing .json umi files
    
    Merge all .json umi files in Umifiles into a unique Merged_UmiFile.json file in Umifiles
    '''
    
    # make a list of umi files
    UmiFiles = [i for i in os.listdir(UmiDir) if i[-5:] == '.json']
    
    # create a dict to hold all umi data
    D = {}
    
    for i in UmiFiles:
        # extract coordinates from file name
        coord = i[:-5]
        # get umi data as dict
        umifile = os.path.join(UmiDir, i)    
        infile = open(umifile)
        data = json.load(infile)
        infile.close()
        # add data in dict
        D[coord] = data
    
    # write merged umi file
    MergedFile = os.path.join(UmiDir, 'Merged_UmiFile.json')
    newfile = open(MergedFile, 'w')
    json.dump(D, newfile, sort_keys = True, indent=4)
    newfile.close()
       

def MergeConsensusFiles(ConsDir):
    '''
    (str) -> None

    :param ConsDir: Directory Consfiles containing .cons consensus files
    
    Merge all .cons consensus files in Consfiles into a unique Merged_ConsensusFile.cons file in Consfiles
    '''
    
    # make a list of consensus files
    ConsFiles = [i for i in os.listdir(ConsDir) if i[-5:] == '.cons' and i.startswith('chr')] 
    # sort files
    L = []
    for i in ConsFiles:
        # extract chromo and start from filename
        chromo, start = i[:i.index(':')], int(i[i.index(':')+1:i.index('-')])
        L.append((chromo, start, i))
    # sort files on chromo and start
    L.sort(key=operator.itemgetter(0, 1))
    
    # make a sorted list of full paths
    S = [os.path.join(ConsDir, i[-1]) for i in L]
    
    # get Header
    infile = open(S[0])
    Header = infile.readline().rstrip().split()
    infile.close()
    
    # write merged consensus file
    MergedFile = os.path.join(ConsDir, 'Merged_ConsensusFile.cons')
    newfile = open(MergedFile, 'w')
    newfile.write('\t'.join(Header) + '\n')
    for i in S:
        infile = open(i)
        # skip header and grab all data
        infile.readline()
        data = infile.read().rstrip()
        infile.close()
        newfile.write(data + '\n')
    newfile.close()
    
    

def MergeDataFiles(DataDir):
    '''
    (str) -> None

    :param DataDir: Directory Datafiles containing .csv data files
    
    Merge all csv datafiles in Datafiles into a unique Merged_DataFile.csv file in Datafiles
    '''
    
    # make a list of datafiles
    DataFiles = [os.path.join(DataDir, i) for i in os.listdir(DataDir) if i.startswith('datafile_') and i[-4:] == '.csv']
    
    Header = ['CHR', 'START', 'END', 'PTU', 'CTU', 'CHILD_NUMS', 'FREQ_PARENTS']
        
    # read each datafile, store in a list
    L = []
    for filename in DataFiles:
        infile = open(filename)
        # skip header and grap data
        infile.readline()
        data = infile.read().rstrip()
        infile.close()
        # get chromosome and start and end positions
        chromo, start = data.split()[0], int(data.split()[1])        
        L.append((chromo, start, data))
    # sort data on chromo and start
    L.sort(key=operator.itemgetter(0, 1))
    
    # write merged datafile
    MergedFile = os.path.join(DataDir, 'Merged_DataFile.csv')
    newfile = open(MergedFile, 'w')
    newfile.write('\t'.join(Header) + '\n')
    for i in L:
        newfile.write('\t'.join(i.split()) + '\n')
    newfile.close()
    

def submit_jobs(bamfile, outdir, reference, famsize, bedfile, countthreshold,
                percentthreshold, distthreshold, postthreshold, refthreshold,
                allthreshold, maxdepth, truncate, ignoreorphans, ignore, merge,
                mydebarcer, mypython, mem, queue):
    '''
    (str, str, str, str, str, int, float, int, int, float, float, int, bool,
    bool, bool, bool, str, str, int, str) -> None
    
    :param bamfile: Path to the bam file
    :param outdir: Directory where .umis, and datafiles are written
    :param reference: Path to the reference genome
    :param famsize: Comma-separated list of minimum umi family size to collapase on
    :param bedfile: Bed file with region coordinates. Chromsome must start with chr and positions are 1-based inclusive  
    :param countthreshold: Base count threshold in pileup column
    :param percentthreshold: Base percent threshold in pileup column
    :param distthreshold: Hamming distance threshold for connecting parent-children umis
    :param postthreshold: Distance threshold in bp for defining families within groups
    :param refthreshold: Reference threshold
    :param allthreshold: Allele threshold
    :param maxdepth: Maximum read depth. Default is 1000000
    :param truncate: Only consider pileup columns in given region. Default is True\
    :param ignoreorphans: Ignore orphans (paired reads that are not in a proper pair). Default is True
    :param ignore: Keep the most abundant family and ignore families at other positions within each group if True. Default is False
    :param merge: Merge datafiles, consensus files and umi files if True
    :param mydebarcer: Path to the debarcer script
    :param mypython: Path to python
    :param mem: Requested memory
    :param queue: Sge queue to submit jobs to 
    
    Submit jobs for Umi Group and Collapse
    '''
    
    # get output directories   
    UmiDir = os.path.join(outdir, 'Umifiles')
    ConsDir = os.path.join(outdir, 'Consfiles')
    QsubDir = os.path.join(outdir, 'Qsubs')
    LogDir = os.path.join(QsubDir, 'Logs')
    DataDir = os.path.join(outdir, 'Datafiles')

    # set up group command
    GroupCmd = '{0} {1} group -o {2} -r \"{3}\" -b {4} -d {5} -p {6} -i {7}'
    # set up collapse cmd
    CollapseCmd = 'sleep 60; {0} {1} collapse -o {2} -b {3} -rf {4} -r \"{5}\" -u {6} -f \"{7}\" -ct {8} -pt {9} -rt {10} -at {11} -p {12} -m {13} -t {14} -i {15}'
    # set qsub command
    QsubCmd1 = 'qsub -b y -cwd -N {0} -o {1} -e {1} -q {2} -l h_vmem={3}g \"bash {4}\"'
    QsubCmd2 = 'qsub -b y -cwd -N {0} -hold_jid {1} -o {2} -e {2} -q {3} -l h_vmem={4}g \"bash {5}\"'
    
    # extract regions from bedfile
    Regions = ExtractRegions(bedfile)

    # create a list of job names
    GroupJobNames, ConsJobNames = [], []

    # loop over regions
    for region in Regions:
        # check region format
        CheckRegionFormat(region)
        
        ### RUN Group and Collapse for a given region ### 
        
        # dump group cmd into a shell script  
        GroupScript = os.path.join(QsubDir, 'UmiGroup_{0}.sh'.format(region.replace(':', '_').replace('-', '_')))
        newfile = open(GroupScript, 'w')
        newfile.write(GroupCmd.format(mypython, mydebarcer, outdir, region, bamfile, str(distthreshold), str(postthreshold), ignore) + '\n')
        newfile.close()
        jobname1 = 'UmiGroup_' + region.replace(':', '_').replace('-', '_')
        # run group umi for region
        subprocess.call(QsubCmd1.format(jobname1, LogDir, queue, str(mem), GroupScript), shell=True)      
        # record jobname
        GroupJobNames.append(jobname1)
        
        # dump collapse cmd into a shell script  
        umifile = os.path.join(UmiDir, '{0}.json'.format(region))
        CollapseScript = os.path.join(QsubDir, 'UmiCollapse_{0}.sh'.format(region.replace(':', '_').replace('-', '_')))
        newfile = open(CollapseScript, 'w')
        newfile.write(CollapseCmd.format(mypython, mydebarcer, outdir, bamfile, reference, region, umifile,
                                         str(famsize), str(countthreshold), str(percentthreshold),
                                         str(refthreshold), str(allthreshold), str(postthreshold),
                                         str(maxdepth), str(truncate), str(ignoreorphans)) +'\n') 
        newfile.close()
        jobname2 = 'UmiCollapse_' + region.replace(':', '_').replace('-', '_')
        # run collapse umi for region
        subprocess.call(QsubCmd2.format(jobname2, jobname1, LogDir, queue, str(mem), CollapseScript), shell=True)      
        # record jobname2
        ConsJobNames.append(jobname2)
        

    if merge  == True:
        
        # submit jobs to merge 
        MergeCmd = '{0} {1} merge -d {2} -dt {3}'

        # collect regions
        Regions = []
        for i in os.listdir(DataDir):
            if i.startswith('datafile_') and i[-4:] == '.csv':
                Regions.append(i[i.index('chr'):-4])
        for i in os.listdir(ConsDir):
            if i[-5:] == '.cons':
                Regions.append(i[:-5])
        for i in os.listdir(UmiDir):
            if i[-5:] == '.umis':
                Regions.append(i[:-5])
        Regions = list(map(lambda x: x.replace(':', '-'), list(set(Regions))))          
        
        # merge datafiles
        MergeScript1 = os.path.join(QsubDir, 'MergeDataFiles.sh')
        newfile = open(MergeScript1, 'w')
        newfile.write(MergeCmd.format(mypython, mydebarcer, DataDir, 'datafiles') + '\n') 
        newfile.close()
        jobname3 = 'MergeDataFiles_' + '_'.join(Regions)
        # run merge datafiles
        subprocess.call(QsubCmd2.format(jobname3, GroupJobNames[-1], LogDir, queue, '10', MergeScript1), shell=True)    
        
        # merge consensus files
        MergeScript2 = os.path.join(QsubDir, 'MergeConsensusFiles.sh')
        newfile = open(MergeScript2, 'w')
        newfile.write(MergeCmd.format(mypython, mydebarcer, ConsDir, 'consensusfiles') + '\n') 
        newfile.close()
        jobname4 = 'MergeConsensusFiles_' + '_'.join(Regions)
        # run merge consensus files
        subprocess.call(QsubCmd2.format(jobname4, ConsJobNames[-1], LogDir, queue, '10', MergeScript2), shell=True)    
        
        # merge umi files     
        MergeScript3 = os.path.join(QsubDir, 'MergeUmiFiles.sh')
        newfile = open(MergeScript3, 'w')
        newfile.write(MergeCmd.format(mypython, mydebarcer, UmiDir, 'umifiles') + '\n')
        newfile.close()
        jobname5 = 'MergeUmiFiles_' + '_'.join(Regions)
        # run merge umi files
        subprocess.call(QsubCmd2.format(jobname5, GroupJobNames[-1], LogDir, queue, '10', MergeScript3), shell=True)
        
        
def merge_umi_datafiles(output_path, id):
    
    path=output_path+"umifiles/"
    merged_file=output_path+"umifiles/merged2_file.csv"
    sorted_merge=output_path+"umifiles/RUNID_"+id+"_Merged_UMI.csv"
    os.chdir(path)
    contig, start, end, total_pumis, total_cumis, child_nums, parent_freq = ([] for i in range(7))
    headers = ['CHR', 'START', 'END', 'PTU', 'CTU', 'CHILD_NUMS', 'FREQ_PARENTS']
    for counter, file in enumerate(glob.glob("datafile_*.csv")):
        f = open(file, "r")
        reader = csv.DictReader(f, delimiter='\t', fieldnames=headers)
        next(reader)
        for row in reader:
            contig.append(row['CHR']); start.append(row['START']); end.append(row['END']); total_pumis.append(row['PTU']); total_cumis.append(row['CTU']); child_nums.append(row['CHILD_NUMS']); parent_freq.append(row['FREQ_PARENTS'])

    data = {'CHR':contig, 'START':start, 'END':end, 'PTU':total_pumis, 'CTU':total_cumis, 'CHILD_NUMS':child_nums, 'FREQ_PARENTS':parent_freq}
    merged_data = pd.DataFrame(data, columns=headers)
    merged_data.set_index('CHR', inplace=True)
    merged_data.to_csv(path_or_buf=merged_file, sep='\t')
    os.system("sort -V "+merged_file+" > "+sorted_merge)
    os.system("rm "+merged_file)

def concat_cons(output_path, config, id):
    path=output_path+"consfiles/"
    vcf_path=output_path+"vcffiles/"
    header_file = output_path+"consfiles/headers.cons"
    sorted_names_file = path+"temp_sorted_filenames.txt"
    os.chdir(path)
    for file in glob.glob("chr*.cons"):
        modify_cons(file, path)
    file_lst = "ls -v "+path+"chr*.cons"
    os.system(file_lst+" > "+sorted_names_file)
    region=""
    first_region=""
    last_region="" 
    with open(sorted_names_file, 'r') as f:
        lines = f.read().splitlines()
        first_line = lines[0]
        last_line = lines[-1]

        first_region = first_line.split('/')[-1]
        first_region = first_region.split('.')[0]
        
        last_region = last_line.split('/')[-1]
        last_region = last_region.split('.')[0]
        
        region=first_region+"_"+last_region
    
    merged_file=output_path+"consfiles/RUNID_"+id+"_Merged_CONS.cons"
    
    file = open(header_file, "w")
    headers = ['INTVL', 'CHROM', 'POS', 'REF', 'A', 'C', 'G', 'T', 'I', 'D', 'N', 'RAWDP', 'CONSDP', 'FAM', 'REF_FREQ', 'MEAN_FAM']
    csv.register_dialect('myDialect', delimiter='\t', quoting=csv.QUOTE_NONE)
    writer = csv.DictWriter(file, dialect='myDialect', fieldnames=headers)
    writer.writeheader()
    file.close()

    #Sort files by region, and concatonate
    mod_file_lst = "$(ls -v "+path+"MOD*.cons)"
    
    os.system("cat "+header_file+" "+mod_file_lst+" > "+merged_file)
    os.system("rm "+path+"MOD*.cons")
    os.system("rm "+header_file)

    print("Running variant call...")
    get_vcf_output(cons_file=merged_file, region_start=first_region, region_end=last_region, output_path=vcf_path, config=config, id=id)


def modify_cons(file_path, output_path):
    output_path = output_path
    file_name = file_path.split('/')[-1]
    region = file_name.split('.')[0]
    with open(file_path) as f:
        num_of_lines = sum(1 for line in f)

    num_of_lines -=1
    interval_file = output_path+"INTVLS_"+file_name
    f = open(interval_file, "w")
    for i in range(num_of_lines):
        write(region+'\n')
    f.close()
    
    #Path to modified cons file, which has no header and has an additional column containing intervals 
    modified_cons = output_path+"MOD_"+file_name
    
    #Create a temporary, un-headered version of the cons file
    unheadered_cons = output_path+"unheadered_temp.cons"
    os.system("sed '1d' "+file_path+" > "+unheadered_cons)
    
    #cmd2 = "awk '{getline l < "+unheadered_cons+"; print $0'\t'l}' "+interval_file+" > "+modified_cons
    cmd = "paste -d '\t' "+interval_file+" "+unheadered_cons+" > "+modified_cons
    os.system(cmd)	
    
    os.system("rm "+interval_file+" "+unheadered_cons)
