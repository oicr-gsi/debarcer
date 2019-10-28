import operator
import os
import json
import subprocess
from src.utilities import CheckRegionFormat, CheckJobs
import uuid


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
        newfile.write('\t'.join(i[-1].split()) + '\n')
    newfile.close()


def name_job(prefix):
    '''
    (str) -> str
    
    :param prefix: Job description (eg UmiCollapse)
    
    Return a unique job name by appending a random string to prefix
    '''
   
    # get a random string
    UID = str(uuid.uuid4())
    jobname = prefix + '_' + UID
    return jobname    

    
def submit_jobs(bamfile, outdir, reference, famsize, bedfile, count_threshold,
                consensus_threshold, dist_threshold, post_threshold, ref_threshold,
                alt_threshold, filter_threshold, maxdepth, truncate, ignoreorphans, ignore, merge,
                plot, report, call, mincov, minratio, minumis, minchildren, extension,
                sample, mydebarcer, mypython, mem, queue):
    '''
    (str, str, str, str, str, int, float, int, int, float, float, int, int, bool,
    bool, bool, bool, bool, bool, bool,
    
    
    str, str, str, str, int, str) -> None
    
    :param bamfile: Path to the bam file
    :param outdir: Directory where .umis, and datafiles are written
    :param reference: Path to the reference genome
    :param famsize: Comma-separated list of minimum umi family size to collapase on
    :param bedfile: Bed file with region coordinates. Chromsome must start with chr and positions are 1-based inclusive  
    :param count_threshold: Base count threshold in pileup column
    :param consensus_threshold: Majority rule consensus threshold in pileup column
    :param dist_threshold: Hamming distance threshold for connecting parent-children umis
    :param post_threshold: Distance threshold in bp for defining families within groups
    :param ref_threshold: Maximum reference frequency (in %) to consider alternative variants
                          (ie. position with ref freq <= ref_threshold is considered variable)
    :param alt_threshold: Minimum allele frequency (in %) to consider an alternative allele at a variable position 
                          (ie. allele freq >= alt_threshold and ref freq <= ref_threshold --> record alternative allele)
    :param filter_threshold: Minimum number of reads to pass alternative variants 
                             (ie. filter = PASS if variant depth >= alt_threshold)
    :param maxdepth: Maximum read depth. Default is 1000000
    :param truncate: Only consider pileup columns in given region. Default is True\
    :param ignoreorphans: Ignore orphans (paired reads that are not in a proper pair). Default is True
    :param ignore: Keep the most abundant family and ignore families at other positions within each group if True. Default is False
    :param merge: Merge datafiles, consensus files and umi files if True
    :param plot: Generate figure plots if True
    :param report: Generate analysis report if True
    :param call: Convert consensus files to VCF if True
    :param mincov: Minimum read depth to label regions
    :param minratio: Minimum ratio to label regions    
    :param minumis: Minimum number of umis to label regions
    :param minchildren: Minimum number of umi children to label regions
    :param extension: Figure file extension
    :param sample: Sample name to appear in report. If empty str, outdir basename is used
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
    GroupCmd = '{0} {1} group -o {2} -r \"{3}\" -b {4} -d {5} -p {6} -i {7} -t {8}'
    # set up collapse cmd
    CollapseCmd = 'sleep 60; {0} {1} collapse -o {2} -b {3} -rf {4} -r \"{5}\" -u {6} -f \"{7}\" -ct {8} -pt {9} -p {10} -m {11} -t {12} -i {13}'
        
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
        newfile.write(GroupCmd.format(mypython, mydebarcer, outdir, region, bamfile, str(dist_threshold), str(post_threshold), ignore, str(truncate)) + '\n')
        newfile.close()
        # get a umique job name
        jobname1 = name_job('UmiGroup' + '_' + region.replace(':', '-'))
        # run group umi for region
        subprocess.call(QsubCmd1.format(jobname1, LogDir, queue, str(mem), GroupScript), shell=True)      
        # record jobname
        GroupJobNames.append(jobname1)
        
        # dump collapse cmd into a shell script  
        umifile = os.path.join(UmiDir, '{0}.json'.format(region))
        CollapseScript = os.path.join(QsubDir, 'UmiCollapse_{0}.sh'.format(region.replace(':', '_').replace('-', '_')))
        newfile = open(CollapseScript, 'w')
        newfile.write(CollapseCmd.format(mypython, mydebarcer, outdir, bamfile, reference, region, umifile,
                                         str(famsize), str(count_threshold), str(consensus_threshold),
                                         str(post_threshold), str(maxdepth), str(truncate), str(ignoreorphans)) +'\n') 
        newfile.close()
        # get a umique job name
        jobname2 = name_job('UmiCollapse' + '_' + region.replace(':', '-'))
        # run collapse umi for region
        subprocess.call(QsubCmd2.format(jobname2, jobname1, LogDir, queue, str(mem), CollapseScript), shell=True)      
        # record jobname2
        ConsJobNames.append(jobname2)
        
    
    # make a list of merge job names
    MergeJobNames = []        
    if merge  == True:
        Regions = list(map(lambda x: x.replace(':', '-'), Regions))
        
        # submit jobs to merge 
        MergeCmd = 'sleep 600; {0} {1} merge -d {2} -dt {3}'

        # check if Goup jobs are still running
        running_group = CheckJobs(GroupJobNames)
        if running_group == False:
            # merge datafiles
            MergeScript1 = os.path.join(QsubDir, 'MergeDataFiles.sh')
            newfile = open(MergeScript1, 'w')
            newfile.write(MergeCmd.format(mypython, mydebarcer, DataDir, 'datafiles') + '\n') 
            newfile.close()
            jobname3 =  name_job('MergeDataFiles')
            MergeJobNames.append(jobname3)
            # run merge datafiles
            #subprocess.call(QsubCmd2.format(jobname3, GroupJobNames[-1], LogDir, queue, '20', MergeScript1), shell=True)    
            subprocess.call(QsubCmd1.format(jobname3, LogDir, queue, '20', MergeScript1), shell=True)    
                
            # merge umi files     
            MergeScript2 = os.path.join(QsubDir, 'MergeUmiFiles.sh')
            newfile = open(MergeScript2, 'w')
            newfile.write(MergeCmd.format(mypython, mydebarcer, UmiDir, 'umifiles') + '\n')
            newfile.close()
            jobname4 = name_job('MergeUmiFiles')
            MergeJobNames.append(jobname4)
            # run merge umi files
            #subprocess.call(QsubCmd2.format(jobname5, GroupJobNames[-1], LogDir, queue, '20', MergeScript3), shell=True)
            subprocess.call(QsubCmd1.format(jobname4, LogDir, queue, '20', MergeScript2), shell=True)
         
        # check if collapse jobs are still running
        running_collapse = CheckJobs(ConsJobNames)
        if running_collapse == False:
            # merge consensus files
            MergeScript3 = os.path.join(QsubDir, 'MergeConsensusFiles.sh')
            newfile = open(MergeScript3, 'w')
            newfile.write(MergeCmd.format(mypython, mydebarcer, ConsDir, 'consensusfiles') + '\n') 
            newfile.close()
            jobname5 = name_job('MergeConsensusFiles')
            MergeJobNames.append(jobname5)
            # run merge consensus files
            #subprocess.call(QsubCmd2.format(jobname4, ConsJobNames[-1], LogDir, queue, '20', MergeScript2), shell=True)    
            subprocess.call(QsubCmd1.format(jobname5, LogDir, queue, '20', MergeScript3), shell=True)    
        
    # make a list of call jobs
    CallJobs = [] 
    if call == True:
        # make a list of jobs, wait until all jobs are done before converting consensus files to VCF 
        Z = GroupJobNames + MergeJobNames
        running_groupmerge = CheckJobs(Z)  
        if running_groupmerge == False:
            # generate VCF from all consensus files 
            # set up vcf command
            VarCallCmd = 'sleep 600; {0} {1} call -o {2} -rf {3} -f \"{4}\" -rt {5} -at {6} -ft {7}'
            CallScript = os.path.join(QsubDir, 'VarCall.sh')
            newfile = open(CallScript, 'w')
            newfile.write(VarCallCmd.format(mypython, mydebarcer, outdir, reference, famsize, ref_threshold, alt_threshold, filter_threshold))
            newfile.close()    
            jobname6 = name_job('Call')
            CallJobs.append(jobname6)
            subprocess.call(QsubCmd1.format(jobname6, LogDir, queue, '20', CallScript), shell=True)    
    
    if plot == True:
        # make a list of jobs. wait until all jobs are done before plotting and reporting
        L = ConsJobNames + GroupJobNames + MergeJobNames + CallJobs
        running_jobs = CheckJobs(L)
        if running_jobs == False:
            # generate plots and report if report is True
            PlotCmd = 'sleep 600; {0} {1} plot -d {2} -e {3} -s {4} -r {5} -mv {0} -mr {7} -mu {8} -mc {9}'
            PlotScript = os.path.join(QsubDir, 'PlotFigures.sh')
            newfile = open(PlotScript, 'w')
            newfile.write(PlotCmd.format(mypython, mydebarcer, outdir, extension, sample, report, mincov, minratio, minumis, minchildren))
            newfile.close()
            jobname7 = name_job('Plot')
            subprocess.call(QsubCmd1.format(jobname7, LogDir, queue, '20', PlotScript), shell=True)    
    