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
from src.generate_vcf import get_vcf_output2


"""
/src/generate_merge.py 
=========================================
Purpose:

Script contains sub-functions to merge umi and cons data files, and functions
used in Debarcer's 'run' sub-process.

Author: Isha Warikoo
Copyright (c) 2018 GSI, Ontario Institute for Cancer Research

"""



def submit_jobs(bamfile, bedfile, output_dir, config, index):
	
	for i in range(index,len(lines)):
		chromosome = lines[i][0]
		pos1 = lines[i][1]
		pos2 = lines[i][2]

		#Create umi scripts
		f = open(output_dir+"umifiles/umigrp_"+chromosome+"_"+pos1+".sh","w")
		os.system("chmod +x "+output_dir+"umifiles/umigrp_"+chromosome+"_"+pos1+".sh")
		f.write("module load /.mounts/labs/PDE/Modules/modulefiles/python-gsi/3.6.4")
		if config:
			f.write("\npython3.6 "+"./debarcer.py group -o "+output_dir+"umifiles/ -r "+chromosome+"\:"+pos1+"-"+pos2+" -b "+bamfile+" -c "+config)
		else:
			f.write("\npython3.6 "+"./debarcer.py group -o "+output_dir+"umifiles/ -r "+chromosome+"\:"+pos1+"-"+pos2+" -b "+bamfile)
		#os.system("qsub -cwd -b y -N UMI_"+str(chromosome)+"_"+str(pos1)+" -e logs -o logs -l h_vmem=10g "+output_dir+"umifiles/umigrp_"+chromosome+"_"+pos1+".sh")


		#Create cons scripts
		f = open(output_dir+"consfiles/cons_"+chromosome+"_"+pos1+".sh","w")
		os.system("chmod +x "+output_dir+"consfiles/cons_"+chromosome+"_"+pos1+".sh")
		f.write("module load /.mounts/labs/PDE/Modules/modulefiles/python-gsi/3.6.4")
		if config:
			f.write("\npython3.6 "+"./debarcer.py collapse -o "+output_dir+"consfiles/ -r "+chromosome+"\:"+pos1+"-"+pos2+" -b "+bamfile+" -u "+output_dir+"umifiles/"+chromosome+"\:"+pos1+"-"+pos2+".umis -c "+config)
		else:
			f.write("\npython3.6 "+"./debarcer.py collapse -o "+output_dir+"consfiles/ -r "+chromosome+"\:"+pos1+"-"+pos2+" -b "+bamfile+" -u "+output_dir+"umifiles/"+chromosome+"\:"+pos1+"-"+pos2+".umis")
		#os.system("qsub -cwd -b y -N CONS_"+str(chromosome)+"_"+str(pos1)+" -e logs -o logs -l h_vmem=10g -hold_jid UMI_* "+output_dir+"consfiles/cons_"+chromosome+"_"+pos1+".sh")

		#Create vcf scripts
		f = open(output_dir+"vcffiles/call_"+chromosome+"_"+pos1+".sh","w")
		os.system("chmod +x "+output_dir+"vcffiles/call_"+chromosome+"_"+pos1+".sh")
		f.write("module load /.mounts/labs/PDE/Modules/modulefiles/python-gsi/3.6.4")
		if config:
			f.write("\npython3.6 "+"./debarcer.py call -o "+output_dir+"vcffiles/ -r "+chromosome+"\:"+pos1+"-"+pos2+" -cf "+output_dir+"consfiles/"+chromosome+"\:"+pos1+"-"+pos2+".cons -f 1,2,5 -c "+config)
		else:
			f.write("\npython3.6 "+"./debarcer.py call -o "+output_dir+"vcffiles/ -r "+chromosome+"\:"+pos1+"-"+pos2+" -cf "+output_dir+"consfiles/"+chromosome+"\:"+pos1+"-"+pos2+".cons")
		#os.system("qsub -cwd -b y -N CALL_"+str(chromosome)+"_"+str(pos1)+" -e logs -o logs -l h_vmem=10g -hold_jid CONS_* "+output_dir+"vcffiles/call_"+chromosome+"_"+pos1+".sh")


def check_merge_flag(var):
	if var == 'umi': 
		print("Pass umi!")
	elif var == 'cons':
		print("Pass cons!")
	else:
		print("ERR: Incorrect flag, "+var+" specified. Input is expected to be 'umi', 'cons' or 'vcf'")


def merge_umi_datafiles(output_path):
	path=output_path+"umifiles/"
	merged_file=output_path+"umifiles/merged_file.csv"
	os.chdir(path)
	#results = pd.DataFrame([])
	results=[]

	file = open(merged_file, "w")

	headers = ['CHR', 'START', 'END', 'PTU', 'CTU', 'CHILD_NUMS', 'FREQ_PARENTS']
	csv.register_dialect('myDialect', delimiter='\t', quoting=csv.QUOTE_NONE)
	
	writer = csv.DictWriter(file, dialect='myDialect', fieldnames=headers)
	writer.writeheader()

	for counter, file in enumerate(glob.glob("*.csv")):
		f = open(file, "r")
		reader = csv.DictReader(f, delimiter='\t', fieldnames=headers)
		next(reader)
		for row in reader:
			contig=row['CHR']; start=row['START']; end=row['END']; total_pumis=row['PTU']; total_cumis=row['CTU']; child_nums=row['CHILD_NUMS']; freq_of_parent_umis=row['FREQ_PARENTS']
			csvrow = {'CHR' : contig, 'START' : start, 'END' : end, 'PTU' : str(total_pumis), 'CTU' : str(total_cumis), 'CHILD_NUMS': child_nums, 'FREQ_PARENTS' : freq_of_parent_umis}
			writer.writerow(csvrow)	


def merge_umi_datafiles2(output_path):
	path=output_path+"umifiles/"
	merged_file=output_path+"umifiles/merged2_file.csv"
	sorted_merge=output_path+"umifiles/merged2_sorted.csv"

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


def concat_cons(output_path, config):

	path=output_path+"consfiles/"
	#merged_file=output_path+"consfiles/SORTED_CONCAT_data.cons"
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

	merged_file=output_path+"consfiles/Merged_"+str(first_region)+"_"+str(last_region)+".cons"

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

	"""
	with open(merged_file, 'r') as f:
		lines = f.read().splitlines()
		first_line = lines[1]
		last_line = lines[-1]

		first_region = first_line.split('\t')[0]
		last_region = last_line.split('\t')[0]

		region=first_region+"_"+last_region		

	with open(sorted_names_file, 'r') as f:
		lines = f.read().splitlines()
		first_line = lines[0]
		last_line = lines[-1]

		first_region = first_line.split('/')[-1]
		first_region = first_line.split('.')[0]

		last_region = last_line.split('/')[-1]
		last_region = last_line.split('.')[1]

		region=first_region+"_"+last_region
	"""

	#os.system("rm "+sorted_names_file)


	print("Running variant call")
	#debarcer_path = config.rsplit('/', 2)[0]
	#debarcer_path = debarcer_path+"/debarcer/debarcer.py"	
	#print(debarcer_path)
	
	get_vcf_output2(cons_file=merged_file, region_start=first_region, region_end=last_region, output_path=path, config=config)


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
		f.write(region+'\n')
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



"""
def check_consfile(cons_file):
	f = open(cons_file, "r")
	#reader = csv.reader(f, delimiter='\t')
	line = f.readline()
	if 'INTVL' in line:
		return True
	else:
		return False
"""



def temp_cons(output_dir, region, row, merged_file):
	file_path=output_path+"/consfiles/temp_"+region+".cons"

	headers = ['INTVL', 'CHROM', 'POS', 'REF', 'A', 'C', 'G', 'T', 'I', 'D', 'N', 'RAWDP', 'CONSDP', 'FAM', 'REF_FREQ', 'MEAN_FAM']
	#f = open(file_path, "w")
	#writer = csv.DictWriter(f, dialect='myDialect' , field=headers)
	
	f = open(file_path, "w")	
	writer = csv.writer(f)
	writer.writerow(headers)

	r_file = open(merged_file, "r")
	#reader = csv.DictReader(r_file, delimiter='\t', fieldnames=headers)
	reader = csv.reader(r_file)
	counter = int(row[0])
	while counter<= int(row[1]):
		writer.writerow = reader[counter]
		counter+=1


	return file_path


