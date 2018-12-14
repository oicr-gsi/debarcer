import os
import sys


"""
/src/create_plots.py 
=========================================
Purpose:

Script contains sub-functions for Debarcer's 'plot' sub-process.

Author: Isha Warikoo
Copyright (c) 2018 GSI, Ontario Institute for Cancer Research

"""


#Umi plots

def check_file(file_name, extension):
	file_exists = file_name.exists()
	name = file_name.split('/')[-1]
	ext = name.split('.')[-1]

	if file_exists and ext == extension:
		return True
	else:
		return False


def create_umi_dfs(file_name):

	f = file_name.split('/')[-1]
	name = f.split('.')[0]

	headers=['CHR', 'START', 'END', 'PTU', 'CTU', 'CHILD_NUMS', 'FREQ_PARENTS']
	df_headers=['INTVL', 'PTU', 'CTU', 'CHILD_NUMS', 'FREQ_PARENTS', 'INTVL_SIZE', 'CP']
	region, total_pumis, total_cumis, child_nums, parent_freq, child_to_parent, size_of_intvl = ([] for i in range(7))

	table = []

	f= open(file_name, "r")
	reader = csv.DictReader(f, delimiter='\t', fieldnames=headers)
	next(reader)
	counter = 0
	for row in reader:
		counter+=1

		#Build Sub-dataframe table
		str_cumi_lst = (row['CHILD_NUMS']).split(','); str_pumi_lst = (row['FREQ_PARENTS']).split(','); cumi_lst = []; pumi_lst = [];
		for i in str_cumi_lst:
			cumi_lst.append(int(i))
		for j in str_pumi_lst:
			pumi_lst.append(int(j))

		table.append(cumi_lst); table.append(pumi_lst)


		intvl_name = row['CHR']+":"+row['START']+"-"+row['END']
		region_len = str((int(row['END']))-(int(row['START'])))
		temp_intvl = row['CHR']+":"+row['START']+"+"+region_len
		c_to_p = (float(row['CTU']))/(float(row['PTU']))

		size_of_intvl.append(int(region_len)); child_to_parent.append(float(round(c_to_p,1))); region.append(temp_intvl); total_pumis.append(int(row['PTU'])); total_cumis.append(int(row['CTU'])); child_nums.append(row['CHILD_NUMS']); parent_freq.append(row['FREQ_PARENTS'])

	line = {'INTVL':region, 'PTU':total_pumis, 'CTU':total_cumis, 'CHILD_NUMS':child_nums, 'FREQ_PARENTS':parent_freq, 'INTVL_SIZE':size_of_intvl, 'CP':child_to_parent}
	df = pd.DataFrame(line, columns=df_headers)
	df.set_index('INTVL', inplace=True)

	headers_subdf=[]
	transp_table = [list(row) for row in itertools.zip_longest(*table, fillvalue=None)]
	col_nums = len(table)
	for i in range(1,col_nums+1):
		headers_subdf.append('col'+str(i))

	subframe = pd.DataFrame(transp_table, columns=headers_subdf)

	return df, subframe, name, col_nums, region

def plot_cp(df, output_path, name):
	#Plot Region vs. Child/Parent Ratio
	fig = plt.figure()
	df.sort_values('CP', ascending=False)['CP'].plot(kind='bar',x='INTVL',y='CP', color='pink', rot=90, title="Interval vs. Children to Parent Umis")
	plt.tight_layout()
	plt.savefig(output_path+"CP_"+name+".png")
	plt.close(fig)

def plot_PTU(df, output_path, name):
	#Plot Region vs. Parent Umi Count
	fig = plt.figure()
	df.sort_values('PTU', ascending=False)['PTU'].plot(kind='bar',x='INTVL',y='PTU', color='red', rot=90, title="Interval vs. Parent Umi Count")
	#plt.gcf().subplots_adjust(bottom=0.15)
	plt.tight_layout()
	plt.savefig(output_path+"PTU_"+name+".png")
	plt.close(fig)


def plot_CTU(df, output_path, name):
	#Plot Region vs. Child Umi Count
	fig = plt.figure()
	df.sort_values('CTU', ascending=False)['CTU'].plot(kind='bar',x='INTVL',y='CTU', color='blue', rot=90, title="Interval vs. Child Umi Count")
	plt.tight_layout()
	plt.savefig(output_path+"CTU_"+name+".png")
	plt.close(fig)

def plot_intvlsize_PTU_CTU(df, output_path, name):
	#Plot Interval size vs. Parent Umi Count & Child Umi Count
	ctu = df.plot(kind='scatter', x='INTVL_SIZE', y='CTU', color='blue', label="Child Umis Count")
	ptu = df.plot(kind='scatter', x='INTVL_SIZE', y='PTU', color = 'red', title="Interval Size vs. PTU and CTU", label="Parent Umi Count", ax=ctu)
	plt.legend()
	plt.tight_layout()
	plt.savefig(output_path+"CTU_PTU_intvlsize_"+name+".png")

def plot_child_pfreq(subframe, output_path, col_nums, regions):
	cnt = 0
	for i in range(0, col_nums, 2):
		#x_col = [subframe.columns[i]]
		#y_col = [subframe.columns[i+1]]

		subframe.plot(kind='scatter', x=i, y=i+1, color='purple', label="Parent Freq", title="No. of Children vs. Parent Freq.")
		plt.tight_layout()
		plt.savefig(output_path+"Children_vs_ParentFreq_"+str(regions[cnt])+".png")
		cnt+=1

def umi_plot(output_path, file_name, umi_flag):
	df, subframe, name, col_nums, regions = create_umi_dfs(file_name)

	if umi_flag == 'rs':
		plot_child_pfreq(subframe, output_path, col_nums, regions)
	elif umi_flag == 'all':
		plot_cp(df, output_path, name)
		plot_PTU(df, output_path, name)
		plot_CTU(df, output_path, name)
		plot_intvlsize_PTU_CTU(df, output_path, name)
		plot_child_pfreq(subframe, output_path, col_nums, regions)

		


#Consensus plots

def create_consdf(consfile):
	df_headers=['INTVL', 'CHROM', 'POS', 'REF', 'A', 'C', 'G', 'T', 'RAWDP', 'CONSDP', 'FAM', 'REF_FREQ', 'MEAN_FAM']
	df_headers2 = ['CHROM', 'POS', 'REF', 'A', 'C', 'G', 'T', 'RAWDP', 'CONSDP', 'FAM', 'REF_FREQ', 'MEAN_FAM']
	#df = pd.read_csv(consfile, sep='\t', columns=df_headers2)

	df = pd.read_csv(consfile, sep='\t')
	df.columns = ['CHROM', 'POS', 'REF', 'A', 'C', 'G', 'T', 'I', 'D', 'N', 'RAWDP', 'CONSDP', 'FAM', 'REF_FREQ', 'MEAN_FAM']
	return df

def plot_depth(df, output_path):
	#df = pd.DataFrame(line, columns=df_headers)
	#df.set_index('INTVL', inplace=True)
	figure(num=None, figsize=(15, 13), dpi=80, facecolor='w', edgecolor='k')
	groups=("zero","one", "two", "five")

	df.plot(kind='scatter', x='POS', y='CONSDP', color=['red', 'green', 'blue', 'purple'], title="Base position vs. Collapsed depth")
	#plt.legend()
	min_pos = min(df['POS'])
	max_pos = max(df['POS'])
	step_pos = (max_pos-min_pos)/3
	min_consdp = min(df['CONSDP'])
	max_consdp = max(df['CONSDP'])
	step_consdp = (max_consdp-min_consdp)/5
	plt.ylim([min_consdp,max_consdp])
	plt.xlim([min_pos, max_pos])
	plt.yticks(np.arange(min_consdp, max_consdp, step_consdp))
	plt.yscale('log')
	plt.xticks(np.arange(min_pos, max_pos, step_pos))
	plt.legend(('red', 'green', 'blue', 'purple'), ('0', '1', '2', '5'))
	plt.savefig(output_path+"base_pos_vs_CONSDP.png")



def cons_plot(output_path, file_name, cons_flag):

	if cons_flag == 'all':
		df = create_consdf(file_name)
		plot_depth(df, output_path)





