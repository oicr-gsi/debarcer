
import pysam
import ggplot as gg
import pandas as pd

def plot_umis(contig, region_start, region_end, cons_file):
	"""TODO"""

	df = pd.read_csv(cons_file, sep='\t')
	df = df[column.startswith('#')]

	##TEST
	print(df.to_string())

	## Plot A
	## ...

	## Plot B
	## ...

##TEST
if __name__=='__main__':
	plot_umis(None, None, None, '/.mounts/labs/gsiprojects/gsi/Debarcer2/debarcer2/data/epic-ds/chr1:16940000-16970000.cons')