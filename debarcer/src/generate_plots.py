
import pysam
import ggplot as gg
import pandas as pd

def plot_umis(contig, region_start, region_end, cons_file):
	"""TODO"""

	cons_df = pd.read_csv(cons_file, sep='\t')

	##TEST
	print(cons_df.to_string())

	## Plot A
	## ...

	## Plot B
	## ...

##TEST
if __name__=='__main__':
	plot_umis(None, None, None, '/.mounts/labs/gsiprojects/gsi/Debarcer2/debarcer2/data/epic-ds/chr1:16959000-16960000.cons')