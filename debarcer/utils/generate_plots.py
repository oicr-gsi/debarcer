
import matplotlib
from matplotlib.ticker import MultipleLocator
matplotlib.use('Agg')

import pandas as pd
from ggplot import *
import argparse

def tally_output(cons_file, output_path, file_stem):
	"""Writes a tally file for easy graphing/analysis.""" 

	tally_file = "{}/{}.tally".format(output_path, file_stem)

	with open(cons_file, "r") as reader, open(tally_file, "w") as writer:

		for line in reader:
			if line.startswith('#'):
				writer.write(line.split('# ')[1])

	return tally_file


def generate_df(tally_file):
	"""Builds a pandas DataFrame from the given tally_file."""

	df = pd.read_csv(tally_file, sep='\t')
	df.columns = ['chr', 'base_pos', 'A', 'C', 'G', 'T', 'I', 'D', 'N', 'raw_dp', 'cons_dp', 'f_size', 'ref_freq']

	return df


def graph_depth(df, output_path, file_stem):
	"""A plot of base position vs. depth, coloured by family size."""

	plot_file = "{}/{}_depth.png".format(output_path, file_stem)
	
	base_pos_min = min(df['base_pos'])
	base_pos_max = max(df['base_pos'])
	base_pos_step = int((base_pos_max - base_pos_min) / 10)

	depth_max = max(df['cons_dp'])
	depth_step = int(depth_max / 10)

	plot = ggplot(aes('base_pos', 'cons_dp'), data=df) + \
		geom_point(aes(colour='f_size')) + \
		xlab("Base Position") + \
		ylim(0, 2000) + \
		scale_x_continuous(breaks=range(base_pos_min, base_pos_max, base_pos_step)) + \
		scale_y_continuous(breaks=range(0, depth_max, depth_step)) + \
		ylab("Depth")

	plot.save(plot_file)


def graph_ref_freq(df, output_path, file_stem):
	"""A plot of base position vs. reference frequency, faceted by family size."""

	plot_file = "{}/{}_freq.png".format(output_path, file_stem)

	base_pos_min = min(df['base_pos'])
	base_pos_max = max(df['base_pos'])
	base_pos_step = int((base_pos_max - base_pos_min) / 10)

	for x in range(base_pos_min, base_pos_max, base_pos_step):
		print(x)

	plot = ggplot(aes('base_pos', 'ref_freq'), data=df) + \
		geom_bar(aes(colour='f_size')) + \
		facet_wrap('f_size', ncol=1) + \
		scale_x_continuous(breaks=range(base_pos_min, base_pos_max, base_pos_step)) + \
		scale_y_continuous(breaks=[0.0, 25.0, 50.0, 75.0, 100.0]) + \
		xlab("Base Position") + \
		ylab("Reference Frequency (%)")

	plot.save(plot_file)


def generate_plots(cons_file, output_path):
	"""(Main) generates plots from the given .cons files."""

	file_stem = cons_file.split('/')[-1].split('.cons')[0]
	tally_file = tally_output(cons_file, output_path, file_stem)
	df = generate_df(tally_file)
	graph_depth(df, output_path, file_stem)


if __name__ == '__main__':

	## Argument + config parsing and error handling
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cons_file', help='Path to your .cons file.')
    parser.add_argument('-o', '--output_path', help='Path to write output files to.')

    args = parser.parse_args()

    generate_plots(args.cons_file, args.output_path)

	


