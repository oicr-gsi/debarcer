# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 12:56:47 2019

@author: rjovelin
"""

# import modules
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import rc
rc('mathtext', default='regular')
import os
import numpy as np
from scipy import stats
import argparse




#### plot coverage. use this script to plot coverage for the different regions ###

def ExtractCoverage(ConsFile):
    '''
    (file) -> float, float
     
    :param ConsFile: Consensus file with raw depth at each position within a given region (ie. not merged)

    Return a tuple with the mean read depth within that interval and the standard error of the mean
    '''
    
    L = []
    infile = open(ConsFile)
    Header = infile.readline().rstrip().split('\t')
    for line in infile:
        if 'chr' in line:
            line = line.rstrip().split('\t')
            fam = line[Header.index('FAM')]
            if fam == '0':
                L.append(float(line[Header.index('RAWDP')]))
    infile.close()
    M = np.mean(L)
    sem = stats.sem(L)
    return M, sem


def ExtractTotalUmis(DataFile):
    '''
    (file) -> tuple    

    :param DataFile: Data file with total umi count for a given region (ie. not merged)

    Return a tuple with interval and total umis
    '''
    
    infile = open(DataFile)
    Header = infile.readline()
    line = infile.readline().strip()
    if line != '':
        line = line.split()
        chromo, start, end, umis = line[Header.index('CHR')], line[Header.index('START')], line[Header.index('END')], int(line[Header.index('PTU')])
    infile.close()
    return (chromo + ':' + start + '-' + end, umis)


def GetSampleUmis(L):
    '''
    (list) -> dict
    
    :param L: A list of full paths to data files with umi count per interval (ie. files not merged)
    
    Returns a dictionary of interval coordinates: total umis value pairs
    '''
    
    D = {}
    for filename in L:
        pos, umis = ExtractTotalUmis(filename)
        D[pos] = umis
    return D

def GetSampleCoverage(L):
    '''
    (list) -> dict
    
    :param L: A list of full paths to consensus files with umi count per interval (ie. files not merged)
    
    Returns a dictionary of interval coordinates with a list with mean and s.e.m. of coverage within the interval
    '''
    
    D = {}
    for filename in L:
        region = os.path.basename(filename)
        region = region[:region.index('.cons')]
        M, sem = ExtractCoverage(filename)
        D[region] = [M, sem]
    return D


def SortPositions(L):
    '''
    (list) -> list
    
    :param L: List of genomic intervals chrN:A-B (1-based inclusive)
    
    Returns a sorted list of intervals, sorted on chromosome and start positions 
    '''
     
    # make a sorted list of chromos
    Chromos = {}
    for i in L:
        i = i.split(':')
        chromo = i[0]
        start = i[1].split('-')[0]
        end = i[1].split('-')[1]
        if chromo in Chromos:
            Chromos[chromo].append((int(start), int(end)))
        else:
            Chromos[chromo] = [(int(start), int(end))]
    # sort positions within chromosomes
    for chromo in Chromos:
        Chromos[chromo].sort()
    
    # make a list of sorted chromosomes
    a = [int(i.replace('chr', '')) for i in Chromos if 'X' not in i]
    a.sort()
    for i in range(len(a)):
        a[i] = 'chr' + str(a[i])
    if 'chrX' in Chromos:
        a.append('chrX')
    
    Positions = []
    for i in a:
        for j in Chromos[i]:
            Positions.append(i + ':' + str(j[0]) + '-' + str(j[1]))
    return Positions                
             
             
             
def CreateCoverageAx(columns, rows, position, figure, data, coordinates, **Options):
    '''
    (int, int, int, figure_object, dict, list, str, str, dict) -> ax object
    
    :param columns: Number of columns
    :param rows: Number of rows
    :param position: Ax position in figure
    :param figure: Figure object opened for writing
    :param data: Dictionary of region: coverage or total umi key, value pairs
    :param coordinates: List of genomic intervals chrN:A-B
    :param Options: Accepted keys are:
                    'firstax': add a 2nd plot sharing axes of the 1st plot
                    'errorbar': error bars for the bar graph 
           
    Return a ax object in figure
    '''
    
    # plot total umi and coverage in a single plot if firstax option is used
    if 'firstax' in Options:
        # plot umi count using axis of 1st graph
        ax = Options['firstax'].twinx()
        ax.scatter([i for i in range(len(coordinates))], [data[i] for i in coordinates], edgecolor = 'black', facecolor = 'black', marker='o', lw = 1, s = 60, alpha = 1)
    else:
        # add a plot coverage to figure (N row, N column, plot N)
        ax = figure.add_subplot(rows, columns, position)
        # plot data
        if 'errorbar' in Options:
            errorbar = Options['errorbar']
            ax.bar([i for i in range(len(coordinates))], [data[i] for i in coordinates], width=0.8, yerr=errorbar,
                    color='white', edgecolor='black', linewidth=0.7, error_kw=dict(elinewidth=0.7, ecolor='black', markeredgewidth=0.7))
        else:
            ax.bar([i for i in range(len(coordinates))], [data[i] for i in coordinates], width=0.8, color='white', edgecolor='black', linewidth=0.7)
        
    # make a list of genomic regions 
    Chromos = []
    for i in coordinates:
        i = i.split(':')
        Chromos.append(i[0] + '\n' + i[1].split('-')[0] + '\n' + i[1].split('-')[1])
    
    # limit y axis
    YMax = [data[i] for i in data]
    YMax = max(YMax)
    YMax = float(YMax + (YMax * 10 /100))
    ax.set_ylim([0, YMax])    
    # set y ticks    
    if YMax <=50:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 10)])
    elif 50 < YMax <=200:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 20)]) 
    elif 200 < YMax <=500:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 50)])
    elif 500 < YMax <=1000:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 100)])  
    elif 1000 < YMax <=2000:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 200)])    
    elif 2000 < YMax <=5000:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 500)])
    elif 5000 < YMax <=10000:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 1000)])
    else:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 2000)])
    
    # set up y axis label and grid
    if 'firstax' not in Options:
        # write label for y axis
        ax.set_ylabel('Mean coverage per region', color = 'black',  size = 14, ha = 'center')
        ax.set_xlabel('Intervals', color = 'black',  size = 14, ha = 'center')
    
        # add a light grey horizontal grid to the plot, semi-transparent, 
        ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.4, linewidth = 0.4)  
        # hide these grids behind plot objects
        ax.set_axisbelow(True)

        # write label for x axis
        xPos = [i for i in range(len(coordinates))]
        #leftLim, rightLim = xPos[0] -1, xPos[-1] +1
        plt.xticks(xPos, Chromos, ha = 'center', rotation = 0, fontsize = 9)
    else:
        # write label for y axis
        ax.set_ylabel('Total umis per region', color = 'black',  size = 14, ha = 'center')
        
    # add space between axis and tick labels
    ax.yaxis.labelpad = 18
    ax.xaxis.labelpad = 18
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)  
       
    # do not show ticks
    plt.tick_params(axis='both', which='both', bottom=True, top=False,
                right=False, left=False, labelbottom=True, colors = 'black',
                labelsize = 12, direction = 'out')  
    
    # add legend
    legend_elements = [Line2D([0], [0], marker='s', markeredgecolor='black', markerfacecolor='white',
                       label='coverage', markersize=8, linestyle='None'),
                       Line2D([0], [0], marker='o', markeredgecolor='black', markerfacecolor='black',
                       label='umis', markersize=8, linestyle='None')]
    ax.legend(handles=legend_elements, frameon=False, ncol = 2, bbox_to_anchor=(0.9, 1.08))
    
    return ax


def PlotCoverage(directory, FigureFileName, extension):
    '''
    (str, str, dict) -> None
    
    :param directory: Directory containaing subdirectories Consfiles and Datafiles
                      respectively with consensus and data files
    :param FigureFileName: Name of the output .png figure file 
    :param extension: Figure format. Accepted values:
                      png, pdf, jpeg, tiff                                         
       
    Generates a plot with mean coverage and total umis per interval
    
    Pre-condition: consensus and data files are not merged (chrN:A-B.cons and chrN:A-B.csv)
    '''
    
    # get the expected subdirectories in directory
    ConsDir = os.path.join(directory, 'Consfiles')
    DataDir = os.path.join(directory, 'Datafiles')
    
    # make a list of consensus files
    ConsFiles = [os.path.join(ConsDir, i) for i in os.listdir(ConsDir) if i.startswith('chr') and i[-5:] == '.cons']
    # make a list of datafiles with umis
    DataFiles = [os.path.join(DataDir, i) for i in os.listdir(DataDir) if (i.startswith('datafile') and 'chr' in i and i[-4:] == '.csv')]
    
    # check that paths to files are valid
    for i in ConsFiles:
        if os.path.isfile == False:
            raise ValueError('ERR: Invalid path to consensus file')
    for i in DataFiles:
        if os.path.isfile == False:
            raise ValueError('ERR: Invalid path to data file')
    
    # get mean coverage per interval
    Coverage = GetSampleCoverage(ConsFiles)
    # get total parent umis for each interval
    Umis = GetSampleUmis(DataFiles)
    # get a sorted list of positions
    Coordinates = SortPositions(list(Coverage.keys()))
    
    # create figure
    figure = plt.figure(1, figsize = (9, 6))
    # create a dict with mean coverage
    M = {}
    for i in Coverage:
        M[i] = Coverage[i][0]
    # create a sorted list with sem
    S = [Coverage[i][1] for i in Coordinates]
    
    # plot data
    ax1 = CreateAx(1, 1, 1, figure, M, Coordinates, 'o', 'white', errorbar=S)
    ax2 = CreateAx(1, 1, 1, figure, Umis, Coordinates, 'o', 'black', firstax=ax1)
        
    # save figure in Figures directory, create directory if doesn't exist    
    FigDir = os.path.join(directory, 'Figures')
    if os.path.isdir(FigDir) == False:
         os.mkdir(FigDir)
    Outputfile = os.path.join(FigDir, FigureFileName + extension)
    plt.tight_layout()
    figure.savefig(Outputfile, bbox_inches = 'tight')

























  
    
#### plot mean fam size ####




def ExtractFamSize(ConsensusFile, Family):
    '''
    (file) -> dict
    Take a consensus file for a given interval and sample and return a dictionary
    of position: mean family size value pairs
    '''
    
    D = {}
    
    infile = open(ConsensusFile)
    Header = infile.readline().rstrip().split('\t')
    for line in infile:
        if 'chr' in line:
            line = line.rstrip().split('\t')
            pos, fam, famsize = line[Header.index('POS')], line[Header.index('FAM')], float(line[Header.index('MEAN_FAM')])
            if Family == fam:
                assert pos not in D
                D[pos] = famsize
    infile.close()
    return D                
                    
  
def CreateAx(Columns, Rows, Position, figure, Data, Color, YLabel, XLabel, Title):
    '''
    (int, int, int, figure_object, list, dict, str, dict) -> ax object
    Return a ax in figure
    '''
    
    # make a sorted list of positions
    pos = list(map(lambda x: int(x), list(Data[0].keys())))
    pos.sort()
    pos = list(map(lambda x: str(x), pos))
    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot data
    for i in range(len(Data)):
        ax.plot([j for j in range(len(pos))], [Data[i][j] for j in pos], color = Color[i], marker='', linewidth=2, linestyle='-', alpha = 1)
    
    # limit y axis
    YMax = []
    for i in range(len(Data)):
        for j in Data[i]:
            YMax.append(Data[i][j])
    YMax = max(YMax)
    if YMax < 10:
        YMax = float(YMax + (YMax * 25 /100))
    else:
        YMax = float(YMax + (YMax * 15 /100))
    
    
    ax.set_ylim([0, YMax])    
        
    # write title   
    ax.set_title(Title, size = 14)
    
    # write label for y and x axis
    ax.set_ylabel(YLabel, color = 'black',  size = 14, ha = 'center')
    ax.set_xlabel(XLabel, color = 'black',  size = 14, ha = 'center')

    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.4, linewidth = 0.4)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)

    # add space between axis and tick labels
    ax.yaxis.labelpad = 18
    ax.xaxis.labelpad = 18
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)  
    # offset the spines
#    for spine in ax.spines.values():
#        spine.set_position(('outward', 10))
    
    # do not show ticks
    plt.tick_params(axis='both', which='both', bottom=True, top=False,
                right=False, left=False, labelbottom=True, colors = 'black',
                labelsize = 12, direction = 'out')  

    # write ticks for y and x axis
    xtickspos = [i for i in range(0, len(pos), 20)]
    xticks = [pos[i] for i in xtickspos]
    plt.xticks(xtickspos, xticks, ha = 'center', rotation = 0, fontsize = 12)
    ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 2)])
    
    # add legend
    legend_elements = [Line2D([0], [0], marker='', label='1', linestyle='-', color = '#4B0082'),
                       Line2D([0], [0], marker='', label='2', linestyle='-', color = '#7B68EE'),
                       Line2D([0], [0], marker='', label='3', linestyle='-', color = '#8A2BE2'),
                       Line2D([0], [0], marker='', label='5', linestyle='-', color = '#BA55D3'),
                       Line2D([0], [0], marker='', label='10', linestyle='-', color = '#DDA0DD')]       
    ax.legend(handles=legend_elements, frameon=False, ncol=5, loc='lower center')
    return ax





# use this function to create a figure for each consensus file for a given sample
def PlotData(args):
    '''
    
    '''
    
    # map sample to pool
    Pools = {'SWID_14058644':'Pool 1-1', 'SWID_14058646':'Pool 1-2', 'SWID_14058648':'Pool 2-1',
    'SWID_14058650':'Pool 2-2', 'SWID_14058652':'Pool 4-1', 'SWID_14058654':'Pool 4-2',
    'SWID_14058656':'Pool 5-1', 'SWID_14058658':'Pool 5-2', 'SWID_14058660':'Pool 6-1',
    'SWID_14058662':'Pool 6-2', 'SWID_14058664':'Pool 7-1', 'SWID_14058666':'Pool 7-2'}

    # get sample directory
    sampledir = os.path.join(args.workingdir, args.sample)
    assert os.path.isdir(sampledir)
    
    # create directory to save figures
    FigDir = os.path.join(sampledir, 'Figures')
    if os.path.isdir(FigDir) == False:
        os.mkdir(FigDir)
    
    ConsFile = os.path.join(sampledir, args.region + '.cons')
    assert os.path.isfile(ConsFile)
    
    # get outputfile name
    OutputFile = args.sample + '_MeanFamilySize_' + args.region.replace(':', '-') + '.png'
    print(OutputFile)
    OutputFile = os.path.join(FigDir, OutputFile)
      
    # extract consesnsus depth
    C = [ExtractFamSize(ConsFile, j) for j in ['1', '2', '3', '5', '10']]
        
    Color = {0: '#4B0082', 1: '#7B68EE', 2: '#8A2BE2', 3: '#BA55D3', 4:'#DDA0DD'}
    figure = plt.figure(1, figsize = (10, 7))
        
    ax = CreateAx(1, 1, 1, figure, C, Color, 'Mean family size', args.region, Pools[args.sample])
    figure.savefig(OutputFile, bbox_inches = 'tight')
        
        
if __name__ == '__main__':

    # create top-level parser
    main_parser = argparse.ArgumentParser(prog = 'PlotMeanFamSize.py', description='Plot mean family size for a given sample and region', add_help=True)

    # form analyses to EGA       
    main_parser.add_argument('-w', '--WorkingDir', dest='workingdir', default='/.mounts/labs/gsiprojects/genomics/CBALL/GroupCollapse/',
                             help='Directory with sample directories containing consensus files. Default is /.mounts/labs/gsiprojects/genomics/CBALL/GroupCollapse/')
    main_parser.add_argument('-s', '--Sample', dest='sample', help='Sample name', required=True)
    main_parser.add_argument('-r', '--Region', dest='region', help='Region', required=True)
    main_parser.set_defaults(func=PlotData)

    # get arguments from the command line
    args = main_parser.parse_args()
    # pass the args to the default function
    args.func(args)




##### plot non-ref #####


# -*- coding: utf-8 -*-
"""
Created on Tue May 28 15:12:00 2019

@author: rjovelin
"""

# use this script to plot non-reference frequency for a given interval and sample

# import modules
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import rc
rc('mathtext', default='regular')
import os
import argparse
import numpy as np


def ExtractNonRefFreq(ConsensusFile, Family):
    '''
    (file) -> dict
    Take a file with mean coverage per interval and return a dictionary
    of interval: coverage value pairs
    '''
    
    D = {}
    
    infile = open(ConsensusFile)
    Header = infile.readline().rstrip().split('\t') 
    for line in infile:
        if 'chr' in line:
            line = line.rstrip().split('\t')
            pos, fam, freq = line[Header.index('POS')], line[Header.index('FAM')], 100 - float(line[Header.index('REF_FREQ')])
            if Family == fam:
                assert pos not in D
                D[pos] = freq
    infile.close()
    return D                
                    
  
def CreateAx(Columns, Rows, Position, figure, Data, Color, **Options):
    '''
    (int, int, int, figure_object, dict, str, str, dict) -> ax object
    Return a ax in figure
    '''
    
    # make a sorted list of positions
    pos = list(map(lambda x: int(x), list(Data.keys())))
    pos.sort()
    pos = list(map(lambda x: str(x), pos))
    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot data
    ax.bar([j for j in range(len(pos))], [Data[j] for j in pos], color=Color, edgecolor=Color, linewidth=0, alpha = 1)
    
    # limit y axis
    YMax = [Data[i] for i in Data]
    YMax = max(YMax)
    YMax = float(YMax + (YMax * 10 /100))
    ax.set_ylim([0, YMax])    
        
    # write y axis ticks
    if YMax <=2:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 0.2)])
    elif 2 < YMax <=6:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 1)])
    elif 6 < YMax <=10:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 2)])
    elif 10 < YMax <=30:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 5)])
    else:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 10)])
    
    # write title   
    if 'title' in Options:
        Title = Options['title']
        ax.set_title(Title, size = 14)
    
    # write label for y and x axis
    if 'YLabel' in Options:
        YLabel=Options['YLabel']
        ax.set_ylabel(YLabel, color = 'black',  size = 14, ha = 'center')
    if 'XLabel' in Options:
        XLabel = Options['XLabel']
        ax.set_xlabel(XLabel, color = 'black',  size = 14, ha = 'center')

    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.4, linewidth = 0.4)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)

    # add space between axis and tick labels
    ax.yaxis.labelpad = 18
    ax.xaxis.labelpad = 18
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)  
    # offset the spines
#    for spine in ax.spines.values():
#        spine.set_position(('outward', 10))
    
    if 'XLabel' in Options:
        # do not show ticks
        plt.tick_params(axis='both', which='both', bottom=True, top=False,
                    right=False, left=False, labelbottom=True, colors = 'black',
                    labelsize = 12, direction = 'out')  

        # write ticks for x axis
        xtickspos = [i for i in range(0, len(pos), 20)]
        xticks = [pos[i] for i in xtickspos]
        plt.xticks(xtickspos, xticks, ha = 'center', rotation = 0, fontsize = 12)
    else:
        # do not show ticks
        plt.tick_params(axis='both', which='both', bottom=True, top=False,
                    right=False, left=False, labelbottom=False, colors = 'black',
                    labelsize = 12, direction = 'out')  
    if 'legend' in Options:
        if Options['legend'] == True:
            # add legend
            legend_elements = [Line2D([0], [0], marker='s', label='0', linestyle='None', color = 'black'),
                               Line2D([0], [0], marker='s', label='1', linestyle='None', color = '#4B0082'),
                               Line2D([0], [0], marker='s', label='2', linestyle='None', color = '#7B68EE'),
                               Line2D([0], [0], marker='s', label='3', linestyle='None', color = '#8A2BE2'),
                               Line2D([0], [0], marker='s', label='5', linestyle='None', color = '#BA55D3'),
                               Line2D([0], [0], marker='s', label='10', linestyle='None', color = '#DDA0DD')]       
            ax.legend(handles=legend_elements, frameon=False, ncol=6, bbox_to_anchor=(0.7, 1.6))
    return ax


# use this function to create a figure for each consensus file for a given sample
def PlotData(args):
    '''
    
    '''
    
    # map sample to pool
    Pools = {'SWID_14058644':'Pool 1-1', 'SWID_14058646':'Pool 1-2', 'SWID_14058648':'Pool 2-1',
    'SWID_14058650':'Pool 2-2', 'SWID_14058652':'Pool 4-1', 'SWID_14058654':'Pool 4-2',
    'SWID_14058656':'Pool 5-1', 'SWID_14058658':'Pool 5-2', 'SWID_14058660':'Pool 6-1',
    'SWID_14058662':'Pool 6-2', 'SWID_14058664':'Pool 7-1', 'SWID_14058666':'Pool 7-2'}

    # get sample directory
    sampledir = os.path.join(args.workingdir, args.sample)
    assert os.path.isdir(sampledir)
    
    # create directory to save figures
    FigDir = os.path.join(sampledir, 'Figures')
    if os.path.isdir(FigDir) == False:
        os.mkdir(FigDir)
    
    ConsFile = os.path.join(sampledir, args.region + '.cons')
    assert os.path.isfile(ConsFile)
    
    # get outputfile name
    OutputFile = args.sample + '_NonRefFreq_' + args.region.replace(':', '-') + '.png'
    OutputFile = os.path.join(FigDir, OutputFile)
      
    # extract non-reference frequency
    Data = [ExtractNonRefFreq(ConsFile, j) for j in ['0', '1', '2', '3', '5', '10']]
    
    # set up colors  
    Color = {0:'black', 1: '#4B0082', 2: '#7B68EE', 3: '#8A2BE2', 4: '#BA55D3', 5:'#DDA0DD'}
    
    # create figure
    figure = plt.figure(1, figsize = (8, 10))
    for i in range(len(Data)):
        if i == 0:
            ax = CreateAx(1, 6, i+1, figure, Data[i], Color[i], title=Pools[args.sample], legend=True)
        elif i == 3:
            ax = CreateAx(1, 6, i+1, figure, Data[i], Color[i], Ylabel='Non ref freq.')
        elif i == 5:
            ax = CreateAx(1, 6, i+1, figure, Data[i], Color[i], XLabel=args.region)
        else:
            ax = CreateAx(1, 6, i+1, figure, Data[i], Color[i])
    plt.tight_layout()
    figure.savefig(OutputFile, bbox_inches = 'tight')
        
        
if __name__ == '__main__':

    # create top-level parser
    main_parser = argparse.ArgumentParser(prog = 'PlotNonRef.py', description='Plot non-reference frequency for a given sample and region', add_help=True)

    # form analyses to EGA       
    main_parser.add_argument('-w', '--WorkingDir', dest='workingdir', default='/.mounts/labs/gsiprojects/genomics/CBALL/GroupCollapse/',
                             help='Directory with sample directories containing consensus files. Default is /.mounts/labs/gsiprojects/genomics/CBALL/GroupCollapse/')
    main_parser.add_argument('-s', '--Sample', dest='sample', help='Sample name', required=True)
    main_parser.add_argument('-r', '--Region', dest='region', help='Region', required=True)
    main_parser.set_defaults(func=PlotData)

    # get arguments from the command line
    args = main_parser.parse_args()
    # pass the args to the default function
    args.func(args)



#### plot rawcons depth ###


# -*- coding: utf-8 -*-
"""
Created on Mon May 27 11:41:50 2019

@author: rjovelin
"""


# use this script to plot coverage for the different regions

# import modules
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import rc
rc('mathtext', default='regular')
import os
import argparse
import numpy as np


def ExtractDepth(ConsensusFile, Family):
    '''
    (file) -> dict
    Take a file with mean coverage per interval and return a dictionary
    of interval: coverage value pairs
    '''
    
    D = {}
    
    infile = open(ConsensusFile)
    Header = infile.readline().rstrip().split('\t')
    for line in infile:
        if 'chr' in line:
            line = line.rstrip().split('\t')
            pos, rawdp = line[Header.index('POS')], int(line[Header.index('RAWDP')])
            consdp, fam = int(line[Header.index('CONSDP')]), line[Header.index('FAM')]
            if Family == fam:
                assert pos not in D
                if Family == '0':
                    D[pos] = rawdp
                else:
                    D[pos] = consdp
    infile.close()
    return D                
                    
  
def CreateAx(Columns, Rows, Position, figure, Data, Color, YLabel, **Options):
    '''
    (int, int, int, figure_object, list, dict, str, dict) -> ax object
    Return a ax in figure
    '''
    
    # make a sorted list of positions
    pos = list(map(lambda x: int(x), list(Data[0].keys())))
    pos.sort()
    pos = list(map(lambda x: str(x), pos))
    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot data
    for i in range(len(Data)):
        ax.plot([j for j in range(len(pos))], [Data[i][j] for j in pos], color = Color[i], marker='', linewidth=2, linestyle='-', alpha = 1)
    
    # limit y axis
    YMax = []
    for i in range(len(Data)):
        for j in Data[i]:
            YMax.append(Data[i][j])
    YMax = max(YMax)
    YMax = float(YMax + (YMax * 10 /100))
    ax.set_ylim([0, YMax])    
        
    if YMax <=50:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 10)])
    elif 50 < YMax <=200:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 20)]) 
    elif 200 < YMax <=500:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 50)])
    elif 500 < YMax <=1000:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 100)])  
    elif 1000 < YMax <=2000:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 200)])    
    elif 2000 < YMax <=5000:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 500)])
    elif 5000 < YMax <=10000:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 1000)])
    else:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 2000)])
    
    # write title   
    if 'title' in Options:
        Title = Options['title']
        ax.set_title(Title, size = 14)
    
    # write label for y and x axis
    ax.set_ylabel(YLabel, color = 'black',  size = 14, ha = 'center')
    if 'XLabel' in Options:
        XLabel = Options['XLabel']
        ax.set_xlabel(XLabel, color = 'black',  size = 14, ha = 'center')

    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.4, linewidth = 0.4)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)

    # add space between axis and tick labels
    ax.yaxis.labelpad = 18
    ax.xaxis.labelpad = 18
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)  
    # offset the spines
#    for spine in ax.spines.values():
#        spine.set_position(('outward', 10))
    
    if 'XLabel' in Options:
        # do not show ticks
        plt.tick_params(axis='both', which='both', bottom=True, top=False,
                    right=False, left=False, labelbottom=True, colors = 'black',
                    labelsize = 12, direction = 'out')  

        # write ticks for y and x axis
        xtickspos = [i for i in range(0, len(pos), 20)]
        xticks = [pos[i] for i in xtickspos]
        plt.xticks(xtickspos, xticks, ha = 'center', rotation = 0, fontsize = 12)
    else:
        # do not show ticks
        plt.tick_params(axis='both', which='both', bottom=True, top=False,
                    right=False, left=False, labelbottom=False, colors = 'black',
                    labelsize = 12, direction = 'out')  
    if 'legend' in Options:
        if Options['legend'] == True:
            # add legend
            legend_elements = [Line2D([0], [0], marker='', label='1', linestyle='-', color = '#4B0082'),
                               Line2D([0], [0], marker='', label='2', linestyle='-', color = '#7B68EE'),
                               Line2D([0], [0], marker='', label='3', linestyle='-', color = '#8A2BE2'),
                               Line2D([0], [0], marker='', label='5', linestyle='-', color = '#BA55D3'),
                               Line2D([0], [0], marker='', label='10', linestyle='-', color = '#DDA0DD')]       
            ax.legend(handles=legend_elements, frameon=False, ncol=5, loc='lower center')
    return ax





# use this function to create a figure for each consensus file for a given sample
def PlotData(args):
    '''
    
    '''
    
    # map sample to pool
    Pools = {'SWID_14058644':'Pool 1-1', 'SWID_14058646':'Pool 1-2', 'SWID_14058648':'Pool 2-1',
    'SWID_14058650':'Pool 2-2', 'SWID_14058652':'Pool 4-1', 'SWID_14058654':'Pool 4-2',
    'SWID_14058656':'Pool 5-1', 'SWID_14058658':'Pool 5-2', 'SWID_14058660':'Pool 6-1',
    'SWID_14058662':'Pool 6-2', 'SWID_14058664':'Pool 7-1', 'SWID_14058666':'Pool 7-2'}

    # get sample directory
    sampledir = os.path.join(args.workingdir, args.sample)
    assert os.path.isdir(sampledir)
    
    # create directory to save figures
    FigDir = os.path.join(sampledir, 'Figures')
    if os.path.isdir(FigDir) == False:
        os.mkdir(FigDir)
    
    ConsFile = os.path.join(sampledir, args.region + '.cons')
    assert os.path.isfile(ConsFile)
    
    # get outputfile name
    OutputFile = args.sample + '_RawConsensusDepth_' + args.region.replace(':', '-') + '.png'
    print(OutputFile)
    OutputFile = os.path.join(FigDir, OutputFile)
      
    # extract rawdepth  data from file
    R = [ExtractDepth(ConsFile, '0')]
    # extract consesnsus depth
    C = [ExtractDepth(ConsFile, j) for j in ['1', '2', '3', '5', '10']]
        
    rColor = {0:'black'}
    cColor = {0: '#4B0082', 1: '#7B68EE', 2: '#8A2BE2', 3: '#BA55D3', 4:'#DDA0DD'}
    figure = plt.figure(1, figsize = (10, 7))
        
    ax1 = CreateAx(1, 2, 1, figure, R, rColor, 'Raw depth', title=Pools[args.sample])
    ax2 = CreateAx(1, 2, 2, figure, C, cColor, 'Consensus depth', XLabel=args.region, legend=True)
    plt.tight_layout()
    figure.savefig(OutputFile, bbox_inches = 'tight')
        
        
if __name__ == '__main__':

    # create top-level parser
    main_parser = argparse.ArgumentParser(prog = 'PlotRawConsDepth.py', description='Plot raw and consensus depth for a given sample and region', add_help=True)

    # form analyses to EGA       
    main_parser.add_argument('-w', '--WorkingDir', dest='workingdir', default='/.mounts/labs/gsiprojects/genomics/CBALL/GroupCollapse/',
                             help='Directory with sample directories containing consensus files. Default is /.mounts/labs/gsiprojects/genomics/CBALL/GroupCollapse/')
    main_parser.add_argument('-s', '--Sample', dest='sample', help='Sample name', required=True)
    main_parser.add_argument('-r', '--Region', dest='region', help='Region', required=True)
    main_parser.set_defaults(func=PlotData)

    # get arguments from the command line
    args = main_parser.parse_args()
    # pass the args to the default function
    args.func(args)



####### existing plots ######


import os
import sys
import matplotlib
matplotlib.use('Agg')
import os
import sys
import pysam
import configparser
import argparse
import operator
import functools
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
import csv
import fnmatch
import itertools
from matplotlib.pyplot import figure

#import plotly.graph_objs as go
#import plotly.offline as off

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
	plt.xlabel('Interval')
	plt.ylabel('Child:Parent Ratio')
	plt.tight_layout()
	plt.savefig(output_path+"CP_"+name+".png")
	plt.close(fig)

def plot_PTU(df, output_path, name):
	#Plot Region vs. Parent Umi Count
	fig = plt.figure()
	df.sort_values('PTU', ascending=False)['PTU'].plot(kind='bar',x='INTVL',y='PTU', color='red', rot=90, title="Interval vs. Parent Umi Count")
	#plt.gcf().subplots_adjust(bottom=0.15)
	plt.xlabel('Interval')
	plt.ylabel('Number of Parent UMIs')
	plt.tight_layout()
	plt.savefig(output_path+"PTU_"+name+".png")
	plt.close(fig)


def plot_CTU(df, output_path, name):
	#Plot Region vs. Child Umi Count
	fig = plt.figure()
	df.sort_values('CTU', ascending=False)['CTU'].plot(kind='bar',x='INTVL',y='CTU', color='blue', rot=90, title="Interval vs. Child Umi Count")
	plt.xlabel('Interval')
	plt.ylabel('Number of Child UMIs')
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
		plt.xlabel('Number of UMI Children')
		plt.ylabel('Frequency of Parents')

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
	figure(num=None, figsize=(15, 13), dpi=80, facecolor='w', edgecolor='k')
	groups=("zero","one", "two", "five")

	colors = ['blue', 'green', 'red', 'purple']
	ax = plt.scatter(x, y, c=label, cmap=matplotlib.colors.ListedColormap(colors))

	#plt.legend()
	plt.yscale('log')
	plt.xlim([min_pos, max_pos])
	plt.yticks([100, 1000, 10000, 100000, 1000000])

	plt.xticks(np.arange(min_pos, max_pos, step_pos))
	plt.ticklabel_format(useOffset=False, style='plain', axis='x')
	plt.xlabel = "Base Position"
	plt.ylabel = "Depth"

	plt.savefig(output_path+"base_pos_vs_CONSDP.png")


def plot_reffeq(df, output_path):

	figure(num=None, figsize=(15, 13), dpi=80, facecolor='w', edgecolor='k')

	groups=("zero","one", "two", "five")

	x = df['POS']
	y = df['REF_FREQ']
	label=df['FAM']

	colors = ['blue', 'green', 'red', 'purple']
	ax = plt.scatter(x, y, c=label, cmap=matplotlib.colors.ListedColormap(colors))

	#plt.legend(label, colors)

	min_pos = min(df['POS'])
	max_pos = max(df['POS'])
	step_pos = (max_pos-min_pos)/5

	min_reffreq = min(df['REF_FREQ'])
	max_reffreq = max(df['REF_FREQ'])
	step_reffreq = (max_reffreq-min_reffreq)/5

	plt.xlim([min_pos, max_pos])
	plt.yticks(np.arange(min_reffreq, max_reffreq, step_reffreq))

	plt.xticks(np.arange(min_pos, max_pos, step_pos))
	plt.ticklabel_format(useOffset=False, style='plain', axis='x')
	plt.xlabel = "Base Position"
	plt.ylabel = "Refrence Frequency"

	plt.savefig(output_path+"base_pos_vs_REFFREQ.png")






def cons_plot(output_path, file_name, cons_flag):

	if cons_flag == 'all':
		df = create_consdf(file_name)
		plot_depth(df, output_path)
		plot_reffeq(df, output_path)




### need to edit specript below to allow commands for each plot



if __name__ == '__main__':

    # create top-level parser
    main_parser = argparse.ArgumentParser(prog = 'PlotCoverage.py', description='Plot mean coverage and total umis per interval', add_help=True)

    # form analyses to EGA       
    main_parser.add_argument('-w', '--WorkingDir', dest='workingdir', default='/.mounts/labs/gsiprojects/genomics/CBALL/GroupCollapse/',
                             help='Directory with sample directories containing consensus files. Default is /.mounts/labs/gsiprojects/genomics/CBALL/GroupCollapse/')
    main_parser.add_argument('-s', '--Sample', dest='sample', help='Sample name', required=True)
    main_parser.set_defaults(func=PlotData)

    # get arguments from the command line
    args = main_parser.parse_args()
    # pass the args to the default function
    args.func(args)

    