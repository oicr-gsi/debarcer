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
from src.utilities import FormatRegion


#### functions for plotting coverage ####

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


def ExtractUmiCounts(DataFile):
    '''
    (file) -> dict    

    :param DataFile: Data file with umi count for a given region (ie. not merged)

    Return a dictionary with umi count for the different umi categories for a given region
    '''
    
    D = {}
    
    infile = open(DataFile)
    Header = infile.readline().strip().split('\t')
    line = infile.readline().strip()
    if line != '':
        line = line.split()
        # get genomic region
        chromo, start, end = line[Header.index('CHR')], line[Header.index('START')], line[Header.index('END')]
        # get total parent umis
        ptu = int(line[Header.index('PTU')])
        # get total child umis
        ctu = int(line[Header.index('CTU')])
        region = chromo + ':' + start + '-' + end
        D[region] = {'PTU': ptu, 'CTU': ctu} 
    infile.close()
    return D


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
        # extract region by reading coordinates in file
        region = FormatRegion(filename)
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


def PlotCoverage(directory, Outputfile):
    '''
    (str, str, dict) -> None
    
    :param directory: Directory containaing subdirectories Consfiles and Datafiles
                      respectively with consensus and data files
    :param Outputfile: Name of the output figure file 
         
    Generates a plot with mean coverage and total umis per interval
    
    Pre-condition: consensus and data files are not merged (chrN:A-B.cons and chrN:A-B.csv)
    '''
    
    # get the expected subdirectories in directory
    ConsDir = os.path.join(directory, 'Consfiles')
    DataDir = os.path.join(directory, 'Datafiles')
    
    if os.path.isdir(ConsDir) == False or os.path.isdir(DataDir) == False:
        raise ValueError('ERR: Invalid Consfiles and/or Datafiles directory')
    
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
    ax1 = CreateCoverageAx(1, 1, 1, figure, M, Coordinates, 'o', 'white', errorbar=S)
    ax2 = CreateCoverageAx(1, 1, 1, figure, Umis, Coordinates, 'o', 'black', firstax=ax1)
        
    plt.tight_layout()
    figure.savefig(Outputfile, bbox_inches = 'tight')


#### functions for plotting mean family size ####

def ExtractFamSize(ConsensusFile):
    '''
    (file) -> dict
    
    :param ConsensusFile: Path to consensus file 
    
    Return a dictionary of position: mean family size value pairs for each family
    for a given interval and sample
    '''
    
    # create a dict {fam: {pos: famsize}}
    D = {}
    
    infile = open(ConsensusFile)
    Header = infile.readline().rstrip().split('\t')
    for line in infile:
        if 'chr' in line:
            line = line.rstrip().split('\t')
            pos, fam, famsize = line[Header.index('POS')], int(line[Header.index('FAM')]), float(line[Header.index('MEAN_FAM')])
            # record data only for defined family size
            if fam != 0:
                # intialize inner dict
                if fam not in D:
                    D[fam] = {}
                # add pos, famsize for given family
                D[fam][pos] = famsize    
    infile.close()
    return D                
                    
  
def CreateMeanFamAx(Columns, Rows, Position, figure, Data, Color, YLabel, XLabel):
    '''
    (int, int, int, figure_object, dict, list, str, str) -> ax object
    
    :param columns: Number of columns
    :param rows: Number of rows
    :param position: Ax position in figure
    :param figure: Figure object opened for writing
    :param Data: mean family size at position for each family size threshold
    :param Color: sorted list of colors for plotting
    :param YLabel: Label of the Y axis
    :param XLabel: Label of the x axis
        
    Return a ax in figure
    '''
    
    # make a sorted list of family sizes
    FamSize = Data.keys()
    FamSize.sort()
    
    # make a sorted list of positions
    pos = list(map(lambda x: int(x), list(Data[FamSize[0]].keys())))
    pos.sort()
    pos = list(map(lambda x: str(x), pos))
    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot data  
    for i in range(len(FamSize)):
        ax.plot([j for j in range(len(pos))], [Data[FamSize[i]][j] for j in pos], color = Color[i], marker='', linewidth=2, linestyle='-', alpha = 1)
    
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
    legend_elements = []
    # loop over family sizes
    for i in range(len(FamSize)):
        legend_elements.append(Line2D([0], [0], marker='', label=str(FamSize[i]), linestyle='-', color = Color[i]))
    ax.legend(handles=legend_elements, frameon=False, ncol=5, loc='lower center')
    return ax


# use this function to create a figure for each consensus file for a given sample
def PlotMeanFamSize(ConsFile, Color, Outputfile):
    '''
    (str, list, str) -> None
    
    :param ConsFile: Path to the consensus file
    :param Color: List with colors for plotting
    :param Outputfile: Name of the output figure file
           
    Pre-condition: consensus file is not merged chrN:A-B.cons 
    '''
    
    # check that file is valid path
    if os.path.isfile(ConsFile) == False:
        raise ValueError('ERR: Invalid path to consensus file')
    
    # extract region from consensus file
    region = FormatRegion(ConsFile)
       
    # extract consensus depth for each family size -> {fam: {pos: meanfamSize}}
    Data = ExtractFamSize(ConsFile)
    
    # create figure object + ax
    figure = plt.figure(1, figsize = (10, 7))
    ax = CreateMeanFamAx(1, 1, 1, figure, Data, Color, 'Mean family size', region)
    
    plt.tight_layout()
    figure.savefig(Outputfile, bbox_inches = 'tight')
    
        

#### plot non-reference frequency ####

def ExtractNonRefFreq(ConsensusFile):
    '''
    (file) -> dict

    :param ConsensusFile: Path to consensus file

    Return a dictionary of interval: coverage value pairs
    
    Pre-condition: consensus file is not merged chrN:A-B.cons
    '''
    
    # create a dict {fam: {pos: non-ref freq}}
    D = {}
    
    infile = open(ConsensusFile)
    Header = infile.readline().rstrip().split('\t') 
    for line in infile:
        if 'chr' in line:
            line = line.rstrip().split('\t')
            pos, fam, freq = line[Header.index('POS')], int(line[Header.index('FAM')]), 100 - float(line[Header.index('REF_FREQ')])
            if fam not in D:
                D[fam] = {}
            D[fam][pos] = freq
    infile.close()
    return D                
                    
 
    
def CreateNonRefFreqAx(Columns, Rows, Position, figure, Data, Color, **Options):
    '''
    (int, int, int, figure_object, dict, str, dict) -> ax object
    
    :param columns: Number of columns
    :param rows: Number of rows
    :param position: Ax position in figure
    :param figure: Figure object opened for writing
    :param Data: Non-reference frequency at each position for a given family size
    :param Color: String color of the data
    :param Options: Accepted keys are:
                    'YLabel': Label of the Y axis
                    'XLabel': Label of the X axis
                    'legend': Add legend (True) or not (False)
                    'fam_size': List of family sizes in consensus file
                    'Colors': List of colors, parallel to fam_size 
                        
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
            FamSize = Options['fam_size']
            Colors = Options['colors']
            # add legend
            legend_elements = []
            for i in range(len(FamSize)):
                legend_elements.append(Line2D([0], [0], marker='s', label=str(FamSize[i]), linestyle='None', color = Colors[i]))
            ax.legend(handles=legend_elements, frameon=False, ncol=6, bbox_to_anchor=(0.7, 1.6))
    return ax


def PlotNonRefFreqData(ConsFile, Color, Outputfile):
    '''
    (str, list, str) -> None
    
    :param ConsFile: Path to the consensus file
    :param Color: List with colors for plotting
    :param Outputfile: Name of the output figure file
           
    Pre-condition: consensus file is not merged chrN:A-B.cons 
    '''
    
    # check that file is valid path
    if os.path.isfile(ConsFile) == False:
        raise ValueError('ERR: Invalid path to consensus file')
    
    # extract region from consensus file
    region = FormatRegion(ConsFile)
       
    # extract non-reference frequency for all family sizes in consensus file
    Data = ExtractNonRefFreq(ConsFile)
    
    # create figure
    figure = plt.figure(1, figsize = (8, 10))
    
    # make a sorted list of family sizes
    FamSize = list(Data.keys())
    FamSize.sort()
    
    # make a list of dicts {pos: non-ref freq}
    L = []
    for i in FamSize:
        d = {}
        for pos in Data[i]:
            d[pos] = Data[i][pos]
        L.append(d)
    
    for i in range(len(L)):
        if i == 0:
            ax = CreateNonRefFreqAx(1, len(L), i+1, figure, L[i], Color[i], legend=True, fam_size=FamSize, colors=Color)
        elif i == len(L) // 2:
            ax = CreateNonRefFreqAx(1, len(L), i+1, figure, L[i], Color[i], Ylabel='Non ref freq.')
        elif i == len(L) - 1:
            ax = CreateNonRefFreqAx(1, len(L), i+1, figure, L[i], Color[i], XLabel= region)
        else:
            ax = CreateNonRefFreqAx(1, len(L), i+1, figure, L[i], Color[i])
    plt.tight_layout()
    figure.savefig(Outputfile, bbox_inches = 'tight')
    


#### plot raw and consensus depth ####


def ExtractDepth(ConsensusFile):
    '''
    (file) -> dict
    
    :param ConsensusFile: Path to the consensus file
    
    Return a dictionary with coverage at each position for each family size
    '''
    
    D = {}
    
    infile = open(ConsensusFile)
    Header = infile.readline().rstrip().split('\t')
    for line in infile:
        if 'chr' in line:
            line = line.rstrip().split('\t')
            pos, rawdp = line[Header.index('POS')], int(line[Header.index('RAWDP')])
            consdp, fam = int(line[Header.index('CONSDP')]), int(line[Header.index('FAM')])
            if fam not in D:
                D[fam] = {}
            if fam == 0:
                 D[fam][pos] = rawdp
            else:
                D[fam][pos] = consdp
    infile.close()
    return D                
                    
  
def CreateConsDepthAx(Columns, Rows, Position, figure, Data, Color, YLabel, **Options):
    '''
    (int, int, int, figure_object, list, list, str, dict) -> ax object
    
    :param columns: Number of columns
    :param rows: Number of rows
    :param position: Ax position in figure
    :param figure: Figure object opened for writing
    :param Data: Depth at each position for each family size
    :param Color: List of colors
    :param Options: Accepted keys are:
                    'XLabel': Label of the X axis
                    'legend': Add legend (True) or not (False)
                    'fam_size': List of family sizes in consensus file
                    'Colors': List of colors, parallel to fam_size 
                        
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
            # get parallel lists of family sizes and colors
            Colors = Options['colors']
            FamSize = Options['fam_size']
            # add legend
            legend_elements = []
            for i in range(len(FamSize)):
                legend_elements.append(Line2D([0], [0], marker='', label=str(FamSize[i]), linestyle='-', color = Colors[i]))
            ax.legend(handles=legend_elements, frameon=False, ncol=5, loc='lower center')
    return ax


def PlotConsDepth(ConsFile, Color, Outputfile):
    '''
    (str, list, str) -> None
    
    :param ConsFile: Path to the consensus file
    :param Color: List with colors for plotting
    :param Outputfile: Name of the output figure file
           
    Pre-condition: consensus file is not merged chrN:A-B.cons 
    '''
    
    # check that file is valid path
    if os.path.isfile(ConsFile) == False:
        raise ValueError('ERR: Invalid path to consensus file')
       
    # extract region from consensus file
    region = FormatRegion(ConsFile)
        
    # extract depth at each position for each family size
    Data = ExtractDepth(ConsFile)
     
    # make a sorted list of dicts {pos: depth} for each family size
    FamSize = Data.keys()
    FamSize.sort()
    L = []
    for i in FamSize:
        d = {}
        for pos in Data[i]:
            d[pos] = Data[i][pos]
        L.append(d)
    
    figure = plt.figure(1, figsize = (10, 7))
    # plot raw depth, family size = 0    
    ax1 = CreateAx(1, 2, 1, figure, L[0:1], Color[0:1], 'Raw depth')
    ax2 = CreateAx(1, 2, 2, figure, L[1:], Color[1:], 'Consensus depth', XLabel=region, legend=True, fam_size=FamSize, colors=Color)
    plt.tight_layout()
    figure.savefig(Outputfile, bbox_inches = 'tight')
        
 
    

####### existing plots ######


def PlotChildParentRatio(directory, Outputfile):
    '''
    (str, str, dict) -> None
    
    :param directory: Directory containaing subdirectories Consfiles and Datafiles
                      respectively with consensus and data files
    :param Outputfile: Name of the output figure file 
         
    Generates a plot with mean children to parent umis count ratios for each region
    
    Pre-condition: consensus and data files are not merged (chrN:A-B.cons and chrN:A-B.csv)
    '''
    
    # get the directory with data files
    DataDir = os.path.join(directory, 'Datafiles')
    if os.path.isdir(DataDir) == False:
        raise ValueError('ERR: Invalid directory: {0}'.format(DataDir))
    
    # make a list of datafiles with umis
    DataFiles = [os.path.join(DataDir, i) for i in os.listdir(DataDir) if (i.startswith('datafile') and 'chr' in i and i[-4:] == '.csv')]
    
    
    print(DataFiles)
    
    
    # check that paths to files are valid
    for i in DataFiles:
        if os.path.isfile == False:
            raise ValueError('ERR: Invalid path to data file')
    
    # extract ptu and ctu for each region
    L = [ExtractUmiCounts(i) for i in DataFiles]
    # compute child/parent umi ratios for each region 
    Data = {}
    for d in L:
        region = list(d.keys())[0]
        ptu, ctu = d[region]['PTU'], d[region]['CTU']
        if ptu != 0:
            Data[region] = ctu/ptu
        
    # get a sorted list of positions
    Coordinates = SortPositions(list(Data.keys()))
    
    # create figure
    figure = plt.figure(1, figsize = (9, 6))
    # add a plot coverage to figure (N row, N column, plot N)
    ax = figure.add_subplot(1, 1, 1)
    # plot ctu/ptu ratio for each region
    ax.scatter([i for i in range(len(Coordinates))], [Data[i] for i in Coordinates], edgecolor = 'black', facecolor = 'pink', marker='o', lw = 1, s = 60, alpha = 1)
    
    # make a list of genomic regions 
    Chromos = []
    for i in Coordinates:
        i = i.split(':')
        Chromos.append(i[0] + '\n' + i[1].split('-')[0] + '\n' + i[1].split('-')[1])
    
    # limit y axis
    YMax = [Data[i] for i in Data]
    YMax = max(YMax)
    YMax = YMax + (YMax * 10/100)
    ax.set_ylim([0, YMax])
    step = round(YMax/10, 2)    
    # set y ticks    
    ax.yaxis.set_ticks([i for i in np.arange(0, YMax, step)])
        
    # write label for y axis
    ax.set_ylabel('Child:Parent Ratio', color = 'black',  size = 14, ha = 'center')
    ax.set_xlabel('Intervals', color = 'black',  size = 14, ha = 'center')
        
    # write title   
    ax.set_title("Interval vs. Children to Parent Umis", size = 14)
    
    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.4, linewidth = 0.4)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)
    
    # write label for x axis
    xPos = [i for i in range(len(Coordinates))]
    plt.xticks(xPos, Chromos, ha = 'center', rotation = 0, fontsize = 9)
               
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
    figure.savefig(Outputfile, bbox_inches = 'tight')




















#############################################################################################






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

    