# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 12:56:47 2019

@author: rjovelin
"""

# import modules
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy import stats
from src.utilities import edit_distance, FormatRegion
from src.umi_error_correct import most_frequent
import networkx as nx
import json
import collections
import pandas as pd
import seaborn as sns
import pygal
from pygal.style import DefaultStyle, Style
import yaml
import math

def SetUpTicks(AxisMax):
    '''
    (num) -> int
    
    :param AxisMax: maximum value of x or y axis
    
    Return the step to set up axis ticks
    '''
    
    if AxisMax <= 10:
        step = 1 
    elif 10 < AxisMax <= 20:
        step = 2
    elif 20 < AxisMax <=50:
        step = 10
    elif 50 < AxisMax <=200:
        step = 20  
    elif 200 < AxisMax <=500:
        step = 50
    elif 500 < AxisMax <=1000:
        step = 100   
    elif 1000 < AxisMax <=2000:
        step = 200     
    elif 2000 < AxisMax <=5000:
        step = 500 
    elif 5000 < AxisMax <=10000:
        step = 1000 
    elif 10000 < AxisMax <= 50000:
        step = 5000
    elif 50000 < AxisMax <= 200000:
        step = 25000
    elif 200000 < AxisMax <= 1000000:
        step = 100000
    elif 1000000 < AxisMax <= 5000000:
        step = 500000
    else:
        step = 1000000
    return step


def GetCoverageRegions(CoverageStats):
    '''
    (str) -> dict
    
    :param CoverageStats: yaml file with mean read depth per region in Stats directory
    
    Returns a dictionary of region: mean read depth 
    '''
    
    with open(CoverageStats) as infile:
        coverage = yaml.safe_load(infile)
    return coverage 
      

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


def ExtractUmiCounts(DataFile):
    '''
    (file) -> dict    

    :param DataFile: Data file with umi count for a given region (ie. not merged and not empty)

    Return a dictionary with umi count for the different umi categories for a given region
    '''
    
    D = {}
    
    # get region from file name
    region = FormatRegion(DataFile)
    infile = open(DataFile)
    Header = infile.readline().strip().split('\t')
    line = infile.readline().strip()
    if line != '':
        line = line.split()
        # get total parent umis
        ptu = int(line[Header.index('PTU')])
        # get total child umis
        ctu = int(line[Header.index('CTU')])
        # get numbers of children
        children = line[Header.index('CHILD_NUMS')]
        parents = line[Header.index('FREQ_PARENTS')]
        D[region] = {'PTU': ptu, 'CTU': ctu, 'children': children, 'parents': parents} 
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
        region = FormatRegion(filename)
        # extract umi counts
        d = ExtractUmiCounts(filename)
        umis = d[region]['PTU']
        D[region] = umis
    return D

def GetSampleCoverage(L):
    '''
    (list) -> dict
    
    :param L: A list of full paths to consensus files with umi count per interval (ie. files not merged and not empty)
    
    Returns a dictionary of interval coordinates with a list with mean and s.e.m. of coverage within the interval
    '''
    
    D = {}
    for filename in L:
        # extract region from filename 
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
        if ':' in i:
            i = i.split(':')
        elif '_' in i:
            i = i.split('_')
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
    contigs = [i.replace('chr', '') for i in Chromos]
    # place non-numeric contigs at beginining
    for i in range(len(contigs)):
        if contigs[i].isnumeric() == False:
            j = contigs.pop(i)
            contigs.insert(0, j)
    # remove non-numeric chromos from contigs and add to new list
    aside = []
    while contigs[0].isnumeric() == False:
        aside.append(contigs.pop(0))
    aside.sort()
    # convert contigs to int and sort
    contigs = sorted(list(map(lambda x:int(x), contigs)))    
    # add back non-numerical chromos
    contigs.extend(aside)
    # add back 'chr' to contigs
    contigs = list(map(lambda x: 'chr' + str(x), contigs))
       
    Positions = []
    for i in contigs:
        for j in Chromos[i]:
            Positions.append(i + ':' + str(j[0]) + '-' + str(j[1]))
    return Positions                
 
             
def PlotDataPerRegion(CoverageStats, DataFiles, **Options):
    '''
    (str, list, dict) -> None or str
    
    :param CoverageStats: yaml file with mean read depth per region in Stats directory
    :param Options: Optional parameters, accepted values:
                    'minval': Minimum value. Values lower are colored in red
                    'outputfile': Name of output SVG file
                    'datatype': Data to be ploted ['umis', coverage', 'children', ratio']
    
    Generates a plot of the datatype per region. Returns a base64 string
    of the svg image by default, unless the outputfile option is used
    '''
    
    # plot coverage by default
    if 'datatype' in Options:
        datatype = Options['datatype']
        if datatype not in ['coverage', 'umis', 'ratio', 'children']:
            datatype = 'coverage'
    else:
        # plot coverage by default
        datatype = 'coverage'
    
    # get extract data. data is a dict with cordinates as key region-data as value 
    if datatype == 'coverage':
        # get mean coverage per interval
        data = GetCoverageRegions(CoverageStats)
        Ytitle='Mean read depth per region'
        Title='Coverage'
    elif datatype in ['children', 'ratio', 'umis']:
        # extract umi counts for each region
        L = [ExtractUmiCounts(i) for i in DataFiles]
        data = {}
        for d in L:
            region = list(d.keys())[0]
            ptu, ctu = d[region]['PTU'], d[region]['CTU']
            if datatype == 'ratio':
                # compute child/parent umi ratios for each region        
                if ptu != 0:
                    data[region] = ctu/ptu
                Ytitle = 'Child:Parent Ratio'
                Title = "Children to Parent UMIs"
            elif datatype == 'umis':
                # get total parent umis for each region
                data[region] = ptu
                Ytitle = 'Total UMI counts'
                Title = "Parent + children UMIs"
            elif datatype == 'children':
                # get children umis for each region
                data[region] = ctu
                Ytitle = 'Children UMI counts'
                Title = "Children UMIs"
     
    # sort coordinates according to data values
    L = list(reversed(sorted([(data[i], i) for i in data])))
    Coordinates = [i[1] for i in L]
    dataVals = [data[i] for i in Coordinates]
        
    # create lists to store data and metadata
    low_data, high_data, all_data = [],[],[]
    if 'minval' in Options:
        minval = Options['minval']
        # set minval to 0 if not proper format
        try:
            float(minval) or int(minval)
        except:
            minval = 0
        for i in range(len(dataVals)):
            if dataVals[i] < minval:
                low_data.append({'value': (i, dataVals[i]), 'color':'red',
                                 'label':Coordinates[i]})
            else:
                high_data.append({'value': (i, dataVals[i]), 'color':'blue',
                                  'label':Coordinates[i]})
    else:
        for i in range(len(dataVals)):
            all_data.append({'value': (i, dataVals[i]), 'color':'blue',
                                  'label':Coordinates[i]})
    # define custom style 
    custom_style = Style(opacity='.3', opacity_hover='.9',
                         value_font_size=15,
                         value_colors=('black',),
                         no_data_font_size=40, no_data_text='No result found',
                         label_font_size=12,
                         background = 'white',
                         font_family='Arial, Verdana, sans-serif',
                         base_style=DefaultStyle)
    # set number of columns for legend
    if len(low_data) != 0:
        ncol=2
    else:
        ncol=1
        
    # use scatter plot
    xy_chart = pygal.XY(stroke=False, show_legend=True,
                        human_readable=True, fill=False,
                        show_x_labels=False,
                        spacing=10, max_scale=10,
                        dots_size=5, truncate_legend=-1,
                        legend_box_size=18,
                        dynamic_print_values=True, style=custom_style,
                        y_title=Ytitle,
                        x_title='Genomic intervals',
                        include_x_axis=True,
                        legend_at_bottom=True,
                        legend_at_bottom_columns=ncol,
                        title=Title,
                        explicit_size=False,
                        width=800,
                        height=400)
    
    if len(low_data) != 0:
        xy_chart.add({'title':'Low (< {0})'.format(minval), 'color':'red'}, low_data)
    if len(high_data) != 0:
        if len(low_data) == 0:
            xy_chart.add({'title':'Read depth'.format(minval), 'color':'blue'}, high_data)
        else:
            xy_chart.add({'title':'High'.format(minval), 'color':'blue'}, high_data)
    if len(all_data) != 0:
        xy_chart.add({'title':'Read depth', 'color':'blue'}, all_data)
    
    # save as svg by default if outputfile provided
    # and return a base64 string otherwise
    if 'outputfile' in Options:
        outputfile = Options['outputfile']
        extension = 'svg'
    else:
        extension = ''
    
    if extension == 'svg':
        xy_chart.render_to_file(outputfile + '.' + extension)
    else:
        return xy_chart.render_data_uri()

    
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
    FamSize = sorted(Data.keys())
        
    # make a sorted list of positions
    # collect all positions across families because some positions may be missing for some family size
    positions = []
    for i in Data:
        positions.extend(list(Data[i].keys()))
    positions = list(map(lambda x: str(x), sorted(list(map(lambda x: int(x), list(set(positions)))))))
    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot data  
    for i in range(len(FamSize)):
        # make a list of values for each position, including missing values
        yvals = [Data[FamSize[i]][j] if j in Data[FamSize[i]] else None for j in positions]
        # convert to numpy array
        yvals = np.array(yvals).astype(np.double)
        # create a mask so that line plots doesn't leave gap between missing values
        ymask = np.isfinite(yvals)
        # create an array with positions
        xvals = np.arange(len(positions))
        # pass the mask to x and y values when plotting
        ax.plot(xvals[ymask], yvals[ymask], color = Color[i], marker='', linewidth=2, linestyle='-', alpha = 1)
        
    # limit y axis
    YMax = []
    for i in Data:
        for j in Data[i]:
            YMax.append(Data[i][j])
    YMax = max(YMax)
    if YMax < 10:
        YMax = float(YMax + (YMax * 25 /100))
    else:
        YMax = float(YMax + (YMax * 15 /100))
    ax.set_ylim([0, YMax])    
    if YMax < 10:
        ystep = 2
    else:
        ystep = SetUpTicks(YMax)    
    ax.yaxis.set_ticks([i for i in np.arange(0, YMax, ystep)])

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
        
    # set up x axis. divide genomic interval in 3
    xtickspos = list(map(lambda x: math.floor(x), [i for i in np.arange(0, len(positions)+1, (len(positions)-1) / 3)]))
    xticks = [positions[i] for i in xtickspos]
    plt.xticks(xtickspos, xticks, ha = 'center', rotation = 0, fontsize = 12)
    
    # do not show y ticks
    plt.tick_params(axis='both', which='both', bottom=True, top=False,
                right=False, left=False, labelbottom=True, colors = 'black',
                labelsize = 12, direction = 'out')  
      
    # add legend
    legend_elements = []
    # loop over family sizes
    for i in range(len(FamSize)):
        legend_elements.append(Line2D([0], [0], marker='', label=str(FamSize[i]), linestyle='-', color = Color[i]))
    ax.legend(handles=legend_elements, frameon=False, bbox_to_anchor=(1.05, 0.95), loc='upper left', borderaxespad=0.)
    return ax


def PlotMeanFamSize(ConsFile, Color, Outputfile, W, H):
    '''
    (str, list, str, int, int) -> None
    
    :param ConsFile: Path to the consensus file
    :param Color: List with colors for plotting
    :param Outputfile: Name of the output figure file
    :param W: Figure width in inches
    :param H: Figure Height in inches
           
    Pre-condition: consensus file is not merged chrN:A-B.cons 
    '''
    
    # extract region from consensus file
    region = FormatRegion(ConsFile)
       
    # extract consensus depth for each family size -> {fam: {pos: meanfamSize}}
    Data = ExtractFamSize(ConsFile)
    
    # create figure object + ax
    figure = plt.figure()
    figure.set_size_inches(W, H)
    ax = CreateMeanFamAx(1, 1, 1, figure, Data, Color, 'Mean family size', region)
       
    plt.tight_layout()
    figure.savefig(Outputfile, bbox_inches = 'tight')
    plt.close()    

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
                    
 
    
def CreateNonRefFreqAx(Columns, Rows, Position, figure, Data, Color, fam_size, **Options):
    '''
    (int, int, int, figure_object, dict, str, int, dict) -> ax object
    
    :param columns: Number of columns
    :param rows: Number of rows
    :param position: Ax position in figure
    :param figure: Figure object opened for writing
    :param Data: Non-reference frequency at each position for a given family size
    :param Color: String color of the data
    :param Options: Accepted keys are:
                    'YLabel': Label of the Y axis
                    'XLabel': Label of the X axis
                    'YLimit': Y axis limit, in variant frequency (0-100) 
                    'title': Title of the subplot
                        
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
    if 'YLimit' in Options:
        try:
            YMax = float(Options['YLimit'])
        except:
            YMax = [Data[i] for i in Data]
            YMax = max(YMax)
    else:
        YMax = [Data[i] for i in Data]
        YMax = max(YMax)
    
    YMax = float(YMax + (YMax * 10 /100))
    #adjust YMax if 0, because of error bottom == top
    if YMax == 0:
        YMax = 1
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
    elif 30 < YMax <=50:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 10)])
    elif YMax > 50:
        ax.yaxis.set_ticks([i for i in np.arange(0, YMax, 20)])
    
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
    
    # set up x axis
    # divide genomic interval in 3
    xtickspos = list(map(lambda x: math.floor(x), [i for i in np.arange(0, len(pos)+1, (len(pos)-1) / 3)]))
    xticks = [pos[i] for i in xtickspos]
    plt.xticks(xtickspos, xticks, ha = 'center', rotation = 0, fontsize = 12)
    
    if 'XLabel' in Options:
        # show ticks
        plt.tick_params(axis='both', which='both', bottom=True, top=False,
                    right=False, left=False, labelbottom=True, colors = 'black',
                    labelsize = 12, direction = 'out')  
    else:
        # do not show ticks
        plt.tick_params(axis='both', which='both', bottom=True, top=False,
                    right=False, left=False, labelbottom=False, colors = 'black',
                    labelsize = 12, direction = 'out')  
    return ax


def PlotNonRefFreqData(ConsFile, Color, Outputfile, W, H, **Options):
    '''
    (str, list, str, int, int, dict) -> None
    
    :param ConsFile: Path to the consensus file
    :param Color: List with colors for plotting
    :param Outputfile: Name of the output figure file
    :param W: Figure width in inches
    :param H: Figure Height in inches
    :param Options: Accepted keys are:
                    'YLimit': Y axis limit, in variant frequency (0-100) 
                    'title': Title of the plot
                               
    Pre-condition: consensus file is not merged chrN:A-B.cons 
    '''
    
    if 'YLimit' in Options:
        try:
            YLimit = float(Options['YLimit'])
        except:
            YLimit = ''
    else:
        YLimit = ''
    
    # extract region from consensus file
    region = FormatRegion(ConsFile)
       
    # extract non-reference frequency for all family sizes in consensus file
    Data = ExtractNonRefFreq(ConsFile)
    
    # handle missing data among family sizes
    # make a list with all positions across all family sizes
    positions =[]
    for i in Data:
        positions.extend(list(Data[i].keys()))
    positions = list(set(positions))    
    # set missing values to 0
    for i in Data:
        for j in positions:
            if j not in Data[i]:
                Data[i][j] = 0
       
    # create figure
    figure = plt.figure()
    figure.set_size_inches(W,H)
    
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
            ax = CreateNonRefFreqAx(1, len(L), i+1, figure, L[i], Color[i], FamSize[i], YLimit=YLimit)
            # add title
            if 'title' in Options:
                title = Options['title']
                ax.set_title(title, size=12, loc='center', ha='center')
            # add legend
            if 'legend' in Options:
                legend_elements = []
                for i in range(len(L)):
                    legend_elements.append(Line2D([0], [0], marker='s', label=str(FamSize[i]), linestyle='None', color = Color[i]))
                ax.legend(handles=legend_elements, frameon=False, ncol=1, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)    
        elif i == len(L) // 2:
            ax = CreateNonRefFreqAx(1, len(L), i+1, figure, L[i], Color[i], FamSize[i], YLimit=YLimit)
            if 'ylabel' in Options:
                ylabel = Options['ylabel']
                # axes coordinates are 0,0 is bottom left and 1,1 is upper right
                ax.text(-0.1, 0.5, ylabel, horizontalalignment='right',
                        verticalalignment='center', rotation='vertical', transform=ax.transAxes, fontsize=12)
        elif i == len(L) - 1:
            ax = CreateNonRefFreqAx(1, len(L), i+1, figure, L[i], Color[i], FamSize[i], XLabel= region, YLimit=YLimit)
        else:
            ax = CreateNonRefFreqAx(1, len(L), i+1, figure, L[i], Color[i], FamSize[i], YLimit=YLimit)

    plt.tight_layout()
    figure.savefig(Outputfile, bbox_inches = 'tight')
    plt.close()

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
                    
  
def CreateConsDepthAx(Columns, Rows, Position, figure, Data, positions, Color, YLabel, **Options):
    '''
    (int, int, int, figure_object, list, list, list, str, dict) -> ax object
    
    :param columns: Number of columns
    :param rows: Number of rows
    :param position: Ax position in figure
    :param figure: Figure object opened for writing
    :param Data: Depth at each position for each family size
    :param positions: Sorted list of positions within genomic intervals
    :param Color: List of colors
    :param Options: Accepted keys are:
                    'XLabel': Label of the X axis
                    'legend': Add legend (True) or not (False)
                    'fam_size': List of family sizes in consensus file
                    'Colors': List of colors, parallel to fam_size 
                        
    Return a ax in figure
    '''
    
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # plot data
    for i in range(len(Data)):
        # add missing values, convert list to numpy array
        yvals = np.array([Data[i][j] if j in Data[i] else None for j in positions]).astype(np.double)
        # create mask so that line plots doesn't leave gaps between missing values
        ymask = np.isfinite(yvals)
        # create array with positions
        xvals = np.arange(len(positions))
        # pass the mask to x and y values when plotting
                
        #ax.plot([j for j in range(len(pos))], [Data[i][j] for j in pos], color = Color[i], marker='', linewidth=2, linestyle='-', alpha = 1)
    
    
        ax.plot(xvals[ymask], yvals[ymask], color = Color[i], marker='', linewidth=2, linestyle='-', alpha = 1)
    
    # limit y axis
    YMax = []
    for i in range(len(Data)):
        for j in Data[i]:
            YMax.append(Data[i][j])
    YMax = max(YMax)
    YMax = float(YMax + (YMax * 10 /100))
    ax.set_ylim([0, YMax])    
    step = SetUpTicks(YMax)    
    ax.yaxis.set_ticks([i for i in np.arange(0, YMax, step)])
     
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

    # set up x axis
    # divide genomic interval in 3
    xtickspos = list(map(lambda x: math.floor(x), [i for i in np.arange(0, len(positions)+1, (len(positions)-1) / 3)]))
    xticks = [str(positions[i]) for i in xtickspos]
    plt.xticks(xtickspos, xticks, ha = 'center', rotation = 0, fontsize = 12)
    
    if 'XLabel' in Options:
        # show ticks
        plt.tick_params(axis='both', which='both', bottom=True, top=False,
                    right=False, left=False, labelbottom=True, colors = 'black',
                    labelsize = 12, direction = 'out')  
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
            ax.legend(handles=legend_elements, frameon=False, bbox_to_anchor=(1.05, 0.95), loc='upper left', borderaxespad=0.)
       
    return ax


def PlotConsDepth(ConsFile, Color, Outputfile, W, H):
    '''
    (str, list, str, int, int) -> None
    
    :param ConsFile: Path to the consensus file
    :param Color: List with colors for plotting
    :param Outputfile: Name of the output figure file
    :param W: Figure width in inches
    :param H: Figure Height in inches       
    
    Pre-condition: consensus file is not merged chrN:A-B.cons and not empty
    '''
    
    # extract region from consensus file
    region = FormatRegion(ConsFile)
        
    # extract depth at each position for each family size
    Data = ExtractDepth(ConsFile)
     
    # make a sorted list of dicts {pos: depth} for each family size
    FamSize = sorted(Data.keys())
    L = []
    for i in FamSize:
        d = {}
        for pos in Data[i]:
            d[pos] = Data[i][pos]
        L.append(d)
    
    # make a sorted list of positions across all family size
    # because some positions my be missing for some family size
    positions = []
    for i in FamSize:
        positions.extend(list(Data[i].keys())) 
    positions = sorted(list(set(positions)))
    
    # create figure
    figure = plt.figure()
    figure.set_size_inches(W, H)
    
    # plot raw depth, family size = 0    
    ax1 = CreateConsDepthAx(1, 2, 1, figure, L[0:1], positions, Color[0:1], 'Raw depth', legend=True, fam_size=FamSize, colors=Color)
    ax2 = CreateConsDepthAx(1, 2, 2, figure, L[1:], positions, Color[1:], 'Consensus depth', XLabel=region)
    
    plt.tight_layout()
    figure.savefig(Outputfile, bbox_inches = 'tight')
    plt.close() 
 
def ComputeIntervalSize(Coordinates):
    '''
    (list) -> list
    
    :param Coordinates: A list of genomic regions in the format chrN:posA-posB
    
    Returns a list of genomic interval size in bp
    '''

    # get the interval size for each region in sorted list of coordinates
    Sizes = []
    for i in range(len(Coordinates)):
        interval = list(map(lambda x: float(x), Coordinates[i][Coordinates[i].index(':')+1:].split('-')))
        Sizes.append(interval[1] - interval[0])
    return Sizes      
    
   
def PlotParentsToChildrenCounts(DataFiles, Outputfile, W, H):
    '''
    (list, str, int, int) -> None
    
    :param DataFiles: List of data files generated after umi grouping 
    :param Outputfile: Name of the output figure file 
    :param Graph: Type of data to plot. Accepted values:
                  'ratio': children to parent umis ratio
                  'parents': total umi count
                  'children': children umi count 
    :param W: Figure width in inches
    :param H: Figure Height in inches
    
    Generates a plot with umi counts (children, parents or children to parents ratio)
    
    Pre-condition: Data files are not merged (datafile_chrN:A-B.csv) and not empty
    '''
    
    # extract umi counts for each region
    L = [ExtractUmiCounts(i) for i in DataFiles]
    
    Data = {}
    for d in L:
        region = list(d.keys())[0]
        ptu, ctu = d[region]['PTU'], d[region]['CTU']
        Data[region] = [ctu, ptu]
    
    # get a sorted list of positions
    Coordinates = SortPositions(list(Data.keys()))
    
    # get the interval size for each region in sorted ist of coordinates
    Sizes = ComputeIntervalSize(Coordinates)
        
    # map colors to interval size          
    cmap = plt.get_cmap('Reds', max(Sizes))
        
    # make parallel lists of children and parent counts
    CTU = [Data[i][0] for i in Coordinates]
    PTU = [Data[i][1] for i in Coordinates]
    
    # create figure
    figure = plt.figure()
    figure.set_size_inches(W, H)
    
    # add a plot coverage to figure (N row, N column, plot N)
    ax = figure.add_subplot(1, 1, 1)
    # plot ctu/ptu ratio for each region
    im = ax.scatter(PTU, CTU, edgecolor = 'lightgrey', clip_on=False, c = Sizes, cmap = cmap, marker='o', lw = 1, s = 160, alpha = 0.7)
        
    # limit y axis to maximum value
    YMax = max(CTU)
    # add 10% to max value
    YMax = YMax + (YMax * 10/100)
    ax.set_ylim([0, YMax])
    # set y ticks    
    step = SetUpTicks(YMax)
    ax.yaxis.set_ticks([i for i in np.arange(0, YMax, step)])
            
    # write label for y axis
    ax.set_ylabel('Number of children UMIs', color = 'black',  size = 14, ha = 'center')
    ax.set_xlabel('Number of parents UMIs', color = 'black',  size = 14, ha = 'center')
        
    # write title   
    ax.set_title('Interval Size vs. PTU and CTU', size = 14)
    
    # limit x axis to maximum value
    XMax = max(PTU)
    # add 10% to max value
    XMax = XMax + (XMax * 10/100)
    ax.set_xlim([0 - (XMax * 10/100)/2, XMax])
    # set y ticks    
    step = SetUpTicks(XMax)
    ax.xaxis.set_ticks([i for i in np.arange(0, XMax, step)])
    
    # add space between axis and tick labels
    ax.yaxis.labelpad = 18
    ax.xaxis.labelpad = 18
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)  
    
    # offset the x axis
    for loc, spine in ax.spines.items():
        spine.set_position(('outward', 10))
        spine.set_smart_bounds(True)
       
    # do not show ticks
    plt.tick_params(axis='both', which='both', bottom=True, top=False,
                right=False, left=False, labelbottom=True, colors = 'black',
                labelsize = 12, direction = 'out')  
        
    # add color bar
    cb = figure.colorbar(im, ax=ax)
    cb.set_label('Interval size', size=14, labelpad=18)
        
    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.4, linewidth = 0.4)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)
    
    figure.savefig(Outputfile, bbox_inches = 'tight')
    plt.close()

def PlotParentFreq(DataFiles, Outputfile, W, H):
    '''
    (list, str, int, int) -> None
    
    :param DataFiles: List of .csv data files generated after umi grouping
    :param Outputfile: Name of the output figure file 
    :param W: Figure width in inches
    :param H: Figure Height in inches 
    
    Generates a plot with parent frequencies vs children count for all regions
    
    Pre-condition: consensus and data files are not merged (chrN:A-B.cons and chrN:A-B.csv)
    '''
    
    # extract umi counts for each region
    L = [ExtractUmiCounts(i) for i in DataFiles]
     
    Data = {}
    for d in L:
        region = list(d.keys())[0]
        # get total count
        ptu = d[region]['PTU']
        # get the counts of children and parents
        children, parents = d[region]['children'].split(','), d[region]['parents'].split(',')
        children = list(map(lambda x: int(x.strip()), children))
        # compute parent frequencies
        parents = list(map(lambda x: int(x.strip()) / ptu, parents))
        # map children to parents   
        k ={children[i]: parents[i] for i in range(len(children))}
        Data[region] = k
    
    # get a sorted list of positions
    Coordinates = SortPositions(list(Data.keys()))
    
    # create a cmap and list of colors using divergening colors PiYG
    cmap = plt.get_cmap('PiYG')
    colors = [cmap(i) for i in np.linspace(0, 1, len(Coordinates))]
        
    # create figure
    figure = plt.figure()
    figure.set_size_inches(W, H)
    
    # add a plot coverage to figure (N row, N column, plot N)
    ax = figure.add_subplot(1, 1, 1)
    # loop over sorted regions
    for i in range(len(Coordinates)):
        # plot parent frequencies vs sorted number of children
        ax.scatter(sorted(Data[Coordinates[i]].keys()), [Data[Coordinates[i]][j] for j in sorted(Data[Coordinates[i]].keys())], edgecolor = colors[i], facecolor = colors[i], marker='o', lw = 1, s = 70, alpha = 0.3, clip_on=False)
    
    # limit y axis to maximum value
    YMax = []
    for i in Data:
        for j in Data[i]:
            YMax.append(Data[i][j])
    YMax = max(YMax)
    # add 10% to max value
    YMax = YMax + (YMax * 10/100)
    ax.set_ylim([0, YMax])
    
    # set Y axis ticks
    step = round(YMax/10, 2)    
    ax.yaxis.set_ticks([i for i in np.arange(0, YMax, step)])
    
    # write label for y axis
    ax.set_ylabel('Parent UMI frequency', color = 'black',  size = 14, ha = 'center')
    ax.set_xlabel('Number of children UMIs', color = 'black',  size = 14, ha = 'center')
        
    # write title   
    ax.set_title('Parent Frequency vs Children UMIs', size = 14)
    
    # write label for x axis
    children = []
    for i in range(len(Coordinates)):
        children.extend(list(Data[Coordinates[i]].keys()))
    children = sorted(list(set(children)))
    minchildren, maxchildren = children[0], children[-1]
    xstep = SetUpTicks(maxchildren)
    children = [i for i in range(0, maxchildren + 1, xstep)]
    xPos = [i for i in children]
    plt.xticks(xPos, list(map(lambda x: str(x), children)), ha = 'center', rotation = 0, fontsize = 9)
               
    # add space between axis and tick labels
    ax.yaxis.labelpad = 18
    ax.xaxis.labelpad = 18
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)  
       
    # offset the x axis
    for loc, spine in ax.spines.items():
        spine.set_position(('outward', 5))
        spine.set_smart_bounds(True)
    
    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.4, linewidth = 0.4)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)
    
    # do not show ticks
    plt.tick_params(axis='both', which='both', bottom=True, top=False,
                right=False, left=False, labelbottom=True, colors = 'black',
                labelsize = 12, direction = 'out')  
    figure.savefig(Outputfile, bbox_inches = 'tight')
    plt.close()

def ParentToChildren(UmiFile):
    '''
    (str) -> list
    
    :param UmiFile: Path to file with umi sequences and child to parent relationships
    
    Return a list of tuples with (parent, child. hamming_distance) umi sequences
    
    Precondition: the umi file is not merged
    '''
    
    # load umi file with umi sequences and chil to parent relations
    infile = open(UmiFile)
    data = json.load(infile)
    infile.close()
    
    # create a dict {parent: [children]}
    D = {}
    for i in data:
        parent = data[i]['parent']
        if parent not in D:
            D[parent] = []
        D[parent].append(i)
    
    # make a list of (parent, child, distance)
    L = []
    for i in D:
        for j in D[i]:
            # compute distance between parent and child umi
            d = edit_distance(i, j)
            L.append((i, j, d))
    return L



def LoadDataToNetwork(G, L):
    '''
    (networkx.classes.digraph.DiGraph, list) -> networkx.classes.digraph.DiGraph
    
    :param G: Directed graph networkx object 
    :param L: List of nodes and edges
    
    Return the network with nodes and edges
    '''
    
    # loop over 
    for i in L:
        # do not add self-edges
        if i[-1] == 0:
            # distance is 0, add node but no edge 
            G.add_node(i[0])
        else:
            # add edge connection parent to child umi
            G.add_edge(i[0], i[1])
    return G
    

def BuildNetwork(UmiFile):
    '''
    (str) -> networkx object
    
    :param UmiFile:  Path to file with umi sequences and child to parent relationships 
    
    Return a networkx object
    '''
    
    # make a list of (parent, child) umi sequences
    L = ParentToChildren(UmiFile)
    # build directed network
    G = nx.DiGraph()
    G = LoadDataToNetwork(G, L)
    return G


def CreateNetworkAx(Columns, Rows, Position, figure, UmiFile):
    '''
    (int, int, int, figure_object, str) -> ax object
    
    :param columns: Number of columns
    :param rows: Number of rows
    :param position: Ax position in figure
    :param figure: Figure object opened for writing
    :param UmiFile: Path to json file with umi parent: children relationships
    
    Return a ax in figure
    '''
   
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    
    # build network
    G = BuildNetwork(UmiFile)
    
    # convert the graph to a dict
    # method differs based on networkx version
    try:
        d = dict(G.adjacency())
    except:
        d = dict(G.adjacency_iter())
    
    # get degree for all nodes
    degree = dict(G.degree())
    
    # make lists of parents, children and nodes without edges 
    singles, parents, children = [], [], []
    for i in d:
        # check if node has edges
        if len(d[i]) == 0:
            # no edge
            singles.append(i)
        else:
            # add parent and children
            parents.append(i)
            for j in d[i]:
                children.append(j)
    
    # make list of all nodes with edges
    AllNodes = list(set(parents + children))
        
    # write title
    ax.set_title('Parent-Child UMI network', size = 14)
       
    # add space between axis and tick labels
    ax.yaxis.labelpad = 18
    ax.xaxis.labelpad = 18
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(False)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)  
       
    # do not show ticks
    plt.tick_params(axis='both', which='both', bottom=False, top=False,
                    right=False, left=False, labelleft=False, labelbottom=False,
                    labelright=False, colors = 'black', labelsize = 12, direction = 'out')  
    
    # set up same network layout for all drawings
    Pos = nx.spring_layout(G)
    # draw edges    
    nx.draw_networkx_edges(G, pos=Pos, width=0.7, edge_color='grey', style='solid',
                           alpha=0.4, ax=ax, arrows=False, node_size=5,
                           nodelist=AllNodes, node_shape='o')
    # draw all nodes, color according to degree
    nodelist = sorted(degree.keys())
    node_color = [degree[i] for i in nodelist]
    
    # limit color map to the highest degree
    cmap=plt.get_cmap('winter', max(node_color)+1)
    
    nodes = nx.draw_networkx_nodes(G, pos=Pos, with_labels=False, node_size=5,
                                   node_color=node_color, node_shape='o', alpha=0.3,
                                   linewidths=0, edgecolors='grey', ax=None,
                                   nodelist=nodelist, cmap=cmap)
    nodes.set_clim(min(node_color), max(node_color)+1) 

    # add discrete color bar for node degree
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=0.05)
    
    
    # make a list of degree 
    degree_vals = sorted(list(set(node_color)))
    cb = figure.colorbar(nodes, cax=cax, orientation = 'horizontal', ticks=[i for i in range(min(degree_vals), max(degree_vals)+2)], use_gridspec=False)
    
    
    #    # write x ticks
    if len(list(set(degree_vals))) < 10:
        step = 1
    elif 10 <= len(list(set(node_color))) < 20:
        step = 2
    elif 20 < len(list(set(node_color))) < 60:
        step = 10
    elif len(list(set(node_color))) >= 60:
        step = 20
    
    xtickslabels = [str(i) if i % step == 0 else '' for i in range(min(degree_vals), max(degree_vals)+2)]
    
    cb.ax.set_xticklabels(xtickslabels)
    
    
    
    #cb = figure.colorbar(nodes, cax=cax, orientation = 'horizontal', ticks=[i for i in range(min(node_color), max(node_color)+2)], use_gridspec=False)
    
    
    
    
#    # write x ticks
#    if len(list(set(node_color))) < 10:
#        step = 1
#    elif 10 <= len(list(set(node_color))) < 20:
#        step = 2
#    elif 20 < len(list(set(node_color))) < 60:
#        step = 10
#    elif len(list(set(node_color))) >= 60:
#        step = 20
#    cb.ax.set_xticklabels([str(i) for i in range(min(node_color), max(node_color)+2, step)])
    cb.set_label('Node degree', size=14, ha='center', color='black', labelpad=18)
            
    return ax


def CreateDegreeAx(Columns, Rows, Position, figure, UmiFile):
    '''
    (int, int, int, figure_object, str) -> ax object
    
    :param columns: Number of columns
    :param rows: Number of rows
    :param position: Ax position in figure
    :param figure: Figure object opened for writing
    :param UmiFile: Path to json file with umi parent: children relationships
    
    Return a ax in figure
    '''

    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    
    # build network
    G = BuildNetwork(UmiFile)
    Degree = dict(G.degree())
    
    # make a list of node degree
    degree_sequence = sorted([d for d in Degree.values()], reverse=True)
    # count nodes with a given degree
    degree_count = collections.Counter(degree_sequence)
    # get node degree 
    degree = sorted(degree_count.keys())
    # compute minimum and maximum degree
    mindegree, maxdegree = degree[0], degree[-1]
    # make parallel lists of degree and count sorted on degree
    # include all degree values between min and max degree
    degree = [i for i in range(mindegree, maxdegree + 1)]
    count = []
    for i in degree:
        if i in degree_count:
            count.append(degree_count[i])
        else:
            count.append(0)
    
    # plot network degree
    ax.bar(degree, count, width=0.4, color='#eaccff', edgecolor=['grey'] * len(degree), linewidth=0.7)
               
    # limit y axis and sset y axis ticks
    YMax = max(count)
    YMax = float(YMax + (YMax * 10 /100))
    ax.set_ylim([0, YMax])    
    step = SetUpTicks(YMax)
    ax.yaxis.set_ticks([i for i in np.arange(0, YMax, step)])
    
    # set up y axis label and grid
    ax.set_ylabel('Count', color = 'black',  size = 14, ha = 'center')
    ax.set_xlabel('Node degree', color = 'black',  size = 14, ha = 'center')
    
    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.4, linewidth = 0.4)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)
    
    # set x ticks
    XMax = float(maxdegree + (maxdegree * 10/100))
    step = SetUpTicks(XMax)
    ax.xaxis.set_ticks([int(i) for i in np.arange(0, XMax, step)])
    
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
    return ax
    
    
def PlotNetworkDegree(UmiFile, Outputfile, W, H):
    '''
    (str, str, int, int) -> None
    
    :param UmiFile: Path to json file with umi parent: children relationships
    :param Outputfile: Name of output figure file
    :param W: Figure width in inches
    :param H: Figure Height in inches
    
    Plot an histogram of network degree and the entire network of umis parent-children relationships    
    '''

    # create figure
    figure = plt.figure()
    figure.set_size_inches(W, H)
    
    # plot network degree
    ax1 = CreateDegreeAx(2, 1, 1, figure, UmiFile)
    # plot the network
    ax2 = CreateNetworkAx(2, 1, 2, figure, UmiFile)
    
    # save figure    
    plt.tight_layout()
    figure.savefig(Outputfile, bbox_inches = 'tight')
    plt.close()

def GetUmiCountFromPreprocessing(Datafile):
    '''
    (str) -> list
    
    :param Datafile: Path to file with UMI counts generated during pre-preprocessing
    
    Returns a list of umi occurence: count of each umi sequeces
    '''
    
    # get the umi count from the data file
    infile = open(Datafile)
    infile.readline()
    L = []
    for line in infile:
        if line.rstrip() != '':
            line = line.rstrip().split()
            count = int(line[1])
            L.append(count)
    infile.close()
    
    return L
    

def GetFamilyReadDepth(UmiFile):
    '''
    (str) -> dict
    
    :param UmiFile: Path to json file with umi-parent relationships and family count after grouping
    
    Returns a dictionary of read depth for each family and position
    '''
    
    infile = open(UmiFile)
    umis = json.load(infile)
    infile.close()

    # group umis by family
    D ={}
    for i in umis:
        parent = umis[i]['parent']
        if parent not in D:
            D[parent] = {}
        for j in umis[i]['positions']:
            # umi count from grouping is already the count of all umis from a same family at a given position
            # grab the count for the first umi of the family, no need to record count of other family members
            D[parent][j] = umis[i]['positions'][j]
    return D


def GetIndividualUmiInfo(UmiFile):
    '''
    (str) -> tuple
    
    :param UmiFile: Path to json file with umi-parent relationships and family
    count before grouping generated by get_umi_families 
    
    Returns a tuple with lists of umi count (=nb reads for given umi and position)
    for all umis and parent umis only
    '''

    # extract dict {parent: {umi:{pos:count}}}
    infile = open(UmiFile)
    umis = json.load(infile)
    infile.close()

    # get the umi occurence at each position for parents only and parents + children
    parent_umis, all_umis  = [], []
    
    for parent in umis:
        for umi in umis[parent]:
            # compute hamming distance
            distance = edit_distance(parent, umi)
            # get count at each position
            count = [umis[parent][umi][pos] for pos in umis[parent][umi]]
            if distance == 0:
                # umi is parent
                parent_umis.extend(count)
                all_umis.extend(count)
            else:
                # umi is children
                all_umis.extend(count)
    
    return all_umis, parent_umis

    
def GetUmiFamilySizeFromGrouping(UmiFile):
    '''
    (str) -> dict
    
    :param UmiFile: Path to json file with umi-parent relationships and family count after grouping
    
    Returns a dictionary of umi occurence: counts
    '''

    infile = open(UmiFile)
    umis = json.load(infile)
    infile.close()

    # group umis by family
    D ={}
    for i in umis:
        parent = umis[i]['parent']
        if parent not in D:
            D[parent] = {}
        for j in umis[i]['positions']:
            # grab all the umi sequences for a given family
            if j in D[parent]:
                D[parent][j].append(i)
            else:
                D[parent][j] = [i]
    
    # get the family size for each family and position
    for parent in D:
        for pos in D[parent]:
            D[parent][pos] = len(list(set(D[parent][pos])))
    return D


def PlotFamSizeReadDepth(UmiFile, Outputfile):
    '''
    (str, str, int, int) -> dict
    
    :param UmiFile: Path to json file with umi-parent relationships and family count after grouping
    :param Outputfile: Name of output figure file
            
    Plot a marginal plot of UMI family size and read depth    
    '''
    
    
    # set seaborn white style
    sns.set(style="white", color_codes=True)
    # restore tick marks
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['ytick.left'] = True
       
    # get the size of each family
    FamSize = GetUmiFamilySizeFromGrouping(UmiFile)
    # get the frequency distribution of read depth for each family
    ReadDepth =  GetFamilyReadDepth(UmiFile)
 
    # make parallel list of family size and read depth
    r, s = {}, {}
    i = 0
    assert FamSize.keys() == ReadDepth.keys()
    for parent in FamSize:
        for pos in FamSize[parent]:
            r[i] = ReadDepth[parent][pos]
            s[i] = FamSize[parent][pos]
            i += 1
    
    # create data frames
    S = pd.DataFrame(list(s.values()), columns=['size'])
    R = pd.DataFrame(list(r.values()), columns=['depth'])
    # join dataframes
    df = S.join(R)    
    
    # seaborn returns an ax object
    ax = sns.jointplot(x='depth', y='size', data=df, kind='scatter', color="pink",
                  space=0, ratio=3, marginal_kws={'bins':15, 'color':'#d9b3ff'}, annot_kws=dict(stat='r'),
                  s=35, linewidth=1).set_axis_labels('Read depth within group', 'Unique UMI seqs per group', size=14, color='black', ha='center')
    
    # use different color for y axis histogram and add edges to bars
    plt.setp(ax.ax_marg_y.patches, color="#80d4ff", lw=1)
    plt.setp(ax.ax_marg_x.patches, color='#d9b3ff', lw=1)
    
    # axes can be accessed using ax_marg_x/y and ax_joint then passing matplotlib axes methods
    # limit y axis and set y ticks
    YMax = max(df['size'])
    YMax = float(YMax + (YMax * 10 /100))
    ax.ax_joint.set_ylim([0, YMax])    
    step = SetUpTicks(YMax)
    ax.ax_joint.yaxis.set_ticks([i for i in np.arange(0, YMax + 1, step)])
    
    # limit x axis and set x ticks
    XMax = max(df['depth'])
    XMax = float(XMax + (XMax * 10 /100))
    ax.ax_joint.set_xlim([0, XMax])    
    step = SetUpTicks(XMax)
    ax.ax_joint.xaxis.set_ticks([i for i in np.arange(0, XMax + 1, step)])
    
    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.ax_joint.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.4, linewidth = 0.4)  
    # hide these grids behind plot objects
    ax.ax_joint.set_axisbelow(True)
    # add space between axis and tick labels
    ax.ax_joint.yaxis.labelpad, ax.ax_joint.xaxis.labelpad = 18, 18
    
    # do not show lines around figure  
    ax.ax_joint.spines["top"].set_visible(False)    
    ax.ax_marg_x.spines['bottom'].set_visible(False)
    ax.ax_joint.spines["bottom"].set_visible(True)    
    ax.ax_joint.spines["right"].set_visible(False)    
    ax.ax_marg_x.spines['right'].set_visible(False)
    ax.ax_joint.spines["left"].set_visible(False)  
    ax.ax_marg_y.spines["left"].set_visible(False)  

    # do not show ticks
    ax.ax_joint.tick_params(axis='both', which='both', bottom=True, top=False,
                right=False, left=True, labelleft=True, labeltop=False, labelbottom=True, colors = 'black',
                labelsize = 12, direction = 'out')  
    ax.ax_marg_x.tick_params(axis='both', which='both', bottom=False, top=False,
                right=False, left=False, labelleft=False, labeltop=False, labelbottom=False, colors = 'black',
                labelsize = 12, direction = 'out')  
    ax.ax_marg_y.tick_params(axis='both', which='both', bottom=False, top=False,
                right=False, left=True, labelleft=True, labeltop=False, labelbottom=False, colors = 'black',
                labelsize = 12, direction = 'out')

    plt.savefig(Outputfile, bbox_inches = 'tight')
    plt.close()

def PlotUMiFrequency(L, Outputfile, Title, overlapping, W, H):
    '''
    (list, str, str, str, bool, int, int) -> None
        
    :param L: List with umi counts
    :param Outputfile: Name of output figure file
    :param YLabel: Label of the Y axis
    :param XLabel: Label of the X axis
    :param overlapping: True or False for 2 distributions. 
                        If True, L is a 2-item list, each being a list of counts
    :param W: Figure width in inches
    :param H: Figure Height in inches        
            
    Plot an histogram of UMI occurence or density plots UMI occurence for multiple distributions
    '''
    
    # create figure
    figure = plt.figure()
    figure.set_size_inches(W, H)
        
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(1, 1, 1)
    
    # check if plotting a single of overlapping histograms
    if overlapping == False:
        # plot data using 20 bins
        counts, bins, patches = ax.hist(L, bins=20, facecolor='pink', lw=1, edgecolor='lightgrey', align='mid')
        # get the maximum values for x and y
        XMax, YMax = max(L), max(counts)
        # limit y axis and set ticks
        YMax = YMax + (10/100 * YMax)
        ax.set_ylim([0, YMax + 1])    
        step = SetUpTicks(YMax)
        ax.set_yticks([i for i in np.arange(0, YMax + 1, step)])    
        XLabel, YLabel = 'UMI occurence', 'Counts'
    
    elif overlapping == True:
        # plot density using seaborn
        ax = sns.kdeplot(L[0], color = '#ff66ff', shade=True, alpha=0.35, ax=ax, lw=2)
        ax = sns.kdeplot(L[1], color = '#00cccc', shade=True, alpha=0.35, ax=ax, lw=2)
        XMax = max(L[0] + L[1])
        # start x axis at 0
        ax.set_xbound(lower=0)
        XLabel, YLabel = 'UMI occurence', 'Density'
        
    if overlapping == False:
        # limit x axis and set x ticks
        ax.set_xlim([0, XMax + 1])    
    
    # customize x ticks     
    step = SetUpTicks(XMax)
    ax.set_xticks([i for i in np.arange(0, XMax + step, step)])
      
    # add title        
    ax.set_title(Title, size = 14)
        
    # set up y axis label and grid
    ax.set_ylabel(YLabel, color = 'black',  size = 14, ha = 'center')
    ax.set_xlabel(XLabel, color = 'black',  size = 14, ha = 'center')
    
    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.4, linewidth = 0.4)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)
    
    # add space between axis and tick labels
    ax.yaxis.labelpad, ax.xaxis.labelpad = 18, 18
    
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
    if overlapping==True:
        legend_elements = []
        legend_elements.append(Patch(facecolor='#ff66ff', edgecolor= '#ff66ff', label='parents + children', alpha=0.35))
        legend_elements.append(Patch(facecolor='#00cccc', edgecolor= '#00cccc', label='parents', alpha=0.35))
        ax.legend(handles=legend_elements, frameon=False, ncol=1, loc='best', prop={'size': 12})
     
    # save figure to file    
    figure.savefig(Outputfile, bbox_inches = 'tight')
    plt.close()

def CreateAxReadDepth(columns, rows, position, figure, data, YLabel):
    
    '''
    (int, int, int, figure_object, list, str) -> ax_object
    
    :param columns: Number of columns
    :param rows: Number of rows
    :param position: Ax position in figure
    :param figure: Figure object opened for writing
    :param data: Values to be plotted
    :param YLabel: Y axis label
        
    Return a ax object in figure
    '''
    
    # create an ax instance in figure
    ax = figure.add_subplot(rows, columns, position)
    # plot density distribution of read depth
    n1, b1, p1 = ax.hist(data[0], bins=15, align='mid', color = '#00cccc', alpha=0.5)
    n2, b2, p2 = ax.hist(data[1], bins=b1, align='mid', color = '#ff66ff', alpha=0.5)
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)  
    # do not show ticks
    plt.tick_params(axis='both', which='both', bottom=True, top=False,
                    right=False, left=False, labelbottom=True, colors = 'black',
                    labelsize = 12, direction = 'out')  
    
    # set x ticks    
    XMax = max(data[0] + data[1])
    step = SetUpTicks(XMax)
    ax.set_xticks([i for i in np.arange(0, XMax + step, step)])
    ax.set_xlabel('Read depth', color = 'black',  size = 14, ha = 'center')    
    
    YMax = max(n1 + n2)
    step = SetUpTicks(YMax)
    ax.set_yticks([i for i in np.arange(0, YMax + 1, step)])
    # set up y axis label and grid
    ax.set_ylabel(YLabel, color = 'black',  size = 14, ha = 'center')
    
    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.4, linewidth = 0.4)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)
    
    # add space between axis and tick labels
    ax.yaxis.labelpad, ax.xaxis.labelpad = 18, 18
    
    # add legend
    legend_elements = []
    legend_elements.append(Patch(facecolor='#00cccc', edgecolor= '#00cccc', label='highest', alpha=1))
    legend_elements.append(Patch(facecolor='#ff66ff', edgecolor='#ff66ff', label='others', alpha=1))
    ax.legend(handles=legend_elements, frameon=False, ncol=2, prop={'size': 12}, 
              bbox_to_anchor=(0.4, 1.1), loc='upper left', borderaxespad=0.)
    return ax
    
    
def PlotReadDepth(UmiFile, Outputfile, W, H):
    '''
    (str, str, int, int) -> None
    
    
    :param UmiFile: Path to json file with umi-parent relationships and family count after grouping
    :param Outputfile: Path to the output figure file
    :param W: Figure width in inches
    :param H: Figure Height in inches
        
    Plot the distribution of read depth at positions of highest and lower umi abundance
    '''
    
    # get read depth for each umi family and position {parent: position: read_depth}
    All = GetFamilyReadDepth(UmiFile)
       
    region = os.path.basename(UmiFile)
    region = region[:-5]
    if '_' in region:
        region = region.replace('_', ':')
    
    # separately record positions with the highest read depth and other positions
    MostFrequent, Others = {}, {}
    for parent in All:
        L = [(int(i.split(':')[1]), All[parent][i]) for i in All[parent]] 
        L.sort()
        # identify pos with most abundant umi within given family <-- (pos, read_depth)
        most_abundant = most_frequent(L)
        # record only positions with the highest read depth within group 
        MostFrequent[parent] = {}
        MostFrequent[parent][most_abundant[0]] = most_abundant[1]
        # record all other positions
        for i in L:
            if (i[0], i[1]) != most_abundant:
                if parent not in Others:
                    Others[parent] = {}
                Others[parent][i[0]] = i[1]
    
    # create figure
    figure = plt.figure()
    figure.set_size_inches(W, H)
        
    # get read depth for all positions and umi groups
    data1 = [list(All[i].values())[0] for i in All]
    # get read depth for positions with highest umi abundance
    data2 = [list(MostFrequent[i].values())[0] for i in MostFrequent]
    # get read depth for positions with lower umi abundance
    data3 = [list(Others[i].values())[0] for i in Others]
    
    # plot histograms of read depth and Umi family count
    ax = CreateAxReadDepth(1, 1, 1, figure, [data2, data3], 'Umi families')    
           
    # save figure to file  
    figure.savefig(Outputfile, bbox_inches = 'tight')
    plt.close()


def PlotIncorrectReads(ReadInfo, Outputfile, datatype, W, H):
    '''
    (str, str, str, int, int) -> None

    :param ReadInfo: Path to json with info about reads generated during pre-processing or grouping
    :param Outputfile: Path to the output figure file
    :param datatype: Step recording reads. Accepted values:
                     preprocessing: record incorrect/correct reads
                     mapping: record unmapped/mapped reads
    :param W: Figure width in inches
    :param H: Figure Height in inches     
        
    Generate a donut graph with proportions of correct/incorrect reads found
    during pre-processing or mapped/unmapped reads found during grouping
    '''

    infile = open(ReadInfo)
    data = json.load(infile)
    infile.close()

    if datatype == 'preprocessing':
        total, good, bad = data['Total'], data['Correct'], data['Incorrect']
        good_name, bad_name = 'correct', 'incorrect'
        Title = 'Pre-processed reads'
        size = [good/total * 100, bad/total * 100]
        s_good, s_bad = format(good, ','), format(bad, ',')
        colors = ['#d9b3ff', '#0073e6']
        
        # sort data 
        data = sorted(list(zip(size, [s_good, s_bad], [good_name, bad_name], colors)), key=lambda x:x[0])
        # get back counts, size, names and colors 
        size = [i[0] for i in data]
        counts = [i[1] for i in data]
        names = [i[2] for i in data]
        colors = [i[3] for i in data]
        # use MathText to highlight substring in bold
        labels = []
        for i in range(len(names)):
            labels.append("{0}\n".format(names[i]) + r"$\bf{" + str(counts[i]) + "}$")
    elif datatype == 'mapping':
        region = list(data.keys())[0]
        mapped, discarded = data[region]['mapped'], data[region]['unmapped'] + data[region]['secondary'] + data[region]['supplementary']
        total = mapped + discarded
        mapped_name, discarded_name = 'mapped', 'discarded'
        Title = 'Filtered reads'
        size = [mapped/total * 100, discarded/total * 100]
        s_mapped, s_discarded = format(mapped, ','), format(discarded, ',')
        colors = ['#ffe680','#990000']
        # sort data
        data = sorted(list(zip(size, [s_mapped, s_discarded], [mapped_name, discarded_name], colors)), key=lambda x:x[0])
        # get back counts and size
        size = [i[0] for i in data]
        counts = [i[1] for i in data]
        names = [i[2] for i in data]
        colors = [i[3] for i in data]
        # use MathText to highlight substring in bold
        labels = []
        for i in range(len(names)):
            labels.append("{0}\n".format(names[i]) + r"$\bf{" + str(counts[i]) + "}$")
           
    # create figure
    figure = plt.figure()
    figure.set_size_inches(W, H)
    
    # create ax instance
    ax = figure.add_subplot(1, 1, 1)
    
    # Create a circle for the center of the plot with radius 0.7
    my_circle=plt.Circle( (0,0), 0.7, color='white')
    
    # plot data as pie chart
    if datatype == 'preprocessing':
        lw = 7
    elif datatype == 'mapping':
        lw = 2
    ax.pie(size, labels=labels, colors=colors,
           textprops={'fontsize':22, 'fontweight':'normal'},
           wedgeprops = { 'linewidth' : lw, 'edgecolor' : 'white' })
    # add circle in the center of the pie to create a donut
    p=plt.gcf()
    p.gca().add_artist(my_circle)
    
    # add title
    ax.set_title(Title, size=18, loc='center', ha='center')
       
    # Equal aspect ratio ensures that pie is drawn as a circle
    ax.axis('equal')  

    plt.tight_layout()
    figure.savefig(Outputfile, bbox_inches = 'tight')
    plt.close()