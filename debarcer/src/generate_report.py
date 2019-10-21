# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 12:05:32 2019

@author: rjovelin
"""

import time
import mistune
import os
import scipy.ndimage
from itertools import zip_longest
import base64
import matplotlib.pyplot as plt
import pygal
from src.generate_plots import PlotDataPerRegion


def ResizeFifure(filename, scaling_factor):
    '''
    (str, float) -> (int, int)
    
    :param filename: Path to figure file
    :param scaling_factor: The factor applied to resize figure 
    
    Return new file size with same proportions as a tuple of height and width
    '''
    # extract the original figure size
    #height, width, channels = matplotlib.pyplot.imread(filename).shape
    height, width, channels = scipy.ndimage.imread(filename).shape
    
    # resize while keeping proportions between height and width
    height, width = list(map(lambda x: x * scaling_factor, [height, width]))
    return height, width        


def EncodeImage(filename):
    '''
    (str)- > str
    
    Return a string representing the encoding of image filename in base64
    '''
    
    # encode image in base64
    with open(filename, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return encoded_string
    

def GetExpectedFigure(FigDir, extension, expected_name):
    '''
    (str, str, str) --> str
    
    :param FigDir: Directory Figures in which figures are located
    :param extension: Extension of the figure files
    :param expected_name: Name of the figure file
    
    Return the full path of the figure expected_name + extension in directory Figures
    or the empty string if file doesn't exist
    '''
    
    figpath = os.path.join(FigDir, expected_name + extension)
    if os.path.isfile(figpath) == True:
        return figpath
    else:
        return ''
     
    
def ListAllRegions(directory, extension):
    '''
    (str, str) -> list
    
    :param directory:
    :param extension: Figure file extension (eg. png)
    
    Returns a list of regions extracted from all figure files located in Figures
    subfolder of directory
    '''

    # make a list if figure files
    FigDir = os.path.join(directory, 'Figures')
    Figures = [i for i in os.listdir(FigDir) if extension in i and os.path.isfile(os.path.join(FigDir, i))]    
    # make a list of regions
    regions = []
    for i in Figures:
        if '_chr' in i:
            region = i[i.index('_chr')+1:i.rindex('.')]
            regions.append(region)
    regions = list(set(regions))
    return regions


def ListExpectedFigures(directory, extension):
    '''
    (str, str) -> dict
    
    :param directory:
    :param extension: Figure file extension (eg. png)
    
    Returns a dictionary with paths to all expected figure files
    Path is the empty string if the expected figure is absent from subfolder Figures
    '''
    
    # make a list of figue files
    FigDir = os.path.join(directory, 'Figures')
          
    # create a dict to store path to each expected figure (can be '')
    D = {}
    
    # map expected figures to figure names for aggregate figures across regions
    N = ['reads', 'preprocessing', 'coverage', 'ratio', 'total', 'children',
         'interval', 'freq']
    L = ['Proportion_correct_reads.', 'UMI_occurence_preprocessing.', 'Coverage_Umi_Count.',
         'Child_Parent_Umis_Ratio.', 'Total_Umis.', 'Children_Umis.', 'PTU_vs_CTU.', 'Children_vs_ParentFreq.']
    
    for i in range(len(N)):
        D[N[i]] = GetExpectedFigure(FigDir, extension, L[i])
    
    # list all regions from figure files
    regions = ListAllRegions(directory, extension)
    
    # map expected figs to names for all region
    SN = ['network', 'marginal', 'depth', 'before_grouping', 'famsize', 'reffreq', 'lowreffreq', 'raw']
    SL = ['UMI_network_degree_{0}.', 'UMI_size_depth_marginal_distribution_{0}.',
          'Read_depth_per_umi_family_{0}.', 'UMI_freq_distribution_{0}.',
          'MeanFamilySize_{0}.', 'NonRefFreq_{0}.', 'NonRefFreq_low_freq_{0}.', 'RawConsensusDepth_{0}.']
          
    for i in regions:
        if i not in D:
            D[i] = {}
        for j in range(len(SN)):
            D[i][SN[j]] = GetExpectedFigure(FigDir, extension, SL[j].format(i))
        
    return D


def AddTitle(L, N, color, font_family, sample, renderer):
    '''
    (list, int, str, str, str, mistune.Markdown) -> None
    
    :param L: List with report strings
    :param N: Number of empty lines following title    
    :param color: Color of the title text
    :param font_family: Comma-separated text fonts
    :param sample: Sample name to appear in title
    :param renderer: markdown renderer
    
    Add title to list L
    '''
    
    # add title
    L.append(renderer('<font size=7><p style="text-align: center; color: {0}; font-family: {1}; font-weight:bold">Report for {2}</p></font>'.format(color, font_family, sample)))
    # add empty lines  
    L.append(renderer('<pre> </pre>' * N))
    

def CountMissingFiles(FigPaths):
    '''
    (dict) -> (int, int)

    :param FigPaths: Dictionary with paths to all expected figures (can be empty str)
    
    Return a tuple with count of valid and missing files
    '''
    
    # count valid and missing files
    missing, valid = 0, 0
   
    # make a list of non-region keys
    for i in FigPaths:
        if not i.startswith('chr'):
            if os.path.isfile(FigPaths[i]) == True:
                valid += 1
            else:
                missing +=1
        else:
            # make a list of non-regions keys
            for j in FigPaths[i]:
                if os.path.isfile(FigPaths[i][j]) == True:
                    valid += 1
                else:
                    missing +=1
    
    return valid, missing


def AddInfo(directory, L, N, color, font_family, FigPaths, renderer):
    '''
    (list, str, int, str, str, dict, mistune.Markdown) -> None
    
    :param directory: Directory with subfolders including Figures 
    :param L: List with report strings
    :param N: Number of empty lines following header
    :param color: Color of the text
    :param font_family: Comma-separated text fonts
    :param FigPaths: Dictionary with paths to all expected figures (can be empty str)
    :param renderer: markdown renderer    
        
    Add information about debarcer, stime stamp and working directory to list.
    Modify list in place
    '''
    
    # get the date
    date = '<b>time stamp:</b> ' + time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    # get debarcer version
    version = '<b>debarcer version:</b> ' + 'xxx'
    # get the directory containing subdirs 
    directory = '<b>directory:</b> ' + directory
    # count the number of valid and missing files
    valid, missing = CountMissingFiles(FigPaths)
    total = valid + missing
    files = '<b>figures:</b> ' + '{0} / {1} missing images'.format(missing, total) 
    text = '<br>'.join([date, version, directory, files])
    L.append(renderer('<pre><font size=3><p style="text-align: left; color: {0}; font-family: {1};">{2}</p></font></pre>'.format(color, font_family, text)))
    L.append(renderer('<pre> </pre>' * N))   


def AddHeader(L, N, color, num, font_family, text, renderer):
    '''
    (list, int, str, str, str, str, mistune.Markdown) -> int
    
    :param L: List with report strings
    :param N: Number of empty lines following header
    :param num: Header number
    :param color: Color of the text
    :param font_family: Comma-separated text fonts
    :param text: Text of the header
    :param renderer: markdown renderer
    
    Add header to list. Modify list in place. And return the header number 
    '''
    
    L.append('<font size=6><p style="text-align: left; color: {0}; font-weight: bold; font-family: {1};">{2}. {3}</p></font>'.format(color, font_family, num, text))
    #L.append('## {0}. {1}').format(num, text)     
    L.append('<pre> </pre>' * N)   
    return num

def AddSubheader(L, N, color, num1, num2, font_family, text):
    '''
    (list, int, str, str, str, str) -> int
    
    :param L: List with report strings
    :param N: Number of empty lines following header
    :param color: Color of the text
    :param num1: Header number
    :param num2: Sub-header number
    :param font_family: Comma-separated text fonts
    :param text: Text of the header
    
    Add header to list. Modify list in place. And return the sub-header number
    '''
    
    L.append('<font size=4><p style="text-align: left; color: {0}; font-weight: bold; font-family: {1};">{2}.{3} {4}</p></font>'.format(color, font_family, num1, num2, text))
    #L.append('## {0}. {1}').format(num, text)     
    L.append('<pre> </pre>' * N)   
    return num2


def AddCoverageFig(L, CovStats, DataFiles, font_family, mincov, extension, FigPaths, figcounter, N):
    '''
    (list, str, str, dict, int)- > int
    
    :param L: List with report strings
    :param font_family: Comma-separated text fonts
    :param extension: Extension of the figure files
    :param FigPaths: Dictionary with paths to all expected figures (can be empty str)
    :param figcounter: Figure number
    :param N: Number of empty lines following last legend
    
    Add figures and legends to list L or a warning if figures don't exist
    and return the number of next figure
    '''

    # scaling factor
    scale = 0.7

    # add a warning if both expected files are missing
    if FigPaths['coverage'] == '':
        missing = 'Coverage_Umi_Count.{0}'.format(extension)
        L.append('<p style="color: Tomato;text-align: left; font-family: Arial, sans-serif; font-weight=bold;">[Warning]<br> Missing expected files:<br>{0} </p>'.format(missing)) 
        L.append('<pre> </pre>')
    
    if FigPaths['coverage'] != '':
        # resize figure
        height, width = ResizeFifure(FigPaths['coverage'], scale)
        # encode base64 image
        encoded_fig = EncodeImage(FigPaths['coverage'])
        # add image and legend
        images = '<img style="padding-right:0px; padding-left:30px" src="data:image/png;base64,{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(encoded_fig, 'coverage', width, height)
        L.append(images)
        # add legend
        legends = '<span style="padding-right: 70px; padding-left:10px; font-family:{0}; font-size:16px"> <b>Figure {1}</b>. Average read depth and umi counts per genomic interval</span>'.format(font_family, figcounter)
        L.append(legends)
        # update figure counter
        figcounter += 1
    # append empty line
    L.append('<pre> </pre>' * N)
    
    
    # add coverage fig
    source_fig = PlotDataPerRegion(CovStats, DataFiles, mincov=mincov, datatype='coverage')
    L.append('<embed type="image/svg+xml" src= {0} width="{1}" height="{2}" />'.format(source_fig, 800, 800))
    
    
    return figcounter


def AddPreprocessingFigs(L, font_family, extension, FigPaths, figcounter, N, renderer):
    '''
    (list, str, str, dict, int, mistune.Markdown)- > int
    
    :param L: List with report strings
    :param font_family: Comma-separated text fonts
    :param extension: Extension of the figure files
    :param FigPaths: Dictionary with paths to all expected figures (can be empty str)
    :param figcounter: Figure number
    :param N: Number of empty lines following last legend
    :param renderer: markdow renderer
    
    Add figures and legends to list L or a warning if figures don't exist
    and return the number of next figure
    '''
    
    # keys to access figure files         
    keys = ['reads', 'preprocessing']
    # scaling factor for resizing figures
    scale = [0.65, 0.7]
    # alternate names on html page
    altfig = ['processed reads', 'umi frequency']
    
    # map expected files to key
    d = {'reads':'Proportion_correct_reads.{0}'.format(extension), 'preprocessing': 'UMI_occurence_preprocessing.{0}'.format(extension)}
    # add a warning if both expected files are missing
    if FigPaths['reads'] == '' or FigPaths['preprocessing'] == '':
        missing = '<br>'.join([d[i] for i in d if FigPaths[i] == ''])
        L.append('<p style="color: Tomato;text-align: left; font-family: Arial, sans-serif; font-weight=bold;">[Warning]<br> Missing expected files:<br>{0} </p>'.format(missing)) 
        L.append('<pre> </pre>')
    
    # add figures inline 
    images = ''    
    fignum = {}
    
    for i in range(len(keys)):
        if FigPaths[keys[i]] != '':
            # resize image
            height, width = ResizeFifure(FigPaths[keys[i]], scale[i])
            # encode base64 image
            encoded_fig = EncodeImage(FigPaths[keys[i]])
        # add image and legend
            if i == 0:
                images += '<img style="padding-right: 100px; padding-left:30px" src="data:image/png;base64,{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(encoded_fig, altfig[i], width, height)
            else:
                images += '<img style="padding-left:30px" src="data:image/png;base64,{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(encoded_fig, altfig[i], width, height)
            #update figure counter
            fignum[keys[i]] = figcounter
            figcounter += 1
    L.append(images)
       
    # add legends
    legends = ''
    for i in range(len(keys)):
        if FigPaths[keys[i]] != '':
            if i == 0:
                legends += '<span style="padding-right: 70px; padding-left:10px; font-family:{0}; font-size:16px"> <b>Figure {1}</b>. Number of reads with correct and incorrect umi-spacer configuration </span>'.format(font_family,fignum[keys[i]])
            else:
                legends += '<span style="padding-left:10px; font-family:{0}; font-size:16px"> <b>Figure {1}</b>. Frequency distribution of umis with correct configurarion </span>'.format(font_family, fignum[keys[i]])
    L.append(renderer(legends))
    
    # append empty line
    L.append(renderer('<pre> </pre>' * N))
    return figcounter
    

def AddSpacerLine(L, renderer):
    '''
    (list) -> None
    
    :param L: List with report strings
    :param renderer: markdown renderer
    
    Add a black line with empty space before and after 
    '''
    # add 1 empty line before and after spacer
    L.append('<pre> </pre>')
    # add spacer line
    L.append(renderer('****'))
    L.append('<pre> </pre>')


def grouper(iterable, n, fillvalue=None):
    '''
    (iterable, int, NoneType) -> iterable
        
    Return an iterable of 2 item tuples, with None if length of iterable is odd 
    '''
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def AddBeforeGroupingSection(L, font_family, extension, FigPaths, figcounter, N):
    '''
    (list, str, str, dict, int)- > int
    
    :param L: List with report strings
    :param font_family: Comma-separated text fonts
    :param extension: Extension of the figure files
    :param FigPaths: Dictionary with paths to all expected figures (can be empty str)
    :param figcounter: Figure number
    :param N: Number of empty lines following last legend
    
    Add figures and legends to list L or a warning if figures don't exist
    and return the number of next figure
    '''
    
    # add description of the figures
    style = 'text-align: justify; text-justify: inter-word; padding-right: 20px;\
    padding-left:10px; font-family:{0}; font-size:18px'.format(font_family)
    L.append('<p style="{0}">Density plots showing the frequency distribution of\
             umi sequences<br>(ie. the number of times a given umi sequence is observed)</p>'.format(style))
    L.append('<pre> </pre>')
    
    # make a sorted list of regions for 'before_grouping' figures
    keys = sorted([i for i in FigPaths.keys() if 'chr' in i])
    
    # add warning for missing files
    missing = '<br>'.join(['UMI_freq_distribution_{0}.{1}'.format(i, extension) for i in keys if FigPaths[i]['before_grouping'] == ''])
    if len(missing) != 0:
        L.append('<p style="color: Tomato;text-align: left; font-family: Arial, sans-serif; font-weight=bold;">[Warning]<br> Missing expected files:<br>{0} </p>'.format(missing)) 
        L.append('<pre> </pre>')
    
    # define scaling factor and set alternate figure name
    scale, altfig = 0.8, 'before grouping'
        
    # make pairs of non-empty figure pairs
    Files = [FigPaths[keys[i]]['before_grouping'] for i in range(len(keys)) if FigPaths[keys[i]]['before_grouping'] != '']
    # group files by 2 
    Files = list(map(lambda x: list(x), list(grouper(Files, 2))))
    for i in range(len(Files)):
        if None in Files[i]:
            Files[i].pop(Files[i].index(None))
    
    # add images and legends for valid files    
    for i in range(len(Files)):
        # map files to regions
        regions = []
        # store images and figure number for given region
        images, fignum = '', {}    
        for j in range(len(Files[i])):
            region = os.path.basename(Files[i][j])
            region = region[region.index('_chr')+ 1: region.rindex('.')]
            regions.append(region)
            # resize image
            height, width = ResizeFifure(Files[i][j], scale)
            # encode base64 image
            encoded_fig = EncodeImage(Files[i][j])
            # add image and legend
            if j == 0:
                images += '<img style="padding-right: 30px; padding-left:10px" src="data:image/png;base64,{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(encoded_fig, altfig, width, height)
            else:
                images += '<img style="padding-left:10px" src="data:image/png;base64,{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(encoded_fig, altfig, width, height)
            #update figure counter
            fignum[region] = figcounter
            figcounter += 1
        L.append(images)
        # add legends
        legends = ''
        for j in range(len(regions)):
            if j == 0:
                legends += '<span style="padding-right: 180px; padding-left:20px; font-family:{0}; font-size:16px"> <b>Figure {1}</b>. Interval {2} </span>'.format(font_family,fignum[regions[j]], regions[j])
            else:
                legends += '<span style="padding-left:20px; font-family:{0}; font-size:16px"> <b>Figure {1}</b>. Interval {2} </span>'.format(font_family, fignum[regions[j]], regions[j])
            #update figure counter
            #figcounter += 1
        L.append(legends)
    
        # append empty line
        L.append('<pre> </pre>' * N)
    
    return figcounter
        


def AddGrouping(L, font_family, extension, FigPaths, figcounter, N, num):
    '''
    (list, str, str, dict, int)- > int
    
    :param L: List with report strings
    :param font_family: Comma-separated text fonts
    :param extension: Extension of the figure files
    :param FigPaths: Dictionary with paths to all expected figures (can be empty str)
    :param figcounter: Figure number
    :param N: Number of empty lines following last legend
    :param num: Number of the header of the Grouping section
    
    Add figures and legends to list L or a warning if figures don't exist
    and return the number of next figure
    '''
    
    # add sub-header
    subnum = AddSubheader(L, 1, 'black', num, 1, font_family, 'QC plots across regions')
            
    # add description of the figures
    L.append('<p style="text-align: left; font-family:{0}; font-size:18px">Number of parent and children umis\
             following family grouping for each genomic interval</p>'.format(font_family))
    L.append('<pre> </pre>')

    # 1. add figures of umi counts across regions

    # keys to access figures in this order
    keys = ['total', 'children', 'ratio', 'interval', 'freq']
    Maps = {'total':['Total_Umis.', 'Umi counts per genomic interval', 0.75, 'total umis'],
            'children':['Children_Umis.', 'Count of children umis per interval', 0.75, 'children umis'],
            'ratio':['Child_Parent_Umis_Ratio.', 'Child to parent umi ratios per interval', 0.75, 'ratio'],
            'interval':['PTU_vs_CTU.', 'Relationship between parent and children umi counts with interval sizes', 0.85, 'interval size'],
            'freq':['Children_vs_ParentFreq.', 'Frequency of parent umis with a given number of children for each interval', 0.85, 'frequency']}
    
    # add warning for missing files
    missing = '<br>'.join([Maps[i][0] for i in Maps if FigPaths[i] == '']).format(extension)
    if len(missing) != 0:
        L.append('<p style="color: Tomato;text-align: left; font-family: Arial, sans-serif; font-weight=bold;">[Warning]<br> Missing expected files:<br>{0} </p>'.format(missing)) 
        L.append('<pre> </pre>')

    # make groups of non-empty figure pairs
    group_fig = 2
    Lgds = [Maps[keys[i]][1] for i in range(len(keys)) if FigPaths[keys[i]] != '']
    Files = [FigPaths[keys[i]] for i in range(len(keys)) if FigPaths[keys[i]] != '']
    ScalingFactors = [Maps[keys[i]][2] for i in range(len(keys)) if FigPaths[keys[i]] != '']
    AltNames = [Maps[keys[i]][3] for i in range(len(keys)) if FigPaths[keys[i]] != '']
    # group files, legends and scaling factors
    Lgds = list(map(lambda x: list(x), list(grouper(Lgds, group_fig))))
    Files = list(map(lambda x: list(x), list(grouper(Files, group_fig))))
    ScalingFactors = list(map(lambda x: list(x), list(grouper(ScalingFactors, group_fig))))
    AltNames = list(map(lambda x: list(x), list(grouper(AltNames, group_fig))))
    for i in range(len(Files)):
        if None in Files[i]:
            Files[i].pop(Files[i].index(None))
        if None in Lgds[i]:
            Lgds[i].pop(Lgds[i].index(None))
        if None in ScalingFactors[i]:
            ScalingFactors[i].pop(ScalingFactors[i].index(None))
        if None in AltNames[i]:
            AltNames[i].pop(AltNames[i].index(None))
        
    # add images and legends for valid files    
    for i in range(len(Files)):
        # store images and figure number for given region
        images, fignum = '', []    
        for j in range(len(Files[i])):
            # resize image
            height, width = ResizeFifure(Files[i][j], ScalingFactors[i][j])
            # encode base64 image
            encoded_fig = EncodeImage(Files[i][j])
            # add images
            if j == 0:
                images += '<img style="padding-right: 30px; padding-left:10px" src="data:image/png;base64,{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(encoded_fig, AltNames[i][j], width, height)
            else:
                images += '<img style="padding-left:10px" src="data:image/png;base64,{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(encoded_fig, AltNames[i][j], width, height)
            #update figure counter
            fignum.append(figcounter)
            figcounter += 1
        L.append(images)
        # add legends
        legends = ''
        for j in range(len(Lgds[i])):
            if j == 0:
                legends += '<span style="padding-right: 190px; padding-left:30px; font-family:{0}; font-size:16px"> <b>Figure {1}</b>. {2}</span>'.format(font_family,fignum[j], Lgds[i][j])
            else:
                legends += '<span style="padding-left:30px; font-family:{0}; font-size:16px"> <b>Figure {1}</b>. {2}</span>'.format(font_family, fignum[j], Lgds[i][j])
        L.append(legends)
        # append empty line
        L.append('<pre> </pre>')
    
    # 2. Add figures specific to each region
    subnum = AddSubheader(L, 1, 'black', num, subnum + 1, font_family, 'Region-specific QC plots')

    intro = ['Degree distribution (left panel) shows the number of edges between umi nodes<br>defined by the hamming distance between umi sequences.\
             Network shows the<br> interaction among umi nodes colored by degree (right panel)',
             'Marginal plots show the relationship between read depth and umis per group',
             'Read depth distribution at positions of highest and lower abundance, shown<br>as proportion of the read depth within family group']
    for i in intro:
        L.append('<ul><li color:black><span style="list-type-position:outside;\
                 list-style-type:circle; display:list-item; text-align: left; padding-right: 10px;\
                 padding-left:12px; font-family:{0}; font-size:18px">{1}</span></li></ul>'.format(font_family, i))
    L.append('<pre> </pre>')

    # make a sorted list of regions for 'grouoing' figures
    regions = sorted([i for i in FigPaths.keys() if 'chr' in i])
    
    # add warning for missing files
    # keys to access figures in this order
    keys = ['network', 'marginal', 'depth']
    # map keys to expected file names, legends, scaling factors and alternate figure names
    Maps = {'network':['UMI_network_degree_{0}.{1}', 'Node degree distribution', 0.6, 'network'],
            'marginal':['UMI_size_depth_marginal_distribution_{0}.{1}', 'Marginal plot', 0.7, 'marginal'],
            'depth':['Read_depth_per_umi_family_{0}.{1}', 'Read depth within group', 0.6, 'read depth']}
        
    
    missing = '<br>'.join([Maps[i][0].format(j, extension) for i in Maps for j in regions if FigPaths[j][i] == ''])
    if len(missing) != 0:
        L.append('<p style="color: Tomato;text-align: left; font-family: Arial, sans-serif; font-weight=bold;">[Warning]<br> Missing expected files:<br>{0} </p>'.format(missing)) 
        L.append('<pre> </pre>')
    
    # make groups of non-empty figure pairs
    # keep groups of files together per region
    Files, Lgds, ScalingFactors, AltNames, Intervals = [], [], [], [], []
    for i in regions:
        f, l, s, a = [], [], [], []
        for j in keys:
            if FigPaths[i][j] != '':
                l.append(Maps[j][1])
                f.append(FigPaths[i][j])
                s.append(Maps[j][2])
                a.append(Maps[j][3])
        if len(f) != 0:
            Files.append(f)
            Lgds.append(l)
            ScalingFactors.append(s)
            AltNames.append(a)
            Intervals.append(i)
       
    # add images and legends for valid files    
    for i in range(len(Files)):
        # add sentence about the genomic interval
        L.append('<pre> </pre>')
        L.append('<ul><li color:black><p style="color:black; display:list-item;\
                 list-style-type:square; text-align: left; font-family: Arial,\
                 sans-serif; font-weight=normal;"><b>Interval {0}</b></p></li></ul>'.format(Intervals[i]))
        # store images and figure number for given region
        images, fignum = '', []    
        for j in range(len(Files[i])):
            # resize image
            height, width = ResizeFifure(Files[i][j], ScalingFactors[i][j])
            # encode base64 image
            encoded_fig = EncodeImage(Files[i][j])
            # add images
            if j == 0:
                images += '<img style="padding-right: 10px; padding-left:10px" src="data:image/png;base64,{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(encoded_fig, AltNames[i][j], width, height)
            else:
                images += '<img style="padding-left:10px" src="data:image/png;base64,{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(encoded_fig, AltNames[i][j], width, height)
            #update figure counter
            fignum.append(figcounter)
            figcounter += 1
        L.append(images)
        # add legends
        legends = ''
        for j in range(len(Lgds[i])):
            if j == 0:
                padding_right,padding_left = 210, 10
            elif j == 1:
                padding_right, padding_left = 180, 10
            else:
                padding_right, padding_left = 0, 10
            legends += '<span style="padding-right: {0}px; padding-left:{1}px; font-family:{2}; font-size:16px"> <b>Figure {3}</b>. {4}</span>'.format(padding_right, padding_left, font_family,fignum[j], Lgds[i][j])
        L.append(legends)
        # append empty line
        L.append('<pre> </pre>')

    return figcounter


def AddCollapsing(L, font_family, extension, FigPaths, figcounter, N, num):
    '''
    (list, str, str, dict, int)- > int
    
    :param L: List with report strings
    :param font_family: Comma-separated text fonts
    :param extension: Extension of the figure files
    :param FigPaths: Dictionary with paths to all expected figures (can be empty str)
    :param figcounter: Figure number
    :param N: Number of empty lines following last legend
        
    Add figures and legends to list L or a warning if figures don't exist
    and return the number of next figure
    '''
    
    # make a sorted list of regions for 'grouoing' figures
    regions = sorted([i for i in FigPaths.keys() if 'chr' in i])
    # keys to access figures in this order
    keys = ['famsize', 'raw', 'reffreq', 'lowreffreq']
    Maps = {'famsize':['MeanFamilySize_{0}.{1}', 'Mean family size', 0.7, 'family size'],
            'reffreq':['NonRefFreq_{0}.{1}', 'Frequency of alternative variants', 0.85, 'alternative variants'],
            'raw':['RawConsensusDepth_{0}.{1}', 'Raw and consensus depth', 0.6, 'read depth'],
            'lowreffreq': ['NonRefFreq_low_freq_{0}.{1}', 'Frequency alternative variants (Y axis cut)', 0.85, 'low frequency']}
    
    # Add figures specific to each region
    subnum = AddSubheader(L, 1, 'black', num, 1, font_family, 'Region-specific QC plots')
        
    intro = ['Average umi family size across a given genomic interval<br>for various minimum size thresholds',
             'Read depth without and with umi collapsing for various<br>minimum family size thresholds',
             'Frequency of alternative variants without and with collapsing<br>on various minimum umi family size thresholds',
             'Focus on low frequency alternative variant']
    
    for i in intro:
        L.append('<ul><li color:black><span style="list-type-position:outside;\
                 list-style-type:circle; display:list-item; text-align: left; padding-right: 10px;\
                 padding-left:12px; font-family:{0}; font-size:18px">{1}</span></li></ul>'.format(font_family, i))
    L.append('<pre> </pre>')
    
    # add warning for missing files
    # map keys to expected file names, legends, scaling factors and alternate figure names
    missing = '<br>'.join([Maps[i][0].format(j, extension) for i in Maps for j in regions if FigPaths[j][i] == ''])
    if len(missing) != 0:
        L.append('<p style="color: Tomato;text-align: left; font-family: Arial, sans-serif; font-weight=bold;">[Warning]<br> Missing expected files:<br>{0} </p>'.format(missing)) 
        L.append('<pre> </pre>')
    
    ## add famsize and raw inline and reffreq and lowreffreq as separate line
    for i in range(len(regions)):
        # check if some files are present
        if len([FigPaths[regions[i]][keys[j]] for j in range(len(keys)) if FigPaths[regions[i]][keys[j]] != '']) != 0:
            # write interval
            L.append('<pre> </pre>')
            L.append('<ul><li color:black><p style="color:black; display:list-item;\
                     list-style-type:square; text-align: left; font-family: Arial,\
                     sans-serif; font-weight=normal;"><b>Interval {0}</b></p></li></ul>'.format(regions[i]))    
            # add famsize and raw figures if they exist
            f, l, s, a = [], [], [], []
            images, fignum = '', []    
            for j in range(len(keys)-2):
                if FigPaths[regions[i]][keys[j]] != '':
                    l.append(Maps[keys[j]][1])
                    f.append(FigPaths[regions[i]][keys[j]])
                    s.append(Maps[keys[j]][2])
                    a.append(Maps[keys[j]][3])
            if len(f) != 0:
                #L.append('<pre> </pre>')
                for j in range(len(f)):
                    # resize image
                    height, width = ResizeFifure(f[j], s[j])
                    # encode base64 image
                    encoded_fig = EncodeImage(f[j])
                    # add images
                    if j == 0:
                        images += '<img style="padding-right: 10px; padding-left:10px" src="data:image/png;base64,{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(encoded_fig, a[j], width, height)
                    else:
                        images += '<img style="padding-left:10px" src="data:image/png;base64,{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(encoded_fig, a[j], width, height)
                    #update figure counter
                    fignum.append(figcounter)
                    figcounter += 1
                L.append(images)
                # add legends
                legends = ''
                for j in range(len(l)):
                    if j == 0:
                        padding_right,padding_left = 440, 10
                    else:
                        padding_right, padding_left = 0, 10
                    legends += '<span style="padding-right: {0}px; padding-left:{1}px; font-family:{2}; font-size:16px"> <b>Figure {3}</b>. {4}</span>'.format(padding_right, padding_left, font_family,fignum[j], l[j])
                L.append(legends)
                # append empty line
                L.append('<pre> </pre>')
            
            
            # add reffreq and lowreffreq if they exist
            f, l, s, a = [], [], [], []
            images, fignum = '', []   
            for j in range(2, len(keys)):
                if FigPaths[regions[i]][keys[j]] != '':
                    l.append(Maps[keys[j]][1])
                    f.append(FigPaths[regions[i]][keys[j]])
                    s.append(Maps[keys[j]][2])
                    a.append(Maps[keys[j]][3])
            if len(f) != 0:
                for j in range(len(f)):
                    # resize image
                    height, width = ResizeFifure(f[j], s[j])
                    # encode base64 image
                    encoded_fig = EncodeImage(f[j])
                    # add images
                    if j == 0:
                        images += '<img style="padding-right: 10px; padding-left:10px" src="data:image/png;base64,{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(encoded_fig, a[j], width, height)
                    else:
                        images += '<img style="padding-left:10px" src="data:image/png;base64,{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(encoded_fig, a[j], width, height)
                    #update figure counter
                    fignum.append(figcounter)
                    figcounter += 1
                L.append(images)
                # add legends
                legends = ''
                for j in range(len(l)):
                    if j  == 0:
                        padding_right,padding_left = 440, 10
                    else:
                        padding_right, padding_left = 0, 10
                    legends += '<span style="padding-right: {0}px; padding-left:{1}px; font-family:{2}; font-size:16px"> <b>Figure {3}</b>. {4}</span>'.format(padding_right, padding_left, font_family,fignum[j], l[j])
                L.append(legends)
                # append empty line
                L.append('<pre> </pre>')
              
    return figcounter



def GetSampleName(directory, **Options):
    '''
    (str, dict) -> str

    :param directory: Directory with subfolders including Figures
    :param Options: Optional parameters. Accepted values: 'sample'

    Return the sample name, either from     
    
    '''
    
    if 'sample' in Options:
        sample = Options['sample']
    else:
        sample = ''
    if sample == '':
        sample = os.path.basename(directory) 
    return sample
    

def WriteReport(directory, CovStats, DataFiles, extension, Outputfile, mincov, renderer=mistune.Markdown(), **Options):
    '''
    (str, str, list, str, str, float, mistune.Markdown, dict) -> None
    
    :param directory: Directory with subfolders including Figures
    :param CovStats: yaml file with mean read depth per region
    :param DataFiles:  List of .csv fles with umi counts 
    :param extension: Extension of the figure files
    :param Outputfile: Name of the html report
    :param mincov: Minimum read depth to label regions
    :param renderer: markdown renderer
    :param Options: Optional parameters. Accepted values: 'sample'
        
    Write an html report of debarcer analysis for a given sample
    '''
    
    # set up font family <- string with multiple values. browser will use values from left to right if not defined  
    font_family = 'Arial, Verdana, sans-serif'
        
    # set up counter for figure files
    figcounter = 1
        
    # get the path of all expected figures (can be '')
    FigPaths = ListExpectedFigures(directory, extension)
       
    # use list to store lines of report
    L = []
    
    ## Title section
    ## get sample name
    sample = GetSampleName(directory, **Options)
    AddTitle(L, 3, 'black', font_family, sample, renderer)
    
    ## Add debarcer info and time stamp
    AddInfo(directory, L, 2, 'black', font_family, FigPaths, renderer)
    
    ## Pre-processing section
    headernum = AddHeader(L, 1, 'black', 1, font_family, 'Pre-processing', renderer)
    # add figures from pre-processing, update figure counter 
    figcounter = AddPreprocessingFigs(L, font_family, 'png', FigPaths, figcounter, 1)
    # add spacer line
    AddSpacerLine(L, renderer)
    
    ## Coverage section
    headernum = AddHeader(L, 1, 'black', headernum+1, font_family, 'Coverage')
    # add figures from pre-processing, update figure counter 
    figcounter = AddCoverageFig(L, CovStats, DataFiles, font_family, mincov, 'png', FigPaths, figcounter, 1) 
    
    # add spacer line
    AddSpacerLine(L, renderer)
    
    
    ########
    
    
    #PlotDataPerRegion(CovStats, DataFiles, outputfile=os.path.join(FigDir, 'Child_Parent_Umis_Ratio'), mincov=args.minratio, datatype='ratio')

    #PlotDataPerRegion(CovStats, DataFiles, outputfile=os.path.join(FigDir, 'Total_Umis'), mincov=args.minumis, datatype='umis')

    #PlotDataPerRegion(CovStats, DataFiles, outputfile=os.path.join(FigDir, 'Children_Umis'), mincov=args.minchildren, datatype='children')

    
    
    ###########
    
    
    ## Pre-grouping section
    headernum = AddHeader(L, 1, 'black', headernum+1, font_family, 'Umi distribution before family grouping')
    # add figures from pre-grouping QC, update figure counter     
    figcounter = AddBeforeGroupingSection(L, font_family, extension, FigPaths, figcounter, 1)
    # add spacer line
    AddSpacerLine(L, renderer)
         
    ## Grouping section
    headernum = AddHeader(L, 1, 'black', headernum+1, font_family, 'Umi Grouping')
    # add figures from Grouping section and update figure counter
    figcounter = AddGrouping(L, font_family, extension, FigPaths, figcounter, 1, headernum)
    # add spacer line
    AddSpacerLine(L, renderer)
    
    ## Collapsing section
    headernum = AddHeader(L, 1, 'black', headernum+1, font_family, 'Umi family Collapsing')
    # add figures from Collapsing section and update figure counter
    figcounter = AddCollapsing(L, font_family, extension, FigPaths, figcounter, 1, headernum)

    # create report string
    S = ''.join(L)
    newfile = open(Outputfile, 'w')
    newfile.write(S)
    newfile.close()
    
