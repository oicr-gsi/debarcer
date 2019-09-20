# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 12:05:32 2019

@author: rjovelin
"""

import mistune
import os
import scipy.ndimage
from xhtml2pdf import pisa 

# using ipython and html to generate report 
# https://plot.ly/python/v3/pdf-reports/


class MyCustomRenderer(mistune.Renderer):
    def header(self, text, level, raw=None):
        return "<h%s class='my-cusom-header-cls'>%s</h%s>" % (level, text, level)
#    def image(self, src, title, alt_text):
#        return "<img src='%s' alt='%s' class='img-responsive'>" % (src, alt_text)
#    
    
    
    def image(self, src, title, text):
        if src.startswith('javascript:'):
            src = ''
        text = mistune.escape(text, quote=True)
        if title:
            title = mistune.escape(title, quote=True)
            html = '<img class="img-responsive center-block" src="%s" alt="%s" title="%s"' % (src, text, title)
        else:
            html = '<img class="img-responsive center-block" src="%s" alt="%s"' % (src, text)
        if self.options.get('use_xhtml'):
            return '%s />' % html
        return '%s>' % html

    
    
    

    
    
    
    
    
    def link(self, link, title, content):
        if link.lower().startswith('http://'):
            link = 'https://' + link[len('http://'):]
        elif link.lower().startswith('https://'):
            pass # we do not need any processing here
        else:
            # so it is a relative link.
            # We can add a domain in front of it.
            link = link.lstrip('/')
            link = 'https://example.com/' + link
        return "<a href='%s' title='%s'>%s</a>" % (link, title, content)




def convert_html_to_pdf(source_html, output_filename):
    # open output file for writing (truncated binary)
    result_file = open(output_filename, "w+b")

    # convert HTML to PDF
    pisa_status = pisa.CreatePDF(
            source_html,                # the HTML to convert
            dest=result_file)           # file handle to recieve result

    # close output file
    result_file.close()                 # close output file

    # return True on success and False on errors
    return pisa_status.err  




class MyCustomRenderer2(mistune.Renderer):
    pass



def GetExpectedFigure(FigDir, extension, expected_name):
    '''
    
    
    
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
    SN = ['network', 'marginal', 'depth', 'before_grouping', 'famsize', 'reffreq', 'raw']
    SL = ['UMI_network_degree_{0}.', 'UMI_size_depth_marginal_distribution_{0}.',
          'Read_depth_per_umi_family_{0}.', 'UMI_freq_distribution_{0}.',
          'MeanFamilySize_{0}.', 'NonRefFreq_{0}.', 'RawConsensusDepth_{0}.']
               
    for i in regions:
        if i not in D:
            D[i] = {}
        for j in range(len(SN)):
            D[i][SN[j]] = GetExpectedFigure(FigDir, extension, SL[j].format(i))
        
    return D



def AddImage(FigPaths, Report, keys, scale, legends, altfig, figcounter):
    '''
    
    
    '''

    for i in range(len(keys)):
        if FigPaths[keys[i]] != '':
            # get original size 
            height, width, channels = scipy.ndimage.imread(FigPaths[keys[i]]).shape
            # rescale
            height, width = list(map(lambda x: x * scale[i], [height, width]))
            # add image and legend
            Report.append('<img src="{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(FigPaths[keys[i]], altfig[i], width, height))
            Report.append(legends[i].format(figcounter))
            #update figure counter
            figcounter += 1
    
    return Report, figcounter










    

def WriteReport(directory, extension, **Options):
    '''
    
    
    
    '''
    
    
    
    ## issues:
    # use specific rendered for lines or blocks
    # center text and image
    # add empty paragraphs
    # add multiple images per line. can add dseveral image separated by empty string
    # but string cannot start with white space
#    truc = ''
#    
#    for i in range(len(keys)):
#        if FigPaths[keys[i]] != '':
#            # get original size 
#            height, width, channels = scipy.ndimage.imread(FigPaths[keys[i]]).shape
#            # rescale
#            height, width = list(map(lambda x: x * scale[i], [height, width]))
#            if i == 0:
#                truc = truc + '<img src="{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(FigPaths[keys[i]], altfig[i], width, height)
#            else:
#                truc = truc + '     ' + '<img src="{0}" alt="{1}" title="{1}" width="{2}" height="{3}" />'.format(FigPaths[keys[i]], altfig[i], width, height)
#            figcounter += 1 
#    L.append(truc)    
#    
    
    
    keys = ['reads',  'preprocessing', 'coverage', 'ratio', 'total', 'children', 'interval', 'freq']
    region_keys = ['network', 'marginal', 'depth', 'before_grouping', 'famsize', 'reffreq', 'raw']
    
    
      
    
    
    
    
    
    
    
    
    
    # set up counter for figure files
    figcounter = 1
    
    
    
    # get the path of all expected figures (can be '')
    FigPaths = ListExpectedFigures(directory, extension)
    
    #renderer = mistune.Renderer(escape=True, hard_wrap=True)
    #renderer = mistune.Renderer()
    renderer = MyCustomRenderer2()
    markdown = mistune.Markdown(renderer=renderer)
    
    # can't use multi-lines string for some unknown reasons. use list to store lines of report
    L = []
    
    # add level 1 title
    # get sample name from directory or options
    if 'sample' in Options:
        sample = Options['sample']
    else:
        sample = os.path.basename(directory)
    L.append('# **debarcer Report - {0}**'.format(sample))
                 
    # add level 2 title
    L.append('## 1. Pre-processing')     
     
    # keys to access figure files         
    keys = ['reads', 'preprocessing']
    # scaling factor
    scale = [0.65, 0.7]
    # legends of the 2 files resulting from pre-proccessing QC
    legends = ['**Figure {0}**. Number of reads with correct and incorrect umi-spacer configuration',
               '**Figure {0}**. Frequency distribution of umis with correct configurarion']
    # alternate names on html page
    altfig = ['processed reads', 'umi frequency']
    # add images and legends  
    L, figcounter = AddImage(FigPaths, L, keys, scale, legends, altfig, figcounter)
    
    
    # add level 2 title
    L.append('## 2. Umi distribution before family grouping')     
    
    keys = sorted([i for i in FigPaths.keys() if 'chr' in i])
    scale = [0.65]
    legends = ['**Figure {0}**. Number of reads with correct and incorrect umi-spacer configuration']
    # alternate names on html page
    altfig = ['before_grouping']
    for i in range(len(keys)):
        L, figcounter = AddImage(FigPaths[keys[i]], L, ['before_grouping'], scale, legends, altfig, figcounter)
    
    
    
    L.append('****')     
    
    
    
    
    # add level 2 title
    L.append('## 3. Umi Grouping')     
    
    # keys to access figure files         
    keys = ['total', 'children', 'ratio', 'interval', 'freq']
    # scaling factor
    scale = [0.63, 0.63, 0.63, 0.68, 0.8] 
    # legends of the 4 files resulting from grouping QC
    legends = ['**Figure {0}**. Frequency distribution of umis with correct configurarion',
               '**Figure {0}**. Frequency distribution of umis with correct configurarion',
               '**Figure {0}**. Frequency distribution of umis with correct configurarion',
               '**Figure {0}**. Frequency distribution of umis with correct configurarion',
               '**Figure {0}**. Frequency distribution of umis with correct configurarion']
    # alternate names on html page
    altfig = ['total umis', 'children umis', 'ratio', 'interval size', 'frequency']
    
    # add images and legends
    L, figcounter = AddImage(FigPaths, L, keys, scale, legends, altfig, figcounter)
    
    keys = sorted([i for i in FigPaths.keys() if 'chr' in i])
    scale = [0.65, 0.85, 0.65]
    legends = ['**Figure {0}**. Number of reads with correct and incorrect umi-spacer configuration'] * 3
    # alternate names on html page
    altfig = ['network', 'marginal', 'depth']
    for i in range(len(keys)):
        L, figcounter = AddImage(FigPaths[keys[i]], L, ['network', 'marginal', 'depth'], scale, legends, altfig, figcounter)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # add level 2 title
    L.append('## 4. Umi family Collapsing')
    
    # keys to access figure files         
    keys = ['coverage']
    # scaling factor
    scale = [0.7]
    # legends of the 4 files resulting from grouping QC
    legends = ['**Figure {0}**. Number of reads with correct and incorrect umi-spacer configuration']
    # asternate names on html page
    altfig = ['coverage']
        
    L, figcounter = AddImage(FigPaths, L, keys, scale, legends, altfig, figcounter)
    
    keys = sorted([i for i in FigPaths.keys() if 'chr' in i])
    scale = [0.65, 0.85, 0.65]
    legends = ['**Figure {0}**. Number of reads with correct and incorrect umi-spacer configuration'] * 3
    # alternate names on html page
    altfig = ['famsize', 'reffreq', 'raw']
    for i in range(len(keys)):
        L, figcounter = AddImage(FigPaths[keys[i]], L, ['famsize', 'reffreq', 'raw'], scale, legends, altfig, figcounter)
    
    # create report string
    S = ''.join([markdown(i) for i in L])
     
    newfile = open('index.html', 'w')
    newfile.write(S)
    newfile.close()
    
    return S


S = WriteReport('./', 'png', sample = 'sample1')  


convert_html_to_pdf(S, 'report.pdf')