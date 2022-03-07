# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 12:40:28 2022

Workbook to compare analysis of angle-dependences of different samples.

Requires BRimmler/analysis-modules library (add to PATH)

@author: rimmler
"""

# ____________________________________________________________________________
# SETTINGS
# Project E2
# ipFileDirNames = [
#     r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-2\220131\D1_0deg\02_angle-dependence\fittingOutput\angleDependence\c-free\fitparams_summary.csv',
#     r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-2\220131\D3_45deg\01_angle-dependence\fittingOutput\angleDependence\c-free\fitparams_summary.csv',
#     r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2960-2\220202\D1_0deg\003_angle-dependence\fittingOutput\angleDependence\c-free\fitparams_summary.csv',
#     r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2960-2\220203\D4_90deg\002_angle-dependence\pos_field\fittingOutput\angleDependence\c-free\fitparams_summary.csv'
#     ]

# sampleNames = ['MA2959-2-D1', 'MA2959-2-D3', 'MA2960-2-D1', 'MA2960-2-D4']

# fit_comps_to_show = ['xyz']
# analysisMode = 1 # Analysis mode used in stfmrAnglePlot for fitting of angle-dependence

# opDir = r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\ProjE-synthesis\220204'

# Project E2: MSN(001) before and after magsats
ipFileDirNames = [
    r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-2\220131\D1_0deg\02_angle-dependence\fittingOutput\angleDependence\c-free\fitparams_summary.csv',
    r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-2\220307\angle\fittingOutput\angleDependence\c-free\fitparams_summary.csv',
    ]

sampleNames = ['MA2959-2-D1', 'MA2959-2-D1-7T']

fit_comps_to_show = ['xyz']
analysisMode = 1 # Analysis mode used in stfmrAnglePlot for fitting of angle-dependence

opDir = r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\ProjE-synthesis\220307_msn001_mag-demag'

# Project F
# ipFileDirNames = [
#     r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2960-2\220202\D1_0deg\003_angle-dependence\fittingOutput\angleDependence\c-free\fitparams_summary.csv',
#     r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2960-2\220203\D4_90deg\002_angle-dependence\pos_field\fittingOutput\angleDependence\c-free\fitparams_summary.csv',
#     r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA3273-1\211209\1_angle-dependence\fittingOutput\angleDependence\c-free\fitparams_summary.csv'
#     ]

# sampleNames = [
#     'MA2960-2-D1', 
#     'MA2960-2-D4', 
#     'MA3273-1-D1'
#     ]

# Project G
# ipFileDirNames = [
    # r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-2\220131\D1_0deg\02_angle-dependence\fittingOutput\angleDependence\c-free\fitparams_summary.csv',
    # r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2960-2\220202\D1_0deg\003_angle-dependence\fittingOutput\angleDependence\c-free\fitparams_summary.csv',
    # r'D:\owncloud\1_Mn3SnN project\ST FMR Other samples\12 nm Mn3Sn\A5_0 deg\Rotation\fittingOutput\angleDependence\c-free\fitparams_summary.csv',
    # r'D:\owncloud\1_Mn3SnN project\ST FMR Other samples\10 nm Mn3Pt 001\Rotation_C5_45deg\fittingOutput\angleDependence\c-free\fitparams_summary.csv',
    # r'D:\owncloud\1_Mn3SnN project\ST FMR Other samples\10 nm Mn3Pt 111\Rotation_C5_45deg\fittingOutput\angleDependence\c-free\fitparams_summary.csv',
    # r'D:\owncloud\1_Mn3SnN project\ST FMR Other samples\12 nm Mn3Sn Cu Py\Rotation A5_0 deg\fittingOutput\angleDependence\c-free\fitparams_summary.csv',
    # r'D:\owncloud\1_Mn3SnN project\ST FMR Other samples\5 nm Pt Py\Rotation\fittingOutput\angleDependence\c-free\fitparams_summary.csv',
    # ]

# sampleNames = [
    # 'MA2959-D1',
    # 'MA2960-D1',
    # 'MA278-A5', 
    # 'MA2483-C5',
    # 'MA2475-C5',
    # 'MA1486-A5',
    # 'MA192-BL2',
    # ]

# fit_comps_to_show = ['xyz']
# analysisMode = 1 # Analysis mode used in stfmrAnglePlot for fitting of angle-dependence


# opDir = r'D:\owncloud\1_Mn3SnN project\ST FMR Other samples\Analysis'

# ____________________________________________________________________________
# MODULES
from files import File 
from helpers.pandas_helpers import add_constant_column, shift_column_to_start
from plots import GenPlot
import pandas as pd
# import matplotlib.pyplot as plt

# ____________________________________________________________________________
# DEFINITIONS
class fittingOutput:
    def __init__(self, ipFile):
        self.ipFile = ipFile
        self.read()
        
    def read(self):
        self.data = pd.read_csv(self.ipFile.fileDirName)
        
    def add_sampleName_to_data(self, name):
        self.data = add_constant_column(self.data, 'sample', name)
        self.data = shift_column_to_start(self.data, 'sample')
            
# ____________________________________________________________________________
# CODE

ipFiles = []
for ipFileDirName in ipFileDirNames:
    ipFiles.append(File(ipFileDirName))
    
fitSummary = pd.DataFrame()
for i, ipFile in enumerate(ipFiles):
    sampleName = sampleNames[i]
    fitOp = fittingOutput(ipFile)
    fitOp.add_sampleName_to_data(sampleName)
    fitSummary = pd.concat([fitSummary, fitOp.data], ignore_index=True)
    
plots = []
for i, row in fitSummary.iterrows():
    if row['fit_comps'] in fit_comps_to_show:
        if analysisMode == 1:
            tau_plt_cols =  ['xAD_norm', 'xFL_norm', 'yAD_norm', 'yFL_norm', 'zAD_norm', 'zFL_norm']
            analysisModeText = 'c-free'
            ylabel = '$\\tau^{norm}$'
        if analysisMode == 2:
            tau_plt_cols =  ['tau_xAD', 'tau_xFL', 'tau_yAD', 'tau_yFL', 'tau_zAD', 'tau_zFL']
            analysisModeText = 'full-quantitative'
            ylabel = '$\tau$ (muT/rad)'
        taus = row[tau_plt_cols]
        
        title = 'Sample: {}, Torques: {}\nAnalysis mode: {}'.format(row['sample'], row['fit_comps'], analysisModeText)
        save_label = '{}_{}'.format(row['sample'], row['fit_comps'])
        plot = GenPlot(title=title, ylabel=ylabel, save_label=save_label)
        plot.bar(tau_plt_cols, taus, label='vals')
        plots.append(plot)
        
        
        # fig, ax = plt.subplots(dpi=300)
        # bar = row[tau_plt_cols].plot.bar(ax=ax)
        # bar.set_ylabel()
        # bar.set_title()
        # bar.ticklabel_format(axis='y', style='sci', useOffset=False)
        # bar.tick_params(direction='in',top=True,right=True)
        
        # bars.append(bar)
    
max_ylims = list(plots[0].ax.get_ylim())
for i in range(2):
    for plot in plots:
        if i == 0:
            ylims = plot.ax.get_ylim()
            if ylims[0] < max_ylims[0]:
                max_ylims[0] = ylims[0]
            if ylims[1] > max_ylims[1]:
                max_ylims[1] = ylims[1]
            
        if i == 1:
            plot.ax.set_ylim(max_ylims)
            plot.ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            plot.report(opDir, saveData=True)
            
            # fig.savefig(opFileFig.fileDirName, bbox_inches='tight')
        





    
    
    
    
    
    
    
    
    
    
    
    