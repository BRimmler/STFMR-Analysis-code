# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 12:40:28 2022

Workbook to compare analysis of angle-dependences of different samples.

Requires BRimmler/analysis-modules library (add to PATH)

@author: rimmler
"""

# ____________________________________________________________________________
# SETTINGS
ipFileDirNames = [
    # r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-2\220131\D1_0deg\02_angle-dependence\fittingOutput\angleDependence\fitparams_summary.csv',
    # r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-2\220131\D3_45deg\01_angle-dependence\fittingOutput\angleDependence\fitparams_summary.csv',
    # r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2960-2\220202\D1_0deg\003_angle-dependence\fittingOutput\angleDependence\fitparams_summary.csv',
    # r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2960-2\220203\D4_90deg\002_angle-dependence\pos_field\fittingOutput\angleDependence\fitparams_summary.csv'
    ]

sampleNames = ['MA2959-2-D1', 'MA2959-2-D3', 'MA2960-2-D1', 'MA2960-2-D4']

fit_comps_to_show = ['y', 'xy', 'yz', 'xyz']


opDir = r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\ProjE-synthesis\220204'

# ____________________________________________________________________________
# MODULES
from files import File 
import helpers as hlp
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
        self.data = hlp.add_constant_column(self.data, 'sample', name)
        self.data = hlp.shift_column_to_start(self.data, 'sample')
            
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
        tau_plt_cols =  ['tau_xAD', 'tau_xFL', 'tau_yAD', 'tau_yFL', 'tau_zAD', 'tau_zFL']
        taus = row[tau_plt_cols]
        
        title = 'Sample: {}, Torques: {}'.format(row['sample'], row['fit_comps'])
        ylabel = 'Val (muT/rad)'
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
        





    
    
    
    
    
    
    
    
    
    
    
    