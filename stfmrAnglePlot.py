# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 10:29:28 2021

@author: rimmler
"""

# ____________________________________________________________________________
# SETTINGS

# Data
'''
"selectFileType"
    How to select input files:
        Mode 0: Select each file seperately through UI
        Mode 1: Select file that specifies all file locations
        Mode 2: Give file locations file in code (need to know what you are doing)    
'''
selectFileType = 0

''' 
"analysisMode":
    Requirements for different modes:
        a) Lineshape analysis (frequency-dependence)
        b) AMR calibration
        c) Irf calibration
    Mode 0:
        Plotting mode. Requires only angle-dependent data
    Mode 1:
        "c-free" fitting. V_amr is a fitting parameter and Vs and Va are fitted
        simulatneously to ensure Vamr is the same for both fits. 
        Requirement: a)
    Mode 2:
        Quantitative fitting. Torques have quantitative meaning. 
        Requirements: a)-c)
'''
analysisMode = 1


voltageMagnitude = 'mu' # V
flipSign = False
fit_phi_offset = True # Only implements for c-free mode
fit_comps_list = ['y', 'yz', 'xyz'] # Select assumed torque components
norm_to = 'yFL' # Only for mode 1. Specify which torque component to normalize to.

plotPhiMode = 1 # 0: raw angle, 1: shifted angle
plotDpi = 600


# ____________________________________________________________________________
# CODE
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import matplotlib.pyplot as plt
from files import File
from plots import GenPlot, BoxText
from helpers.file_handling import read_csv_Series
import numpy as np
import modules.stfmrAnglePlotFitting as apf
from modules.stfmrAnglePlotFittingCFree import angleDepFittingCFree, get_norm_torques
import helpers.stfmrAnglePlotFitHelpers as aph
from units import rad2deg
from helpers.stfmrAnglePlotUIHelper import get_ipFileLocationsFilesFromUI


if selectFileType == 0:
    ipFileLocationsFiles = [get_ipFileLocationsFilesFromUI(analysisMode)]

elif selectFileType == 1:
    root = tk.Tk()
    root.withdraw()
    ipFileLocationsFiles = [File(filedialog.askopenfilename(parent=root, 
                                                            title='Choose .csv file with input files locations'))]

elif selectFileType == 2:
    ipFileLocationsFiles = [
        # File(r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-2\220131\D1_0deg\02_angle-dependence\fittingOutput\angleDependence\MA2959-2-D1_angleDep_input_files.csv'),
        # File(r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-2\220131\D3_45deg\01_angle-dependence\fittingOutput\angleDependence\MA2959-2-D3_angleDep_input_files.csv'),
        # File(r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2960-2\220202\D1_0deg\003_angle-dependence\fittingOutput\angleDependence\MA2960-2-D1_angleDep_input_files.csv'),
        # File(r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2960-2\220203\D4_90deg\002_angle-dependence\pos_field\fittingOutput\angleDependence\MA2960-2-D4_angleDep_input_files.csv')
        
        File(r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\Reference_Binoy\A5_0deg_demag\fittingOutput\angleDependence\refBinoy_angleDep_input_files.csv')
        ]
else:
    raise ValueError(f'Select files type "{selectFileType}" not defined')

inputFiles = []
ipFileLocations = []
for ipFileLocationsFile in ipFileLocationsFiles:
    
    ipFileLocations = read_csv_Series(ipFileLocationsFile.fileDirName)
    ipAngleDepFittingSummaryFile = File(ipFileLocations['angle dependence fitting summary'])
    inputData = pd.read_csv(ipAngleDepFittingSummaryFile.fileDirName,index_col=False)
    
    if voltageMagnitude == 'mu':
        y_label = 'V ($\mu$V)'
        voltageDivider = 1e-6
    
        if plotPhiMode == 0:
            try:
                x =  inputData['Angle (deg)']
            except:
                try:
                    x = inputData['fieldAngle (deg)']
                except:
                    raise ValueError
            x_label = '$\phi$ (deg)'
    
            Vs = inputData['Vsym (V)']
            Vas = inputData['Vas (V)']
    
        elif plotPhiMode == 1:
            x =  inputData.sort_values(by='fieldAngle (deg)')['fieldAngle (deg)']
            x_label = '$\phi$ (deg)'
    
            Vs = inputData.sort_values(by='fieldAngle (deg)')['Vsym (V)']
            Vas = inputData.sort_values(by='fieldAngle (deg)')['Vas (V)']
    
    I = float(inputData['Current (mA)'][0])
    P = float(inputData['rf Power (dBm)'][0])
    f = float(inputData['Frequency (GHz)'][0])
    
    if flipSign == True:
        Vs *= -1
        Vas *= -1
    
    # _________________________________________________________________________
    # ANALYSIS MODE 0
    if analysisMode == 0:
        # Simple data plotting without fit
        fig, ax = plt.subplots()
        ax.scatter(x, Vs, label='Vs')
        ax.scatter(x, Vas, label='Vas')
        plt.plot(x, Vs)
        plt.plot(x, Vas)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.legend()
        ax.set_xticks(np.arange(0, 361, 60))
        ax.set_title('I = {} mA, f = {} GHz, P = {} dBm'.format(I, f, P))
        outputFileSubdir = ipAngleDepFittingSummaryFile.fileDir + '/angleDependence/plot-only'
        outputFile = File(outputFileSubdir, ipAngleDepFittingSummaryFile.fileNameWOExt + '_anglePlot.png')
        outputFile.makeDirIfNotExist()
        fig.savefig(outputFile.fileDirName, bbox_inches="tight", dpi=plotDpi)
        
    # _________________________________________________________________________
    # ANALYSIS MODE 1
    elif analysisMode == 1:
        ''' c-free fitting '''
        opFileDir = ipAngleDepFittingSummaryFile.fileDir + '/angleDependence/c-free'
        opFileParams = File(opFileDir, 'fitparams_summary.csv')
        opParamsSum = pd.DataFrame()
        for fit_comps in fit_comps_list:
            title = 'I = {} mA, f = {} GHz, P = {} dBm \nAssumed components: {}'.format(I, f, P, fit_comps)
            phiDepPlt = GenPlot(title=title, xlabel=x_label, ylabel=y_label, dpi=plotDpi)
            phiDepPlt.ax.set_xticks(np.arange(0, 361, 60))
            phiDepPlt.scatter(x, Vs/voltageDivider, label='Vs_data')
            phiDepPlt.scatter(x, Vas/voltageDivider, label='Va_data')
    
            x_plt = np.linspace(0, 360, 100)
            
            cps = aph.get_cps(1, ipFileLocationsFile)
            fitting_output = angleDepFittingCFree(x, x_plt, Vs, Vas, cps, fit_comps, fit_phi_offset, do_check_fit=False)
            
            params, params_dict, Vs_fit, Vs_plt, Va_fit, Va_plt = fitting_output
            torques, torques_norm = get_norm_torques(params, norm_to)
            if not params_dict['Vamr_s'] == params_dict['Vamr_a']:
                raise # They are forced to be the same
            Vamr = params_dict['Vamr_s']
            if not params_dict['phi0_s'] == params_dict['phi0_a']:
                raise 
            Vamr = params_dict['Vamr_s']
            phi0 = rad2deg(params_dict['phi0_s'])
    
            phiDepPlt.plot(x_plt, Vs_plt/voltageDivider, label=f'Vs_fit_{fit_comps}')
            phiDepPlt.plot(x_plt, Va_plt/voltageDivider, label=f'Va_fit_{fit_comps}')
            
            box = BoxText(1.03, 1)
            box.add_text('Fitted params:')
            box.add_empty_line()
            box.add_param('Vamr', Vamr, rep='e')
            box.add_param('phi0', phi0)
            for key, param in torques.items():
                box.add_param(key, param)              
            for key, param in torques_norm.items():
                box.add_param(key, param)
                
            phiDepPlt.make_boxtext(box)
    
            opFileFig = File(opFileDir, 'plt_'+fit_comps+'.png')
            opFileFig.makeDirIfNotExist()
            phiDepPlt.report(opFileFig.fileDir, opFileFig.fileName, saveData=True)
            
            opParams = pd.Series(params_dict|torques_norm)
            opParams['fit_comps'] = fit_comps
    
            opParamsSum = opParamsSum.append(opParams, ignore_index=True)
    
        opParamsSum = opParamsSum.set_index('fit_comps')
        opParamsSum.to_csv(opFileParams.fileDirName, index=True)
    
    # _________________________________________________________________________
    # ANALYSIS MODE 0
    elif analysisMode == 2:
        ''' Quantitative fitting of angle-dependent data '''
        opFileDir = ipAngleDepFittingSummaryFile.fileDir + '/angleDependence/full-quantitative'
        opFileParams = File(opFileDir, 'fitparams_summary.csv')
        opParamsSum = pd.DataFrame()
        for fit_comps in fit_comps_list:
            fig, ax = plt.subplots()
            ax.scatter(x, Vs/voltageDivider, label='Vs')
            ax.scatter(x, Vas/voltageDivider, label='Vas')
    
            x_plt = np.linspace(0, 360, 100)
            
            cps = aph.get_cps(2, ipFileLocationsFile)
            params, Vs_fit, Vs_plt, Va_fit, Va_plt = apf.opt_V_ana_full(fit_comps, x, Vs, Vas, x_plt, cps)
            sotr = apf.get_sotr(params, cps) # spin torque ratios
    
            def calc_r2(y, y_fit):
                ss_res = np.sum((y - y_fit) ** 2) # residual sum of squares
                ss_tot = np.sum((y - np.mean(y)) ** 2) # total sum of squares
                return 1 - (ss_res / ss_tot) # r-squared (coefficient of determination)
    
            Vs_r2 = calc_r2(Vs, Vs_fit)
            Va_r2 = calc_r2(Vas, Va_fit)
    
            ax.plot(x_plt, Vs_plt/voltageDivider, label='Vs fit ('+fit_comps+', $R^2=${:.3f})'.format(Vs_r2))
            ax.plot(x_plt, Va_plt/voltageDivider, label='Vas fit ('+fit_comps+', $R^2=${:.3f})'.format(Va_r2))
    
            # if norm_torques == True:
            #     params_norm = aph.norm_torques(params)
            #     boxtext = 'Torques (norm): \n\n'
            #     params = params_norm
            # else:
            boxtext = 'Torques: \n\n'
    
            for key in params:
                comp = key.split('_')[1]
                boxtext += comp
                # if norm_torques is True:
                #     boxtext += ' = {:.2f}'.format(params[key])
                # else:
                boxtext += ' = {:.1f} $\mu$T/rad'.format(params[key]*1e6)
                boxtext += '\n'
            boxtext = boxtext[:-1]
    
            props = dict(boxstyle='round', facecolor='white', alpha=0.5)
            ax.text(1.03, 1, boxtext, verticalalignment='top',
                    transform=ax.transAxes, bbox=props, fontsize=10)
    
            ax.set_title('I = {} mA, f = {} GHz, P = {} dBm \nAssumed components: {}'.format(I, f, P, fit_comps))
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.legend()
            ax.set_xticks(np.arange(0, 361, 60))
    
            opFileFig = File(opFileDir, 'plt_'+fit_comps+'.png')
            opFileFig.makeDirIfNotExist()
            fig.savefig(opFileFig.fileDirName, bbox_inches="tight", dpi=plotDpi)
    
            opFileCurves = File(opFileDir,'curve_'+fit_comps+'.csv')
            opFileCurves.makeDirIfNotExist()
            opCurves = pd.DataFrame()
            opCurves['phi_plt (deg)'] = x_plt
            opCurves['Vs_plt (muV)'] = Vs_plt
            opCurves['Va_plt (muV)'] = Va_plt
            opCurves.to_csv(opFileCurves.fileDirName, index=False)
    
            opParams = pd.Series({**params, **sotr})
            opParams['fit_comps'] = fit_comps
            opParams['Vs_r2'] = Vs_r2
            opParams['Va_r2'] = Va_r2
    
            opParamsSum = opParamsSum.append(opParams, ignore_index=True)
    
        opParamsSum = opParamsSum.set_index('fit_comps')
        opParamsSum.to_csv(opFileParams.fileDirName, index=True)


















