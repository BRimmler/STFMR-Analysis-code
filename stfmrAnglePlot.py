# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 10:29:28 2021

@author: rimmler
"""

# ____________________________________________________________________________
# SETTINGS
voltageMagnitude = 'mu' # V
plotPhiMode = 1 # 0: raw angle, 1: shifted angle
flipSign = False
do_fit = True
c_free = False # True if prefactor in fit is a fit parameter.
norm_torques = False # Normalize torques
select_files_ui = False

fit_comps_list = ['y', 'xy', 'yz', 'xyz']

opFileSubdir = '/angleDependence/'
plotDpi = 600


# ____________________________________________________________________________
# CODE
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import matplotlib.pyplot as plt
from helpers.stfmrHelpers import File, File2, read_csv_Series
import numpy as np
import modules.stfmrAnglePlotFitting as apf
import helpers.stfmrAnglePlotFitHelpers as aph

inputFiles = []
ipFileLocations = []
if select_files_ui is True:
    root = tk.Tk()
    root.withdraw()
    ipFileLocationsFiles = [File(filedialog.askopenfilename(parent=root, 
                                                            title='Choose .csv file with input file locations'))]
else:
    ipFileLocationsFiles = [
        # File(r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-2\220131\D1_0deg\02_angle-dependence\fittingOutput\angleDependence\MA2959-2-D1_angleDep_input_files.csv'),
        # File(r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-2\220131\D3_45deg\01_angle-dependence\fittingOutput\angleDependence\MA2959-2-D3_angleDep_input_files.csv'),
        # File(r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2960-2\220202\D1_0deg\003_angle-dependence\fittingOutput\angleDependence\MA2960-2-D1_angleDep_input_files.csv'),
        # File(r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2960-2\220203\D4_90deg\002_angle-dependence\pos_field\fittingOutput\angleDependence\MA2960-2-D4_angleDep_input_files.csv')
        ]

for ipFileLocationsFile in ipFileLocationsFiles:
    
    ipFileLocations = read_csv_Series(ipFileLocationsFile.file_fulldir)
    ipAngleDepFittingSummaryFile = File(ipFileLocations['angle dependence fitting summary'])
    inputData = pd.read_csv(ipAngleDepFittingSummaryFile.file_fulldir,index_col=False)
    
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
    
    if do_fit == False:
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
        outputFile = File(ipAngleDepFittingSummaryFile.filedir + '/' + ipAngleDepFittingSummaryFile.filename_wo_ext + '_anglePlot.png')
        fig.savefig(outputFile.file_fulldir, bbox_inches="tight", dpi=plotDpi)
    
    else:
        if norm_torques == True:
            suffix = '_norm'
        else:
            suffix = ''
        opFileDir = ipAngleDepFittingSummaryFile.filedir + opFileSubdir
        opFileParams = File2(opFileDir, 'fitparams_summary'+suffix+'.csv')
        # opParamsSum = pd.DataFrame(columns=['fit_comps',
        #                                     'tau_xAD', 'tau_xFL', 'tau_yAD', 'tau_yFL', 'tau_zAD', 'tau_zFL',
        #                                     'Vs_r2', 'Va_r2'])
        opParamsSum = pd.DataFrame()
        for fit_comps in fit_comps_list:
            fig, ax = plt.subplots()
            ax.scatter(x, Vs/voltageDivider, label='Vs')
            ax.scatter(x, Vas/voltageDivider, label='Vas')
    
            x_plt = np.linspace(0, 360, 100)
            if c_free == True:
                Vs_fit_opt, Vs_fit_cov, Vs_fit, Vs_plt = apf.opt_Vs_ana_free(fit_comps, x, Vs, x_plt)
                Va_fit_opt, Va_fit_cov, Va_fit, Va_plt = apf.opt_Va_ana_free(fit_comps, x, Vas, x_plt)
                params = aph.get_torques(fit_comps, Vs_fit_opt, Va_fit_opt)
            else:
                cps = aph.get_cps(ipFileLocationsFile)
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
    
            if norm_torques == True:
                params_norm = aph.norm_torques(params)
                boxtext = 'Torques (norm): \n\n'
                params = params_norm
            else:
                boxtext = 'Torques: \n\n'
    
            for key in params:
                comp = key.split('_')[1]
                boxtext += comp
                if norm_torques is True:
                    boxtext += ' = {:.2f}'.format(params[key])
                else:
                    boxtext += ' = {:.1f} $\mu$T/rad'.format(params[key]*1e6)
                boxtext += '\n'
            boxtext = boxtext[:-1]
    
            props = dict(boxstyle='round', facecolor='white', alpha=0.5)
            ax.text(1.03, 1, boxtext, verticalalignment='top',
                    transform=ax.transAxes, bbox=props, fontsize=10)
    
            ax.set_title('I = {} mA, f = {} GHz, P = {} dBm \nAssumed components: {}, Normalized: {}'.format(I, f, P, fit_comps, norm_torques))
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.legend()
            ax.set_xticks(np.arange(0, 361, 60))
    
            opFileFig = File(opFileDir + 'plt_'+fit_comps+suffix+'.png')
            opFileFig.makeDirIfNotExist()
            fig.savefig(opFileFig.file_fulldir, bbox_inches="tight", dpi=plotDpi)
    
            opFileCurves = File2(opFileDir, 'curve_'+fit_comps+suffix+'.csv')
            opCurves = pd.DataFrame()
            opCurves['phi_plt (deg)'] = x_plt
            opCurves['Vs_plt (muV)'] = Vs_plt
            opCurves['Va_plt (muV)'] = Va_plt
            opCurves.to_csv(opFileCurves.fileDirName, index=False)
    
            if c_free is True:
                opParams = pd.Series({**params})
            else:
                opParams = pd.Series({**params, **sotr})
            opParams['fit_comps'] = fit_comps
            opParams['Vs_r2'] = Vs_r2
            opParams['Va_r2'] = Va_r2
    
    
            opParamsSum = opParamsSum.append(opParams, ignore_index=True)
    
        opParamsSum = opParamsSum.set_index('fit_comps')
        opParamsSum.to_csv(opFileParams.fileDirName, index=True)


















