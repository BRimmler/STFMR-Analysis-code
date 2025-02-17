# -*- coding: utf-8 -*-
'''
Analysis module for analysis of angle-dependence

Author: 
    Berthold Rimmler, 
    Max Planck Institute of Microstructure Physics, Halle
    Weinberg 2
    06120 Halle
    berthold.rimmler@mpi-halle.mpg.de

'''

''' Input zone '''
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
selectFileType = 2

''' 
"analysisMode":
    Requirements for different modes:
        a) Lineshape analysis (frequency-dependence)
        b) AMR calibration
        c) Irf calibration
        d) PHE and AHE calibration
    
    Mode 0:
        Plotting mode. Requires only angle-dependent data
    Mode 1:
        "c-free" fitting. V_amr is a fitting parameter and Vs and Va are fitted
        simulatneously to ensure Vamr is the same for both fits. 
        Requirement: a)
    Mode 2:
        Quantitative fitting. Torques have quantitative meaning. 
        Requirements: a)-c)
    Mode 3:
        Semi-quantitative fitting with generalized Karimeddiny artifact description.
        Requirements: a)-c)
    Mode 4:
        Semi-quantitative fitting with generalized Karimeddiny artifact descirption 
        in XX and XY direction.
        Requirements: a)-d)
'''
analysisMode = 4

'''
"Vset_mode":
    Only for analysisMode 4.
    Specify which data to use for fitting.
    0: Vsxx, Vaxx, Vsxy
    1: Vsxx, Vaxx, Vaxy
'''
Vset_mode = 0
    


voltageMagnitude = 'mu' # V
flipSign = False
fit_phi_offset = False # Only implements for c-free mode
fit_comps_list = ['xyz'] # Select assumed torque components
assume_arts = True
norm_to = 'yFL' # Only for mode 1. Specify which torque component to normalize to.

plotPhiMode = 1 # 0: raw angle, 1: shifted angle
delta_phi = 45 # distance between angle tick values (deg)
plotDpi = 600

saveData = True

''' Input zone ends here. '''
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
from modules.stfmrKarimeddinyFitting import V_Karimeddiny_fitting, get_norm_torques_karimed, calc_Ru
from modules.stfmrKarimeddinyHallFitting import V_Karimeddiny_Hall_fitting, get_norm_torques_karimed, calc_Ru
import stfmrHelpers.stfmrAnglePlotFitHelpers as aph
from units import rad2deg
from stfmrHelpers.stfmrAnglePlotUIHelper import get_ipFileLocationsFilesFromUI


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
        
        File(r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-2\220131\D1_0deg\02_angle-dependence\fittingOutput\angleDependence\MA2959-2-D1_angleDep_input_files.csv')
        ]
else:
    raise ValueError(f'Select files type "{selectFileType}" not defined')

inputFiles = []
ipFileLocations = []
for ipFileLocationsFile in ipFileLocationsFiles:
    
    # Get input file locations
    ipFileLocations = read_csv_Series(ipFileLocationsFile.fileDirName)
    ipAngleDepFittingSummaryFile = File(ipFileLocations['angle dependence fitting summary'])
    # Get input data
    inputData = pd.read_csv(ipAngleDepFittingSummaryFile.fileDirName,index_col=False)
    if analysisMode == 4:
        # Get additional data from XY measurement
        ipAngleDepFittingXYSummaryFile = File(ipFileLocations['angle dependence fitting summary transversal'])
        inputDataXY = pd.read_csv(ipAngleDepFittingXYSummaryFile.fileDirName,index_col=False)
    
    # Extract important collumns
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
            if analysisMode == 4:
                Vsxx = Vs
                Vaxx = Vas
                Vsxy = inputDataXY['Vsym (V)']
                Vaxy = inputDataXY['Vas (V)']
    
        elif plotPhiMode == 1:
            x =  inputData.sort_values(by='fieldAngle (deg)')['fieldAngle (deg)']
            x_label = '$\phi$ (deg)'
    
            Vs = inputData.sort_values(by='fieldAngle (deg)')['Vsym (V)']
            Vas = inputData.sort_values(by='fieldAngle (deg)')['Vas (V)']
    
    # Extract fixed parameters
    I = float(inputData['Current (mA)'][0])
    P = float(inputData['rf Power (dBm)'][0])
    f = float(inputData['Frequency (GHz)'][0])
    
    # Flip sign if defined
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
        ax.set_xticks(np.arange(0, 361, delta_phi))
        ax.set_title('I = {} mA, f = {} GHz, P = {} dBm'.format(I, f, P))
        outputFileSubdir = ipAngleDepFittingSummaryFile.fileDir + '/angleDependence/plot-only'
        outputFile = File(outputFileSubdir, ipAngleDepFittingSummaryFile.fileNameWOExt + '_anglePlot.png')
        outputFile.makeDirIfNotExist()
        if saveData is True:
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
            phiDepPlt.ax.set_xticks(np.arange(0, 361, delta_phi))
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
        if saveData is True:
            opParamsSum.to_csv(opFileParams.fileDirName, index=True)
    
    # _________________________________________________________________________
    # ANALYSIS MODE 2
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
            ax.set_xticks(np.arange(0, 361, delta_phi))
    
            opFileFig = File(opFileDir, 'plt_'+fit_comps+'.png')
            opFileFig.makeDirIfNotExist()
            
            if saveData is True:
                fig.savefig(opFileFig.fileDirName, bbox_inches="tight", dpi=plotDpi)
    
            opFileCurves = File(opFileDir,'curve_'+fit_comps+'.csv')
            opFileCurves.makeDirIfNotExist()
            opCurves = pd.DataFrame()
            opCurves['phi_plt (deg)'] = x_plt
            opCurves['Vs_plt (muV)'] = Vs_plt
            opCurves['Va_plt (muV)'] = Va_plt
            
            if saveData is True:
                opCurves.to_csv(opFileCurves.fileDirName, index=False)
    
            opParams = pd.Series({**params, **sotr})
            opParams['fit_comps'] = fit_comps
            opParams['Vs_r2'] = Vs_r2
            opParams['Va_r2'] = Va_r2
    
            opParamsSum = opParamsSum.append(opParams, ignore_index=True)
    
        opParamsSum = opParamsSum.set_index('fit_comps')
        
        if saveData is True:
            opParamsSum.to_csv(opFileParams.fileDirName, index=True)

    # _________________________________________________________________________
    # ANALYSIS MODE 3
    elif analysisMode == 3:
        ''' Semi-quantitative fitting with generalized Karimeddiny artifact description '''
        
        opFileDir = ipAngleDepFittingSummaryFile.fileDir + '/angleDependence/karimeddiny'
        opFileParams = File(opFileDir, 'fitparams_summary.csv')
        opParamsSum = pd.DataFrame()
        for fit_comps in fit_comps_list:
            title = 'I = {} mA, f = {} GHz, P = {} dBm \nAssumed components: {}, assume artifacts: {}'.format(I, f, P, fit_comps, assume_arts)
            phiDepPlt = GenPlot(title=title, xlabel=x_label, ylabel=y_label, dpi=plotDpi)
            phiDepPlt.ax.set_xticks(np.arange(0, 361, delta_phi))
            phiDepPlt.scatter(x, Vs/voltageDivider, label='Vs_data')
            phiDepPlt.scatter(x, Vas/voltageDivider, label='Va_data')
            
            x_plt = np.linspace(0, 360, 100)
            
            # Get constant parameters
            cps = aph.get_cps(3, ipFileLocationsFile)
            
            # Fit
            fitting_output = V_Karimeddiny_fitting(fit_comps, x, Vs, Vas, x_plt, cps, assume_arts=assume_arts)
            params, params_dict, Vs_fit, Vs_plt, Va_fit, Va_plt = fitting_output
            torques, torques_norm = get_norm_torques_karimed(params, norm_to)
            
            # Fit quality:
            Ru_s = calc_Ru(Vs, Vs_fit)
            Ru_a = calc_Ru(Vas, Va_fit)
            
            # Plot
            phiDepPlt.plot(x_plt, Vs_plt/voltageDivider, label=f'Vs_fit_{fit_comps}')
            phiDepPlt.plot(x_plt, Va_plt/voltageDivider, label=f'Va_fit_{fit_comps}')
            
            box = BoxText(1.03, 1)
            box.add_text('Fitted params:')
            box.add_empty_line()
            box.add_param('Tart', params_dict['Tart'], rep='e')
            
            for key, param in torques.items():
                box.add_param(key, param)              
            for key, param in torques_norm.items():
                box.add_param(key, param)
                
            box.add_empty_line()
            box.add_param('Ru_s', Ru_s*100, unit=' %', rep='f')
            box.add_param('Ru_a', Ru_a*100, unit=' %', rep='f')
                
            phiDepPlt.make_boxtext(box)
    
            opFileFig = File(opFileDir, f'plt_{fit_comps}_arts={assume_arts}.png')
            opFileFig.makeDirIfNotExist()
            phiDepPlt.report(opFileFig.fileDir, opFileFig.fileName, saveData=True)
            
            opParams = pd.Series(params_dict|torques_norm)
            opParams['fit_comps'] = fit_comps
            opParams['assume_arts'] = assume_arts
            opParams['Ru_s'] = Ru_s
            opParams['Ru_a'] = Ru_a
    
            opParamsSum = opParamsSum.append(opParams, ignore_index=True)
    
        opParamsSum = opParamsSum.set_index('fit_comps')
        if saveData is True:
            opParamsSum.to_csv(opFileParams.fileDirName, index=True)
            
    # _________________________________________________________________________
    # ANALYSIS MODE 4
    elif analysisMode == 3:
        ''' Semi-quantitative fitting with generalized Karimeddiny artifact description in XX and XY direction '''
        
        opFileDir = ipAngleDepFittingSummaryFile.fileDir + '/angleDependence/karimeddiny'
        opFileParams = File(opFileDir, 'fitparams_summary.csv')
        opParamsSum = pd.DataFrame()
        
        
        for fit_comps in fit_comps_list:
            title = 'I = {} mA, f = {} GHz, P = {} dBm \nAssumed components: {}, assume artifacts: {}'.format(I, f, P, fit_comps, assume_arts)
            phiDepPlt = GenPlot(mode='vstack-share-x', title=title, xlabel=x_label, dpi=plotDpi)
            phiDepPlt.ax[0].set_xticks(np.arange(0, 361, delta_phi))
            phiDepPlt.scatter(x, Vsxx/voltageDivider, axis=0, label='Vsxx_data')
            phiDepPlt.scatter(x, Vaxx/voltageDivider,  axis=0, label='Vaxx_data')
            phiDepPlt.scatter(x, Vsxy/voltageDivider, axis=1, label='Vsxy_data')
            phiDepPlt.scatter(x, Vaxy/voltageDivider,  axis=1, label='Vay_data')
            
            x_plt = np.linspace(0, 360, 100)
            
            # Get constant parameters
            cps = aph.get_cps(4, ipFileLocationsFile)
            
            # Fit
            fitting_output = V_Karimeddiny_Hall_fitting(fit_comps, x, Vset_mode, Vs, Vas, x_plt, cps, assume_arts=assume_arts)
            params, params_dict, Vsxx_fit, Vsxx_plt, Vaxx_fit, Vaxx_plt, Vsxy_fit, Vsxy_plt, Vaxy_fit, Vaxy_plt = fitting_output
            torques, torques_norm = get_norm_torques_karimed(params, norm_to)
            
            # Fit quality:
            Ru_sxx = calc_Ru(Vsxx, Vsxx_fit)
            Ru_axx = calc_Ru(Vaxx, Vaxx_fit)
            Ru_sxy = calc_Ru(Vsxy, Vsxy_fit)
            Ru_axy = calc_Ru(Vaxy, Vaxy_fit)
            
            # Plot
            phiDepPlt.plot(x_plt, Vsxx_plt/voltageDivider, axis=0, label=f'Vsxx_fit_{fit_comps}')
            phiDepPlt.plot(x_plt, Vaxx_plt/voltageDivider, axis=0, label=f'Vaxx_fit_{fit_comps}')
            phiDepPlt.plot(x_plt, Vsxy_plt/voltageDivider, axis=1, label=f'Vsxy_fit_{fit_comps}')
            phiDepPlt.plot(x_plt, Vaxy_plt/voltageDivider, axis=1, label=f'Vaxy_fit_{fit_comps}')
            
            box = BoxText(1.03, 1)
            box.add_text('Fitted params:')
            box.add_empty_line()
            box.add_param('Tart', params_dict['Tart'], rep='e')
            
            for key, param in torques.items():
                box.add_param(key, param)              
            for key, param in torques_norm.items():
                box.add_param(key, param)
                
            box.add_empty_line()
            box.add_param('Ru_sxx', Ru_sxx*100, unit=' %', rep='f')
            box.add_param('Ru_axx', Ru_axx*100, unit=' %', rep='f')
            box.add_param('Ru_sxy', Ru_sxy*100, unit=' %', rep='f')
            box.add_param('Ru_axy', Ru_axy*100, unit=' %', rep='f')
                
            phiDepPlt.make_boxtext(box)
    
            opFileFig = File(opFileDir, f'plt_{fit_comps}_arts={assume_arts}.png')
            opFileFig.makeDirIfNotExist()
            phiDepPlt.report(opFileFig.fileDir, opFileFig.fileName, saveData=True)
            
            opParams = pd.Series(params_dict|torques_norm)
            opParams['fit_comps'] = fit_comps
            opParams['assume_arts'] = assume_arts
            opParams['Vset_mode'] = Vset_mode
            opParams['Ru_sxx'] = Ru_sxx
            opParams['Ru_axx'] = Ru_axx
            opParams['Ru_sxy'] = Ru_sxy
            opParams['Ru_axy'] = Ru_axy
    
            opParamsSum = opParamsSum.append(opParams, ignore_index=True)
    
        opParamsSum = opParamsSum.set_index('fit_comps')
        if saveData is True:
            opParamsSum.to_csv(opFileParams.fileDirName, index=True)


























