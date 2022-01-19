# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 14:43:54 2021

@author: rimmler
"""

#_____________________________________________________________________________
# INPUT
corrOffset = True # Correct for measurement offset between RF and DC measurement

Pmode = 'dBm' # power in mW or power level in dBm for plot
plotDpi = 300 # Resolution


#_____________________________________________________________________________
# CODE

import tkinter as tk
from tkinter import filedialog
import pandas as pd
import matplotlib.pyplot as plt
from stfmrHelpers import File, File2
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress
from datetime import datetime


#_____________________________________________________________________________
# Data

root = tk.Tk()
root.withdraw()
inputFileDC = File(filedialog.askopenfilename(parent=root, title='Idc')).to_File2()
inputFileRF = File(filedialog.askopenfilename(parent=root, title='Irf')).to_File2()

inputFileDCs = [inputFileDC]
inputFileRFs = [inputFileRF]

# inputFileDCs = [
#     File(r'D:\owncloud\0_Personal\DATA\Mn3SnN\ST-FMR\MA2959-2\211126/001_Idc_1mA_15mA.txt').to_File2(),
#     File(r'D:\owncloud\0_Personal\DATA\Mn3SnN\ST-FMR\MA2960-2\211126/007_Idc_1mA_15mA_30p.txt').to_File2(),
#     File(r'D:\owncloud\0_Personal\DATA\Mn3SnN\ST-FMR\MA1487\211126/001_MA1487-left-B4_Idc_1mA_15mA.txt').to_File2()
#     ]
# inputFileRFs = [
#     File(r'D:\owncloud\0_Personal\DATA\Mn3SnN\ST-FMR\MA2959-2\211126/003_Irf_3dBm_15dBm_DC_1mA_above_10dBm.txt').to_File2(),
#     File(r'D:\owncloud\0_Personal\DATA\Mn3SnN\ST-FMR\MA2960-2\211126/008_Irf_7dBm_9dBm_1mA_9GHz_above_10dBm.txt').to_File2(),
#     File(r'D:\owncloud\0_Personal\DATA\Mn3SnN\ST-FMR\MA1487\211126/002_MA1487-left-B4_Irf_7dBm_9dBm_1mA_9GHz_above_10dBm.txt').to_File2()
#     ]




#_____________________________________________________________________________
# Functions

def get_dc_data(ipFile):
    dcData = pd.read_csv(ipFile.fileDirName)
    Idc = dcData['I(A)'].to_numpy()
    Rdc = dcData[' R(ohm)'].to_numpy()
    return dcData, Idc, Rdc

def get_rf_data(ipFile):
    rfData = pd.read_csv(ipFile.fileDirName)
    # print(rfData)
    LPrf = rfData['Power(dBm)'].to_numpy()
    Rrf = rfData[' R(ohm)'].to_numpy()
    Imeas = rfData[' I(A)'].to_numpy()[0]
    frf = rfData[' f(Hz)'].to_numpy()[0]
    return rfData, LPrf, Rrf, Imeas, frf

def corr_Rrf_offset(Prf_mW, Rrf, Imeas, a_opt, b_opt):
    ''' There may be an arbitrary resistance offset between the DC and RF
    measurements. This is subtracted by extrapolating Rrf as f(P) to P=0,
    where the real resistance should overlap with that of Rdc(I=Imeas) '''
    slope, intercept, r, p, se = linregress(Prf_mW, Rrf)
    Rmeas = fit_Rdc(Imeas, a_opt, b_opt)
    offset = intercept - Rmeas
    Rrf_offset = Rrf - offset

    ''' In addition, the data must be corrected by the differential resistance
    caused by the DC measurmeent current. Note that b_opt corresponds to the
    theoretical resistance at 0 DC current'''
    Rrf_corr = Rrf_offset - (Rmeas-b_opt)
    return Rrf_corr

def conv_dBm_mW(LP):
    ''' Convert Leistungspegel LP in dBm to power in mW '''
    return 10**(LP/10) # mW

def fit_Rdc(Idc, a, b):
    return a * Idc**2 + b

def fit_Idc(Rdc, a, b):
    ''' Fit function for measurement Rdc vs. Idc '''
    return np.sqrt((Rdc - b) / a)

def fit_Irf(P_sqrt, m):
    ''' P_sqrt = sqrt(P) '''
    return m * P_sqrt

def calc_Idc_equiv_rf(Rrf, a, b):
    ''' The equivalent dc current that would cause the same Joule heating as
    the applied rf current '''
    return np.sqrt(abs(Rrf - b) / a)

def calc_Irf_peak(Idc_equiv_rf, Imeas):
    ''' Peak rf current from equivalent dc current '''
    # if subImeas == True:
    #     Irf_peak = np.sqrt(2) * (Idc_equiv_rf - Imeas)
    # elif subImeas == False:
    Irf_peak = np.sqrt(2) * (Idc_equiv_rf)
    return Irf_peak

def dict2df(d):
    ''' Converts a dictionary with different lengths to a pd.DataFrame '''
    return pd.DataFrame.from_dict(d, orient='index').transpose()

#_____________________________________________________________________________
# Main

for i in range(len(inputFileDCs)):
    # Get data
    inputFileDC = inputFileDCs[i]
    inputFileRF = inputFileRFs[i]

    dcData, Idc, Rdc = get_dc_data(inputFileDC)
    rfData, PLrf_dBm, Rrf, Imeas, frf = get_rf_data(inputFileRF)

    Prf_mW = conv_dBm_mW(PLrf_dBm)

    # Fit Rdc vs. Idc to fit function to get a, b
    # Idc_fit_opt, Idc_fit_cov = curve_fit(fit_Idc, Idc, Rdc)
    Rdc_fit_opt, Idc_fit_cov = curve_fit(fit_Rdc, Idc, Rdc)

    a_opt = Rdc_fit_opt[0]
    b_opt = Rdc_fit_opt[1]


    # Get Rrf measurement offset and subtract from Rrf data
    if corrOffset is True:
        Rrf = corr_Rrf_offset(Prf_mW, Rrf, Imeas, a_opt, b_opt)

    # Calculate the dc current equivalent to applied rf current
    Idc_equiv_rf = calc_Idc_equiv_rf(Rrf, a_opt, b_opt)

    # Calculate resulting rf peak current
    Irf_peak = calc_Irf_peak(Idc_equiv_rf, Imeas)

    # Fit rf peak current vs. sqrt(rf power in mW)
    Prf_mW_sqrt = np.sqrt(Prf_mW)
    PLrf_dBm_sqrt = np.sqrt(PLrf_dBm)
    Irf_fit_opt, Irf_fit_cov = curve_fit(fit_Irf, Prf_mW_sqrt, Irf_peak)

    m_opt = Irf_fit_opt[0] # unit depends on Pmode


    #_____________________________________________________________________________
    # Plots
    # Rdc vs. Idc and Rrf vs. P: Data and fit
    Idc_plt = np.linspace(Idc[0], Idc[-1], 100)
    Rdc_plt = fit_Rdc(Idc_plt, a_opt, b_opt)

    Prf_mw_sqrt_plt = np.linspace(Prf_mW_sqrt[0], Prf_mW_sqrt[-1], 100)
    PLrf_dBm_sqrt_plt = np.linspace(PLrf_dBm_sqrt[0], PLrf_dBm_sqrt[-1], 100)

    Irf_peak_plt = fit_Irf(Prf_mw_sqrt_plt, m_opt)


    fig, axs = plt.subplots(1 ,2, figsize=(10,4))
    ax=axs[0]
    ln1=ax.scatter(Idc*1e3, Rdc, label='Rdc Data')
    ln2=ax.plot(Idc_plt*1e3, Rdc_plt, c='r', label='Rdc fit')
    ax2 = ax.twiny()
    if corrOffset is True:
        label = 'Rrf data, corrected'
    else:
        label = 'Rrf data, raw'


    if Pmode == 'mW':
        ln3=ax2.scatter(Prf_mW, Rrf, facecolors='none', edgecolors='g', label=label)
        ax2.set_xlabel('$P_{rf}$ (mW)')
        ax2.set_xlim(left=0)
    elif Pmode == 'dBm':
        ln3=ax2.scatter(PLrf_dBm, Rrf, facecolors='none', edgecolors='g', label=label)
        ax2.set_xlabel('$PL_{rf}$ (dBm)')

    ax.set_xlabel('$I_{dc}$ (mA)')
    ax.set_ylabel('$R$ ($\Omega$)')

    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=0)

    ax3=axs[1]
    ax3.scatter(Prf_mW_sqrt, Irf_peak*1e3, label='$I_{rf}^{peak}$ data')
    ax3.plot(Prf_mw_sqrt_plt, Irf_peak_plt*1e3, c='r', label='$I_{rf}^{peak}$ fit')
    ax3.set_xlabel('$(P_{rf})^{0.5}$ (mW$^{1/2}$)')
    ax3.set_ylabel('$I_{rf}^{peak}$ (mA)')
    ax3.set_xlim(left=0)
    ax3.set_ylim(bottom=0)
    ax3.legend()

    boxtext = '$I_{dc}^{meas}$'+' = {:.2f} mA\n'.format(Imeas*1e3)
    boxtext += '$I_{rf}^{peak}$ = m*$(P [mW])^{1/2}$\n'
    boxtext += 'm = {:.3f} '.format(m_opt*1e3)+'mA/mW$^{0.5}$'
    boxprops = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax3.text(0.5, 0.22, boxtext, verticalalignment='top',
            transform=ax3.transAxes, bbox=boxprops, fontsize=10)

    plt.show()

    # Irf_peak vs. sqrt(P): Data and fit


    #_____________________________________________________________________________
    # Report
    today = datetime.today().strftime('%y%m%d')
    opFileDir = inputFileDC.fileDir + f'/irf_calib_output_{today}'

    # Plot:
    opPlotFile = File2(opFileDir, 'Irf_calib_plot.png')
    opPlotFile.makeDirIfNotExist()

    fig.savefig(opPlotFile.fileDirName, bbox_inches="tight", dpi=plotDpi)

    # Curves:
    opCurvesFile = File2(opFileDir, 'Irf_calib_curves.csv')

    opCurves_dict = {'Idc (A)': Idc,
                              'Idc_plt (A)': Idc_plt,
                              'Rdc (Ohm)': Rdc,
                              'Rdc_plt (Ohm)': Rdc_plt,
                              'Prf_mW (mW)': Prf_mW,
                              'PLrf_dBm (dBm)': PLrf_dBm,
                              'Rrf (Ohm)': Rrf,
                              'Prf_mW_sqrt (sqrt(mW))': Prf_mW_sqrt,
                              'Prf_mw_sqrt_plt (sqrt(mW))': Prf_mw_sqrt_plt,
                              'Irf_peak (A)': Irf_peak,
                              'Irf_peak_plt (A)': Irf_peak_plt
                              }

    opCurves = dict2df(opCurves_dict)
    opCurves.to_csv(opCurvesFile.fileDirName, index=False)


    # Parameters:
    opParamsFile = File2(opFileDir, 'Irf_calib_fit_params.csv')
    opParamsFile.makeDirIfNotExist()

    opParams = pd.Series({''
                        'm (A/sqrt(mW))': m_opt,
                          'corrOffset': corrOffset,
                          'inputFileDC': inputFileDC.fileDirName,
                          'inputFileRF': inputFileRF.fileDirName
                          })

    opParams.to_csv(opParamsFile.fileDirName, header=False)














