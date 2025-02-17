# -*- coding: utf-8 -*-
"""
Created on Fri Aug 20 16:06:58 2021

@author: rimmler
"""

#_____________________________________________________________________________
# INPUT
'''
IP Data must be PPMS .dat file of measurement Rxx/Rxy vs. IP/OP angle

Units:
'''
sampleID = 'MA2959-2-D4'
effect = 'amr'

#_____________________________________________________________________________
# CODE

import tkinter as tk
from tkinter import filedialog
import pandas as pd
import matplotlib.pyplot as plt
from files import File
import numpy as np
from numpy import sin, cos, pi
from scipy.optimize import curve_fit


#_____________________________________________________________________________
# Data
root = tk.Tk()
root.withdraw()
ipFile = File(filedialog.askopenfilename(parent=root, title='Select input file'))
# ipFile = File(r'D:/DATA/Mn3SnN/PPMS/MA2959-1-D4/210830/003_IP_13-14_MA2959-1-D4_ave_10_Rxxphi_H_1000Oe__I_1e-3A.dat').to_File2()

#_____________________________________________________________________________
# Functions
def deg2rad(x):
    return x * pi / 180

def rad2deg(x):
    return x * 180 / pi

def fit_Ramr(phi, R0, phi0, DeltaR):
    return R0 + DeltaR * (cos(phi-phi0))**2

def calc_dRamrdphi(phi, phi0, DeltaR):
    ''' Derivative dR/dphi in Ohm/rad '''
    return -DeltaR * sin(2*(phi-phi0))

def fit_Rphe(phi, R0, phi0, DeltaR):
    return 

def calc_dRphedphi(phi, phi0, DeltaR):
    return 

def fit_Rahe(phi, R0, phi0, DeltaR):
    return

def calc_dRahedphi(phi, phi0, DeltaR):
    return

#_____________________________________________________________________________
# Main

# Get data
ipData = pd.read_csv(ipFile.fileDirName, sep='\t')
for col in ipData.columns:
    if 'Cur' in col.split('_'):
        cur_col_name = col
    elif 'Volt' in col.split('_'):
        volt_col_name = col
    elif 'Res' in col.split('_'):
        res_col_name = col

T_fix = round(np.average(ipData['Temperature (K)']))
H_fix = round(np.average(ipData['Field (Oe)']))
I_fix = np.average(ipData[cur_col_name])

phi_deg = ipData['Angle']
phi_rad = deg2rad(phi_deg)
R = ipData[res_col_name]

# Fit R(phi)
if effect in ['AMR', 'amr']:
    effect = 'AMR'
    fit_R = fit_Ramr
    calc_dRdphi = calc_dRamrdphi
    phi_label = '$\phi$'
elif effect in ['PHE', 'phe']:
    effect = 'PHE'
    fit_R = fit_Rphe
    calc_dRdphi = calc_dRphedphi
    phi_label = '$\phi$'
elif effect in ['AHE', 'ahe']:
    effect = 'AHE'
    fit_R = fit_Rahe
    calc_dRdphi = calc_dRahedphi
    phi_label = '$\theta$'

R_fit_opt, R_fit_cov = curve_fit(fit_R, phi_rad, R)
R0_fit = R_fit_opt[0]
phi0_fit = R_fit_opt[1]
phi0_fit_deg = rad2deg(phi0_fit)
DeltaR_fit = R_fit_opt[2]


#_____________________________________________________________________________
# Plots
# R and dR/dphi vs. phi
phi_plt = np.linspace(0, 360, len(phi_deg))
phi_plt_rad = deg2rad(phi_plt)
R_plt = fit_R(phi_plt_rad, R0_fit, phi0_fit, DeltaR_fit)



dRdphi_plt = calc_dRdphi(phi_plt_rad, phi0_fit, DeltaR_fit)



fig, ax = plt.subplots(dpi=600)
ax.scatter(phi_deg, R, label='R data')
ax.plot(phi_plt, R_plt, c='g', label='R fit')
ax2 = ax.twinx()
ax2.plot(phi_plt, dRdphi_plt, c='r', label='dR/dphi fit')

lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc=0)

boxtext = '$\phi_0$ = {:.0f}°\nR0 = {:.2f} $\Omega$\nDeltaR = {:.3f} $\Omega$'.format(phi0_fit_deg, R0_fit, DeltaR_fit)
boxprops = dict(boxstyle='round', facecolor='white', alpha=0.8)
ax2.text(0.7, 0.2, boxtext, verticalalignment='top',
        transform=ax.transAxes, bbox=boxprops, fontsize=10)

ax.set_title('{}, {}: $\phi$={}, T={} K, H={} Oe, I={:.1f} mA'.format(sampleID, effect, phi_label, T_fix, H_fix, I_fix*1e3))
ax.set_xlabel('$\phi$ (deg)')
ax.set_ylabel('R ($\Omega$)')
ax2.set_ylabel('dR/d$\phi$ ($\Omega$/rad)')
# ax.legend()
ax.set_xticks(np.arange(0, 361, 90))

#_____________________________________________________________________________
# Report
opFileDir = ipFile.fileDir + f'/{effect}_fitting_output'

# Plot
opFigFile = File(opFileDir, f'/{effect}_fit_plt.png')
opFigFile.makeDirIfNotExist()
fig.savefig(opFigFile.fileDirName, bbox_inches = 'tight')

# Curves: R vs. phi, fit R, dR/dphi vs. plt phi
opCurvesFile = File(opFileDir, f'/{effect}_fit_curves.csv')
opCurves = pd.DataFrame({'phi (deg)': phi_deg.to_numpy(),
                         'R_data (Ohm)': R.to_numpy(),
                         'phi_plt (deg)': phi_plt,
                         'R_fit (Ohm)': R_plt,
                         'dRdphi (Ohm/rad)': dRdphi_plt
                         })

opCurves.to_csv(opCurvesFile.fileDirName, index=False)


# Parameters: R0_fit, DeltaR_fit
opParamsFile = File(opFileDir, f'/{effect}_fit_params.csv')

opParams = pd.Series({'Input file name': ipFile.fileName,
                      'T (K)': T_fix,
                      'H (Oe)': H_fix,
                      'I (A)': I_fix,
                      'R0_fit (Ohm)': R0_fit,
                      'phi0_fit (deg)': rad2deg(phi0_fit),
                      'DeltaR_fit (Ohm)': DeltaR_fit
                      })

opParams.to_csv(opParamsFile.fileDirName)

















