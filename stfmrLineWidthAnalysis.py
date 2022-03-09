# -*- coding: utf-8 -*-
'''
Analysis module for analysis of dc-bias dependence ("line width analysis")

Author: 
    Berthold Rimmler, 
    Max Planck Institute of Microstructure Physics, Halle
    Weinberg 2
    06120 Halle
    berthold.rimmler@mpi-halle.mpg.de

'''

''' Input zone '''
# ____________________________________________________________________________
# INPUT

# Data
selectInputUI = False       # True to select file through UI. Otherwise specify in code below.

currentLimit = 1.           # Maximum dc current used for linear fit in mA

saveOutput = True
plotDpi = 600               # Plot resolution [default: 600]


''' Input zone ends here. '''
# ____________________________________________________________________________
# CODE
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import numpy as np
from helpers.file_handling import read_csv_Series
from scipy.optimize import curve_fit
from files import File
from plots import GenPlot, BoxText
from scipy import constants
import outputs as op

# Get data and sample
if selectInputUI is True:
    root = tk.Tk()
    ipFileLocDirName = filedialog.askopenfilename(parent=root, title='Choose .csv file with input file locations')
    root.destroy()
else:
    ipFileLocDirName = r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-2\220307\dcbias\ana\fittingOutput\lineWidthAna\linewidth_input_files.csv'


ipFileLocationsFile = File(ipFileLocDirName)
ipFileLocations = read_csv_Series(ipFileLocationsFile.fileDirName)

# Get fitting output of dc bias measurements and lineshape analysis output
ipFittingOutputFile = File(ipFileLocations['linewidth raw fitting output'])
ipLineshapeAnaOutputFile = File(ipFileLocations['lineshape analysis output'])
ipLineshapeFittingOutputFile = File(ipFileLocations['lineshape raw fitting output'])
ipResistivitiesFile = File(ipFileLocations['resistivities'])
ipDeviceDimsFile = File(ipFileLocations['device dimensions'])


ipFittingOutput = pd.read_csv(ipFittingOutputFile.fileDirName, index_col='Index')

# Get frequency and angle
f_data = ipFittingOutput['Frequency (GHz)']
if len(f_data.unique()) > 1:
    raise ValueError('Frenquency is not constant.')
f = f_data.unique()[0]*1e9 # to Hz

phi_data = ipFittingOutput['fieldAngle (deg)']
if len(phi_data.unique()) > 1:
    raise ValueError('fieldAngle is not constant.')

C = {}
C['f (Hz)'] = f
C['phi (deg)'] = phi_data.unique()[0]

# Limit to currents below currentLimit
ipFittingOutput = ipFittingOutput[abs(ipFittingOutput['Current (mA)']) <= currentLimit]

# Split in positive and negative field
fittingData = {'pos': ipFittingOutput[ipFittingOutput['Hres (Oe)']>0],
               'neg': ipFittingOutput[ipFittingOutput['Hres (Oe)']<0]}

# Get resonance fields
for key, data in fittingData.items():
    C[f'H0_{key} (Oe)'] = data['Hres (Oe)'].mean()
    C[f'H0Error_{key} (Oe)'] = data['HresError (Oe)'].mean()
    
    

ipLineshapeFittingOutput = pd.read_csv(ipLineshapeFittingOutputFile.fileDirName, index_col='Index')

ipLineshapeAnaOutput = read_csv_Series(ipLineshapeAnaOutputFile.fileDirName)
C = {**C, **ipLineshapeAnaOutput}

resistivities = read_csv_Series(ipResistivitiesFile.fileDirName).to_dict()
C = {**C, **resistivities}

device_dimensions = read_csv_Series(ipDeviceDimsFile.fileDirName).to_dict()
C = {**C, **device_dimensions}


# Linear fit dDelta over dIdc
def fDelta(Idc, m, b):
    return m * Idc + b

m_Delta_linfit = {}
for key in ['neg', 'pos']:
    popt, pcov = curve_fit(fDelta, fittingData[key]['Current (mA)'], fittingData[key]['DeltaH (Oe)'])
    perr = np.sqrt(np.diag(pcov))
    m_Delta_cgs = popt[0] # Oe/mA
    m_Delta_err_cgs = perr[0] # Oe/mA
    C[f'm_Delta_{key} (Oe/mA)'] = m_Delta_cgs
    C[f'm_Delta_{key}_err (Oe/mA)'] = m_Delta_err_cgs
    m_Delta_linfit[key] = {'m': popt[0], 'b': popt[1]}
    

# Clean unit --> Convert everything to SI
C['Meffopt (A/m)'] = C['Meffopt (emu/cm3)'] * 1e3
C['MeffoptError (A/m)'] = C['MeffoptError (emu/cm3)'] * 1e3
C['Ms (A/m)'] = C['Ms (emu/cm3)'] * 1e3
C['MsError (A/m)'] = C['MsError (emu/cm3)'] * 1e3

for key in ['neg', 'pos']:
    C[f'H0_{key} (A/m)'] = C[f'H0_{key} (Oe)'] / (4*np.pi * 1e-3)
    C[f'H0Error_{key} (A/m)'] = C[f'H0Error_{key} (Oe)'] / (4*np.pi * 1e-3)
    C[f'm_Delta_{key} (1/m)'] = C[f'm_Delta_{key} (Oe/mA)'] / (4*np.pi*1e-6)
    C[f'm_Delta_{key}_err (1/m)'] = C[f'm_Delta_{key}_err (Oe/mA)'] / (4*np.pi*1e-6)


# Calculate SHA

# m_halpha from parameters from previous fitting and input parameters
def calc_m_alpha(f, g, phi, Ms, t, H0, Meff, MsErr, tErr, H0Err, MeffErr):
    hbar = constants.hbar
    e = constants.e
    mu_0 = constants.mu_0
    c1 = 2*np.pi*f/g*hbar/(2*e)
    sin_phi = np.sin(phi*np.pi/180)
    c2 = mu_0*Ms*t*(H0 + 0.5*Meff)
    c2_err = mu_0*np.sqrt(
        (MsErr*t*(H0 + 0.5*Meff))**2
        +(Ms*tErr*(H0 + 0.5*Meff))**2
        +(Ms*t**H0Err)**2
        +(Ms*t*0.5*MeffErr)**2)
    m_alpha = c1 * sin_phi / c2
    m_alpha_err = c1 * sin_phi * c2_err / (c2**2)
    return m_alpha, m_alpha_err

for key in ['neg', 'pos']:
    Ms = C['Ms (A/m)']
    Meff = C['Meffopt (A/m)']
    if key == 'neg':
        Ms *= -1    # Ms and Meff become negative (matters for term in bracket)
        Meff *= -1
        
    m_alpha, m_alpha_err = calc_m_alpha(C['f (Hz)'], C['g'], C['phi (deg)'], Ms, C['t (m)'],
                                    C[f'H0_{key} (A/m)'], Meff, C['MsError (A/m)'], 
                                    C['tError (m)'], C[f'H0Error_{key} (A/m)'], C['MeffoptError (A/m)'])
    
    C[f'm_alpha_{key}'] = m_alpha
    C[f'm_alpha_{key}_err'] = m_alpha_err

# Resistance ratio from resistivities input file
def calc_RR(rho_FM, rho_SHM, rho_FM_err, rho_SHM_err):
    rr = 1 + rho_SHM / rho_FM
    rr_err = np.sqrt((rho_SHM_err/rho_FM)**2+(rho_SHM*rho_FM_err/rho_FM**2)**2)
    return rr, rr_err

rr, rr_err = calc_RR(C['rho_fm (muOhmcm)'], C['rho_shm (muOhmcm)'], 
                     C['rho_fm_err (muOhmcm)'], C['rho_shm_err (muOhmcm)'])

C['rr'] = rr
C['rr_err'] = rr_err

# Device cross section from device_dimensions input file
def calc_A_dev(w, h, wErr, hErr):
    ''' device width w in um and error, device height h in nm and error '''
    w *= 1e-6 # in m
    h *= 1e-9 # in m
    wErr *= 1e-6
    hErr *= 1e-9
    A_dev = w*h # m2
    A_dev_err = np.sqrt((wErr*h)**2 + (w*hErr)**2)
    return A_dev, A_dev_err

A_dev, A_dev_err = calc_A_dev(C['device_width (um)'], C['device_height (nm)'], 
                        C['device_width_error (um)'], C['device_height_error (nm)'])

C['A_dev'] = A_dev
C['A_dev_err'] = A_dev_err


# Spin-Hall angle
def calc_SHA(m_Delta, m_alpha, rr, A_dev, m_Delta_err, m_alpha_err, rr_err, A_dev_err):
    SHA = m_Delta / m_alpha * rr * A_dev
    SHA_err = np.sqrt(
        (m_Delta_err / m_alpha * rr * A_dev)**2 + (m_Delta *m_alpha_err / m_alpha**2 * rr * A_dev)**2
        + (m_Delta / m_alpha * rr_err * A_dev)**2 + (m_Delta / m_alpha * rr * A_dev_err)**2)
    return SHA, SHA_err

SHA = {}
SHA_err = {}
for key in ['neg', 'pos']:
    sha, sha_err = calc_SHA(C[f'm_Delta_{key} (1/m)'], C[f'm_alpha_{key}'], C['rr'], C['A_dev'],
                            C[f'm_Delta_{key}_err (1/m)'], C[f'm_alpha_{key}_err'], C['rr_err'], C['A_dev_err'])
    SHA[key] = sha
    SHA_err[key] = sha_err
    C[f'SHA_{key}'] = sha
    C[f'SHA_{key}_err'] = sha_err
    

# print('SHA = {:.2e} +- {:.2e}'.format(SHA, SHA_err))
    

# Plotting
lwPlot = GenPlot(xlabel='$I_{dc}$ (mA)', ylabel='$\Delta$ (Oe)', dpi=plotDpi)
lwBox = BoxText(1.025, 1)
for key in ['neg', 'pos']:
    # Data
    I_dc = fittingData[key]['Current (mA)']
    lwPlot.errorbar_scatter(I_dc, fittingData[key]['DeltaH (Oe)'],
                            yerr=fittingData[key]['DeltaHError (Oe)'], label=f'{key}: data')
    
    # Linear fit
    lwPlot.plot(I_dc, fDelta(I_dc, m_Delta_linfit[key]['m'], m_Delta_linfit[key]['b']), label=f'{key}: lin fit')
    
    lwBox.add_param(f'SHA_{key}', SHA[key], rep='e', error=SHA_err[key])
lwPlot.make_boxtext(lwBox)

if saveOutput is True:
    # Plot
    opDir = ipFileLocationsFile.fileDir + '/output'
    lwPlot_opFile = File(opDir, 'lineWidthAnaFitting.png' )
    lwPlot.report(lwPlot_opFile.fileDir, opName=lwPlot_opFile.fileNameWOExt, saveData=True)
    
    # Parameter
    params_output = op.SeriesOutput(C, opDir, 'lineWidthAnaParams.csv')
    params_output.save()
    






















