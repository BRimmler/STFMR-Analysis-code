# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 15:09:08 2022

@author: rimmler
"""

# ____________________________________________________________________________
# INPUT

# Data
inputUI = False # not working, give file dir and name
ipFileLocDirName = r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-2\210805\005_linewidth_D1\fittingOutput\lineWidthAna\linewidth_input_files.csv'

plotDpi = 300

# ____________________________________________________________________________
# CODE
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from helpers.file_handling import read_csv_Series
from scipy.optimize import curve_fit
from files import File
from plots import GenPlot, BoxText
from scipy import constants

# Get input data
if inputUI is True:
    pass # Not currently working
    root = tk.Tk()
    # root.withdraw()
    ipFileLocDirName = filedialog.askopenfilename(parent=root, title='Choose .csv file with input file locations')

ipFileLocationsFile = File(ipFileLocDirName)
ipFileLocations = read_csv_Series(ipFileLocationsFile.fileDirName)

# Get fitting output of dc bias measurements and lineshape analysis output
ipFittingOutputFile = File(ipFileLocations['linewidth raw fitting output'])
ipLineshapeAnaOutputFile = File(ipFileLocations['lineshape analysis output'])
ipLineshapeFittingOutputFile = File(ipFileLocations['lineshape raw fitting output'])
ipResistivitiesFile = File(ipFileLocations['resistivities'])
ipDeviceDimsFile = File(ipFileLocations['device dimensions'])


ipFittingOutput = pd.read_csv(ipFittingOutputFile.fileDirName, index_col='Index')
f_data = ipFittingOutput['Frequency (GHz)']
if len(f_data.unique()) > 1:
    raise ValueError('Frenquency is not constant.')
f = f_data[0]*1e9 # to Hz

phi_data = ipFittingOutput['fieldAngle (deg)']
if len(phi_data.unique()) > 1:
    raise ValueError('fieldAngle is not constant.')

Const = {}
Const['f (Hz)'] = f
Const['phi (deg)'] = phi_data[0]

ipLineshapeFittingOutput = pd.read_csv(ipLineshapeFittingOutputFile.fileDirName, index_col='Index')
Const['H0 (Oe)'] = ipLineshapeFittingOutput['Hres (Oe)'][ipLineshapeFittingOutput['Frequency (GHz)']==f*1e-9].iloc[0]
Const['H0Error (Oe)'] = ipLineshapeFittingOutput['HresError (Oe)'][ipLineshapeFittingOutput['Frequency (GHz)']==f*1e-9].iloc[0]

ipLineshapeAnaOutput = read_csv_Series(ipLineshapeAnaOutputFile.fileDirName)
Const = {**Const, **ipLineshapeAnaOutput}

resistivities = read_csv_Series(ipResistivitiesFile.fileDirName).to_dict()
Const = {**Const, **resistivities}

device_dimensions = read_csv_Series(ipDeviceDimsFile.fileDirName).to_dict()
Const = {**Const, **device_dimensions}



# Linear fit dDelta over dIdc
Idc = ipFittingOutput['Current (mA)']
Delta = ipFittingOutput['DeltaH (Oe)']
def fDelta(Idc, m, b):
    return m * Idc + b
popt, pcov = curve_fit(fDelta, Idc, Delta)
perr = np.sqrt(np.diag(pcov))
m_Delta_cgs = popt[0] # Oe/mA
m_Delta_err_cgs = perr[0] # Oe/mA

# Clean unit --> Convert everything to SI
Const['Meffopt (A/m)'] = Const['Meffopt (emu/cm3)'] * 1e3
Const['MeffoptError (A/m)'] = Const['MeffoptError (emu/cm3)'] * 1e3
Const['Ms (A/m)'] = Const['Ms (emu/cm3)'] * 1e3
Const['MsError (A/m)'] = Const['MsError (emu/cm3)'] * 1e3
Const['H0 (A/m)'] = Const['H0 (Oe)'] / (4*np.pi * 1e-3)
Const['H0Error (A/m)'] = Const['H0Error (Oe)'] / (4*np.pi * 1e-3)

m_Delta_SI = m_Delta_cgs / (4*np.pi) # 1/m
m_Delta_err_SI = m_Delta_err_cgs / (4*np.pi) # 1/m

Params = {'m_Delta_SI (1/m)': m_Delta_SI,
          'm_Delta_err_SI': m_Delta_err_SI}


# Calculate SHA
def calc_m_alpha(f, g, phi, Ms, t, H0, Meff, MsErr, tErr, H0Err, MeffErr):
    hbar = constants.hbar
    e = constants.e
    mu_0 = constants.mu_0
    c1 = 2*np.pi*f/g*hbar/(2*e)
    sin_phi = np.sin(phi*np.pi/180)
    c2 = mu_0*Ms*t*(H0+0.5*Meff)
    c2_err = mu_0*np.sqrt((MsErr*t*(H0+0.5*Meff))**2+(Ms*tErr*(H0+0.5*Meff))**2+(Ms*t**H0Err)**2+(Ms*t*0.5*MeffErr)**2)
    m_alpha = c1 * sin_phi / c2
    m_alpha_err = c1 * sin_phi * c2_err / (c2**2)
    return m_alpha, m_alpha_err


def calc_RR(rho_FM, rho_SHM, rho_FM_err, rho_SHM_err):
    rr = 1 + rho_SHM / rho_FM
    rr_err = np.sqrt((rho_SHM_err/rho_FM)**2+(rho_SHM*rho_FM_err/rho_FM**2)**2)
    return rr, rr_err

def SHA(m_Delta, m_alpha, RR, A_C):
    return m_Delta / m_alpha * RR * A_C

m_alpha, m_alpha_err = calc_m_alpha(f, Const['g'], Const['phi (deg)'], Const['Ms (A/m)'], Const['t (m)'],
                                    Const['H0 (A/m)'], Const['Meffopt (A/m)'], Const['MsError (A/m)'], 
                                    Const['tError (m)'], Const['H0Error (A/m)'], Const['MeffoptError (A/m)'])

rr, rr_err = calc_RR(Const['rho_fm (muOhmcm)'], Const['rho_shm (muOhmcm)'], 
                     Const['rho_fm_err (muOhmcm)'], Const['rho_shm_err (muOhmcm)'])

A_C = Const['device_width (um)']*1e-6*Const['device_height (nm)']*1e-9 # m2

Params['m_alpha'] = m_alpha
Params['m_alpha_err'] = m_alpha_err
Params['rr'] = rr
Params['rr_err'] = rr_err

# !!! Continue here: Error of A_C, calc SHE


# Plotting























