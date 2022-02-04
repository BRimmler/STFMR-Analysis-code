# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 15:09:08 2022

@author: rimmler
"""

# ____________________________________________________________________________
# INPUT
# Directory where to find the input file with file locations and additional parameters
ipDir = r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\MA2959-1\210805\005_linewidth_D1\linewidthAnalysis'
ipFileLocationsFileName = r'linewidth_input_files.csv'
ipParamsFileName = r'linewidth_input_params.csv'

plotDpi = 600

# ____________________________________________________________________________
# CODE
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import helpers.stfmrLineWidthAnalysisHelpers as hlp
from scipy.optimize import curve_fit
from base_modules.files import File
import base_modules.helpers as ut
from base_modules.plots import GenPlot, BoxText
from scipy import constants

# Get input data

# root = tk.Tk()
# root.withdraw()
# inputFile = File(filedialog.askopenfilename(parent=root, title='Choose .csv file with fitting summary'))
ipFileLocationsFile = File(ipDir, ipFileLocationsFileName)
ipParamsFile = File(ipDir, ipParamsFileName)

ipParams = hlp.read_csv_Series(ipParamsFile)
Const = ipParams.to_dict()

ipFileLocations = hlp.read_csv_Series(ipFileLocationsFile)

# Get fitting output of dc bias measurements and lineshape analysis output
ipFittingOutputFile = File(ipFileLocations['linewidth raw fitting output'])
ipLineshapeAnaOutputFile = File(ipFileLocations['lineshape analysis output'])

ipFittingOutput = pd.read_csv(ipFittingOutputFile.fileDirName, index_col='Index')
f_data = ipFittingOutput['Frequency (GHz)']
if len(f_data.unique()) > 1:
    raise ValueError('Frenquency is not constant.')
f = f_data[0]*1e9 # to Hz
Const['f (Hz)'] = f

ipLineshapeAnaOutput = hlp.read_csv_Series(ipLineshapeAnaOutputFile)
Const = {**Const, **ipLineshapeAnaOutput}

# Linear fit dDelta over dIdc
Idc = ipFittingOutput['Current (mA)']
Delta = ipFittingOutput['DeltaH (Oe)']
def fDelta(Idc, m, b):
    return m * Idc + b
popt, pcov = curve_fit(fDelta, Idc, Delta)
m_Delta = popt[0] # Oe/mA

# Clean unit --> Convert everything to SI
Const['l_C (m)'] = Const['l_C (mum)'] * 1e-6
Const['Meffopt (emu/cm3)']

# CLEAN UNIT HERE


# Calculate SHA
def calc_m_alpha(f, g, phi, Ms, t, H0, Meff):
    hbar = constants.hbar
    e = constants.e
    mu_0 = constants.mu_0
    c1 = 2*np.pi*f/g*hbar/(2*e**2) # e**2 or not??!!
    sin_phi = np.sin(phi*np.pi/180)
    c2 = mu_0*Ms*t*(H0+0.5*Meff)
    m_alpha = c1 * sin_phi / c2
    return m_alpha

def calc_RR(R_FM, R_SHM):
    return (R_FM + R_SHM) / R_FM

def calc_A_C(t, d, l_C, w_C):
    return (t + d) * l_C * w_C

def SHA(m_Delta, m_alpha, RR, A_C):
    return m_Delta / m_alpha * RR * A_C


# Plotting























