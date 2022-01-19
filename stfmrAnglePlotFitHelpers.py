# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 11:46:02 2021

@author: rimmler
"""

import tkinter as tk
from tkinter import filedialog
from stfmrHelpers import File
import numpy as np
import pandas as pd

def conv_dBm_mW(LP):
    ''' Convert Leistungspegel LP in dBm to power in mW '''
    return 10**(LP/10) # mW

def print_dict(dc, title=''):
    print(title)
    for key, val in dc.items():
        print(f'{key}: {val}')

def get_torques(fit_comps, Vs_fit_opt, Va_fit_opt):
    params = {}
    if fit_comps == 'x':
        params['tau_xAD'] = Vs_fit_opt[0]
        params['tau_xFL'] = Va_fit_opt[0]
    elif fit_comps == 'y':
        params['tau_yAD'] = Vs_fit_opt[0]
        params['tau_yFL'] = Va_fit_opt[0]
    elif fit_comps == 'z':
        params['tau_zFL'] = Vs_fit_opt[0]
        params['tau_zAD'] = Va_fit_opt[0]
    elif fit_comps == 'xy':
        params['tau_xAD'] = Vs_fit_opt[0]
        params['tau_yAD'] = Vs_fit_opt[1]
        params['tau_xFL'] = Va_fit_opt[0]
        params['tau_yFL'] = Va_fit_opt[1]
    elif fit_comps == 'xz':
        params['tau_xAD'] = Vs_fit_opt[0]
        params['tau_zFL'] = Vs_fit_opt[1]
        params['tau_xFL'] = Va_fit_opt[0]
        params['tau_zAD'] = Va_fit_opt[1]
    elif fit_comps == 'yz':
        params['tau_yAD'] = Vs_fit_opt[0]
        params['tau_zFL'] = Vs_fit_opt[1]
        params['tau_yFL'] = Va_fit_opt[0]
        params['tau_zAD'] = Va_fit_opt[1]
    elif fit_comps == 'xyz':
        params['tau_xAD'] = Vs_fit_opt[0]
        params['tau_yAD'] = Vs_fit_opt[1]
        params['tau_zFL'] = Vs_fit_opt[2]
        params['tau_xFL'] = Va_fit_opt[0]
        params['tau_yFL'] = Va_fit_opt[1]
        params['tau_zAD'] = Va_fit_opt[2]

    return params

def norm_torques(params):
    params_norm = {}
    if 'tau_yAD' in params.keys():

        params_norm['tau_yAD_norm'] = 1
        params_norm['tau_yFL_norm'] = 1

        if 'tau_xAD' in params.keys():
            params_norm['tau_xAD_norm'] = params['tau_xAD'] / params['tau_yAD']
            params_norm['tau_xFL_norm'] = params['tau_xFL'] / params['tau_yFL']
        if 'tau_zFL' in params.keys():
            params_norm['tau_zFL_norm'] = params['tau_zFL'] / params['tau_yAD']
            params_norm['tau_zAD_norm'] = params['tau_zAD'] / params['tau_yFL']

    return params_norm

def get_cps(ipFileLocationsFile=None, print_extracted_params=False):
    ''' Get constant parameters '''
    def get_P(PhiDepData):
        P_dBm_array = PhiDepData['rf Power (dBm)'].to_numpy()
        P_dBm = np.average(P_dBm_array)
        for p in P_dBm_array:
            if abs(p-P_dBm) > 0:
                raise ValueError('Prf is not constant throughout the angle-dependent measurement')
        return P_dBm

    if ipFileLocationsFile == None:
        root = tk.Tk()
        root.withdraw()
        ipFileLocationsFile = File(filedialog.askopenfilename(parent=root, title='Choose .csv file with locations of required output files'))

    def _read_csv_Series(file):
        return pd.read_csv(file, header = None, index_col = 0, squeeze = True)
    ipFileLocations = _read_csv_Series(ipFileLocationsFile.file_fulldir)

    ipFileLS = File(ipFileLocations['lineshape analysis output'])
    ipFilePhiDep = File(ipFileLocations['angle dependence fitting summary'])
    ipFileIrf = File(ipFileLocations['Irf calibration fitting output'])
    ipFileAMR = File(ipFileLocations['AMR measurement fitting output'])

    LSData = _read_csv_Series(ipFileLS.file_fulldir)
    PhiDepData = pd.read_csv(ipFilePhiDep.file_fulldir)
    # return PhiDepData
    IrfData = _read_csv_Series(ipFileIrf.file_fulldir)
    AMRData = _read_csv_Series(ipFileAMR.file_fulldir)

    alpha = LSData['alphaopt']
    Meff = LSData['Meffopt (emu/cm3)'] # emu/cm3
    Meff_SI = Meff * 1e3 # A/m

    H0 = np.average(PhiDepData['Hres (Oe)'].to_numpy()) # Oe
    H0_SI = H0*1e3/(4*np.pi) # A/m

    t = LSData['t (m)']
    d = LSData['d (m)']
    Ms = LSData['Ms (emu/cm3)']
    Ms_SI = Ms * 1e3 # A/m

    m = float(IrfData['m (A/sqrt(mW))'])

    Prf_dBm = get_P(PhiDepData)
    Prf_mW = conv_dBm_mW(Prf_dBm)
    Irf = m * np.sqrt(Prf_mW) # A

    DeltaR = AMRData['DeltaR_fit (Ohm)']

    cps = {
        'alpha': alpha,
        't (m)': t,
        'd (m)': d,
        'Ms (emu/cm3)': Ms,
        'Ms_SI (A/m)': Ms_SI,
        'Meff (emu/cm3)': Meff,
        'Meff_SI (A/m)': Meff_SI,
        'H0 (Oe)': H0,
        'H0_SI (A/m)': H0_SI,
        'PL_rf (dBm)': Prf_dBm,
        'P_rf (mW)': Prf_mW,
        'm (A/sqrt(mW))': m,
        'Irf (A)': Irf,
        'DeltaR (Ohm/rad)': float(DeltaR)
        }
    if print_extracted_params is True:
        print_dict(cps, title='Extracted parameters for angle-dependence fitting:')
    return cps












