# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 15:10:22 2022

@author: rimmler
"""

import pandas as pd


def read_csv_Series(file):
    return pd.read_csv(file.fileDirName, header=None, index_col=0, squeeze=True)

# ipFileLocations = _read_csv_Series(ipFileLocationsFile.file_fulldir)

# ipFileLS = File(ipFileLocations['lineshape analysis output'])
# ipFilePhiDep = File(ipFileLocations['angle dependence fitting summary'])
# ipFileIrf = File(ipFileLocations['Irf calibration fitting output'])
# ipFileAMR = File(ipFileLocations['AMR measurement fitting output'])

# LSData = _read_csv_Series(ipFileLS.file_fulldir)
# PhiDepData = pd.read_csv(ipFilePhiDep.file_fulldir)
# # return PhiDepData
# IrfData = _read_csv_Series(ipFileIrf.file_fulldir)
# AMRData = _read_csv_Series(ipFileAMR.file_fulldir)

# alpha = LSData['alphaopt']
# Meff = LSData['Meffopt (emu/cm3)'] # emu/cm3
# Meff_SI = Meff * 1e3 # A/m

# H0 = np.average(PhiDepData['Hres (Oe)'].to_numpy()) # Oe
# H0_SI = H0*1e3/(4*np.pi) # A/m

# t = LSData['t (m)']
# d = LSData['d (m)']
# Ms = LSData['Ms (emu/cm3)']
# Ms_SI = Ms * 1e3 # A/m

# m = float(IrfData['m (A/sqrt(mW))'])

# DeltaR = AMRData['DeltaR_fit (Ohm)']

# cps = {
#     'alpha': alpha,
#     't (m)': t,
#     'd (m)': d,
#     'Ms (emu/cm3)': Ms,
#     'Ms_SI (A/m)': Ms_SI,
#     'Meff (emu/cm3)': Meff,
#     'Meff_SI (A/m)': Meff_SI,
#     'H0 (Oe)': H0,
#     'H0_SI (A/m)': H0_SI,
#     'PL_rf (dBm)': Prf_dBm,
#     'P_rf (mW)': Prf_mW,
#     'm (A/sqrt(mW))': m,
#     'Irf (A)': Irf,
#     'DeltaR (Ohm/rad)': float(DeltaR)
#     }
# if print_extracted_params is True:
#     print_dict(cps, title='Extracted parameters for angle-dependence fitting:')
# return cps












