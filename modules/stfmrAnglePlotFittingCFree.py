# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 10:25:58 2022

@author: rimmler
"""
import matplotlib.pyplot as plt
import numpy as np

from lmfit import Parameters, minimize, report_fit
from modules.functions.stfmrAnglePlotMixingEffs import comp_Vs, comp_Va
from helpers.maths_helpers import deg2rad

default_torques = ['xAD', 'xFL', 'yAD', 'yFL', 'zAD', 'zFL']

def angleDepFittingCFree(phi_data, phi_plt, Vs_data, Va_data, 
                            cps, comps, fit_phi_offset=False, do_check_fit=False, do_lmfit_report=False):
    ''' Alternative fitting procedure for fitting angle dependence. Herein, 
    the AMR voltage is taken as a fitting parameter in the two equations of
    Vs and Va '''
    
    def Vs_dataset(params, phi):
        Vs = comp_Vs(phi, params['phi0_s'], params['Vamr_s'], 
                     params['alpha'], params['H0'], params['Meff'], 
                     params['xAD'], params['yAD'], params['zFL'])
        return Vs
    
    def Va_dataset(params, phi):
        Va = comp_Va(phi, params['phi0_a'],params['Vamr_a'], 
                     params['alpha'], params['H0'], params['Meff'], 
                     params['xFL'], params['yFL'], params['zAD'])
        return Va
    
    def objective(params, phi, Vs_data, Va_data):
        '''Calculate total residual for fits.'''
        data = np.array([Vs_data, Va_data])
        resid = 0.*data[:]
    
        # make residual per data set
        resid[0, :] = data[0, :] - Vs_dataset(params, phi)
        resid[1, :] = data[1, :] - Va_dataset(params, phi)
        
        # now flatten this to a 1D array, as minimize() needs
        return resid.flatten()
    
    # Get fixed parameters from cps
    alpha = cps['alpha']
    H0 = cps['H0_SI (A/m)']
    Meff = cps['Meff_SI (A/m)']
    
    # _________________________________________________________________________
    # Convert angle to rad
    phi_data = deg2rad(phi_data)
    phi_plt = deg2rad(phi_plt)
    
    fit_params = Parameters()
    
    # Fixed parameters:
    fit_params.add('alpha', value=alpha)
    fit_params.add('H0', value=H0)
    fit_params.add('Meff', value=Meff)
    for param in ['alpha', 'H0', 'Meff']:
        fit_params[param].vary=False
    
    # Fitting parameters:
    fit_params.add('Vamr_s', value=-1e-3, max=0, min=-1e-2) # AMR voltage in equation for Vs
    fit_params.add('Vamr_a', expr='Vamr_s')
    fit_params.add('phi0_s', value=0, max=deg2rad(10), min=deg2rad(-10))
    fit_params.add('phi0_a', expr='phi0_s')
    
    for comp in ['x', 'y', 'z']:
        if comp in comps: 
            # Case: If the component is in the fitting components (e.g. x is in xyz),
            # set its initial value to 1 and let it vary.
            fit_params.add(f'{comp}AD', value=1)
            fit_params.add(f'{comp}FL', value=1)
        else:
            # Else, fix its value at 0
            fit_params.add(f'{comp}AD', value=0)
            fit_params.add(f'{comp}FL', value=0)
            fit_params[f'{comp}AD'].vary=False
            fit_params[f'{comp}FL'].vary=False
    
    # _________________________________________________________________________
    # Run the global fit and show the fitting result
    out = minimize(objective, fit_params, args=(phi_data, Vs_data, Va_data))
    if do_lmfit_report is True:
        report_fit(out.params)
    
    # _________________________________________________________________________
    # Plot the data sets and fits
    Vs_fit = Vs_dataset(out.params, phi_data)
    Va_fit = Va_dataset(out.params, phi_data)
    Vs_plt = Vs_dataset(out.params, phi_plt)
    Va_plt = Va_dataset(out.params, phi_plt)
    
    if do_check_fit is True:
        plt.figure()
        plt.plot(phi_data, Vs_data, 'o', phi_data, Vs_fit, '-')
        plt.plot(phi_data, Va_data, 'o', phi_data, Va_fit, '-')
        
    # Get parameters in dict
    params_dict = {}
    for key in out.params.keys():
        params_dict[key] = out.params[key].value
        
    torques = {}
    for key in default_torques:
        torques[key] = out.params[key].value
        
    return out.params, params_dict, Vs_fit, Vs_plt, Va_fit, Va_plt



def get_norm_torques(params, norm_to):
    ''' Get normalized torques from fitting parameter output '''
    if not params['Vamr_s'].value == params['Vamr_a'].value:
        raise
        
    torques = {}
    torques_norm = {}
    for key in default_torques:
        tau = params[key].value
        tau_norm = tau / params[norm_to].value
        torques[key] = tau
        torques_norm[f'{key}_norm'] = tau_norm
        
    return torques, torques_norm
        
        
        
    
        


















