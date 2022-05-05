# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 13:00:04 2022

@author: rimmler
"""

from lmfit import Parameters, minimize
import modules.functions.stfmrAnglePlotFitFuncKarimeddiny as fk
from units import deg2rad
import numpy as np

default_torques = ['tau_xDL', 'tau_xFL', 'tau_yDL', 'tau_yFL', 'tau_zDL', 'tau_zFL']

def check_comps(comps):
    if comps not in ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz']:
        raise ValueError
        
def V_Karimeddiny_Hall_fitting(comps, phi_data_deg, Vset_mode,
                               Vsxx_data, Vaxx_data, Vsxy_data, Vaxy_data, 
                               phi_fit_deg, cps,
                               assume_arts=True, do_lmfit_report=False, do_check_fit=False):
    
    def Vsxx_dataset(params, phi):
        return fk.Vsxx(phi, params['Irf'], params['Ramr'], params['Tart'], params['alpha'],
                     params['w1'], params['w2'], params['wp'], params['W (m)'],
                     params['tau_xDL'], params['tau_xFL'], params['tau_yDL'],
                     params['tau_yFL'], params['tau_zDL'], params['tau_zFL'])
    
    def Vaxx_dataset(params, phi):
        return fk.Vaxx(phi, params['Irf'], params['Ramr'], params['alpha'],
                     params['w1'], params['w2'], params['wp'],
                     params['tau_xDL'], params['tau_xFL'], params['tau_yDL'],
                     params['tau_yFL'], params['tau_zDL'], params['tau_zFL'])
    
    def Vsxy_dataset(params, phi):
        return fk.Vsxy(phi, params['Irf'], params['Rphe'], params['Rahe'], params['Tart'], params['alpha'],
                     params['w1'], params['w2'], params['wp'], params['L (m)'],
                     params['tau_xDL'], params['tau_xFL'], params['tau_yDL'],
                     params['tau_yFL'], params['tau_zDL'], params['tau_zFL'])
    
    def Vaxy_dataset(params, phi):
        return fk.Vaxy(phi, params['Irf'], params['Rphe'], params['Rahe'], params['alpha'],
                     params['w1'], params['w2'], params['wp'],
                     params['tau_xDL'], params['tau_xFL'], params['tau_yDL'],
                     params['tau_yFL'], params['tau_zDL'], params['tau_zFL'])

    def objective(params, phi, Vset_mode, Vsxx_data, Vaxx_data, Vsxy_data, Vaxy_data):
        '''Calculate total residual for fits.'''
        V1 = Vsxx_data
        V2 = Vaxx_data
        if Vset_mode == 0:
            V3 = Vsxy_data
        elif Vset_mode == 1:
            V3 = Vaxy_data
        
        data = np.array([V1, V2, V3])
        resid = 0.*data[:]
        
        # make residual per data set
        resid[0, :] = data[0, :] - Vsxx_dataset(params, phi)
        resid[1, :] = data[1, :] - Vaxx_dataset(params, phi)
        if Vset_mode == 0:
            resid[2, :] = data[2, :] - Vsxy_dataset(params, phi)
        elif Vset_mode == 1:
            resid[2, :] = data[2, :] - Vaxy_dataset(params, phi)
        
        # now flatten this to a 1D array, as minimize() needs
        return resid.flatten()
    
    # Get fixed parameters from cps
    alpha = cps['alpha']
    H0 = cps['H0_SI (A/m)']
    B0 = fk.B0(H0)
    Meff = cps['Meff_SI (A/m)']
    gamma = fk.gamma(cps['g_e'])
    w1 = fk.w1(gamma, B0)
    w2 = fk.w2(gamma, B0, Meff)
    wp = fk.wp(w1, w2)
    Irf = cps['Irf (mA)']*1e-3 # to A
    Ramr = cps['Ramr (Ohm)']
    Rphe = cps['Rphe (Ohm)']
    Rahe = cps['Rahe (Ohm)']
    W = cps['W (m)']
    L = cps['L (m)']
    
    # Convert angle to rad
    check_comps(comps)
    phi_data = deg2rad(phi_data_deg)
    phi_plt = deg2rad(phi_fit_deg)
    
    # Fixed parameters
    fit_params = Parameters()
    
    fit_params.add('Irf', value=Irf)
    fit_params.add('Ramr', value=Ramr)
    fit_params.add('Rphe', value=Rphe)
    fit_params.add('Rahe', value=Rahe)
    fit_params.add('W', value=W)
    fit_params.add('L', value=L)
    
    fit_params.add('alpha', value=alpha)
    fit_params.add('w1', value=w1)
    fit_params.add('w2', value=w2)
    fit_params.add('wp', value=wp)
    
    for param in ['Irf', 'Ramr', 'Rphe', 'Rahe', 'W', 'L', 'alpha', 'w1', 'w2', 'wp']:
        fit_params[param].vary=False
    
    # Varying parameters (fit parameters)
    
    for comp in ['x', 'y', 'z']:
        if comp in comps: 
            # Case: If the component is in the fitting components (e.g. x is in xyz),
            # set its initial value to 1 and let it vary.
            fit_params.add(f'tau_{comp}DL', value=1)
            fit_params.add(f'tau_{comp}FL', value=1)
        else:
            # Else, fix its value at 0
            fit_params.add(f'tau_{comp}DL', value=0)
            fit_params.add(f'tau_{comp}FL', value=0)
            fit_params[f'tau_{comp}DL'].vary=False
            fit_params[f'tau_{comp}FL'].vary=False
            
    ### 
    # Iterative fitting to include artifacts
    if assume_arts is True:
        iterats = 2
    else:
        iterats = 1
        
    for i in range(iterats):
        if i == 0:
            # First iteration: No artifacts
            fit_params.add('Tart', value=0)
            fit_params['Tart'].vary=False
            
            out = minimize(objective, fit_params, args=(phi_data, Vsxx_data, Vaxx_data, Vsxy_data, Vaxy_data))
            # if do_lmfit_report is True:
            #     report_fit(out.params)
            
        elif i == 1:
            # Second iteration: if artifacts assumed: start from output params 
            # of previous fit and let Tart vary
            fit_params['Tart'].value = -1
            fit_params['Tart'].min = -1e3
            fit_params['Tart'].max = 0 
            fit_params['Tart'].vary = True
            
            for comp in ['x', 'y', 'z']:
                if comp in comps: 
                    # Case: If the component is in the fitting components (e.g. x is in xyz),
                    # set its initial value to 1 and let it vary.
                    fit_params[f'tau_{comp}DL'] = out.params[f'tau_{comp}DL']
                    fit_params[f'tau_{comp}FL'] = out.params[f'tau_{comp}FL']
                    
            out = minimize(objective, fit_params, args=(phi_data, Vset_mode, Vsxx_data, Vaxx_data, Vsxy_data, Vaxy_data))
    
    # _________________________________________________________________________
    # Run the global fit and show the fitting result
            
    # _________________________________________________________________________
    # Plot the data sets and fits
    Vsxx_fit = Vsxx_dataset(out.params, phi_data)
    Vaxx_fit = Vaxx_dataset(out.params, phi_data)
    Vsxx_plt = Vsxx_dataset(out.params, phi_plt)
    Vaxx_plt = Vaxx_dataset(out.params, phi_plt)
    
    Vsxy_fit = Vsxy_dataset(out.params, phi_data)
    Vaxy_fit = Vaxy_dataset(out.params, phi_data)
    Vsxy_plt = Vsxy_dataset(out.params, phi_plt)
    Vaxy_plt = Vaxy_dataset(out.params, phi_plt)
    
    # if do_check_fit is True:
    #     plt.figure()
    #     plt.plot(phi_data, Vs_data, 'o', phi_data, Vs_fit, '-')
    #     plt.plot(phi_data, Va_data, 'o', phi_data, Va_fit, '-')
    
    # Get parameters in dict
    params_dict = {}
    for key in out.params.keys():
        params_dict[key] = out.params[key].value
        
    return out.params, params_dict, Vsxx_fit, Vsxx_plt, Vaxx_fit, Vaxx_plt, Vsxy_fit, Vsxy_plt, Vaxy_fit, Vaxy_plt
    


def get_norm_torques_karimed(params, norm_to):
    ''' Get normalized torques from fitting parameter output '''
    if norm_to.split('_')[0] != 'tau':
        norm_to = 'tau_'+norm_to
    
    torques = {}
    torques_norm = {}
    for key in default_torques:
        tau = params[key].value
        tau_norm = tau / params[norm_to].value
        torques[key] = tau
        torques_norm[f'{key}_norm'] = tau_norm
        
    return torques, torques_norm
        

def calc_Ru(yobs, ycalc):
    return (np.sum(np.abs(yobs)**2) - np.sum(np.abs(ycalc)**2)) / np.sum(np.abs(yobs)**2)
        



