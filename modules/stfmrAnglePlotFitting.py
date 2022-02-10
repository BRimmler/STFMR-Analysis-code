# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 10:12:47 2022

@author: rimmler
"""

from lmfit import Model
from scipy import constants
import modules.functions.stfmrAnglePlotFitFunc as f
from helpers.maths_helpers import deg2rad

mu0 = constants.mu_0
e = constants.e
hbar = constants.Planck

def check_comps(comps):
    if comps not in ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz']:
        raise ValueError

# ____________________________________________________________________________
def opt_V_ana_full(comps, phi_data_deg, Vs_data, Va_data, phi_fit_deg, cps):
    ''' comps: components as "xy", "xyz" etc.
        cps: parameters required to calculate prefactor c'''
    check_comps(comps)
    phi_data = deg2rad(phi_data_deg)
    phi_plt = deg2rad(phi_fit_deg)

    # Construct tau based on given comps:
    C_s = f.calc_c_s(cps)
    C_a = f.calc_c_a(cps)
    Vs_model = Model(f.Vs_xyz)	
    Va_model = Model(f.Va_xyz)
    
    Vs_model.set_param_hint('c_s', value=C_s, vary=False)
    Va_model.set_param_hint('c_a', value=C_a, vary=False)
    
    if 'x' in comps:
        Vs_model.set_param_hint('tau_xAD', value=1)
        Va_model.set_param_hint('tau_xFL', value=1)
    else:
        Vs_model.set_param_hint('tau_xAD', value=0, vary=False)
        Va_model.set_param_hint('tau_xFL', value=0, vary=False)
        
    if 'y' in comps:
        Vs_model.set_param_hint('tau_yAD', value=1)
        Va_model.set_param_hint('tau_yFL', value=1)
    else:
        Vs_model.set_param_hint('tau_yAD', value=0, vary=False)
        Va_model.set_param_hint('tau_yFL', value=0, vary=False)
        
    if 'z' in comps:
        Vs_model.set_param_hint('tau_zFL', value=1)
        Va_model.set_param_hint('tau_zAD', value=1)
    else:
        Vs_model.set_param_hint('tau_zFL', value=0, vary=False)
        Va_model.set_param_hint('tau_zAD', value=0, vary=False)
        
    Vs_params = Vs_model.make_params()
    Va_params = Va_model.make_params()
            
    Vs_results = Vs_model.fit(Vs_data, Vs_params, phi=phi_data)
    Va_results = Va_model.fit(Va_data, Va_params, phi=phi_data)

    torques = {}
    for torque_key in ['tau_xAD', 'tau_yAD', 'tau_zFL']:
        if torque_key in Vs_results.params.keys():
            torques[torque_key] = Vs_results.params[torque_key].value
        else:
            torques[torque_key] = 0
    for torque_key in ['tau_xFL', 'tau_yFL', 'tau_zAD']:
        if torque_key in Va_results.params.keys():
            torques[torque_key] = Va_results.params[torque_key].value
        else:
            torques[torque_key] = 0

    Vs_fit = f.Vs_xyz(phi_data, torques['tau_xAD'], torques['tau_yAD'],torques['tau_zFL'], C_s)
    Vs_plt = f.Vs_xyz(phi_plt, torques['tau_xAD'], torques['tau_yAD'],torques['tau_zFL'], C_s)
    Va_fit = f.Va_xyz(phi_data, torques['tau_xFL'], torques['tau_yFL'],torques['tau_zAD'], C_a)
    Va_plt = f.Va_xyz(phi_plt, torques['tau_xFL'], torques['tau_yFL'],torques['tau_zAD'], C_a)

    return torques, Vs_fit, Vs_plt, Va_fit, Va_plt


def get_sotr(torques, cps):
    Ms = cps['Ms_SI (A/m)']
    d_Py = cps['d (m)']
    t_film = cps['t (m)']

    a = e * mu0 * Ms * d_Py * t_film / hbar
    tau_xAD = torques['tau_xAD']
    tau_yAD = torques['tau_yAD']
    tau_yFL = torques['tau_yFL']
    tau_zAD = torques['tau_zAD']

    theta_x = a * tau_xAD / tau_yFL
    theta_y = a * tau_yAD / tau_yFL
    theta_z = a * tau_zAD / tau_yFL

    sotr = {'theta_x': theta_x,
            'theta_y': theta_y,
            'theta_z': theta_z}

    return sotr

























