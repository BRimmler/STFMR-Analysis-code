# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 13:45:57 2021

@author: rimmler
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 13:07:03 2021

@author: rimmler
"""

import numpy as np
import scipy as sp


def opt_Vs_ana(comps, phi_data_deg, Vs_data, phi_fit_deg):
    ''' Give components as "xy", "xyz" etc. '''
    if comps not in ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz']:
        raise ValueError

    def deg2rad(x):
        return x * np.pi / 180
    phi_data = deg2rad(phi_data_deg)
    phi_fit = deg2rad(phi_fit_deg)

    def pref(phi, c):
        return c * np.sin(2*phi)

    def t_sx(tau_xAD, phi):
        return tau_xAD * np.sin(phi)
    def t_sy(tau_yAD, phi):
        return tau_yAD * np.cos(phi)
    def t_sz(tau_zFL, phi):
        return tau_zFL

    def Vs_x(phi, tau_xAD, c):
        return pref(phi, c) * t_sx(tau_xAD, phi)
    def Vs_y(phi, tau_yAD, c):
        return pref(phi, c) * t_sy(tau_yAD, phi)
    def Vs_z(phi, tau_zFL, c):
        return pref(phi, c) * t_sz(tau_zFL, phi)

    def Vs_xy(phi, tau_xAD, tau_yAD, c):
        return pref(phi, c) * (Vs_x(phi, tau_xAD, c) + Vs_y(phi, tau_yAD, c))
    def Vs_xz(phi, tau_xAD, tau_zFL, c):
        return pref(phi, c) * (Vs_x(phi, tau_xAD, c) + Vs_z(phi, tau_zFL, c))
    def Vs_yz(phi, tau_yAD, tau_zFL, c):
        return pref(phi, c) * (Vs_y(phi, tau_yAD, c) + Vs_z(phi, tau_zFL, c))

    def Vs_xyz(phi, tau_xAD, tau_yAD, tau_zFL, c):
        return pref(phi, c) * (Vs_x(phi, tau_xAD, c) + Vs_y(phi, tau_yAD, c) + Vs_z(phi, tau_zFL, c))

    if comps == 'x':
        Vs_fit_opt, Vs_fit_cov = sp.optimize.curve_fit(Vs_x, phi_data, Vs_data)
        Vs_fit = Vs_x(phi_fit, Vs_fit_opt[0], Vs_fit_opt[-1])
    elif comps == 'y':
        Vs_fit_opt, Vs_fit_cov = sp.optimize.curve_fit(Vs_y, phi_data, Vs_data)
        Vs_fit = Vs_y(phi_fit, Vs_fit_opt[0], Vs_fit_opt[-1])
    elif comps == 'z':
        Vs_fit_opt, Vs_fit_cov = sp.optimize.curve_fit(Vs_z, phi_data, Vs_data)
        Vs_fit = Vs_z(phi_fit, Vs_fit_opt[0], Vs_fit_opt[-1])
    elif comps == 'xy':
        Vs_fit_opt, Vs_fit_cov = sp.optimize.curve_fit(Vs_xy, phi_data, Vs_data)
        Vs_fit = Vs_xy(phi_fit, Vs_fit_opt[0], Vs_fit_opt[1], Vs_fit_opt[-1])
    elif comps == 'xz':
        Vs_fit_opt, Vs_fit_cov = sp.optimize.curve_fit(Vs_xz, phi_data, Vs_data)
        Vs_fit = Vs_xz(phi_fit, Vs_fit_opt[0], Vs_fit_opt[1], Vs_fit_opt[-1])
    elif comps == 'yz':
        Vs_fit_opt, Vs_fit_cov = sp.optimize.curve_fit(Vs_yz, phi_data, Vs_data)
        Vs_fit = Vs_yz(phi_fit, Vs_fit_opt[0], Vs_fit_opt[1], Vs_fit_opt[-1])
    elif comps == 'xyz':
        Vs_fit_opt, Vs_fit_cov = sp.optimize.curve_fit(Vs_xyz, phi_data, Vs_data)
        Vs_fit = Vs_xyz(phi_fit, Vs_fit_opt[0], Vs_fit_opt[1], Vs_fit_opt[2], Vs_fit_opt[-1])

    return Vs_fit_opt, Vs_fit_cov, Vs_fit

def opt_Va_ana(comps, phi_data_deg, Va_data, phi_fit_deg):
    ''' Give components as "xy", "xyz" etc. '''
    if comps not in ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz']:
        raise ValueError

    def deg2rad(x):
        return x * np.pi / 180
    phi_data = deg2rad(phi_data_deg)
    phi_fit = deg2rad(phi_fit_deg)

    def pref(phi, c):
        return c * np.sin(2*phi)

    def Va_x(phi, tau_xFL, c):
        return pref(phi, c) * tau_xFL * np.sin(phi)
    def Va_y(phi, tau_yFL, c):
        return pref(phi, c) * tau_yFL * np.cos(phi)
    def Va_z(phi, tau_zAD, c):
        return pref(phi, c) * tau_zAD

    def Va_xy(phi, tau_xFL, tau_yFL, c):
        return pref(phi, c) * (Va_x(phi, tau_xFL, c) + Va_y(phi, tau_yFL, c))
    def Va_xz(phi, tau_xFL, tau_zAD, c):
        return pref(phi, c) * (Va_x(phi, tau_xFL, c) + Va_z(phi, tau_zAD, c))
    def Va_yz(phi, tau_yFL, tau_zAD, c):
        return pref(phi, c) * (Va_y(phi, tau_yFL, c) + Va_z(phi, tau_zAD, c))

    def Va_xyz(phi, tau_xFL, tau_yFL, tau_zAD, c):
        return pref(phi, c) * (Va_x(phi, tau_xFL, c) + Va_y(phi, tau_yFL, c) + Va_z(phi, tau_zAD, c))

    if comps == 'x':
        Va_fit_opt, Va_fit_cov = sp.optimize.curve_fit(Va_x, phi_data, Va_data)
        Va_fit = Va_x(phi_fit, Va_fit_opt[0], Va_fit_opt[-1])
    elif comps == 'y':
        Va_fit_opt, Va_fit_cov = sp.optimize.curve_fit(Va_y, phi_data, Va_data)
        Va_fit = Va_y(phi_fit, Va_fit_opt[0], Va_fit_opt[-1])
    elif comps == 'z':
        Va_fit_opt, Va_fit_cov = sp.optimize.curve_fit(Va_z, phi_data, Va_data)
        Va_fit = Va_z(phi_fit, Va_fit_opt[0], Va_fit_opt[-1])
    elif comps == 'xy':
        Va_fit_opt, Va_fit_cov = sp.optimize.curve_fit(Va_xy, phi_data, Va_data)
        Va_fit = Va_xy(phi_fit, Va_fit_opt[0], Va_fit_opt[1], Va_fit_opt[-1])
    elif comps == 'xz':
        Va_fit_opt, Va_fit_cov = sp.optimize.curve_fit(Va_xz, phi_data, Va_data)
        Va_fit = Va_xz(phi_fit, Va_fit_opt[0], Va_fit_opt[1], Va_fit_opt[-1])
    elif comps == 'yz':
        Va_fit_opt, Va_fit_cov = sp.optimize.curve_fit(Va_yz, phi_data, Va_data)
        Va_fit = Va_yz(phi_fit, Va_fit_opt[0], Va_fit_opt[1], Va_fit_opt[-1])
    elif comps == 'xyz':
        Va_fit_opt, Va_fit_cov = sp.optimize.curve_fit(Va_xyz, phi_data, Va_data)
        Va_fit = Va_xyz(phi_fit, Va_fit_opt[0], Va_fit_opt[1], Va_fit_opt[2], Va_fit_opt[-1])

    return Va_fit_opt, Va_fit_cov, Va_fit