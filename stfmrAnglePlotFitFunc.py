# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 13:07:03 2021

@author: rimmler
"""

import numpy as np
from scipy.optimize import curve_fit
from lmfit import Model
from scipy import constants

mu0 = constants.mu_0
e = constants.e
hbar = constants.Planck

def check_comps(comps):
    if comps not in ['x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz']:
        raise ValueError

def deg2rad(x):
    return x * np.pi / 180

def opt_Vs_ana_free(comps, phi_data_deg, Vs_data, phi_fit_deg):
    ''' Give components as "xy", "xyz" etc. '''
    check_comps(comps)
    phi_data = deg2rad(phi_data_deg)
    phi_plt = deg2rad(phi_fit_deg)

    def pref(phi, c):
        return -c * np.sin(2*phi)

    def Vs_x(phi, tau_xAD, c):
        return pref(phi, c) * tau_xAD * np.sin(phi)
    def Vs_y(phi, tau_yAD, c):
        return pref(phi, c) * tau_yAD * np.cos(phi)
    def Vs_z(phi, tau_zFL, c):
        return pref(phi, c) * tau_zFL

    def Vs_xy(phi, tau_xAD, tau_yAD, c):
        return Vs_x(phi, tau_xAD, c) + Vs_y(phi, tau_yAD, c)
    def Vs_xz(phi, tau_xAD, tau_zFL, c):
        return Vs_x(phi, tau_xAD, c) + Vs_z(phi, tau_zFL, c)
    def Vs_yz(phi, tau_yAD, tau_zFL, c):
        return Vs_y(phi, tau_yAD, c) + Vs_z(phi, tau_zFL, c)

    def Vs_xyz(phi, tau_xAD, tau_yAD, tau_zFL, c):
        return Vs_x(phi, tau_xAD, c) + Vs_y(phi, tau_yAD, c) + Vs_z(phi, tau_zFL, c)

    if comps == 'x':
        Vs_fit_opt, Vs_fit_cov = curve_fit(Vs_x, phi_data, Vs_data)
        Vs_fit = Vs_x(phi_data, Vs_fit_opt[0], Vs_fit_opt[-1])
        Vs_plt = Vs_x(phi_plt, Vs_fit_opt[0], Vs_fit_opt[-1])
    elif comps == 'y':
        Vs_fit_opt, Vs_fit_cov = curve_fit(Vs_y, phi_data, Vs_data)
        Vs_fit = Vs_y(phi_data, Vs_fit_opt[0], Vs_fit_opt[-1])
        Vs_plt = Vs_y(phi_plt, Vs_fit_opt[0], Vs_fit_opt[-1])
    elif comps == 'z':
        Vs_fit_opt, Vs_fit_cov = curve_fit(Vs_z, phi_data, Vs_data)
        Vs_fit = Vs_z(phi_data, Vs_fit_opt[0], Vs_fit_opt[-1])
        Vs_plt = Vs_z(phi_plt, Vs_fit_opt[0], Vs_fit_opt[-1])
    elif comps == 'xy':
        Vs_fit_opt, Vs_fit_cov = curve_fit(Vs_xy, phi_data, Vs_data)
        Vs_fit = Vs_xy(phi_data, Vs_fit_opt[0], Vs_fit_opt[1], Vs_fit_opt[-1])
        Vs_plt = Vs_xy(phi_plt, Vs_fit_opt[0], Vs_fit_opt[1], Vs_fit_opt[-1])
    elif comps == 'xz':
        Vs_fit_opt, Vs_fit_cov = curve_fit(Vs_xz, phi_data, Vs_data)
        Vs_fit = Vs_xz(phi_data, Vs_fit_opt[0], Vs_fit_opt[1], Vs_fit_opt[-1])
        Vs_plt = Vs_xz(phi_plt, Vs_fit_opt[0], Vs_fit_opt[1], Vs_fit_opt[-1])
    elif comps == 'yz':
        Vs_fit_opt, Vs_fit_cov = curve_fit(Vs_yz, phi_data, Vs_data)
        Vs_fit = Vs_yz(phi_data, Vs_fit_opt[0], Vs_fit_opt[1], Vs_fit_opt[-1])
        Vs_plt = Vs_yz(phi_plt, Vs_fit_opt[0], Vs_fit_opt[1], Vs_fit_opt[-1])
    elif comps == 'xyz':
        Vs_fit_opt, Vs_fit_cov = curve_fit(Vs_xyz, phi_data, Vs_data)
        Vs_fit = Vs_xyz(phi_data, Vs_fit_opt[0], Vs_fit_opt[1], Vs_fit_opt[2], Vs_fit_opt[-1])
        Vs_plt = Vs_xyz(phi_plt, Vs_fit_opt[0], Vs_fit_opt[1], Vs_fit_opt[2], Vs_fit_opt[-1])

    return Vs_fit_opt, Vs_fit_cov, Vs_fit, Vs_plt


def opt_Va_ana_free(comps, phi_data_deg, Va_data, phi_fit_deg):
    ''' Give components as "xy", "xyz" etc. '''
    check_comps(comps)
    phi_data = deg2rad(phi_data_deg)
    phi_plt = deg2rad(phi_fit_deg)

    def pref(phi, c):
        return -c * np.sin(2*phi)

    def Va_x(phi, tau_xFL, c):
        return pref(phi, c) * tau_xFL * np.sin(phi)
    def Va_y(phi, tau_yFL, c):
        return pref(phi, c) * tau_yFL * np.cos(phi)
    def Va_z(phi, tau_zAD, c):
        return pref(phi, c) * tau_zAD

    def Va_xy(phi, tau_xFL, tau_yFL, c):
        return Va_x(phi, tau_xFL, c) + Va_y(phi, tau_yFL, c)
    def Va_xz(phi, tau_xFL, tau_zAD, c):
        return Va_x(phi, tau_xFL, c) + Va_z(phi, tau_zAD, c)
    def Va_yz(phi, tau_yFL, tau_zAD, c):
        return Va_y(phi, tau_yFL, c) + Va_z(phi, tau_zAD, c)

    def Va_xyz(phi, tau_xFL, tau_yFL, tau_zAD, c):
        return Va_x(phi, tau_xFL, c) + Va_y(phi, tau_yFL, c) + Va_z(phi, tau_zAD, c)

    if comps == 'x':
        Va_fit_opt, Va_fit_cov = curve_fit(Va_x, phi_data, Va_data)
        Va_fit = Va_x(phi_data, Va_fit_opt[0], Va_fit_opt[-1])
        Va_plt = Va_x(phi_plt, Va_fit_opt[0], Va_fit_opt[-1])
    elif comps == 'y':
        Va_fit_opt, Va_fit_cov = curve_fit(Va_y, phi_data, Va_data)
        Va_fit = Va_y(phi_data, Va_fit_opt[0], Va_fit_opt[-1])
        Va_plt = Va_y(phi_plt, Va_fit_opt[0], Va_fit_opt[-1])
    elif comps == 'z':
        Va_fit_opt, Va_fit_cov = curve_fit(Va_z, phi_data, Va_data)
        Va_fit = Va_z(phi_data, Va_fit_opt[0], Va_fit_opt[-1])
        Va_plt = Va_z(phi_plt, Va_fit_opt[0], Va_fit_opt[-1])
    elif comps == 'xy':
        Va_fit_opt, Va_fit_cov = curve_fit(Va_xy, phi_data, Va_data)
        Va_fit = Va_xy(phi_data, Va_fit_opt[0], Va_fit_opt[1], Va_fit_opt[-1])
        Va_plt = Va_xy(phi_plt, Va_fit_opt[0], Va_fit_opt[1], Va_fit_opt[-1])
    elif comps == 'xz':
        Va_fit_opt, Va_fit_cov = curve_fit(Va_xz, phi_data, Va_data)
        Va_fit = Va_xz(phi_data, Va_fit_opt[0], Va_fit_opt[1], Va_fit_opt[-1])
        Va_plt = Va_xz(phi_plt, Va_fit_opt[0], Va_fit_opt[1], Va_fit_opt[-1])
    elif comps == 'yz':
        Va_fit_opt, Va_fit_cov = curve_fit(Va_yz, phi_data, Va_data)
        Va_fit = Va_yz(phi_data, Va_fit_opt[0], Va_fit_opt[1], Va_fit_opt[-1])
        Va_plt = Va_yz(phi_plt, Va_fit_opt[0], Va_fit_opt[1], Va_fit_opt[-1])
    elif comps == 'xyz':
        Va_fit_opt, Va_fit_cov = curve_fit(Va_xyz, phi_data, Va_data)
        Va_fit = Va_xyz(phi_data, Va_fit_opt[0], Va_fit_opt[1], Va_fit_opt[2], Va_fit_opt[-1])
        Va_plt = Va_xyz(phi_plt, Va_fit_opt[0], Va_fit_opt[1], Va_fit_opt[2], Va_fit_opt[-1])

    return Va_fit_opt, Va_fit_cov, Va_fit, Va_plt

# ____________________________________________________________________________
def opt_V_ana_full(comps, phi_data_deg, Vs_data, Va_data, phi_fit_deg, cps):
    ''' comps: components as "xy", "xyz" etc.
        cps: parameters required to calculate prefactor c'''
    check_comps(comps)
    phi_data = deg2rad(phi_data_deg)
    phi_plt = deg2rad(phi_fit_deg)

    def c_s(cps):
        Irf = cps['Irf (A)']
        DeltaR = cps['DeltaR (Ohm/rad)']
        alpha = cps['alpha'] # 1
        H0 = cps['H0_SI (A/m)']
        Meff = cps['Meff_SI (A/m)']
        return -Irf / 2. * DeltaR / ( alpha*mu0*(2*H0 + Meff))

    def c_a(cps):
        H0 = cps['H0_SI (A/m)']
        Meff = cps['Meff_SI (A/m)']
        return c_s(cps) * np.sqrt(1 + Meff/H0)

    def Vs_x(phi, tau_xAD, c):
        return c * np.sin(2*phi) * (tau_xAD * np.sin(phi))
    def Vs_y(phi, tau_yAD, c):
        return c * np.sin(2*phi) * (tau_yAD * np.cos(phi))
    def Vs_z(phi, tau_zFL, c):
        return c * np.sin(2*phi) * (tau_zFL)

    def Vs_xy(phi, tau_xAD, tau_yAD, c):
        return Vs_x(phi, tau_xAD, c) + Vs_y(phi, tau_yAD, c)
    def Vs_xz(phi, tau_xAD, tau_zFL, c):
        return Vs_x(phi, tau_xAD, c) + Vs_z(phi, tau_zFL, c)
    def Vs_yz(phi, tau_yAD, tau_zFL, c):
        return Vs_y(phi, tau_yAD, c) + Vs_z(phi, tau_zFL, c)

    def Vs_xyz(phi, tau_xAD, tau_yAD, tau_zFL, c):
        return Vs_x(phi, tau_xAD, c) + Vs_y(phi, tau_yAD, c) + Vs_z(phi, tau_zFL, c)

    #____________________
    def Va_x(phi, tau_xFL, c):
        return c * np.sin(2*phi) * (tau_xFL * np.sin(phi))
    def Va_y(phi, tau_yFL, c):
        return c * np.sin(2*phi) * (tau_yFL * np.cos(phi))
    def Va_z(phi, tau_zAD, c):
        return c * np.sin(2*phi) * (tau_zAD)

    def Va_xy(phi, tau_xFL, tau_yFL, c):
        return Va_x(phi, tau_xFL, c) + Va_y(phi, tau_yFL, c)
    def Va_xz(phi, tau_xFL, tau_zAD, c):
        return Va_x(phi, tau_xFL, c) + Va_z(phi, tau_zAD, c)
    def Va_yz(phi, tau_yFL, tau_zAD, c):
        return Va_y(phi, tau_yFL, c) + Va_z(phi, tau_zAD, c)

    def Va_xyz(phi, tau_xFL, tau_yFL, tau_zAD, c):
        return Va_x(phi, tau_xFL, c) + Va_y(phi, tau_yFL, c) + Va_z(phi, tau_zAD, c)

    # Construct tau based on given comps:
    C_s = c_s(cps)
    C_a = c_a(cps)
    if comps == 'x':
        Vs_model = Model(Vs_x)
        Vs_params = Vs_model.make_params(tau_xAD=1, c=C_s)
        Va_model = Model(Va_x)
        Va_params = Va_model.make_params(tau_xFL=1, c=C_a)
    elif comps == 'y':
        Vs_model = Model(Vs_y)
        Vs_params = Vs_model.make_params(tau_yAD=1, c=C_s)
        Va_model = Model(Va_y)
        Va_params = Va_model.make_params(tau_yFL=1, c=C_a)
    elif comps == 'z':
        Vs_model = Model(Vs_z)
        Vs_params = Vs_model.make_params(tau_zFL=1, c=C_s)
        Va_model = Model(Va_z)
        Va_params = Va_model.make_params(tau_zAD=1, c=C_a)
    elif comps == 'xy':
        Vs_model = Model(Vs_xy)
        Vs_params = Vs_model.make_params(tau_xAD=1, tau_yAD=1, c=C_s)
        Va_model = Model(Va_xy)
        Va_params = Va_model.make_params(tau_xFL=1, tau_yFL=1, c=C_a)
    elif comps == 'xz':
        Vs_model = Model(Vs_xz)
        Vs_params = Vs_model.make_params(tau_xAD=1, tau_zFL=1, c=C_s)
        Va_model = Model(Va_xz)
        Va_params = Va_model.make_params(tau_xFL=1, tau_zAD=1, c=C_a)
    elif comps == 'yz':
        Vs_model = Model(Vs_yz)
        Vs_params = Vs_model.make_params(tau_yAD=1, tau_zFL=1, c=C_s)
        Va_model = Model(Va_yz)
        Va_params = Va_model.make_params(tau_yFL=1, tau_zAD=1, c=C_a)
    elif comps == 'xyz':
        Vs_model = Model(Vs_xyz)
        Vs_params = Vs_model.make_params(tau_xAD=1, tau_yAD=1, tau_zFL=1, c=C_s)
        Va_model = Model(Va_xyz)
        Va_params = Va_model.make_params(tau_xFL=1, tau_yFL=1, tau_zAD=1, c=C_a)

    Vs_params['c'].vary = False
    Va_params['c'].vary = False
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

    Vs_fit = Vs_xyz(phi_data, torques['tau_xAD'], torques['tau_yAD'],torques['tau_zFL'], C_s)
    Vs_plt = Vs_xyz(phi_plt, torques['tau_xAD'], torques['tau_yAD'],torques['tau_zFL'], C_s)
    Va_fit = Va_xyz(phi_data, torques['tau_xFL'], torques['tau_yFL'],torques['tau_zAD'], C_a)
    Va_plt = Va_xyz(phi_plt, torques['tau_xFL'], torques['tau_yFL'],torques['tau_zAD'], C_a)


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

























