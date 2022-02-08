# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 13:07:03 2021

@author: rimmler
"""

import numpy as np
from scipy import constants

mu0 = constants.mu_0

def calc_a(cps):
    Irf = cps['Irf (A)']
    DeltaR = cps['DeltaR (Ohm/rad)']
    alpha = cps['alpha'] # 1
    H0 = cps['H0_SI (A/m)']
    Meff = cps['Meff_SI (A/m)']
    return -Irf / 2. * DeltaR / ( alpha*mu0*(2*H0 + Meff))

def calc_c_s(cps):
    return calc_a(cps)

def calc_c_a(cps):
    H0 = cps['H0_SI (A/m)']
    Meff = cps['Meff_SI (A/m)']
    return calc_a(cps) * np.sqrt(1 + Meff/H0)

def Vs_x(phi, tau_xAD, c_s):
    return c_s * np.sin(2*phi) * (tau_xAD * np.sin(phi))
def Vs_y(phi, tau_yAD, c_s):
    return c_s * np.sin(2*phi) * (tau_yAD * np.cos(phi))
def Vs_z(phi, tau_zFL, c_s):
    return c_s * np.sin(2*phi) * (tau_zFL)

def Vs_xy(phi, tau_xAD, tau_yAD, c_s):
    return Vs_x(phi, tau_xAD, c_s) + Vs_y(phi, tau_yAD, c_s)
def Vs_xz(phi, tau_xAD, tau_zFL, c_s):
    return Vs_x(phi, tau_xAD, c_s) + Vs_z(phi, tau_zFL, c_s)
def Vs_yz(phi, tau_yAD, tau_zFL, c_s):
    return Vs_y(phi, tau_yAD, c_s) + Vs_z(phi, tau_zFL, c_s)

def Vs_xyz(phi, tau_xAD, tau_yAD, tau_zFL, c_s):
    return Vs_x(phi, tau_xAD, c_s) + Vs_y(phi, tau_yAD, c_s) + Vs_z(phi, tau_zFL, c_s)

def Va_x(phi, tau_xFL, c_a):
    return c_a * np.sin(2*phi) * (tau_xFL * np.sin(phi))
def Va_y(phi, tau_yFL, c_a):
    return c_a * np.sin(2*phi) * (tau_yFL * np.cos(phi))
def Va_z(phi, tau_zAD, c_a):
    return c_a * np.sin(2*phi) * (tau_zAD)

def Va_xy(phi, tau_xFL, tau_yFL, c_a):
    return Va_x(phi, tau_xFL, c_a) + Va_y(phi, tau_yFL, c_a)
def Va_xz(phi, tau_xFL, tau_zAD, c_a):
    return Va_x(phi, tau_xFL, c_a) + Va_z(phi, tau_zAD, c_a)
def Va_yz(phi, tau_yFL, tau_zAD, c_a):
    return Va_y(phi, tau_yFL, c_a) + Va_z(phi, tau_zAD, c_a)

def Va_xyz(phi, tau_xFL, tau_yFL, tau_zAD, c_a):
    return Va_x(phi, tau_xFL, c_a) + Va_y(phi, tau_yFL, c_a) + Va_z(phi, tau_zAD, c_a)

