# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 08:25:11 2022

@author: rimmler
"""

# ____________________________________________________________________________
# Karimeddiny description

from numpy import sin, cos, sqrt
from scipy.constants import mu_0, physical_constants, hbar

mu_B = physical_constants['Bohr magneton'][0]

def calc_Vamr(cps):
    return cps['Irf (A)'] * cps['DeltaR (Ohm/rad)']

def gamma(g_e):
    return g_e * mu_B / hbar

def B0(H0):
    return mu_0 * H0

def w1(gamma, B0):
    return gamma*B0

def w2(gamma, B0, Meff):
    return gamma*(B0+mu_0*Meff)

def wp(w1, w2):
    return w1+w2

def tau_x(phi, tau_xDL, tau_yDL, tau_zFL):
    return tau_xDL*sin(phi)-tau_yDL*cos(phi)-tau_zFL

def tau_z(phi, tau_xFL, tau_yFL, tau_zDL):
    return tau_xFL*sin(phi)-tau_yFL*cos(phi)+tau_zDL

def Vart(Csp, Clssene, wp, w1, w2, alpha, gamma, tau_yDL, tau_yFL):
    Vart = (Csp + Clssene*wp/(alpha*gamma))*1/(alpha**2*wp**2)*(w1*tau_yDL**2 + w2*tau_yFL**2)
    return Vart

### Vs, Va
def Vs(phi, Vamr, Tart, alpha, w1, w2, wp, 
            tau_xDL, tau_xFL, tau_yDL, tau_yFL, tau_zDL, tau_zFL):
    tau_x_val = tau_x(phi, tau_xDL, tau_yDL, tau_zFL)
    tau_z_val = tau_z(phi, tau_xFL, tau_yFL, tau_zDL)
    Vs = Vamr/(2*alpha*wp)*sin(2*phi)*tau_x_val + Tart/((alpha*wp)**2)*sin(phi)*(w1*tau_x_val**2+w2*tau_z_val**2)
    return Vs

def Va(phi, Vamr, alpha, w1, w2, wp,
            tau_xDL, tau_xFL, tau_yDL, tau_yFL, tau_zDL, tau_zFL):
    tau_z_val = tau_z(phi, tau_xFL, tau_yFL, tau_zDL)
    Va = Vamr/(2*alpha*wp)*sin(2*phi)*sqrt(w2/w1)*tau_z_val
    return Va























