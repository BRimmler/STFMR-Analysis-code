# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 10:27:24 2022

@author: rimmler
"""

from numpy import sin, cos, sqrt
from scipy.constants import mu_0

def phi_offset(phi, phi0):
    return phi - phi0

def Vs(Vamr, Ts):
    return Vamr * Ts

def Va(Vamr, Ta):
    return Vamr * Ta

def Ts(ts, tau_par):
    return ts * tau_par

def Ta(ta, tau_perp):
    return ta * tau_perp

def ts(alpha, H0, Meff):
    return 1/(mu_0*alpha*(2*H0+Meff))

def ta(alpha, H0, Meff):
    return sqrt(1+Meff/H0)/(mu_0*alpha*(2*H0+Meff))

def tau_par(phi, phi0, xAD, yAD, zFL):
    return (xAD*sin(phi_offset(phi, phi0)) + yAD*cos(phi_offset(phi, phi0)) + zFL)*sin(2*phi_offset(phi, phi0))

def tau_perp(phi, phi0, xFL, yFL, zAD):
    return (xFL*sin(phi_offset(phi, phi0)) + yFL*cos(phi_offset(phi, phi0)) + zAD)*sin(2*phi_offset(phi, phi0))

def comp_Vs(phi, phi0, Vamr, alpha, H0, Meff, xAD, yAD, zFL):
    return Vs(Vamr, Ts(ts(alpha, H0, Meff), tau_par(phi, phi0, xAD, yAD, zFL)))

def comp_Va(phi, phi0, Vamr, alpha, H0, Meff, xFL, yFL, zAD):
    return Va(Vamr, Ta(ta(alpha, H0, Meff), tau_perp(phi, phi0,xFL, yFL, zAD)))