# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 13:23:00 2021

@author: rimmler
"""

tau_xAD = 0
tau_yAD = 1
tau_zAD = 1

tau_xFL = 0
tau_yFL = 0
tau_zFL = 0

c = 1

import numpy as np
import matplotlib.pyplot as plt


def deg2rad(x):
    return x * np.pi / 180

def pref(phi, c):
    return c * np.sin(2*phi)

def t_sx(tau_xAD, phi):
    return tau_xAD * np.sin(phi)
def t_sy(tau_yAD, phi):
    return tau_yAD * np.cos(phi)
def t_sz(tau_zFL, phi):
    return tau_zFL

def t_ax(tau_xFL, phi):
    return tau_xFL * np.sin(phi)
def t_ay(tau_yFL, phi):
    return tau_yFL * np.cos(phi)
def t_az(tau_zAD, phi):
    return tau_zAD

def Vs_x(phi, tau_xAD, c):
    return pref(phi, c) * t_sx(tau_xAD, phi)
def Vs_y(phi, tau_yAD, c):
    return pref(phi, c) * t_sy(tau_yAD, phi)
def Vs_z(phi, tau_zFL, c):
    return pref(phi, c) * t_sz(tau_zFL, phi)

def Va_x(phi, tau_xFL, c):
    return pref(phi, c) * t_ax(tau_xAD, phi)
def Va_y(phi, tau_yFL, c):
    return pref(phi, c) * t_ay(tau_yAD, phi)
def Va_z(phi, tau_zAD, c):
    return pref(phi, c) * t_az(tau_zFL, phi)

def Vs_xyz(phi, tau_xAD, tau_yAD, tau_zFL, c):
    return pref(phi, c) * (t_sx(tau_xAD, phi) + t_sy(tau_yAD, phi) + t_sz(tau_zFL, phi))
def Va_xyz(phi, tau_xFL, tau_yFL, tau_zAD, c):
    return pref(phi, c) * (t_ax(tau_xFL, phi) + t_ay(tau_yFL, phi) + t_az(tau_zAD, phi))


phi_deg = np.linspace(0, 360, 100)
phi_rad = deg2rad(phi_deg)

plt.plot(phi_deg, Vs_xyz(phi_rad, tau_xAD, tau_yAD, tau_zFL, c))
plt.plot(phi_deg, Va_xyz(phi_rad, tau_xFL, tau_yFL, tau_zAD, c))
plt.xticks(np.arange(0, 361, 60))