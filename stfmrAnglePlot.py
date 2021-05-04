# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 10:29:28 2021

@author: rimmler
"""

# ____________________________________________________________________________
# SETTINGS
voltageMagnitude = 'mu' # V
plotPhiMode = 1 # 0: raw angle, 1: shifted angle

plotDpi = 300


# ____________________________________________________________________________
# CODE
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import matplotlib.pyplot as plt


class File:
    def __init__(self, file):
        self.file_fulldir = file
        self.filedir = '/'.join(file.split('/')[:-1])
        self.filename = file.split('/')[-1]
        self.fileext = '.' + self.filename.split('.')[-1]
        self.filename_wo_ext = '.'.join(self.filename.split('.')[:-1])

root = tk.Tk()
root.withdraw()
inputFile = File(filedialog.askopenfilename(parent=root, title='Choose .csv file with fitting summary'))
# inputFile = File(r'D:\DATA\Mn3SnN\ST-FMR\MA2427-1\210331\004_angle_dependence\fittingOutput\000_fittingSummary.csv')
outputFile = File(inputFile.filedir + '/' + inputFile.filename_wo_ext + '_anglePlot.png')

inputData = pd.read_csv(inputFile.file_fulldir,index_col=False )


if voltageMagnitude == 'mu':
    y_label = 'V ($\mu$V)'
    voltageDivider = 1e-6

    if plotPhiMode == 0:
        x =  inputData['Angle (deg)']
        x_label = '$\Theta$ (deg)'

        Vs = inputData['Vsym (V)'] / voltageDivider
        Vas = inputData['Vas (V)'] / voltageDivider

    elif plotPhiMode == 1:
        x =  inputData.sort_values(by='AngleShifted (deg)')['AngleShifted (deg)']
        x_label = '$\phi$ (deg)'

        Vs = inputData.sort_values(by='AngleShifted (deg)')['Vsym (V)'] / voltageDivider
        Vas = inputData.sort_values(by='AngleShifted (deg)')['Vas (V)'] / voltageDivider

I = float(inputData['Current (mA)'][0])
P = float(inputData['rf Power (dBm)'][0])
f = float(inputData['Frequency (GHz)'][0])



plt.scatter(x, Vs, label='Vs')
plt.plot(x, Vs)
plt.scatter(x, Vas, label='Vas')
plt.plot(x, Vas)
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.legend()
plt.title('I = {} mA, f = {} GHz, P = {} dBm'.format(I, f, P))
plt.show
plt.savefig(outputFile.file_fulldir, bbox_inches="tight", dpi=plotDpi)











