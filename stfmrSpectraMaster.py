# -*- coding: utf-8 -*-
'''
Analysis module for fitting of  ST-FMR spectra taken with RT ST-FMR setup at MPI Halle.

Author: 
    Berthold Rimmler, 
    Max Planck Institute of Microstructure Physics, Halle
    Weinberg 2
    06120 Halle
    berthold.rimmler@mpi-halle.mpg.de

'''

''' Input zone '''

# MAIN INFO
# Definition of angles:
measMode = 0 # 0: normal, 1: angle-dependence (specify stageAngle in separate file)
deviceAngle = 0 # Angle of device with respect to sample/crystal axes
stageAngle = 45 # For angle-dependence defined in separate file
offsetField = -16 # Oe, offset field of the magnet
baseVoltageMultiplier = 1 # Multiplies the voltageArray. For angle dependent measurements multiplied by value specified in freq-angle-corresp. file

system = 'Berthold' # use to specify file system in stfmrAnalysis

# ______________
# OTHER SETTINGS
plotAndCheck = 0 # 1: Check each file for selecting the right range, 0: for analysis without checking
plotAllTogether = 0 # 0: Separate files for each curve, 1: Plot all curves together, 2: Separate files as well as together
legendMode = 3 # 0: Frequency, 1: Current, 2: rf power, 3: Current and Frequency, Any other Number: All three written together
numberOFHeaderLines = 4

# File Format Inputs
hAxis = 0 # Unit Oe
vAxis = 1 # Unit V
mSize = 3 # Size of markers in plot

''' Input zone ends here. '''
#______________________________________________________________________________
import os
import tkinter as tk
from tkinter import filedialog
from modules.stfmrAnalysis import stfmrAnalysis
from files import File

def ui_get_IPFiles():
    root = tk.Tk()
    root.withdraw()
    inputFolder= filedialog.askdirectory(parent=root, title='Choose directory with input data files.')
    inputFileNames = os.listdir(inputFolder)
    root.destroy()
    IPFiles = []
    for fileName in inputFileNames:
        fullInputFileName = inputFolder + '/' + fileName
        file = File(fullInputFileName)
        if file.fileExt == '.txt':
            IPFiles.append(file)
    return IPFiles

def ui_get_facFile():
    ''' Get the .csv file with file-angle-correspondence '''
    root = tk.Tk()
    root.withdraw()
    facfile = filedialog.askopenfilename(parent=root, title='Choose .csv file with file-angle-correspondence')
    return File(facfile)

IPFiles = ui_get_IPFiles()
if measMode == 1:
    facFile = ui_get_facFile()
else:
    facFile = None


stfmrAna = stfmrAnalysis(measMode, deviceAngle, stageAngle, offsetField, IPFiles, facFile, plotAndCheck,
              plotAllTogether, legendMode, numberOFHeaderLines,
              hAxis, vAxis, baseVoltageMultiplier, mSize, system)
stfmrAna.do()








