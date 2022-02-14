# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 10:57:20 2022

@author: rimmler
"""
import tkinter as tk
from tkinter import filedialog
import pandas as pd
import os
from files import File

def get_file_from_ui(title):
    root = tk.Tk()
    file = filedialog.askopenfilename(parent=root, title=title)
    
    root.withdraw()
    return file

def get_ipFileLocationsFilesFromUI(analysisMode):
    inputFileLocations = {}
    
    # Get angle-dependent data
    angle_dep_loc = get_file_from_ui('Select fitting summary of angle-dependent data.')
    inputFileLocations['angle dependence fitting summary'] = angle_dep_loc
    
    if analysisMode in [1, 2]:
        # Lineshape analysis        
        lineshape_loc = get_file_from_ui('Select fitting summary of lineshape analysis.')
        inputFileLocations['lineshape analysis output'] = lineshape_loc
        
    if analysisMode in [2]:
        # Irf
        irf_loc = get_file_from_ui('Select Irf calibration fitting output.')
        inputFileLocations['Irf calibration fitting output'] = irf_loc
        # AMR
        amr_loc = get_file_from_ui('Select AMR calibration fitting output.')
        inputFileLocations['AMR measurement fitting output'] = amr_loc
    
    inputFileLocations = pd.Series(inputFileLocations)
    current_path = File(os.path.realpath(__file__)).fileDir
    temp_path = current_path + '/temp'
    temp_file = File(temp_path, 'inputFileLocationsTemp.csv')
    temp_file.makeDirIfNotExist()
    inputFileLocations.to_csv(temp_file.fileDirName, header=False)
    return temp_file

        