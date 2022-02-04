# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 14:53:43 2020

@author: apandeya
"""


# import sys
import numpy as np
import pandas as pd
from itertools import zip_longest
import os
import tkinter as tk
from lineWidth import lineWidth
import matplotlib.pyplot as plt
# import matplotlib.font_manager as font_manager

'''
Function to check if a value is float
'''


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


# Function to write a list into a file
def changeListToString(list1):
    string1 = ""
    for item in list1:
        string1 = string1 + item + ","
    string1 = string1[:-1] + "\n"
    return string1

root = tk.Tk()
root.lift()
root.withdraw()
# path = tk.filedialog.askdirectory(parent=root, title='Choose a folder containing a fit parameter file')
path = tk.filedialog.askopenfilename(parent=root, title='Select input file')
number_of_header_lines = 2
currentAxis = 2
frequencyAxis = 4
delataHAxis = 5
delataHErrorAxis = 6
hResAxis = 7
error_in_hRes = 50
file_pre_slice = 5
file_post_slice = 43
#path = "toProcess"
nameFileArray = tuple(path + '/' + file for file in os.listdir(path) if file.endswith(".dat") or file.endswith("csv"))


lineWidthAnalysisArray = []
couldNotAnalyzeArray = []
def listToCsvString(*args):
    outstring = ""
    for item in args:
        outstring = outstring + str(item) + ","
    outstring += "\n"
    return outstring

def formatField (Field, x, xerr):
    """


    Parameters
    ----------
    Field : The name of the field
    x : Value of the field
    xerr : Error in the field

    Returns
    -------
    Str
        Field: x +/- xerr

    """
    return str(Field) + " " + str(x) + " +/- " + str(xerr) + "\n"
def formatWithComma (Field, x, xerr):
    """


    Parameters
    ----------
    Field : The name of the field
    x : Value of the field
    xerr : Error in the field

    Returns
    -------
    Str
        Field: x +/- xerr

    """
    return str(Field) + ", " + str(x) + ", +/-, " + str(xerr) + "\n"

for inputFileName in nameFileArray:
    if inputFileName[-3:] == 'dat':
            delim = '\t'
    elif inputFileName[-3:] == 'csv':
        delim = ','
    fileArray = np.genfromtxt(inputFileName,  skip_header=2, delimiter = delim)
    hres = round(fileArray[0, hResAxis])
    currentArray = np.array(fileArray[:, currentAxis])
    deltaHArray = np.array(fileArray[:, delataHAxis])
    deltaHErrorArray = np.array(fileArray[:, delataHErrorAxis])
    frequency = fileArray[0, frequencyAxis]
    hRes = fileArray[0, hResAxis]
    fileName = inputFileName[len(path)+file_pre_slice: len(path) + file_post_slice]
    try:
        lineWidthAnalysisArray.append(lineWidth(currentArray, deltaHArray, deltaHErrorArray, frequency, hRes, fileName))
    except:
        couldNotAnalyzeArray.append(inputFileName)

path_out = path + '/output'
if os.path.isdir(path_out) == False:
    os.makedirs(path_out)
np.savetxt(path_out+"/NotAnalyzed.csv",couldNotAnalyzeArray, fmt='%s', delimiter= ","  )
with open(path_out + '/FittingParam.csv', "w") as outFile1:
    outFile1.write("FileName,DeviceName,Frequency,Positive Slope, Positive Slope Error, Negative Slope, Negative Slope Error, Average Slope, Average Slope Error, Positive DeltaH, Positive DeltaH Error, Negative DeltaH, Negative DeltaHError\n")
    outFile1.write(",,Ghz,Oe/mA,Oe/mA,Oe/mA,Oe/mA,Oe/mA,Oe/mA,Oe,Oe,Oe,Oe\n")
for lw_item1 in lineWidthAnalysisArray:
    if lw_item1.hRes > 0:
        for lw_item2 in lineWidthAnalysisArray:
            if lw_item2.hRes<0 and lw_item2.frequency == lw_item1.frequency:
                fig, ax = plt.subplots(1, 1, constrained_layout = True)
                fig.set_size_inches(8,6)
                lw_item1.fitLine()
                lw_item2.fitLine()

                #print (lw_item1.fileName)

                outputFileName = path + "/output/" + lw_item1.fileName + '_Hres' + str(lw_item1.hRes) + ".csv"
                header = "#" + formatWithComma("Positive Slope", lw_item1.slope, lw_item1.slopeError)
                header = header + "#" + formatWithComma("Negative Slope", lw_item2.slope, lw_item2.slopeError)
                header = header +  "#" + formatWithComma("Positive Delta H", lw_item1.deltaHAtZero, lw_item1.deltaHErrorAtZero)
                header = header +"#" +  formatWithComma("Negative Delta H", lw_item2.deltaHAtZero, lw_item2.deltaHErrorAtZero)
                header = header + "#" + "Frequency: " + str( lw_item1.frequency) + " GHz" + " \n"
                header = header + "#" +  lw_item1.fileName + " \n"
                header = header + 'Current, DHp, DHpError, DHpFitting DHn, DHnError, DHnFitting \n'
                header = header + 'mA, Oe, Oe, Oe, Oe, Oe, Oe \n'
                with open(outputFileName, "a") as outFile1:
                    outFile1.write(header)
                for (current, deltaH1, deltaH1Error, yArray1, deltaH2, deltaH2Error, yArray2) in\
                    zip_longest(lw_item1.currentArray,lw_item1.deltaHArray, lw_item1.deltaHErrorArray, lw_item1.yArray,\
                        lw_item2.deltaHArray, lw_item2.deltaHErrorArray, lw_item2.yArray, fillvalue=''):
                        with open(outputFileName, "a") as outFile1:
                            outFile1.write(listToCsvString(current, deltaH1, deltaH1Error, yArray1, deltaH2, deltaH2Error, yArray2))

                #np.savetxt(outputFileName, output_array, header=header,  delimiter = ',')

                #Plot the graphs

                ax.errorbar(lw_item1.currentArray, lw_item1.deltaHArray, yerr=lw_item1.deltaHErrorArray, fmt="or", label =str(lw_item1.frequency) + ' GHz,' + " +H")
                ax.errorbar(lw_item2.currentArray, lw_item2.deltaHArray, yerr=lw_item2.deltaHErrorArray, fmt="ob", label =str(lw_item1.frequency) + ' GHz,' + " -H")
                ax.plot(lw_item1.currentArray, lw_item1.yArray, "-r", label= "Slope: {:.3f} $\pm$ {:.3f}".format(lw_item1.slope, lw_item1.slopeError))
                ax.plot(lw_item2.currentArray, lw_item2.yArray, "-b", label= "Slope: {:.3f} $\pm$ {:.3f}".format(lw_item2.slope, lw_item2.slopeError))
                ax.set_ylabel(r"$\delta \Delta$H (Oe)")
                ax.set_xlabel(r"Current (mA)")

                ax.legend()

                plt.show()

                fig.savefig(outputFileName[:-4] + ".pdf" + ".pdf", dpi = 600, bbox_inches="tight")
                with open(path + '/output/FittingParam.csv', "a") as outFile1:
                    outFile1.write(str(lw_item1.fileName) + ',' + str(lw_item1.fileName[-9:-6]) + ',' + str(lw_item1.frequency) + ',' +
                                   str(lw_item1.slope) + ','  + str(lw_item1.slopeError)+ ',' +
                                   str(lw_item2.slope) + ',' + str(lw_item2.slopeError)+ ',' +
                                   str((abs(lw_item1.slope) + abs(lw_item2.slope))/2) + ','+ str((abs(lw_item1.slopeError) + abs(lw_item2.slopeError))/2)+ ',' +
                                   str(lw_item1.deltaHAtZero) + ','+ str(lw_item1.deltaHErrorAtZero)+ ',' +
                                   str(lw_item2.deltaHAtZero) + ',' + str(lw_item2.deltaHErrorAtZero)+ '\n')

                plt.close()







