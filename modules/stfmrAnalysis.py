''' Credit: Avanindra Kumar Pandeya, MPI Halle, Python God <3 '''

# import sys
# sys.path.append("..")
import os
from modules.stfmrSpectraFitting import stfmrSpectraFitting
import numpy as np
import matplotlib.pyplot as plt
import modules.stfmrRanges as stfmrRanges
from matplotlib import cm
import pandas as pd
# import tkinter as tk
# from tkinter import filedialog

class stfmrAnalysis:
    def __init__(self, measMode, deviceAngle, stageAngle, offsetField, IPFiles, facFile, plotANdCheck, 
                 plotAllTogether, legendMode, numberOFHeaderLines, hAxis, vAxis, 
                 baseVoltageMultiplier, mSize, system='default', dpi=600):
        self.measMode = measMode
        self.deviceAngle = deviceAngle
        self.stageAngle = stageAngle
        self.offsetField = offsetField
        self.IPFiles = IPFiles
        self.facFile = facFile
        self.plotANdCheck = plotANdCheck
        self.plotAllTogether = plotAllTogether
        self.legendMode = legendMode
        self.numberOFHeaderLines = numberOFHeaderLines
        self.hAxis = hAxis
        self.vAxis = vAxis
        self.baseVoltageMultiplier = baseVoltageMultiplier
        self.mSize = mSize
        self.system = system       
        self.dpi = dpi
        
    # ____________________________________________________________________________
    # Functions
    
    def parametersFromFileName(self, inputFileName, inputFileDir):
        underScoreNumber = 0

        i = len(inputFileName) - 1
        positionArray = np.array([])

        # Iterate till 1st element and keep on decrementing i
        while i >= 0:
            if inputFileName[i] == '_':
                underScoreNumber += 1
                positionArray = np.append(positionArray, i)
            i -= 1
        current =  float(inputFileName[int(positionArray[2])+1:int(positionArray[1])-2])
        dBm = float(inputFileName[int(positionArray[3])+1:int(positionArray[2])-3])
        frequency = float(inputFileName[int(positionArray[4])+1:int(positionArray[3])-3])
        name =  inputFileName[int(positionArray[-1])+1:int(positionArray[4])-3]
        number = float((inputFileName[len(inputFileDir)+1:len(inputFileDir)+4]))
        return current, dBm, frequency, name, number

    def get_angleShifted(self, angle, deviceAngle):
        angleDiff = angle - deviceAngle
        if angleDiff >= 0:
            return angleDiff
        else:
            return 360 - abs(angleDiff)

    def get_fieldAngle(self, deviceAngle, stageAngle):
        phi = stageAngle-deviceAngle
        phi_mod = phi % 360 # Module operator to map angle in range 0 to 360 deg
        return phi_mod

    # ____________________________________________________________________________
    # Do analysis
    def do(self):
        outputSubfolder = '/fittingOutput'
        if self.system == 'Berthold':
            OPFolder = self.IPFiles[0].fileDir.replace('DATA', 'ANALYSIS')+outputSubfolder
        else:
            OPFolder = self.IPFiles[0].fileDir + outputSubfolder
    
        if os.path.isdir(OPFolder) == False:
            os.makedirs(OPFolder)
    
        OPSummaryFileDir = OPFolder + '/000_fittingSummary.csv'
    
        if self.plotAllTogether != 0:
            fig2, ax2 = plt.subplots()
    
        fileParameter = {}
        sortedInputFileArray = np.array([])
        colorNumber = np.linspace(0, 1, len(self.IPFiles))
        colors = [ cm.gist_rainbow(x) for x in colorNumber ]
    
        if self.measMode == 1:
            fac = pd.read_csv(self.facFile.fileDirName)
            # print(fac)
    
    
        for IPFile in self.IPFiles:
            current, dBm, frequency, name, number = self.parametersFromFileName(IPFile.fileDirName, IPFile.fileDir)
            fileParameter[IPFile.fileDirName] = {"i" : current, "rf": dBm, "f": frequency, "name": name, "number" : number}
    
            if self.measMode == 1:
                stageAngle = float(fac.loc[fac['file_number']==number]['angle'].to_numpy())
                fieldAngle = self.get_fieldAngle(self.deviceAngle, stageAngle)
                voltageMultiplier = float(fac.loc[fac['file_number']==number]['V_polarity'].to_numpy()) * self.baseVoltageMultiplier
            else:
                stageAngle = self.stageAngle
                fieldAngle = self.get_fieldAngle(self.deviceAngle, self.stageAngle)
                voltageMultiplier = self.baseVoltageMultiplier
    
            fileParameter[IPFile.fileDirName]['deviceAngle'] = self.deviceAngle
            fileParameter[IPFile.fileDirName]['stageAngle'] = stageAngle
            fileParameter[IPFile.fileDirName]['fieldAngle'] = fieldAngle
            fileParameter[IPFile.fileDirName]['voltageMultiplier'] = voltageMultiplier
    
            positionFound = 0
    
            for index, fileName in enumerate(sortedInputFileArray):
    
                if fileParameter[fileName]["f"] < frequency and not(positionFound) :
    
                    positionFound = 1
    
                    leftArray = sortedInputFileArray[:index]
                    rightArray = sortedInputFileArray[index:]
                    sortedInputFileArray = np.append(leftArray, IPFile.fileDirName)
                    sortedInputFileArray = np.append(sortedInputFileArray, rightArray)
    
            if not(positionFound):
                sortedInputFileArray = np.append(sortedInputFileArray, IPFile.fileDirName)
    
        outFile1Cont = []
        for index, inputFileName in enumerate(sortedInputFileArray):
            def progress(i, array):
                print(f'Progress: {i+1}/{len(array)}')
    
            progress(index, sortedInputFileArray)
            print(f'InputFileName: {inputFileName}')
            voltageMultiplier = fileParameter[inputFileName]['voltageMultiplier']
            stageAngle = fileParameter[inputFileName]['stageAngle']
    
            IPFileObject = None
            for file in self.IPFiles:
                if file.fileDirName == inputFileName:
                    IPFileObject = file
    
            fieldArray = np.array([])
            amplitudeArray = np.array([])
            fileArray = np.genfromtxt(inputFileName, skip_header=self.numberOFHeaderLines, delimiter = '\t')
            hArray = fileArray[:, self.hAxis]
    
            vArray = fileArray[:, self.vAxis]*voltageMultiplier
    
            if self.plotANdCheck:
                #plot and check
                fig, ax = plt.subplots()
                legend = str(fileParameter[inputFileName]["i"]) + " mA " + str(fileParameter[inputFileName]["f"]) + " GHz"
                ax.plot( hArray, vArray, "ro", label = legend)
    
                ax.set_xlabel("Magnetic Field (Oe)")
                ax.set_ylabel(r"$V_{mix}$ (V)")
                ax.legend()
    
                plt.show()
                fig.savefig(inputFileName + "_rawData.png", bbox_inches="tight", dpi=self.dpi)
    
                minRange = float(input ("Type the MINIMUM field value for fitting\n"))
                maxRange = float(input ("Type the MAXIMUM field value for fitting\n"))
    
            else:
                minRange = stfmrRanges.minRange[fileParameter[inputFileName]["f"]]
                maxRange = stfmrRanges.maxRange[fileParameter[inputFileName]["f"]]
    
            for field, voltage in zip(hArray, vArray):
                if abs(field) >= abs(minRange) and abs(field) <= abs(maxRange):
                    fieldArray = np.append(fieldArray, field)
                    amplitudeArray = np.append(amplitudeArray, voltage)
    
            # Correct for offset field
            fieldArray -= self.offsetField
    
            # Do fitting:
            stfmrFit = stfmrSpectraFitting(fieldArray, amplitudeArray, fileParameter[inputFileName]["f"], inputFileName)
            
            # For summary file
            outFile1ContLine = {
                'Index': fileParameter[inputFileName]["number"],
                'stageAngle (deg)': fileParameter[inputFileName]["stageAngle"],
                'fieldAngle (deg)': fileParameter[inputFileName]["fieldAngle"],
                'Current (mA)': fileParameter[inputFileName]["i"],
                'rf Power (dBm)': fileParameter[inputFileName]["rf"],
                'Frequency (GHz)': fileParameter[inputFileName]["f"],
                'offsetField (Oe)': self.offsetField,
                'DeltaH (Oe)': stfmrFit.deltaH,
                'DeltaHError (Oe)': stfmrFit.deltaHErr,
                'Hres (Oe)': stfmrFit.Hres,
                'HresError (Oe)': stfmrFit.HresErr,
                'V0 (V)': stfmrFit.V0,
                'V0Error (V)': stfmrFit.V0Err,
                'V1 (V)': stfmrFit.V1,
                'V1Error (V)': stfmrFit.V1Err,
                'Vsym (V)': stfmrFit.Vsym,
                'VsymError (V)': stfmrFit.VsymErr,
                'Vas (V)': stfmrFit.Vas,
                'VasError (V)': stfmrFit.VasErr,
                }
            outFile1Cont.append(outFile1ContLine)
            
            
            # Lines for each file
            with open(OPFolder + '/' + IPFileObject.fileName + '_fit.csv', "w") as outFile2:
                outFile2.write("stageAngle: %s mA \n" %fileParameter[inputFileName]["i"])
                outFile2.write("Current: %s deg \n" %fileParameter[inputFileName]["stageAngle"])
                outFile2.write("Frequency: %.2f GHz \n" %fileParameter[inputFileName]["f"])
                outFile2.write("offsetField: %.2f Oe \n" %self.offsetField)
                outFile2.write("Field,Vmix,Vsym,Vas,V0,V1\n")
                outFile2.write("Oe,V,V,V,V,V\n")
                for field, vmix, vsym, vas, v0, v1 in zip(stfmrFit.h, stfmrFit.v, stfmrFit.symmetric, stfmrFit.asymmetric, stfmrFit.v0Array, stfmrFit.v1Array):
                    outFile2.write(str(field) + ',' + str(vmix) + ',')
                    outFile2.write(str(vsym) + ',' + str(vas) + ',')
                    outFile2.write(str(v0) + ',' + str(v1) + '\n')
    
    
    
    
            if self.legendMode == 0:
                legend = str(fileParameter[inputFileName]["f"]) + " GHz"
            elif self.legendMode == 1:
                legend = str(fileParameter[inputFileName]["i"]) + " mA"
            elif self.legendMode == 2:
                legend = str(fileParameter[inputFileName]["rf"]) + " dBm"
            elif self.legendMode == 3:
                legend = str(fileParameter[inputFileName]["f"]) + " GHz  " + str(fileParameter[inputFileName]["i"]) + " mA"
            else:
                legend = str(current) + " mA, " + str(frequency) + " GHz, " + str(dBm) + " dBm"
    
            legendFitting = r"$\Delta$H: %.2e $\pm$ %.1e" % (stfmrFit.deltaH, stfmrFit.deltaHErr) + "\n"
            legendFitting = legendFitting + r"$H_{res}$: %.2e $\pm$ %.1e" % (stfmrFit.Hres, stfmrFit.HresErr)
    
            if self.plotAllTogether != 1 :
                # print('cekc')
                # print(stfmrFit.h)
                fig, ax = plt.subplots()
                ax.plot(fieldArray, amplitudeArray, "ro", label = legend, markersize = self.mSize)
                ax.plot(stfmrFit.h, stfmrFit.v, label = legendFitting)
                ax.plot(stfmrFit.h, stfmrFit.Hres*stfmrFit.V1+ stfmrFit.V0 + stfmrFit.symmetric, label = r"$V_{sym}$: %.2e $\pm$  %.1e" % (stfmrFit.Vsym, stfmrFit.VsymErr))
                ax.plot(stfmrFit.h, stfmrFit.Hres*stfmrFit.V1+ stfmrFit.V0 + stfmrFit.asymmetric, label = r"$V_{as}$: %.2e $\pm$  %.1e" % (stfmrFit.Vas, stfmrFit.VasErr))
    
                ax.ticklabel_format(style = "sci", scilimits = (0,0), axis="y", useMathText = True)
    
                ax.ticklabel_format(style = "sci", scilimits = (0,0), axis="y", useMathText = True)
    
                ax.set_xlabel("Magnetic Field (Oe)")
                ax.set_ylabel(r"$V_{mix}$ (V)")
                ax.legend()
    
                fig.savefig(OPFolder + '/' + IPFileObject.fileNameWOExt + '_fit.png', bbox_inches="tight", dpi=self.dpi)
    
            if self.plotAllTogether != 0 :
                ax2.plot( fieldArray, amplitudeArray - stfmrFit.V0, "o", label = legend, color = colors[index], markersize = self.mSize)
                ax2.plot(stfmrFit.h, stfmrFit.v - stfmrFit.V0, "-b")
    
                ax2.ticklabel_format(style = "sci", scilimits = (0,0), axis="y", useMathText = True)
    
                ax2.set_xlabel("Magnetic Field (Oe)")
                ax2.set_ylabel(r"$V_{mix}$ (V)")
                ax2.legend()
                fig2.savefig( "toProcess/" + fileParameter[inputFileName]["name"] + "fitting_together.png", bbox_inches="tight", dpi=self.dpi)
                fig2.savefig( "toProcess/" + fileParameter[inputFileName]["name"] + "fitting_together.eps", bbox_inches="tight", dpi=self.dpi)
    
        self.outFile1Cont = pd.DataFrame(outFile1Cont)
        self.outFile1Cont.to_csv(OPSummaryFileDir, index=False)
        
        self.stfmrFit = stfmrFit
        plt.show()

























