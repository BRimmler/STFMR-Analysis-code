''' Credit: Avanindra Pandeya, MPI Halle, Python God <3 '''

import os
from Stfmr import Stfmr
import numpy as np
import matplotlib.pyplot as plt
import stfmrRange
from matplotlib import cm
import pandas as pd
# import tkinter as tk
# from tkinter import filedialog


def stfmrAnalysis(measMode, deviceAngle, IPFiles, facFile, plotANdCheck, plotAllTogether, legendMode, numberOFHeaderLines, hAxis, vAxis, baseVoltageMultiplier, mSize):

# ____________________________________________________________________________
# Subfunctions

    def parametersFromFileName(inputFileName, inputFileDir):
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

    def get_angleShifted(angle, deviceAngle):
        angleDiff = angle - deviceAngle
        if angleDiff >= 0:
            return angleDiff
        else:
            return 360 - abs(angleDiff)


# ____________________________________________________________________________



    outputSubfolder = '/fittingOutput'
    OPFolder = IPFiles[0].filedir + outputSubfolder
    if os.path.isdir(OPFolder) == False:
        os.makedirs(OPFolder)

    OPSummaryFileDir = OPFolder + '/000_fittingSummary.csv'


    with open(OPSummaryFileDir, "a") as outFile1:
        outFile1.write("Number,Angle (deg),AngleShifted (deg),Current (mA),rf Power (dBm),Frequency (GHz),DeltaH (Oe),DeltaHError (Oe),Hres (Oe),HresError (Oe),V0 (V),V0Error (V),V1 (V),V1Errror (V),Vsym (V),VsymError (V),Vas (V),VasError (V)\n")
        # outFile1.write(",deg,deg,mA,dBm,GHz,Oe,Oe,Oe,Oe,V,V,V,V,V,V,V,V\n")

    if plotAllTogether != 0:
        fig2, ax2 = plt.subplots()

    fileParameter = {}
    sortedInputFileArray = np.array([])
    colorNumber = np.linspace(0, 1, len(IPFiles))
    colors = [ cm.gist_rainbow(x) for x in colorNumber ]

    if measMode == 2:
        fac = pd.read_csv(facFile.file_fulldir)
        print(fac)


    for IPFile in IPFiles:
        current, dBm, frequency, name, number = parametersFromFileName(IPFile.file_fulldir, IPFile.filedir)
        fileParameter[IPFile.file_fulldir] = {"i" : current, "rf": dBm, "f": frequency, "name": name, "number" : number}

        if measMode == 2:
            angle = float(fac.loc[fac['file_number']==number]['angle'].to_numpy())
            angleShifted = get_angleShifted(angle, deviceAngle)
            voltageMultiplier = float(fac.loc[fac['file_number']==number]['V_polarity'].to_numpy()) * baseVoltageMultiplier
        else:
            angle = None
            angleShifted = None
            voltageMultiplier = baseVoltageMultiplier

        fileParameter[IPFile.file_fulldir]['angle'] = angle
        fileParameter[IPFile.file_fulldir]['angleShifted'] = angleShifted
        fileParameter[IPFile.file_fulldir]['voltageMultiplier'] = voltageMultiplier


        positionFound = 0

        for index, fileName in enumerate(sortedInputFileArray):

            if fileParameter[fileName]["f"] < frequency and not(positionFound) :

                positionFound = 1

                leftArray = sortedInputFileArray[:index]
                rightArray = sortedInputFileArray[index:]
                sortedInputFileArray = np.append(leftArray, IPFile.file_fulldir)
                sortedInputFileArray = np.append(sortedInputFileArray, rightArray)

        if not(positionFound):
            sortedInputFileArray = np.append(sortedInputFileArray, IPFile.file_fulldir)


    for index, inputFileName in enumerate(sortedInputFileArray):
        voltageMultiplier = fileParameter[inputFileName]['voltageMultiplier']
        angle = fileParameter[inputFileName]['angle']
        # print('voltageMultiplier = ', voltageMultiplier)
        # print('angle = ', angle)

        IPFileObject = None
        for file in IPFiles:
            if file.file_fulldir == inputFileName:
                IPFileObject = file

        fieldArray = np.array([])
        amplitudeArray = np.array([])
        fileArray = np.genfromtxt(inputFileName,  skip_header=numberOFHeaderLines, delimiter = '\t')
        hArray = fileArray[:, hAxis]

        vArray = fileArray[:, vAxis]*voltageMultiplier

        if plotANdCheck:
            #plot and check
            fig, ax = plt.subplots()
            legend = str(fileParameter[inputFileName]["i"]) + " mA " + str(fileParameter[inputFileName]["f"]) + " GHz"
            ax.plot( hArray, vArray, "ro", label = legend)

            ax.set_xlabel("Magnetic Field (Oe)")
            ax.set_ylabel(r"$V_{mix}$ (V)")
            ax.legend()

            plt.show()
            fig.savefig(inputFileName + "_rawData.png", bbox_inches="tight", dpi=600)

            minRange = float(input ("Type the MINIMUM field value for fitting\n"))
            maxRange = float(input ("Type the MAXIMUM field value for fitting\n"))

        else:
            minRange = stfmrRange.minRange[fileParameter[inputFileName]["f"]]
            maxRange = stfmrRange.maxRange[fileParameter[inputFileName]["f"]]

        for field, voltage in zip(hArray, vArray):
            if abs(field) >= abs(minRange) and abs(field) <= abs(maxRange):
                fieldArray = np.append(fieldArray, field)
                amplitudeArray = np.append(amplitudeArray, voltage)



        STFMR = Stfmr(fieldArray, amplitudeArray, fileParameter[inputFileName]["f"], inputFileName)
        with open(OPSummaryFileDir, "a") as outFile1:
            outFile1.write(str(fileParameter[inputFileName]["number"])+ ',')
            outFile1.write(str(fileParameter[inputFileName]["angle"])+ ',')
            outFile1.write(str(fileParameter[inputFileName]["angleShifted"])+ ',')
            outFile1.write( str(fileParameter[inputFileName]["i"]) +',' + str(fileParameter[inputFileName]["rf"]) +',' + str(fileParameter[inputFileName]["f"]) + ',')
            outFile1.write( str(STFMR.deltaH) +',' +  str(STFMR.deltaHErr) +',')
            outFile1.write( str(STFMR.Hres) +',' +  str(STFMR.HresErr) +',')
            outFile1.write( str(STFMR.V0) +',' +  str(STFMR.V0Err) +',')
            outFile1.write( str(STFMR.V1) +',' +  str(STFMR.V1Err) +',')
            outFile1.write(str(STFMR.Vsym) +',' +  str(STFMR.VsymErr) + ',')
            outFile1.write(str(STFMR.Vas) +',' +  str(STFMR.VasErr) + ',\n')
        with open(OPFolder + '/' + IPFileObject.filename_wo_ext + '_fit.csv', "w") as outFile2:
            outFile2.write("Angle: %s mA \n" %fileParameter[inputFileName]["i"])
            outFile2.write("Current: %s deg \n" %fileParameter[inputFileName]["angle"])
            outFile2.write("Frequency: %.2f GHz \n" %fileParameter[inputFileName]["f"])
            outFile2.write("Field,Vmix,Vsym,Vas,V0,V1\n")
            outFile2.write("Oe,V,V,V,V,V\n")
            for field, vmix, vsym, vas, v0, v1 in zip(STFMR.h, STFMR.v, STFMR.symmetric, STFMR.asymmetric, STFMR.v0Array, STFMR.v1Array):
                outFile2.write(str(field) + ',' + str(vmix) + ',')
                outFile2.write(str(vsym) + ',' + str(vas) + ',')
                outFile2.write(str(v0) + ',' + str(v1) + '\n')




        if legendMode == 0:
            legend = str(fileParameter[inputFileName]["f"]) + " GHz"
        elif legendMode == 1:
            legend = str(fileParameter[inputFileName]["i"]) + " mA"
        elif legendMode == 2:
            legend = str(fileParameter[inputFileName]["rf"]) + " dBm"
        elif legendMode == 3:
            legend = str(fileParameter[inputFileName]["f"]) + " GHz  " + str(fileParameter[inputFileName]["i"]) + " mA"
        else:
            legend = str(current) + " mA, " + str(frequency) + " GHz, " + str(dBm) + " dBm"

        legendFitting = r"$\Delta$H: %.2e $\pm$ %.1e" % (STFMR.deltaH, STFMR.deltaHErr) + "\n"
        legendFitting = legendFitting + r"$H_{res}$: %.2e $\pm$ %.1e" % (STFMR.Hres, STFMR.HresErr)

        if plotAllTogether != 1 :
            fig, ax = plt.subplots()
            ax.plot( fieldArray, amplitudeArray, "ro", label = legend, markersize = mSize)
            ax.plot(STFMR.h, STFMR.v, label = legendFitting)
            ax.plot(STFMR.h, STFMR.Hres*STFMR.V1+ STFMR.V0 + STFMR.symmetric, label = r"$V_{sym}$: %.2e $\pm$  %.1e" % (STFMR.Vsym, STFMR.VsymErr))
            ax.plot(STFMR.h, STFMR.Hres*STFMR.V1+ STFMR.V0 + STFMR.asymmetric, label = r"$V_{as}$: %.2e $\pm$  %.1e" % (STFMR.Vas, STFMR.VasErr))

            ax.ticklabel_format(style = "sci", scilimits = (0,0), axis="y", useMathText = True)

            ax.ticklabel_format(style = "sci", scilimits = (0,0), axis="y", useMathText = True)

            ax.set_xlabel("Magnetic Field (Oe)")
            ax.set_ylabel(r"$V_{mix}$ (V)")
            ax.legend()

            fig.savefig(OPFolder + '/' + IPFileObject.filename_wo_ext + '_fit.png', bbox_inches="tight", dpi=600)

        if plotAllTogether != 0 :
            ax2.plot( fieldArray, amplitudeArray - STFMR.V0, "o", label = legend, color = colors[index], markersize = mSize)
            ax2.plot(STFMR.h, STFMR.v - STFMR.V0, "-b")

            ax2.ticklabel_format(style = "sci", scilimits = (0,0), axis="y", useMathText = True)

            ax2.set_xlabel("Magnetic Field (Oe)")
            ax2.set_ylabel(r"$V_{mix}$ (V)")
            ax2.legend()
            fig2.savefig( "toProcess/" + fileParameter[inputFileName]["name"] + "fitting_together.png", bbox_inches="tight", dpi=600)
            fig2.savefig( "toProcess/" + fileParameter[inputFileName]["name"] + "fitting_together.eps", bbox_inches="tight", dpi=600)

    plt.show()

























