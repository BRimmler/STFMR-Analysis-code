import os
from Stfmr import Stfmr
import numpy as np
import matplotlib.pyplot as plt
import stfmrRange
from matplotlib import cm
import tkinter.filedialog
import tkinter as tk

""" This is the input zone. Here we will enter everything which needs to be entered manually """

# Are you not sure about the range of fitting. Set plotAndCheck to True to get a chance to review
# Tip: Put only a few representativie file in /toProcess folder before doing it. Unless you want to keep typing the whole day.
plotANdCheck = 0 # 1: Check each file for selcting the right range, 0: for analysis without checking

plotAllTogether = 0 # 0: Separate files for each curve, 1: Plot all curves together, 2: Separate files as well as together

legendMode = 3 # 0: Frequency, 1: Current, 2: rf power, 3: Current and Frequency, Anyother Number: All three written together

numberOFHeaderLines = 5

#File Format Inputs
hAxis = 0 #Unit Oe
vAxis = 1 #Unit V

voltageMmultiplier = 1 # It multiplies the voltageArray to, for example, flip the voltage

mSize = 3 # Size of markers in plot

root = tk.Tk()
root.withdraw()
path = tkinter.filedialog.askdirectory(parent=root, title='Choose a sample parameter file')


# path = "toProcess" # Path where all the files to be analyzed are kept

""" Input zone is over! Careful before changing anything unless you know what you are doing! """

nameFileArray = tuple(path + '/' + file for file in os.listdir(path) if file.endswith(".txt"))



def parametersFromFileName(inputFileName):
    """
    

    Parameters
    ----------
    inputFileName :String
        File Name

    Returns
    -------
    current : flaot
       Value of current from filename.
    dBm : flaot
        Value of power from filename..
    frequency : float
        Value of frequency from filename..

    """
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
    number = float((inputFileName[len(path)+1:len(path)+4]))
    return current, dBm, frequency, name, number

#Create a file to log all the fitting parameters
if not os.path.exists(path + '/parameter'):
    os.mkdir(path + '/parameter')
with open(path + '/parameter/FittingParam' + " .csv", "a") as outFile1:
    outFile1.write("Current, rf Power, Frequency,DeltaH,DeltaHError,Hres,HresError,V0,V0Error,V1,V1Errror,Vsym,VsymError,Vas,VasError,Number\n")
    outFile1.write("mA,dBm,GHz,Oe,Oe,Oe,Oe,V,V,V,V,V,V,V,V,\n")

#If all the files are to be plotted together, generate figure.
if plotAllTogether != 0:
    fig2, ax2 = plt.subplots()

#Variable to store all the parameters in the file
fileParameter = {}

#Variable to store the names of all the files in sorted order
sortedInputFileArray = np.array([])

#Generate number of colors of different plots based on number of files
colorNumber = np.linspace(0, 1, len(nameFileArray)) 

#Generate a color pallete 
colors = [ cm.gist_rainbow(x) for x in colorNumber ]

#Create a sorted array of all the files and store their parameters in the variable declared above
for inputFileName in nameFileArray:
    current, dBm, frequency, name, number = parametersFromFileName(inputFileName)
    fileParameter[inputFileName] = {"i" : current, "rf": dBm, "f": frequency, "name": name, "number" : number}
    
    positionFound = 0

    
    for index, fileName in enumerate(sortedInputFileArray):
        
        if fileParameter[fileName]["f"] < frequency and not(positionFound) :
            
            positionFound = 1
            
            leftArray = sortedInputFileArray[:index]
            rightArray = sortedInputFileArray[index:]
            sortedInputFileArray = np.append(leftArray, inputFileName)
            sortedInputFileArray = np.append(sortedInputFileArray, rightArray)
            
    if not(positionFound):
        sortedInputFileArray = np.append(sortedInputFileArray, inputFileName)

#Perform fittings and generate all the graphs
for index, inputFileName in enumerate(sortedInputFileArray):
    fieldArray = np.array([])
    amplitudeArray = np.array([])
    fileArray = np.genfromtxt(inputFileName,  skip_header=numberOFHeaderLines, delimiter = '\t')
    hArray = fileArray[:, hAxis]
    vArray = fileArray[:, vAxis]*voltageMmultiplier
    
    if plotANdCheck:
        #plot and check
        fig, ax = plt.subplots()
        legend = str(fileParameter[inputFileName]["i"]) + " mA " + str(fileParameter[inputFileName]["f"]) + " GHz" 
        ax.plot( hArray, vArray, "ro", label = legend)
        
        ax.set_xlabel("Magnetic Field (Oe)")
        ax.set_ylabel(r"$V_{mix}$ (V)")
        ax.legend()
        
        plt.show()
        fig.savefig(inputFileName + "_rawData.pdf", bbox_inches="tight", dpi=600)
        
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
    with open(path + '/parameter/FittingParam' + " .csv", "a") as outFile1:
        outFile1.write( str(fileParameter[inputFileName]["i"]) +',' + str(fileParameter[inputFileName]["rf"]) +',' + str(fileParameter[inputFileName]["f"]) + ',')
        outFile1.write( str(STFMR.deltaH) +',' +  str(STFMR.deltaHErr) +',')
        outFile1.write( str(STFMR.Hres) +',' +  str(STFMR.HresErr) +',')
        outFile1.write( str(STFMR.V0) +',' +  str(STFMR.V0Err) +',')
        outFile1.write( str(STFMR.V1) +',' +  str(STFMR.V1Err) +',')
        outFile1.write(str(STFMR.Vsym) +',' +  str(STFMR.VsymErr) + ',')
        outFile1.write(str(STFMR.Vas) +',' +  str(STFMR.VasErr) + ',')
        outFile1.write(str(fileParameter[inputFileName]["number"])+ '\n')
    with open(inputFileName + 'fititng.csv', "w") as outFile2:
        outFile2.write("Cuurent: %s mA \n" %fileParameter[inputFileName]["i"])
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
    
        fig.savefig(inputFileName + "fitting.pdf", bbox_inches="tight", dpi=600)
        
    if plotAllTogether != 0 :
        ax2.plot( fieldArray, amplitudeArray - STFMR.V0, "o", label = legend, color = colors[index], markersize = mSize)
        ax2.plot(STFMR.h, STFMR.v - STFMR.V0, "-b")
        
        ax2.ticklabel_format(style = "sci", scilimits = (0,0), axis="y", useMathText = True)
        
        ax2.set_xlabel("Magnetic Field (Oe)")
        ax2.set_ylabel(r"$V_{mix}$ (V)")
        ax2.legend()
        fig2.savefig( "toProcess/" + fileParameter[inputFileName]["name"] + "fitting_together.pdf", bbox_inches="tight", dpi=600)
        fig2.savefig( "toProcess/" + fileParameter[inputFileName]["name"] + "fitting_together.eps", bbox_inches="tight", dpi=600)

plt.show()

        
    
