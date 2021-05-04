import numpy as np
from scipy.optimize import curve_fit
import logging



class Stfmr():
    def __init__(self, hArray, vArray, frequency, inputFileName):
        self.fieldArray = hArray
        self.amplitudeArray = vArray
        #self.fieldArray, self.amplitudeArray  = self.assignFieldAndAmpArray(hArray, vArray, minRange, maxRange)
        self.deltaH = self.findDh(self.fieldArray, self.amplitudeArray) #Units: T
        self.Hres = self.fieldArray[int(len(self.fieldArray)/2)] #Units: T
        self.frequency = frequency #Units : Hz
        self.V0 = (self.amplitudeArray[0] + self.amplitudeArray[-1])/2
        self.V1 = 0
        self.Vsym = np.amax(self.amplitudeArray) + np.amin(self.amplitudeArray) - 2*self.V0
        self.Vas = np.amax(self.amplitudeArray) - np.amin(self.amplitudeArray)
        
        self.current, self.dBm, self.frequency = self.parametersFromFileName(inputFileName)
        self.h = np.array([])
        self.v = np.array([])
        self.symmetric = np.array([])
        self.asymmetric = np.array([])
        self.v1Array = np.array([])
        self.v0Array = np.array([])

        logging.basicConfig(filename= inputFileName + 'fitting.log',level=logging.INFO)

        self.params = [self.deltaH, self.Hres, self.V0, self.V1, self.Vsym,self.Vas ]
        logging.info("------------------------------------")
        logging.info(inputFileName)
        logging.info("------------------------------------")
        logging.info("Initialization of Parameters")
        logging.info("------------------------------------")
        self.logparameters(self.params)
        self.fitCurve()
        self.correctVasSign()
        self.fitFunction()

    def __del__(self):
        logging.shutdown()
    
    def correctVasSign(self):
        if self.deltaH < 0:
            self.Vas = -self.Vas
            self.deltaH = -self.deltaH
            
    def assignFieldAndAmpArray(self, hArray, vArray, minRange, maxRange):
        fieldArray = np.array([])
        amplitudeArray = np.array([])
        for field, voltage in zip(hArray, vArray):
            if abs(field) >= abs(minRange) and abs(field) <= abs(maxRange):
                fieldArray = np.append(fieldArray, field)
                amplitudeArray = np.append(amplitudeArray, voltage)
        return fieldArray, amplitudeArray
                
    def denominator(self, H, deltaH, Hres):
        return np.square(deltaH) + np.square(H - Hres)

    def symmetricFunction(self, H, deltaH, Hres):
        return np.square(deltaH)/self.denominator(H, deltaH, Hres)

    def antiSymmetricFunction(self, H, deltaH, Hres):
        return (deltaH*(H - Hres))/self.denominator(H, deltaH, Hres)

    def vMix(self, H, deltaH, Hres, V0, V1, Vsym, Vas):
        return V0 + V1*H + Vsym*self.symmetricFunction(H, deltaH, Hres) + Vas*self.antiSymmetricFunction(H, deltaH, Hres)

    def fitFunction(self):
        h = self.fieldArray[-1]
        while (abs(h) < abs(self.fieldArray[0])):
            self.h = np.append(self.h, h)
            self.v = np.append(self.v, self.vMix(h, self.deltaH, self.Hres, self.V0, self.V1, self.Vsym,self.Vas))
            self.symmetric = np.append(self.symmetric, self.Vsym*self.symmetricFunction(h, self.deltaH, self.Hres))
            self.asymmetric = np.append(self.asymmetric, self.Vas*self.antiSymmetricFunction(h, self.deltaH, self.Hres))
            self.v0Array = np.append(self.v0Array, self.V0)
            self.v1Array = np.append(self.v1Array, h*self.V1)
            h += (self.fieldArray[0]-self.fieldArray[-1])/abs(self.fieldArray[0]-self.fieldArray[-1])
        
        logging.debug("self.h: " + str(self.h))

        logging.debug("self.v: " + str(self.v))
    
    def findDh(self, x, y):
        for a, b in zip(x, y):
            if b == np.amax(y):
                amax = a

            if b == np.amin(y):
                amin  = a

        return amax - amin
    
#    def vMixFun(self, field):
#        
#        return self.vMix(field, self.deltaH, self.Hres, self.V0, self.V1, self.Vsym, self.Vas)
#    
#    def symmFun(self, field):
#        
#        return self.Vsym*self.symmetricFunction(field, self.deltaH, self.Hres)
#    
#    def antisymmFun(self, field):
#        
#        return self.Vas*self.antiSymmetricFunction(field, self.deltaH, self.Hres)



    def fitCurve(self):
        self.fitParams, self.fitConv = curve_fit(self.vMix, self.fieldArray, self.amplitudeArray, self.params, method = 'lm', ftol = 1e-10)
        self.deltaH = self.fitParams[0]
        self.Hres = self.fitParams[1]
        self.V0 = self.fitParams[2]
        self.V1 = self.fitParams[3]
        self.Vsym = self.fitParams[4]
        self.Vas = self.fitParams[5]
        
        
        self.deltaHErr = np.sqrt(self.fitConv[0,0])
        self.HresErr = np.sqrt(self.fitConv[1,1])
        self.V0Err = np.sqrt(self.fitConv[2,2])
        self.V1Err = np.sqrt(self.fitConv[3,3])
        self.VsymErr = np.sqrt(self.fitConv[4,4])
        self.VasErr = np.sqrt(self.fitConv[5,5])
        
        
        
        
        logging.info("------------------------------------")
        logging.info("After Fitting")
        logging.info("------------------------------------")
        self.logparameters(self.fitParams)
        logging.info("---------------------THE END--------------------------------")
        logging.info("\n")

    def logparameters(self, arrayx):
        logging.info(self.addStringFloat("DeltaH: ", arrayx[0]))
        logging.info(self.addStringFloat("Hres: ", arrayx[1]))
        logging.info(self.addStringFloat("V0: ", arrayx[2]))
        logging.info(self.addStringFloat("V1: ", arrayx[3]))
        logging.info(self.addStringFloat("Vsym: ", arrayx[4]))
        logging.info(self.addStringFloat("Vas: ", arrayx[5]))

    
    def addStringFloat(self, stringvalue, floatValue):
        return stringvalue + str(floatValue)
    
    def parametersFromFileName(self, inputFileName):
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
        return current, dBm, frequency