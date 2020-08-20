import numpy as np
from scipy.optimize import curve_fit
import logging



class Stfmr():
    def __init__(self, hArray, vArray, frequency, inputFileName):
        """
        

        Parameters
        ----------
        hArray : Array of magnetic Field
        vArray : Array of voltages
        frequency : Float
        inputFileName : String
        Returns
        -------
        None.

        """
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
        """
        Negative value of DeltaH does not make physicsl sense.
        This function can flip the the sign of delta H in case fitting gives a negative value. 

        Returns
        -------
        None.

        """
        if self.deltaH < 0:
            self.Vas = -self.Vas
            self.deltaH = -self.deltaH
            
    def assignFieldAndAmpArray(self, hArray, vArray, minRange, maxRange):
        """
        

        Parameters
        ----------
        hArray : Array
            Magnetic Field
        vArray : Array
            Volatge.
        minRange : float
            Minimum value of magnetic field for which fitting needs to be done.
        maxRange : float
            Maximum value of magnetic field for which fitting needs to be done.

        Returns
        -------
        fieldArray : Array
            Magneitc Field.
        amplitudeArray : Array
            Volatge.
        This fuinction would select data in the given range

        """
        fieldArray = np.array([])
        amplitudeArray = np.array([])
        for field, voltage in zip(hArray, vArray):
            if abs(field) >= abs(minRange) and abs(field) <= abs(maxRange):
                fieldArray = np.append(fieldArray, field)
                amplitudeArray = np.append(amplitudeArray, voltage)
        return fieldArray, amplitudeArray
                
    def denominator(self, H, deltaH, Hres):
        """

        Parameters
        ----------
        H : flaot
            Magnetic Field.
        deltaH : float
            Linewidth.
        Hres : float
            Resonance Field.

        Returns
        -------
        TYPE
            Denominator of Lorentzian function .

        """
        return np.square(deltaH) + np.square(H - Hres)

    def symmetricFunction(self, H, deltaH, Hres):
        """
    
        Parameters
        ----------
        H : flaot
            Magnetic Field.
        deltaH : float
            Linewidth.
        Hres : float
            Resonance Field.

        Returns
        -------
        float
            Symmetirc Part of Lorentzian function .

        """
        return np.square(deltaH)/self.denominator(H, deltaH, Hres)

    def antiSymmetricFunction(self, H, deltaH, Hres):
        """
    
        Parameters
        ----------
        H : flaot
            Magnetic Field.
        deltaH : float
            Linewidth.
        Hres : float
            Resonance Field.

        Returns
        -------
        float
            Anti Symmetirc Part of Lorentzian function .

        """
        return (deltaH*(H - Hres))/self.denominator(H, deltaH, Hres)

    def vMix(self, H, deltaH, Hres, V0, V1, Vsym, Vas):
        """
        

        Parameters
        ----------
        H : flaot
            Magnetic Field.
        deltaH : float
            Linewidth.
        Hres : float
            Resonance Field.
        V0 : float
            DC part of the mixing voltage.
        V1 : float
            Linear part of the mixing voltage
        Vsym : float
            Symmetric component of the mixing voltage.
        Vas : float
            Anti Symmetric component of the mixing voltage..

        Returns
        -------
        float
        Mixing Voltage

        """
        return V0 + V1*H + Vsym*self.symmetricFunction(H, deltaH, Hres) + Vas*self.antiSymmetricFunction(H, deltaH, Hres)

    def fitFunction(self):
        """
        Create different arrays after fitting

        Returns
        -------
        None.

        """
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
        """
        

        Parameters
        ----------
        x : Array
            Array of dependent variables (magnetic field in this case).
        y : Array
            Array of independent variables (mixing voltage in this case).

        Returns
        -------
        TYPE
            float: Linewidth.

        """
        for a, b in zip(x, y):
            if b == np.amax(y):
                amax = a

            if b == np.amin(y):
                amin  = a

        return amax - amin
    

    def fitCurve(self):
        """
        This function performs actual fitting

        Returns
        -------
        None.

        """
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
        """
        

        Parameters
        ----------
        arrayx : Array
            Logs all the fitting parameters given in the array in a file.

        Returns
        -------
        None.

        """
        logging.info(self.addStringFloat("DeltaH: ", arrayx[0]))
        logging.info(self.addStringFloat("Hres: ", arrayx[1]))
        logging.info(self.addStringFloat("V0: ", arrayx[2]))
        logging.info(self.addStringFloat("V1: ", arrayx[3]))
        logging.info(self.addStringFloat("Vsym: ", arrayx[4]))
        logging.info(self.addStringFloat("Vas: ", arrayx[5]))

    
    def addStringFloat(self, stringvalue, floatValue):
        """
        

        Parameters
        ----------
        stringvalue : string
        floatValue : float

        Returns
        -------
        String
            Concatenates the two values.

        """
        return stringvalue + str(floatValue)
    
    def parametersFromFileName(self, inputFileName):
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
        return current, dBm, frequency