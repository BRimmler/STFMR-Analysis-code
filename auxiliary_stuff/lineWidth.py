# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 15:38:45 2020

@author: apandeya
"""
import numpy as np
from scipy.optimize import curve_fit

def line(x, a, b):
        """

        Parameters
        ----------
        x : Array.
        a : slope.
        b : Intercept.

        Returns
        -------
        TYPE
            Array.

        """
        return a * x + b

class lineWidth():
    def __init__(self,currentArray, deltaHArray, deltaHErrorArray, frequency, hRes, fileName):
        """
        

        Parameters
        ----------
        currentArray : Array
        deltaHArray : Array
        deltaHErrorArray : Array
        frequency : float
        hRes : float
        fileName :string

        Returns
        -------
        None.

        """
        self.currentArray = np.array(currentArray)
        self.deltaHArray = np.array(deltaHArray)
        self.deltaHErrorArray = np.array(deltaHErrorArray)
        self.frequency = frequency
        self.hRes = hRes
        self.fileName = fileName
        
        for current, deltaH, deltaHError in zip(self.currentArray, self.deltaHArray, self.deltaHErrorArray ):
            if current == 0:
                self.deltaHAtZero = deltaH
                self.deltaHErrorAtZero = deltaHError
        
        self.deltaHArray = self.deltaHArray - self.deltaHAtZero
        
    def __repr__(self):
        return (self.fileName + "_f:" + str(self.frequency) + "_H:" + str(self.hRes))

    def fitLine(self):
        """
        
        Returns
        -------
        None. Creates two more array xArray and yArray for the object and calculates slope and intecept for the given inputs

        """
        popt, pcov = curve_fit(line, self.currentArray , self.deltaHArray, sigma=self.deltaHErrorArray)
        
        self.slope = popt[0]
        self.intercept = popt[1]
        
        self.slopeError = np.sqrt(pcov[0,0])
        self.interceptError = np.sqrt(pcov[1,1])
        
        
        
        self.yArray = line(self.currentArray, self.slope, self.intercept)
        
    
        
        