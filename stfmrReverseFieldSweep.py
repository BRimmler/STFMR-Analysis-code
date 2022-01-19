# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 09:49:14 2021

@author: rimmler

Used to artificially reverse the field sweep direction of ST-FMR measurements.
This is required because part of the main analysis code does not work for field
sweep from low field to high field.

"""

inputDir = r'D:\DATA\Mn3SnN\ST-FMR\MA2959-2\210805\003_angle-dependence_D1'
inputExt = '.txt'
outputDir = inputDir + r'\reversed_field_sweep'




import os

class File:
    def __init__(self, fileDir, fileName):
        self.fileDir = fileDir
        self.fileName = fileName
        self.fileDirName = fileDir + '/' + fileName
        self.fileExt = '.' + self.fileName.split('.')[-1]
        self.fileNameWOExt = '.'.join(self.fileName.split('.')[:-1])

    def makeDirIfNotExist(self):
        if os.path.isdir(self.fileDir) == False:
            os.makedirs(self.fileDir)
            print(self.fileDir+' created.')

ipFileNames = os.listdir(inputDir)
ipFiles = []
for ipFileName in ipFileNames:
    if os.path.isfile(os.path.join(inputDir, ipFileName)):
        ipFiles.append(File(inputDir, ipFileName))

for ipFile in ipFiles:
    if ipFile.fileExt == inputExt:
        with open(ipFile.fileDirName) as f:
            ipContent = f.readlines()
        opContentWOHead = ipContent[3:]
        opContentWOHead.reverse()
        opContent = ipContent[:3] + opContentWOHead

        opFile = File(outputDir, ipFile.fileName)
        opFile.makeDirIfNotExist()

        with open(opFile.fileDirName, 'w') as f:
            opContent = f.writelines(opContent)

