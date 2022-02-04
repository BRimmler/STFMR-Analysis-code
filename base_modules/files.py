# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 11:15:55 2021

@author: rimmler
"""

import os

# Generic file class
class File:
    def __init__(self, *fileLoc):
        if len(fileLoc) == 1: # Single string file directory and name
            fileLoc = fileLoc[0].replace('\\', '/')
            self.fileDirName = fileLoc
            self.fileDir = '/'.join(fileLoc.split('/')[:-1])
            self.fileName = fileLoc.split('/')[-1]
        elif len(fileLoc) == 2: # File directory and file name
            self.fileDir = fileLoc[0]
            self.fileName = fileLoc[1]
            self.fileDirName = self.fileDir + '/' + self.fileName
        else:
            raise SyntaxError('fileLoc must be of length 1 or 2.')

        self.fileExt = '.' + self.fileName.split('.')[-1]
        self.fileNameWOExt = '.'.join(self.fileName.split('.')[:-1])

    def makeDirIfNotExist(self):
        if os.path.isdir(self.fileDir) == False:
            os.makedirs(self.fileDir)
            print(self.fileDir+' created.')

    def create_subDir(self, name):
        subDir = self.firDir + '/' + name
        if os.path.isdir(subDir) is False:
            os.makedirs(subDir)
            print(subDir+' created.')
        return subDir




