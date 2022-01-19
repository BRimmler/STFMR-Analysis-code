# -*- coding: utf-8 -*-
"""
Created on Fri May 14 11:13:20 2021

@author: rimmler
"""
import os

class File:
    def __init__(self, file):
        file = file.replace('\\', '/')
        self.file_fulldir = file
        self.filedir = '/'.join(file.split('/')[:-1])
        self.filename = file.split('/')[-1]
        self.fileext = '.' + self.filename.split('.')[-1]
        self.filename_wo_ext = '.'.join(self.filename.split('.')[:-1])


    def makeDirIfNotExist(self):
        if os.path.isdir(self.filedir) == False:
            os.makedirs(self.filedir)
            print(self.filedir+' created.')

    def to_File2(self):
        return File2(self.filedir, self.filename)

def make_file(filedir, filename, fileext=''):
    file_fulldir = filedir + '/' + filename + fileext
    return File(file_fulldir)

class File2:
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

