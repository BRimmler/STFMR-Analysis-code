# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 10:09:41 2022

@author: rimmler
"""

import shutil, errno

src = r'D:\owncloud\0_Personal\ANALYSIS\python\analysis-modules'
fileNames = ['files.py', 'helpers.py', 'parameters.py', 'plots.py', 'plots_settings.py']
trg = r'D:\owncloud\0_Personal\ANALYSIS\Mn3SnN\ST-FMR\Code\stfmr\base_modules'

def get(src, files, dst):
    for fileName in fileNames:
        fileDirName = src + '/' + fileName
        try:
            shutil.copy(fileDirName, dst)
        except OSError as exc: # python >2.5
            if exc.errno in (errno.ENOTDIR, errno.EINVAL):
                shutil.copy(fileDirName, dst)
            else: raise

print(f'Copy base modules from from {src} to {trg}')
get(src, fileNames, trg)
print('Done.')