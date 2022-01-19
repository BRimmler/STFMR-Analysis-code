# -*- coding: utf-8 -*-
"""
Created on Fri May 14 11:11:03 2021

@author: rimmler
"""

import sys
from PIL import Image
import stfmrHelpers as aux
import os
import tkinter as tk
from tkinter import filedialog

aspect_ratio = (5, 3)

# def ui_get_spectra():
root = tk.Tk()
root.withdraw()
spectra = filedialog.askopenfilenames(parent=root, title='Chose spectra .png files')

spectra_files = []
for spec in spectra:
    spectra_files.append(aux.File(spec))

N = len(spectra)

images = [Image.open(x) for x in spectra]
widths, heights = zip(*(i.size for i in images))

total_width = sum(widths)
max_height = max(heights)

new_im = Image.new('RGB', (total_width, max_height))

x_offset = 0
for im in images:
  new_im.paste(im, (x_offset,0))
  x_offset += im.size[0]

op_file = aux.make_file(spectra_files[0].filedir, '000SpectraMerged.jpg')

new_im.save(op_file.file_fulldir)