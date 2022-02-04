# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 15:03:11 2022

@author: Berthold Rimmler
"""

class Param:
    def __init__(self, arg, val, unit=None, err=None):
        if not isinstance(arg, str):
            raise TypeError

        self.arg = arg
        self.val = float(val)
        if unit == 1 or unit is None:
            self.unit_rep = ''
        else:
            self.unit_rep = str(unit)
        self.err = err
        self.unit = unit

    def __str__(self):
        if self.err is None:
            rep = f'Param {self.arg} = {self.val} {self.unit_rep}'
        else:
            rep = f'Param {self.arg} = ({self.val}$\pm${self.err} {self.unit_rep}'
        return rep

# a = Param('a', 1009000)
# print(a)