# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 16:24:19 2021

@author: rimmler
"""

import pandas as pd
from openpyxl import load_workbook

# _____________________________________________________________________________
# File handling
def excel_to_pandas(excel_file, sheet_name='Sheet1'):
    '''


    Parameters
    ----------
    excel_file : file.File
        An input excel file to read.
    sheet_name : str, optional
        The name of the excel sheet to read from. The default is 'Sheet1'.

    Returns
    -------
    df : pd.DataFrame
        The content of the excel sheet as pandas DataFrame.

    '''
    wb = load_workbook(excel_file.fileDirName)
    sheet = wb[sheet_name]
    df = pd.DataFrame(sheet.values)
    df.columns = df.iloc[0]
    df = df[1:]
    df.head()
    return df


# _____________________________________________________________________________
# File content handling
def clean_str(string):
    '''


    Parameters
    ----------
    string : str
        String with or without undesired "\n".

    Returns
    -------
    out : str
        Input without "\n" if present.

    '''
    if string[-1:] == '\n':
        out = string[:-1]
    else:
        out = string
    return out

def clean(obj):
    '''


    Parameters
    ----------
    obj : str, dict, array-type
        An str or a type of array of strings with or without undesired "\n".

    Returns
    -------
    out : as input type
        Input without "\n" present.

    '''
    if isinstance(obj, str):
        out = clean_str(obj)
    elif isinstance(obj, dict):
        out = {}
        for key in obj:
            value_in = obj[key]
            value_out = clean_str(value_in)
            out[key] = value_out
    elif isinstance(obj, list):
        out = []
        for item in obj:
            out.append(clean_str(item))
    return out

# _____________________________________________________________________________
# Python object handling
def get_lst_or_float(obj, i):
    '''


    Parameters
    ----------
    obj : float or list or tuple of floats
        A float or a list or tuple of floats.
    i : int
        Index of the float to be returned.

    Returns
    -------
    out : float
        If input is float --> returns input. Else returns value at given index.

    '''
    if isinstance(obj, list) or isinstance(obj, tuple):
        out = obj[i]
    else:
        out = obj
    return out


def make_lst(obj):
    '''


    Parameters
    ----------
    obj : An array-type object or an object
        The generic input is made into a list, regardless of the input type.

    Returns
    -------
    out : list
        List of inputs.

    '''
    if isinstance(obj, list) or isinstance(obj, tuple):
        out = obj
    elif isinstance(obj, str) or isinstance(obj, int) or isinstance(obj, float):
        out = [obj]
    return out


def pd_find_closest(df, val):
    '''


    Parameters
    ----------
    df : pd.Series
        A pandas Series.
    val : float
        A target value.

    Returns
    -------
    index : int
        The index of the element in df closest to val.
    val : TYPE
        The value of element in df closest to val.

    '''
    sort = df.iloc[(df-val).abs().argsort()[:1]]
    index = sort.index.tolist()[0]
    val = sort.tolist()[0]
    return index, val


# _____________________________________________________________________________
# Pandas
class Param:
    def __init__(self, arg, val, unit, verbose=True):
        if not isinstance(arg, str):
            if verbose is True: print(f'arg "{arg}" converted to str')
            self.arg = str(arg)
        else:
            self.arg = arg
        if not isinstance(unit, str):
            if verbose is True: print(f'unit "{unit}" converted to str')
            self.unit = str(unit)
        else:
            self.unit = unit
        if not isinstance(val, float):
            if verbose is True: print(f'val "{val}" converted to float')
            self.val = float(val)
        else:
            self.val = val

def item_to_param(key, val):
    key_split = key.split(' ')
    if len(key_split) > 2:
        raise ValueError('Simple conversion requires key as type "parameter (unit)"')
    arg, unit = key_split
    return Param(arg, float(val), unit)

def add_constant_column(df, name, val):
    # adding column with constant value
    df[name] = pd.Series([val for x in range(len(df.index))])
    return df
    
def shift_column_to_pos(df, col_name, loc):
    col = df[col_name]
    df = df.drop(columns=[col_name])
    df.insert(loc=loc, column=col_name, value=col)
    return df
    
def shift_column_to_start(df, col_name):
    df = shift_column_to_pos(df, col_name, 0)
    return df

# _____________________________________________________________________________
# SI and CGS system (who the fuck invented this shit??!!)





















