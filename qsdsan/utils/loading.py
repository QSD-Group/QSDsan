#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

import os
path = os.path.dirname(os.path.realpath(__file__))
data_path = path[:-6] + '/data/'
del os

import pandas as pd

__all__ = ('load_data', 'data_path')
 

def load_data(path=None, sheet=None):
    if path[-4:] == 'xlsx' or path[-4:] == '.xls':
        try: data = pd.read_excel(path, sheet_name=sheet, index_col=0, engine='openpyxl')
        except: data = pd.read_excel(path, index_col=0, engine='openpyxl')
    elif path[-4:] == '.csv':
        data = pd.read_csv(path, index_col=0)
    else:
        raise ValueError('Only csv or xlsx files can be loaded.')
    return data