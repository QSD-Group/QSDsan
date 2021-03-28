#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import os
path = os.path.dirname(os.path.realpath(__file__))
data_path = path[:-6] + '/data/'
del os

import pandas as pd

__all__ = ('load_data', 'data_path')
 

def load_data(path=None, sheet=None, **kwargs):
    kwargs.setdefault('index_col', 0)
    last_4 = path[-4:]
    if last_4 == '.tsv':
        data = pd.read_csv(path, sep='\t', **kwargs)
    elif last_4 == '.csv':
        data = pd.read_csv(path, **kwargs)
    elif last_4 in ('xlsx', '.xls'):
        try: data = pd.read_excel(path, sheet_name=sheet, engine='openpyxl', **kwargs)
        except: data = pd.read_excel(path, engine='openpyxl', **kwargs)
    else:
        raise ValueError('Only tsv, csv, or xlsx files can be loaded.')
    return data