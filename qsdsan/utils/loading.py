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
from .. import _pk

__all__ = ('load_data', 'data_path', 'dct_from_str', 'save_pickle', 'load_pickle')


def load_data(path=None, sheet=None, index_col=0, **kwargs):
    '''For data importing.'''
    last_4 = path[-4:]
    if last_4 in ('.tsv', '.txt'):
        data = pd.read_csv(path, sep='\t', index_col=index_col, **kwargs)
    elif last_4 == '.csv':
        data = pd.read_csv(path, index_col=index_col, **kwargs)
    elif last_4 in ('xlsx', '.xls'):
        try: data = pd.read_excel(path, sheet_name=sheet, engine='openpyxl',
                                  index_col=index_col, **kwargs)
        except: data = pd.read_excel(path, engine='openpyxl',
                                     index_col=index_col, **kwargs)
    else:
        raise ValueError('Only tsv, csv, or xlsx files can be loaded.')
    return data


def dct_from_str(dct_str, sep=',', dtype='float'):
    '''
    Use to parse str into a dict,
    the str should be written in `k1=v1, k2=v2, ..., kn=vn`
    (separated by comma or other symbols defined by `sep`).
    '''
    splitted = [i.split('=') for i in dct_str.replace(' ', '').split(sep)]

    if dtype == 'float':
        return {k:float(v) for k, v in splitted}

    elif dtype == 'int':
        return {k:int(v) for k, v in splitted}

    else:
        return {k:v for k, v in splitted}


def save_pickle(obj, path):
    '''Save object as a pickle file using Pickle Protocol 5.'''
    if _pk is None:
        raise RuntimeError('Current environment does not support Pickle Protocol 5, '
                           'cannot save pickle files.')
    f = open(path, 'wb')
    _pk.dump(obj, f, protocol=5)
    f.close()


def load_pickle(path):
    '''Load object saved as a pickle file using Pickle Protocol 5.'''
    if _pk is None:
        raise RuntimeError('Current environment does not support Pickle Protocol 5, '
                           'cannot load pickle files.')
    f = open(path, 'rb')
    obj = _pk.load(f)
    f.close()
    return obj