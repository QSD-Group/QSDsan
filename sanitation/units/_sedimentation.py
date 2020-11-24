#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Sanitation Explorer: Sustainable design of non-sewered sanitation technologies
Copyright (C) 2020, Sanitation Explorer Development Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-for-WaSH/sanitation/blob/master/LICENSE.txt
for license details.

Ref:
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.

'''


# %%

import numpy as np
from .. import SanUnit
from ..utils.loading import load_data, data_path

__all__ = ('Sedimentation',)

data_path += 'unit_data/Sedimentation.csv'

#!!! Potentially make these a Process, or make a Treatment class
from ._toilet import Toilet
_allocate_N_reduction = Toilet._allocate_N_reduction


class Sedimentation(SanUnit):
    '''Sedimentation of wastes into liquid and solid phases.'''
    
    def __init__(self, ID='', ins=None, outs=(),
                 if_CH4_captured=True, if_N_degradation=True,
                 **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs)

        data = load_data(path=data_path)
        for para in data.index:
            if para == 'split':
                value = eval(data.loc[para]['expected'])
                setattr(self, para, value)
            else:
                value = float(data.loc[para]['expected'])
                setattr(self, '_'+para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)

    @property
    def tau(self):
        '''[float] Residence time, [d].'''
        return self._tau
    @tau.setter
    def tau(self, i):
        self._tau = float(i)

    @property
    def split(self):
        '''
        [float] or [dict] Fractions of material retention in the settled solids
        before degradation. If a single number is provided, then it is assumed
        that retentions of all Components in the WasteStream are the same.
        
        Note
        ----
            Set state variable values (e.g., COD) will be retained if the retention
            ratio is a single number (treated like the loss stream is split
            from the original stream), but not when the ratio is a dict.

        '''
        return self._split
    @split.setter
    def split(self, i):
        try:
            self._split = float(i)
            self._split_type = 'float'
        except:
            if isinstance(i, dict):
                self._split = i
                self._split = 'dict'
            else:
                raise TypeError(f'Only float or dict allowed, not {type(i).__name__}.')

















