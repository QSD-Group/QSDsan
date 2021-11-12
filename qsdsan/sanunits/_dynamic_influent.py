# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from .. import SanUnit, WasteStream
from ..utils import ospath, load_data, data_path
import numpy as np, pandas as pd

__all__ = ('DynamicInfluent',)
dynamic_inf_path = ospath.join(data_path, 'sanunit_data/_inf_dry_2006.tsv')

class DynamicInfluent(SanUnit):
    
    _N_ins = 0
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), data_file=None, thermo=None,
                 init_with='WasteStream', isdynamic=True, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, isdynamic=isdynamic)
        self._init_from_file(dynamic_inf_path)
        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_from_file(self, file_path=None):
        path = file_path or dynamic_inf_path
        df = load_data(path, index_col=None)
        self._data = df
    
    #!!! use natural cubic spline to fit the data?
    
    def _run(self):
        pass