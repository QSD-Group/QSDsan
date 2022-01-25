#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 09:53:41 2021

@author: lewisrowles
"""

"""
Created on Tue Mar 23 09:07:43 2021

@author: lewisrowles
"""

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Lewis Rowles <stetsonsc@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''
# %%

from collections.abc import Iterable
from biosteam._graphics import UnitGraphics
import numpy as np
from .. import SanUnit, Construction
from ..utils import load_data, data_path

__all__ = ('ConstructionBiogenicRefinery',)

#path to csv with all the inputs

#data_path = '/Users/lewisrowles/opt/anaconda3/lib/python3.8/site-packages/exposan/biogenic_refinery/_housing_biogenic_refinery.csv'


### 
class ConstructionBiogenicRefinery(SanUnit, isabstract=True):
    '''
    ConstructionBiogenicRefinery. 

    
    Reference documents
    -------------------
    N/A
    
    Parameters
    ----------
    ins : none
    outs : none

        
    References
    ----------
    .. N/A
    
    '''
    
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, 
                 init_with='WasteStream',const_wage=15, const_person_days=10000000000,
                 **kwargs):
        # if isinstance(ins, str) or (not isinstance(ins, Iterable)):
        #     self._N_outs = self._N_ins = 1
        # else:
        #     self._N_outs = self._N_ins = len(ins)
        self._graphics = UnitGraphics.box(self._N_ins, self._N_outs)
        SanUnit.__init__(self, ID, ins, outs)
        self.const_wage = const_wage
        self.const_person_days = const_person_days

        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    # def _run(self):
    #     waste = self.ins[0]
    #     treated = self.outs[0]
    #     treated.copy_like(self.ins[0])
    _N_ins=0
    _N_outs=0

        
        
 
    #_cost based on amount of steel and stainless plus individual components
    def _cost(self):
        
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        self.baseline_purchase_costs['Labor'] = (self.const_wage*self.const_person_days) 

        self._BM = dict.fromkeys(self.baseline_purchase_costs.keys(), 1)


    





       
