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

__all__ = ('HousingBiogenicRefinery',)

#path to csv with all the inputs

import os
su_data_path = os.path.join(data_path, 'sanunit_data/_housing_biogenic_refinery.tsv')

### 
class HousingBiogenicRefinery(SanUnit):
    '''
    Housing Biogenic Refinery is composed of the casing around the system, 
    containers, and the concrete slab. 

    
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
    

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 const_wage=15, const_person_days=100,
                 **kwargs):
        if isinstance(ins, str) or (not isinstance(ins, Iterable)):
            self._N_outs = self._N_ins = 1
        else:
            self._N_outs = self._N_ins = len(ins)
        self._graphics = UnitGraphics.box(self._N_ins, self._N_outs)
        SanUnit.__init__(self, ID, ins, outs, F_BM_default=1)
        self.price_ratio = 1
        self.const_wage = const_wage
        self.const_person_days = const_person_days

# load data from csv each name will be self.name    
        data = load_data(path=su_data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    def _run(self):
        waste = self.ins[0]
        treated = self.outs[0]
        treated.copy_like(self.ins[0])

        
    #_design will include all the construction or captial impacts  
    def _design(self):
        design = self.design_results
        
        # defining the quantities of materials/items
        # note that these items to be to be in the _impacts_items.xlsx
        design['Steel'] = S_quant = (2000 + 4000)
        design['StainlessSteelSheet'] = SSS_quant = (4.88 * 2 * 3 * 4.5)
        design['Concrete'] = Con_quant = self.concrete_thickness * self.footprint
        
        self.construction = (
            Construction(item='Steel', quantity = S_quant, quantity_unit = 'kg'),
            Construction(item='StainlessSteelSheet', quantity = SSS_quant, quantity_unit = 'kg'),
            Construction(item='Concrete', quantity = Con_quant, quantity_unit = 'm3'),
            )
        self.add_construction()
        
 
    #_cost based on amount of steel and stainless plus individual components
    def _cost(self):
        
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        C = self.baseline_purchase_costs
        C['Containers'] = (self.container20ft_cost + self.container40ft_cost)
        C['Equip Housing'] = (self.design_results['StainlessSteelSheet'] / 4.88 * self.stainless_steel_housing)  
        C['Concrete'] = (self.design_results['Concrete'] * self.concrete_cost) 
        C['Labor'] = (self.const_wage*self.const_person_days) 
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio 


    





       
