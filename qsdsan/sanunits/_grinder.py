
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 19:55:59 2021

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


import numpy as np
from qsdsan import SanUnit, Construction
from qsdsan.utils.loading import load_data, data_path
import os

__all__ = ('Grinder',)

#path to csv with all the inputs
#data_path = '/Users/lewisrowles/opt/anaconda3/lib/python3.8/site-packages/exposan/biogenic_refinery/_grinder.csv'
data_path += 'sanunit_data/_grinder.tsv'

### 
class Grinder(SanUnit):
    '''
    Grinder is used to break up solids.

    
    Reference documents
    -------------------
    N/A
    
    Parameters
    ----------
    ins : WasteStream
        Solids
    outs : WasteStream 
        Solids

        
    References
    ----------
    .. N/A
    
    '''
    

    def __init__(self, ID='', ins=None, outs=(), **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, F_BM_default=1)

# load data from csv each name will be self.name    
        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)




        
# define the number of influent and effluent streams    
    _N_ins = 1
    _N_outs = 1

# in _run: define influent and effluent streams and treatment processes 
    def _run(self):
        waste_in = self.ins[0]
        waste_out = self.outs[0]
        waste_out.copy_like(self.ins[0])

        moisture_content = (waste_in.imass['H2O'] / waste_in.F_mass) # fraction
        TS_in = waste_in.F_mass - waste_in.imass['H2O'] # kg TS dry/hr
        
        self.waste_sol_flow = waste_in.imass['OtherSS']

        # set necessary moisture content of effluent as 35%
        waste_out.imass['H2O'] = (0.65/0.35) * TS_in # fraction
        waste_out._COD = (waste_in.COD * (waste_in.F_vol/waste_out.F_vol))
        
    def _design(self):
        design = self.design_results
        
        # defining the quantities of materials/items
        # note that these items to be to be in the _impacts_items.xlsx
        # add function to calculate the number of screw presses required based on 
        # influent flowrate 
        design['Steel'] = S_quant = (self.grinder_steel)

        
        self.construction = (
            Construction(item='Steel', quantity = S_quant, quantity_unit = 'kg'),
            )
        self.add_construction(add_cost=False)
        

    #_cost based on amount of steel and stainless plus individual components
    def _cost(self):
        
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        self.baseline_purchase_costs['Grinder'] = (self.grinder)
        
        self.F_BM = dict.fromkeys(self.baseline_purchase_costs.keys(), 1)


       
        power_demand = (self.grinder_electricity * (self.waste_sol_flow)) #kW (all items are per hour)

        self.power_utility(power_demand)

        

