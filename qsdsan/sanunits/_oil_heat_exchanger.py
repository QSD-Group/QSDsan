#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Mar  9 13:41:03 2021

@author: lewisrowles stetson@gmail.com
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
__all__ = ('OilHeatExchanger',)

su_data_path = os.path.join(data_path,'sanunit_data/_oil_heat_exchanger.tsv')

### 
class OilHeatExchanger(SanUnit):
    '''
    Oil heat exhanger utilizes an organic Rankin cycle. This CHP system is 
    used to generate additional electricity that the refinery and/or 
    facility can use to decrease the units electrical demand on the 
    electrical grid. This type of system is required for ISO 31800 
    certification as the treatment unity needs to be energy independent 
    when processing faecal sludge. 
    
    Reference documents
    -------------------
    N/A
    
    Parameters
    ----------
    ins : WasteStream
        Hot gas.
    outs : WasteStream 
        Hot gas.
            
            
    References
    ----------
    .. N/A
    
    '''
    

    def __init__(self, ID='', ins=None, outs=(), **kwargs):
        
        
        SanUnit.__init__(self, ID, ins, outs)

        self.price_ratio = 1 
    
        data = load_data(path=su_data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    

        
    
    _N_ins = 1
    _N_outs = 1
    
    def _run(self):
        heat_in = self.ins[0]
        heat_out = self.outs[0]
        heat_out.copy_like(heat_in)
        heat_out.phase = 'g'
        
        # set temperature
        heat_out.T = self.ohx_temp

        # calculate the power that was delivered to the ORC
        self.power_delivery_orc = self.oil_flowrate * ((273.15 + self.oil_temp_out) - (273.15 + self.oil_temp_in)) * self.oil_density * self.oil_specific_heat * (60/1000)  # MJ/hr
        # calculate losses through pipe
        self.pipe_heat_loss = ((273.15 + self.oil_temp_out) - (273.15 + self.amb_temp)) / (self.oil_r_pipe_k_k_w + self.oil_r_insulation_k_k_w) * 3.6  # MJ/hr

    
    #_design will include all the construction or captial impacts  
    def _design(self):
        design = self.design_results
        
        design['OilHeatExchanger'] = OHX_quant = (4/200)
        design['Pump'] = Pump_quant = (2.834/2.27)
        
        
        self.construction = (
            Construction(item='OilHeatExchanger', quantity = OHX_quant, quantity_unit = 'ea'),
            Construction(item='Pump', quantity = Pump_quant, quantity_unit = 'ea'),
            )
        self.add_construction()

    
    
    #_cost based on amount of steel and stainless plus individual components
    def _cost(self):
    
        self.baseline_purchase_costs['Oil Heat Exchanger'] = (self.orc_cost) * self.price_ratio
        
        
        power_demand = (self.oil_pump_power - self.oil_electrical_energy_generated) #kWh
        self.power_utility(power_demand) #kWh
        #breakpoint()
        
        



       
