#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 18:43:10 2021

@author: lewisrowles
"""

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

__all__ = ('HydronicHeatExchanger',)

#data_path = '/Users/lewisrowles/opt/anaconda3/lib/python3.8/site-packages/exposan/biogenic_refinery/_hydronic_heat_exchanger.csv'

#data_path = os.path.abspath(os.path.dirname('_hydronic_heat_exchanger.csv'))
data_path += '/sanunit_data/_hydronic_heat_exchanger.tsv'

### 
class HydronicHeatExchanger(SanUnit):
    '''
    Hydronic heat exchanger is used for applications that require drying 
    of the feedstock before the refinery is capable of processing the 
    material. The heat is exchanged between the exhaust gas and water, 
    which is then pumped into radiators connected to a dryer. 
    The refinery monitors the temperature of the water to ensure that 
    the feedstock is being sufficiently dried before entering the refinery.


    
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
    

    def __init__(self, ID='', ins=None, outs=(),init_with='WasteStream', **kwargs):
        
        
        SanUnit.__init__(self, ID, ins, outs, init_with=init_with)

    
        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    

        
    
    _N_ins = 1
    _N_outs = 1
    
    def _run(self):
        hot_gas_in = self.ins[0]
        hot_gas_out = self.outs[0]
        hot_gas_out.phase = 'g'
        
        # set temperature
        hot_gas_out.T = self.hhx_temp
        # calculate the heat that was delivered to the HHX
        self.heat_output_water = self.water_flowrate *  ((273.15 + self.water_out_temp) - (273.15 + self.water_in_temp)) * self.water_density_kg_m_3 * self.water_heat_capacity_k_j_kg_k * (60/1000)  # MJ/hr
        # calculate losses through water pipe
        self.heat_loss_water_pipe = ((273.15 + self.water_out_temp) - (273.15 + self.ambient_temp)) / (self.water_r_pipe_k_k_w + self.water_r_insulation_k_k_w) * 3.6 # MJ/hr

        
    #_design will include all the construction or captial impacts  
    def _design(self):
        design = self.design_results
        
        # defining the quantities of materials/items
        # note that these items to be to be in the _impacts_items.xlsx
        design['StainlessSteel'] = SS_quant = (self.heat_exchanger_hydronic_stainless)
        design['Steel'] = S_quant = (self.heat_exchanger_hydronic_steel)
        design['HydronicHeatExchanger'] =  HHX_quant = 1
        design['Pump'] = Pump_quant = (17.2/2.72)
        
        
        
        self.construction = (
            Construction(item='StainlessSteel', quantity = SS_quant, quantity_unit = 'kg'),
            Construction(item='Steel', quantity = S_quant, quantity_unit = 'kg'),
            Construction(item='HydronicHeatExchanger', quantity = HHX_quant, quantity_unit = 'ea'),
            Construction(item='Pump', quantity = Pump_quant, quantity_unit = 'ea'),
            )
        self.add_construction()
        
    
    #_cost based on amount of steel and stainless plus individual components
    def _cost(self):
        
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        self.baseline_purchase_costs['Stainless steel'] = (self.stainless_steel_cost 
                            * self.design_results['StainlessSteel'])
        self.baseline_purchase_costs['Steel'] = (self.steel_cost
                            * self.design_results['Steel'])
        self.baseline_purchase_costs['Misc. parts'] = (self.hhx_stack 
                                              + self.hhx_stack_thermocouple 
                                              + self.hhx_oxygen_sensor 
                                              + self.hhx_inducer_fan 
                                              + self.hhx_flow_meter 
                                              + self.hhx_pump 
                                              + self.hhx_water_in_thermistor 
                                              + self.hhx_water_out_thermistor 
                                              + self.hhx_load_tank 
                                              + self.hhx_expansion_tank 
                                              + self.hhx_heat_exchanger 
                                              + self.hhx_values 
                                              + self.hhx_thermal_well 
                                              + self.hhx_hot_water_tank 
                                              + self.hhx_overflow_tank)
       
        self.F_BM = dict.fromkeys(self.baseline_purchase_costs.keys(), 1)
        
        
        #certain parts need to be replaced based on an expected lifefime
        #the cost of these parts is considered along with the cost of the labor to replace them
        cb_replacement_parts_annual_cost = (0) # USD/yr only accounts for time running
        
        cb_annual_maintenance = ((self.service_team_adjustdoor_hhx / 60 * 12 * self.service_team_wages) 
                                 + (self.service_team_replacewaterpump_hhx / 60 * (1 / self.frequency_corrective_maintenance) * self.service_team_wages)
                                 + (self.service_team_purgewaterloop_hhx / 60 * (1 / self.frequency_corrective_maintenance) * self.service_team_wages)) #USD/yr only accounts for time running
        
        #need to be a cost per hour
        self.add_OPEX =  (cb_replacement_parts_annual_cost + cb_annual_maintenance) / (365 * 24) # USD/hr (all items are per hour)
        
        power_demand = (self.water_pump_power + self.hhx_inducer_fan_power) #kW
        self.power_utility(power_demand) #kW
        #breakpoint()
        
        
        
        
        
