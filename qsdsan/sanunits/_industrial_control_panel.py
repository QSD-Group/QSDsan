#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
from qsdsan import SanUnit, Construction
from qsdsan.utils.loading import load_data, data_path
import os
__all__ = ('IndustrialControlPanel',)

#path to csv with all the inputs
#data_path = '/Users/lewisrowles/opt/anaconda3/lib/python3.8/site-packages/exposan/biogenic_refinery/_industrial_control_panel.csv'
#data_path = os.path.abspath(os.path.dirname('_industrial_control_panel.csv'))
data_path += 'sanunit_data/_industrial_control_panel.tsv'

### 
class IndustrialControlPanel(SanUnit):
    '''
    Industrial Controll Panel is composed of the electronic controller and 
    its components. 

    
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
                 **kwargs):
        if isinstance(ins, str) or (not isinstance(ins, Iterable)):
            self._N_outs = self._N_ins = 1
        else:
            self._N_outs = self._N_ins = len(ins)
        self._graphics = UnitGraphics.box(self._N_ins, self._N_outs)
        SanUnit.__init__(self, ID, ins, outs)


# load data from csv each name will be self.name    
        data = load_data(path=data_path)
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
        design['Electronics'] = Elect_quant = 2
        design['ElectricConnectors'] = ElectCon_quant = 0.5
        design['ElectricCables'] = ElectCables_quant = 3
        
        self.construction = (
            Construction(item='Electronics', quantity = Elect_quant, quantity_unit = 'kg'),
            Construction(item='ElectricConnectors', quantity = ElectCon_quant, quantity_unit = 'kg'),
            Construction(item='ElectricCables', quantity = ElectCables_quant, quantity_unit = 'm'),
            )
        self.add_construction()
    
    #_cost based on amount of steel and stainless plus individual components
    def _cost(self):
        
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        self.baseline_purchase_costs['Misc. parts'] = (self.icp_controller_board + 
                                              self.icp_variable_frequence_drives + 
                                              self.icp_power_meter + 
                                              self.icp_line_filter + 
                                              self.icp_transformer + 
                                              self.icp_power_meter_transformer + 
                                              self.icp_AC_to_DC + 
                                              self.icp_DC_to_AC + 
                                              self.icp_touch_screen)    

        
        self.F_BM = dict.fromkeys(self.baseline_purchase_costs.keys(), 1)


        
        #certain parts need to be replaced based on an expected lifefime
        #the cost of these parts is considered along with the cost of the labor to replace them
        cb_replacement_parts_annual_cost = (0) # USD/yr only accounts for time running
        
        cb_annual_maintenance = ((((self.electrician_replacecables_icp + self.electrician_replacewires_icp) / 60 * self.certified_electrician_wages)
                                  + (self.service_team_replacetouchscreen_icp / 60 * self.service_team_wages)
                                  + (self.facility_manager_configurevariable_icp / 60 * self.facility_manager_wages)
                                  + ((self.biomass_controls_replaceboard_icp + self.biomass_controls_codemalfunctioning_icp) / 60 * self.biomass_controls_wages)
                                  / self.frequency_corrective_maintenance)) #USD/yr only accounts for time running
        
        #need to be a cost per hour
        self.add_OPEX =  (cb_replacement_parts_annual_cost + cb_annual_maintenance) / (365 * 24) # USD/hr (all items are per hour)

       
        power_demand = (self.icp_controller_board_power + self.icp_variable_frequence_drives_power) #kWh/hr (all items are per hour)
        #needs to be power demand per hour
        self.power_utility(power_demand)
