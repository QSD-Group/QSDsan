#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 13:41:55 2021

@author: torimorgan
"""

import numpy as np
from warnings import warn
from qsdsan import SanUnit, Construction
#from ._decay import Decay
from qsdsan.utils.loading import load_data, data_path


__all__ = ('SystemReclaimer',)

data_path += 'sanunit_data/_system_reclaimer.csv'

# R = 1 #number of reclaimers



class SystemReclaimer(SanUnit):
    '''
    Cost and life cycle impacts of the system components for the Reclaimer 2.0
    
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', 
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, F_BM_default=1)
        self.price_ratio = 1

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

    def _design(self):
        #find rough value for FRP for tank 
        design = self.design_results
        #!!! Consider doing later design['Brass'] = brass_quant = self.aluminum_weight
        design['Steel'] = steel_quant = self.steel_weight * 1
        self.construction = (
                            (Construction(item='Steel', quantity = steel_quant, quantity_unit = 'kg')))
        self.add_construction(add_cost=False)
        
 
    def _cost(self):
        C = self.baseline_purchase_costs

        C['System'] = ((self.T_nut + self.die_cast_hinge + self.SLS_locks + self.DC_round_key
                                         + self.handle_rod + self.eight_mm_bolt + self.button_headed_nut
                                         + self.twelve_mm_bolt + self.ten_mm_CSK + self.sixteen_mm_bolt 
                                         + self.coupling_brass + self.socket + self.onehalf_tank_nipple + self.onehalf_in_coupling_brass
                                         + self.onehalf_in_fitting + self.plate + self.pump + self.three_way_valve + self.lofted_tank)) 
                                         # + ((self.T_nut + self.die_cast_hinge + self.SLS_locks + self.DC_round_key
                                         # + self.handle_rod + self.eight_mm_bolt + self.button_headed_nut
                                         # + self.twelve_mm_bolt + self.ten_mm_CSK + self.sixteen_mm_bolt 
                                         # + self.coupling_brass + self.socket + self.onehalf_tank_nipple + self.onehalf_in_coupling_brass
                                         # + self.onehalf_in_fitting + self.plate + self.pump + self.three_way_valve + self.lofted_tank) * R * .05))
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        #self._BM = dict.fromkeys(self.baseline_purchase_costs.keys(), 1)
        
        # #certain parts need to be replaced based on an expected lifefime
        # #the cost of these parts is considered along with the cost of the labor to replace them
        # replacement_parts_annual_cost = (self.T_nut + self.die_cast_hinge + self.SLS_locks + self.DC_round_key
        #                                  + self.handle_rod + self.eight_mm_bolt + self.button_headed_nut
        #                                  + self.twelve_mm_bolt + self.ten_mm_CSK + self.sixteen_mm_bolt 
        #                                  + self.coupling_brass + self.socket + self.onehalf_tank_nipple + self.onehalf_in_coupling_brass
        #                                  + self.onehalf_in_fitting + self.plate + self.pump + self.three_way_valve + self.lofted_tank) * .1
     
        #self.power_utility(self.power_demand / 1000) #kW
        self.power_utility(self.power_demand * 0)
        
    def _calc_replacement_cost(self):
        controls_replacement_cost = (self.replacement_costs) / 20 * self.price_ratio #USD/yr
        return controls_replacement_cost/ (365 * 24) # USD/hr (all items are per hour)
        


