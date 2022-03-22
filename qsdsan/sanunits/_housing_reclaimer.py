#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 10:56:00 2021

@author: torimorgan
"""

from qsdsan import SanUnit, Construction
#from ._decay import Decay
from qsdsan.utils.loading import load_data, data_path


__all__ = ('HousingReclaimer',)

data_path += 'sanunit_data/_housing_reclaimer.csv'

#To scale the system from 0 - 120 users based on units corresponding to 30 users: 
ppl = 120

if (ppl <=30): R = 1 
elif (ppl >30 and ppl <=60): R = 2
elif (ppl >60 and ppl <=90): R = 3
elif (ppl >90 and ppl <=120): R = 4
else: R = 5

#To scale the system from 0-120 users based on units corresponding to 25 users:
if (ppl <=25): D = 1 
elif (ppl >25 and ppl <=50): R = 2
elif (ppl >50 and ppl <=75): R = 3
elif (ppl >75 and ppl <=100): R = 4
else: R = 5

class HousingReclaimer(SanUnit):
    '''
    Cost and life cycle impacts of the housing for the Reclaimer 2.0
    
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', 
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, F_BM_default=1)
        self.price_ratio = 1
  
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
        design['Steel'] = steel_quant = (self.steel_weight + (self.framework_weight/4) + self.fittings_weight)
        self.construction = ((Construction(item='Steel', quantity = steel_quant, quantity_unit = 'kg')))
        self.add_construction(add_cost=False)
        
 
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Housing'] = ((self.frame + self.extrusion + self.angle_frame + self.angle + self.door_sheet + self.plate_valve + self.powder)
        + ((self.frame + self.extrusion + self.angle_frame + self.angle + self.door_sheet + self.plate_valve + self.powder) * .1 * R)
        + (self.portable_toilet * (D)))
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        

