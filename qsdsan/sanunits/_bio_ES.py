#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 12:03:02 2021

@author: torimorgan
"""
import numpy as np
from warnings import warn
from qsdsan import SanUnit, Construction
from ._decay import Decay
from ..utils import load_data, data_path

__all__ = ('BioES',)


data_path += 'sanunit_data/_bio_ES_1.tsv'




class BioES(SanUnit):
    '''
    All of the costs and impacts associated with biological pretreatment for the Eco-san system

    Parameters
    ----------
    ins : WasteStream
        Waste for treatment.
    outs : WasteStream
        Treated waste, fugitive CH4, and fugitive N2O.
   
        
    References
    ----------
    .. 2019.06 Technical report for BMGF V3 _ CC 2019.06.13.pdf
    See Also
    --------
    :ref:`qsdsan.sanunits.Decay <sanunits_Decay>`
    
    '''
    
    #no replacement parts for the anaerobic tank and cleaning performed 
#throughout whole system is considered in TEA 

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', 
                 **kwargs):
        # if isinstance(ins, str) or (not isinstance(ins, Iterable)):
        #     self._N_outs = self._N_ins = 1
        # else:
        #     self._N_outs = self._N_ins = len(ins)
        # self._graphics = UnitGraphics.box(self._N_ins, self._N_outs)
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
        
    def _design(self):
        #find rough value for FRP for tank 
        design = self.design_results
        design['FRP'] = FRP_quant = self.FRP
        design['Pump'] = pump_quant = self.pump_ammount
        self.construction = (
                             Construction(item='Pump', quantity = pump_quant, quantity_unit = 'ea'),
                             Construction(item='FRP', quantity = FRP_quant, quantity_unit = 'kg'))
        self.add_construction()
 
    def _cost(self):
        
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        self.purchase_costs['Tanks'] = (self.tank_cost)
        self.purchase_costs['Misc. parts'] = (self.middle_tank_pump_cost + 
                                              self.fan_cost + 
                                              self.level_guage_cost + 
                                              self.pipe_cost +
                                              self.regulating_tank_pump_cost)
        
        self._BM = dict.fromkeys(self.purchase_costs.keys(), 1)
        
        #certain parts need to be replaced based on an expected lifefime
        #the cost of these parts is considered along with the cost of the labor to replace them
        bio_replacement_parts_annual_cost = (self.guage_life * self.level_guage_cost) 
        + (self.pump_life * self.pump_replacement_cost)
        + (self.bio_tube_diffuser_replacement_cost * self.diffuser_life)
        # USD/yr only accounts for time running
                                                  
        self.add_OPEX =  (bio_replacement_parts_annual_cost) / (365 * 24) # USD/hr (all items are per hour)
        
        #self.power_utility(self.power_demand * self.working_time)
        
