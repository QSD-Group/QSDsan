#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tues Jan 4 13:37:22 2022

@author: torimorgan
"""


from qsdsan import SanUnit, Construction
from qsdsan.utils.loading  import load_data, data_path

__all__ = ('SolarReclaimer',)

data_path += 'sanunit_data/_solar_reclaimer.csv'


class SolarReclaimer(SanUnit):
    '''
    All of the costs related to the solar energy system

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
        design['Solar'] = solar_quant = self.solar_capacity
        design['Battery'] = battery_quant = self.battery_kg
        self.construction = (
                             Construction(item='Solar', quantity = solar_quant, quantity_unit = 'm2'))
        self.construction = (
                             Construction(item='Battery', quantity = battery_quant, quantity_unit = 'kg'))
        self.add_construction(add_cost=False)

    def _calc_labor_cost(self):
        labor_cost = self.wages * self.pannel_cleaning
        return labor_cost/(365 * 24)  # USD/hr (all items are per hour)
    
    def _calc_replacement_cost(self):
        solar_replacement_parts_annual_cost = (self.solar_replacement * self.solar_cost)  
        return solar_replacement_parts_annual_cost/ (365 * 24) * self.price_ratio # USD/hr (all items are per hour)
                  
    def _calc_maintenance_labor_cost(self):
        solar_maintenance_labor = (self.pannel_cleaning * self.wages)
        return solar_maintenance_labor / (365 * 24) # USD/hr (all items are per hour)
    
    F_BM = {'Battery System': 1, 'Solar Cost':1}
    def _cost(self):
        #purchase_costs is used for capital costs
        C=self.baseline_purchase_costs
        C['Battery System'] = (self.battery_storage_cost + self.battery_holder_cost)
        C['Solar Cost'] = ((self.solar_cost * self.power_demand_30users) + self.solar_module_system 
                                                      + self.inverter_cost)
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        add_OPEX = self._calc_replacement_cost() 
        self._add_OPEX = {'Additional OPEX': add_OPEX}
        
  
        
        
        
        
        
        
        
        
        
        
        
