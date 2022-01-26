#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 17:33:27 2021

@author: torimorgan
"""

import numpy as np
from warnings import warn
from qsdsan import SanUnit, Construction
from ._decay import Decay
from ..utils import load_data, data_path

__all__ = ('MBRECR',)

data_path += 'sanunit_data/_MBR_ECR.tsv'


class MBRECR(SanUnit, Decay):
    '''
   Electrochemical treatment with chlorine dosing 

    Parameters
    ----------
    ins : WasteStream
        Waste for treatment, salt. 
    outs : Treated
        Treated waste, fugitive CH4, and fugitive N2O.
   
        
    References
    ----------
    .. 2019.06 Technical report for BMGF V3 _ CC 2019.06.13.pdf
    See Also
    --------
    :ref:`qsdsan.sanunits.Decay <sanunits_Decay>`
    
    '''
    
    def __init__(self, ID='', ins=None, outs=(), **kwargs):
        SanUnit.__init__(self, ID, ins, outs)
        
        self.price_ratio = 1
        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    
    _N_ins = 3
    _N_outs = 1

    
    def _run(self):
        waste, salt, HCl_acid = self.ins
        treated = self.outs[0]
        treated.copy_like(self.ins[0])
    
       
        
        salt.imass['NaCl'] = self.salt_dosing/7 #salt demand per day
        
        
        HCl_density = 1.2 #g/ml
        HCl_acid.imass['HCl'] = self.HCl_life/52/24/7 * HCl_density * 1000 #kg/h
        
        # COD removal
        COD_deg = treated._COD*treated.F_vol/1e3*self.COD_removal # kg/hr
        treated._COD *= (1-self.COD_removal) #mg/L
    


    def _design(self):
        #find rough value for FRP for tank 
        design = self.design_results
        design['Metal'] = electrode_quant = self.electrode
        design['Pump'] = pump_quant = self.pump
        self.construction = ((Construction(item='Pump', quantity = pump_quant, quantity_unit = 'ea')),
                             (Construction(item='Metal', quantity = electrode_quant, quantity_unit = 'kg')))
        self.add_construction()
 
    def _cost(self):
        C = self.baseline_purchase_costs
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        C['level_guage_cost'] = (self.level_guage_cost)
        C['pump_cost'] = (self.pump_cost)
        C['fan_cost'] = (self.fan_cost)
        C['salt_dosing_device_cost'] = (self.salt_dosing_device_cost)
        C['UPVC_pipe_cost'] = (self.UPVC_pipe_cost)
        C['UPVC_electric_ball_cost'] = (self.UPVC_electric_ball_cost)
        C['GAC_cost'] = (self.GAC_cost)
        C['electrode_cost'] = (self.electrode_cost)
        C['ECR_reactor_cost'] = (self.ECR_reactor_cost)
        C['power_supply_cost'] = (self.power_supply_cost)
        C['plastic_spraying_cabinent_cost'] = (self.plastic_spraying_cabinent_cost)
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        # self._BM = dict.fromkeys(self.purchase_costs.keys(), 1)
        
        #need to be a cost per hour
        
        ECR_replacement_parts_annual_cost = ((self.electrode_replacement_cost)
                                             +(self.HCl_replacement_cost * self.HCl_life) 
                                               + (self.salt_replacement_cost * 52 )) * self.price_ratio
     
        self.add_OPEX =  (ECR_replacement_parts_annual_cost) / (365 * 24) # USD/hr (all items are per hour)

       
        self.power_utility(self.power_demand * self.working_time) #kWh/hr (all items are per hour)

