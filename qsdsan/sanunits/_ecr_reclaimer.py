#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 14:31:47 2021

@author: torimorgan
"""
import numpy as np
from warnings import warn
from qsdsan import SanUnit, Construction
#from ._decay import Decay
from qsdsan.utils.loading import load_data, data_path

__all__ = ('ECR_Reclaimer',)

data_path += 'sanunit_data/_ECR_Reclaimer.csv'

X = 4 #number of reclaimers

class ECR_Reclaimer(SanUnit):
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
        SanUnit.__init__(self, ID, ins, outs, F_BM_default=1)
        
        self.price_ratio = 1
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
        waste = self.ins[0]
        treated = self.outs[0]
        treated.copy_like(self.ins[0])

              

    def _design(self):
        design = self.design_results
        design['Titanium'] = electrode_quant = self.Titanium_weight * X
        self.construction = ((Construction(item='Titanium', quantity = electrode_quant, quantity_unit = 'kg')))
        self.add_construction(add_cost=False)
 
    def _cost(self):
        
        self.EC_brush_scaled = (self.EC_brush) * X
        self.EC_cell_scaled = (self.EC_cell) * X
        C = self.baseline_purchase_costs

        C['EC_brush'] = self.EC_brush_scaled
        C['EC_cell'] = self.EC_cell_scaled
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        #self._BM = dict.fromkeys(self.baseline_purchase_costs.keys(), 1)
        
        self.power_utility(self.power_demand * X / 1000) #kW
        # self.power_utility(self.power_demand * 0)
    
    def _calc_replacement_cost(self):
        ecr_replacement_cost = ((self.EC_cell * (20/self.EC_cell_lifetime)) + 
        (self.EC_brush * (20/self.EC_brush_lifetime))) * X #USD/yr
        return ecr_replacement_cost/ (365 * 24) * self.price_ratio # USD/hr (all items are per hour)

        
       
        
