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
        design['Titanium'] = electrode_quant = self.Titanium_weight * 4
        self.construction = ((Construction(item='Titanium', quantity = electrode_quant, quantity_unit = 'kg')))
        self.add_construction(add_cost=False)
 
    def _cost(self):
        
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        self.baseline_purchase_costs['EC_brush'] = (self.EC_brush * 4)
        self.baseline_purchase_costs['EC_cell'] = (self.EC_cell * 4)
        
        self._BM = dict.fromkeys(self.baseline_purchase_costs.keys(), 1)
        
        self.power_utility(self.power_demand * 4 / 1000) #kW
        #self.power_utility(self.power_demand * self.working_time)
    
    def _calc_replacement_cost(self):
        ecr_replacement_cost = ((self.EC_cell * (20/self.EC_cell_lifetime)) + 
        (self.EC_brush * (20/self.EC_brush_lifetime))) * 4 #USD/yr
        return ecr_replacement_cost/ (365 * 24) # USD/hr (all items are per hour)

        
       
        
