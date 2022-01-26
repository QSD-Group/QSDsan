#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 11:06:37 2021

@author: torimorgan
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 14:31:47 2021

@author: torimorgan
"""
import numpy as np
from warnings import warn
from qsdsan import SanUnit, Construction
from ._decay import Decay
from ..utils import load_data, data_path

__all__ = ('MBR',)

data_path += 'sanunit_data/_MBR.tsv'


class MBR(SanUnit, Decay):
    '''
    Aerobic digestion of waste
    Parameters with a membrane bioreactor 
    ----------
    ins : WasteStream
        Waste for treatment.
    outs : WasteStream
        Treated waste, and fugitive N2O.
   
        
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
    
    
    _N_ins = 1
    _N_outs = 3

#look up literature values for percentages typically removed by anaerobic/follow yalin/john's assumptions \
#Anaerobic is 60-80%     
    def _run(self):
        waste = self.ins[0]
        treated, N2O, CH4 = self.outs
        treated.copy_like(self.ins[0])
        CH4.phase = N2O.phase = 'g'
        
        
        # COD removal
        COD_deg = treated._COD*treated.F_vol/1e3*self.COD_removal # kg/hr
        treated._COD *= (1-self.COD_removal) #mg/L
        
        CH4_prcd = COD_deg*self.MCF_decay*self.max_CH4_emission
        CH4.imass['CH4'] = CH4_prcd 
            
        N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.tau/365,
                                            max_decay=self.N_max_decay)
        N_loss_tot = N_loss*waste.TN/1e3*waste.F_vol #work through the 
            #conversions to make sure they are applicable
        NH3_rmd, NonNH3_rmd = \
            self.allocate_N_removal(N_loss_tot, waste.imass['NH3'])
        treated.imass ['NH3'] = waste.imass['NH3'] - NH3_rmd
        treated.imass['NonNH3'] = waste.imass['NonNH3'] - NonNH3_rmd
        N2O.imass['N2O'] = N_loss_tot*self.N2O_EF_decay*44/28 #check units
                
#no replacement parts for the anaerobic tank and cleaning performed 
#throughout whole system is considered in TEA 


    def _design(self):
        #find rough value for FRP for tank 
        design = self.design_results
        design['FRP'] = FRP_quant = self.FRP_per_tank     
        self.construction = (Construction(item='FRP', quantity = FRP_quant, quantity_unit = 'kg'),)
        self.add_construction()
 
    def _cost(self):
        C = self.baseline_purchase_costs

    
        #purchase_costs is used for capital costs
        #can use quantities from above (e.g., self.design_results['StainlessSteel'])
        #can be broken down as specific items within purchase_costs or grouped (e.g., 'Misc. parts')
        C['MBRTank'] = (self.MBR_cost)
        C['BioTank'] = (self.tank_cost)
        C['PrimaryTank'] = (self.primary_tank_cost)
        C['RegulatingTank'] = (self.regulating_pump_cost)
        C['fancost'] = (self.fan_cost)
        C['guagecost'] = (self.guage_cost)
        C['pipecost'] = (self.pipe_cost)
        C['midtankpump'] = (self.mid_tank_pump)
        # self._BM = dict.fromkeys(self.purchase_costs.keys(), 1)
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
         #need to be a cost per hour
        MBR_replacement_parts_annual_cost = ((self.MBR_cost * self.MBR_replacement)
                    + (self.guage_cost * self.guage_life)) 
        
        self.add_OPEX =  (MBR_replacement_parts_annual_cost) / (365 * 24) * self.price_ratio # USD/hr (all items are per hour)
        #self.power_utility(self.power_demand * self.working_time)
        
       #power is negligible in the aerobic Tank 
    