#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 23 14:26:28 2021

@author: torimorgan
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 16 15:31:50 2021

@author: torimorgan
"""

import numpy as np
from warnings import warn
from qsdsan import SanUnit, Construction
from ._decay import Decay
from ..utils import load_data, data_path

__all__ = ('PrimaryMBR',)

data_path += 'sanunit_data/_primary_MBR.tsv'





class PrimaryMBR(SanUnit, Decay):
    '''
    Anaerobic digestion of wastes in primary treatment of Ecosan

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
    
    
    _N_ins = 2
    _N_outs = 5

#Anaerobic is 60-80%     
    def _run(self):
        waste, MgOH2 = self.ins
        treated, CH4, N2O, sludge, struvite = self.outs
        treated.copy_like(self.ins[0])
        CH4.phase = N2O.phase = 'g'
        MgOH2.phase = struvite.phase = 's'
        
        MgOH2_mass = (self.orthoP_post_AnMBR / self.MW_P 
                    * self.Mg_dose * self.MW_MgOH2) * waste.F_vol /1e3  # kg Mg(OH)2 per hr;
        MgOH2.imass['MagnesiumHydroxide'] = MgOH2_mass
        
        #MgOH2_cost_day = MgOH2.imass['MagnesiumHydroxide'] * self.MgOH2_cost * 24
        
        P_precipitated = self.orthoP_post_AnMBR * self.P_recovery * waste.F_vol / 1e3  # kg precipitated-P/hr
        N_precipitated = (self.orthoP_post_AnMBR * self.P_recovery / 
                          self.MW_P * self.N_P_ratio_struvite * self.MW_N) *  waste.F_vol /1e3 # kg precipitated-N/hr
        
        
        struvite_production_time = P_precipitated / self.MW_P * self.MW_struvite  # kg (NH4)MgPO4•6(H2O) / hr
        struvite.imass['Struvite'] = struvite_production_time
        
        # COD removal
        COD_deg = treated.COD*treated.F_vol/1e3*self.COD_removal # kg/hr
        #add conversion of COD to C here 
        treated._COD = waste.COD * (1-self.COD_removal)
        
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
        
        N_loss_tot = N_loss*waste.TN/1e3*waste.F_vol #work through the 
            #conversions to make sure they are applicable
        NH3_rmd, NonNH3_rmd = \
            self.allocate_N_removal(N_loss_tot, waste.imass['NH3'])
        treated.imass ['NH3'] = waste.imass['NH3'] - NH3_rmd
        treated.imass['NonNH3'] = waste.imass['NonNH3'] - NonNH3_rmd
        N2O.imass['N2O'] = N_loss_tot*self.N2O_EF_decay*44/28 #check units
        
        #sludge production
        self.sludge_TN = waste.imass['N'] * self.N_max_decay
        self.sludge_TN_F_mass = self.sludge_TN * waste.F_vol * 1e-3
        sludge.imass['N'] = self.sludge_TN
        
        self.sludge_TP = waste.imass['P'] * self.TP_removal
        self.sludge_TP_F_mass = self.sludge_TP * waste.F_vol * 1e-3 
        sludge.imass['P'] = self.sludge_TP 
        
        self.sludge_COD = waste._COD * self.COD_removal 
        self.sludge_COD_F_mass = self.sludge_COD * waste.F_vol * 1e-3
        sludge._COD = self.sludge_COD 
        
        self.sludge_prcd = (self.reactor_volume * waste.imass['OtherSS']
                            * self.VSS_TSS_ratio)
        sludge.imass['OtherSS'] = (self.sludge_prcd - self.sludge_TN_F_mass 
         - self.sludge_TP_F_mass - self.sludge_COD_F_mass)
        
        sludge.imass['H2O'] = sludge.F_mass * 0.5        
#no replacement parts for the anaerobic tank and cleaning performed 
#throughout whole system is considered in TEA 


    def _design(self):
        #find rough value for FRP for tank 
        design = self.design_results
        design['FRP'] = FRP_quant = self.FRP_per_tank  
        self.construction = (Construction(item='FRP', quantity = FRP_quant, quantity_unit = 'kg'))
        self.add_construction()
 
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tanks'] = (self.FRP_tank_cost)
        # self._BM = dict.fromkeys(self.purchase_costs.keys(), 1)
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
         #certain parts need to be replaced based on an expected lifefime
        #the cost of these parts is considered along with the cost of the labor to replace the
        PrimaryMBR_replacement_costs = (self.Mg_dose * self.MgOH2_cost / 365 / 24) * self.price_ratio
       
        self.add_OPEX =  (PrimaryMBR_replacement_costs) 
        
        self.power_utility(self.power_demand)
