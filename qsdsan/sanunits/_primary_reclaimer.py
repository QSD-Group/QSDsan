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

# __all__ = ('PrimaryES',)
__all__ = ('PrimaryReclaimer',)

data_path += 'sanunit_data/_primary_reclaimer.csv'

P = 4/4

class PrimaryReclaimer(SanUnit, Decay):
    '''
    Anaerobic digestion of wastes in primary treatment of Reclaimer

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
        SanUnit.__init__(self, ID, ins, outs, F_BM_default=1)
        
        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    
    _N_ins = 2
    _N_outs = 5

     
    def _run(self):
        waste, MgOH2 = self.ins
        self.ins[1].imass['MagnesiumHydroxide'] = (self.orthoP_post / self.MW_P / self.Mg_dose * self.MW_MgOH2)/24  # kg Mg(OH)2 per hr;

        treated, CH4, N2O, sludge, struvite = self.outs
        treated.copy_like(self.ins[0])
        CH4.phase = N2O.phase = 'g'
        
        
        P_precipitated = self.orthoP_post * self.P_recovery  # kg precipitated-P/hr
        N_precipitated = (self.orthoP_post * self.P_recovery / self.MW_P * self.N_P_ratio_struvite * self.MW_N)  # kg precipitated-N/hr
        
        
        struvite_production_time = P_precipitated / self.MW_P * self.MW_struvite  # kg (NH4)MgPO4•6(H2O) / hr
        struvite.imass['Struvite'] = struvite_production_time
        
        # COD removal
        COD_deg = treated.COD*treated.F_vol/1e3*self.COD_removal # kg/hr
        #add conversion of COD to C here 
        treated._COD = waste.COD * (1-self.COD_removal)
       
        
        #do this for all sanunits
        
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
        
        #sludge production
        self.sludge_TN = waste.imass['N'] * self.N_max_decay
        self.sludge_TN_F_mass = self.sludge_TN * waste.F_vol * 1e-3
        sludge.imass['N'] = self.sludge_TN
        
        
        self.sludge_COD = waste._COD * self.COD_removal 
        self.sludge_COD_F_mass = self.sludge_COD * waste.F_vol * 1e-3
        sludge._COD = self.sludge_COD 
        

        sludge.imass['H2O'] = sludge.F_mass * 0.5

    def _design(self):
        #find rough value for FRP for tank 
        design = self.design_results
        design['FRP'] = FRP_quant = self.FRP_per_tank 
        design['Pump'] = pump_quant = self.pump_lca * 2
        self.construction = (Construction(item='FRP', quantity = FRP_quant, quantity_unit = 'kg'),
                             Construction(item='Pump', quantity = pump_quant, quantity_unit = 'each'))
        self.add_construction(add_cost = False)
 
    def _cost(self):
        
        self.baseline_purchase_costs['Tanks'] = ((self.FRP_tank_cost * P) + self.pump)
        self._BM = dict.fromkeys(self.baseline_purchase_costs.keys(), 1)
    
    def _calc_replacement_cost(self):
        Mg_cost = (self.orthoP_post / self.MW_P / self.Mg_dose * self.MW_MgOH2 * self.MgOH2_cost) * P
        return Mg_cost # USD/hr (all items are per hour)
        

