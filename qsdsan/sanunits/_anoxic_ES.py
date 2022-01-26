#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  6 09:16:06 2021

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

__all__ = ('AnoxicES',)


data_path += 'sanunit_data/_anoxic_ES.tsv'


class AnoxicES(SanUnit, Decay):
    '''
    Anoxic digestion of wastes.
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
        treated, CH4, N2O = self.outs
        treated.copy_like(self.ins[0])
        CH4.phase = N2O.phase = 'g'
        
        # COD removal
        COD_deg = treated.COD*treated.F_vol/1e3*self.COD_removal # kg/hr
        treated._COD = waste.COD * (1-self.COD_removal)
        
        #fix this and look at IPCC for anoxic 
        # CH4_prcd = COD_deg*self.MCF_decay*self.max_CH4_emission
        # CH4.imass['CH4'] = CH4_prcd  
    
        CH4_prcd = COD_deg*self.MCF_decay*self.max_CH4_emission
        CH4.imass['CH4'] = CH4_prcd       
        N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.tau/365,
                                            max_decay=self.N_max_decay)    

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
#the costs of the tanks and construction impacts are accounted for in bio files

