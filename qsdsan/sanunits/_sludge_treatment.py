# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 11:15:25 2022

@author: raisa
"""
from .. import SanUnit, WasteStream
import numpy as np

class Thickener(SanUnit):
    
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, 
                  init_with='WasteStream', F_BM_default=None, thickner_perc=7, TSS_removal_perc=98, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic, init_with=init_with, 
                         F_BM_default=F_BM_default)
        
        @property
        def thickner_perc(self):
            '''tp is the percentage of Suspended Sludge in the underflow of the thickener'''
            return self._tp

        @thickner_perc.setter
        def thickner_perc(self, tp):
            if tp is not None:
                if tp>=100 or tp<=0:
                    raise ValueError(f'should be between 0 and 100 not {tp}')
                self._tp = tp
            else: 
                raise ValueError('percentage of SS in the underflow of the thickener expected from user')
                
        @property
        def TSS_removal_perc(self):
            '''The percentage of suspended solids removed in the thickner'''
            return self._TSS_rmv

        @TSS_removal_perc.setter
        def TSS_removal_perc(self, TSS_rmv):
            if TSS_rmv is not None:
                if TSS_rmv>=100 or TSS_rmv<=0:
                    raise ValueError(f'should be between 0 and 100 not {TSS_rmv}')
                self._TSS_rmv = TSS_rmv
            else: 
                raise ValueError('percentage of suspended solids removed in the thickner expected from user')
                
        @property
        def thickner_factor(self):
            inf, = self.ins[0]
            if not self.ins: return
            elif inf.isempty(): return
            else: 
                TSS_in = inf.get_TSS()
                if TSS_in > 0:
                    thickner_factor = self._tp*10000/self.ins[0].get_TSS()
                    if thickner_factor<1:
                        thickner_factor=1
                    return thickner_factor
                else: return None
        
        @property
        def thinning_factor(self):
            thickner_factor = self.thickner_factor
            if thickner_factor<1:
                thinning_factor=0
            else:
                Qu_factor = self._TSS_rmv/(100*thickner_factor)
                thinning_factor = (1 - (self._TSS_rmv/100))/(1 - Qu_factor)
            return thinning_factor
        
        def _run(self):
            
            # self.inlet = WasteStream('inlet')
            # self.inlet = self.ins
            inf, = self.ins
            uf, of = self.outs
            
            split = 
            
            inf.split_to(uf, of, split)
            
            
            cmps = self.components
            
          
                
                
            
            
            
            
            
        
                
        
                
        
                
        
        
        
        
        
