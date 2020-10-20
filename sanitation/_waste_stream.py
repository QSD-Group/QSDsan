#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 07:48:41 2020

@author: yalinli_cabbi
"""


# import thermosteam as tmo
from thermosteam import Stream


__all__ = ('WasteStream',)


# %%

# =============================================================================
# Define the WasteStream class
# =============================================================================

class WasteStream(Stream):
    '''A subclass of the Stream object in the thermosteam package with additional attributes and methods for waste treatment    '''
    
    def show(self, T=None, P=None, flow=None, composition=None, N=None,
             stream_info=True):
        '''Show WasteStream information'''
        
        # Stream-related specifications
        if stream_info:
            super().show(T, P, flow, composition, N)
        else:
            print(self._basic_info())
            display_units = self.display_units
            T_units = T or display_units.T
            P_units = P or display_units.P
            print(self._info_phaseTP(self.phase, T_units, P_units))
        
        # Component-related properties
        info = '\n Component-specific properties:\n'
        info += f'  charge: {self.charge} mol/hr\n'
        
        print(info)
        
    _ipython_display_ = show
    
    @property
    def charge(self):
        '''[float] Total charge of the stream in + mol/hr'''
        return (self.chemicals.i_charge * self.mol).sum()
    
    
    
    