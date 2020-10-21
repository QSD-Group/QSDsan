#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 07:48:41 2020

@author: yalinli_cabbi
"""


# import thermosteam as tmo
from thermosteam import Stream, utils


__all__ = ('WasteStream',)


# %%

# =============================================================================
# Define the WasteStream class
# =============================================================================

@utils.registered(ticket_name='ws')
class WasteStream(Stream):
    '''A subclass of the Stream object in the thermosteam package with additional attributes and methods for waste treatment    '''
    
    def show(self, T=None, P=None, flow='kg/hr', composition=None, N=None,
             stream_info=True):
        '''Show WasteStream information'''        
        info = ''

        # Stream-related specifications
        if stream_info:
            super().show(T, P, flow, composition, N)
        else:
            info += self._basic_info()
            display_units = self.display_units
            T_units = T or display_units.T
            P_units = P or display_units.P
            info += self._info_phaseTP(self.phase, T_units, P_units)
        
        # Component-related properties
        info += '\n Component-specific properties:\n'
        info += f'  TC      : {self.TC:.1f} g C/hr\n'
        info += f'  TN      : {self.TN:.1f} g N/hr\n'
        info += f'  TP      : {self.TP:.1f} g P/hr\n'
        info += f'  TK      : {self.TK:.1f} g K/hr\n'
        #!!! This display is definitely weired
        info += f'  charge  : {self.charge:.1f} mol/hr\n'
        info += f'  COD     : {self.COD:.1f} g O2/hr\n'
        info += f'  BOD5    : {self.BOD5:.1f} g O2/hr\n'
        info += f'  uBOD    : {self.uBOD:.1f} g O2/hr\n'
        info += f'  Totmass : {self.Totmass:.1f} g/hr\n'
        info += f'  Vmass   : {self.Vmass:.1f} g/hr\n'
        
        print(info)
        
    _ipython_display_ = show
    
    @property
    def components(self):
        return self._thermo.chemicals

    #!!! Suspect many of the units below aren't correct, just putting here as placeholders
    # double-check using self.mass or self.mol
    @property
    def TC(self):
        '''[float] Total carbon content of the stream in g C/hr'''
        return (self._thermo.chemicals.i_C * self.mass).sum()
    
    @property
    def TN(self):
        '''[float] Total nitrogen content of the stream in g N/hr'''
        return (self._thermo.chemicals.i_N * self.mass).sum()

    @property
    def TP(self):
        '''[float] Total phosphorus content of the stream in g P/hr'''
        return (self._thermo.chemicals.i_P * self.mass).sum()

    @property
    def TK(self):
        '''[float] Total potassium content of the stream in g K/hr'''
        return (self._thermo.chemicals.i_K * self.mass).sum()

    @property
    def charge(self):
        '''[float] Total charge of the stream in + mol/hr'''
        return (self._thermo.chemicals.i_charge * self.mol).sum()
    
    @property
    def COD(self):
        '''[float] Chemical oxygen demand of the stream in g O2/hr'''
        #!!! This is definitely not correct...
        return self.F_mass * 1000

    @property
    def BOD5(self):
        '''[float] Five-day biological oxygen demand of the stream in g O2/hr'''
        return (self._thermo.chemicals.f_BOD5_COD * self.COD).sum()

    @property
    def uBOD(self):
        '''[float] Ultimate biological oxygen demand of the stream in g O2/hr'''
        return (self._thermo.chemicals.f_uBOD_COD * self.COD).sum()
    

    #!!! How is this different from self.F_mass?
    @property
    def Totmass(self):
        '''[float] Total mass of the stream in g /hr'''
        return (self._thermo.chemicals.i_mass * self.mass).sum()
    
    @property
    def Vmass(self):
        '''[float] Total volatile mass in g /hr'''
        return (self._thermo.chemicals.f_Vmass_Totmass * self.Totmass).sum()














    
    