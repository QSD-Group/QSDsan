#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 17 07:48:41 2020

@author: yalinli_cabbi
"""

import thermosteam as tmo

__all__ = ('WasteStream',)


# %%

# =============================================================================
# Define the WasteStream class
# =============================================================================

class WasteStream(tmo.Stream):
    '''
    A subclass of the Stream object in the thermosteam package with additional attributes and methods for waste treatment
    '''

    #!!! Below are just copied from Component, haven't updated yet
    def __new__(cls, ID, formula=None, phase='l', i_C=0, i_N=0, i_P=0, i_K=0, i_mass=0,
                i_charge=0, f_BOD5_COD=0, f_BODu_COD=0, f_Vmass_Totmass=0,
                description=None, particle_size='Soluble', degradability='Undegradable',
                organic=False, measure_unit=None):

        self = super().__new__(cls, ID=ID, search_db=False)
        self._ID = ID
        if formula:
            self._formula = formula
        tmo._chemical.lock_phase(self, phase)
        self._i_C = i_C
        self._i_N = i_N
        self._i_P = i_P
        self._i_K = i_K
        self._i_mass = i_mass
        self._i_charge = i_charge
        self._f_BOD5_COD = f_BOD5_COD
        self._f_BODu_COD = f_BODu_COD
        self._f_Vmass_Totmass = f_Vmass_Totmass
        self._description = description
        assert particle_size in ('Dissolved gas', 'Soluble', 'Colloidal', 'Particulate'), \
            "particle_size must be 'Dissolved gas', 'Soluble', 'Colloidal', or 'Particulate'"
        self._particle_size = particle_size
        assert degradability in ('Biological', 'Chemical', 'Undegradable'), \
            "degradability must be 'Biological', 'Chemical', or 'Undegradable'"
        self._degradability = degradability
        assert organic.__class__ is bool, "organic must be 'True' or 'False'"
        self._organic = organic
        self._measure_unit = measure_unit
        return self