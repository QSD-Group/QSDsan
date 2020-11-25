#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Sanitation Explorer: Sustainable design of non-sewered sanitation technologies
Copyright (C) 2020, Sanitation Explorer Development Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-for-WaSH/sanitation/blob/master/LICENSE.txt
for license details.

Ref:
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.

'''


# %%

import numpy as np
from warnings import warn
from .. import SanUnit
from ._decay import Decay
from ..utils.loading import load_data, data_path

__all__ = ('Lagoon',)


class Lagoon(SanUnit, Decay):
    '''Anaerobic and facultive lagoon treatment.'''
    
    def __init__(self, ID='', ins=None, outs=(), design_type='anaerobic',
                 if_N2O_emission=True, **kwargs):
        
        '''

        Parameters
        ----------
        ins : WasteStream
            Waste for treatment.
        outs : WasteStream
            Treated waste, fugitive CH4, and fugitive N2O.
        design_type : [str]
            Can be 'anaerobic' or 'facultive'.
        if_N2O_emission : [bool]
            If consider N2O emission from N degradation the process.

        '''        
        
        SanUnit.__init__(self, ID, ins, outs)
        self._tau = None
        self._P_removal = 0.
        
        anaerobic_path = data_path + 'unit_data/AnaerobicLagoon.csv'
        self._anaerobic_defaults = load_data(path=anaerobic_path)
        facultive_path = data_path + 'unit_data/FacultiveLagoon.csv'
        self._facultive_defaults = load_data(path=facultive_path)
        
        self._design_type = None
        self.design_type = design_type
        self.if_N2O_emission = if_N2O_emission
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    __doc__ += __init__.__doc__
    __init__.__doc__ = __doc__
        
    
    _N_ins = 1
    _N_outs = 3

    def _run(self):
        waste = self.ins[0]
        treated, CH4, N2O = self.outs
        CH4.phase = N2O.phase = 'g'

        treated.copy_like(waste)
        removed_frac = self.COD_removal*self.COD_degradation
        treated._COD *= 1 - self.COD_removal
        treated.imass['OtherSS'] *= 1 - self.COD_removal
        
        CH4.imass['CH4'] = \
            removed_frac*waste.F_vol/1e3*self.MCF_decay*self.max_CH4_emission
        
        if self.if_N2O_emission:
            N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.tau/365,
                                            max_removal=self.N_max_removal)
            N_loss_tot = N_loss*waste.TN/1e3*waste.F_vol
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(N_loss_tot, waste.imass['NH3'])
            treated.imass ['NH3'] -=  NH3_rmd
            treated.imass['NonNH3'] -= NonNH3_rmd
            N2O.imass['N2O'] = N_loss_tot*self.N2O_EF_decay*44/28
        else:
            N2O.empty()

        treated.imass['P'] *= 1 - self.P_removal

    def _cost(self):
        design = self.design_results
        design['Lagoon number'] = N = self.N_lagoon
        design['Single lagoon volume'] = V = self.lagoon_V
        design['Lagoon length'] = L = self.lagoon_L
        design['Lagoon width'] = W = self.lagoon_W
        design['Lagoon depth'] = depth = V / (L*W)
        design['Total excavation'] = V * N
        liner_A = (L*W + 2*depth*(L+W)) * N
        design['Total liner'] = liner_A * self.liner_unit_mass
    
    @property
    def design_type(self):
        '''[str] Lagoon type, can be either 'anaerobic' or 'facultive'.'''
        return self._design_type
    @design_type.setter
    def design_type(self, i):
        if i == self._design_type: pass
        else:
            if i == 'anaerobic':
                data = self._anaerobic_defaults
            elif i == 'facultive':
                data = self._facultive_defaults
            else:
                raise ValueError("design_type can only be 'anaerobic' or 'facultive',"
                                 f'not {i}.')
            for para in data.index:
                value = float(data.loc[para]['expected'])
                setattr(self, para, value)
        self._design_type = i

    @property
    def COD_removal(self):
        '''[float] Fraction of COD removed during treatment.'''
        return self._COD_removal
    @COD_removal.setter
    def COD_removal(self, i):
        self._COD_removal = float(i)
    
    @property
    def COD_degradation(self):
        '''[float] Fraction of removed COD that degrades.'''
        return self._COD_degradation
    @COD_degradation.setter
    def COD_degradation(self, i):
        self._COD_degradation = float(i)
        
    @property
    def P_removal(self):
        '''[float] Fraction of P removed during treatment.'''
        return self._P_removal
    @P_removal.setter
    def P_removal(self, i):
        self._P_removal = float(i)
        
    @property
    def N_lagoon(self):
        '''[int] Number of lagoons, float will be converted to the smallest integer.'''
        return self._N_lagoon
    @N_lagoon.setter
    def N_lagoon(self, i):
        self._N_lagoon = np.ceil(i)    

    @property
    def tau(self):
        '''[float] Residence time, [d].'''
        if self._lagoon_V:
            return self._lagoon_V*self.N_lagoon/(self.F_mass_in*24)
        else:
            return self._tau
    @tau.setter
    def tau(self, i):
        if self._lagoon_V:
            msg = f'Residence time set, the original lagoon volume of {self._lagoon_V} m3 is ignored'
            warn(msg, source=self)
            self._lagoon_V = None
        self._tau = float(i)

    @property
    def lagoon_V(self):
        '''[float] Volume of the lagoon, [m3].'''
        if self._tau:
            return self._tau*self.F_mass_in*24/self.N_lagoon
        else:
            return self._lagoon_V
    @lagoon_V.setter
    def lagoon_V(self, i):
        if self._tau:
            msg = f'Lagoon volume set, the original residence time of {self._tau} d is ignored'
            warn(msg, source=self)
            self._tau = None
        self._lagoon_V = float(i)

    @property
    def lagoon_L(self):
        '''[float] Length of the lagoon, [m].'''
        return self._lagoon_L
    @lagoon_L.setter
    def lagoon_L(self, i):
        self._lagoon_L = float(i)

    @property
    def lagoon_W(self):
        '''[float] Width of the lagoon, [m].'''
        return self._lagoon_W
    @lagoon_W.setter
    def lagoon_W(self, i):
        self._lagoon_W = float(i)

    @property
    def liner_unit_mass(self):
        '''[float] Unit mass of the lagoon liner, [kg/m2].'''
        return self._liner_unit_mass
    @liner_unit_mass.setter
    def liner_unit_mass(self, i):
        self._liner_unit_mass = float(i)











