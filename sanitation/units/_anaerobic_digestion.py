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

TODO:
    [1] Incorporate ADM, or change this to SimpleAD or something

'''


# %%

import numpy as np
from .. import SanUnit
from ..utils.loading import load_data, data_path

__all__ = ('AnaerobicDigestion',)

data_path += 'unit_data/AnaerobicDigestion.csv'

#!!! Potentially make these a Process, or make a Treatment class
from ._toilet import Toilet
_allocate_N_reduction = Toilet._allocate_N_reduction


# @cost(basis, CE, cost, n)
class AnaerobicDigestion(SanUnit):
    '''Anaerobic digestion of wastes with production of biogas.'''
    
    def __init__(self, ID='', ins=None, outs=(),
                 if_CH4_captured=True, if_N_degradation=True,
                 **kwargs):
        '''

        Parameters
        ----------
        if_CH4_captured : [bool]
            If the generated CH4 is captured.
        if_N_degradation : [bool]
            If N degradation and N2O emission occur during treatment.

        '''
        SanUnit.__init__(self, ID, ins, outs)
        self.if_CH4_captured = if_CH4_captured
        self.if_N_degradation = if_N_degradation
        self._tau_previous = 0.
        self._max_CH4_emission = 0.25
    
        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, '_'+para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    _N_ins = 1
    _N_outs = 3

    def _run(self):
        waste = self.ins[0]
        treated, CH4, N2O = self.outs
        treated.copy_like(self.ins[0])
        CH4.phase, N2O.phase = 'g', 'g'
        
        # COD removal
        COD_deg = treated._COD*treated.F_vol/1e3*self.COD_removal # kg/hr
        treated._COD *= (1-self.COD_removal)
        treated.mass *= (1-self.COD_removal)
        CH4_prd = COD_deg*self.MCF*self.max_CH4_emission
        if self.if_CH4_captured:
            CH4.empty()
        else:
            CH4.imass['CH4'] = CH4_prd

        #!!! Check with Hannah about this algorithm, previous storage time,
        # additional storage time, etc.
        # Maybe allow tau_previous = 0, otherwise add ValueError
        if self.if_N_degradation:
            N_deg = waste.TN*self.N_max_removal
            k = self.decay_k_N
            t0 = self.tau_previous/365
            t = self.tau/365
            N0 = N_deg*k*t0/(1-np.exp(-k*t0)) #!!! This definitely isn't correct
            N1 = N_deg*np.exp(-k*(t0+t))
            N_loss = N0 - N1
            N2O.imass['N2O'] = N_loss*waste.F_vol/1e3*self.N2O_EF*44/28
            NH3_rmd, NonNH3_rmd = \
                _allocate_N_reduction(N_loss*waste.F_vol/1e3, waste.imass['NH3'])
            treated.imass ['NH3'] = waste.imass['NH3'] - NH3_rmd
            treated.imass['NonNH3'] = waste.imass['NonNH3'] - NonNH3_rmd
        else:
            N2O.empty()

    def _design(self):
        D = self.design_results
        D['Flow rate'] = Q = self.ins[0].F_vol
        D['Reactor number'] = N = self.N_reactor
        V_tot = Q * self.tau*24
        # one extra as a backup
        D['Single reactor volume'] = V_single = V_tot/(1-self.headspace_frac)/(N-1)
        # Rx modeled as a cylinder
        D['Reactor diameter'] = dia = (4*V_single*self.aspect_ratio/np.pi)**(1/3)
        D['Reactor height'] = h = self.aspect_ratio * dia
        D['Reactor thickness'] = thick = self.rx_thickness
        D['Total concrete volume'] = V_concrete = N*thick*(2*np.pi/4*(dia**2)+np.pi*dia*h)
        D['Total excavation volume'] = V_tot
        #!!! Excavation is an activity, not material, maybe change the name to "construction"?
        
    # #!!! No opex assumption in ref [1]
    # # Use the Material/Construction class
    # def _cost(self):
    #     self.purchase_cost['Concrete'] = self.design_results['Total concrete volume']
        
        
        

    @property
    def COD_removal(self):
        '''[float] Fraction of COD removed during anaerobic digestion.'''
        return self._COD_removal
    @COD_removal.setter
    def COD_removal(self, i):
        self._COD_removal = float(i)

    @property
    def max_CH4_emission(self):
        '''[float] Maximum methane emssion as a fraction of degraded COD, [g CH4/g COD].'''
        return self._max_CH4_emission
    @max_CH4_emission.setter
    def max_CH4_emission(self, i):
        self._max_CH4_emission = float(i)

    @property
    def MCF(self):
        '''[float] Methane correction factor for COD during treatment.'''
        return self._MCF
    @MCF.setter
    def MCF(self, i):
        self._MCF = float(i)

    @property
    def N_max_removal(self):
        '''[float] Maximum fraction of N removed through denitrification during storage given sufficient time.'''
        return self._N_max_removal
    @N_max_removal.setter
    def N_max_removal(self, i):
        self._N_max_removal = float(i)

    @property
    def tau_previous(self):
        '''[float] Time between the waste production and anaerobic digestion, [d].'''
        return self._tau_previous
    @tau_previous.setter
    def tau_previous(self, i):
        self._tau_previous = float(i)

    @property
    def tau(self):
        '''[float] Storage time of the waste prior to anaerobic digestion, [d].'''
        return self._tau
    @tau.setter
    def tau(self, i):
        self._tau = float(i)

    @property
    def N2O_EF(self):
        '''[float] Fraction of N emitted as N2O.'''
        return self._N2O_EF
    @N2O_EF.setter
    def N2O_EF(self, i):
        self._N2O_EF = float(i)

    @property
    def decay_k_N(self):
        '''[float] Rate constant for N decay, [/yr].'''
        return self._decay_k_N
    @decay_k_N.setter
    def decay_k_N(self, i):
        self._decay_k_N = float(i)

    @property
    def N_reactor(self):
        '''[int] Number of reactors, non-integer will be converted to the smallest integer.'''
        return self._N_reactor
    @N_reactor.setter
    def N_reactor(self, i):
        self._N_reactor = np.ceil(i)

    @property
    def aspect_ratio(self):
        '''[float] Diameter-to-height ratio of the reactor.'''
        return self._aspect_ratio
    @aspect_ratio.setter
    def aspect_ratio(self, i):
        self._aspect_ratio = float(i)

    @property
    def headspace_frac(self):
        '''[float] Fraction of the reactor volume for headspace gas.'''
        return self._headspace_frac
    @headspace_frac.setter
    def headspace_frac(self, i):
        self._headspace_frac = float(i)

    @property
    def rx_thickness(self):
        '''[float] Wall thickness of the concrete reactor.'''
        return self._rx_thickness
    @rx_thickness.setter
    def rx_thickness(self, i):
        self._rx_thickness = float(i)











