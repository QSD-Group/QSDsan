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
from ._decay import Decay
from ..utils.loading import load_data, data_path

__all__ = ('AnaerobicDigestion',)

data_path += 'unit_data/_anaerobic_digestion.csv'


# @cost(basis, CE, cost, n)
class AnaerobicDigestion(SanUnit, Decay):
    '''Anaerobic digestion of wastes with production of biogas.'''
    
    def __init__(self, ID='', ins=None, outs=(),
                 if_CH4_captured=True, if_N_degradation=True,
                 **kwargs):
        
        '''

        Parameters
        ----------
        ins : WasteStream
            Waste for treatment.
        outs : WasteStream
            Treated waste, fugitive CH4, and fugitive N2O.
        if_CH4_captured : [bool]
            If the generated CH4 is captured.
        if_N_degradation : [bool]
            If N degradation and N2O emission occur during treatment.

        '''
        
        SanUnit.__init__(self, ID, ins, outs)
        self.if_CH4_captured = if_CH4_captured
        self.if_N_degradation = if_N_degradation
        self._tau_previous = 0.
    
        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, '_'+para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    __doc__ += __init__.__doc__
    __init__.__doc__ = __doc__
    
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
        CH4_prcd = COD_deg*self.MCF_decay*self.max_CH4_emission
        if self.if_CH4_captured:
            CH4.empty()
        else:
            CH4.imass['CH4'] = CH4_prcd

        #!!! Check with Hannah about this algorithm, previous storage time,
        # additional storage time, etc.
        if self.if_N_degradation:
            N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.tau/365,
                                            max_removal=self.N_max_removal)
            N_loss_tot = N_loss*waste.TN/1e3*waste.F_vol
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(N_loss_tot, waste.imass['NH3'])
            treated.imass ['NH3'] = waste.imass['NH3'] - NH3_rmd
            treated.imass['NonNH3'] = waste.imass['NonNH3'] - NonNH3_rmd
            N2O.imass['N2O'] = N_loss_tot*self.N2O_EF_decay*44/28
        else:
            N2O.empty()

    def _design(self):
        design = self.design_results
        design['Flow rate'] = Q = self.ins[0].F_vol
        design['Reactor number'] = N = self.N_reactor
        V_tot = Q * self.tau*24
        # one extra as a backup
        design['Single reactor volume'] = V_single = V_tot/(1-self.headspace_frac)/(N-1)
        # Rx modeled as a cylinder
        design['Reactor diameter'] = D = (4*V_single*self.aspect_ratio/np.pi)**(1/3)
        design['Reactor height'] = H = self.aspect_ratio * D
        design['Reactor thickness'] = thick = self.rx_thickness
        design['Total concrete volume'] = N*thick*(2*np.pi/4*(D**2)+np.pi*D*H)
        design['Total excavation volume'] = V_tot
        #!!! Excavation is an activity, not material, maybe change the name to "construction"?
        
    #!!! No opex assumption in ref [1]
    # Use the Material/Construction class
    def _cost(self):
        pass
        # self.purchase_cost['Concrete'] = self.design_results['Total concrete volume']
        
    
    #!!! Maybe do not need this tau_previous
    @property
    def tau_previous(self):
        '''[float] Time between the waste production and anaerobic digestion, [d].'''
        return self._tau_previous
    @tau_previous.setter
    def tau_previous(self, i):
        self._tau_previous = float(i)

    @property
    def tau(self):
        '''[float] Residence time, [d].'''
        return self._tau
    @tau.setter
    def tau(self, i):
        self._tau = float(i)

    @property
    def COD_removal(self):
        '''[float] Fraction of COD removed during treatment.'''
        return self._COD_removal
    @COD_removal.setter
    def COD_removal(self, i):
        self._COD_removal = float(i)

    @property
    def N_reactor(self):
        '''[int] Number of reactors, float will be converted to the smallest integer.'''
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











