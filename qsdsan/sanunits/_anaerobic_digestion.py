#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.

TODO (maybe):
    [1] Incorporating ADM as a Process, or change this to SimpleAD or similar

'''


# %%

import numpy as np
from .. import SanUnit, Construction
from ._decay import Decay
from ..utils.loading import load_data, data_path

__all__ = ('AnaerobicDigestion',)

data_path += 'sanunit_data/_anaerobic_digestion.tsv'


class AnaerobicDigestion(SanUnit, Decay):
    '''
    Anaerobic digestion of wastes with the production of biogas based on Trimmer et al. [1]_

    Parameters
    ----------
    ins : WasteStream
        Waste for treatment.
    outs : WasteStream
        Treated waste, biogas, and fugitive N2O.
    flow_rate : float
        Total flow rate through the reactor (for sizing purpose), [m3/d].
        If not provided, will use F_vol_in.
    if_capture_biogas : bool
        If produced biogas will be captured, otherwise it will be treated
        as fugitive CH4.
    if_N2O_emission : bool
        If consider N2O emission from N degradation the process.
        
    References
    ----------
    .. [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.
        
    See Also
    --------
    :ref:`qsdsan.sanunits.Decay <sanunits_Decay>`
    
    '''
    
    def __init__(self, ID='', ins=None, outs=(), flow_rate=None,
                 if_capture_biogas=True, if_N2O_emission=False, **kwargs):
        SanUnit.__init__(self, ID, ins, outs)
        self._flow_rate = flow_rate
        self.if_capture_biogas = if_capture_biogas
        self.if_N2O_emission = if_N2O_emission
    
        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, '_'+para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    _N_ins = 1
    _N_outs = 4

    def _run(self):
        waste = self.ins[0]
        treated, biogas, CH4, N2O = self.outs
        treated.copy_like(self.ins[0])
        biogas.phase = CH4.phase = N2O.phase = 'g'
        
        # COD removal
        COD_deg = treated._COD*treated.F_vol/1e3*self.COD_removal # kg/hr
        treated._COD *= (1-self.COD_removal)
        
        #!!! Which assumption is better?
        treated.imass['OtherSS'] *= (1-self.COD_removal)
        # treated.mass *= (1-self.COD_removal)
        
        CH4_prcd = COD_deg*self.MCF_decay*self.max_CH4_emission
        if self.if_capture_biogas:
            biogas.imass['CH4'] = CH4_prcd
            CH4.empty()
        else:
            CH4.imass['CH4'] = CH4_prcd
            biogas.empty()

        if self.if_N2O_emission:
            N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.tau/365,
                                            max_decay=self.N_max_decay)
            N_loss_tot = N_loss*waste.TN/1e3*waste.F_vol
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(N_loss_tot, waste.imass['NH3'])
            treated.imass ['NH3'] = waste.imass['NH3'] - NH3_rmd
            treated.imass['NonNH3'] = waste.imass['NonNH3'] - NonNH3_rmd
            N2O.imass['N2O'] = N_loss_tot*self.N2O_EF_decay*44/28
        else:
            N2O.empty()

    _units = {
        'Volumetric flow rate': 'm3/hr',
        'Residence time': 'd',
        'Single reactor volume': 'm3',
        'Reactor diameter': 'm',
        'Reactor height': 'm'
        }

    def _design(self):
        design = self.design_results
        design['Volumetric flow rate'] = Q = self.flow_rate
        design['Residence time'] = tau = self.tau
        design['Reactor number'] = N = self.N_reactor
        V_tot = Q * tau*24
        # one extra as a backup
        design['Single reactor volume'] = V_single = V_tot/(1-self.headspace_frac)/(N-1)
        # Rx modeled as a cylinder
        design['Reactor diameter'] = D = (4*V_single*self.aspect_ratio/np.pi)**(1/3)
        design['Reactor height'] = H = self.aspect_ratio * D
        concrete =  N*self.concrete_thickness*(2*np.pi/4*(D**2)+np.pi*D*H)
        self.construction = (
            Construction(item='Concrete', quantity=concrete, quantity_unit='m3'),
            Construction(item='Excavation', quantity=V_tot, quantity_unit='m3'),
            )
        self.add_construction()
        

    @property
    def flow_rate(self):
        '''
        [float] Total flow rate through the reactor (for sizing purpose), [m3/d].
        If not provided, will calculate based on F_vol_in.
        '''
        return self._flow_rate if self._flow_rate else self.F_vol_in*24
    @flow_rate.setter
    def flow_rate(self, i):
        self._flow_rate = float(i)

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
        self._N_reactor = int(np.ceil(i))

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
    def concrete_thickness(self):
        '''[float] Thickness of the concrete wall.'''
        return self._concrete_thickness
    @concrete_thickness.setter
    def concrete_thickness(self, i):
        self._concrete_thickness = float(i)











