#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


# %%

import numpy as np
from .. import SanUnit, Construction
from ._decay import Decay
from ..utils import load_data, data_path

__all__ = ('AnaerobicBaffledReactor',)

data_path += 'sanunit_data/_anaerobic_baffled_reactor.tsv'


class AnaerobicBaffledReactor(SanUnit, Decay):
    '''
    Anaerobic baffled reactor with the production of biogas based on
    `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_

    Parameters
    ----------
    ins : WasteStream
        Waste for treatment.
    outs : WasteStream
        Treated waste, biogas, fugitive CH4, and fugitive N2O.
    degraded_components : tuple
        IDs of components that will degrade (at the same removal as `COD_removal`).
    if_capture_biogas : bool
        If produced biogas will be captured, otherwise it will be treated
        as fugitive CH4.
    if_N2O_emission : bool
        If considering fugitive N2O generated from the degraded N.

    Examples
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_

    References
    ----------
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
    Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
    Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
    https://doi.org/10.1021/acs.est.0c03296.

    See Also
    --------
    :ref:`qsdsan.sanunits.Decay <sanunits_Decay>`
    '''

    gravel_density = 1600

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), if_capture_biogas=True,
                 if_N2O_emission=False, **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=1)
        self.degraded_components = tuple(degraded_components)
        self.if_capture_biogas = if_capture_biogas
        self.if_N2O_emission = if_N2O_emission

        self.construction = (
            Construction('concrete', item='Concrete', quantity_unit='m3'),
            Construction('gravel', item='Gravel', quantity_unit='kg'),
            Construction('excavation', item='Excavation', quantity_unit='m3'),
            )

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
        COD_deg = waste._COD*waste.F_vol/1e3*self.COD_removal # kg/hr
        treated._COD *= (1-self.COD_removal)
        treated.imass[self.degraded_components] *= (1-self.COD_removal)

        CH4_prcd = COD_deg*self.MCF_decay*self.max_CH4_emission
        if self.if_capture_biogas:
            biogas.imass['CH4'] = CH4_prcd
            CH4.empty()
        else:
            CH4.imass['CH4'] = CH4_prcd
            biogas.empty()

        N_tot = waste.TN/1e3 * waste.F_vol
        N_loss_tot = N_tot * self.N_removal
        NH3_rmd, NonNH3_rmd = \
            self.allocate_N_removal(N_loss_tot, waste.imass['NH3'])
        treated.imass ['NH3'] = waste.imass['NH3'] - NH3_rmd
        treated.imass['NonNH3'] = waste.imass['NonNH3'] - NonNH3_rmd

        if self.if_N2O_emission:
            N2O.imass['N2O'] = N_loss_tot*self.N_max_decay*self.N2O_EF_decay*44/28
        else:
            N2O.empty()

    _units = {
        'Residence time': 'd',
        'Reactor length': 'm',
        'Reactor width': 'm',
        'Reactor height': 'm',
        'Single reactor volume': 'm3'
        }

    def _design(self):
        design = self.design_results
        design['Residence time'] = self.tau
        design['Reactor number'] = N = self.N_reactor
        design['Baffle number'] = N_b = self.N_baffle
        design['Reactor length'] = L = self.reactor_L
        design['Reactor width'] = W = self.reactor_W
        design['Reactor height'] = H = self.reactor_H
        design['Single reactor volume'] = V = L*W*H

        constr = self.construction
        concrete = N*self.concrete_thickness*(2*L*W+2*L*H+(2+N_b)*W*H)*self.add_concrete
        constr[0].quantity = concrete
        constr[1].quantity = N*V/(N_b+1) * self.gravel_density
        constr[2].quantity = N * V # excavation

        self.add_construction()


    @property
    def tau(self):
        '''[float] Residence time, [d].'''
        return self._tau
    @tau.setter
    def tau(self, i):
        self._tau = i

    @property
    def COD_removal(self):
        '''[float] Fraction of COD removed during treatment.'''
        return self._COD_removal
    @COD_removal.setter
    def COD_removal(self, i):
        self._COD_removal = i

    @property
    def N_removal(self):
        '''[float] Fraction of N removed during treatment.'''
        return self._N_removal
    @N_removal.setter
    def N_removal(self, i):
        self._N_removal = i

    @property
    def N_reactor(self):
        '''[int] Number of reactors, float will be converted to the smallest integer.'''
        return self._N_reactor
    @N_reactor.setter
    def N_reactor(self, i):
        self._N_reactor = int(np.ceil(i))

    @property
    def reactor_L(self):
        '''[float] Reactor length, [m].'''
        return self._reactor_L
    @reactor_L.setter
    def reactor_L(self, i):
        self._reactor_L = i

    @property
    def reactor_W(self):
        '''[float] Reactor width, [m].'''
        return self._reactor_W
    @reactor_W.setter
    def reactor_W(self, i):
        self._reactor_W = i

    @property
    def reactor_H(self):
        '''[float] Reactor height, [m].'''
        return self._reactor_H
    @reactor_H.setter
    def reactor_H(self, i):
        self._reactor_H = i

    @property
    def N_baffle(self):
        '''[int] Number of reactors, float will be converted to the smallest integer.'''
        return self._N_baffle
    @N_baffle.setter
    def N_baffle(self, i):
        self._N_baffle = int(np.ceil(i))

    @property
    def add_concrete(self):
        '''
        [float] Additional concrete as a fraction of the reactor concrete usage
        to account for receiving basin and biogas tank.
        '''
        return self._add_concrete
    @add_concrete.setter
    def add_concrete(self, i):
        self._add_concrete = i

    @property
    def concrete_thickness(self):
        '''[float] Thickness of the concrete wall.'''
        return self._concrete_thickness
    @concrete_thickness.setter
    def concrete_thickness(self, i):
        self._concrete_thickness = i