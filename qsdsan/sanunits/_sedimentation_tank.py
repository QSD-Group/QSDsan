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
from .. import Construction
from ._decay import Decay
from ._sludge_separator import SludgeSeparator
from ..utils import load_data, data_path

__all__ = ('SedimentationTank',)

data_path += 'sanunit_data/_sedimentation_tank.tsv'


class SedimentationTank(SludgeSeparator, Decay):
    '''
    Sedimentation of wastes into liquid and solid phases based on
    `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_

    Parameters
    ----------
    ins : WasteStream
        Waste for treatment.
    outs : WasteStream
        Liquid, settled solids, fugitive CH4, and fugitive N2O.
    degraded_components : tuple
        IDs of components that will degrade (simulated by first-order decay).
    if_N2O_emission : bool
        If consider N2O emission from N degradation in the process.

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

    def __init__(self, ID='', ins=None, outs=(),thermo=None, init_with='WasteStream',
                 split=None, settled_frac=None,
                 degraded_components=('OtherSS',), if_N2O_emission=False, **kwargs):

        SludgeSeparator.__init__(self, ID, ins, outs, thermo, init_with,
                                 split, settled_frac, F_BM_default=1)
        self.degraded_components = tuple(degraded_components)
        self.if_N2O_emission = if_N2O_emission

        self.construction = (
            Construction('concrete', item='Concrete', quantity_unit='m3'),
            Construction('steel', item='Steel', quantity_unit='kg'),
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
        liq, sol, CH4, N2O = self.outs
        CH4.phase = N2O.phase = 'g'

        # Retention in the settled solids
        SludgeSeparator._run(self)

        # COD degradation in settled solids
        COD_loss = self.first_order_decay(k=self.decay_k_COD,
                                          t=self.tau/365,
                                          max_decay=self.COD_max_decay)

        tot_COD_kg = sol._COD * sol.F_vol / 1e3
        sol.imass[self.degraded_components] *= 1 - COD_loss

        # Adjust total mass of of the settled solids by changing water content
        liq, sol = self._adjust_solid_water(waste, liq, sol)

        COD_loss_kg = tot_COD_kg * COD_loss
        CH4.imass['CH4'] = COD_loss_kg * self.max_CH4_emission * self.MCF_decay
        sol._COD = tot_COD_kg*(1-COD_loss)/sol.F_vol*1e3

        # N degradation
        if self.if_N2O_emission:
            N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.tau/365,
                                            max_decay=self.N_max_decay)
            N_loss_tot = N_loss*sol.TN/1e3*sol.F_vol
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(N_loss_tot, sol.imass['NH3'])
            sol.imass ['NH3'] -=  NH3_rmd
            sol.imass['NonNH3'] -= NonNH3_rmd
            N2O.imass['N2O'] = N_loss_tot*self.N2O_EF_decay*44/28
        else:
            N2O.empty()

    _units = {
        'Single tank volume': 'm3',
        'Single tank height': 'm',
        'Single tank width': 'm',
        'Single tank length': 'm',
        'Single roof area': 'm2'
        }

    def _design(self):
        design = self.design_results
        # `tau` not used, might be that working volume fraction not known
        design['Tank number'] = N = self.N_tank
        design['Single tank volume'] = V_single = self.tank_V
        L2W = self.tank_L_to_W
        W2H = self.tank_W_to_H
        design['Single tank height'] = H = (V_single/(L2W*(W2H**2)))**(1/3)
        design['Single tank width'] = W = H * W2H
        design['Single tank length'] = L = W * L2W
        design['Single roof area'] = N*L*W/(np.cos(self.roof_slope/180*np.pi))
        side_area = N*2*(L*H + W*H)

        # Concrete
        thick = self.concrete_thickness
        side_concrete = N*thick*(L*W+2*W*H+2*L*H)
        column_concrete = N*(thick**2)*H*self.column_per_side*2

        constr = self.construction
        constr[0].quantity = side_concrete + column_concrete
        constr[1].quantity = (design['Single roof area']+side_area) * self.roof_unit_mass # steel

        self.add_construction()

    @property
    def tau(self):
        '''[float] Residence time, [d].'''
        return self._tau
    @tau.setter
    def tau(self, i):
        self._tau = i

    @property
    def tank_V(self):
        '''[float] Volume of the sedimentation tank.'''
        return self._tank_V
    @tank_V.setter
    def tank_V(self, i):
        self._tank_V = i

    @property
    def tank_L_to_W(self):
        '''[float] Length-to-width ratio of the sedimentation tank.'''
        return self._tank_L_to_W
    @tank_L_to_W.setter
    def tank_L_to_W(self, i):
        self._tank_L_to_W = i

    @property
    def tank_W_to_H(self):
        '''[float] Width-to-height ratio of the sedimentation tank.'''
        return self._tank_W_to_H
    @tank_W_to_H.setter
    def tank_W_to_H(self, i):
        self._tank_W_to_H = i

    @property
    def N_tank(self):
        '''[int] Number of sedimentation tanks, float will be converted to the smallest integer.'''
        return self._N_tank
    @N_tank.setter
    def N_tank(self, i):
        self._N_tank = int(np.ceil(i))

    @property
    def column_per_side(self):
        '''[int] Number of columns per side of sedimentation tanks, float will be converted to the smallest integer.'''
        return self._column_per_side
    @column_per_side.setter
    def column_per_side(self, i):
        self._column_per_side = int(np.ceil(i))

    @property
    def concrete_thickness(self):
        '''[float] Thickness of the concrete wall.'''
        return self._concrete_thickness
    @concrete_thickness.setter
    def concrete_thickness(self, i):
        self._concrete_thickness = i

    @property
    def roof_slope(self):
        '''[float] Slope of the tank roof, [°].'''
        return self._roof_slope
    @roof_slope.setter
    def roof_slope(self, i):
        self._roof_slope = i

    @property
    def roof_unit_mass(self):
        '''[float] Unit mass of the tank roof, [kg/m2].'''
        return self._roof_unit_mass
    @roof_unit_mass.setter
    def roof_unit_mass(self, i):
        self._roof_unit_mass = i