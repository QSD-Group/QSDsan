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

from math import ceil
from warnings import warn
from ... import SanUnit
from .._decay import Decay
from ...utils import ospath, load_data, data_path

__all__ = ('Toilet',)

toilet_path = ospath.join(data_path, 'sanunit_data/_toilet.tsv')


# %%

class Toilet(SanUnit, Decay, isabstract=True):
    '''
    Abstract class containing common parameters and design algorithms for toilets
    based on `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_

    Parameters
    ----------
    degraded_components : tuple
        IDs of components that will degrade (simulated by first-order decay).
    N_user : int, float
        Number of people per toilet.
        Note that this number can be a float when calculated from `N_tot_user` and `N_toilet`.
    N_toilet : int
        Number of parallel toilets.
        In calculation, `N_toilet` will be calculated as `ceil(N_tot_user/N_user)`.
    N_tot_user : int
        Total number of users.

        .. note::

            If `N_tot_user` is provided (i.e., not "None"),
            then updating `N_user` will recalculate `N_toilet`, and vice versa.

    if_toilet_paper : bool
        If toilet paper is used.
    if_flushing : bool
        If water is used for flushing.
    if_cleansing : bool
        If water is used for cleansing.
    if_desiccant : bool
        If desiccant is used for moisture and odor control.
    if_air_emission : bool
        If emission to air occurs
        (i.e., if the pit is completely sealed off from the atmosphere).
    if_ideal_emptying : bool
        If the toilet appropriately emptied to avoid contamination to the
        environmental.
    CAPEX : float
        Capital cost of a single toilet.
    OPEX_over_CAPEX : float
        Fraction of annual operating cost over total capital cost.

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

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), N_user=1, N_toilet=1, N_tot_user=None,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                 CAPEX=None, OPEX_over_CAPEX=None):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=1)
        self.degraded_components = tuple(degraded_components)
        self._N_user = self._N_toilet = self._N_tot_user = None
        self.N_user = N_user
        self.N_toilet = N_toilet
        self.N_tot_user = N_tot_user
        self.if_toilet_paper = if_toilet_paper
        self.if_flushing = if_flushing
        self.if_cleansing = if_cleansing
        self.if_desiccant = if_desiccant
        self.if_air_emission = if_air_emission
        self.if_ideal_emptying = if_ideal_emptying
        self.CAPEX = CAPEX
        self.OPEX_over_CAPEX = OPEX_over_CAPEX

        data = load_data(path=toilet_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            if para in ('desiccant_V', 'desiccant_rho'):
                setattr(self, para, value)
            else:
                setattr(self, '_'+para, value)
        del data

        self._empty_ratio = 0.59

    _N_ins = 6
    _outs_size_is_fixed = False


    def _run(self):
        ur, fec, tp, fw, cw, des = self.ins
        tp.imass['Tissue'] = int(self.if_toilet_paper)*self.toilet_paper
        fw.imass['H2O'] = int(self.if_flushing)*self.flushing_water
        cw.imass['H2O'] = int(self.if_cleansing)*self.cleansing_water
        des.imass['WoodAsh'] = int(self.if_desiccant)*self.desiccant

    density_dct = {
        'Sand': 1442,
        'Gravel': 1600,
        'Brick': 1750,
        'Plastic': 0.63,
        'Steel': 7900,
        'StainlessSteelSheet': 2.64
        }

    def _cost(self):
        self.baseline_purchase_costs['Total toilets'] = self.CAPEX * self.N_toilet
        add_OPEX = self.baseline_purchase_costs['Total toilets']*self.OPEX_over_CAPEX/365/24
        self._add_OPEX = {'Additional OPEX': add_OPEX}


    @staticmethod
    def get_emptying_emission(waste, CH4, N2O, empty_ratio, CH4_factor, N2O_factor):
        '''
        Calculate emissions due to non-ideal emptying based on
        `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_,

        Parameters
        ----------
        stream : WasteStream
            Excreta stream that is not appropriately emptied (before emptying).
        CH4 : WasteStream
            Fugitive CH4 gas (before emptying).
        N2O : WasteStream
            Fugitive N2O gas (before emptying).
        empty_ratio : float
            Fraction of excreta that is appropriately emptied..
        CH4_factor : float
            Factor to convert COD removal to CH4 emission.
        N2O_factor : float
            Factor to convert COD removal to N2O emission.

        Returns
        -------
        stream : WasteStream
            Excreta stream that is not appropriately emptied (after emptying).
        CH4 : WasteStream
            Fugitive CH4 gas (after emptying).
        N2O : WasteStream
            Fugitive N2O gas (after emptying).

        '''
        COD_rmd = waste.COD*(1-empty_ratio)/1e3*waste.F_vol
        CH4.imass['CH4'] += COD_rmd * CH4_factor
        N2O.imass['N2O'] += COD_rmd * N2O_factor
        waste.mass *= empty_ratio
        return waste, CH4, N2O

    @property
    def N_user(self):
        '''[int, float] Number of people per toilet.'''
        return self._N_user or self.N_tot_user/self.N_toilet
    @N_user.setter
    def N_user(self, i):
        if i is not None:
            N_user = self._N_user = int(i)
            old_toilet = self._N_toilet
            if old_toilet and self.N_tot_user:
                new_toilet = ceil(self.N_tot_user/N_user)
                warn(f'With the provided `N_user`, the previous `N_toilet` of {old_toilet} '
                     f'is recalculated from `N_tot_user` and `N_user` as {new_toilet}.')
                self._N_toilet = None
        else:
            self._N_user = i

    @property
    def N_toilet(self):
        '''[int] Number of parallel toilets.'''
        return self._N_toilet or ceil(self.N_tot_user/self.N_user)
    @N_toilet.setter
    def N_toilet(self, i):
        if i is not None:
            N_toilet = self._N_toilet = ceil(i)
            old_user = self._N_user
            if old_user and self.N_tot_user:
                new_user = self.N_tot_user/N_toilet
                warn(f'With the provided `N_toilet`, the previous `N_user` of {old_user} '
                     f'is recalculated from `N_tot_user` and `N_toilet` as {new_user}.')
                self._N_user = None
        else:
            self._N_toilet = i

    @property
    def N_tot_user(self):
        '''[int] Number of total users.'''
        return self._N_tot_user
    @N_tot_user.setter
    def N_tot_user(self, i):
        if i is not None:
            self._N_tot_user = int(i)
        else:
            self._N_tot_user = None

    @property
    def toilet_paper(self):
        '''
        [float] Amount of toilet paper used
        (if ``if_toilet_paper`` is True), [kg/cap/hr].
        '''
        return self._toilet_paper
    @toilet_paper.setter
    def toilet_paper(self, i):
        self._toilet_paper = i

    @property
    def flushing_water(self):
        '''
        [float] Amount of water used for flushing
        (if ``if_flushing_water`` is True), [kg/cap/hr].
        '''
        return self._flushing_water
    @flushing_water.setter
    def flushing_water(self, i):
        self._flushing_water = i

    @property
    def cleansing_water(self):
        '''
        [float] Amount of water used for cleansing
        (if ``if_cleansing_water`` is True), [kg/cap/hr].
        '''
        return self._cleansing_water
    @cleansing_water.setter
    def cleansing_water(self, i):
        self._cleansing_water = i

    @property
    def desiccant(self):
        '''
        [float] Amount of desiccant used (if ``if_desiccant`` is True), [kg/cap/hr].

        .. note::

            Value set by ``desiccant_V`` and ``desiccant_rho``.

        '''
        return self.desiccant_V*self.desiccant_rho

    @property
    def N_volatilization(self):
        '''
        [float] Fraction of input N that volatilizes to the air
        (if ``if_air_emission`` is True).
        '''
        return self._N_volatilization
    @N_volatilization.setter
    def N_volatilization(self, i):
        self._N_volatilization = i

    @property
    def empty_ratio(self):
        '''
        [float] Fraction of excreta that is appropriately emptied.

        .. note::

            Will be 1 (i.e., 100%) if ``if_ideal_emptying`` is True.

        '''
        if self.if_ideal_emptying:
            return 1.
        return self._empty_ratio
    @empty_ratio.setter
    def empty_ratio(self, i):
        if self.if_ideal_emptying:
            warn(f'`if_ideal_emptying` is True, the set value {i} is ignored.')
        self._empty_ratio = i

    @property
    def MCF_aq(self):
        '''[float] Methane correction factor for COD lost due to inappropriate emptying.'''
        return self._MCF_aq
    @MCF_aq.setter
    def MCF_aq(self, i):
        self._MCF_aq = i

    @property
    def N2O_EF_aq(self):
        '''[float] Fraction of N emitted as N2O due to inappropriate emptying.'''
        return self._N2O_EF_aq
    @N2O_EF_aq.setter
    def N2O_EF_aq(self, i):
        self._N2O_EF_aq = i