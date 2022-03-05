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
from math import ceil
from warnings import warn
from . import Decay
from .. import SanUnit, WasteStream, Construction
from ..utils import ospath, load_data, data_path, dct_from_str

__all__ = ('Toilet', 'PitLatrine', 'UDDT',)


# %%

toilet_path = ospath.join(data_path, 'sanunit_data/_toilet.tsv')

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


# %%

pit_path = ospath.join(data_path, 'sanunit_data/_pit_latrine.tsv')

class PitLatrine(Toilet):
    '''
    Single pit latrine based on `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_,
    a subclass of :class:`qsdsan.sanunits.Toilet`.

    The following components should be included in system thermo object for simulation:
    Tissue, WoodAsh, H2O, NH3, NonNH3, P, K, CH4, N2O.

    The following impact items should be pre-constructed for life cycle assessment:
    Cement, Sand, Gravel, Brick, Plastic, Steel, Wood, Excavation.

    Parameters
    ----------
    ins : stream obj
        Excreta.
    outs : Iterable(stream obj)
        Recyclable mixed excreta, stream leached to soil, fugitive CH4, and fugitive N2O.
    lifetime : int
        Lifetime of this pit latrine, [yr].
    if_leaching : bool
        If infiltration to soil occurs
        (i.e., if the pit walls and floors are permeable).
    if_shared : bool
        If the toilet is shared.
    if_pit_above_water_table : bool
        If the pit is above local water table.

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
    :ref:`qsdsan.sanunits.Toilet <sanunits_toilets>`
    '''

    # Legacy code to add checkers
    # _P_leaching = Frac_D(name='P_leaching')

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), N_user=1, N_toilet=1, N_tot_user=None,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                 CAPEX=449, OPEX_over_CAPEX=0.05,
                 lifetime=8, if_leaching=True, if_shared=True,
                 if_pit_above_water_table=True, **kwargs):

        Toilet.__init__(self, ID, ins, outs, thermo, init_with,
                        degraded_components, N_user, N_toilet, N_tot_user,
                        if_toilet_paper, if_flushing, if_cleansing, if_desiccant,
                        if_air_emission, if_ideal_emptying, CAPEX, OPEX_over_CAPEX)

        self.lifetime = lifetime
        self.if_leaching = if_leaching
        self.if_pit_above_water_table = if_pit_above_water_table
        self.if_shared = if_shared
        self._pit_depth = 4.57 # m
        self._pit_area = 0.8 # m2
        self._liq_leaching = None
        self._mixed = WasteStream(f'{self.ID}_mixed')

        self.construction = (
            Construction('cement', linked_unit=self, item='Cement', quantity_unit='kg'),
            Construction('sand', linked_unit=self, item='Sand', quantity_unit='kg'),
            Construction('gravel', linked_unit=self, item='Gravel', quantity_unit='kg'),
            Construction('brick', linked_unit=self, item='Brick', quantity_unit='kg'),
            Construction('liner', linked_unit=self, item='Plastic', quantity_unit='kg'),
            Construction('steel', linked_unit=self, item='Steel', quantity_unit='kg'),
            Construction('wood', linked_unit=self, item='Wood', quantity_unit='m3'),
            Construction('excavation', linked_unit=self, item='Excavation', quantity_unit='m3'),
            )

        data = load_data(path=pit_path)
        for para in data.index:
            if para in ('MCF_decay', 'N2O_EF_decay'):
                value = dct_from_str(data.loc[para]['expected'], dtype='float')
            else:
                value = float(data.loc[para]['expected'])
            setattr(self, '_'+para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    _N_outs = 4

    def _run(self):
        Toilet._run(self)
        waste, leachate, CH4, N2O = self.outs
        CH4.phase = N2O.phase = 'g'

        mixed = self._mixed
        mixed.mix_from(self.ins)
        tot_COD_kg = sum(float(getattr(i, '_COD') or getattr(i, 'COD'))*i.F_vol for i in self.ins)/1e3

        # All composite variables in mg/L
        # Leaching
        # Here COD change due to leaching is not considered
        if self.if_leaching:
            # Additional assumption not in ref [1]
            leachate.imass['H2O'] = mixed.imass['H2O'] * self.liq_leaching
            leachate.imass['NH3'], leachate.imass['NonNH3'] = \
                self.allocate_N_removal(mixed.TN/1e3*mixed.F_vol*self.N_leaching,
                                        mixed.imass['NH3'])
            leachate.imass['P'] = mixed.imass['P'] * self.P_leaching
            leachate.imass['K'] = mixed.imass['K'] * self.K_leaching
            mixed.mass -= leachate.mass

        # Air emission
        if self.if_air_emission:
            # N loss due to ammonia volatilization
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(mixed.TN/1e3*mixed.F_vol*self.N_volatilization,
                                        mixed.imass['NH3'])
            mixed.imass ['NH3'] -= NH3_rmd
            mixed.imass['NonNH3'] -= NonNH3_rmd
            # Energy/N loss due to degradation
            COD_loss = self.first_order_decay(k=self.decay_k_COD,
                                              t=self.emptying_period,
                                              max_decay=self.COD_max_decay)
            COD_loss_kg = tot_COD_kg * COD_loss
            CH4.imass['CH4'] = COD_loss_kg * self.max_CH4_emission * self.MCF_decay
            mixed.imass[self.degraded_components] *= 1 - COD_loss

            N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.emptying_period,
                                            max_decay=self.N_max_decay)
            N_loss_tot = mixed.TN/1e3*mixed.F_vol*N_loss
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(N_loss_tot,
                                        mixed.imass['NH3'])
            mixed.imass ['NH3'] -= NH3_rmd
            mixed.imass['NonNH3'] -= NonNH3_rmd
            N2O.imass['N2O'] = N_loss_tot * self.N2O_EF_decay * 44/28
            mixed._COD = (tot_COD_kg-COD_loss_kg)*1e3/mixed.F_vol
        else:
            CH4.empty()
            N2O.empty()

        # Aquatic emission when not ideally emptied
        if not self.if_ideal_emptying:
            mixed, CH4, N2O = self.get_emptying_emission(
                waste=mixed, CH4=CH4, N2O=N2O,
                empty_ratio=self.empty_ratio,
                CH4_factor=self.COD_max_decay*self.MCF_aq*self.max_CH4_emission,
                N2O_factor=self.N2O_EF_decay*44/28)

        # Drain extra water, assume density of water to be 1 kg/L
        sludge = self.sludge_accum_rate/(365*24)
        diff = mixed.F_mass - sludge
        if diff > 0:
            mixed_COD = mixed._COD * mixed.F_vol
            mixed.imass['H2O'] -= diff
            mixed.imass['H2O'] = max(0, mixed.imass['H2O'])
            mixed._COD = mixed_COD / mixed.F_vol

        waste.copy_like(mixed)

        # Scale up the effluent based on the number of user per toilet and
        # toilet number
        N_tot_user = self.N_tot_user or self.N_toilet*self.N_user
        for i in self.outs:
            if not i.F_mass == 0:
                i.F_mass *= N_tot_user

    _units = {
        'Emptying period': 'yr',
        'Single pit volume': 'm3',
        'Single pit area': 'm2',
        'Single pit depth': 'm'
        }

    def _design(self):
        design = self.design_results
        design['Number of users per toilet'] = self.N_user
        design['Parallel toilets'] = N = self.N_toilet
        design['Emptying period'] = self.emptying_period
        design['Single pit volume'] = self.pit_V
        design['Single pit area'] = self.pit_area
        design['Single pit depth'] = self.pit_depth

        density = self.density_dct
        constr = self.construction
        constr[0].quantity = 700 * N # cement
        constr[1].quantity = 2.2 * density['Sand'] * N
        constr[2].quantity = 0.8 * density['Gravel'] * N
        constr[3].quantity = 54*0.0024 * density['Brick'] * N
        constr[4].quantity = 16 * density['Plastic'] * N
        constr[5].quantity = 0.00425 * density['Steel'] * N
        constr[6].quantity = 0.19 * N # wood
        constr[7].quantity = self.pit_V * N # excavation

        self.add_construction(add_cost=False)


    @property
    def pit_depth(self):
        '''[float] Depth of the pit, [m].'''
        return self._pit_depth
    @pit_depth.setter
    def pit_depth(self, i):
        self._pit_depth = i

    @property
    def pit_area(self):
        '''[float] Area of the pit, [m2].'''
        return self._pit_area
    @pit_area.setter
    def pit_area(self, i):
        self._pit_area = i

    # With baseline assumptions, this is about the same as the total volume of the excreta
    @property
    def pit_V(self):
        '''[float] Volume of the pit, [m3].'''
        return self.pit_area*self.pit_depth

    @property
    def emptying_period(self):
        '''[float] Time interval between pit emptying, [yr].'''
        return self._emptying_period
    @emptying_period.setter
    def emptying_period(self, i):
        self._emptying_period = i

    @property
    def sludge_accum_rate(self):
        '''[float] Sludge accumulation rate, [L/cap/yr].'''
        return self._sludge_accum_rate
    @sludge_accum_rate.setter
    def sludge_accum_rate(self, i):
        self._sludge_accum_rate = i

    @property
    def liq_leaching(self):
        '''
        [float] Fraction of input water that leaches to the soil
        (if if_leaching is True). If not set, then return the maximum of
        fraction of N, P, K leaching
        '''
        return self._liq_leaching or \
            max(self.N_leaching, self.P_leaching, self.K_leaching)
    @liq_leaching.setter
    def liq_leaching(self, i):
        self._liq_leaching = i

    @property
    def N_leaching(self):
        '''
        [float] Fraction of input N that leaches to the soil
        (if if_leaching is True).
        '''
        return self._N_leaching
    @N_leaching.setter
    def N_leaching(self, i):
        self._N_leaching = i
        # Legacy code to add checkers
        # @Frac_C(self)
        # def N_leaching(): return i

    @property
    def P_leaching(self):
        '''
        [float] Fraction of input P that leaches to the soil
        (if if_leaching is True).
        '''
        return self._P_leaching
    @P_leaching.setter
    def P_leaching(self, i):
        self._P_leaching = i

    @property
    def K_leaching(self):
        '''
        [float] Fraction of input K that leaches to the soil
        (if if_leaching is True).
        '''
        return self._K_leaching
    @K_leaching.setter
    # This is faster than using descriptors or decorators, but (I think) less elegant
    def K_leaching(self, i):
        if i < 0:
            raise ValueError('Value for K_leaching cannot be negative')
        self._K_leaching = i

    def _return_MCF_EF(self):
        # self._MCF and self._N2O_EF are dict for
        # single_above_water, communal_above_water, below_water
        if self.if_pit_above_water_table:
            if not self.if_shared:
                return 'single_above_water'
            else:
                return 'communal_above_water'
        else:
            return 'below_water'

    @property
    def MCF_decay(self):
        '''[float] Methane correction factor for COD degraded during storage.'''
        return float(self._MCF_decay[self._return_MCF_EF()])
    @MCF_decay.setter
    def MCF_decay(self, i):
        self._MCF_decay[self._return_MCF_EF()] = i

    @property
    def N2O_EF_decay(self):
        '''[float] Fraction of N emitted as N2O during storage.'''
        return float(self._N2O_EF_decay[self._return_MCF_EF()])
    @N2O_EF_decay.setter
    def N2O_EF_decay(self, i):
        self._N2O_EF_decay[self._return_MCF_EF()] = i


# %%

uddt_path = ospath.join(data_path, 'sanunit_data/_uddt.tsv')

class UDDT(Toilet):
    '''
    Urine-diverting dry toilet with liquid storage tank and dehydration vault
    for urine and feces storage, respectively, based on
    `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_,
    a subclass of :class:`qsdsan.sanunits.Toilet`.

    The following components should be included in system thermo object for simulation:
    Tissue, WoodAsh, H2O, NH3, NonNH3, P, K, Mg, CH4, N2O.

    The following impact items should be pre-constructed for life cycle assessment:
    Cement, Sand, Gravel, Brick, Plastic, Steel, StainlessSteelSheet, Wood.

    Parameters
    ----------
    ins : stream obj
        Excreta stream.
    outs : Iterable(stream obj)
        Recyclable liquid urine, recyclable solid feces, struvite scaling (irrecoverable),
        hydroxyapatite scaling (irrecoverable), fugitive CH4, and fugitive N2O.
    lifetime : int
        Lifetime of this UDDT, [yr].
    T : float
        Temperature, [K].
    safety_factor : float
        Safety factor for pathogen removal during onsite treatment,
        must be larger than 1.
    if_treatment : bool
        If has onsite treatment.

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
    :ref:`qsdsan.sanunits.Toilet <sanunits_toilets>`
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), N_user=1, N_toilet=1, N_tot_user=None,
                 if_toilet_paper=True, if_flushing=False, if_cleansing=False,
                 if_desiccant=True, if_air_emission=True, if_ideal_emptying=True,
                 CAPEX=553, OPEX_over_CAPEX=0.1, lifetime=8,
                 T=273.15+25, safety_factor=1, if_prep_loss=True, if_treatment=False,
                 **kwargs):

        Toilet.__init__(self, ID, ins, outs, thermo, init_with,
                        degraded_components, N_user, N_toilet, N_tot_user,
                        if_toilet_paper, if_flushing, if_cleansing, if_desiccant,
                        if_air_emission, if_ideal_emptying, CAPEX, OPEX_over_CAPEX)
        self.lifetime = lifetime
        self.T = T
        self._safety_factor = safety_factor
        self.if_prep_loss = if_prep_loss
        self.if_treatment = if_treatment

        self.construction = (
            Construction('cement', linked_unit=self, item='Cement', quantity_unit='kg'),
            Construction('sand', linked_unit=self, item='Sand', quantity_unit='kg'),
            Construction('gravel', linked_unit=self, item='Gravel', quantity_unit='kg'),
            Construction('brick', linked_unit=self, item='Brick', quantity_unit='kg'),
            Construction('liner', linked_unit=self, item='Plastic', quantity_unit='kg'),
            Construction('steel', linked_unit=self, item='Steel', quantity_unit='kg'),
            Construction('ss_sheet', linked_unit=self, item='StainlessSteelSheet', quantity_unit='kg'),
            Construction('wood', linked_unit=self, item='Wood', quantity_unit='m3'),
            )

        data = load_data(path=uddt_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, '_'+para, value)
        del data

        self._tank_V = 60/1e3 # m3
        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_outs = 6

    def _run(self):
        Toilet._run(self)
        liq, sol, struvite, HAP, CH4, N2O = self.outs
        liq.copy_like(self.ins[0])
        # Assume all of additives (toilet paper and desiccant) and water
        # (cleansing and flushing) go to feces
        _COD = self.ins[1]._COD or self.ins[1].COD
        sol_COD = _COD*self.ins[1].F_vol/1e3
        sol.mix_from(self.ins[1:])
        sol._COD = sol_COD*1e3/sol.F_vol
        struvite.phase = HAP.phase = 's'
        CH4.phase = N2O.phase = 'g'

        # Modified from ref [1], assume this only happens when air emission occurs
        # (to be consistent with pit latrine)
        if self.if_air_emission:
            # N loss due to ammonia volatilization
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(liq.TN/1e3*liq.F_vol*self.N_volatilization,
                                        liq.imass['NH3'])
            liq.imass ['NH3'] -= NH3_rmd
            liq.imass['NonNH3'] -= NonNH3_rmd
            # Energy/N loss due to degradation
            COD_loss = self.first_order_decay(k=self.decay_k_COD,
                                              t=self.collection_period/365,
                                              max_decay=self.COD_max_decay)
            CH4.imass['CH4'] = sol._COD/1e3*sol.F_vol*COD_loss * \
                self.max_CH4_emission*self.MCF_decay # COD in mg/L (g/m3)
            sol._COD *= 1 - COD_loss
            sol.imass[self.degraded_components] *= 1 - COD_loss

            N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.collection_period/365,
                                            max_decay=self.N_max_decay)
            N_loss_tot = sol.TN/1e3*sol.F_vol*N_loss
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(N_loss_tot, sol.imass['NH3'])
            sol.imass ['NH3'] -= NH3_rmd
            sol.imass['NonNH3'] -= NonNH3_rmd
            N2O.imass['N2O'] = N_loss_tot * self.N2O_EF_decay * 44/28
        else:
            CH4.empty()
            N2O.empty()

        # N and P losses due to struvite and hydroxyapatite (HAp)
        if self.if_prep_loss:
            # Struvite
            NH3_mol = liq.imol['NH3']
            P_mol = liq.imol['P']
            Mg_mol = liq.imol['Mg']
            Ksp = 10**(-self.struvite_pKsp)
            # Ksp = (initial N - struvite)(initial P - struvite)(initial Mg - struvite)
            coeff = [1, -(NH3_mol+P_mol+Mg_mol),
                     (NH3_mol*P_mol + P_mol*Mg_mol + Mg_mol*NH3_mol),
                     (Ksp - NH3_mol*P_mol*Mg_mol)]
            struvite_mol = 0
            for i in np.roots(coeff):
                if 0 < i < min(NH3_mol, P_mol, Mg_mol):
                    struvite_mol = i
            struvite.imol['Struvite'] = \
                max(0, min(NH3_mol, P_mol, Mg_mol, struvite_mol))
            liq.imol['NH3'] -= struvite_mol
            liq.imol['P'] -= struvite_mol
            liq.imol['Mg'] -= struvite_mol
            # HAP
            left_P = liq.imol['P'] - 3*(liq.imol['Ca']/5)
            # Remaining P enough to precipitate all Ca as HAP
            if left_P > 0:
                HAP.imol['HAP'] = liq.imol['Ca']/5
                liq.imol['P'] = left_P
                liq.imol['Ca'] = 0
            else:
                HAP.imol['HAP'] = liq.imol['P']/3
                liq.imol['Ca'] -= 5*(liq.imol['P']/3)
                liq.imol['P'] = 0
        else:
            struvite.empty()
            HAP.empty()

        # Onsite treatment
        if self.if_treatment:
            NH3_mmol = liq.imol['NH3'] * 1e3
            ur_DM = 1 - liq.imass['H2O']/liq.F_mass
            pKa = 0.09018 + (2729.92/self.T)
            f_NH3_Emerson = 1 / (10**(pKa-self.ur_pH)+1)
            alpha = 0.82 - 0.011*np.sqrt(NH3_mmol+1700*ur_DM)
            beta = 1.17 + 0.02 * np.sqrt(NH3_mmol+1100*ur_DM)
            f_NH3_Pitzer = f_NH3_Emerson * \
                (alpha + ((1-alpha)*(f_NH3_Emerson**beta)))
            NH3_conc = NH3_mmol * f_NH3_Pitzer

            # Time (in days) to reach desired inactivation level
            self.treatment_tau = ((3.2 + self.log_removal) \
                             / (10**(-3.7+0.062*(self.T-273.15)) * (NH3_conc**0.7))) \
                        * 1.14*self.safety_factor
            # Total volume in m3
            self.treatment_V = self.treatment_tau * liq.F_vol*24
        else:
            self.treatment_tau = self.treatment_V = 0

        # Feces water loss if desiccant is added
        if self.if_desiccant:
            sol_COD = sol._COD*sol.F_vol/1e3
            MC_min = self.fec_moi_min
            r = self.fec_moi_red_rate
            t = self.collection_period
            mixed = sol.copy()
            mixed.mix_from(self.ins[1:])
            fec_moi_int = mixed.imass['H2O']/mixed.F_mass
            fec_moi = MC_min + (fec_moi_int-MC_min)/(r*t)*(1-np.exp(-r*t))
            dry_mass = sol.F_mass-sol.imass['H2O']
            tot_mass = (sol.F_mass-sol.imass['H2O']) / (1-fec_moi)
            sol.imass['H2O'] = tot_mass - dry_mass
            sol._COD = sol_COD*1e3/sol.F_vol

        self._vault_V = sol.F_vol*self.collection_period*24 # in day

        # Non-ideal emptying
        if not self.if_ideal_emptying:
            liq, CH4, N2O = self.get_emptying_emission(
                waste=liq, CH4=CH4, N2O=N2O,
                empty_ratio=self.empty_ratio,
                CH4_factor=self.COD_max_decay*self.MCF_aq*self.max_CH4_emission,
                N2O_factor=self.N2O_EF_decay*44/28)
            sol, CH4, N2O = self.get_emptying_emission(
                waste=sol, CH4=CH4, N2O=N2O,
                empty_ratio=self.empty_ratio,
                CH4_factor=self.COD_max_decay*self.MCF_aq*self.max_CH4_emission,
                N2O_factor=self.N2O_EF_decay*44/28)

        # Scale up the effluent based on the number of user per toilet and
        # toilet number
        N_tot_user = self.N_tot_user or self.N_toilet*self.N_user
        for i in self.outs:
            if not i.F_mass == 0:
                i.F_mass *= N_tot_user


    _units = {
        'Collection period': 'd',
        'Single tank volume': 'm3',
        'Single vault volume': 'm3',
        'Treatment time': 'd',
        'Treatment volume': 'm3'
        }

    def _design(self):
        design = self.design_results
        design['Number of users per toilet'] = self.N_user
        design['Parallel toilets'] = N = self.N_toilet
        design['Collection period'] = self.collection_period
        design['Single tank volume'] = self.tank_V
        design['Single vault volume'] = self.vault_V
        design['Treatment time'] = self.treatment_tau
        design['Treatment volume'] = self.treatment_V

        density = self.density_dct
        constr = self.construction
        constr[0].quantity = 200 * N # cement
        constr[1].quantity = 0.6 * density['Sand'] * N
        constr[2].quantity = 0.2 * density['Gravel'] * N
        constr[3].quantity = 682*0.0024 * density['Brick'] * N
        constr[4].quantity = 4 * density['Plastic'] * N
        constr[5].quantity = 0.00351 * density['Steel'] * N
        constr[6].quantity = 28.05*density['StainlessSteelSheet'] * N
        constr[7].quantity = 0.222 * N # wood

        self.add_construction(add_cost=False)


    @property
    def safety_factor(self):
        return self._safety_factor
    @safety_factor.setter
    def safety_factor(self, i):
        if i < 1:
            raise ValueError(f'safety_factor must be larger than 1, not {i}.')
        self._safety_factor = i

    @property
    def collection_period(self):
        '''[float] Time interval between storage tank collection, [d].'''
        return self._collection_period
    @collection_period.setter
    def collection_period(self, i):
        self._collection_period = i

    @property
    def treatment_tau(self):
        '''[float] Time for onsite treatment (if treating), [d].'''
        return self._treatment_tau
    @treatment_tau.setter
    def treatment_tau(self, i):
        self._treatment_tau = i

    @property
    def treatment_V(self):
        '''[float] Volume needed to achieve treatment target (if treating), [d].'''
        return self._treatment_V
    @treatment_V.setter
    def treatment_V(self, i):
        self._treatment_V = i

    @property
    def tank_V(self):
        '''[float] Volume of the urine storage tank, [m3].'''
        return self._tank_V
    @tank_V.setter
    def tank_V(self, i):
        self._tank_V = i

    @property
    def vault_V(self):
        '''[float] Volume of the feces dehydration vault, [m3].'''
        return self._vault_V

    @property
    def struvite_pKsp(self):
        '''[float] Precipitation constant of struvite.'''
        return self._struvite_pKsp
    @struvite_pKsp.setter
    def struvite_pKsp(self, i):
        self._struvite_pKsp = i

    @property
    def prep_sludge(self):
        '''
        [float]
        Fraction of total precipitate appearing as sludge that can
        settle and be removed.
        '''
        return self._prep_sludge
    @prep_sludge.setter
    def prep_sludge(self, i):
        self._prep_sludge = i

    @property
    def log_removal(self):
        '''Desired level of pathogen inactivation.'''
        return self._log_removal
    @log_removal.setter
    def log_removal(self, i):
        self._log_removal = i

    @property
    def ur_pH(self):
        '''Urine pH.'''
        return self._ur_pH
    @ur_pH.setter
    def ur_pH(self, i):
        self._ur_pH = i

    @property
    def fec_moi_min(self):
        '''[float] Minimum moisture content of feces.'''
        return self._fec_moi_min
    @fec_moi_min.setter
    def fec_moi_min(self, i):
        self._fec_moi_min = i

    @property
    def fec_moi_red_rate(self):
        '''[float] Exponential reduction rate of feces moisture.'''
        return self._fec_moi_red_rate
    @fec_moi_red_rate.setter
    def fec_moi_red_rate(self, i):
        self._fec_moi_red_rate = i