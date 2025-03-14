#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''

This module is developed by:
    Yuyao Huang <yuyaoh20@gmail.com>
    Siqi Tang <siqit@outlook.com>
for Enviroloo (EL) Clear Toilet system

'''
# %%
import numpy as np
from math import ceil
from warnings import warn
from qsdsan import WasteStream
from qsdsan import SanUnit, Construction
from qsdsan.sanunits import IdealClarifier
from qsdsan.sanunits._tank import StorageTank
from qsdsan.processes._decay import Decay
from qsdsan.utils import ospath, load_data, data_path, price_ratio
# %% This callable file will be reposited to qsdsan.SanUnit subbranch with the name of _enviroloo
__all__ = (
    'EL_Excretion', # excretion
    'EL_Toilet', # toilet
    'EL_MURT', # toilet
    'EL_CT', # Collection tank
    'EL_PC', # Primary clarifier
    'EL_Anoxic', # Anoxic tank
    'EL_Aerobic', # Aerobic tank
    'EL_MBR', # Membrane filter
    'EL_CWT', # Clear water tank
    'EL_PT', # Pressure tank
    'EL_blower', # blower
    'EL_System', # System-level summary
    'EL_Housing', # Housing of EL_System, such as equipment's armor
    )

EL_su_data_path = ospath.join(data_path, 'sanunit_data/el')

# %%
excretion_path = ospath.join(EL_su_data_path, '_EL_excretion.tsv')

class EL_Excretion(SanUnit):
    '''
    Estimation of N, P, K, and COD in urine and feces based on dietary intake
    for one person based on `Trimmer et al. <https://doi.org/10.1021/acs.est.0c03296>`_

    Parameters
    ----------
    waste_ratio : float
        A ratio in [0, 1] to indicate the amount of intake calories and nutrients
        (N, P, K) that is wasted.

    Examples
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_

    References
    ----------
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
    Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
    Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
    https://doi.org/10.1021/acs.est.0c03296
    '''

    _N_ins = 0
    _N_outs = 2

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 waste_ratio=0, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.waste_ratio = waste_ratio

        data = load_data(path=excretion_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            # value = float(data.loc[para]['low'])
            # value = float(data.loc[para]['high'])
            setattr(self, '_'+para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _run(self):
        ur, fec = self.outs
        ur.empty()
        fec.empty()

        not_wasted = 1 - self.waste_ratio
        factor = 24 * 1e3 # from g per person per day to kg per hour

        ur_N = (self.p_veg+self.p_anim)/factor*self.N_prot \
           * self.N_exc*self.N_ur*not_wasted
        ur.imass['NH3'] = ur_N * self.N_ur_NH3
        ur.imass['NonNH3'] = ur_N - ur.imass['NH3']

        ur.imass['P'] = (self.p_veg*self.P_prot_v+self.p_anim*self.P_prot_a)/factor \
            * self.P_exc*self.P_ur*not_wasted

        e_cal = self.e_cal / 24 * not_wasted
        ur.imass['K'] = e_cal/1e3 * self.K_cal/1e3 * self.K_exc*self.K_ur
        ur.imass['Mg'] = self.Mg_ur / factor
        ur.imass['Ca'] = self.Ca_ur / factor

        ur_exc = self.ur_exc / factor
        ur.imass['H2O'] = self.ur_moi * ur_exc
        ur.imass['OtherSS'] = ur_exc - ur.F_mass

        fec_exc = self.fec_exc / factor
        fec_N = (1-self.N_ur)/self.N_ur * ur_N
        fec.imass['NH3'] = fec_N * self.N_fec_NH3
        fec.imass['NonNH3'] = fec_N - fec.imass['NH3']
        fec.imass['P'] = (1-self.P_ur)/self.P_ur * ur.imass['P']
        fec.imass['K'] = (1-self.K_ur)/self.K_ur * ur.imass['K']
        fec.imass['Mg'] = self.Mg_fec / factor
        fec.imass['Ca'] = self.Ca_fec / factor
        fec.imass['H2O'] = self.fec_moi * fec_exc
        fec.imass['OtherSS'] = fec_exc - fec.F_mass

        # 14 kJ/g COD, the average lower heating value of excreta
        tot_COD = e_cal*self.e_exc*4.184/14/1e3 # in kg COD/hr
        ur._COD = tot_COD*(1-self.e_fec) / (ur.F_vol/1e3) # in mg/L
        fec._COD = tot_COD*self.e_fec / (fec.F_vol/1e3) # in mg/L

    @property
    def e_cal(self):
        '''[float] Caloric intake, [kcal/cap/d].'''
        return self._e_cal
    @e_cal.setter
    def e_cal(self, i):
        self._e_cal = i

    @property
    def p_veg(self):
        '''[float] Vegetal protein intake, [g/cap/d].'''
        return self._p_veg
    @p_veg.setter
    def p_veg(self, i):
        self._p_veg = i

    @property
    def p_anim(self):
        '''[float] Animal protein intake, [g/cap/d].'''
        return self._p_anim
    @p_anim.setter
    def p_anim(self, i):
        self._p_anim = i

    @property
    def N_prot(self):
        '''[float] Nitrogen content in protein, [wt%].'''
        return self._N_prot
    @N_prot.setter
    def N_prot(self, i):
        self._N_prot = i

    @property
    def P_prot_v(self):
        '''[float] Phosphorus content in vegetal protein, [wt%].'''
        return self._P_prot_v
    @P_prot_v.setter
    def P_prot_v(self, i):
        self._P_prot_v = i

    @property
    def P_prot_a(self):
        '''[float] Phosphorus content in animal protein, [wt%].'''
        return self._P_prot_a
    @P_prot_a.setter
    def P_prot_a(self, i):
        self._P_prot_a = i

    @property
    def K_cal(self):
        '''[float] Potassium intake relative to caloric intake, [g K/1000 kcal].'''
        return self._K_cal
    @K_cal.setter
    def K_cal(self, i):
        self._K_cal = i

    @property
    def N_exc(self):
        '''[float] Nitrogen excretion factor, [% of intake].'''
        return self._N_exc
    @N_exc.setter
    def N_exc(self, i):
        self._N_exc = i

    @property
    def P_exc(self):
        '''[float] Phosphorus excretion factor, [% of intake].'''
        return self._P_exc
    @P_exc.setter
    def P_exc(self, i):
        self._P_exc = i

    @property
    def K_exc(self):
        '''[float] Potassium excretion factor, [% of intake].'''
        return self._K_exc
    @K_exc.setter
    def K_exc(self, i):
        self._K_exc = i

    @property
    def e_exc(self):
        '''[float] Energy excretion factor, [% of intake].'''
        return self._e_exc
    @e_exc.setter
    def e_exc(self, i):
        self._e_exc = i

    @property
    def N_ur(self):
        '''[float] Nitrogen recovered in urine, [wt%].'''
        return self._N_ur
    @N_ur.setter
    def N_ur(self, i):
        self._N_ur = i

    @property
    def P_ur(self):
        '''[float] Phosphorus recovered in urine, [wt%].'''
        return self._P_ur
    @P_ur.setter
    def P_ur(self, i):
        self._P_ur = i

    @property
    def K_ur(self):
        '''[float] Potassium recovered in urine, [wt%].'''
        return self._K_ur
    @K_ur.setter
    def K_ur(self, i):
        self._K_ur = i

    @property
    def e_fec(self):
        '''[float] Percent of excreted energy in feces, [%].'''
        return self._e_fec
    @e_fec.setter
    def e_fec(self, i):
        self._e_fec = i

    @property
    def N_ur_NH3(self):
        '''[float] Reduced inorganic nitrogen in urine, modeled as NH3, [% of total urine N].'''
        return self._N_ur_NH3
    @N_ur_NH3.setter
    def N_ur_NH3(self, i):
        self._N_ur_NH3 = i

    @property
    def N_fec_NH3(self):
        '''[float] Reduced inorganic nitrogen in feces, modeled as NH3, [% of total feces N].'''
        return self._N_fec_NH3
    @N_fec_NH3.setter
    def N_fec_NH3(self, i):
        self._N_fec_NH3 = i

    @property
    def ur_exc(self):
        '''[float] Urine generated per day, [g/cap/d].'''
        return self._ur_exc
    @ur_exc.setter
    def ur_exc(self, i):
        self._ur_exc = i

    @property
    def fec_exc(self):
        '''[float] Feces generated per day, [g/cap/d].'''
        return self._fec_exc
    @fec_exc.setter
    def fec_exc(self, i):
        self._fec_exc = i

    @property
    def ur_moi(self):
        '''[float] Moisture (water) content of urine, [wt%].'''
        return self._ur_moi
    @ur_moi.setter
    def ur_moi(self, i):
        self._ur_moi = i

    @property
    def fec_moi(self):
        '''[float] Moisture (water) content of feces, [wt%].'''
        return self._fec_moi
    @fec_moi.setter
    def fec_moi(self, i):
        self._fec_moi = i

    @property
    def Mg_ur(self):
        '''[float] Magnesium excreted in urine, [g Mg/cap/d].'''
        return self._Mg_ur
    @Mg_ur.setter
    def Mg_ur(self, i):
        self._Mg_ur = i

    @property
    def Mg_fec(self):
        '''[float] Magnesium excreted in feces, [g Mg/cap/d].'''
        return self._Mg_fec
    @Mg_fec.setter
    def Mg_fec(self, i):
        self._Mg_fec = i

    @property
    def Ca_ur(self):
        '''[float] Calcium excreted in urine, [g Ca/cap/d].'''
        return self._Ca_ur
    @Ca_ur.setter
    def Ca_ur(self, i):
        self._Ca_ur = i

    @property
    def Ca_fec(self):
        '''[float] Calcium excreted in feces, [g Ca/cap/d].'''
        return self._Ca_fec
    @Ca_fec.setter
    def Ca_fec(self, i):
        self._Ca_fec = i

    @property
    def waste_ratio(self):
        '''
        [float] The amount of intake calories and nutrients
        (N, P, K) that is wasted.

        .. note::
            Not considered for Mg and Ca.
        '''
        return self._waste_ratio
    @waste_ratio.setter
    def waste_ratio(self, i):
        self._waste_ratio = i

# %%
toilet_path = ospath.join(EL_su_data_path, '_EL_toilet.tsv')

class EL_Toilet(SanUnit, Decay, isabstract=True):
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
    price_ratio : float
        Calculated capital cost will be multiplied by this number
        to consider the effect in cost difference from different locations.

    References
    ----------
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
    Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
    Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
    https://doi.org/10.1021/acs.est.0c03296.

    See Also
    --------
    :ref:`qsdsan.processes.Decay <processes_Decay>`

    '''
    _N_ins = 6
    _outs_size_is_fixed = False
    density_dct = {
        'Sand': 1442,
        'Gravel': 1600,
        'Brick': 1750,
        'Plastic': 0.63,
        'Steel': 7900,
        'StainlessSteelSheet': 2.64
        }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), N_user=1, N_toilet=1, N_tot_user=None,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                 CAPEX=None, OPEX_over_CAPEX=None, price_ratio=1.):

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
        self.price_ratio = price_ratio

        data = load_data(path=toilet_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            if para in ('desiccant_V', 'desiccant_rho'):
                setattr(self, para, value)
            else:
                setattr(self, '_'+para, value)
        del data

        # self._empty_ratio = 0.59 default
        self._empty_ratio = 0


    def _run(self):
        ur, fec, tp, fw, cw, des = self.ins
        tp.imass['Tissue'] = int(self.if_toilet_paper)*self.toilet_paper
        fw.imass['H2O'] = int(self.if_flushing)*self.flushing_water
        cw.imass['H2O'] = int(self.if_cleansing)*self.cleansing_water
        des.imass['WoodAsh'] = int(self.if_desiccant)*self.desiccant

    def _scale_up_outs(self):
        '''
        Scale up the effluent based on the number of user per toilet and
        toilet number.
        '''
        N_tot_user = self.N_tot_user or self.N_toilet*self.N_user
        for i in self.outs:
            if not i.F_mass == 0:
                i.F_mass *= N_tot_user


    def _cost(self):
        self.baseline_purchase_costs['Total toilets'] = self.CAPEX * self.N_toilet * self.price_ratio
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
        COD_rmvd = waste.COD*(1-empty_ratio)/1e3*waste.F_vol
        CH4.imass['CH4'] += COD_rmvd * CH4_factor
        N2O.imass['N2O'] += COD_rmvd * N2O_factor
        waste.mass *= empty_ratio

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

    @property
    def if_N2O_emission(self):
        '''[bool] Whether to consider N degradation and fugitive N2O emission.'''
        return self.if_air_emission
    @if_N2O_emission.setter
    def if_N2O_emission(self, i):
        raise ValueError('Setting `if_N2O_emission` for `PitLatrine` is not supported, '
                         'please set `if_air_emission` instead.')

# %%
murt_path = ospath.join(EL_su_data_path, '_EL_murt.tsv')

@price_ratio()
class EL_MURT(EL_Toilet):
    '''
    Multi-unit reinvented toilet.

    The following components should be included in system thermo object for simulation:
    Tissue, WoodAsh, H2O, NH3, NonNH3, P, K, Mg, CH4, N2O.

    The following impact items should be pre-constructed for life cycle assessment:
    Ceramic, Fan.

    Parameters
    ----------
    ins : Iterable(stream)
        waste_in: mixed excreta.
    Outs : Iterable(stream)
        waste_out: degraded mixed excreta.
        CH4: fugitive CH4.
        N2O: fugitive N2O.
    N_squatting_pan_per_toilet : int
        The number of squatting pan per toilet.
    N_urinal_per_toilet : int
        The number of urinals per toilet.
    if_include_front_end : bool
        If False, will not consider the capital and operating costs of this unit.

    See Also
    --------
    :ref:`qsdsan.sanunits.Toilet <sanunits_toilet>`
    '''
    _N_outs = 3
    _units = {
        'Collection period': 'd',
        }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), N_user=1, N_tot_user=1,
                 N_toilet=None, if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=True, if_air_emission=True, if_ideal_emptying=True,
                 CAPEX=0, OPEX_over_CAPEX=0, lifetime=10,
                 N_squatting_pan_per_toilet=1, N_urinal_per_toilet=1,
                 if_include_front_end=True, **kwargs):

        EL_Toilet.__init__(
            self, ID, ins, outs, thermo=thermo, init_with=init_with,
            degraded_components=degraded_components,
            N_user=N_user, N_tot_user=N_tot_user, N_toilet=N_toilet,
            if_toilet_paper=if_toilet_paper, if_flushing=if_flushing,
            if_cleansing=if_cleansing, if_desiccant=if_desiccant,
            if_air_emission=if_air_emission, if_ideal_emptying=if_ideal_emptying,
            CAPEX=CAPEX, OPEX_over_CAPEX=OPEX_over_CAPEX
            )
        self.N_squatting_pan_per_toilet = N_squatting_pan_per_toilet
        self.N_urinal_per_toilet = N_urinal_per_toilet
        self.if_include_front_end = if_include_front_end
        self._mixed_in = WasteStream(f'{self.ID}_mixed_in')

        data = load_data(path=murt_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [
            Construction(item='Ceramic', linked_unit=self, quantity_unit='kg'),
            Construction(item='Fan', linked_unit=self, quantity_unit='ea'),
        ]

    def _run(self):
        EL_Toilet._run(self)
        mixed_out, CH4, N2O = self.outs
        CH4.phase = N2O.phase = 'g'

        mixed_in = self._mixed_in
        mixed_in.mix_from(self.ins)
        #tot_COD_kg = sum(float(getattr(i, 'COD', 0)) * i.F_vol for i in self.ins) / 1e3
        tot_COD_kg = sum(float(getattr(i, 'COD')) * i.F_vol for i in self.ins) / 1e3
        
        # Air emission
        if self.if_air_emission:
            # N loss due to ammonia volatilization
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(mixed_in.TN/1e3*mixed_in.F_vol*self.N_volatilization,
                                        mixed_in.imass['NH3'])
            mixed_in.imass ['NH3'] -= NH3_rmd
            mixed_in.imass['NonNH3'] -= NonNH3_rmd
            
            # Energy/N loss due to degradation
            mixed_in._COD = tot_COD_kg * 1e3 / mixed_in.F_vol # accounting for COD loss in leachate
            Decay._first_order_run(self, waste=mixed_in, treated=mixed_out, CH4=CH4, N2O=N2O)
        else:
            mixed_out.copy_like(mixed_in)
            CH4.empty()
            N2O.empty()
            
        # Aquatic emission when not ideally emptied
        if not self.if_ideal_emptying:
           self.get_emptying_emission(
                waste=mixed_out, CH4=CH4, N2O=N2O,
                empty_ratio=self.empty_ratio,
                CH4_factor=self.COD_max_decay*self.MCF_aq*self.max_CH4_emission,
                N2O_factor=self.N2O_EF_decay*44/28)
        
        self._scale_up_outs()

    def _design(self):
        design = self.design_results
        constr = self.construction
        if self.if_include_front_end:
            design['Number of users per toilet'] = self.N_user
            design['Parallel toilets'] = N = self.N_toilet
            design['Collection period'] = self.collection_period
            design['Ceramic'] = Ceramic_quant = (
                self.squatting_pan_weight * self.N_squatting_pan_per_toilet+
                self.urinal_weight * self.N_urinal_per_toilet
                )
            design['Fan'] = Fan_quant = 1  # assume fan quantity is 1
            constr[0].quantity = Ceramic_quant * N
            constr[1].quantity = Fan_quant * N
            self.add_construction(add_cost=False)
        else:
            design.clear()
            for i in constr: i.quantity = 0

    def _cost(self):
        C = self.baseline_purchase_costs
        if self.if_include_front_end:
            N_toilet = self.N_toilet
            C['Ceramic Toilets'] = (
                self.squatting_pan_cost * self.N_squatting_pan_per_toilet +
                self.urinal_cost * self.N_urinal_per_toilet
                ) * N_toilet
            C['Fan'] = self.fan_cost * N_toilet
            C['Misc. parts'] = (
                self.led_cost +
                self.anticor_floor_cost +
                self.circuit_change_cost +
                self.pipe_cost
                ) * N_toilet

            ratio = self.price_ratio
            for equipment, cost in C.items():
                C[equipment] = cost * ratio
        else:
            self.baseline_purchase_costs.clear()

        sum_purchase_costs = sum(v for v in C.values())
        self.add_OPEX = (
            self._calc_replacement_cost() +
            self._calc_maintenance_labor_cost() +
            sum_purchase_costs * self.OPEX_over_CAPEX / (365 * 24)
            )

    def _calc_replacement_cost(self):
        return 0

    def _calc_maintenance_labor_cost(self):
        return 0

    @property
    def collection_period(self):
        '''[float] Time interval between storage tank collection, [d].'''
        return self._collection_period
    @collection_period.setter
    def collection_period(self, i):
        self._collection_period = float(i)
        
    @property
    def tau(self):
        '''[float] Retention time of the unit, same as `collection_period`.'''
        return self.collection_period
    @tau.setter
    def tau(self, i):
        self.collection_period = i

# %%
CollectionTank_path = ospath.join(EL_su_data_path, '_EL_CT.tsv')

@price_ratio()
class EL_CT(StorageTank):
    
    '''
    Name
    ----
    Collection tank in the Enviroloo (EL) Clear Toilet system.
    
    Parameters
    ----------
    Ins: 
    (1) Mixed wastewater
    (2) Primary clarifier effluent spill
    (3) Clear water tank effluent spill 
    (4) Primary clarifier sludge return

    Outs:
    (1) Treated water
    (2) Methane (CH4)
    (3) Nitrous oxide (N2O)

    Attributes
    ----------
    length_to_diameter : float
        Ratio of the tank length to diameter.
    vessel_material : str
        Material used for constructing the vessel.
    sludge_moisture_content : float
        Moisture content of the sludge (mass of water/total mass).
    COD_removal : float
        Fraction of COD removed in the collection tank.

    References
    ----------
    Similar to the :class:`biosteam.units.MixTank`, but can calculate material usage.

    See Also
    ----------
    class:`biosteam.units.StorageTank`
    '''
    
    _N_ins = 4 
    _N_outs = 1
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    exponent_scale = 0.1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, ppl=None, baseline_ppl=None,
                 vessel_type= 'Field erected', tau=24, V_wf=None, vessel_material='Stainless steel', kW_per_m3=0.1,
                 init_with='WasteStream', F_BM_default=1,
                 include_construction=True, **kwargs):
        StorageTank.__init__(self, 
                      # Basic parameters
                      ID=ID, # The unique identifier of the tank
                      ins=ins, # The input stream to the tank
                      outs=outs, # The output streams from the tank
                      thermo=thermo, # The thermodynamic property package for simulating physical or chemical processes
                      # Other control parameters
                      init_with=init_with, # The method to initialize the tank contents.
                      F_BM_default=F_BM_default, # The default bare module factor for the tank cost estimation
                      include_construction=include_construction, # A Boolean value indicating whether the tank's construction material is considered in the cost analysis or life cycle assessment.
                      # Design parameters
                      vessel_type=vessel_type, # The type of the tank
                      tau=tau, # The retention time of the tank contents, the important parameters involved in mixing, reaction or separation processes.
                      V_wf=V_wf, # The volume working fraction of the tank, the ratio of the volume of the tank contents to the total volume of the tank.
                      vessel_material=vessel_material, # The material of the tank
                      kW_per_m3=kW_per_m3, # The power consumption per unit volume of the tank
                      )
        self.ppl = ppl
        self.baseline_ppl = baseline_ppl
        self._mixed_in = WasteStream(f'{self.ID}_mixed_in')

        data = load_data(path=CollectionTank_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]

    def _run(self):
        # Input stream
        WasteWater = self.ins[0]
        sludge_return = self.ins[1]
        PC_spill_return = self.ins[2]
        CWT_spill_return = self.ins[3]

        mixed_in = self._mixed_in
        # Define input streams
        input_streams = [WasteWater, sludge_return, PC_spill_return, CWT_spill_return]
        mixed_in.mix_from(input_streams)

        # Output stream
        TreatedWater = self.outs[0]
        
        # Copy the mixed result to the outflow
        TreatedWater.copy_like(mixed_in)

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tank'] = self.collection_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
    
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.power_demand_CT
        self.power_utility(power_demand)
    
    def _calc_replacement_cost(self):
        scale  = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        CT_replacement_cost = (
            self.collection_tank_cost / self.collection_tank_lifetime +               
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime) * scale
        CT_replacement_cost = CT_replacement_cost / (365 * 24)  # convert to USD/hr
        return CT_replacement_cost

# %%
PrimaryClarifier_path = ospath.join(EL_su_data_path, '_EL_PC.tsv')

@price_ratio()
class EL_PC(SanUnit):
    """
    Primary clarifier in the Enviroloo (EL) Clear Toilet system for COD and suspended solids removal.

    Parameters
    ----------
    ID : str
        Unique identifier for the unit.
    ins : tuple
        Input streams: (0) Wastewater from lift pump, (1) Nitrate return from membrane tank.
    outs : tuple
        Output streams: (0) Treated water, (1) Spill return, (2) Sludge return.
    sludge_flow_rate : float
        Sludge flow rate (m³/d). Required for sludge return calculation.
    solids_removal_efficiency : float
        Fraction of suspended solids removed (0 to 1). Default is 0.5.
    max_overflow : float
        Maximum allowable overflow rate (m³/h). Default is 15.
    ppl : float
        Current population served.
    baseline_ppl : float
        Baseline population for scaling design and cost.

    Notes
    -----
    - COD and suspended solids are tracked via the `components` object.
    - Spill return occurs if treated water flow exceeds `max_overflow`.
    - Inherits from `SanUnit`, not `IdealClarifier`, for flexibility.
    """

    _N_ins = 2
    _N_outs = 3
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    exponent_scale = 0.1

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 sludge_flow_rate=None, solids_removal_efficiency=0.5, max_overflow=15,
                 F_BM=1, transportation=[],
                 ppl=None, baseline_ppl=None, **kwargs):
        SanUnit.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo, init_with=init_with, 
                         F_BM_default=F_BM, transportation=transportation)
        self.sludge_flow_rate = sludge_flow_rate  # m³/d
        self.solids_removal_efficiency = solids_removal_efficiency  # 0 to 1
        self.max_overflow = max_overflow  # m³/h
        self.ppl = ppl
        self.baseline_ppl = baseline_ppl

        # Load default parameters from data file
        data = load_data(path=PrimaryClarifier_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items(): 
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]

    def _run(self):
        """Simulate the primary clarifier operation."""
        # Inputs
        wastewater, nitrate_return = self.ins
        # Outputs
        treated_water, sludge_return, spill_return = self.outs

        # Mix input streams
        mixed = wastewater.copy()  # Temporary stream for mixing
        mixed.mix_from([wastewater, nitrate_return])

        # Get total inflow (m³/d) and TSS (mg/L)
        Q_in = mixed.F_vol * 24  # Convert m³/h to m³/d
        TSS_in = sludge_return.F_mass # assume all return flow only including 
                                      # sewage sludge that was treated as TSS
        TSS_in = mixed.get_TSS()  # Total suspended solids in mg/L

        # Handle case with no solids
        if TSS_in <= 0:
            treated_water.copy_like(mixed)
            sludge_return.empty()
            spill_return.empty()
            return

        # Validate parameters
        if not self.sludge_flow_rate or not self.solids_removal_efficiency:
            raise ValueError("Must specify 'sludge_flow_rate' and 'solids_removal_efficiency'.")

        # Calculate sludge split
        Qs = self.sludge_flow_rate  # m³/d
        e_rmv = self.solids_removal_efficiency
        f_Qu = Qs / Q_in  # Fraction of flow to sludge
        f_Xu = e_rmv + (1 - e_rmv) * f_Qu  # Fraction of solids to sludge

        # Component-wise split (simplified for COD and SS)
        cmps = self.components
        split_to_sludge = np.zeros(len(cmps))
        for i, cmp in enumerate(cmps):
            if cmp.ID == 'OtherSS':  # Suspended solids
                split_to_sludge[i] = f_Xu
            elif cmp.ID == 'S_COD':  # Soluble COD, assume some removal with solids
                split_to_sludge[i] = f_Qu * 0.1  # Arbitrary small fraction
            else:
                split_to_sludge[i] = f_Qu  # Other components follow flow split
        split_to_sludge = np.clip(split_to_sludge, 0, 1)

        # Split mixed stream
        mixed.split_to(sludge_return, treated_water, split_to_sludge)

        # Handle spill return
        # if self.max_overflow is not None:
        #     if treated_water.F_vol > self.max_overflow:
        #         # Spill return exists
        #         spill_vol = treated_water.F_vol - self.max_overflow  
        #         if not hasattr(self, '_f_spill'):
        #             self._f_spill = None
        #         if not hasattr(self, '_f_overflow'):
        #             self._f_overflow = None
        #             self._f_spill = spill_vol / treated_water.F_vol  
        #             self._f_overflow = 1 - self._f_spill  

        #             spill_return.copy_like(treated_water)  
        #             spill_return.F_mass *= self._f_spill  

        #             TreatedWater.F_mass *= self._f_overflow  
        #     else:
        #         # max_overflow is not none, but TreatedWater < max_overflow
        #         spill_return.empty()
        #         if hasattr(self, '_f_spill'):
        #             del self._f_spill  
        #         if hasattr(self, '_f_overflow'):
        #             del self._f_overflow  
        # else:
        #     # max_overflow is none, no spill return
        #     spill_return.empty()
        
        max_overflow_m3d = self.max_overflow * 24  # Convert m³/h to m³/d
        Q_treated = treated_water.F_vol * 24  # m³/d
        if Q_treated > max_overflow_m3d:
            spill_vol = Q_treated - max_overflow_m3d
            f_spill = spill_vol / Q_treated
            spill_return.copy_like(treated_water)
            spill_return.F_mass *= f_spill
            treated_water.F_mass *= (1 - f_spill)
        else:
            spill_return.empty()

    def _design(self):
        """Calculate design parameters."""
        self.design_results['StainlessSteel'] = (
            self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
        )
        self.construction = [
            Construction(item='StainlessSteel', quantity=self.design_results['StainlessSteel'], quantity_unit='kg')
        ]
        self.add_construction(add_cost=False)

    def _cost(self):
        """Calculate capital and operating costs."""
        C = self.baseline_purchase_costs
        C['Tank'] = self.PC_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings

        ratio = self.price_ratio  # Assume price_ratio decorator sets this
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost()
        power_demand = self.power_demand_PC  # Default to 0 if not set
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        """Calculate replacement cost in USD/hr."""
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        replacement_cost = (
            self.PC_tank_cost / self.PC_tank_lifetime +
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime
        ) * scale
        return replacement_cost / (365 * 24)  # Convert to USD/hr



# class EL_PC(IdealClarifier):
    
#     """
#     Name
#     ----
#     Primary clarifier in the Enviroloo (EL) Clear Toilet system.

#     Introduction
#     ------------
#     The primary treatment of the EL system uses anaerobic digestion
#     to treat wastes (similar to a septic tank).

#     It can be used in conjunction with a membrane bioreactor (MBR)
#     to recovery N and P as struvite.

#     The following impact items should be pre-constructed for life cycle assessment:
#     FRP.

#     Parameters
#     ----------
#     Ins:
#     (1) influent of treated wastetwater from lift pump
#     (2) nitrate return flow from membrane tank

#     Outs:
#     (1) effluent of treated wastewater
#     (2) sludge return flow to collection tank
#     (3) spill flow to collection tank
#     (4) fugitive CH4 emission
#     (5) fugitive N2O emission

#     Attributes
#     ----------

    
#     References
#     ----------
#      refer to the qsdsan.sanunits.PrimaryClarifier module

#      """
    
#     _N_ins = 2
#     _N_outs = 3
#     _ins_size_is_fixed = True
#     _outs_size_is_fixed = True
#     exponent_scale = 0.6
    
#     def __init__(self, ID='', ins=None, outs=(), thermo=None, isdynamic=False, max_overflow=15,
#                  ppl = None, baseline_ppl = None, sludge_flow_rate=None,
#                  solids_removal_efficiency=None,
#                  F_BM_default=1, init_with='WasteStream', **kwargs):
#         """
#         Initialize the primary clarifier with default parameters:
#         - sludge_flow_rate: Default to  m3/d
#         - solids_removal_efficiency: Default to 50% (0.5)
#         """
#         IdealClarifier.__init__(self, ID=ID, ins=ins, outs=outs, thermo=thermo, sludge_flow_rate=sludge_flow_rate,
#                                   solids_removal_efficiency=solids_removal_efficiency,
#                                   isdynamic=isdynamic, init_with=init_with, F_BM_default=F_BM_default
#                                   )
        
#         self.ppl = ppl
#         self.baseline_ppl = baseline_ppl

#         self.max_overflow = max_overflow  # m^3/hr
#         # self.if_with_MBR = if_with_MBR

#         data = load_data(path=PrimaryClarifier_path)
#         for para in data.index:
#             value = float(data.loc[para]['expected'])
#             setattr(self, para, value)
#         del data

#         for attr, value in kwargs.items(): 
#             setattr(self, attr, value)

#     def _init_lca(self):
#         self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]
    
#     def _run(self):
        
#           # Input stream
#           WasteWater, MT_nitrate_return = self.ins
#           MT_nitrate_return.mass = WasteWater.mass * 0.1
#           #MT_nitrate_return = self.ins[1]  # Nitrate from membrane tank over return pump
        
#           # Output stream
#           TreatedWater = self.outs[0]
#           PC_spill_return = self.outs[1]  # Spill water to collection tank
#           PC_sludge_return = self.outs[2]  # Sludge to collection tank over return pump
          
#           # Inherited input stream properties
#           #TreatedWater.copy_like(WasteWater)
#           self._mixed.mix_from([WasteWater, MT_nitrate_return])
#           TreatedWater.copy_like(self._mixed)
#           #PC_sludge_return.F_mass = TreatedWater.F_mass * 0.1

#           # Sludge with water removal
#           PC_sludge_return.empty()
#           PC_sludge_return.mass = TreatedWater.mass * 0.1  # sludge return ratio, temeporary assumption
#           PC_sludge_return.copy_flow(TreatedWater, ('Mg', 'Ca', 'OtherSS', 'Tissue', 'WoodAsh'), remove=True)
#           PC_sludge_return.imass['OtherSS'] = WasteWater.F_mass * self.solids_removal_efficiency
#           PC_sludge_return.imass['H2O'] = PC_sludge_return.imass['OtherSS']/(1-self.solids_moisture_content) - PC_sludge_return.imass['OtherSS']
          
#           # Spill water return
#           if self.max_overflow is not None:
#               if TreatedWater.F_vol > self.max_overflow:
                      
#                   # Spill return exists
#                   spill_vol = TreatedWater.F_vol - self.max_overflow  

#                   if not hasattr(self, '_f_spill'):
#                       self._f_spill = None
#                   if not hasattr(self, '_f_overflow'):
#                       self._f_overflow = None

#                   self._f_spill = spill_vol / TreatedWater.F_vol  
#                   self._f_overflow = 1 - self._f_spill  

#                   PC_spill_return.copy_like(TreatedWater)  
#                   PC_spill_return.F_mass *= self._f_spill  

#                   TreatedWater.F_mass *= self._f_overflow  
#               else:
#                   # max_overflow is not none, but TreatedWater < max_overflow
#                   PC_spill_return.empty()
#                   if hasattr(self, '_f_spill'):
#                       del self._f_spill  
#                   if hasattr(self, '_f_overflow'):
#                       del self._f_overflow  
#           else:
#               # max_overflow is none, no spill return
#               PC_spill_return.empty()
          
#     def _design(self):
#         design = self.design_results
#         constr = self.construction
#         design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
#         self.add_construction(add_cost=False)

#     def _cost(self):  
#         C = self.baseline_purchase_costs
#         C['Tank'] = self.PC_tank_cost
#         C['Pipes'] = self.pipeline_connectors
#         C['Fittings'] = self.weld_female_adapter_fittings
    
#         ratio = self.price_ratio 
#         for equipment, cost in C.items():
#             C[equipment] = cost * ratio
        
#         self.add_OPEX = self._calc_replacement_cost()
        
#         power_demand = self.power_demand_PC
#         self.power_utility(power_demand)
    
#     def _calc_replacement_cost(self):
#         scale  = (self.ppl / self.baseline_ppl) ** self.exponent_scale
#         PC_replacement_cost = (
#             self.PC_tank_cost / self.PC_tank_lifetime +               
#             self.pipeline_connectors / self.pipeline_connectors_lifetime +
#             self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime) * scale
#         PC_replacement_cost = PC_replacement_cost / (365 * 24)  # convert to USD/hr
#         return PC_replacement_cost

# %%
Anoxic_path = ospath.join(EL_su_data_path, '_EL_Anoxic.tsv')

@price_ratio()
class EL_Anoxic(SanUnit, Decay):
    '''
    Anoxic treatment unit in the EL system.

    Parameters
    ----------
    Ins:
    (1) influent of treated wastetwater from primary clarifier
    (2) nitrate return flow from membrane tank
    (3) glucose addition
    (4) agitation pump

    Outs:
    (1) effluent of treated wastewater
    (2) fugitive CH4 emission
    (3) fugtivie N2O emission

    Attributes
    ----------

    
    References
    ----------
     refer to the qsdsan.sanunits.PrimaryClarifier module
    
    
    '''
    _N_ins = 4
    _N_outs = 3
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    baseline_ppl = 30
    exponent_scale = 0.1

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', methane_density=0.657,
                 degraded_components=('OtherSS',),  F_BM_default=1, ppl = None, baseline_ppl = None,
                 glucose_density=1560, if_capture_biogas=False, if_N2O_emission=True, **kwargs):
        Decay.__init__(self, ID, ins, outs, thermo=thermo,
                       init_with=init_with, F_BM_default=F_BM_default,
                       degraded_components=degraded_components,
                       if_capture_biogas=if_capture_biogas,
                       if_N2O_emission=if_N2O_emission,)
        
        self.ppl = ppl
        self.baseline_ppl = baseline_ppl
        self.methane_density = methane_density  # kg/m^3
        self.glucose_density = glucose_density  # kg/m^3
        self._mixed_in = WasteStream(f'{self.ID}_mixed_in')

        data = load_data(path=Anoxic_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]      
        
    def _run(self):

          # Input stream
          WasteWater = self.ins[0]
          nitrate_return = self.ins[1]  # Nitrate from membrane tank over return pump
          
          glucose = self.ins[2]  # Extra carbon source
          glucose_consumed = WasteWater.F_mass * self.chemical_glucose_dosage 
          #TODO: this dosage is the same as PAC. Check data source
          glucose.imass['Glucose'] = glucose_consumed #kg/h
          # glucose.mass = glucose_consumed
          
          agitation = self.ins[3]  # Agitation pump works here, but it does not attend mass flow balance calculation
          agitation.empty()
          
          mixed_in = self._mixed_in
          # Define input streams
          input_streams = [WasteWater, nitrate_return, glucose, agitation]
          mixed_in.mix_from(input_streams)
          
          # Output stream
          TreatedWater = self.outs[0]
          CH4_emission = self.outs[1]
          N2O_emission = self.outs[2]

          # Output stream
          # TreatedWater = self.outs[0]
          # Copy the mixed result to the outflow
          TreatedWater.copy_like(mixed_in)
          # Inherited input stream
          # TreatedWater.copy_like(WasteWater)
          
          # Manually add return nitrate (NO3)
          TreatedWater.imass['NO3'] += nitrate_return.imass['NO3']
          
          # Glucose consumption
          # glucose_consumed = WasteWater.F_vol * self.chemical_glucose_dosage * self.glucose_density
          # glucose.imass['Glucose'] = glucose_consumed
          TreatedWater.imass['Glucose'] += glucose.imass['Glucose']
          
          # glucose.mass = glucose_consumed

          # COD removal
          COD_removal = self.EL_anoT_COD_removal
          removed_COD = (WasteWater.COD + glucose.COD) / 1e3 * WasteWater.F_vol * COD_removal  # kg/hr
          
          # Now we explicitly remove Glucose from TreatedWater
          TreatedWater.imass['Glucose'] = 0  # All glucose is consumed in reaction
          glucose.empty()
          
          # Sludge produced
          sludge_prcd = self.EL_anoT_sludge_yield * removed_COD  # produced biomass
          
          for component in ('Mg', 'Ca', 'OtherSS', 'Tissue', 'WoodAsh'):
              TreatedWater.imass[component] += sludge_prcd  # New sludge produced
          
          # CH4 produced
          CH4_prcd = removed_COD * self.EL_anoT_methane_yield * self.methane_density  # kg CH4 produced/hr
          CH4_soluble = self.EL_anoT_soluble_methane_fraction * CH4_prcd
          CH4_emitted = CH4_prcd - CH4_soluble
          CH4_emission.imass['SolubleCH4'] = CH4_emitted
          TreatedWater.imass['SolubleCH4'] = CH4_soluble
          
          # N2O produced
          N_removal = self.EL_anoT_TN_removal
          if self.if_N2O_emission:
          # Assume the removal part covers both degradation loss
          # and other unspecified removal as well
              N_loss = self.first_order_decay(k = self.decay_k_N,
                                    t = self.EL_anoT_tau / 365,
                                    max_decay = self.N_max_decay)
              if N_loss > N_removal:
                  warn(f'Nitrogen degradation loss ({N_loss:.2f}) larger than the given removal ({N_removal:.2f})), '
                        'the degradation loss is assumed to be the same as the removal.')
                  N_loss = N_removal
            
              # N2O only from the degraded part
              N_loss_tot = N_loss * WasteWater.TN / 1e3 * WasteWater.F_vol
              N2O_emission.imass['N2O'] = N_loss_tot * self.N2O_EF_decay * 44 / 28
          else:
              N2O_emission.empty()
              N2O_emission.imass['N2O'] = 0  # All N2O is emitted
              
          # Assume all NO3 is consumed and does not appear in TreatedWater
          TreatedWater.imass['NO3']  = 0  
          
          # NH3 & NonNH3, P, K calculation
          total_solubles = np.array([
              CH4_soluble,
              WasteWater.imass['NH3'] * (1 - N_loss),  # In water or sludge
              WasteWater.imass['NonNH3'] * (1 - N_loss),  # In water or sludge
              0,  # NO3 tolly consumed
              WasteWater.imass['P'],  # In water or sludge
              WasteWater.imass['K'],  # In water or sludge
              ])
          
          # Removed solubles in the sludge, assume minimal used for growth
          sludge_solubles = total_solubles * np.array([
                            0,  # CH₄ will not go into sludge
                            N_removal,  # Part of it in water and in sludge
                            N_removal,  # Part of it in water and in sludge
                            0,  # NO3 tolly consumed
                            self.EL_anoT_TP_removal,  # Part of it in water and in sludge
                            0,  # K will not go into sludge
                            ])
          
          # Calculate solutes entering TreatedWater
          liquid_solubles = total_solubles - sludge_solubles
          
          # Assign to outputs (TreatedWater now includes sludge)
          TreatedWater.imass['SolubleCH4', 'NH3', 'NonNH3', 'NO3', 'P', 'K'] = liquid_solubles
          TreatedWater.imass['SolubleCH4', 'NH3', 'NonNH3', 'NO3', 'P', 'K'] += sludge_solubles  # Here is only one output stream (TreatedWater + sludge)
          
          # Final COD
          TreatedWater._COD = WasteWater.COD * (1 - self.EL_anoT_COD_removal)
          
    # def _run(self):
        
    #     # Input stream
    #     WasteWater = self.ins[0]
    #     sludge_return = self.ins[1]
    #     glucose = self.ins[2]
    #     agiation = self.ins[3]
        
    #     # Output stream
    #     TreatedWater = self.outs[0]
    #     CH4_emission = self.outs[1]
    #     N2O_emission = self.outs[2]
        
    #     input_streams = [WasteWater, sludge_return, glucose]
            
    #     # Mix all inputs into a single stream
    #     self._mixed.empty()
    #     self._mixed.mix_from(input_streams)

    #     # Copy the mixed result to the outflow
    #     TreatedWater.copy_like(self._mixed)
    
   
    # def _run(self):
        
    #     # Input stream
    #     WasteWater = self.ins[0]
    #     sludge_return = self.ins[1]
    #     glucose = self.ins[2]

    #     # Output stream
    #     TreatedWater = self.outs[0]
    #     CH4_emission = self.outs[1]
    #     N2O_emission = self.outs[2]

    #     # Combined sludge return and influent
    #     WasteWater.F_mass += sludge_return.F_mass
    #     WasteWater.imass += sludge_return.imass

    #     # Glucose in
    #     WasteWater.imass['Glucose'] += glucose.imass['Glucose']
    #     WasteWater.F_mass += glucose.F_mass

    #     # Simulate complete glucose consumption
    #     consumed_glucose = WasteWater.imass['Glucose']
    #     WasteWater.imass['Glucose'] -= consumed_glucose
    #     WasteWater.F_mass -= consumed_glucose

    #     # Adjust N2O emission factor based on reduction ratio
    #     original_N2O_EF = self.N2O_EF_decay
    #     self.N2O_EF_decay *= 0.25

    #     # Call Decay._first_order_run with the updated N2O_EF_decay
    #     super()._first_order_run(waste=WasteWater,
    #                              treated=TreatedWater,
    #                              CH4=CH4_emission,
    #                              N2O=N2O_emission
    #                              )

    #     # Restore the original N2O emission factor
    #     self.N2O_EF_decay = original_N2O_EF

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # assume linear scale
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        # massflow_anoxic = self.ins[0].mass
        C['Tank'] = self.anoxic_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
        # C['Chemcial_glucose'] = self.chemical_glucose_dosage * massflow_anoxic * self.chemical_glucose_price  # make sense the unit of treated water flow
        # Glucose cost is already accounted for in the WasteStream
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.power_demand_AnoxicTank
        self.power_utility(power_demand)
    
    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        Anoxic_tank_replacement_cost = (self.anoxic_tank_cost /self.anoxic_tank_lifetime +
                                        self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
                                        self.pipeline_connectors / self.pipeline_connectors_lifetime) * scale
        Anoxic_tank_replacement_cost = Anoxic_tank_replacement_cost / (365 * 24)  # convert to USD/hr
        return Anoxic_tank_replacement_cost

# %%
Aerobic_path = ospath.join(EL_su_data_path, '_EL_Aerobic.tsv')

@price_ratio()
class EL_Aerobic(SanUnit, Decay):

    '''
    Aerobic treatment unit in the EL system.

    Parameters
    ----------
    Ins:
    (1) influent of treated wastetwater from anoxic tank
    (2) PAC addition
    (3) blower

    Outs:
    (1) effluent of treated wastewater
    (2) fugitive CH4 emission
    (3) fugtivie N2O emission

    Attributes
    ----------

    
    References
    ----------
     refer to the qsdsan.sanunits.PrimaryClarifier module
    '''

    _N_ins = 3  # treated water, PAC, blower
    _N_outs = 3  # treated water, CH4, N2O
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    exponent_scale = 0.1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), ppl = None, baseline_ppl = None, F_BM_default=1,
                 if_capture_biogas=False, if_N2O_emission=True,**kwargs):
        Decay.__init__(self, ID, ins=ins, outs=outs, thermo=thermo,
                       init_with=init_with, F_BM_default=F_BM_default,
                       degraded_components=degraded_components,
                       if_capture_biogas=if_capture_biogas,
                       if_N2O_emission=if_N2O_emission,)

        self.ppl = ppl
        self.baseline_ppl = baseline_ppl

        data = load_data(path=Aerobic_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]  
    
    def _run(self):

        # Input streams
        WasteWater = self.ins[0]
        PAC = self.ins[1]
        air = self.ins[2]
        
        # Output streams
        TreatedWater = self.outs[0]
        CH4_emission = self.outs[1] 
        N2O_emission = self.outs[2]
        
        # Inherited input stream
        TreatedWater.copy_like(WasteWater)
        
        # Manually add return nitrate (NO3)
        # TreatedWater.imass['PAC'] += PAC.imass['PAC']
        PAC.imass['PAC'] = WasteWater.imass['H2O'] * self.chemical_PAC_dosage/100
        
        # COD removal
        COD_removal = self.EL_aeroT_COD_removal
        removed_COD = WasteWater.COD / 1e3 * WasteWater.F_vol * COD_removal  # kg/hr
        
        # Sludge produced
        sludge_prcd = self.EL_aeroT_sludge_yield * removed_COD  # produced biomass
          
        for component in ('Mg', 'Ca', 'OtherSS', 'Tissue', 'WoodAsh'):
            TreatedWater.imass[component] += sludge_prcd  # New sludge produced
            
        # CH4 emission
        CH4_emission.imass['CH4'] += TreatedWater.imass['SolubleCH4']  # Let all soluble CH4 transfer from solution phase to gas phase
        TreatedWater.imass['SolubleCH4'] = 0  # Ensure that treated water will not include soluble CH4
        
        # N2O produced
        N_removal = self.EL_aeroT_TN_removal
        if self.if_N2O_emission:
          # Assume the removal part covers both degradation loss
          # and other unspecified removal as well
              N_loss = self.first_order_decay(k = self.decay_k_N,
                                    t = self.EL_aeroT_tau / 365,
                                    max_decay = self.N_max_decay)
              if N_loss > N_removal:
                  warn(f'Nitrogen degradation loss ({N_loss:.2f}) larger than the given removal ({N_removal:.2f})), '
                        'the degradation loss is assumed to be the same as the removal.')
                  N_loss = N_removal
            
              # N2O only from the degraded part
              N_loss_tot = N_loss * WasteWater.TN / 1e3 * WasteWater.F_vol
              N2O_emission.imass['N2O'] = N_loss_tot * self.N2O_EF_decay * 44 / 28
        else:
              N2O_emission.empty()
              N2O_emission.imass['N2O'] = 0  # All N2O is emitted
        
        # NO3 conversion
        NH3_mass = TreatedWater.imass['NH3']  # Inherite NH4 property from anoxic tank
        NO3_mass_generated = NH3_mass * self.NO3_produced_ratio
        TreatedWater.imass['NO3'] += NO3_mass_generated
        TreatedWater.imass['NH3'] = 0
        
        # N2O emission
        N2O_mass_generated = NH3_mass * (1 - self.NO3_produced_ratio)
        N2O_emission.imass['N2O'] += N2O_mass_generated
    
    # def _run(self):

    #     # Input streams
    #     WasteWater = self.ins[0]
    #     PAC = self.ins[1]
    #     air = self.ins[2]

    #     # Output stream
    #     TreatedWater = self.outs[0] 
    #     CH4_emission = self.outs[1] 
    #     N2O_emission = self.outs[2]
        
    #     input_streams = [WasteWater, PAC, air]
            
    #     # Mix all inputs into a single stream
    #     self._mixed.empty()
    #     self._mixed.mix_from(input_streams)

    #     # Copy the mixed result to the outflow
    #     TreatedWater.copy_like(self._mixed)


    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        # massflow_aerobic = self.ins[0].mass #kg/h
        C['Tank'] = self.aerobic_tank_cost
        C['Pipes'] = self.pipeline_connectors
        C['Fittings'] = self.weld_female_adapter_fittings
        # C['Chemical_PAC'] = self.chemical_PAC_dosage * massflow_aerobic * self.chemical_PAC_price
        #PAC cost is already accounted in WasteStream

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.power_demand_AerobicTank
        self.power_utility(power_demand) # kWh
    
    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) * self.exponent_scale
        Aerobic_tank_replacement_cost = (self.aerobic_tank_cost / self.aerobic_tank_lifetime +
                                        self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
                                        self.pipeline_connectors / self.pipeline_connectors_lifetime) * scale
        Aerobic_tank_replacement_cost = Aerobic_tank_replacement_cost / (365 * 24)  # convert to USD/hr
        return Aerobic_tank_replacement_cost

# %%
MBR_path = ospath.join(EL_su_data_path, '_EL_MBR.tsv')

@price_ratio()
class EL_MBR(SanUnit, Decay):

    '''
    Aerobic treatment unit in the EL system.

    Parameters
    ----------
    Ins:
    (1) influent of treated wastetwater from aerobic tank
    (2) blower

    Outs:
    (1) effluent of treated wastewater
    (2) nitrate return flow to anoxic tank
    (3) nitrate return flow to primary clarifier
    (4) fugitive CH4 emission
    (5) fugtivie N2O emission

    Attributes
    ----------

    
    References
    ----------
     refer to the exposan.eco-san.MBR module
    '''
    _N_ins = 2  # treated water, Blower
    _N_outs = 6  # treated water, CH4, N2O, Nitrate return to Primary Clarifier, Nitrate return to Anoxic Tank
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    exponent_scale = 0.1

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 degraded_components=('OtherSS',), ppl = None, baseline_ppl = None, F_BM_default=1,
                 if_capture_biogas=False, if_N2O_emission=True, **kwargs):
        Decay.__init__(self, ID, ins=ins, outs=outs, thermo=thermo,
                       init_with=init_with, F_BM_default=F_BM_default,
                       degraded_components=degraded_components,
                       if_capture_biogas=if_capture_biogas,
                       if_N2O_emission=if_N2O_emission,)

        self.ppl = ppl
        self.baseline_ppl = baseline_ppl

        data = load_data(path=MBR_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),
                             Construction(item ='PVDF_membrane', linked_unit=self, quantity_unit='kg'),]
    
    def _run(self):
        
        # Input stream
        WasteWater = self.ins[0]
        air = self.ins[1]
        
        # Output stream
        TreatedWater = self.outs[0]
        PC_nitrate_return = self.outs[1]
        AnoT_nitrate_return = self.outs[2]
        CH4_emission = self.outs[3]
        N2O_emission = self.outs[4]
        sludge = self.outs[5]
        
        # Inherited input stream
        TreatedWater.copy_like(WasteWater)
        # TreatedWater.empty()
    
        # # Transfer 99% of components to sludge
        sludge.empty()
        for component in ('P','K','NH3','NonNH3','Mg', 'Ca', 'OtherSS', 'Tissue', 'WoodAsh'):
            mass_in_treated = TreatedWater.imass[component]  # Obtain every components' property
            mass_to_sludge = 0.99 * mass_in_treated          # Transfer 99% components content to sludge
            sludge.imass[component] = mass_to_sludge
            TreatedWater.imass[component] -= mass_to_sludge  # The last components content to treated water
        
        # COD removal
        COD_removal = self.EL_mbrT_COD_removal
        removed_COD = WasteWater.COD / 1e3 * WasteWater.F_vol * COD_removal  # kg/hr
        
        # Sludge produced
        sludge_prcd = removed_COD * self.EL_mbrT_sludge_yield  # Typical yield: 0.5 kg VSS/kg COD
            
        # Sludge handling
        sludge.copy_flow(TreatedWater, ('P', 'K', 'NH3','NonNH3', 'Mg', 'Ca', 'OtherSS', 'Tissue', 'WoodAsh'), remove=True)
        sludge.imass['OtherSS'] += sludge_prcd
        sludge.imass['H2O'] = sludge.imass['OtherSS']/(1-self.sludge_moisture_content) - sludge.imass['OtherSS']
        
        # Water balance
        TreatedWater.imass['H2O'] -= sludge.imass['H2O']
        if TreatedWater.imass['H2O'] <= 0:
            sludge.imass['H2O'] = WasteWater.imass['H2O']
            TreatedWater.imass['H2O'] = 0
        
        # CH4 emission
        CH4_emission.imass['CH4'] = WasteWater.imass['SolubleCH4']  # Let all soluble CH4 transfer from solution phase to gas phase
        TreatedWater.imass['SolubleCH4'] = 0  # Ensure that treated water will not include soluble CH4
        
        # N2O produced
        N_removal = self.EL_mbrT_TN_removal
        if self.if_N2O_emission:
          # Assume the removal part covers both degradation loss
          # and other unspecified removal as well
              N_loss = self.first_order_decay(k = self.decay_k_N,
                                    t = self.EL_mbrT_tau / 365,
                                    max_decay = self.N_max_decay)
              if N_loss > N_removal:
                  warn(f'Nitrogen degradation loss ({N_loss:.2f}) larger than the given removal ({N_removal:.2f})), '
                        'the degradation loss is assumed to be the same as the removal.')
                  N_loss = N_removal
            
              # N2O only from the degraded part
              N_loss_tot = N_loss * WasteWater.TN / 1e3 * WasteWater.F_vol
              N2O_emission.imass['N2O'] = N_loss_tot * self.N2O_EF_decay * 44 / 28
        else:
              N2O_emission.empty()
        
        # NO3 conversion
        NH3_mass = TreatedWater.imass['NH3']  # Inherite NH4 property from anoxic tank
        NO3_mass_generated = NH3_mass * self.NO3_produced_ratio
        TreatedWater.imass['NO3'] += NO3_mass_generated
        TreatedWater.imass['NH3'] = 0
        
        # N2O emission
        N2O_mass_generated = NH3_mass * (1 - self.NO3_produced_ratio)
        N2O_emission.imass['N2O'] += N2O_mass_generated
        
        # Split NO3 into PC_nitrate_return and AnoT_nitrate_return
        NO3_return_ratio = self.NO3_split_ratio  # Ratio of NO3 going to PC vs. AnoT
        PC_nitrate_return.imass['NO3'] = TreatedWater.imass['NO3'] * NO3_return_ratio
        AnoT_nitrate_return.imass['NO3'] = TreatedWater.imass['NO3'] * (1 - NO3_return_ratio)

        # Remove NO3 from TreatedWater since it is returned to the previous tanks
        TreatedWater.imass['NO3'] -= (PC_nitrate_return.imass['NO3'] + AnoT_nitrate_return.imass['NO3'])
        
        # Update COD
        TreatedWater._COD = WasteWater.COD * (1 - self.EL_mbrT_COD_removal)
        sludge._COD = WasteWater.COD * self.EL_mbrT_COD_removal
    
    # def _run(self):

    #     # Input streams
    #     WasteWater = self.ins[0]
    #     air = self.ins[1]

    #     # Output stream
    #     TreatedWater = self.outs[0] 
    #     PC_sludge_return = self.outs[1]
    #     AnoT_sludge_return = self.outs[2]
    #     CH4_emission = self.outs[3] 
    #     N2O_emission = self.outs[4]
        
    #     input_streams = [WasteWater, air]
            
    #     # Mix all inputs into a single stream
    #     self._mixed.empty()
    #     self._mixed.mix_from(input_streams)

    #     # Copy the mixed result to the outflow
    #     TreatedWater.copy_like(self._mixed)
    
    # def _run(self):
        
    #     # Input streams
    #     WasteWater = self.ins[0]
    #     air = self.ins[1]

    #     # Output streams
    #     TreatedWater = self.outs[0]
    #     PC_sludge_return = self.outs[1]
    #     AnoT_sludge_return = self.outs[2]
    #     CH4_emission = self.outs[3]
    #     N2O_emission = self.outs[4]

    #     # Air does not affect mass balance
    #     # Air is used for aeration, so we do not modify WasteWater mass

    #     # Step 1: Copy WasteWater to TreatedWater
    #     TreatedWater.copy_like(WasteWater)

    #     # Step 2: Split sludge for returns (1% to PC, 99% to AnoT)
    #     total_f_sludge = 0.05  # Total sludge accounts for 5% of the wastewater mass
    #     f_PC_sludge = 0.01     # Primary clarifier accounts for 1% of the total sludge
    #     f_AnoT_sludge = 0.99   # Anoxic tank accounts for 99% of the total sludge

    #     # Sludge mass
    #     total_sludge_mass = WasteWater.F_mass * total_f_sludge
    #     PC_sludge_return.mass = total_sludge_mass * f_PC_sludge
    #     AnoT_sludge_return.mass = total_sludge_mass * f_AnoT_sludge

    #     # Update TreatedWater after sludge removal
    #     TreatedWater.F_mass -= total_sludge_mass
    #     for component in TreatedWater.components:
    #         TreatedWater.imass[component] -= (
    #             PC_sludge_return.imass[component] + AnoT_sludge_return.imass[component]
    #         )

    #     # Step 3: Calculate emissions using Decay._first_order_run
    #     super()._first_order_run(waste=WasteWater,
    #                              treated=TreatedWater,
    #                              biogas=None,  # Assume no biogas is captured
    #                              CH4=CH4_emission,
    #                              N2O=N2O_emission
    #                              )

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # assume linear scaling
        design['PVDF_membrane'] = constr[1].quantity = self.membrane_material_weight
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
        C['MBR_tank'] = self.MBR_tank_cost
        C['pipeline'] = self.pipeline_connectors
        C['fittings'] = self.weld_female_adapter_fittings
        C['Membrane_material'] = self.membrane_material_price * self.membrane_material_weight
        C['Membrane_cleaning'] = self.membrane_cleaning_fee

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio 
        
        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.power_demand_MBR #TODO: power_demand unit should be kW, but the tsv file is in kWh/day
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        MBR_replacement_cost = (
            self.MBR_tank_cost / self.MBR_tank_lifetime +
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
            self.membrane_material_price * self.membrane_material_weight / self.membrane_material_lifetime
            ) * scale
        MBR_replacement_cost = MBR_replacement_cost / (365 * 24) * self.price_ratio # USD/hr
        return MBR_replacement_cost

# %%
ClearWaterTank_path = ospath.join(EL_su_data_path, '_EL_CWT.tsv')

@price_ratio()
class EL_CWT(StorageTank):

    '''
    Introduction
    ------------
    To only collect the treated water

    Parameters
    ----------
    Ins:
    (1) influent of treated wastetwater from self-priming pump
    (2) O3 dosing
    (3) air-dissolve influent

    Outs:
    (1) effluent of treated wastewater (clear water)
    (2) spill flow to collection tank


    Attributes
    ----------

    
    References
    ----------
     refer to the qsdsan.sanunits.storagetank module

    '''
    _N_ins = 3  # number of ins
    _N_outs = 2  # number of outs
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    exponent_scale = 0.1

    def __init__(self, ID = '', ins = None, outs = (), thermo = None, ppl = None, baseline_ppl = None,
                 vessel_type='Field erected',tau = 24, V_wf = None, vessel_material='Stainless steel', kW_per_m3 = 0.1,
                 init_with = 'WasteStream', F_BM_default = 1, max_overflow=15, include_construction = True, **kwargs):
        StorageTank.__init__(self, ID=ID, ins=ins, outs=outs, thermo = thermo, init_with = init_with, 
                             F_BM_default = F_BM_default, include_construction = include_construction,
                             kW_per_m3= kW_per_m3, vessel_type= vessel_type, tau= tau,
                             vessel_material= vessel_material, V_wf= V_wf)

        self.ppl = ppl
        self.baseline_ppl = baseline_ppl
        self.max_overflow = max_overflow

        data = load_data(path = ClearWaterTank_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]
    
    def _run(self):
        
        # Input streams
        WasteWater = self.ins[0]
        ozone = self.ins[1]
        air = self.ins[2]
        
        # Output streams
        TreatedWater = self.outs[0]
        spill_return = self.outs[1]
        # CWT_cycle = self.outs[2]
        
        # Inherited input stream
        TreatedWater.copy_like(WasteWater)
        
        # Spill water return
        max_overflow_m3d = self.max_overflow * 24  # Convert m³/h to m³/d
        Q_treated = TreatedWater.F_vol * 24  # m³/d
        if Q_treated > max_overflow_m3d:
            spill_vol = Q_treated - max_overflow_m3d
            f_spill = spill_vol / Q_treated
            spill_return.copy_like(TreatedWater)
            spill_return.F_mass *= f_spill
            TreatedWater.F_mass *= (1 - f_spill)
        else:
            # max_overflow is none, no spill return
            spill_return.empty()
        

        # # Input streams
        # WasteWater = self.ins[0]
        # ozone = self.ins[1]
        # air = self.ins[2]

        # # Output streams
        # TreatedWater = self.outs[0]
        # spill_return = self.outs[1]
        # CWT_cycle = self.outs[2]
        
        # input_streams = [WasteWater, ozone, air]
            
        # # Mix all inputs into a single stream
        # self._mixed.empty()
        # self._mixed.mix_from(input_streams)
        

        # # Copy the mixed result to the outflow
        # TreatedWater.copy_like(self._mixed)
        
    
    
    # def _run(self):
        
    #     # Input streams
    #     WasteWater = self.ins[0]
    #     ozone = self.ins[1]
    #     air = self.ins[2]

    #     # Output streams
    #     TreatedWater = self.outs[0]
    #     CT_spill_return = self.outs[1]
    #     PC_spill_return = self.outs[2]
 
    #     # Constants for overflow conditions
    #     CT_spill_fraction = 0.6
    #     PC_spill_fraction = 0.4

    #     # Step 1: Process inputs
    #     # Air and ozone are balanced, not affecting mass balance
    #     TreatedWater.copy_like(WasteWater)

    #     # Step 2: Handle overflow if necessary
    #     if TreatedWater.F_vol > self.max_overflow:
    #         # Calculate excess volume
    #         spill_vol = TreatedWater.F_vol - self.max_overflow

    #         # Calculate spill fractions
    #         self._f_spill = spill_vol / TreatedWater.F_vol
    #         self._f_treated = 1 - self._f_spill

    #         # Assign spill fractions to CT_spill_return and PC_spill_return
    #         CT_spill_return.F_mass[:] = TreatedWater.F_mass[:] * self._f_spill * CT_spill_fraction
    #         PC_spill_return.F_mass[:] = TreatedWater.F_mass[:] * self._f_spill * PC_spill_fraction

    #         # Adjust the normal treated water to within max capacity
    #         TreatedWater.F_mass[:] *= self._f_treated
    #     else:
    #         # No spill return if overflow is within capacity
    #         CT_spill_return.empty()
    #         PC_spill_return.empty()
    #         self._f_spill = 0.0
    #         self._f_treated = 1.0

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # to be defined in .tsv file
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
        C['Clear water tank'] = self.clear_water_tank_cost
        C['pipeline'] = self.pipeline_connectors
        C['fittings'] = self.weld_female_adapter_fittings
        C['O3 generator'] = self.O3_generation_machine_cost

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.power_demand_CWT
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        CWR_replacement_cost = (
            self.clear_water_tank_cost / self.clear_water_tank_lifetime +
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime +
            self.O3_generation_machine_cost / self.O3_generation_machine_lifetime
            ) * scale
        CWR_replacement_cost = CWR_replacement_cost / (365 * 24) # convert to USD/hr
        return CWR_replacement_cost

# %%
PressureTank_path = ospath.join(EL_su_data_path, '_EL_PT.tsv')

@price_ratio()
class EL_PT(StorageTank):
    
    '''
    Introduction
    ------------
    To only collect the clear water within pressure tank

    Parameters
    ----------
    Ins:
    (1) influent of pressurized clear water from clear water tank

    Outs:
    (1) recycle to flush toilet


    Attributes
    ----------

    
    References
    ----------
     refer to the qsdsan.sanunits.storagetank module

    '''
    _N_ins = 1  # number of ins
    _N_outs = 2  # number of outs
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    exponent_scale = 0.1

    def __init__(self, ID = '', ins = None, outs = (), thermo = None, ppl = None, baseline_ppl = None,
                 vessel_type='Field erected',tau = 24, V_wf = None, vessel_material='Stainless steel', kW_per_m3 = 0.1,
                 init_with = 'WasteStream', F_BM_default = 1, include_construction = True, **kwargs):
        StorageTank.__init__(self, ID=ID, ins=ins, outs=outs, thermo = thermo, init_with = init_with, 
                             F_BM_default = F_BM_default, include_construction = include_construction,
                             kW_per_m3= kW_per_m3, vessel_type= vessel_type, tau= tau,
                             vessel_material= vessel_material, V_wf= V_wf)

        self.ppl = ppl
        self.baseline_ppl = baseline_ppl

        self._mixed = WasteStream()

        data = load_data(path = PressureTank_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]
    
    def _run(self):

        # Input streams
        ClearWater = self.ins[0]

        # Output streams
        TreatedWater = self.outs[0]
        nullConnection = self.outs[1]
        
        # Inherited input stream
        TreatedWater.copy_like(ClearWater)
        # TreatedWater.empty()
        nullConnection.empty()
        
        
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.tank_steel_volume * self.steel_density * (self.ppl / self.baseline_ppl)  # to be defined in .tsv file
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
        C['Pressure water tank'] = self.pressure_water_tank_cost
        C['pipeline'] = self.pipeline_connectors
        C['fittings'] = self.weld_female_adapter_fittings

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost()

        power_demand = self.power_demand_PT
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        PW_replacement_cost = (
            self.pressure_water_tank_cost / self.pressure_water_tank_lifetime +
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.weld_female_adapter_fittings_lifetime) * scale
        PW_replacement_cost = PW_replacement_cost / (365 * 24) # convert to USD/hr
        return PW_replacement_cost

# %%
blower_path = ospath.join(EL_su_data_path, '_EL_blower.tsv')

@price_ratio()
class EL_blower(SanUnit):
    '''
    Introduction
    ------------
    To areate air for aerobic tank and membrane tank

    Parameters
    ----------
    Ins:
    (1) air

    Outs:
    (1) air


    Attributes
    ----------

    
    References
    ----------
     refer to the qsdsan.equipments.Blower module

    '''
    _N_ins = 1  # number of ins
    _N_outs = 1  # number of outs
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    exponent_scale = 0.1

    def __init__(self, ID = '', ins = None, outs = (), init_with = 'WasteStream',
                 # F_BM={
                 #     'Blowers': 2.22,
                 #     'Blower piping': 1,
                 #     'Blower building': 1.11,
                 #     },
                 F_BM = 2.22,
                 lifetime=15, lifetime_unit='yr',
                 # units={
                 #     'Total gas flow': 'CFM',
                 #     'Blower capacity': 'CFM',
                 #     'Number of blowers': '',
                 #     'Total blower power': 'kW',
                 #     },
                 # N_reactor=2, # the number of the reactors where the gas sparging modules will be installed
                 # gas_demand_per_reactor=1, # gas demand per reactor
                 # TDH=6, # total dynamic head for rhe blower, in psi
                 # eff_blower=0.7, # efficiency of the blower in fraction
                 # eff_motor=0.7, # efficiency of the motor in fraction
                 # AFF=3.33, # air flow fraction
                 # building_unit_cost=9, # unit cost of the building, in USD/ft2
                 thermo = None, ppl = None, baseline_ppl = None, **kwargs):
        # super().__init__(ID=ID, lifetime = lifetime, lifetime_unit = lifetime_unit, F_BM=F_BM,
        #                 units=units, N_reactor=N_reactor, gas_demand_per_reactor=gas_demand_per_reactor,
        #                 TDH=TDH, eff_blower=eff_blower, eff_motor=eff_motor, AFF=AFF, building_unit_cost=building_unit_cost,)
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=F_BM)

        self.ppl = ppl
        self.baseline_ppl = baseline_ppl
        self.lifetime = lifetime
        self.lifetime_unit = lifetime_unit

        data = load_data(path = blower_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [Construction(item='StainlessSteel', linked_unit=self, quantity_unit='kg'),]
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.blower_steel_weight  # to be defined in .tsv file
        self.add_construction(add_cost=False)
    

    def _cost(self):
        C = self.baseline_purchase_costs # the below items need to be defined in .tsv file
        C['Blower'] = self.blower_cost
        C['pipeline'] = self.pipeline_connectors
        C['fittings'] = self.weld_female_adapter_fittings

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost()
        
        power_demand = self.power_demand_blower
        self.power_utility(power_demand)
    
    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        Blower_replacement_cost = (
            self.blower_cost * self.blower_lifetime +
            self.pipeline_connectors / self.pipeline_connectors_lifetime +
            self.weld_female_adapter_fittings / self.fittings_lifetime) * scale
        Blower_replacement_cost = Blower_replacement_cost / (365 * 24) # convert to USD/hr
        return Blower_replacement_cost

# %%
housing_path = ospath.join(EL_su_data_path, '_EL_housing.tsv')

#@price_ratio()
class EL_Housing(SanUnit):
    '''
     non_reactive unit for the Enviroloo Clear system
    '''
    _N_ins = 1  # number of ins
    _N_outs = 1  # number of outs
    _ins_size_is_fixed = True
    _outs_size_is_fixed = True
    ppl_per_MURT = 30  # number of people per MURT

    def __init__(self, ID = '', ins = None, outs = (), thermo = None, init_with = None,
                 price_ratio=0.9,
                 ppl = 1000, baseline_ppl = 30, F_BM_default= 1, **kwargs):
        init_with = init_with or {}
        super().__init__(ID=ID, ins=ins, outs=outs, thermo = thermo, 
                         init_with = init_with, F_BM_default=F_BM_default)
        
        self.ppl = ppl
        self.baseline_ppl = baseline_ppl
        self.price_ratio = price_ratio

        data = load_data(path = housing_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self): # replace the actual materials used in the EL
        self.construction = [
            Construction(item = 'StainlessSteel', linked_unit= self, quantity_unit= 'kg'),
            Construction(item = 'Plastic', linked_unit= self, quantity_unit= 'kg'),]

    def _design(self): # replace the actual materials used in the EL
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = (self.steel_weight + self.steel_framework_weight + self.steel_fittings_weight) * (self.ppl / self.baseline_ppl)  # assume linear scaling
        design['Plastic'] = constr[1].quantity = (self.LLDPE_weight) * (self.ppl / self.baseline_ppl)   # assume linear scaling
        self.add_construction(add_cost= False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Housing'] = (self.frame + self.extrusion + 
                        self.angle_frame + self.angle +
                        self.door_sheet + self.plate +
                        self.powder_coating) * (1 + 0.1 * (self.N_EL -1))
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
    
    @property
    def N_EL(self): # determine the number of EL system needed
        return ceil(self.ppl / self.baseline_ppl)
    
    @property
    def N_toilets(self): # determine the number of toilets needed
        return ceil(self.ppl / self.ppl_per_MURT)

# %%
system_path = ospath.join(EL_su_data_path, '_EL_system.tsv')

@price_ratio()
class EL_System(SanUnit, isabstract=True):
    '''
    Relate to connection components in the EL system
    '''
    _N_ins = 1
    _N_outs = 0
    exponent_scale = 0.1

    def __init__(self, ID='', ins=(), outs=None, thermo = None, 
                 init_with = 'WasteStream',
                 # init_with=None,
                 if_gridtied = True, ppl = None, baseline_ppl = None, F_BM_default = 1, **kwargs):
        SanUnit.__init__(self, ID=ID, ins=ins, outs=outs, thermo = thermo, init_with = init_with,  F_BM_default = F_BM_default)
        
        self.ppl = ppl
        self.baseline_ppl = baseline_ppl
        self.if_gridtied = if_gridtied

        data = load_data(path = system_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [
            Construction(item = 'PVC_generic', linked_unit= self, quantity_unit= 'kg'),
            Construction(item = 'HDPE', linked_unit= self, quantity_unit= 'kg'),
            ]

    def _design(self):
        design = self.design_results
        design['PVC_generic'] = self.construction[0].quantity = self.PVC_weight * (self.ppl / self.baseline_ppl)
        design['HDPE'] = self.construction[1].quantity = self.HDPE_weight * (self.ppl / self.baseline_ppl)
        self.add_construction(add_cost= False)
    
    def _cost(self): # replace these items below by that listed in the _EL_system.tsv
        C = self.baseline_purchase_costs
        C['System'] = (
            self.membrane_filters_M +
            self.membrane_filters_size +
            self.membrane_filters_pause_size +
            self.membrane_filters_chassis_M +
            self.membrane_filters_air_diffuser +
            self.membrane_filters_air_diffuser_chassis +
            self.overflow_membrane2collection +
            self.overflow_clear_water2collection +
            self.overflow_primary_clarifier2anoixc +
            self.overflow_anoxic2aerobic +
            self.overflow_aerobic2membrane +
            self.ozone_pipeline_200m +
            self.pipeline_fittings_32mm +
            self.pipeline_fittings_40mm +
            self.pipeline_fittings_60mm +
            self.ball_valves_50mm +
            self.aerobic_air_diffuser)

        ratio = self.price_ratio # ratio of the price of the new system to the baseline system
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX = self._calc_replacement_cost()  # add the cost of replacement

        if self.if_gridtied:
            power_demand = (self.power_demand_system / 1000) * self.N_EL  # in W/d
        else:
            power_demand = 0
        
        self.power_utility(power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        system_replacement_cost = (
            self.membrane_filters_M / self.lifetime_membrane_filters_M +
            self.membrane_filters_size / self.lifetime_membrane_filters_size +
            # self.aerobic_basin / self.lifetime_aerobic_basin +
            self.membrane_filters_air_diffuser / self.lifetime_membrane_filters_air_diffuser +
            self.membrane_filters_air_diffuser_chassis / self.lifetime_membrane_filters_air_diffuser_chassis +
            self.overflow_membrane2collection / self.lifetime_overflow_membrane2collection +
            self.overflow_clear_water2collection / self.lifetime_overflow_clear_water2collection +
            self.overflow_primary_clarifier2anoixc / self.lifetime_overflow_primary_clarifier2anoixc +
            self.overflow_anoxic2aerobic / self.lifetime_overflow_anoxic2aerobic +
            self.overflow_aerobic2membrane / self.lifetime_overflow_aerobic2membrane +
            self.ozone_pipeline_200m / self.lifetime_200m_ozone_pipeline +
            self.pipeline_fittings_32mm / self.lifetime_32mm_pipeline_fittings +
            self.pipeline_fittings_40mm / self.lifetime_40mm_pipeline_fittings +
            self.pipeline_fittings_60mm / self.lifetime_60mm_pipeline_fittings +
            self.ball_valves_50mm / self.lifetime_50mm_ball_valves) * scale
        system_replacement_cost = system_replacement_cost / (365 * 24)   # convert from USD/year to USD/hour
        return system_replacement_cost
    @property
    def N_EL(self): # determine the number of EL system needed
       return ceil(self.ppl / self.baseline_ppl)
