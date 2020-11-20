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
from ..utils.loading import load_data, data_path

__all__ = ('Toilet',)

data_path += 'unit_data/Toilet.csv'


# %%

class Toilet(SanUnit, isabstract=True):
    '''
    Abstract class to hold common parameters for PitLatrine and UDDT
    (urine-diverting dry toliet).
    '''
    
    def __init__(self, ID='', ins=None, outs=(), N_user=1, life_time=8,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                 OPEX_over_CAPEX=None):
        '''
        
        Parameters
        ----------
        N_user : [float]
            Number of people that use the toilet per hour.
        life_time : [float]
            Life time of the toilet in year.
        if_toilet_paper : [bool]
            If toilet paper is used.
        if_flushing : [bool]
            If water is used for flushing.
        if_cleansing : [bool]
            If water is used for cleansing.
        if_desiccant : [bool]
            If desiccant is used for moisture and odor control.
        if_air_emission : [bool]
            If emission to air occurs
            (i.e., if the pit is completely sealed off from the atmosphere).
        if_ideal_emptying : [bool]
            If the toilet appropriately emptied to avoid contamination to the
            environmental.
        OPEX_over_CAPEX : [float]
            Fraction of annual operating cost over total capital cost.
        '''

        SanUnit.__init__(self, ID, ins, outs)
        self.N_user = N_user
        self.life_time = life_time       
        self.if_toilet_paper = if_toilet_paper       
        self.if_flushing = if_flushing       
        self.if_cleansing = if_cleansing
        self.if_desiccant = if_desiccant
        self.if_air_emission = if_air_emission
        self.if_ideal_emptying = if_ideal_emptying
        self.OPEX_over_CAPEX = OPEX_over_CAPEX

        data = load_data(path=data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            if para in ('desiccant_V', 'desiccant_rho'):
                setattr(self, para, value)
            else:
                setattr(self, '_'+para, value)
        del data
        
        self._empty_ratio = 0.59
        # Assuming tau_deg is 2 yr and log_deg is 3
        self._decay_k = (-1/2)*np.log(10**-3)
        self._max_CH4_emission = 0.25
        
    _N_ins = 6
    _outs_size_is_fixed = False
    _units = {
        'Cement': 'kg',
        'Sand': 'kg',
        'Gravel': 'kg',
        'Bricks': 'kg',
        'Plastic': 'kg',
        'Steel': 'kg',
        'Stainless steel sheet': 'kg',
        'Wood': 'm3',
        'Excavation': 'm3'
        }
    _BM = {'Toilet': 1}

    def _run(self):
        ur, fecs, tp, fw, cw, des = self.ins
        
        N_user = self.N_user
        tp.imass['Tissue'] = int(self.if_toilet_paper)*self.toilet_paper * N_user
        fw.imass['H2O'] = int(self.if_flushing)*self.flushing_water * N_user
        cw.imass['H2O'] = int(self.if_cleansing)*self.cleansing_water * N_user
        des.imass['WoodAsh'] = int(self.if_desiccant)*self.desiccant * N_user
      
    @staticmethod
    def _allocate_N_reduction(tot_red, NH3):
        '''
        Allocate the total amount of N removal to NH3 and non-NH3 Components.
        NH3 will be firstly removed before non-NH3.

        Parameters
        ----------
        tot_red : [float]
            Total amount of N to be removed.
        NH3 : [float]
            Current NH3 content.

        Returns
        -------
        [float]
            Amount of NH3 to be removed.
        [float]
            Amount of non-NH3 to be removed.

        '''
        if not NH3 > 0:
            return 0, tot_red
        elif NH3 > tot_red:
            return tot_red, 0
        else:
            return NH3, tot_red-NH3  
      
    @staticmethod
    def get_degradation_loss(k, t, max_removal, tot=1):
        '''
        To calculate first-order degradation loss.

        Parameters
        ----------
        k : [float]
            Degradation rate constant.
        t : [float]
            Degradation time.
        max_removal : [float]
            Maximum removal ratio.
        tot : [float], optional
            Total degradable amount.
            If set to 1 (default), the return is the relative ratio (i.e., loss/tot).

        Returns
        -------
        loss : [float]
            Amount lost due to degradation.

        '''

        max_deg = tot * max_removal
        after = max_deg/(k*t) * (1-np.exp(-k*t))
        loss = max_deg - after
        return loss

    @staticmethod
    def get_emptying_emission(waste, CH4, N2O, app_ratio, CH4_factor, N2O_factor):
        '''
        Calculate emissions due to non-ideal emptying.

        Parameters
        ----------
        stream : WasteStream
            Excreta stream that is not appropriately empited (before emptying).
        CH4 : WasteStream
            Fugitive CH4 gas (before emptying).
        N2O : WasteStream
            Fugitive N2O gas (before emptying).
        app_ratio : [float]
            Fraction of excreta that is appropriately emptied..
        CH4_factor : [float]
            Factor to convert COD removal to CH4 emission.
        N2O_factor : [float]
            Factor to convert COD removal to N2O emission.

        Returns
        -------
        stream : WasteStream
            Excreta stream that is not appropriately empited (before emptying).
        CH4 : WasteStream
            Fugitive CH4 gas (before emptying).
        N2O : WasteStream
            Fugitive N2O gas (before emptying).

        '''
        COD_rmd = waste.COD*(1-app_ratio)/1e3*waste.F_vol
        CH4.imass['CH4'] += COD_rmd * CH4_factor
        waste._COD *= app_ratio
        N2O.imass['N2O'] += COD_rmd * N2O_factor
        waste.mass *= app_ratio
        return waste, CH4, N2O

    @property
    def toilet_paper(self):
        '''
        [float] Amount of toilet paper used
        (if if_toilet_paper is True), [kg/cap/hr].
        '''
        return self._toilet_paper
    @toilet_paper.setter
    def toilet_paper(self, i):
        self._toilet_paper = float(i)
        
    @property
    def flushing_water(self):
        '''
        [float] Amount of water used for flushing
        (if if_flushing_water is True), [kg/cap/hr].
        '''
        return self._flushing_water
    @flushing_water.setter
    def flushing_water(self, i):
        self._flushing_water = float(i)
    
    @property
    def cleansing_water(self):
        '''
        [float] Amount of water used for cleansing 
        (if if_cleansing_water is True), [kg/cap/hr].
        '''
        return self._cleansing_water
    @cleansing_water.setter
    def cleansing_water(self, i):
        self._cleansing_water = float(i)
        
    @property
    def desiccant(self):
        '''
        [float] Amount of desiccant used (if if_desiccant is True), [kg/cap/hr].

        Note
        ----
            Value set by desiccant_V and desiccant_rho.

        '''
        return self.desiccant_V*self.desiccant_rho

    @property
    def N_vol(self):
        '''
        [float] Fraction of input N that volatizes to the air
        (if if_air_emission is True).
        '''
        return self._N_vol
    @N_vol.setter
    def N_vol(self, i):
        self._N_vol = float(i)

    @property
    def empty_ratio(self):
        '''
        [float] Fraction of excreta that is appropriately emptied.

        Note
        ----
            Will be 1 (i.e., 100%) if if_ideal_emptying is True.

        '''
        if self.if_ideal_emptying:
            return 1.
        return self._empty_ratio
    @empty_ratio.setter
    def empty_ratio(self, i):
        if self.if_ideal_emptying:
            msg = f'if_ideal_emptying is True, the set value {i} is ignored.'
            warn(msg, source=self)
        self._empty_ratio = float(i)

    @property
    def COD_max_removal(self):
        '''[float] Maximum raction of COD removed during storage given sufficient time.'''
        return self._COD_max_removal
    @COD_max_removal.setter
    def COD_max_removal(self, i):
        self._COD_max_removal = float(i)

    @property
    def N_max_removal(self):
        '''[float] Maximumraction of N removed through denitrification during storage given sufficient time.'''
        return self._N_max_removal
    @N_max_removal.setter
    def N_max_removal(self, i):
        self._N_max_removal = float(i)

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
    def decay_k(self):
        '''[float] Rate constant for COD and N decay, [/yr].'''
        return self._decay_k
    @decay_k.setter
    def decay_k(self, i):
        self._decay_k = i

    @property
    def max_CH4_emission(self):
        '''[float] Maximum methane emssion as a fraction of degraded COD.'''
        return self._max_CH4_emission
    @max_CH4_emission.setter
    def max_CH4_emission(self, i):
        self._max_CH4_emission = i




