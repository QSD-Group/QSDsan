#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.

Ref:
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.


'''


# %%

from warnings import warn
from .. import SanUnit
from ._decay import Decay
from ..utils.loading import load_data, data_path

__all__ = ('Toilet',)

data_path += 'unit_data/_toilet.csv'


# %%

class Toilet(SanUnit, Decay, isabstract=True):
    '''Abstract class containing common parameters and design algorithms for toilets.'''
    
    def __init__(self, ID='', ins=None, outs=(), N_user=1, N_toilet=1, life_time=8,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                 OPEX_over_CAPEX=None):
        '''
        
        Parameters
        ----------
        N_user : [float]
            Number of people that share this toilet.
        N_toilet : [float]
            Number of paralle toilets.
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
        self._N_user = 1
        self._N_toilet = 1
        self.N_user = N_user
        self.N_toilet = N_toilet
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
        
    _N_ins = 6
    _outs_size_is_fixed = False

    def _run(self):
        ur, fec, tp, fw, cw, des = self.ins
        
        N_user = self.N_user
        tp.imass['Tissue'] = int(self.if_toilet_paper)*self.toilet_paper * N_user
        fw.imass['H2O'] = int(self.if_flushing)*self.flushing_water * N_user
        cw.imass['H2O'] = int(self.if_cleansing)*self.cleansing_water * N_user
        des.imass['WoodAsh'] = int(self.if_desiccant)*self.desiccant * N_user
      

    @staticmethod
    def get_emptying_emission(waste, CH4, N2O, empty_ratio, CH4_factor, N2O_factor):
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
        empty_ratio : [float]
            Fraction of excreta that is appropriately emptied..
        CH4_factor : [float]
            Factor to convert COD removal to CH4 emission.
        N2O_factor : [float]
            Factor to convert COD removal to N2O emission.

        Returns
        -------
        stream : WasteStream
            Excreta stream that is not appropriately empited (after emptying).
        CH4 : WasteStream
            Fugitive CH4 gas (after emptying).
        N2O : WasteStream
            Fugitive N2O gas (after emptying).

        '''
        COD_rmd = waste.COD*(1-empty_ratio)/1e3*waste.F_vol
        CH4.imass['CH4'] += COD_rmd * CH4_factor
        waste._COD *= empty_ratio
        N2O.imass['N2O'] += COD_rmd * N2O_factor
        waste.mass *= empty_ratio
        return waste, CH4, N2O

    @property
    def N_user(self):
        '''[float] Number of people that use the toilet per hour.'''
        return self._N_user
    @N_user.setter
    def N_user(self, i):
        self._N_user = float(i)
        
    @property
    def N_toilet(self):
        '''[float] Number of parallel toilets.'''
        return self._N_toilet
    @N_toilet.setter
    def N_toilet(self, i):
        self._N_toilet = float(i)

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
    def N_volatilization(self):
        '''
        [float] Fraction of input N that volatizes to the air
        (if if_air_emission is True).
        '''
        return self._N_volatilization
    @N_volatilization.setter
    def N_volatilization(self, i):
        self._N_volatilization = float(i)

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
    def MCF_aq(self):
        '''[float] Methane correction factor for COD lost due to inappropriate emptying.'''
        return self._MCF_aq
    @MCF_aq.setter
    def MCF_aq(self, i):
        self._MCF_aq = float(i)

    @property
    def N2O_EF_aq(self):
        '''[float] Fraction of N emitted as N2O due to inappropriate emptying.'''
        return self._N2O_EF_aq
    @N2O_EF_aq.setter
    def N2O_EF_aq(self, i):
        self._N2O_EF_aq = float(i)






