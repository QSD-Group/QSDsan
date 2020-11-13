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
'''

import os
path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                    'inputs/units.xlsx')
del os

import pandas as pd
import thermosteam as tmo
from warnings import warn
from biosteam.units.decorators import cost
from sanitation import SanUnit
from bwaise._utils import load_data

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction


class Excretion(SanUnit):
    '''
    Estimation of N, P, K, and COD in urine and feces based on dietary intake.

    Parameters
    ----------
    N_user : [float]
        Number of people that uses the toilet per hour.

    '''
    
    _N_ins = 0
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), N_user=1, **kwargs):                
        SanUnit.__init__(self, ID, ins, outs)
        self.N_user = N_user
        data = load_data(path=path, sheet='Excretion')
        for para in data.index:
            setattr(self, '_'+para, data.loc[para]['expected'])
        del data
        for attr, value in kwargs:
            setattr(self, attr, value)

    def _run(self):
        
        ur, fec = self.outs
        ur.empty()
        fec.empty()
        # From g per person per day to kg per hour        
        factor = self.N_user / 24 / 1e3
        e_cal = self.e_cal * self.N_user / 24
        ur_exc = self.ur_exc * factor
        ur_N = (self.p_veg+self.p_anim)*factor*self.N_prot \
           * self.N_exc*self.N_ur
        ur.imass['NH3'] = ur_N * self.N_ur_NH3
        ur.imass['NonNH3'] = ur_N - ur.imass['NH3']
        ur.imass['P'] = (self.p_veg*self.P_prot_v+self.p_anim*self.P_prot_a)*factor \
            * self.P_exc*self.P_ur
        ur.imass['K'] = e_cal/1e3 * self.K_cal/1e3 * self.K_exc*self.K_ur
        ur.imass['Mg'] = self.Mg_ur * factor
        ur.imass['Ca'] = self.Ca_ur * factor
        ur.imass['H2O'] = self.ur_moi * ur_exc
        ur.imass['Other_SS'] = ur_exc - ur.F_mass
        
        fec_exc = self.fec_exc * factor
        fec_N = (1-self.N_ur)/self.N_ur * ur_N
        fec.imass['NH3'] = fec_N * self.N_fec_NH3   
        fec.imass['NonNH3'] = fec_N - fec.imass['NH3']
        fec.imass['P'] = (1-self.P_ur)/self.P_ur * ur.imass['P']
        fec.imass['K'] = (1-self.K_ur)/self.K_ur * ur.imass['K']
        fec.imass['Mg'] = self.Mg_fec * factor
        fec.imass['Ca'] = self.Ca_fec * factor
        fec.imass['H2O'] = self.fec_moi * fec_exc
        fec.imass['Other_SS'] = fec_exc - fec.F_mass
        
        # 14 kJ/g COD, the average lower heating value of excreta,
        # 239 to convert it to kcal/kg COD
        tot_COD = e_cal*self.e_exc / (14*239) # in kg COD/hr
        ur._COD = tot_COD*(1-self.e_fec) / (ur.F_vol/1e3) # in mg/L
        fec._COD = tot_COD*self.e_exc / (fec.F_vol/1e3) # in mg/L
        

    @property
    def e_cal(self):
        '''[float] Caloric intake, [kcal/cap/d].'''
        return self._e_cal
    @e_cal.setter
    def e_cal(self, i):
        self._e_cal = float(i)
    
    @property
    def p_veg(self):
        '''[float] Vegetal protein intake, [g/cap/d].'''
        return self._p_veg
    @p_veg.setter
    def p_veg(self, i):
        self._p_veg = float(i)

    @property
    def p_anim(self):
        '''[float] Animal protein intake, [g/cap/d].'''
        return self._p_anim
    @p_anim.setter
    def p_anim(self, i):
        self._p_anim = float(i)

    @property
    def N_prot(self):
        '''[float] Nitrogen content in protein, [wt%].'''
        return self._N_prot
    @N_prot.setter
    def N_prot(self, i):
        self._N_prot = float(i)

    @property
    def P_prot_v(self):
        '''[float] Phosphorus content in vegetal protein, [wt%].'''
        return self._P_prot_v
    @P_prot_v.setter
    def P_prot_v(self, i):
        self._P_prot_v = float(i)

    @property
    def P_prot_a(self):
        '''[float] Phosphorus content in animal protein, [wt%].'''
        return self._P_prot_a
    @P_prot_a.setter
    def P_prot_a(self, i):
        self._P_prot_a = float(i)

    @property
    def K_cal(self):
        '''[float] Potassium intake relative to caloric intake, [g K/1000 kcal].'''
        return self._K_cal
    @K_cal.setter
    def K_cal(self, i):
        self._K_cal = float(i)

    @property
    def N_exc(self):
        '''[float] Nitrogen excretion factor, [% of intake].'''
        return self._N_exc
    @N_exc.setter
    def N_exc(self, i):
        self._N_exc = float(i)

    @property
    def P_exc(self):
        '''[float] Phosphorus excretion factor, [% of intake].'''
        return self._P_exc
    @P_exc.setter
    def P_exc(self, i):
        self._P_exc = float(i)

    @property
    def K_exc(self):
        '''[float] Potassium excretion factor, [% of intake].'''
        return self._K_exc
    @K_exc.setter
    def K_exc(self, i):
        self._K_exc = float(i)

    @property
    def e_exc(self):
        '''[float] Energy excretion factor, [% of intake].'''
        return self._e_exc
    @e_exc.setter
    def e_exc(self, i):
        self._e_exc = float(i)

    @property
    def N_ur(self):
        '''[float] Nitrogen content of urine, [wt%].'''
        return self._N_ur
    @N_ur.setter
    def N_ur(self, i):
        self._N_ur = float(i)

    @property
    def P_ur(self):
        '''[float] Phosphorus content of urine, [wt%].'''
        return self._P_ur
    @P_ur.setter
    def P_ur(self, i):
        self._P_ur = float(i)

    @property
    def K_ur(self):
        '''[float] Potassium content of urine, [wt%].'''
        return self._K_ur
    @K_ur.setter
    def K_ur(self, i):
        self._K_ur = float(i)

    @property
    def e_fec(self):
        '''[float] Percent of excreted energy in feces, [%].'''
        return self._e_fec
    @e_fec.setter
    def e_fec(self, i):
        self._e_fec = float(i)

    @property
    def N_ur_NH3(self):
        '''[float] Reduced inorganic nitrogen in urine, modeled as NH3, [% of total urine N].'''
        return self._N_ur_NH3
    @N_ur_NH3.setter
    def N_ur_NH3(self, i):
        self._N_ur_NH3 = float(i)

    @property
    def N_fec_NH3(self):
        '''[float] Reduced inorganic nitrogen in feces, modeled as NH3, [% of total feces N].'''
        return self._N_fec_NH3
    @N_fec_NH3.setter
    def N_fec_NH3(self, i):
        self._N_fec_NH3 = float(i)

    @property
    def ur_exc(self):
        '''[float] Urine generated per day, [g/cap/d].'''
        return self._ur_exc
    @ur_exc.setter
    def ur_exc(self, i):
        self._ur_exc = float(i)

    @property
    def fec_exc(self):
        '''[float] Feces generated per day, [g/cap/d].'''
        return self._fec_exc
    @fec_exc.setter
    def fec_exc(self, i):
        self._fec_exc = float(i)

    @property
    def ur_moi(self):
        '''[float] Moisture (water) content of urine, [wt%].'''
        return self._ur_moi
    @ur_moi.setter
    def ur_moi(self, i):
        self._ur_moi = float(i)

    @property
    def fec_moi(self):
        '''[float] Moisture (water) content of feces, [wt%].'''
        return self._fec_moi
    @fec_moi.setter
    def fec_moi(self, i):
        self._fec_moi = float(i)

    @property
    def Mg_ur(self):
        '''[float] Magnesium excreted in urine, [g Mg/cap/d].'''
        return self._Mg_ur
    @Mg_ur.setter
    def Mg_ur(self, i):
        self._Mg_ur = float(i)

    @property
    def Mg_fec(self):
        '''[float] Magnesium excreted in feces, [g Mg/cap/d].'''
        return self._Mg_fec
    @Mg_fec.setter
    def Mg_fec(self, i):
        self._Mg_fec = float(i)

    @property
    def Ca_ur(self):
        '''[float] Calcium excreted in urine, [g Ca/cap/d].'''
        return self._Ca_ur
    @Ca_ur.setter
    def Ca_ur(self, i):
        self._Ca_ur = float(i)

    @property
    def Ca_fec(self):
        '''[float] Calcium excreted in feces, [g Ca/cap/d].'''
        return self._Ca_fec
    @Ca_fec.setter
    def Ca_fec(self, i):
        self._Ca_fec = float(i)





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
        data = load_data(path=path, sheet='Toilet')
        for para in data.index:
            if para in ('desiccant_V', 'desiccant_rho'):
                setattr(self, para, data.loc[para]['expected'])
            else:
                setattr(self, '_'+para, data.loc[para]['expected'])
        del data
        
        self._empty_ratio = 0.59
        
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
            # raise ValueError('if_ideal_emptying is True, empty_ratio should be 1,'
            #                  f'the set value {i} is ignored.')
        self._empty_ratio = float(i)

    @property
    def COD_removal(self):
        '''[float] Fraction of COD removed during storage.'''
        return self._COD_removal
    @COD_removal.setter
    def COD_removal(self, i):
        self._COD_removal = float(i)

    @property
    def N_removal(self):
        '''[float] Fraction of N removed through denitrification during storage.'''
        return self._N_removal
    @N_removal.setter
    def N_removal(self, i):
        self._N_removal = float(i)

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
        '''[float] Fraction of N emitted as N2O due to inappropriate emptying.'''
        self._N2O_EF_aq = i




class PitLatrine(Toilet):
    '''Single pit latrine.'''

    def __init__(self, ID='', ins=None, outs=(), N_user=1, life_time=8,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=False, if_air_emission=True, if_ideal_emptying=True, 
                 OPEX_over_CAPEX=0.05,
                 if_infiltration=True, if_shared=True,
                 if_pit_above_water_table=True, **kwargs):

        '''

        if_infiltration : [bool]
            If infiltration to soil occurs
            (i.e., if the pit walls and floors are permeable).
        if_pit_above_water_table : [bool]
            If the pit is above local water table.
            
        '''

        Toilet.__init__(self, ID, ins, outs, N_user, life_time,
                        if_toilet_paper, if_flushing, if_cleansing, if_desiccant,
                        if_air_emission, if_ideal_emptying, OPEX_over_CAPEX)
    
        self.if_infiltration = if_infiltration
        self.if_pit_above_water_table = if_pit_above_water_table
        self.if_shared = if_shared
        data = load_data(path=path, sheet='PitLatrine')
        for para in data.index:
            if para in ('MCF', 'N2O_EF'):
                setattr(self, '_'+para, eval(data.loc[para]['expected']))
            setattr(self, '_'+para, data.loc[para]['expected'])
        del data
        self._pit_depth = 4.57 # m
        self._pit_area = 0.8 # m2
        for attr, value in kwargs:
            setattr(self, attr, value)
        
    __init__.__doc__ = __doc__ + Toilet.__init__.__doc__ + __init__.__doc__
    __doc__ = __init__.__doc__

    _N_outs = 1

    def _run(self):
        Toilet._run(self)

        N_leaching = int(self.if_infiltration)*self.N_leaching
        

    def _design(self):
        design = self.design_results
        design['Cement'] = 700
        design['Sand'] = 2.2 * 1442
        design['Gravel'] = 0.8 * 1600
        design['Bricks'] = 54 * 0.0024 * 1750
        design['Plastic'] = 16 * 0.63
        design['Steel'] = 0.00425  * 7900
        design['Wood'] = 0.19
        design['Excavation'] = self.pit_depth * self.pit_area
        
    def _cost(self):
        self.purchase_costs['Toilet'] = 449
        self._OPEX = self.purchase_costs['Toilet']*self.OPEX_over_CAPEX/365/24

    @property
    def pit_depth(self):
        '''[float] Depth of the pit, [m].'''
        return self._pit_depth
    @pit_depth.setter
    def pit_depth(self, i):
        self._pit_depth = float(i)
        
    @property
    def pit_area(self):
        '''[float] Area of the pit, [m2].'''
        return self._pit_area
    @pit_area.setter
    def pit_area(self, i):
        self._pit_area = float(i)

    @property
    def emptying_period(self):
        '''[float] Time interval between pit emptying, [yr].'''
        return self._emptying_period
    @emptying_period.setter
    def emptying_period(self, i):
        self._emptying_period = float(i)

    @property
    def sludge_accum_rate(self):
        '''[float] Sludge accumulation rate, [L/cap/yr].'''
        return self._sludge_accum_rate
    @sludge_accum_rate.setter
    def sludge_accum_rate(self, i):
        self._sludge_accum_rate = float(i)

    @property
    def N_leaching(self):
        '''
        [float] Fraction of input N that leaches to the soil
        (if if_infiltration is True).
        '''
        return self._N_leaching
    @N_leaching.setter
    def N_leaching(self, i):
        self._N_leaching = float(i)

    @property
    def P_leaching(self):
        '''
        [float] Fraction of input P that leaches to the soil
        (if if_infiltration is True).
        '''
        return self._P_leaching
    @P_leaching.setter
    def P_leaching(self, i):
        self._P_leaching = float(i)

    @property
    def K_leaching(self):
        '''
        [float] Fraction of input K that leaches to the soil
        (if if_infiltration is True).
        '''
        return self._K_leaching
    @K_leaching.setter
    def K_leaching(self, i):
        self._K_leaching = float(i)

    def _return_EF_num(self):
        # self._MCF and self._N2O_EF are tuples for
        # single_above_water, communal_above_water, below_water
        if not self.if_pit_above_water_table:
            return 2
        elif self.if_shared:
            return 1
        else:
            return 0

    @property
    def MCF_loss(self):
        '''[float] Methane correction factor for COD degraded during storage.'''
        return self._MCF_loss[self._return_EF_num()]
    @MCF_loss.setter
    def MCF_loss(self, i):
        self._MCF_loss[self._return_EF_num()]= float(i)

    @property
    def N2O_EF_loss(self):
        '''[float] Fraction of N emitted as N2O during storage.'''
        return self._N2O_EF_loss[self._return_EF_num()]
    @N2O_EF_loss.setter
    def N2O_EF_loss(self, i):
        self._N2O_EF_loss[self._return_EF_num()]= float(i)



class UDDT(Toilet):
    '''
    Urine-diverting dry toilet with liquid storage tank and dehydration vault
    for urine and feces storage, respectively.
    '''
    
    def __init__(self, ID='', ins=None, outs=(), N_user=1, life_time=8,
                 if_toilet_paper=True, if_flushing=True, if_cleansing=False,
                 if_desiccant=False, if_air_emission=True, if_ideal_emptying=True,
                 OPEX_over_CAPEX=0.1,
                 T=273.15+24, safety_factor=1, if_prep_loss=True, if_treatment=False,
                 **kwargs):

        '''

        T : [float]
            Temperature, [K].
        safety_factor : [float]
            Safety factor for pathogen removal during onsite treatment,
            must be larger than 1.            
        if_treatment : [bool]
            If has onsite treatment.
        if_pit_above_water_table : [bool]
            If the pit is above local water table.
            
        '''

        Toilet.__init__(self, ID, ins, outs, N_user, life_time,
                        if_toilet_paper, if_flushing, if_cleansing, if_desiccant,
                        if_air_emission, if_ideal_emptying, OPEX_over_CAPEX)
    
        self.T = T
        self._safety_factor = safety_factor
        self.if_prep_loss = if_prep_loss
        self.if_treatment = if_treatment
        data = load_data(path=path, sheet='UDDT')
        for para in data.index:
            setattr(self, '_'+para, data.loc[para]['expected'])
        del data
        self._tank_V = 60/1e3 # m3
        for attr, value in kwargs:
            setattr(self, attr, value)
    
    __init__.__doc__ = __doc__ + Toilet.__init__.__doc__ + __init__.__doc__
    __doc__ = __init__.__doc__
    
    _N_outs = 2

    def _run(self):
        ur, fec = self.outs
        Toilet._run(self)
        self.outs[0].mix_from(self.ins)

    def _design(self):
        design = self.design_results
        design['Cement'] = 200
        design['Sand'] = 0.6 * 1442
        design['Gravel'] = 0.2 * 1600
        design['Bricks'] = 682 * 0.0024 * 1750
        design['Plastic'] = 4 * 0.63
        design['Steel'] = 0.00351 * 7900
        design['Stainless steel sheet'] = 28.05 * 2.64
        design['Wood'] = 0.222
        
    def _cost(self):
        self.purchase_costs['Toilet'] = 553
        #!!! What is operating hours is different, maybe better to make this in TEA
        self._OPEX = self.purchase_costs['Toilet']*self.OPEX_over_CAPEX/365/24

    @property
    def safety_factor(self):
        return self._safety_factor
    @safety_factor.setter
    def safety_factor(self, i):
        if i < 1:
            raise ValueError(f'safety_factor must be larger than 1, not {i}')
        self._safety_factor = float(i)

    @property
    def collection_period(self):
        '''[float] Time interval between storage tank collection, [d].'''
        return self._collection_period
    @collection_period.setter
    def collection_period(self, i):
        self._collection_period = float(i)

    @property
    def tank_V(self):
        '''[float] Tank volume, [m3].'''
        return self._tank_V
    @tank_V.setter
    def tank_V(self, i):
        self._tank_V = float(i)

    @property
    def struvite_pKsp(self):
        '''[float] Precipitation constant of struvite.'''
        return self._struvite_pKsp
    @struvite_pKsp.setter
    def struvite_pKsp(self, i):
        self._struvite_pKsp = float(i)

    @property
    def prep_sludge(self):
        '''
        [float] Fraction of total precipitate appearing as sludge that can
        settle and be removed.
        '''
        return self._prep_sludge
    @prep_sludge.setter
    def prep_sludge(self, i):
        self._prep_sludge = float(i)

    @property
    def log_removal(self):
        '''Desired level of pathogen inactivation.'''
        return self._log_removal
    @log_removal.setter
    def log_removal(self, i):
        self._log_removal = float(i)

    @property
    def ur_pH(self):
        '''Urine pH.'''
        return self._ur_pH
    @ur_pH.setter
    def ur_pH(self, i):
        self._ur_pH = float(i)

    @property
    def MCF_loss(self):
        '''[float] Methane correction factor for COD degraded during storage.'''
        return self._MCF_loss
    @MCF_loss.setter
    def MCF_loss(self, i):
        self._MCF_loss = float(i)

    @property
    def N2O_EF_loss(self):
        '''[float] Fraction of N emitted as N2O during storage.'''
        return self._N2O_EF_loss
    @N2O_EF_loss.setter
    def N2O_EF_loss(self, i):
        self._N2O_EF_loss= float(i)

    @property
    def moi_min(self):
        '''[float] Minimum moisture content of feces.'''
        return self._moi_min
    @moi_min.setter
    def moi_min(self, i):
        self._moi_min= float(i)

    @property
    def moi_red_rate(self):
        '''[float] Exponential reduction rate of feces moisture.'''
        return self._moi_red_rate
    @moi_red_rate.setter
    def moi_red_rate(self, i):
        self._moi_red_rate= float(i)




























