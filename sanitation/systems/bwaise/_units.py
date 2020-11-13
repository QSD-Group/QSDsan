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
import pandas as pd
import thermosteam as tmo
from biosteam.units.decorators import cost
from sanitation import SanUnit

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
    
    def __init__(self, ID='', ins=None, outs=(), load_defaults=True, N_user=1, **kwargs):                
        SanUnit.__init__(self, ID, ins, outs)
        if load_defaults:
            self.load_data()
        self.N_user = N_user
        for attr, value in kwargs:
            setattr(self, attr, value)
            
    def load_data(self, path=None):
        if not path:
            path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'inputs/defaults.xlsx')
            data = pd.read_excel(path, sheet_name='Excretion', index_col=0)
        elif path[-4:] == 'xlsx' or path[-4:] == '.xls':
            try: data = pd.read_excel(path, sheet_name='Excretion', index_col=0)
            except: data = pd.read_excel(path, index_col=0)
        elif path[-4:] == '.csv':
            data = pd.read_csv(path, index_col=0)
        else:
            raise ValueError('Only csv or xlsx files can be loaded.')
        for para in data.index:
            setattr(self, para, data.loc[para]['expected'])
    
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



class PitLatrine(SanUnit):
    
    pit_depth = 4.57 # m
    pit_area = 0.8 # m2
    toilet_paper = 12.4*545/1e6 / 24 # kg/cap/hr
    flushing_water = 10 / 24 # kg/cap/hr
    cleaning_water = 1 / 24 # kg/cap/hr
    
    def __init__(self, ID='', ins=None, outs=(), N_user=1, life_time=8,
                 if_toilet_paper=True, if_flushing=True, if_cleaning=False,
                 OPEX_over_CAPEX=0.05):
        '''
        Pit latrine.

        Parameters
        ----------
        N_user : [float]
            Number of people that use the toilet per hour.
        life_time : [float]
            Life time of the toilet in year.
        if_toilet_paper : [bool]
        if_flushing : [bool]
        if_cleaning : [bool]
        opex_over_capex : [float]
            Fraction of annual operating cost over total capital cost

        '''

        SanUnit.__init__(self, ID, ins, outs)
        self.N_user = N_user
        self.life_time = life_time
        self.if_toilet_paper = if_toilet_paper
        self.if_flushing = if_flushing
        self.if_cleaning = if_cleaning
        self.OPEX_over_CAPEX = OPEX_over_CAPEX
    
    _N_ins = 5
    _N_outs = 1
    _units = {
        'Cement': 'kg',
        'Sand': 'kg',
        'Gravel': 'kg',
        'Bricks': 'kg',
        'Plastic': 'kg',
        'Steel': 'kg',        
        'Wood': 'm3',
        'Excavation': 'm3'
        }
    _BM = {'Toilet': 1}

    def _run(self):
        ur, fecs, tp, fw, cw = self.ins
        mixture, = self.outs
        
        N_user = self.N_user
        tp.imass['Tissue'] = int(self.if_toilet_paper)*self.toilet_paper * N_user
        fw.imass['H2O'] = int(self.if_flushing)*self.flushing_water * N_user
        cw.imass['H2O'] = int(self.if_cleaning)*self.cleaning_water * N_user
        
        mixture.mix_from(self.ins)        

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
        self._OPEX = self.purchase_costs['Toilet'] * self.OPEX_over_CAPEX



















