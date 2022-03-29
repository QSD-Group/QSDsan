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
from qsdsan import SanUnit, Construction
from qsdsan.utils.loading import load_data, data_path

__all__ = ('SludgePasteurization',)

import os
data_path = os.path.join(data_path, 'sanunit_data/_sludge_pasteurization.tsv')


#To scale the system from 0 - 125 users based on sludge pasteurization units corresponding to 25 users: 
ppl = 120

if (ppl <=25): P = 1/4
elif (ppl >25 and ppl <=50): P = 2/4
elif (ppl >50 and ppl <=75): P = 3/4
elif (ppl >75 and ppl <=100): P = 4/4
else: P = 5/4

class SludgePasteurization(SanUnit):
    '''
    SLUDGE PASTEURIZATION FOR SLUDGE FROM SEPTIC TANK
    Parameters
    ----------
    if_combustion : bool
         If include combusion reaction during simulation.
    temp_pasteurization : float
        Pasteurization temperature (Kelvin)
    Examples
        temp_pasteurization : float
        Pasteurization temperature is 70C or 343.15 Kelvin
    sludge_temp : float
        Temperature of sludge is 10C or 283.15K
    target_MC : float
        Target moisture content is 10%
    heat_loss : float
        Heat loss during pasteurization process is assumed to be 10%
    lhv_lpg : float
        Lower heating value of Liquid Petroleum Gas at 25C/298.15K is 46-51 MJ/kg from World Nuclear Organization
        The lower heating value (also known gross calorific value or gross energy) of a fuel is defined as the amount of heat released 
        by a specified quantity (initially at 25°C) once it is combusted and the products remain evaporated and into atmosphere.
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', 
                 heat_loss=0.1, target_MC = 0.1, sludge_temp = 10 + 273.15, 
                 temp_pasteurization= 70 + 273.15, lhv_lpg = 48.5, **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=1)
        self._heat_loss = heat_loss  
        self._target_MC = target_MC      
        self._sludge_temp = sludge_temp
        self._temp_pasteurization = temp_pasteurization
        self.lhv_lpg = lhv_lpg  
        self.price_ratio = 1
        
        data = load_data(path=data_path)
        
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    _N_ins = 3
    _N_outs = 1

    def _run(self):
        air, sludge, lpg = self.ins
        lpg.phase = 'l'
        sludge.phase = 's'
        treated_sludge = self.outs[0] 
        treated_sludge.copy_like(sludge) 
        
        
        #Constants
        # Specific Heat capacity of water 
        self.Cp_w = 4.184 #kJ kg^-1 C^-1 
        # Specific Heat capacity of dry matter (sludge)
        self.Cp_dm = 1.231 #kJ kg^-1 C^-1
        # Specific latent heat of vaporization of water
        self.l_w = 2260 # kJ kg^-1
        
        # Mass calculations
        # total amount of water in sludge
        self.M_w = sludge.imass['H2O'] #kg/hrs
        # total amount of dry matter in sludge
        self.M_dm = (sludge.F_mass - sludge.imass['H2O']) #kg/hr
        # total amount of water to be evaporated to reach target dry matter content
        self.M_we = (sludge.imass['H2O'] - sludge.F_mass * self.target_MC) #kg/hr
        
        # Overall heat required for pasteurization
        self.Q_d = (self.M_w * self.Cp_w * (self.temp_pasteurization - self.sludge_temp) + self.M_dm * self.Cp_dm * (self.temp_pasteurization - self.sludge_temp) + self.M_we * self.l_w)*(1 + self.heat_loss) #kJ/hr
        
        self.lpg_vol_reqd = ((self.Q_d/1000) / self.lhv_lpg)# kg/hr 
        lpg.imass['LPG'] = self.lpg_vol_reqd
           
        treated_sludge.imass['H2O'] = sludge.F_mass * self.target_MC
        treated_sludge.imass['OtherSS'] = sludge.F_mass - self.M_we        
        treated_sludge.imass['P'] = sludge.imass['P']     
        treated_sludge.imass['N'] = sludge.imass['N']
        
    def _design(self):
        design = self.design_results
        design['Steel'] = S_quant = self.sludge_dryer_weight + self.sludge_barrel_weight      
        self.construction = (
        Construction(item='Steel', quantity = S_quant, quantity_unit = 'kg'),
            )
        self.add_construction(add_cost=False)
        
    def _cost(self):
        C= self.baseline_purchase_costs
        C['Dryer'] = 0 
        C['Barrel'] = 0 
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.add_OPEX =  self._calc_replacement_cost() + self._calc_labor_cost() # USD/hr (all items are per hour)
        self.power_demand = 0
        self.power_utility(self.power_demand)

    def _calc_replacement_cost(self):
        sludge_replacement_cost = (self.sludge_dryer + self.sludge_barrel) / 10 * P #USD/yr
        return sludge_replacement_cost/ (20 * 365 * 24) * self.price_ratio# USD/hr (all items are per hour)


    def _calc_labor_cost(self):
        labor_cost = (self.wages * self.sludge_labor_maintenance) * P
        return labor_cost/(365 * 24)  # USD/hr (all items are per hour)

    @property
    def heat_loss(self):
        return self._heat_loss
    @heat_loss.setter
    def heat_loss(self, i):
        self._heat_loss = i

    @property
    def target_MC(self):
        return self._target_MC
    @target_MC.setter
    def target_MC(self, i):
        self._target_MC = i
        
    @property
    def sludge_temp(self):
        return self._sludge_temp
    @sludge_temp.setter
    def sludge_temp(self, i):
        self._sludge_temp = i        
        
    @property
    def temp_pasteurization(self):
        return self._temp_pasteurization
    @temp_pasteurization.setter
    def temp_pasteurization(self, i):
        self._temp_pasteurization = i
        
        
        
        
        
        