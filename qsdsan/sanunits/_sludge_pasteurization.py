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

from .. import SanUnit

__all__ = ('SludgePasteurization',)


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
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_
    '''
    
    #TELL SHION ABOUT ADDING 273.15, THE TWO MISTAKES CANCEL EACH OTHER ABOUT BUT FOR CONSISTENCY AND ACCURACY YOU PROBABLY WANT TO FIX IT 
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', 
                 heat_loss=0.1, target_MC = 0.1, sludge_temp = 10 + 273.15, 
                 temp_pasteurization= 70 + 273.15, lhv_lpg = 50125):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=1)
        self._heat_loss = heat_loss  
        self._target_MC = target_MC      
        self._sludge_temp = sludge_temp
        self._temp_pasteurization = temp_pasteurization
        
    _N_ins = 3
    _N_outs = 4

    def _run(self):
        air, sludge, lpg = self.ins
        lpg.phase = 'l'
        sludge.phase = 's'
        struvite, treated_sludge = self.outs 
        treated_sludge.copy_like(sludge) 
        treated_sludge = treated_sludge.mass
        
        #Constants
        # Specific Heat capacity of water 
        self.Cp_w = 4.184 #kJ kg^-1 C^-1 
        # Specific Heat capacity of dry matter (sludge)
        self.Cp_dm = 1.231 #kJ kg^-1 C^-1
        # Specific latent heat of vaporization of water
        self.l_w = 2260 # kJ kg^-1
        
        #
        
        #DON'T DIVIDE BY 1000 IF YOU WANT KG/HR 
        # Mass calculations
        # total amount of water in sludge
        self.M_w = sludge.imass['H2O'] #kg/hr
        # total amount of dry matter in sludge
        self.M_dm = (sludge.F_mass - sludge.imass['H2O']) #kg/hr
        # total amount of water to be evaporated to reach target dry matter content
        self.M_we = (sludge.imass['H2O'] - sludge.F_mass * self.target_MC) #kg/hr
        
        #WHAT IS THE TEMP OF EVAPORATED WATER, WHERE DOES THE EVAPORATED WATER GO, RECOVER HEAT, 
        #AND FINAL TEMP OF SLUDGE, AT WHAT TEMP IS WATER REMOVED. 
        
        #IF IT SHOULD BE THAT WET CONTENT IS 10% AFTER SLUDGE PASTERUZIATION 
        #TELL SHION THAT HEAT LOSS IS IN .1% IF SHE DIVIDES BY 100
        # Overall heat required for pasteurization
        self.Q_d = ((self.M_w * self.Cp_w * (self.temp_pasteurization - self.sludge_temp) 
                    + self.M_dm * self.Cp_dm * (self.temp_pasteurization - self.sludge_temp) 
                    + self.M_we * self.l_w)*(1 + self.heat_loss)) #kJ/hr
        

        lpg_vol_reqd = self.Q_d / self.lhv_lpg # kg/hr 
        lpg.imass['CH4'] = lpg_vol_reqd
           
        treated_sludge.imass['H2O'] = sludge.F_mass * self.target_MC
        treated_sludge.imass['OtherSS'] = sludge.F_mass - self.M_we        
        treated_sludge.imass['P'] = sludge.imass['P']     
        treated_sludge.imass['N'] = sludge.imass['N']
        

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
            
        
        
        