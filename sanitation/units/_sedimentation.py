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
from ._decay import Decay
from ..utils.loading import load_data, data_path

__all__ = ('Sedimentation',)

data_path += 'unit_data/Sedimentation.csv'


class Sedimentation(SanUnit, Decay):
    '''Sedimentation of wastes into liquid and solid phases.'''
    
    def __init__(self, ID='', ins=None, outs=(), if_N_degradation=True,
                 **kwargs):
        
        '''

        Parameters
        ----------
        ins : WasteStream
            Waste for treatment.
        outs : WasteStream
            Drained liquid, settled solids, fugitive CH4, and fugitive N2O.
        if_N_degradation : [bool]
            If N degradation and N2O emission occur during treatment.

        '''        
        
        SanUnit.__init__(self, ID, ins, outs)
        self.if_N_degradation = if_N_degradation

        data = load_data(path=data_path)
        for para in data.index:
            if para == 'split':
                value = eval(data.loc[para]['expected'])
                setattr(self, para, value)
            else:
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
        split = self.split
        if self._split_type == 'float':
            liq.copy_like(waste)
            sol.copy_like(waste)
            sol.mass *= self.split
            liq.mass -= sol.mass
        else:
            for var in self.split.keys():
                #!!! In the future this should be best by changing the state variable
                if var == 'TS':
                    sol.imass['OtherSS'] = split[var] * waste.imass['OtherSS']
                elif var == 'COD':
                    sol._COD = split[var] * waste._COD
                    liq._COD = waste._COD - sol._COD
                elif var == 'N':
                    N_sol = split[var]*(waste.imass['NH3']+waste.imass['NonNH3'])
                    NonNH3_rmd, NH3_rmd = \
                        self.allocate_N_removal(N_sol, waste.imass['NonNH3'])
                    sol.imass ['NonNH3'] = NonNH3_rmd
                    sol.imass ['NH3'] = NH3_rmd
                else:
                    sol.imass[var] = split[var] * waste.imass[var]
            liq.mass = waste.mass - sol.mass
            
        # COD degradation in settled solids
        COD_loss = self.first_order_decay(k=self.decay_k_COD,
                                          t=self.tau/365,
                                          max_removal=self.COD_max_removal)

        sol._COD *= 1 - COD_loss
        sol.imass['OtherSS'] *= 1 - COD_loss
        sol.imass['H2O'] = waste.F_mass * self.setted_frac - sol.F_mass
        if sol.imass['H2O'] < 0:
            sol.imass['H2O'] = 0
            msg = 'Negative water content calcualted for settled solids' \
                'try smaller split or larger settled_frac.'
            warn(msg, source=self)
        liq.imass['H2O'] = waste.imass['H2O'] - sol.imass['H2O']
        CH4.imass['CH4'] = sol.COD/1e3*sol.F_vol*COD_loss * \
            self.max_CH4_emission*self.MCF_decay # COD in mg/L (g/m3)

        # N degradation
        if self.if_N_degradation:
            N_loss = self.first_order_decay(k=self.decay_k_N,
                                            t=self.tau/365,
                                            max_removal=self.N_max_removal)
            N_loss_tot = N_loss*sol.TN/1e3*sol.F_vol
            NH3_rmd, NonNH3_rmd = \
                self.allocate_N_removal(N_loss_tot, sol.imass['NH3'])
            sol.imass ['NH3'] -=  NH3_rmd
            sol.imass['NonNH3'] -= NonNH3_rmd

        N2O.imass['N2O'] = N_loss_tot*self.N2O_EF_decay*44/28
    
    
    def _design(self):
        design = self.design_results
        #!!! Why isn't tau used?
        design['Tank number'] = N = self.N_tank
        design['Single tank volume'] = V_single = self.tank_V
        L2W = self.tank_L_to_W
        W2H = self.tank_W_to_H
        design['Tank height'] = H = (V_single/(L2W*(W2H**2)))**(1/3)
        design['Tank width'] = W = H * W2H
        design['Tank length'] = L = W * L2W
        # Concrete
        thick = self.tank_thickness
        side_concrete = N*thick*(L*W+2*W*H+2*L*H)
        column_concrete = N*(thick**2)*H*self.column_per_side*2
        design['Concrete volume'] = side_concrete + column_concrete
        # Steel
        design['Roof area'] = N*L*W/(np.cos(self.roof_slope/180*np.pi))
        side_area = N*2*(L*H + W*H)
        design['Steel'] = design['Roof area'] + side_area
        
    
    def _cost(self):
        pass
    
    

    @property
    def tau(self):
        '''[float] Residence time, [d].'''
        return self._tau
    @tau.setter
    def tau(self, i):
        self._tau = float(i)

    @property
    def split(self):
        '''
        [float] or [dict] Fractions of material retention in the settled solids
        before degradation. If a single number is provided, then it is assumed
        that retentions of all Components in the WasteStream are the same.
        
        Note
        ----
            Set state variable values (e.g., COD) will be retained if the retention
            ratio is a single number (treated like the loss stream is split
            from the original stream), but not when the ratio is a dict.

        '''
        return self._split
    @split.setter
    def split(self, i):
        try:
            self._split = float(i)
            self._split_type = 'float'
        except:
            if isinstance(i, dict):
                self._split = i
                self._split_type = 'dict'
            else:
                raise TypeError(f'Only float or dict allowed, not {type(i).__name__}.')

    @property
    def setted_frac(self):
        '''[float] Fraction of influent that settles as solids.'''
        return self._setted_frac
    @setted_frac.setter
    def setted_frac(self, i):
        self._setted_frac = float(i)

    @property
    def tank_V(self):
        '''[float] Volume of the sedimentation tank.'''
        return self._tank_V
    @tank_V.setter
    def tank_V(self, i):
        self._tank_V = float(i)

    @property
    def tank_L_to_W(self):
        '''[float] Length-to-width ratio of the sedimentation tank.'''
        return self._tank_L_to_W
    @tank_L_to_W.setter
    def tank_L_to_W(self, i):
        self._tank_L_to_W = float(i)

    @property
    def tank_W_to_H(self):
        '''[float] Width-to-height ratio of the sedimentation tank.'''
        return self._tank_W_to_H
    @tank_W_to_H.setter
    def tank_W_to_H(self, i):
        self._tank_W_to_H = float(i)

    @property
    def N_tank(self):
        '''[int] Number of sedimentation tanks, float will be converted to the smallest integer.'''
        return self._N_tank
    @N_tank.setter
    def N_tank(self, i):
        self._N_tank = np.ceil(i)

    @property
    def column_per_side(self):
        '''[int] Number of columns per side of sedimentation tanks, float will be converted to the smallest integer.'''
        return self._column_per_side
    @column_per_side.setter
    def column_per_side(self, i):
        self._column_per_side = np.ceil(i)

    @property
    def tank_thickness(self):
        '''[float] Wall thickness of the concrete tank.'''
        return self._tank_thickness
    @tank_thickness.setter
    def tank_thickness(self, i):
        self._tank_thickness = float(i)

    @property
    def roof_slope(self):
        '''[float] Slope of the tank roof, [°].'''
        return self._roof_slope
    @roof_slope.setter
    def roof_slope(self, i):
        self._roof_slope = float(i)

    @property
    def roof_mass(self):
        '''[float] Mass of the tank roof, [kg/m2].'''
        return self._roof_mass
    @roof_mass.setter
    def roof_mass(self, i):
        self._roof_mass = float(i)








