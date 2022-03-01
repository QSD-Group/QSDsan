# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

Part of this module is based on the BioSTEAM package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from .. import Equipment
from math import ceil

__all__ = ('VerticalMixer',)

class VerticalMixer(Equipment):
    
    def __init__(self, linked_unit=None, ID=None, 
                 units={
                     'power': 'kW',
                     'mixing_intensity': 's^(-1)',
                     'unit_price': 'USD',
                     'Required mixing power': 'kW',
                     'N mixer': ''
                     },
                 F_BM=1., F_D=1., F_P=1., F_M=1., 
                 lifetime=20, lifetime_unit='yr', 
                 power=3.7, mixing_intensity=250, 
                 unit_price=10200, **kwargs):
        
        Equipment.__init__(self, linked_unit, ID, units, F_BM, F_D, F_P, F_M,
                           lifetime, lifetime_unit, **kwargs)
        self.power = power
        self.mixing_intensity = mixing_intensity
        self.unit_price = unit_price
    
    def _design(self):
        return {'Required mixing power': self.P_mix,
                'N mixer': self.N_mix}
        
    def _cost(self):
        return self.unit_price * self.N_mix
    
    @property
    def power(self):
        return self._P
    @power.setter
    def power(self, P):
        if P <= 0: raise ValueError(f'Mixer power [kW] must be positive, not {P}')
        self._P = P
    
    @property
    def mixing_intensity(self):
        return self._G
    @mixing_intensity.setter
    def mixing_intensity(self, g):
        if g <= 0: raise ValueError(f'Mixing intensity [s^-1] must be positive, not {g}')
        self._G = g
    
    @property
    def unit_price(self):
        return self._price
    @unit_price.setter
    def unit_price(self, p):
        self._price = p
    
    @property
    def mu(self):
        for ws in self.linked_unit.outs: 
            if ws.phase == 'l': 
                mu = ws.mu  # kg/hr
                break
        return mu
    
    @property
    def V(self):
        return self.linked_unit.design_results['Reactor volume']
    
    @property
    def P_mix(self):
        return self._G ** 2 * self.mu * self.V / 1000 # kW
    
    @property
    def N_mix(self):
        return ceil(self.P_mix / self.power)