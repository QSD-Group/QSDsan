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
from math import log

__all__ = ('VacuumPump',)

class VacuumPump(Equipment):
    
    def __init__(self, linked_unit=None, ID=None, 
                 units={'P_suction': 'torr', 'Vacuum power': 'kW'},
                 F_BM=1., F_D=1., F_P=1., F_M=1., 
                 lifetime=20, lifetime_unit='yr', 
                 P_suction=100, **kwargs):
        
        Equipment.__init__(self, linked_unit, ID, units, F_BM, F_D, F_P, F_M,
                           lifetime, lifetime_unit, **kwargs)
        self.P_suction = P_suction
    
    @property
    def P_vacuum(self):
        P_break = 21.4 * self._SF()**0.924
        e_motor = 0.8 + 0.0319*log(P_break, 10) - 0.00182 * log(P_break, 10)**2
        return P_break/e_motor # in kW

    def _SF(self, for_cost=False):
        # P_sunction in torr, Qm in kg/hr
        for ws in self.linked_unit.outs: 
            if ws.phase == 'g': 
                Qm = ws.F_mass  # kg/hr
                break
        Ps = self.P_suction  # torr
        if for_cost: 
            return Qm/Ps*2.2  # convert Qm to lb/hr
        else:
            return min(16, max(0.2, Qm/Ps))
    
    def _design(self):
        return {'Vacuum power': self.P_vacuum}
        
    def _cost(self):
        SF = self._SF(True)
        return 1915 * SF**0.41
    
    @property
    def P_suction(self):
        return self._Ps
    
    @P_suction.setter
    def P_suction(self, P):
        if P <= 0: raise ValueError(f'Suction pressure [torr] must be positive, not {P}')
        self._Ps = P