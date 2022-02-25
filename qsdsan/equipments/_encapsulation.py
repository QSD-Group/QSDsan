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

__all__ = ('Beads',)

class Beads(Equipment):
    
    def __init__(self, linked_unit=None, ID=None, 
                 units={'Diameter': 'mm',
                        'Density': 'kg/m3',
                        'Unit cost': 'USD/kg',
                        'Bead total mass': 'kg'},
                 F_BM=1.15, F_D=1., F_P=1., F_M=1., 
                 lifetime=0.5, lifetime_unit='yr', 
                 d_bead=1.0, rho_bead=265, p_bead=1440,
                 **kwargs):
        
        Equipment.__init__(self, linked_unit, ID, units, F_BM, F_D, F_P, F_M,
                           lifetime, lifetime_unit, **kwargs)
        self.d_bead = d_bead
        self.rho_bead = rho_bead
        self.p_bead = p_bead
    
    def _design(self):
        V_bead = self.linked_unit.design_results['Bead total volume']
        return {'Bead total mass': V_bead * self.rho_bead}
        
    def _cost(self):
        m = self._design()['Bead total mass']
        return m * self.p_bead