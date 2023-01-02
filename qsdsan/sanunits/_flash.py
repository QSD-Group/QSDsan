#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import biosteam as bst, qsdsan as qs

__all__ = ('Flash',)

_lb_to_kg = qs.utils.auom('lb').conversion_factor('kg')

class Flash(bst.units.Flash):
    '''
    Similar to biosteam.units.Flash, but includes construction.
    
    See Also
    --------
    :class:`biosteam.units.Flash`
    '''
    
    def _design(self):
        super()._design()
        D = self.design_results
        
        if self.outs[0].F_mass != 0:
            construction = getattr(self, 'construction', [])
            if construction: construction[0].quantity = D['Weight']*_lb_to_kg
            # N = 1, don't have to * numbers of the reactor
            else:
                self.construction = [
                    qs.Construction('carbon_steel', linked_unit=self, item='Carbon_steel', 
                                    quantity=D['Weight']*_lb_to_kg, quantity_unit='kg'),
                    ]