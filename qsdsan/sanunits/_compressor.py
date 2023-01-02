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
from math import ceil, floor

__all__ = ('IsothermalCompressor',)

class IsothermalCompressor(bst.units.IsothermalCompressor):
    '''
    Similar to biosteam.units.IsothermalCompressor, but can calculate number of units.
    '''
    
    def _design(self):
        super()._design()
        D = self.design_results
        power = D['Ideal power']/D['Driver efficiency']
        D['Number of 300 kW unit'] = floor(power/300)
        D['Number of 4 kW unit'] = 0
        if (power - D['Number of 300 kW unit']*300) <= 60:
            # according to Ecoinvent 3: the impact of at most 15 4 kW unit is smaller than 1 300 kW unit
            # therefore, if the rest of power is smaller than 60 kW, use multiple small units
            # else, add one large unit
            D['Number of 4 kW unit'] = ceil((power - D['Number of 300 kW unit']*300)/4)
        else:
            D['Number of 300 kW unit'] += 1
        
        construction = getattr(self, 'construction', [])
        if construction:
            construction[0].quantity = D['Number of 4 kW unit']
            construction[1].quantity = D['Number of 300 kW unit']
        else:
            self.construction = [
                qs.Construction('compressor_4kW', linked_unit=self, item='Compressor_4kW', quantity_unit='ea', quantity=D['Number of 4 kW unit']),
                qs.Construction('compressor_300kW', linked_unit=self, item='Compressor_300kW', quantity_unit='ea', quantity=D['Number of 4 kW unit'])
                ]