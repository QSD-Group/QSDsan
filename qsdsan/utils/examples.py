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

from .. import Component, Components

__all__ = ('load_example_cmps',)


# %%

# =============================================================================
# Examplary components to be used in documentations
# =============================================================================

def load_example_cmps():
    '''Load some default components for documentation purpose.'''

    H2O = Component('H2O', search_ID='H2O', particle_size='Soluble',
                     degradability='Undegradable', organic=False)
    
    CO2 = Component('CO2', search_ID='CO2', particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)  

    N2O = Component('N2O', search_ID='N2O', particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)

    NaCl = Component('NaCl', search_ID='NaCl', particle_size='Soluble',
                     degradability='Undegradable', organic=False)

    H2SO4 = Component('H2SO4', search_ID='H2SO4', particle_size='Soluble',
                      degradability='Undegradable', organic=False)
    
    CH4 = Component('CH4', search_ID='CH4', particle_size='Dissolved gas',
                     degradability='Readily', organic=True)

    Methanol = Component('Methanol', search_ID='Methanol',
                         particle_size='Soluble',
                         degradability='Readily', organic=True)

    Ethanol = Component('Ethanol', search_ID='Ethanol',
                         particle_size='Soluble',
                         degradability='Readily', organic=True)
    
    cmps = Components((H2O, CO2, N2O, NaCl, H2SO4, CH4, Methanol, Ethanol))

    cmps.compile()
    cmps.set_synonym('H2O', 'Water')
    cmps.set_synonym('CH4', 'Methane')
    
    return cmps