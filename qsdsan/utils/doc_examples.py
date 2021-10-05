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

from .. import (
    Component, Components, SanStream, System,
    set_thermo,
    sanunits as su
    )

__all__ = ('load_example_cmps', 'load_example_sys',)


# %%

# =============================================================================
# Examplary components
# =============================================================================

def load_example_cmps():
    '''Load some components for documentation purpose.'''

    H2O = Component('H2O', search_ID='H2O', particle_size='Soluble',
                     degradability='Undegradable', organic=False)
    
    CO2 = Component('CO2', search_ID='CO2', phase='g',
                    particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)  

    N2O = Component('N2O', search_ID='N2O', phase='g',
                    particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)

    NaCl = Component('NaCl', search_ID='NaCl', phase='s', particle_size='Soluble',
                     degradability='Undegradable', organic=False)

    H2SO4 = Component('H2SO4', search_ID='H2SO4', phase='s', particle_size='Soluble',
                      degradability='Undegradable', organic=False)
    
    CH4 = Component('CH4', search_ID='CH4', phase='g', particle_size='Dissolved gas',
                     degradability='Readily', organic=True)

    Methanol = Component('Methanol', search_ID='Methanol', phase='l',
                         particle_size='Soluble',
                         degradability='Readily', organic=True)

    Ethanol = Component('Ethanol', search_ID='Ethanol', phase='l',
                         particle_size='Soluble',
                         degradability='Readily', organic=True)
    
    cmps = Components((H2O, CO2, N2O, NaCl, H2SO4, CH4, Methanol, Ethanol))
    for cmp in cmps:
        cmp.default()

    cmps.compile()
    cmps.set_synonym('H2O', 'Water')
    cmps.set_synonym('CH4', 'Methane')
    
    return cmps

# %%

# =============================================================================
# Examplary systems
# =============================================================================

def load_example_sys(cmps=None):
    '''Load a pre-constructed system for documentation purpose.'''
    
    if cmps:
        set_thermo(cmps)
    
    salt_water = SanStream('salt_water', Water=2000, NaCl=50, units='kg/hr')
    methanol = SanStream('methanol', Methanol=20, units='kg/hr')
    ethanol = SanStream('ethanol', Ethanol=10, units='kg/hr')

    M1 = su.MixTank('M1', ins=(salt_water, 'recycled_brine', methanol, ethanol))
    P1 = su.Pump('P1', ins=M1-0)
    H1 = su.HXutility('H1', ins=P1-0, T=350)
    S1 = su.ComponentSplitter('S1', ins=H1-0, split_keys=('Methanol', 'Ethanol'))
    M2 = su.Mixer('M2', ins=(S1-0, S1-1), outs='alcohols')
    S2 = su.Splitter('S2', ins=S1-2, outs=(1-M1, 'waste_brine'), split=0.2)
    sys = System('sys', path=(M1, P1, H1, S1, M2, S2))
    
    return sys
    













