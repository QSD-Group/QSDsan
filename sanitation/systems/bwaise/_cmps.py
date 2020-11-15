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
'''

import thermosteam as tmo
from sanitation import Component, Components

__all__ = ('cmps', )

NH3 = Component.from_chemical('NH3', tmo.Chemical('NH3'), measured_as='N',
                              phase='l', i_N=1, i_mass=1, particle_size='Soluble',
                              degradability='Undegradable', organic=False)

NonNH3 = Component.from_chemical('NonNH3', tmo.Chemical('N'), measured_as='N',
                                 phase='l', i_N=1, particle_size='Soluble',
                                 degradability='Undegradable', organic=False,
                                 description='Non-NH3 nitrogen')

P = Component.from_chemical('P', tmo.Chemical('P'),
                            phase='l', i_P=1, particle_size='Soluble',
                            degradability='Undegradable', organic=False)

K = Component.from_chemical('K', tmo.Chemical('K'),
                            phase='l', i_K=1, particle_size='Soluble',
                            degradability='Undegradable', organic=False)

Mg = Component.from_chemical('Mg', tmo.Chemical('Mg'),
                             phase='l', i_Mg=1, particle_size='Soluble',
                             degradability='Undegradable', organic=False)

Ca = Component.from_chemical('Ca', tmo.Chemical('Ca'),
                             phase='l', i_Ca=0, particle_size='Soluble',
                             degradability='Undegradable', organic=False)

H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'),
                              phase='l', particle_size='Soluble',
                              degradability='Undegradable', organic=False)

OtherSS = Component('OtherSS', phase='l', particle_size='Soluble',
                    degradability='Undegradable', organic=False,
                    description='Unspecified soluble solids')

N2O = Component.from_chemical('N2O', tmo.Chemical('N2O'),
                              phase='g', i_N=28/44, particle_size='Dissolved gas',
                              degradability='Undegradable', organic=False)

CH4 = Component.from_chemical('CH4', tmo.Chemical('CH4'),
                              phase='g', i_C=12/16, particle_size='Dissolved gas',
                              degradability='Biological', organic=True)

for cmp in (NonNH3, P, K, Mg, Ca, OtherSS):
    cmp.default()
    cmp.copy_models_from(H2O, ('sigma', 'epsilon', 'kappa', 'V', 'Cn', 'mu'))

Tissue = Component('Tissue', MW=1, phase='s', particle_size='Particulate',
                    degradability='Undegradable', organic=False,
                    description='Tissue for toilet paper')
# 375 kg/m3 is the average of 250-500 for tissue from
# https://paperonweb.com/density.htm (accessed 2020-11-12)
V_model = tmo.functional.rho_to_V(375, Tissue.MW)
Tissue.V.add_model(V_model)

WoodAsh = Component('WoodAsh', MW=1, phase='s', i_Mg=0.0224, i_Ca=0.3034,
                    particle_size='Particulate', degradability='Undegradable',
                    organic=False, description='Wood ash for desiccant')
V_model = tmo.functional.rho_to_V(760, WoodAsh.MW)
WoodAsh.V.add_model(V_model)

for i in (Tissue, WoodAsh):
    i.copy_models_from(tmo.Chemical('Glucose'), ('Cn', 'mu'))
    
    

cmps = Components((NH3, NonNH3, P, K, Mg, Ca, H2O, OtherSS, N2O, CH4,
                   Tissue, WoodAsh))
cmps.compile()

cmps.set_synonym('H2O', 'Water')




