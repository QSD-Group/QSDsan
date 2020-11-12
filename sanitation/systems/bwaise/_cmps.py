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

import biosteam as bst
from sanitation import Component, Components

__all__ = ('cmps', )

kcal = Component('kcal', HHV=4184, # 1 kcal is 4184 J
                 phase='l', i_C=0, i_N=0, i_P=0, i_K=0,
                 i_mass=1, i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                 f_Vmass_Totmass=0, particle_size='Soluble',
                 degradability='Biological', organic=True)


NH3 = Component.from_chemical('NH3', bst.Chemical('NH3'), measured_as='N',
                              phase='l', i_C=0, i_N=1, i_P=0, i_K=0,
                              i_mass=1, i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                              f_Vmass_Totmass=0, particle_size='Soluble',
                              degradability='Undegradable', organic=False)

nonNH3 = Component.from_chemical('nonNH3', bst.Chemical('N'), measured_as='N',
                                 phase='l', i_C=0, i_N=1, i_P=0, i_K=0,
                                 i_mass=1, i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                                 f_Vmass_Totmass=0, particle_size='Soluble',
                                 degradability='Undegradable', organic=False,
                                 description='Non-NH3 nitrogen')

P = Component.from_chemical('P', bst.Chemical('P'),
                            phase='l', i_C=0, i_N=0, i_P=1, i_K=0,
                            i_mass=1, i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                            f_Vmass_Totmass=0, particle_size='Soluble',
                            degradability='Undegradable', organic=False)

K = Component.from_chemical('K', bst.Chemical('K'),
                            phase='l', i_C=0, i_N=0, i_P=0, i_K=1,
                            i_mass=1, i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                            f_Vmass_Totmass=0, particle_size='Soluble',
                            degradability='Undegradable', organic=False)

Mg = Component.from_chemical('Mg', bst.Chemical('Mg'),
                             phase='l', i_C=0, i_N=0, i_P=0, i_K=0,
                             i_mass=1, i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                             f_Vmass_Totmass=0, particle_size='Soluble',
                             degradability='Undegradable', organic=False)

Ca = Component.from_chemical('Ca', bst.Chemical('Ca'),
                             phase='l', i_C=0, i_N=0, i_P=0, i_K=0,
                             i_mass=1, i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                             f_Vmass_Totmass=0, particle_size='Soluble',
                             degradability='Undegradable', organic=False)

H2O = Component.from_chemical('H2O', bst.Chemical('H2O'),
                              phase='l', i_C=0, i_N=0, i_P=0, i_K=0,
                              i_mass=1, i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                              f_Vmass_Totmass=0, particle_size='Soluble',
                              degradability='Undegradable', organic=False)


for cmp in (kcal, nonNH3, P, K, Mg, Ca):
    cmp.default()
    cmp.copy_models_from(H2O, ('sigma', 'epsilon', 'kappa', 'V', 'Cn', 'mu'))

others = kcal.copy('others', HHV=0, degradability='Undegradable', organic=False,
                    description='Fillter Component for unspecified mass')

cmps = Components((kcal, NH3, nonNH3, P, K, Mg, Ca, H2O, others))
cmps.compile()






