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

import pytest

def test_component():
    import thermosteam as tmo
    from sanitation import Component, Components
    components = Components.load_default(default_compile=False)
    with pytest.raises(TypeError):
        H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'))
    H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'),
                                  i_C=0, i_N=0, i_P=0, i_K=0, i_mass=1,
                                  i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                                  f_Vmass_Totmass=0,
                                  particle_size='Soluble',
                                  degradability='Undegradable', organic=False)
    with pytest.raises(ValueError):
        components.append(H2O)
    components = Components.load_default()
    tmo.settings.set_thermo(components)
    
# This just means that if pytest runs this module, it calls the test_component function
if __name__ == '__main__':
    test_component()