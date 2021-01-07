#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

import pytest

def test_component():
    import thermosteam as tmo
    from qsdsan import Component, Components
    components = Components.load_default(default_compile=False)
    with pytest.raises(AssertionError):
        H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'))
    H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'),
                                  particle_size='Soluble',
                                  degradability='Undegradable', organic=False)
    with pytest.raises(ValueError):
        components.append(H2O)
    components = Components.load_default()
    tmo.settings.set_thermo(components)
    
# This just means that if pytest runs this module, it calls the test_component function
if __name__ == '__main__':
    test_component()