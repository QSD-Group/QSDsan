#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 09:20:36 2020

@author: yalinli_cabbi
"""

import pytest

def test_component():
    import thermosteam as tmo
    from sanitation import Component, utils
    components = utils.load_components_from_excel()
    with pytest.raises(TypeError):
        H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'))
    H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'),
                                  i_C=0, i_N=0, i_P=0, i_K=0, i_mass=1,
                                  i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                                  f_Vmass_Totmass=0,
                                  particle_size='Soluble',
                                  degradability='Undegradable', organic=False)

    components.append(H2O)
    for i in components:
        i.default()
        i.copy_models_from(components.H2O, ['sigma', 'epsilon', 'kappa', 'V', 'Cn', 'mu'])
    
    components.compile()
    tmo.settings.set_thermo(components)
    
# This just means that if pytest runs this module, it calls the test_component function
if __name__ == '__main__':
    test_component()