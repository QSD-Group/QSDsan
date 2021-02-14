#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Cheung <joycheung1994@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

import pytest

def test_component():
    import thermosteam as tmo
    from qsdsan import Component, Components
    from chemicals.elements import molecular_weight
    from math import isclose
    
    SNH4 = Component('SNH4', formula='NH4+', measured_as='N', 
                 f_BOD5_COD=0, f_uBOD_COD=0, f_Vmass_Totmass=0,
                 description="Ammonium", particle_size="Soluble",
                 degradability="Undegradable", organic=False)
    assert SNH4.i_N == 1
    assert SNH4.i_NOD == molecular_weight({'O':4})/molecular_weight({'N':1})
    SNH4.measured_as = None
    assert SNH4.i_mass == 1

    SAc = Component('SAc', formula='CH3COO-', measured_as='COD', f_BOD5_COD=0.717, 
                    f_uBOD_COD=0.863, f_Vmass_Totmass=1,
                    description="Acetate", particle_size="Soluble",
                    degradability="Readily", organic=True) 
    assert SAc.i_COD == 1
    SAc.measured_as = None
    assert SAc.i_mass == 1
    assert SAc.i_COD == molecular_weight({'O':4})/molecular_weight({'C':2, 'H':3, 'O':2})
    
    SHS = Component.from_chemical('SHS', tmo.Chemical('Hydrosulfide'), 
                                  particle_size="Soluble",
                                  degradability="Undegradable", organic=False)
    assert SHS.i_charge < 0
    SHS.measured_as = 'S'
    assert SHS.i_mass > 1    
    
    components = Components.load_default(default_compile=False)
    with pytest.raises(AssertionError):
        H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'))
    H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'),
                                  particle_size='Soluble',
                                  degradability='Undegradable', organic=False)
    with pytest.raises(ValueError):
        components.append(H2O)
    components = Components.load_default()
    assert components.SH2.measured_as == 'COD'
    assert components.SH2.i_COD == 1
    assert isclose(components.SN2.i_COD, - molecular_weight({'O':1.5})/molecular_weight({'N':1}), rel_tol=1e-3)
    assert isclose(components.SNO3.i_COD, - molecular_weight({'O':4})/molecular_weight({'N':1}), rel_tol=1e-3)
    tmo.settings.set_thermo(components)
    
# This just means that if pytest runs this module, it calls the test_component function
if __name__ == '__main__':
    test_component()
