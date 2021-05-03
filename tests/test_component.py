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

__all__ = ('test_component',)

def test_component():
    import pytest
    import thermosteam as tmo
    from qsdsan import Component, Components
    from chemicals.elements import molecular_weight
    from math import isclose
    
    S_NH4 = Component('S_NH4', formula='NH4+', measured_as='N', 
                 f_BOD5_COD=0, f_uBOD_COD=0, f_Vmass_Totmass=0,
                 description="Ammonium", particle_size="Soluble",
                 degradability="Undegradable", organic=False)
    assert S_NH4.i_N == 1
    assert S_NH4.i_NOD == molecular_weight({'O':4})/molecular_weight({'N':1})
    S_NH4.measured_as = None
    assert S_NH4.i_mass == 1

    S_Ac = Component('S_Ac', formula='CH3COO-', measured_as='COD', f_BOD5_COD=0.717, 
                    f_uBOD_COD=0.863, f_Vmass_Totmass=1,
                    description="Acetate", particle_size="Soluble",
                    degradability="Readily", organic=True) 
    assert S_Ac.i_COD == 1
    S_Ac.measured_as = None
    assert S_Ac.i_mass == 1
    assert S_Ac.i_COD == molecular_weight({'O':4})/molecular_weight({'C':2, 'H':3, 'O':2})
    
    S_HS = Component.from_chemical('S_HS', tmo.Chemical('Hydrosulfide'), 
                                  particle_size="Soluble",
                                  degradability="Undegradable", organic=False)
    assert S_HS.i_charge < 0
    S_HS.measured_as = 'S'
    assert S_HS.i_mass > 1    
    
    components = Components.load_default(default_compile=False)

    #!!! Should we allow None for particle_size, degradability, and organic?
    # with pytest.raises(AssertionError):
    #     H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'))
    
    H2O = Component.from_chemical('H2O', tmo.Chemical('H2O'),
                                  particle_size='Soluble',
                                  degradability='Undegradable', organic=False)
    with pytest.raises(ValueError):
        components.append(H2O)
    components = Components.load_default()
    assert components.S_H2.measured_as == 'COD'
    assert components.S_H2.i_COD == 1
    assert isclose(components.S_N2.i_COD, - molecular_weight({'O':1.5})/molecular_weight({'N':1}), rel_tol=1e-3)
    assert isclose(components.S_NO3.i_COD, - molecular_weight({'O':4})/molecular_weight({'N':1}), rel_tol=1e-3)
    tmo.settings.set_thermo(components)
