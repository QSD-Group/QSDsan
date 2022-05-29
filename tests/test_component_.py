#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

__all__ = ('test_component',)

def test_component():
    import pytest
    from qsdsan import Chemical, Component, Components, set_thermo, \
        _waste_stream as ws_module
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

    S_HS = Component.from_chemical('S_HS', Chemical('Hydrosulfide'),
                                  particle_size="Soluble",
                                  degradability="Undegradable", organic=False)
    assert S_HS.i_charge < 0
    S_HS.measured_as = 'S'
    assert S_HS.i_mass > 1

    # Check default components
    cmps1 = Components.load_default(default_compile=False)
    H2O_chemical = Chemical('H2O')
    H2O = Component.from_chemical('H2O', H2O_chemical)
    with pytest.raises(ValueError): # H2O already in default components
        cmps1.append(H2O)
    with pytest.raises(RuntimeError): # key chemical-related properties missing
        cmps1.compile()
    # Can compile with default-filling those missing properties
    cmps1.default_compile(lock_state_at='', particulate_ref='NaCl')

    cmps2 = Components((cmp for cmp in cmps1 if cmp.ID != 'H2O'))
    H2O = Component.from_chemical('H2O', Chemical('H2O'),
                                  particle_size='Soluble',
                                  degradability='Undegradable', organic=False)
    cmps2.append(H2O)
    cmps2.default_compile(lock_state_at='', particulate_ref='NaCl')

    cmps3 = Components.load_default()
    assert cmps3.S_H2.measured_as == 'COD'
    assert cmps3.S_H2.i_COD == 1
    assert isclose(cmps3.S_N2.i_COD, - molecular_weight({'O':1.5})/molecular_weight({'N':1}), rel_tol=1e-3)
    assert isclose(cmps3.S_NO3.i_COD, - molecular_weight({'O':4})/molecular_weight({'N':1}), rel_tol=1e-3)
    set_thermo(cmps3)

    # Check if the default groups are up-to-date
    cached_cmp_IDs = ws_module._default_cmp_IDs
    cached_cmp_groups = ws_module._specific_groups
    assert set(cmps3.IDs) == cached_cmp_IDs
    get_IDs = lambda attr: {cmp.ID for cmp in getattr(cmps3, attr)}
    for attr, IDs in cached_cmp_groups.items():
        assert IDs == get_IDs(attr)


if __name__ == '__main__':
    test_component()