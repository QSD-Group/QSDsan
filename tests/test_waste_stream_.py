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

__all__ = ('test_waste_stream',)

def test_waste_stream():
    import pytest, numpy as np
    from numpy.testing import assert_allclose
    from math import isclose
    from qsdsan import set_thermo, Components, WasteStream

    components = Components.load_default()
    set_thermo(components)

    ws1 = WasteStream.codstates_inf_model('ws1', 1e5)
    ws2 = WasteStream.codstates_inf_model('ws2', 1e5*24/1e3, units=('m3/d', 'g/m3'))
    assert isclose(ws1.COD, 430, rel_tol=1e-2)
    assert isclose(ws1.TKN, 40, rel_tol=1e-2)
    assert isclose(ws1.TP, 10, rel_tol=1e-2)
    assert isclose(ws1.F_vol, ws2.F_vol)

    ws3 = WasteStream(S_Ac=5, H2O=1000, units='kg/hr')
    ws4 = WasteStream(X_NOO=10, H2O=1000, units='kg/hr')
    ws5 = WasteStream()
    ws5.mix_from((ws3, ws4))
    assert_allclose(ws5.F_mass, 2015.0)
    # TODO: After updating the default component properties,
    # add in tests here to make sure COD, etc. are calculated correctly
    assert_allclose(ws5.COD, 7414.267796, rtol=1e-2)

    # Make sure below attributes are calculated based on flow info, cannot be set
    with pytest.raises(AttributeError):
        ws5.COD = 5

    # Concentration calclation
    ws6 = WasteStream(X_CaCO3=1, H2O=1000, units='kg/hr')
    assert_allclose(np.abs(ws6.conc.value-ws6.mass/ws6.F_vol*1e3).sum(), 0, atol=1e-6)
    ws6.imass['X_B_Subst', 'X_GAO_PHA'] = (100, 1)
    ws7 = WasteStream(X_CaCO3=1, X_B_Subst=100, X_GAO_PHA=1, H2O=1000, units='kg/hr')
    assert_allclose(np.abs(ws6.conc.value-ws7.mass/ws7.F_vol*1e3).sum(), 0, atol=1e-6)
    ws6.mass[:] = 1e-3
    ws6.imass['H2O'] = 1e3
    diff = ws6.conc.value - np.ones_like(ws6.conc.value)
    diff[components.index('H2O')] = 0
    assert_allclose(np.max(np.abs(diff)), 0, atol=1e-2)


if __name__ == '__main__':
    test_waste_stream()