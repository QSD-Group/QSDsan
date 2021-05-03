#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

__all__ = ('test_sanunit',)

def test_sanunit():
    from numpy.testing import assert_allclose
    import qsdsan as qs
    components = qs.Components.load_default()
    qs.set_thermo(components)

    # Test unit initiation
    ws1 = qs.WasteStream(S_Ac=5, H2O=1000, units='kg/hr')
    ws2 = qs.WasteStream(X_NOO=10, H2O=1000, units='kg/hr')
    M1 = qs.sanunits.Mixer('M1', ins=(ws1, ws2, ''), outs='mixture',
                           init_with='WasteStream')
    M1.show()
    assert type(M1.ins[0]).__name__ == 'WasteStream'

    S1 = qs.sanunits.Splitter('S1', ins=M1-0, outs=('', ''), split=0.2)
    ins3 = qs.WasteStream(S_CH3OH=7, H2O=1000, units='kg/hr')
    P1 = qs.sanunits.Pump('P1', ins=ins3)
    M2 = qs.sanunits.MixTank('M2', ins=(S1-0, P1-0), tau=2)
    sys = qs.System('sys', path=(M1, S1, P1, M2))
    sys.simulate()

    assert_allclose(M2.installed_cost, 41808.1524967656, rtol=1e-3)

    # Test mixing of different classes of streams
    ss1 = qs.SanStream(H2O=100)
    ss2 = ss1.copy()

    M3 = qs.sanunits.MixTank('M3', ins=ss1, init_with='WasteStream')
    M3.F_BM['Tanks'] = 1
    M3.simulate()
    M3.show()
    assert type(M3.ins[0]).__name__ == 'SanStream'
    assert type(M3.outs[0]).__name__ == 'WasteStream'
    assert_allclose(M3.installed_cost, 4386.336513753271, rtol=1e-3)

    M4 = qs.sanunits.MixTank('M4', ins=ss2, init_with='Stream')
    assert type(M4.outs[0]).__name__ == 'Stream'
    M4.simulate()
    M4.show()
    assert_allclose(M4.installed_cost, 7237.455247692897, rtol=1e-3)