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

from numpy.testing import assert_allclose

def test_sanunit():
    import biosteam as bst
    from qsdsan import Components, WasteStream, sanunits
    components = Components.load_default()
    bst.settings.set_thermo(components)
    ws1 = WasteStream(S_Ac=5, H2O=1000, units='kg/hr')
    ws2 = WasteStream(X_NOO=10, H2O=1000, units='kg/hr')
    M1 = sanunits.Mixer('M1', ins=(ws1, ws2, ''), outs='mixture')
    M1.show()
    assert type(M1.ins[0]).__name__ == 'WasteStream'
    
    S1 = sanunits.Splitter('S1', ins=M1-0, outs=('', ''), split=0.2)
    ins3 = WasteStream(S_CH3OH=7, H2O=1000, units='kg/hr')
    P1 = sanunits.Pump('P1', ins=ins3)    
    M2 = sanunits.MixTank('M2', ins=(S1-0, P1-0), tau=2)
    M2-0-2-M1
    System = bst.System('System', path=(M1, S1, P1, M2), recycle=M2-0)
    System.simulate()
    assert_allclose(M2.installed_cost, 65519.00446342958, rtol=1e-3)


if __name__ == '__main__':
    test_sanunit()