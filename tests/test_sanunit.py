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

from numpy.testing import assert_allclose

def test_sanunit():
    import biosteam as bst
    from sanitation import Components, WasteStream, units
    components = Components.load_default()
    bst.settings.set_thermo(components)
    ws1 = WasteStream(SAc=5, H2O=1000, units='kg/hr')
    ws2 = WasteStream(XNOO=10, H2O=1000, units='kg/hr')
    M1 = units.Mixer('M1', ins=(ws1, ws2, ''), outs='mixture')
    assert type(M1.ins[0]).__name__ == 'WasteStream'
    
    S1 = units.Splitter('S1', ins=M1-0, outs=('', ''), split=0.2)
    ins3 = WasteStream(SCH3OH=7, H2O=1000, units='kg/hr')
    P1 = units.Pump('P1', ins=ins3)    
    M2 = units.MixTank('M2', ins=(S1-0, P1-0), tau=2)
    M2-0-2-M1
    System = bst.System('System', path=(M1, S1, P1, M2), recycle=M2-0)
    System.simulate()
    assert_allclose(M2.installed_cost, 65735.21240525639, rtol=1e-3)


# This just means that if pytest runs this module, it calls the test_sanunit function
if __name__ == '__main__':
    test_sanunit()