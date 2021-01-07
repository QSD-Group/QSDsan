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
from numpy.testing import assert_allclose

def test_waste_stream():
    import thermosteam as tmo
    from qsdsan import Components, WasteStream
    components = Components.load_default()
    tmo.settings.set_thermo(components)
    ws1 = WasteStream(SAc=5, H2O=1000, units='kg/hr')
    ws2 = WasteStream(XNOO=10, H2O=1000, units='kg/hr')
    ws3 = WasteStream()
    ws3.mix_from((ws1, ws2))
    assert_allclose(ws3.F_mass, 2015.0)
    # TODO: After updating the default component properties,
    # add in tests here to make sure COD, etc. are calculated correctly
    assert_allclose(ws3.COD, 7424.606289711915, rtol=1e-3)
    
    # Make sure below attributes are calculated based on flow info, cannot be set
    with pytest.raises(AttributeError):
        ws3.COD = 5
    
# This just means that if pytest runs this module, it calls the test_waste_stream function
if __name__ == '__main__':
    test_waste_stream()