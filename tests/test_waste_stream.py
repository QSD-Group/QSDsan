#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 10:55:35 2020

@author: yalinli_cabbi
"""

import pytest
from numpy.testing import assert_allclose

def test_waste_stream():
    import thermosteam as tmo
    from sanitation import utils, WasteStream
    components = utils.load_default_components()
    tmo.settings.set_thermo(components)
    ws1 = WasteStream(SAc=5, H2O=1000, units='kg/hr')
    ws2 = WasteStream(SHAc=10, H2O=1000, units='kg/hr')
    ws3 = WasteStream()
    ws3.mix_from((ws1, ws2))
    assert_allclose(ws3.F_mass, 2015.0)
    # TODO: After updating the default component properties,
    # add in tests here to make sure COD, etc. are calculated correctly
    # assert_allclose(ws3.F_mass, 2015.0)
    
    # Make sure below attributes are calculated based on flow info, cannot be set
    with pytest.raises(AttributeError):
        ws3.COD = 5
    
# This just means that if pytest runs this module, it calls the test_waste_stream function
if __name__ == '__main__':
    test_waste_stream()