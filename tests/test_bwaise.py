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

__all__ = ('test_bwaise',)

def test_bwaise():
    from numpy.testing import assert_allclose
    from exposan import bwaise as bw

    assert_allclose(bw.teaA.NPV, -24526872.85395059, rtol=1e-3)
    assert_allclose(bw.teaB.NPV, -2231960.1805307646, rtol=1e-3)
    assert_allclose(bw.teaC.NPV, -97235311.17219843, rtol=1e-3)

    assert_allclose(bw.lcaA.total_impacts['GlobalWarming'], 128530982.2469723, rtol=1e-3)
    assert_allclose(bw.lcaB.total_impacts['GlobalWarming'], 11987512.603749081, rtol=1e-3)
    assert_allclose(bw.lcaC.total_impacts['GlobalWarming'], 69592120.06238395, rtol=1e-3)