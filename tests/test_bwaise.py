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

__all__ = ('test_bwaise',)

def test_bwaise():
    from exposan import bwaise as bw
    
    assert_allclose(bw.teaA.NPV, -22732728.213841617, rtol=1e-3)
    assert_allclose(bw.teaB.NPV, -2231960.180530765, rtol=1e-3)
    assert_allclose(bw.teaC.NPV, -94340457.26937833, rtol=1e-3)

    assert_allclose(bw.lcaA.total_impacts['GlobalWarming'], 146386354.786746, rtol=1e-3)
    assert_allclose(bw.lcaB.total_impacts['GlobalWarming'], 8794967.822499081, rtol=1e-3)
    assert_allclose(bw.lcaC.total_impacts['GlobalWarming'], 56832681.43125322, rtol=1e-3)
    