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


import numpy as np
import biosteam as bst
import qsdsan as qs
from numpy.testing import assert_allclose

__all__ = ('test_Splitter', 'test_Pump', 'test_MixTank', 'test_StorageTank',
           'test_HXutility', 'test_HXprocess',)


bst.default_utilities()
chems = bst.Chemicals(('Methanol', 'Ethanol'))

ws_data = {
    'particle_size': 'Soluble',
    'degradability': 'Readily',
    'organic': True
    }
cmps = qs.Components((qs.Component('Methanol', search_ID='Methanol', **ws_data),
                      qs.Component('Ethanol', search_ID='Ethanol', **ws_data)))


def create_streams(num):
    bst.settings.set_thermo(chems)
    bst_s = []
    for n in range(num):
        s = bst.Stream(Methanol=100*(n+1), Ethanol=100*(n+1), units='kg/hr')
        bst_s.append(s)

    qs.set_thermo(cmps)
    qs_ws = []
    for n in range(num):
        ws = qs.WasteStream(Methanol=100*(n+1), Ethanol=100*(n+1), units='kg/hr')
        qs_ws.append(ws)
    
    return bst_s, qs_ws
    

def check_results(bst_unit, qs_unit):
    bst_unit.simulate()
    qs_unit.simulate()
    
    bst_s = bst_unit.ins + qs_unit.outs
    qs_ws = qs_unit.ins + qs_unit.outs
    for n, s in enumerate(bst_s):
        assert_allclose(np.abs(s.mol-qs_ws[n].mol).sum(), 0, atol=1e-6)
        
    assert_allclose(bst_unit.installed_cost, qs_unit.installed_cost, atol=1e-6)
    assert_allclose(bst_unit.utility_cost, qs_unit.utility_cost, atol=1e-6)
    assert_allclose(bst_unit.power_utility.rate, qs_unit.power_utility.rate, atol=1e-6)

    

# %%

# =============================================================================
# Testing functions
# =============================================================================

def test_Splitter():
    bst_s, qs_ws = create_streams(1)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.Splitter(ins=bst_s, split=0.1)
    
    qs.set_thermo(cmps)
    qs_unit = qs.sanunits.Splitter(ins=qs_ws, split=0.1)
    
    check_results(bst_unit, qs_unit)


def test_Pump():
    bst_s, qs_ws = create_streams(1)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.Pump(ins=bst_s)
    
    qs.set_thermo(cmps)
    qs_unit = qs.sanunits.Pump(ins=qs_ws)
    
    check_results(bst_unit, qs_unit)


def test_MixTank():
    bst_s, qs_ws = create_streams(2)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.MixTank(ins=bst_s)
    
    qs.set_thermo(cmps)
    qs_unit = qs.sanunits.MixTank(ins=qs_ws)
    
    check_results(bst_unit, qs_unit)


def test_StorageTank():
    bst_s, qs_ws = create_streams(1)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.StorageTank(ins=bst_s)
    
    qs.set_thermo(cmps)
    qs_unit = qs.sanunits.StorageTank(ins=qs_ws)
    
    check_results(bst_unit, qs_unit)


def test_HXutility():
    bst_s, qs_ws = create_streams(1)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.HXutility(ins=bst_s, T=400, rigorous=False) #!!! Try True
    
    qs.set_thermo(cmps)
    qs_unit = qs.sanunits.HXutility(ins=qs_ws, T=400, rigorous=False) #!!! Try True
    
    check_results(bst_unit, qs_unit)

def test_HXprocess():
    bst_s, qs_ws = create_streams(2)
    bst_s[0].T = qs_ws[0].T = 400
    
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.HXprocess(ins=bst_s, phase0='l', phase1='l')
    
    qs.set_thermo(cmps)
    qs_unit = qs.sanunits.HXprocess(ins=qs_ws, phase0='l', phase1='l')
    
    check_results(bst_unit, qs_unit)
