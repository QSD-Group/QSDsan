#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, biosteam as bst, qsdsan as qs
from numpy.testing import assert_allclose

__all__ = (
    'test_default',
    'test_BinaryDistillation',
    'test_Flash',
    'test_HXprocess',
    'test_HXutility',
    'test_IsothermalCompressor',
    'test_MixTank',
    'test_Pump',
    'test_Splitter',
    'test_StorageTank',
    )


bst.default()
chems = bst.Chemicals(('Water', 'Methanol', 'Ethanol', 'Glycerol', 'H2'))

ws_data = {
    'particle_size': 'Soluble',
    'degradability': 'Readily',
    'organic': True
    }
cmps = qs.Components([
    qs.Component('Water',
                 particle_size='Soluble',
                 degradability='Undegradable',
                 organic=False),
    qs.Component('Methanol', **ws_data),
    qs.Component('Ethanol', **ws_data),
    qs.Component('Glycerol', **ws_data),
    qs.Component('H2',
                 particle_size='Dissolved gas',
                 degradability='Readily',
                 organic=True),
    ])


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


def check_results(biosteam_unit, qsdsan_unit):
    global bst_unit, qs_unit
    bst_unit, qs_unit = biosteam_unit, qsdsan_unit
    bst_unit.simulate()
    qs_unit.simulate()

    bst_s = bst_unit.ins + bst_unit.outs
    qs_ws = qs_unit.ins + qs_unit.outs
    for n, s in enumerate(bst_s):
        assert_allclose(np.abs(s.mol-qs_ws[n].mol).sum(), 0, atol=1e-6)

    assert_allclose(bst_unit.installed_cost, qs_unit.installed_cost, atol=1e-6)
    bst_unit_ucost = 0 if bst_unit.utility_cost is None else bst_unit.utility_cost
    qs_unit_ucost = 0 if qs_unit.utility_cost is None else qs_unit.utility_cost   
    assert_allclose(bst_unit_ucost, qs_unit_ucost, atol=1e-6)
    assert_allclose(bst_unit.power_utility.rate, qs_unit.power_utility.rate, atol=1e-6)



# %%

# =============================================================================
# Testing functions
# =============================================================================

def test_default():
    qs.default() # default everything


def test_BinaryDistillation():
    bst.settings.set_thermo(chems)
    stream_kwargs = dict(Water=80, Methanol=100, Glycerol=25, units='kmol/hr')
    bst_s = bst.Stream(**stream_kwargs)
    bst_s.T = T = bst_s.bubble_point_at_P().T # Feed at bubble point T
    unit_kwargs = dict(LHK=('Methanol', 'Water'), y_top=0.99, x_bot=0.01, k=2, is_divided=True)
    bst_unit = bst.units.BinaryDistillation(ins=bst_s, **unit_kwargs)
    bst_unit.simulate() # need extra simulation for better convergence
    
    qs.set_thermo(cmps)
    qs_s = qs.WasteStream(T=T, **stream_kwargs)
    qs_unit = qs.sanunits.BinaryDistillation(ins=qs_s, **unit_kwargs)
    qs_unit.simulate()
    
    check_results(bst_unit, qs_unit)


def test_Flash():
    bst.settings.set_thermo(chems)
    stream_kwargs = dict(Glycerol=300, Water=1000, units='kmol/hr')
    bst_s = bst.Stream(**stream_kwargs)
    bst_s.T = T = bst_s.bubble_point_at_P().T # Feed at bubble point T
    unit_kwargs = dict(P=101325, T=410.15)
    bst_unit = bst.units.Flash(ins=bst_s, **unit_kwargs)
    
    qs.set_thermo(cmps)
    qs_s = qs.WasteStream(T=T, **stream_kwargs)
    qs_unit = qs.sanunits.Flash(ins=qs_s, **unit_kwargs)
    
    check_results(bst_unit, qs_unit)


def test_HXprocess():
    bst_s, qs_ws = create_streams(2)
    bst_s[0].T = qs_ws[0].T = 400

    bst.settings.set_thermo(chems)
    bst_unit = bst.units.HXprocess(ins=bst_s, phase0='l', phase1='l')

    qs.set_thermo(cmps)
    qs_unit = qs.sanunits.HXprocess(ins=qs_ws, phase0='l', phase1='l')

    check_results(bst_unit, qs_unit)


def test_HXutility():
    bst_s, qs_ws = create_streams(1)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.HXutility(ins=bst_s, T=400, rigorous=True)

    qs.set_thermo(cmps)
    qs_unit = qs.sanunits.HXutility(ins=qs_ws, T=400, rigorous=True)

    check_results(bst_unit, qs_unit)


def test_IsothermalCompressor():
    bst.settings.set_thermo(chems)
    stream_kwargs = dict(H2=1, T=298.15, P=20e5, units='kmol/hr', phase='g')
    bst_s = bst.Stream(**stream_kwargs)
    unit_kwargs = dict(P=350e5, eta=1)
    bst_unit = bst.units.IsothermalCompressor(ins=bst_s, **unit_kwargs)
    
    qs.set_thermo(cmps)
    qs_s = qs.WasteStream(**stream_kwargs)
    qs_unit = qs.sanunits.IsothermalCompressor(ins=qs_s, **unit_kwargs)
    
    check_results(bst_unit, qs_unit)


def test_MixTank():
    bst_s, qs_ws = create_streams(2)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.MixTank(ins=bst_s)

    qs.set_thermo(cmps)
    qs_unit = qs.sanunits.MixTank(ins=qs_ws)

    check_results(bst_unit, qs_unit)


def test_Pump():
    bst_s, qs_ws = create_streams(1)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.Pump(ins=bst_s)

    qs.set_thermo(cmps)
    qs_unit = qs.sanunits.Pump(ins=qs_ws)

    check_results(bst_unit, qs_unit)


def test_Splitter():
    bst_s, qs_ws = create_streams(1)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.Splitter(ins=bst_s, split=0.1)

    qs.set_thermo(cmps)
    qs_unit = qs.sanunits.Splitter(ins=qs_ws, split=0.1)

    check_results(bst_unit, qs_unit)


def test_StorageTank():
    bst_s, qs_ws = create_streams(1)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.StorageTank(ins=bst_s)

    qs.set_thermo(cmps)
    qs_unit = qs.sanunits.StorageTank(ins=qs_ws)

    check_results(bst_unit, qs_unit)


if __name__ == '__main__':
    test_default()
    test_BinaryDistillation()
    test_Flash()
    test_HXprocess()
    test_HXutility()
    test_IsothermalCompressor()
    test_MixTank()
    test_Pump()
    test_Splitter()
    test_StorageTank()