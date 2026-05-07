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

__all__ = ('test_units_of_measure',)


def test_units_of_measure():
    import qsdsan as qs
    import qsdsan.units_of_measure as uom
    from qsdsan.utils import parse_unit

    assert qs.units_of_measure is uom
    assert qs.utils.auom is uom.auom

    assert uom.auom('MGD').convert(1, 'm3/d') > 3000
    assert uom.auom('mmscfd').convert(1, 'm3/d') > 20000
    assert uom.auom('point').convert(1, 'points') == 1

    unit, remaining = uom.parse_lca_unit('kg CO2-eq')
    assert str(unit) == 'kg'
    assert remaining == 'CO2-eq'

    unit, remaining = uom.parse_lca_unit('MJ-eq')
    assert str(unit) == 'MJ'
    assert remaining == 'eq'

    unit, remaining = uom.parse_lca_unit('kg*km')
    assert str(unit) == 'kg*km'
    assert remaining == ''

    assert parse_unit('kg CO2-eq') == uom.parse_lca_unit('kg CO2-eq')
