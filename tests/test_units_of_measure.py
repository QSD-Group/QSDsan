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

    # 'cu_in' must resolve to a usable cubic-inch unit (pint's default), and
    # 'CFM'/'CFS' aliases must stay registered for QSDsan design units.
    assert abs(uom.auom('cu_in').convert(1, 'm3') - 1.6387064e-5) < 1e-10
    assert uom.auom('CFM').convert(1, 'm3/s') == uom.auom('cfm').convert(1, 'm3/s')
    assert uom.auom('CFS').convert(1, 'm3/s') == uom.auom('cfs').convert(1, 'm3/s')

    # 'each'/'ea' are QSDsan-added dimensionless count units used as the
    # quantity_unit for many Construction items in QSDsan and EXPOsan.
    assert uom.auom('each').convert(1, 'ea') == 1
    assert str(uom.auom('each').dimensionality) == 'dimensionless'

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


def test_units_definition_file_no_redefining_warnings(caplog):
    '''Loading QSDsan's units_of_measure.txt on top of pint + thermosteam
    must not emit any 'Redefining' warnings. Guards against contributors
    adding a definition whose canonical name overlaps a pint or thermosteam
    built-in (which is the bug class this module was restructured to
    eliminate, e.g. the historic 'yr = year = y' triggering pint's
    'Redefining: yr' warning).
    '''
    import logging
    import os

    import pint
    import thermosteam  # noqa: F401  (imported so its definitions are findable)
    import qsdsan.units_of_measure as uom

    qsdsan_dir = os.path.dirname(os.path.realpath(uom.__file__))
    qsdsan_txt = os.path.join(qsdsan_dir, 'units_of_measure.txt')
    thermosteam_txt = os.path.join(
        os.path.dirname(os.path.realpath(thermosteam.__file__)),
        'units_of_measure.txt',
    )

    # Build a clean registry equivalent to what QSDsan loads into at import:
    # pint defaults + thermosteam additions. Then capture warnings only
    # during the QSDsan load, so the assertion is specific to this file.
    fresh = pint.UnitRegistry()
    fresh._on_redefinition = 'warn'
    fresh.load_definitions(thermosteam_txt)

    with caplog.at_level(logging.WARNING, logger='pint.util'):
        fresh.load_definitions(qsdsan_txt)

    redefs = [r.getMessage() for r in caplog.records
              if 'Redefining' in r.getMessage()]
    assert not redefs, f'unexpected pint Redefining warnings: {redefs}'
