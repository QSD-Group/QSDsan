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

from thermosteam.units_of_measure import (
    ureg,
    AbsoluteUnitsOfMeasure as auom,
    RelativeUnitsOfMeasure as ruom,
    UnitsOfMeasure,
    DisplayUnits,
    DisplayNotation,
    Quantity,
    convert,
    format_units,
    format_plot_units,
    parse_units_notation,
    reformat_units,
    get_dimensionality,
    chemical_units_of_measure,
    stream_units_of_measure,
    power_utility_units_of_measure,
    heat_utility_units_of_measure,
    )

__all__ = (
    'ureg',
    'auom',
    'ruom',
    'UnitsOfMeasure',
    'DisplayUnits',
    'DisplayNotation',
    'Quantity',
    'convert',
    'format_units',
    'format_plot_units',
    'parse_units_notation',
    'reformat_units',
    'get_dimensionality',
    'chemical_units_of_measure',
    'stream_units_of_measure',
    'power_utility_units_of_measure',
    'heat_utility_units_of_measure',
    'component_units_of_measure',
    'waste_stream_units_of_measure',
    'lca_units_of_measure',
    'construction_units_of_measure',
    'transportation_units_of_measure',
    'parse_lca_unit',
    )


def _unit_is_defined(name):
    return name in ureg._units


def _define_unit(definition):
    names = tuple(i.strip().split()[0] for i in definition.split('='))
    if all(_unit_is_defined(i) for i in names):
        return
    ureg.define(definition)


for _definition in (
        'sq_m = m2',
        'cu_m = m3',
        'sq_cm = cm2',
        'cu_cm = cm3',
        'sq_ft = ft2',
        'cu_ft = ft3',
        'cu_in = in3',
        'yd3 = yard**3 = cu_yd',
        'cfm = cf/minute = CFM',
        'cfs = cf/second = CFS',
        'yr = year = y',
        'hr = hour = h',
        'd = day',
        'each = count = ea',
        'unit = count',
        'point = count = points',
        'MGD = 1e6 * gallon / day',
        'mgd = MGD',
        'mmscfd = 1e6 * ft3 / day = MMSCfd = MMSCFD',
        ):
    _define_unit(_definition)
del _definition


component_units_of_measure = {
    'i_C': auom('g'),
    'i_N': auom('g'),
    'i_P': auom('g'),
    'i_K': auom('g'),
    'i_Mg': auom('g'),
    'i_Ca': auom('g'),
    'i_mass': auom('g'),
    'i_charge': auom('mol'),
    'i_COD': auom('g'),
    'i_NOD': auom('g'),
}

waste_stream_units_of_measure = {
    'flow': auom('kg/hr'),
    'molar_flow': auom('kmol/hr'),
    'concentration': auom('mg/L'),
    'charge_concentration': auom('mmol/L'),
    'volumetric_flow': auom('L/hr'),
}

lca_units_of_measure = {
    'mass': auom('kg'),
    'energy': auom('MJ'),
    'transportation': auom('kg*km'),
    'volume_transportation': auom('m3*km'),
    'points': auom('point'),
}

construction_units_of_measure = {
    'quantity': None,
    'lifetime': auom('yr'),
}

transportation_units_of_measure = {
    'distance': auom('km'),
    'interval': auom('hr'),
    'mass_load': auom('kg'),
    'volume_load': auom('m3'),
    'mass_quantity': auom('kg*km'),
    'volume_quantity': auom('m3*km'),
}


def parse_lca_unit(value):
    '''Parse impact assessment units into a pint-compatible unit and suffix.'''
    if not value:
        return None, value

    value = str(value).strip()

    parts = value.split()
    if len(parts) > 1:
        unit = parts[0]
        remaining = ' '.join(parts[1:])
        try:
            return auom(unit), remaining
        except Exception:
            pass

    unit, sep, remaining = value.partition('-')
    if sep:
        try:
            return auom(unit), remaining
        except Exception:
            pass

    try:
        return auom(value), ''
    except Exception:
        return None, value
