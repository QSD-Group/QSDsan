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

import os, pint

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


def _names_in_definition_file(path):
    '''Yield every name added to the registry by a pint definition file.

    Parses ``name = value [= alias ...]`` and ``@alias canonical = alias ...``
    lines, ignoring comments, blank lines, and other ``@`` directives. Used
    to derive the post-load validation list from the file itself so the two
    stay in sync automatically as definitions are added or removed.
    '''
    with open(path) as f:
        for raw_line in f:
            line = raw_line.split('#', 1)[0].strip()
            if not line:
                continue
            if line.startswith('@alias'):
                # '@alias canonical = alias1 [= alias2 ...]' — the canonical
                # name is already registered; only the aliases are new.
                tokens = [p.strip().split()[0]
                          for p in line[len('@alias'):].split('=')]
                yield from tokens[1:]
            elif line.startswith('@'):
                # Other pint directives (@import, @system, @context) don't
                # add named units.
                continue
            else:
                # 'name = value [= alias1 [= alias2 ...]]'
                tokens = [p.strip().split()[0] for p in line.split('=')]
                yield tokens[0]
                yield from tokens[2:]   # skip the value expression


# Load QSDsan's additions to the shared pint registry from an external
# definition file, mirroring the pattern in thermosteam/units_of_measure.py.
# The reload guard prevents re-parsing the file (and a stream of pint
# 'Redefining' warnings) if qsdsan is reloaded in the same Python session.
if not getattr(pint, 'QSDsan_units_loaded', False):
    _txt_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             'units_of_measure.txt')
    ureg.load_definitions(_txt_path)
    pint.QSDsan_units_loaded = True

    # Fail loudly at import time if any name QSDsan just registered does
    # not resolve to a usable unit. pint defers RHS resolution until first
    # use, so a typo or keyword collision (e.g. the historic 'cu_in = in3'
    # bug, where 'in' is a Python keyword) would otherwise silently produce
    # a broken unit. The list is derived from the .txt file itself so it
    # stays in sync as definitions are added or removed.
    for _name in _names_in_definition_file(_txt_path):
        auom(_name)
    del _name, _txt_path

del os, pint


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
