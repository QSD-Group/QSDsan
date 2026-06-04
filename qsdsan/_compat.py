#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Compatibility aliases for BioSTEAM/Thermosteam APIs used by QSDsan.

This module centralizes upstream attributes that are version-sensitive or not
part of the most stable public import surface.
'''

import thermosteam as tmo

from biosteam._graphics import UnitGraphics
from biosteam._unit import ProcessSpecification

try:
    from biosteam.utils import NotImplementedMethod
except ImportError:
    from biosteam.utils import AbstractMethod as NotImplementedMethod

try:
    from biosteam.utils import Timer
except ImportError:
    from biosteam.utils import TicToc as Timer

chemical_fields = tmo._chemical._chemical_fields
checked_properties = tmo._chemical._checked_properties
lock_phase = tmo._chemical.lock_phase
display_asfunctor = tmo._chemical.display_asfunctor
get_chemical_data = tmo._chemical.get_chemical_data
prepare_chemicals = tmo._chemicals.prepare
chemical_data_array = tmo._chemicals.chemical_data_array
_reset_constant = getattr(tmo._chemical, 'reset_constant', None)


def set_chemical_MW(chemical, MW):
    '''
    Set a chemical's molecular weight in a way that works across Thermosteam
    versions.

    Older releases expose a writable ``Chemical.MW`` property; newer releases
    make it read-only (the setter raises ``AttributeError``). Try the public
    setter first (so versions where it works keep using the public API), and
    fall back to Thermosteam's internal constant-reset helper, which is what the
    writable setter itself uses under the hood.
    '''
    MW = float(MW)
    try:
        chemical.MW = MW
    except AttributeError:
        if _reset_constant is not None:
            _reset_constant(chemical, 'MW', MW)
        else:  # pragma: no cover - extremely unlikely with supported versions
            chemical._MW = MW


__all__ = (
    'UnitGraphics',
    'ProcessSpecification',
    'NotImplementedMethod',
    'Timer',
    'chemical_fields',
    'checked_properties',
    'lock_phase',
    'display_asfunctor',
    'get_chemical_data',
    'prepare_chemicals',
    'chemical_data_array',
    'set_chemical_MW',
    )
