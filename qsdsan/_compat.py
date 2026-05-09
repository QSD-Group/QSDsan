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
    )
