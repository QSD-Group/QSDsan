#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import biosteam as bst
import thermosteam as tmo

from qsdsan import _compat


def test_biosteam_compatibility_aliases():
    assert _compat.ProcessSpecification is bst._unit.ProcessSpecification
    assert _compat.UnitGraphics is bst._graphics.UnitGraphics
    expected_timer = getattr(bst.utils, 'Timer', None)
    if expected_timer is None:
        expected_timer = bst.utils.TicToc
    assert _compat.Timer is expected_timer
    assert hasattr(_compat.NotImplementedMethod, '__call__')


def test_thermosteam_chemical_compatibility_aliases():
    assert _compat.chemical_fields is tmo._chemical._chemical_fields
    assert _compat.checked_properties is tmo._chemical._checked_properties
    assert _compat.lock_phase is tmo._chemical.lock_phase
    assert _compat.display_asfunctor is tmo._chemical.display_asfunctor
    assert _compat.get_chemical_data is tmo._chemical.get_chemical_data
    assert _compat.prepare_chemicals is tmo._chemicals.prepare
    assert _compat.chemical_data_array is tmo._chemicals.chemical_data_array
