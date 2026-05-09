#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = (
    'test_formatting_helpers_cover_common_numeric_ranges',
    'test_load_data_supports_text_csv_npy_and_pickle',
    'test_copy_attr_respects_skip_and_same_options',
    'test_register_with_prefix_handles_explicit_and_ticket_ids',
    'test_sum_system_utility_handles_power_heat_and_invalid_utility',
    'test_create_example_components',
    'test_create_example_system',
    'test_create_example_model',
    )

from pathlib import Path
from shutil import rmtree
from uuid import uuid4

import pytest


def test_formatting_helpers_cover_common_numeric_ranges():
    from qsdsan.utils import format_number, format_str

    assert format_number(None) == 'None'
    assert format_number(0) == 'None'
    assert format_number(float('nan')) == '0'
    assert format_number(2_000_000) == '2.00E+06'
    assert format_number(0.00001) == '1.00E-05'
    assert format_number(0.00543) == '0.0054'
    assert format_number(12.0) == '12'
    assert format_number(12.345) == '12.35'
    assert format_str('mixed string-name') == 'mixed_string_name'


def test_load_data_supports_text_csv_npy_and_pickle():
    import numpy as np
    import pandas as pd
    from qsdsan.utils import load_data, load_pickle, save_pickle

    frame = pd.DataFrame({'value': [1.5, 2.5]}, index=[1., 2.])
    scratch = Path('tmps') / f'test_utils_{uuid4().hex}'
    scratch.mkdir(parents=True)
    try:
        tsv = scratch / 'data.tsv'
        csv = scratch / 'data.csv'
        npy = scratch / 'data.npy'
        pickle = scratch / 'data.pkl'

        frame.to_csv(tsv, sep='\t')
        frame.to_csv(csv)
        np.save(npy, np.array([[1., 1.5], [2., 2.5]]))

        assert load_data(str(tsv)).equals(frame)
        assert load_data(str(csv)).equals(frame)
        assert load_data(str(npy), columns=['value']).equals(frame)

        save_pickle({'loaded': True}, pickle)
        assert load_pickle(pickle) == {'loaded': True}

        with pytest.raises(ValueError, match='Only tab deliminted'):
            load_data(str(scratch / 'data.json'))
    finally:
        rmtree(scratch, ignore_errors=True)


def test_copy_attr_respects_skip_and_same_options():
    from qsdsan.utils import copy_attr

    class Slots:
        __slots__ = ('kept', 'skipped', 'shared')

    original = Slots()
    original.kept = [1]
    original.skipped = ['original']
    original.shared = {'same': True}

    new = Slots()
    new.kept = []
    new.skipped = ['new']
    new.shared = {}

    copy_attr(new, original, skip=('skipped',), same=('shared',))

    assert new.kept == [1]
    assert new.kept is not original.kept
    assert new.skipped == ['new']
    assert new.shared is original.shared


def test_register_with_prefix_handles_explicit_and_ticket_ids():
    from qsdsan.utils import register_with_prefix

    class Registry:
        def __init__(self):
            self.data = {}

        def register(self, ID, obj):
            self.data[ID] = obj

        def register_safely(self, ID, obj):
            self.register(ID, obj)

    class Registered:
        registry = Registry()
        tickets = iter(('1', '2'))

        def _take_ticket(self):
            return next(self.tickets)

    first = Registered()
    register_with_prefix(first, 'unit', '')
    assert first.registry.data['unit_1'] is first

    second = Registered()
    register_with_prefix(second, 'unit', 'pump')
    assert second.registry.data['unit_pump'] is second


def test_sum_system_utility_handles_power_heat_and_invalid_utility():
    from qsdsan.utils import sum_system_utility

    class PowerUtility:
        def __init__(self, consumption, rate):
            self.consumption = consumption
            self.rate = rate

    class HeatUtility:
        def __init__(self, ID, duty, flow):
            self.ID = ID
            self.duty = duty
            self.flow = flow

    class Unit:
        def __init__(self, power_utility=None, heat_utilities=()):
            self.power_utility = power_utility
            self.heat_utilities = list(heat_utilities)

    class System:
        operating_hours = 2

        def __init__(self, units):
            self.units = units

    steam = HeatUtility('low_pressure_steam', 3e6, 1)
    cooling = HeatUtility('cooling_water', -1e6, 1)
    unit = Unit(PowerUtility(consumption=4, rate=-1), (steam, cooling))
    excluded = Unit(PowerUtility(consumption=10, rate=10), ())
    system = System([unit, excluded])

    assert sum_system_utility(system, exclude_units=excluded) == 8
    assert sum_system_utility(
        system,
        exclude_units=excluded,
        utility='electricity',
        calculate_net_utility=True,
        ) == -2
    assert sum_system_utility(system, utility='steam') == 6
    assert sum_system_utility(system, utility='heating') == 6
    assert sum_system_utility(system, utility='cooling') == -2

    with pytest.raises(ValueError, match='utility'):
        sum_system_utility(system, utility='water')


def test_create_example_components():
    from qsdsan.utils import create_example_components
    cmps = create_example_components()
    assert 'H2O' in cmps.IDs
    assert 'CH4' in cmps.IDs
    assert 'Ethanol' in cmps.IDs
    assert len(cmps) == 8
    assert cmps.H2O.particle_size == 'Soluble'

    cmps2 = create_example_components(set_thermo=False)
    assert len(cmps2) == 8


def test_create_example_system():
    from qsdsan.utils import create_example_system
    sys = create_example_system()
    unit_ids = [u.ID for u in sys.path]
    assert len(unit_ids) == 6
    assert 'M1' in unit_ids
    assert 'P1' in unit_ids
    assert 'H1' in unit_ids


def test_create_example_model():
    from qsdsan.utils import create_example_model
    model = create_example_model(evaluate=False)
    assert len(model.parameters) == 6
    assert len(model.metrics) == 4

    model_eval = create_example_model(evaluate=True, N=20, seed=554)
    assert model_eval.table is not None
    assert len(model_eval.table) == 20
