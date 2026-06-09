# -*- coding: utf-8 -*-
"""
Tests for QSDsan's self-contained public API surface.

The goal is that users rarely need to ``import biosteam``/``import thermosteam``
directly: the common surface is reachable through ``qsdsan``. These tests assert
that each re-exported name resolves and *is the same object* as its
BioSTEAM/Thermosteam source. The identity (``is``) checks are deliberate: if
BioSTEAM rebinds or renames one of these globals in a future release, the
re-export silently goes stale, and these tests fail loudly (in CI, before a
release) instead of a user hitting it.
"""

import biosteam as bst
import thermosteam as tmo
import qsdsan as qs


# %% Top-level plain re-exports (objects / classes / functions)

def test_toplevel_reexports_match_biosteam():
    # New re-exports added to make qsdsan self-contained.
    assert qs.settings is bst.settings
    assert qs.preferences is bst.preferences
    assert qs.stream_utility_prices is bst.stream_utility_prices
    assert qs.Thermo is bst.Thermo
    assert qs.UtilityAgent is bst.UtilityAgent
    assert qs.Facility is bst.Facility
    assert qs.AgileSystem is bst.AgileSystem
    assert qs.get_OSBL is bst.get_OSBL
    assert qs.MissingStream is bst.utils.MissingStream


def test_existing_toplevel_reexports_match_biosteam():
    # Pre-existing re-exports; covered here so a BioSTEAM rebind also trips a test.
    assert qs.Stream is bst.Stream
    assert qs.MultiStream is bst.MultiStream
    assert qs.Chemical is bst.Chemical
    assert qs.Chemicals is bst.Chemicals
    assert qs.HeatUtility is bst.HeatUtility
    assert qs.PowerUtility is bst.PowerUtility
    assert qs.Unit is bst.Unit
    assert qs.System is bst.System
    assert qs.Model is bst.Model
    assert qs.Metric is bst.Metric
    assert qs.Parameter is bst.Parameter
    assert qs.Reaction is bst.Reaction


def test_settings_is_mutated_in_place_not_rebound():
    # Plain re-export is correct only because callers mutate the singleton in
    # place (e.g. ``settings.thermo.ideal()``) rather than rebinding the name.
    # If this identity ever breaks, switch to a read-only accessor property.
    assert qs.settings is bst.settings
    assert qs.preferences is bst.preferences
    assert qs.stream_utility_prices is bst.stream_utility_prices


# %% Mutable-global behavior

def test_cepci_is_live_view_of_bst_CE():
    # qs.CEPCI is a settable property that writes through to bst.CE.
    old = bst.CE
    try:
        bst.CE = 123.0
        assert qs.CEPCI == 123.0           # read-through
        qs.CEPCI = 456.0
        assert bst.CE == 456.0             # write-through
    finally:
        bst.CE = old


def test_powerutility_price_writes_through_via_class():
    # No accessor needed: qs.PowerUtility is the same class object as
    # bst.PowerUtility, so setting the class attribute is shared global state.
    assert qs.PowerUtility is bst.PowerUtility
    old = bst.PowerUtility.price
    try:
        qs.PowerUtility.price = 0.123
        assert bst.PowerUtility.price == 0.123
    finally:
        bst.PowerUtility.price = old


def test_heatutility_agents_accessible_via_class():
    # HeatUtility classmethods/attrs work through the re-exported class.
    assert qs.HeatUtility is bst.HeatUtility
    assert qs.HeatUtility.get_agent('cooling_water') is bst.HeatUtility.get_agent('cooling_water')


# %% Lazy submodule re-export

def test_report_submodule_is_lazy_biosteam_report():
    assert qs.report is bst.report
    assert qs.report.voc_table is bst.report.voc_table
    assert qs.report.FOCTableBuilder is bst.report.FOCTableBuilder


# %% Raw BioSTEAM units surfaced through qsdsan.unit_operations.bst

def test_bst_namespace_reexports_raw_biosteam_units():
    from qsdsan.unit_operations import bst as su_bst

    assert su_bst.IsenthalpicValve is bst.IsenthalpicValve
    assert su_bst.Stripper is bst.Stripper
    assert su_bst.MolecularSieve is bst.MolecularSieve
    assert su_bst.BatchBioreactor is bst.BatchBioreactor
    assert su_bst.VacuumSystem is bst.VacuumSystem
    assert su_bst.Boiler is bst.Boiler
    assert su_bst.BoilerTurbogenerator is bst.BoilerTurbogenerator
    assert su_bst.ChilledWaterPackage is bst.ChilledWaterPackage
    assert su_bst.CoolingTower is bst.CoolingTower
    assert su_bst.SolidsCentrifuge is bst.SolidsCentrifuge


def test_raw_units_listed_in_namespace_all():
    from qsdsan.unit_operations import bst as su_bst

    for name in ('IsenthalpicValve', 'Stripper', 'MolecularSieve', 'BatchBioreactor',
                 'VacuumSystem', 'Boiler', 'BoilerTurbogenerator', 'ChilledWaterPackage',
                 'CoolingTower', 'SolidsCentrifuge'):
        assert name in su_bst.__all__


# %% Helper functions in qsdsan.utils

def test_utils_functional_helpers_match_thermosteam():
    assert qs.utils.rho_to_V is tmo.functional.rho_to_V
    assert qs.utils.V_to_rho is tmo.functional.V_to_rho


def test_utils_cost_decorator_matches_biosteam():
    # qsdsan thinly wraps BioSTEAM's @cost to also accept ``CEPCI`` as an alias for
    # ``CE``; ``__wrapped__`` still points at the BioSTEAM original, so an upstream
    # rename surfaces here (and at import time).
    from biosteam.units.decorators import cost as bst_cost
    assert qs.utils.cost.__wrapped__ is bst_cost


def test_utils_var_helpers_match_biosteam():
    from biosteam.evaluation._utils import var_columns, var_indices
    assert qs.utils.var_columns is var_columns
    assert qs.utils.var_indices is var_indices
