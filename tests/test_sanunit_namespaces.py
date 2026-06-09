# -*- coding: utf-8 -*-
"""Tests for public sanunit namespace groups."""

import qsdsan as qs


def test_sanunit_namespaces_are_public():
    from qsdsan.unit_operations import bst, static, dynamic

    assert qs.unit_operations.bst is bst
    assert qs.unit_operations.static is static
    assert qs.unit_operations.dynamic is dynamic


def test_bst_namespace_reexports_compatible_wrappers():
    from qsdsan.unit_operations import bst

    assert bst.Mixer is qs.unit_operations.Mixer
    assert bst.Splitter is qs.unit_operations.Splitter
    assert bst.Pump is qs.unit_operations.Pump
    assert bst.Flash is qs.unit_operations.Flash
    assert bst.BinaryDistillation is qs.unit_operations.BinaryDistillation
    assert bst.HXutility is qs.unit_operations.HXutility
    assert bst.StorageTank is qs.unit_operations.StorageTank
    assert bst.IsothermalCompressor is qs.unit_operations.IsothermalCompressor


def test_static_namespace_reexports_qsdsan_static_units():
    from qsdsan.unit_operations import static

    assert static.Excretion is qs.unit_operations.Excretion
    assert static.SepticTank is qs.unit_operations.SepticTank
    assert static.Sedimentation is qs.unit_operations.Sedimentation
    assert static.Screening is qs.unit_operations.Screening
    assert static.SludgePasteurization is qs.unit_operations.SludgePasteurization


def test_dynamic_namespace_reexports_dynamic_units():
    from qsdsan.unit_operations import dynamic

    assert dynamic.DynamicInfluent is qs.unit_operations.DynamicInfluent
    assert dynamic.HydraulicDelay is qs.unit_operations.HydraulicDelay
    assert dynamic.CSTR is qs.unit_operations.CSTR
    assert dynamic.PFR is qs.unit_operations.PFR
