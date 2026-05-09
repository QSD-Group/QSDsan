# -*- coding: utf-8 -*-
"""Tests for public sanunit namespace groups."""

import qsdsan as qs


def test_sanunit_namespaces_are_public():
    from qsdsan.sanunits import bst, static, dynamic

    assert qs.sanunits.bst is bst
    assert qs.sanunits.static is static
    assert qs.sanunits.dynamic is dynamic


def test_bst_namespace_reexports_compatible_wrappers():
    from qsdsan.sanunits import bst

    assert bst.Mixer is qs.sanunits.Mixer
    assert bst.Splitter is qs.sanunits.Splitter
    assert bst.Pump is qs.sanunits.Pump
    assert bst.Flash is qs.sanunits.Flash
    assert bst.BinaryDistillation is qs.sanunits.BinaryDistillation
    assert bst.HXutility is qs.sanunits.HXutility
    assert bst.StorageTank is qs.sanunits.StorageTank
    assert bst.IsothermalCompressor is qs.sanunits.IsothermalCompressor


def test_static_namespace_reexports_qsdsan_static_units():
    from qsdsan.sanunits import static

    assert static.Excretion is qs.sanunits.Excretion
    assert static.SepticTank is qs.sanunits.SepticTank
    assert static.Sedimentation is qs.sanunits.Sedimentation
    assert static.Screening is qs.sanunits.Screening
    assert static.SludgePasteurization is qs.sanunits.SludgePasteurization


def test_dynamic_namespace_reexports_dynamic_units():
    from qsdsan.sanunits import dynamic

    assert dynamic.DynamicInfluent is qs.sanunits.DynamicInfluent
    assert dynamic.HydraulicDelay is qs.sanunits.HydraulicDelay
    assert dynamic.CSTR is qs.sanunits.CSTR
    assert dynamic.PFR is qs.sanunits.PFR
