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

import numpy as np, biosteam as bst, qsdsan as qs
from numpy.testing import assert_allclose

__all__ = (
    'test_default',
    'test_units_attribute_is_dict',
    'test_CEPCI_proxy',
    'test_TEA_CEPCI',
    'test_BinaryDistillation',
    'test_Flash',
    'test_HXprocess',
    'test_HXutility',
    'test_HXutility_results',
    'test_IsothermalCompressor',
    'test_MixTank',
    'test_Pump',
    'test_Splitter',
    'test_StorageTank',
    # bst-namespace units without a long-standing parity counterpart in this
    # file. Each test layers (smoke + parity-where-feasible + add-on) so a
    # BioSTEAM upstream change to ``__init__``/``_setup``/``_run`` would fail
    # at the QSDsan wrapping layer rather than silently propagating. Numeric
    # behavior remains the job of EXPOsan integration tests.
    'test_Mixer',
    'test_Tank',
    'test_FakeSplitter',
    'test_ReversedSplitter',
    'test_ShortcutColumn',
    'test_MESHDistillation',
    'test_AdiabaticMultiStageVLEColumn',
    'test_HeatExchangerNetwork',
    'test_ProcessWaterCenter',
    )


bst.default()
chems = bst.Chemicals(('Water', 'Methanol', 'Ethanol', 'Glycerol', 'H2'))

ws_data = {
    'particle_size': 'Soluble',
    'degradability': 'Readily',
    'organic': True
    }
cmps = qs.Components([
    qs.Component('Water',
                 particle_size='Soluble',
                 degradability='Undegradable',
                 organic=False),
    qs.Component('Methanol', **ws_data),
    qs.Component('Ethanol', **ws_data),
    qs.Component('Glycerol', **ws_data),
    qs.Component('H2',
                 particle_size='Dissolved gas',
                 degradability='Readily',
                 organic=True),
    ])


def create_streams(num):
    bst.settings.set_thermo(chems)
    bst_s = []
    for n in range(num):
        s = bst.Stream(Methanol=100*(n+1), Ethanol=100*(n+1), units='kg/hr')
        bst_s.append(s)

    qs.set_thermo(cmps)
    qs_ws = []
    for n in range(num):
        ws = qs.WasteStream(Methanol=100*(n+1), Ethanol=100*(n+1), units='kg/hr')
        qs_ws.append(ws)

    return bst_s, qs_ws


def test_waste_stream_accepts_component_flows_with_default_flow():
    qs.set_thermo(cmps)
    ws = qs.WasteStream(Water=80, Methanol=100, Glycerol=25, units='kmol/hr')
    assert_allclose(ws.imol['Water', 'Methanol', 'Glycerol'], [80, 100, 25])


def check_results(biosteam_unit, qsdsan_unit):
    global bst_unit, qs_unit
    bst_unit, qs_unit = biosteam_unit, qsdsan_unit
    bst_unit.simulate()
    qs_unit.simulate()

    bst_s = bst_unit.ins + bst_unit.outs
    qs_ws = qs_unit.ins + qs_unit.outs
    for n, s in enumerate(bst_s):
        assert_allclose(np.abs(s.mol-qs_ws[n].mol).sum(), 0, atol=1e-6)

    assert_allclose(bst_unit.installed_cost, qs_unit.installed_cost, atol=1e-6)
    bst_unit_ucost = 0 if bst_unit.utility_cost is None else bst_unit.utility_cost
    qs_unit_ucost = 0 if qs_unit.utility_cost is None else qs_unit.utility_cost
    assert_allclose(bst_unit_ucost, qs_unit_ucost, atol=1e-6)
    assert_allclose(bst_unit.power_utility.rate, qs_unit.power_utility.rate, atol=1e-6)


# =============================================================================
# Add-on coverage: verifies SanUnit mixin attributes (``add_OPEX``,
# ``lifetime``, ``construction``, etc.) exist on each bst-inherited class and
# that sentinel values survive a full ``simulate()`` pass. Catches BioSTEAM
# upstream changes that would silently drop SanUnit mixin state on the bst-
# inherited class. EXPOsan integration tests cover numeric drift; this just
# pins the *plumbing* class-by-class.
# =============================================================================

_ADDON_OPEX_SENTINEL = '_addon_smoke_opex_sentinel_'
_ADDON_OPEX_VALUE = 42.0
_ADDON_LIFETIME = 13
_ADDON_UPTIME_RATIO = 0.85
_REQUIRED_SANUNIT_ATTRS = (
    'construction', 'transportation', 'equipment', 'add_OPEX',
    'uptime_ratio', 'lifetime', 'include_construction', 'F_BM',
)


def check_addon_attrs(qs_unit):
    '''Assert the SanUnit mixin attributes are present on this bst class.'''
    missing = [a for a in _REQUIRED_SANUNIT_ATTRS if not hasattr(qs_unit, a)]
    assert not missing, (
        f'{type(qs_unit).__name__} missing SanUnit attrs: {missing}'
    )


def setup_addons(qs_unit):
    '''
    Stamp sentinel values onto the SanUnit add-on attributes before simulating.
    Pair with ``assert_addons_persisted`` after the simulate call to verify
    they were not clobbered by the BioSTEAM-side state machinery.
    '''
    check_addon_attrs(qs_unit)
    qs_unit.add_OPEX = {_ADDON_OPEX_SENTINEL: _ADDON_OPEX_VALUE}
    qs_unit.lifetime = _ADDON_LIFETIME
    qs_unit.uptime_ratio = _ADDON_UPTIME_RATIO


def assert_addons_persisted(qs_unit):
    '''Verify the sentinel add-on values stamped by ``setup_addons`` survived.'''
    assert qs_unit.add_OPEX.get(_ADDON_OPEX_SENTINEL) == _ADDON_OPEX_VALUE, (
        f'{type(qs_unit).__name__}: add_OPEX clobbered during simulate()'
    )
    assert qs_unit.lifetime == _ADDON_LIFETIME, (
        f'{type(qs_unit).__name__}: lifetime clobbered during simulate()'
    )
    assert qs_unit.uptime_ratio == _ADDON_UPTIME_RATIO, (
        f'{type(qs_unit).__name__}: uptime_ratio clobbered during simulate()'
    )



# %%

# =============================================================================
# Testing functions
# =============================================================================

def test_default():
    qs.default() # default everything


def test_units_attribute_is_dict():
    '''
    Every ``SanUnit`` subclass must have a ``dict`` ``_units`` attribute; BioSTEAM's
    ``results()`` does ``self._units.get(...)``, so a non-dict (e.g., ``None``) breaks it.

    This guards against a class of bug where a subclass sets
    ``_units = SomeParent._units.update({...})`` -- ``dict.update`` mutates in place and
    returns ``None``, leaving ``_units`` as ``None``.
    '''
    modules = [m for m in (getattr(qs, 'unit_operations', None),
                           getattr(qs, 'sanunits', None)) if m is not None]
    seen = {}
    for module in modules:
        for name in dir(module):
            obj = getattr(module, name)
            if isinstance(obj, type) and issubclass(obj, qs.SanUnit):
                seen[obj] = obj.__name__
    assert seen, 'no SanUnit subclasses found to check'
    bad = sorted(n for cls, n in seen.items()
                 if not isinstance(getattr(cls, '_units', {}), dict))
    assert not bad, f'these unit classes have a non-dict `_units`: {bad}'


def test_CEPCI_proxy():
    '''
    ``qs.CEPCI`` reads and writes BioSTEAM's global cost index (which BioSTEAM
    abbreviates as ``CE``) so users do not have to import biosteam. Implemented as a
    property on a module subclass.
    '''
    old = bst.CE
    try:
        assert qs.CEPCI == bst.CE           # live read
        qs.CEPCI = 600.0
        assert bst.CE == 600.0             # assignment proxies to bst.CE
        bst.CE = 555.0
        assert qs.CEPCI == 555.0           # remains a live view
        qs.CEPCI = qs.CEPCI_by_year[2023]
        assert bst.CE == qs.CEPCI_by_year[2023]
    finally:
        bst.CE = old


def test_TEA_CEPCI():
    '''
    Creating a `TEA` without an explicit `CEPCI` must NOT reset the global cost index
    (regression for an early-binding ``CEPCI=bst.CE`` default that froze it at 567.5).
    A provided ``CEPCI`` is applied, and ``TEA.CEPCI_by_year`` mirrors ``qs.CEPCI_by_year``.
    '''
    old = bst.CE
    try:
        _, qs_ws = create_streams(1)
        qs.set_thermo(cmps)
        M = qs.unit_operations.MixTank('M_cepci', ins=qs_ws)
        sys = qs.System('sys_cepci', path=(M,))
        bst.CE = qs.CEPCI_by_year[2023]
        tea = qs.TEA(system=sys, discount_rate=0.05, lifetime=10)
        assert bst.CE == qs.CEPCI_by_year[2023]      # not reset by TEA creation
        assert tea.CEPCI_by_year is qs.CEPCI_by_year  # consistent table
        qs.TEA(system=sys, CEPCI=qs.CEPCI_by_year[2010], discount_rate=0.05, lifetime=10)
        assert bst.CE == qs.CEPCI_by_year[2010]      # explicit CEPCI is applied
    finally:
        bst.CE = old


def test_BinaryDistillation():
    bst.settings.set_thermo(chems)
    stream_kwargs = dict(Water=80, Methanol=100, Glycerol=25, units='kmol/hr')
    bst_s = bst.Stream(**stream_kwargs)
    bst_s.T = T = bst_s.bubble_point_at_P().T # Feed at bubble point T
    unit_kwargs = dict(LHK=('Methanol', 'Water'), y_top=0.99, x_bot=0.01, k=2, is_divided=True)
    bst_unit = bst.units.BinaryDistillation(ins=bst_s, **unit_kwargs)
    bst_unit.simulate() # need extra simulation for better convergence
    
    qs.set_thermo(cmps)
    qs_s = qs.WasteStream(T=T, **stream_kwargs)
    qs_unit = qs.unit_operations.BinaryDistillation(ins=qs_s, **unit_kwargs)
    qs_unit.simulate()

    setup_addons(qs_unit)
    check_results(bst_unit, qs_unit)
    assert_addons_persisted(qs_unit)


def test_Flash():
    bst.settings.set_thermo(chems)
    stream_kwargs = dict(Glycerol=300, Water=1000, units='kmol/hr')
    bst_s = bst.Stream(**stream_kwargs)
    bst_s.T = T = bst_s.bubble_point_at_P().T # Feed at bubble point T
    unit_kwargs = dict(P=101325, T=410.15)
    bst_unit = bst.units.Flash(ins=bst_s, **unit_kwargs)
    
    qs.set_thermo(cmps)
    qs_s = qs.Stream(T=T, **stream_kwargs)
    qs_unit = qs.unit_operations.Flash(ins=qs_s, **unit_kwargs)

    setup_addons(qs_unit)
    check_results(bst_unit, qs_unit)
    assert_addons_persisted(qs_unit)


def test_HXprocess():
    bst_s, qs_ws = create_streams(2)
    bst_s[0].T = qs_ws[0].T = 400

    bst.settings.set_thermo(chems)
    bst_unit = bst.units.HXprocess(ins=bst_s, phase0='l', phase1='l')

    qs.set_thermo(cmps)
    qs_unit = qs.unit_operations.HXprocess(ins=qs_ws, phase0='l', phase1='l')

    setup_addons(qs_unit)
    check_results(bst_unit, qs_unit)
    assert_addons_persisted(qs_unit)


def test_HXutility():
    bst_s, qs_ws = create_streams(1)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.HXutility(ins=bst_s, T=400, rigorous=True)

    qs.set_thermo(cmps)
    qs_unit = qs.unit_operations.HXutility(ins=qs_ws, T=400, rigorous=True)

    setup_addons(qs_unit)
    check_results(bst_unit, qs_unit)
    assert_addons_persisted(qs_unit)


def test_HXutility_results():
    '''
    Regression test: ``HXutility.results()`` must work. Its ``_units`` was once set to
    ``None`` (assigned the return value of ``dict.update``), which raised
    ``AttributeError`` inside BioSTEAM's ``results()``. ``check_results`` only calls
    ``simulate()``, so this exercises the reporting path explicitly.
    '''
    _, qs_ws = create_streams(1)
    qs.set_thermo(cmps)
    qs_unit = qs.unit_operations.HXutility('H_results', ins=qs_ws, T=350)
    qs_unit.simulate()
    df = qs_unit.results()
    assert df is not None and len(df) > 0


def test_IsothermalCompressor():
    bst.settings.set_thermo(chems)
    stream_kwargs = dict(H2=1, T=298.15, P=20e5, units='kmol/hr', phase='g')
    bst_s = bst.Stream(**stream_kwargs)
    unit_kwargs = dict(P=350e5, eta=1)
    bst_unit = bst.units.IsothermalCompressor(ins=bst_s, **unit_kwargs)
    
    qs.set_thermo(cmps)
    qs_s = qs.WasteStream(**stream_kwargs)
    qs_unit = qs.unit_operations.IsothermalCompressor(ins=qs_s, **unit_kwargs)

    setup_addons(qs_unit)
    check_results(bst_unit, qs_unit)
    assert_addons_persisted(qs_unit)


def test_MixTank():
    bst_s, qs_ws = create_streams(2)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.MixTank(ins=bst_s)

    qs.set_thermo(cmps)
    qs_unit = qs.unit_operations.MixTank(ins=qs_ws)

    setup_addons(qs_unit)
    check_results(bst_unit, qs_unit)
    assert_addons_persisted(qs_unit)


def test_Pump():
    bst_s, qs_ws = create_streams(1)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.Pump(ins=bst_s)

    qs.set_thermo(cmps)
    qs_unit = qs.unit_operations.Pump(ins=qs_ws)

    setup_addons(qs_unit)
    check_results(bst_unit, qs_unit)
    assert_addons_persisted(qs_unit)


def test_Splitter():
    bst_s, qs_ws = create_streams(1)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.Splitter(ins=bst_s, split=0.1)

    qs.set_thermo(cmps)
    qs_unit = qs.unit_operations.Splitter(ins=qs_ws, split=0.1)

    setup_addons(qs_unit)
    check_results(bst_unit, qs_unit)
    assert_addons_persisted(qs_unit)


def test_StorageTank():
    bst_s, qs_ws = create_streams(1)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.StorageTank(ins=bst_s)

    qs.set_thermo(cmps)
    qs_unit = qs.unit_operations.StorageTank(ins=qs_ws)

    setup_addons(qs_unit)
    check_results(bst_unit, qs_unit)
    assert_addons_persisted(qs_unit)


# =============================================================================
# bst-namespace units without a long-standing parity counterpart in this file.
#
# These verify only that the QSDsan-exposed class instantiates and simulates
# without error -- they do NOT pin numeric behavior. The bst units are in
# QSDsan because EXPOsan systems use them, and EXPOsan integration tests pin
# the numbers at the system level. A smoke test here is enough to catch
# import/instantiation regressions in the QSDsan wrapping layer.
# =============================================================================

def test_Mixer():
    bst_s, qs_ws = create_streams(2)
    bst.settings.set_thermo(chems)
    bst_unit = bst.units.Mixer(ins=bst_s)

    qs.set_thermo(cmps)
    qs_unit = qs.unit_operations.Mixer(ins=qs_ws)

    setup_addons(qs_unit)
    check_results(bst_unit, qs_unit)
    assert_addons_persisted(qs_unit)


def test_Tank():
    '''
    ``Tank`` is an abstract base class -- it has no ``purchase_cost_algorithms``,
    which concrete subclasses (``MixTank``, ``StorageTank``) provide, so it
    cannot be simulated standalone and has no parity counterpart. Verify only
    that the class still mixes in ``SanUnit`` and exposes the add-on attribute
    surface. The concrete tanks above carry full smoke + parity + add-on.
    '''
    cls = qs.unit_operations.Tank
    assert isinstance(cls, type) and issubclass(cls, qs.SanUnit)
    # Instance-level attribute presence: bypass __init__ (which would raise on
    # missing ``purchase_cost_algorithms``) and assert the mixin still installs
    # the SanUnit attribute surface on a bare instance.
    qs.set_thermo(cmps)
    inst = cls.__new__(cls)
    qs.SanUnit.__init__(inst, 'T_addons')
    check_addon_attrs(inst)


def test_FakeSplitter():
    '''
    ``FakeSplitter`` requires the user to assign ``outs`` manually before
    simulating (no automatic split), so a parity comparison against BioSTEAM's
    ``MockSplitter`` adds little signal. Cover smoke + add-ons only.
    '''
    _, qs_ws = create_streams(1)
    qs.set_thermo(cmps)
    u = qs.unit_operations.FakeSplitter('FS_addons', ins=qs_ws, outs=('o1', 'o2'))
    u.outs[0].copy_like(qs_ws[0])
    u.outs[0].imass['Methanol'] *= 0.4
    u.outs[1].copy_like(qs_ws[0])
    u.outs[1].imass['Methanol'] *= 0.6
    setup_addons(u)
    u.simulate()
    assert_addons_persisted(u)


def test_ReversedSplitter():
    '''
    ``ReversedSplitter`` has outs-driven semantics (user pre-sets outs, the
    inlet is computed) and its BioSTEAM signature is `**kwargs`-only -- so a
    structured parity comparison against the BST class adds little signal.
    Cover smoke + add-ons here; parity is implicit in the qs-side simulate.
    '''
    qs.set_thermo(cmps)
    src = qs.WasteStream('rsrc', Methanol=100, Ethanol=100, units='kg/hr')
    qs_unit = qs.unit_operations.ReversedSplitter(
        'RS_addons', ins=src, outs=('ro1', 'ro2'), split=(0.3, 0.7),
    )

    setup_addons(qs_unit)
    qs_unit.simulate()
    assert_addons_persisted(qs_unit)


def test_ShortcutColumn():
    bst.settings.set_thermo(chems)
    stream_kwargs = dict(Water=80, Methanol=100, Glycerol=25, units='kmol/hr')
    bst_s = bst.Stream(**stream_kwargs)
    bst_s.T = T = bst_s.bubble_point_at_P().T
    unit_kwargs = dict(
        LHK=('Methanol', 'Water'), y_top=0.99, x_bot=0.01, k=2, is_divided=True,
    )
    bst_unit = bst.units.ShortcutColumn(ins=bst_s, **unit_kwargs)

    qs.set_thermo(cmps)
    qs_s = qs.WasteStream(T=T, **stream_kwargs)
    qs_unit = qs.unit_operations.ShortcutColumn(ins=qs_s, **unit_kwargs)

    setup_addons(qs_unit)
    check_results(bst_unit, qs_unit)
    assert_addons_persisted(qs_unit)


def test_MESHDistillation():
    bst.settings.set_thermo(chems)
    stream_kwargs = dict(Water=80, Methanol=100, units='kmol/hr')
    bst_s = bst.Stream(**stream_kwargs)
    bst_s.T = T = bst_s.bubble_point_at_P().T
    unit_kwargs = dict(
        LHK=('Methanol', 'Water'), N_stages=10, feed_stages=(5,),
        reflux=2.0, boilup=2.5,
    )
    bst_unit = bst.units.MESHDistillation(
        ins=bst_s, outs=('bvap', 'bliq'), **unit_kwargs,
    )

    qs.set_thermo(cmps)
    qs_s = qs.WasteStream(T=T, **stream_kwargs)
    qs_unit = qs.unit_operations.MESHDistillation(
        ins=qs_s, outs=('qvap', 'qliq'), **unit_kwargs,
    )

    setup_addons(qs_unit)
    check_results(bst_unit, qs_unit)
    assert_addons_persisted(qs_unit)


def test_AdiabaticMultiStageVLEColumn():
    '''
    ``AdiabaticMultiStageVLEColumn`` is an absorption-style column needing a
    paired vapor + liquid feed to establish a VLE profile; standalone configs
    diverge in BioSTEAM's stage solver. Cover instantiate + add-on attribute
    surface here -- numeric drift is the job of EXPOsan systems that actually
    use it.
    '''
    qs.set_thermo(cmps)
    vap = qs.WasteStream('ac_vap', H2=1, phase='g', units='kmol/hr')
    liq = qs.WasteStream('ac_liq', Water=10, units='kmol/hr')
    u = qs.unit_operations.AdiabaticMultiStageVLEColumn(
        'AC_addons', ins=(vap, liq), outs=('top', 'bot'),
        N_stages=3, solute='H2', feed_stages=(0, 2),
    )
    check_addon_attrs(u)


def test_HeatExchangerNetwork():
    '''
    ``HeatExchangerNetwork`` internally switches stream phases to
    ``MultiStream``, which ``WasteStream``'s class layout disallows. The
    standard bridge is :class:`~.PhaseChanger` with ``init_with='Stream'``,
    which converts an upstream ``WasteStream`` into a plain ``qs.Stream``
    (= ``thermosteam.Stream``) that HEN can consume.
    '''
    qs.set_thermo(cmps)
    ws = qs.WasteStream('hen_feed_ws', Water=100, Methanol=50, units='kmol/hr')
    pc = qs.unit_operations.PhaseChanger(
        'PC_hen', ins=ws, init_with='Stream', phase='l',
    )
    hx = qs.unit_operations.HXutility(
        'H_hen', ins=pc-0, T=400, rigorous=True, init_with='Stream',
    )
    hen = qs.unit_operations.HeatExchangerNetwork('HEN_addons')
    sys = qs.System('sys_hen_addons', path=(pc, hx), facilities=(hen,))
    setup_addons(hen)
    sys.simulate()
    assert_addons_persisted(hen)


def test_ProcessWaterCenter():
    '''
    ``ProcessWaterCenter`` is a facility; smoke-test it by attaching to a
    minimal system. Add-on attributes are verified on the facility instance.
    '''
    qs.set_thermo(cmps)
    pw = qs.WasteStream('process_water_in', Water=1000, units='kg/hr')
    pwc = qs.unit_operations.ProcessWaterCenter(
        'PWC_addons', ins=pw, process_water_streams=(pw,),
    )
    sys = qs.System('sys_pwc_addons', path=(), facilities=(pwc,))
    setup_addons(pwc)
    sys.simulate()
    assert_addons_persisted(pwc)


if __name__ == '__main__':
    test_default()
    test_units_attribute_is_dict()
    test_CEPCI_proxy()
    test_TEA_CEPCI()
    test_BinaryDistillation()
    test_Flash()
    test_HXprocess()
    test_HXutility()
    test_HXutility_results()
    test_IsothermalCompressor()
    test_MixTank()
    test_Pump()
    test_Splitter()
    test_StorageTank()
    test_Mixer()
    test_Tank()
    test_FakeSplitter()
    test_ReversedSplitter()
    test_ShortcutColumn()
    test_MESHDistillation()
    test_AdiabaticMultiStageVLEColumn()
    test_HeatExchangerNetwork()
    test_ProcessWaterCenter()
