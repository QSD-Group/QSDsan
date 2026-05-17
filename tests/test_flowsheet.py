#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ('test_flowsheet',)


def test_flowsheet():
    import qsdsan as qs
    from qsdsan._sanflowsheet import SanFlowsheet, SanMainFlowsheet

    # ── 1. SanFlowsheet has all four LCA registries ───────────────────────────
    fs = qs.Flowsheet('_test_lca_1')
    assert hasattr(fs, 'indicator'),      "flowsheet missing 'indicator' registry"
    assert hasattr(fs, 'item'),           "flowsheet missing 'item' registry"
    assert hasattr(fs, 'construction'),   "flowsheet missing 'construction' registry"
    assert hasattr(fs, 'transportation'), "flowsheet missing 'transportation' registry"

    # ── 2. qs.Flowsheet IS SanFlowsheet ──────────────────────────────────────
    assert qs.Flowsheet is SanFlowsheet, "qs.Flowsheet should be SanFlowsheet"

    # ── 3. qs.main_flowsheet IS SanMainFlowsheet ─────────────────────────────
    assert isinstance(qs.main_flowsheet, SanMainFlowsheet), \
        "qs.main_flowsheet should be a SanMainFlowsheet instance"

    # ── 4. Item registry is isolated per flowsheet ───────────────────────────
    with qs.Flowsheet('_test_sys_a') as fs_a:
        GWP_a = qs.ImpactIndicator('GlobalWarming_t', unit='kg CO2-eq')
        steel_a = qs.ImpactItem('Steel_t', 'kg', GlobalWarming_t=2.55)
        assert qs.ImpactItem.get_item('Steel_t') is steel_a

    # After exiting sys_a, we are back to the previous flowsheet
    assert qs.ImpactItem.get_item('Steel_t') is None, \
        "Steel_t should not be visible outside its flowsheet"

    with qs.Flowsheet('_test_sys_b') as fs_b:
        assert qs.ImpactItem.get_item('Steel_t') is None, \
            "Steel_t from sys_a should not leak into sys_b"
        GWP_b = qs.ImpactIndicator('GlobalWarming_t', unit='kg CO2-eq')
        steel_b = qs.ImpactItem('Steel_t', 'kg', GlobalWarming_t=3.0)
        assert qs.ImpactItem.get_item('Steel_t') is steel_b

    # ── 5. ImpactIndicators are isolated per flowsheet ───────────────────────
    with qs.Flowsheet('_test_sys_c') as fs_c:
        assert qs.ImpactIndicator.get_indicator('GlobalWarming_t') is None, \
            "ImpactIndicator from another flowsheet should not be visible"

    assert qs.ImpactIndicator.get_indicator('GlobalWarming_t') is None, \
        "ImpactIndicator should not persist after leaving its flowsheet"

    # ── 6. clear() detaches StreamImpactItem from its stream ─────────────────
    components = qs.Components.load_default()
    qs.set_thermo(components)

    with qs.Flowsheet('_test_sys_d') as fs_d:
        GWP_d = qs.ImpactIndicator('GlobalWarming_t', unit='kg CO2-eq')
        s = qs.SanStream('_test_stream_d', H2O=100, units='kg/hr')
        item_d = qs.StreamImpactItem(linked_stream=s, GlobalWarming_t=1.5)
        assert s.stream_impact_item is item_d
        fs_d.clear()
        assert s.stream_impact_item is None, \
            "clear() should unlink StreamImpactItem from stream"
        assert qs.ImpactItem.get_item(item_d.ID) is None, \
            "clear() should remove item from registry"

    # ── 7. clear() detaches Construction from its unit ───────────────────────
    with qs.Flowsheet('_test_sys_e') as fs_e:
        GWP_e = qs.ImpactIndicator('GlobalWarming_t', unit='kg CO2-eq')
        concrete = qs.ImpactItem('Concrete_t', 'kg', GlobalWarming_t=0.1)
        M = qs.unit_operations.MixTank(
            '_test_M_e', ins=qs.WasteStream(H2O=100, units='kg/hr')
        )
        c = qs.Construction(item=concrete, quantity=50, quantity_unit='kg',
                            linked_unit=M)
        M.construction = [c]
        assert len(M.construction) == 1
        fs_e.clear()
        assert M.construction == [], \
            "clear() should empty unit.construction"
        assert qs.ImpactItem.get_item('Concrete_t') is None


def test_construction_specs():
    import qsdsan as qs
    from qsdsan.unit_operations import MixTank

    components = qs.Components.load_default()
    qs.set_thermo(components)

    # ── Define a unit subclass with default construction specs ────────────────
    class SpecTank(MixTank):
        # Each spec: dict with keys item, quantity, quantity_unit,
        # and optionally lifetime, lifetime_unit.
        _construction_specs = (
            dict(item='SpecSteel', quantity=10., quantity_unit='kg',
                 lifetime=20., lifetime_unit='yr'),
            dict(item='SpecConcrete', quantity=5., quantity_unit='m3'),
        )

    # ── 1. LCA raises clearly when required items are not loaded ─────────────
    with qs.Flowsheet('_test_spec_a') as fs_a:
        feed = qs.WasteStream('_spec_feed_a', H2O=1000, units='kg/hr')
        tank = SpecTank('_spec_tank_a', ins=feed)
        sys = qs.System('_spec_sys_a', path=(tank,))
        sys.simulate()

        GWP = qs.ImpactIndicator('GWP_spec', unit='kg CO2-eq')
        # Item NOT loaded — LCA should raise
        import pytest
        with pytest.raises(RuntimeError, match="SpecSteel"):
            lca = qs.LCA(sys, lifetime=20, indicators=(GWP,),
                         simulate_system=False)

    # ── 2. LCA resolves specs when items are loaded ───────────────────────────
    with qs.Flowsheet('_test_spec_b') as fs_b:
        feed = qs.WasteStream('_spec_feed_b', H2O=1000, units='kg/hr')
        tank = SpecTank('_spec_tank_b', ins=feed)
        sys = qs.System('_spec_sys_b', path=(tank,))
        sys.simulate()

        GWP     = qs.ImpactIndicator('GWP_spec', unit='kg CO2-eq')
        steel   = qs.ImpactItem('SpecSteel',    'kg', GWP_spec=2.55)
        concrete = qs.ImpactItem('SpecConcrete', 'm3', GWP_spec=100.)

        lca = qs.LCA(sys, lifetime=20, indicators=(GWP,),
                     simulate_system=False)

        # Spec-derived constructions appear in LCA construction units
        assert tank in lca.construction_units, \
            "unit with _construction_specs should appear in lca.construction_units"

        impacts = lca.get_construction_impacts()
        # 10 kg steel × 2.55 kg CO2-eq/kg = 25.5; 5 m3 concrete × 100 = 500
        from numpy.testing import assert_allclose
        assert_allclose(impacts['GWP_spec'], 25.5 + 500., rtol=1e-6)

    # ── 3. Explicit unit.construction overrides spec for same item ────────────
    with qs.Flowsheet('_test_spec_c') as fs_c:
        feed = qs.WasteStream('_spec_feed_c', H2O=1000, units='kg/hr')
        tank = SpecTank('_spec_tank_c', ins=feed)

        GWP    = qs.ImpactIndicator('GWP_spec', unit='kg CO2-eq')
        steel  = qs.ImpactItem('SpecSteel',    'kg', GWP_spec=2.55)
        conc   = qs.ImpactItem('SpecConcrete', 'm3', GWP_spec=100.)

        # User provides explicit construction for SpecSteel → spec is skipped
        explicit = qs.Construction(item=steel, quantity=99., quantity_unit='kg',
                                   linked_unit=tank)
        tank.construction = [explicit]

        sys = qs.System('_spec_sys_c', path=(tank,))
        sys.simulate()
        lca = qs.LCA(sys, lifetime=20, indicators=(GWP,), simulate_system=False)

        impacts = lca.get_construction_impacts()
        # 99 kg steel × 2.55 = 252.45; 5 m3 concrete × 100 = 500
        assert_allclose(impacts['GWP_spec'], 252.45 + 500., rtol=1e-6)

    # ── 4. _construction_specs defaults to empty tuple on base SanUnit ────────
    assert qs.SanUnit._construction_specs == (), \
        "SanUnit._construction_specs should default to empty tuple"


def test_two_system_switch():
    """Verify that switching between two LCA-enabled systems produces no
    cross-contamination in registries or unit state."""
    import qsdsan as qs
    from numpy.testing import assert_allclose

    components = qs.Components.load_default()
    qs.set_thermo(components)

    # ── System A ──────────────────────────────────────────────────────────────
    with qs.Flowsheet('_switch_sys_a') as fs_a:
        GWP = qs.ImpactIndicator('GWP_switch', unit='kg CO2-eq')

        feed_a = qs.WasteStream('_sw_feed_a', H2O=1000, units='kg/hr')
        M_a    = qs.unit_operations.MixTank('_sw_M_a', ins=feed_a)
        sys_a  = qs.System('_sw_sys_a', path=(M_a,))
        sys_a.simulate()

        steel_a  = qs.ImpactItem('SwitchSteel', 'kg', GWP_switch=2.0)
        constr_a = qs.Construction(item=steel_a, quantity=10., quantity_unit='kg',
                                   linked_unit=M_a)
        M_a.construction = [constr_a]

        lca_a = qs.LCA(sys_a, lifetime=10, indicators=(GWP,),
                       simulate_system=False)
        impacts_a = lca_a.get_construction_impacts()
        # 10 kg × 2.0 = 20
        assert_allclose(impacts_a['GWP_switch'], 20., rtol=1e-6)

    # ── After exiting sys_a: SwitchSteel should not be visible ───────────────
    assert qs.ImpactItem.get_item('SwitchSteel') is None, \
        "SwitchSteel from sys_a should not be visible after exiting its flowsheet"

    # ── System B (same item name, different CF) ───────────────────────────────
    with qs.Flowsheet('_switch_sys_b') as fs_b:
        GWP = qs.ImpactIndicator('GWP_switch', unit='kg CO2-eq')

        assert qs.ImpactItem.get_item('SwitchSteel') is None, \
            "SwitchSteel should not leak from sys_a into sys_b"

        feed_b = qs.WasteStream('_sw_feed_b', H2O=500, units='kg/hr')
        M_b    = qs.unit_operations.MixTank('_sw_M_b', ins=feed_b)
        sys_b  = qs.System('_sw_sys_b', path=(M_b,))
        sys_b.simulate()

        # Same name, different CF — no conflict warning expected
        steel_b  = qs.ImpactItem('SwitchSteel', 'kg', GWP_switch=5.0)
        constr_b = qs.Construction(item=steel_b, quantity=3., quantity_unit='kg',
                                   linked_unit=M_b)
        M_b.construction = [constr_b]

        lca_b = qs.LCA(sys_b, lifetime=10, indicators=(GWP,),
                       simulate_system=False)
        impacts_b = lca_b.get_construction_impacts()
        # 3 kg × 5.0 = 15
        assert_allclose(impacts_b['GWP_switch'], 15., rtol=1e-6)

    # ── sys_a impacts are unchanged ───────────────────────────────────────────
    with qs.Flowsheet('_switch_sys_a'):
        assert_allclose(
            lca_a.get_construction_impacts()['GWP_switch'], 20., rtol=1e-6,
            err_msg="sys_a LCA impacts should be unchanged after sys_b was created"
        )


if __name__ == '__main__':
    test_flowsheet()
    test_construction_specs()
    test_two_system_switch()
