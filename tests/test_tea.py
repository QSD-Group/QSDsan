#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = (
    'test_tea_metric_relationships',
    'test_capital_hierarchy_collapsed_by_default',
    'test_capital_hierarchy_with_subclass_markups',
    'test_FOC_components',
    'test_material_cost_and_sales_from_stream_prices',
    'test_annualized_equipment_cost_equals_capex_without_replacement',
    'test_annualized_equipment_cost_int_lifetime',
    'test_annualized_equipment_cost_dict_lifetime',
    'test_depreciation_no_effect_without_taxable_income',
    'test_depreciation_changes_npv_with_tax_and_profit',
    'test_depreciation_schedule_must_fit_lifetime',
    'test_CEPCI_accessors',
    'test_simpletea_is_deprecated_alias',
    'test_agile_system_npv_charges_equipment_replacement',
    )

import pytest


def _A(n, r):
    '''Annuity (present-worth) factor.'''
    return (1 - (1 + r) ** -n) / r


def _example_tea(lifetime=10, **kwargs):
    import qsdsan as qs
    from qsdsan.utils import create_example_system
    sys = create_example_system()
    sys.simulate()
    return qs.TEA(sys, discount_rate=0.05, lifetime=lifetime, simulate_system=False, **kwargs), sys


def _two_equipment_system(flowsheet_ID):
    '''A controlled 1-in/1-out unit with two equipment items (no numba / heat exchange).'''
    import qsdsan as qs
    from qsdsan import SanUnit, WasteStream, Component, Components
    cmps = Components([Component('Water', search_ID='Water', particle_size='Soluble',
                                 degradability='Undegradable', organic=False)])
    cmps.compile()
    qs.set_thermo(cmps)
    qs.main_flowsheet.set_flowsheet(qs.Flowsheet(flowsheet_ID))

    class TwoEquip(SanUnit):
        _N_ins = 1
        _N_outs = 1
        def _run(self):
            self.outs[0].copy_like(self.ins[0])
        def _cost(self):
            self.baseline_purchase_costs['A'] = 1e6   # F_BM defaults to 1 -> installed == purchase
            self.baseline_purchase_costs['B'] = 2e6

    ws = WasteStream('ws', Water=1000)
    u = TwoEquip('u', ins=ws, outs='e')
    sys = qs.System('sys', path=(u,))
    sys.simulate()
    return sys, u


def _example_agile_system(flowsheet_ID, equipment_lifetime=0):
    '''The two-equipment system wrapped in an AgileSystem of two operation modes.'''
    import qsdsan as qs
    sys, u = _two_equipment_system(flowsheet_ID)
    u.equipment_lifetime = equipment_lifetime
    sys.simulate()
    ag = qs.AgileSystem()

    @ag.operation_parameter
    def set_flow(flow):
        u.ins[0].imass['Water'] = flow

    ag.operation_mode(sys, operating_hours=4000., flow=1000.)
    ag.operation_mode(sys, operating_hours=4760., flow=500.)
    ag.simulate()
    return ag


# %%

def test_tea_metric_relationships():
    '''AOC/VOC/EAC/CAPEX identities and IRR<->discount_rate alias.'''
    tea, sys = _example_tea()
    assert tea.IRR == tea.discount_rate == 0.05
    assert tea.VOC == pytest.approx(tea.material_cost + tea.utility_cost)
    assert tea.AOC == pytest.approx(tea.FOC + tea.VOC)
    assert tea.EAC == pytest.approx(tea.annualized_CAPEX + tea.AOC)
    assert tea.CAPEX == pytest.approx(tea.TCI)
    # a pure-cost system has a negative NPV
    assert tea.NPV < 0


def test_capital_hierarchy_collapsed_by_default():
    '''With no markups and no working capital, the whole capital chain is one number.'''
    tea, sys = _example_tea()
    ic = tea.installed_equipment_cost
    assert ic > 0
    for attr in ('DPI', 'TDC', 'FCI', 'TCI', 'CAPEX'):
        assert getattr(tea, attr) == pytest.approx(ic)


def test_capital_hierarchy_with_subclass_markups():
    '''A subclass can reopen the collapsed hierarchy by overriding the build-up methods.'''
    import qsdsan as qs
    from qsdsan.utils import create_example_system
    sys = create_example_system(); sys.simulate()

    class CustomTEA(qs.TEA):
        def _DPI(self, installed_equipment_cost):
            return installed_equipment_cost * 1.30
        def _TDC(self, DPI):
            return DPI * 1.10

    tea = CustomTEA(sys, discount_rate=0.05, lifetime=10, simulate_system=False)
    ic = tea.installed_equipment_cost
    assert tea.DPI == pytest.approx(ic * 1.30)
    assert tea.TDC == pytest.approx(ic * 1.30 * 1.10)
    assert tea.FCI == pytest.approx(tea.TDC)            # _FCI defaults to TDC
    assert tea.CAPEX == pytest.approx(tea.TCI)


def test_FOC_components():
    '''FOC = FCI*annual_maintenance + annual_labor + total_add_OPEX.'''
    tea, sys = _example_tea(annual_maintenance=0.02, annual_labor=1e5)
    expected = tea.FCI * 0.02 + 1e5 + tea.total_add_OPEX
    assert tea.FOC == pytest.approx(expected)
    # maintenance is the only FOC term tied to capital
    tea.annual_maintenance = 0.
    tea.annual_labor = 0.
    assert tea.FOC == pytest.approx(tea.total_add_OPEX)


def test_material_cost_and_sales_from_stream_prices():
    '''Priced feeds -> material_cost (VOC); priced products -> sales.'''
    import qsdsan as qs
    from qsdsan.utils import create_example_system
    sys = create_example_system(); sys.simulate()
    feed = [s for s in sys.feeds if s.ID == 'methanol'][0]
    product = [s for s in sys.products if s.ID == 'alcohols'][0]
    feed.price = 0.5
    product.price = 2.0
    tea = qs.TEA(sys, discount_rate=0.05, lifetime=10, simulate_system=False)
    assert tea.material_cost > 0
    assert tea.sales > 0
    assert tea.VOC == pytest.approx(tea.material_cost + tea.utility_cost)


def test_annualized_equipment_cost_equals_capex_without_replacement():
    '''Equipment assumed to last the project -> the two annualized metrics coincide.'''
    tea, sys = _example_tea(lifetime=10)
    assert tea.annualized_equipment_cost == pytest.approx(tea.annualized_CAPEX)


def test_annualized_equipment_cost_int_lifetime():
    '''An integer equipment lifetime annualizes the installed cost over that lifetime.'''
    sys, u = _two_equipment_system('test_tea_int_life')
    import qsdsan as qs
    u.equipment_lifetime = 8
    sys.simulate()
    tea = qs.TEA(sys, discount_rate=0.05, lifetime=20, simulate_system=False)
    installed = u.installed_cost
    assert tea.annualized_equipment_cost == pytest.approx(installed / _A(8, 0.05))
    # a replacement is charged in the cash flow, so annualized_CAPEX differs (is higher)
    assert tea.annualized_CAPEX > tea.annualized_equipment_cost


def test_annualized_equipment_cost_dict_lifetime():
    '''Per-equipment lifetime dict annualizes each item over its own lifetime.

    Regression test for the variable-shadowing bug in
    ``get_unit_annualized_equipment_cost``'s dict branch.
    '''
    import qsdsan as qs
    sys, u = _two_equipment_system('test_tea_dict_life')
    r = 0.05
    A_cost, B_cost = u.installed_costs['A'], u.installed_costs['B']

    # different lifetimes per equipment
    u.equipment_lifetime = {'A': 10, 'B': 20}
    sys.simulate()
    tea = qs.TEA(sys, discount_rate=r, lifetime=20, simulate_system=False)
    assert tea.annualized_equipment_cost == pytest.approx(
        A_cost / _A(10, r) + B_cost / _A(20, r))

    # an item missing from the dict falls back to the TEA lifetime
    u.equipment_lifetime = {'A': 10}
    sys.simulate()
    tea = qs.TEA(sys, discount_rate=r, lifetime=20, simulate_system=False)
    assert tea.annualized_equipment_cost == pytest.approx(
        A_cost / _A(10, r) + B_cost / _A(20, r))

    # a uniform dict must match the equivalent integer lifetime
    u.equipment_lifetime = {'A': 10, 'B': 10}
    sys.simulate()
    tea_dict = qs.TEA(sys, discount_rate=r, lifetime=20, simulate_system=False)
    u.equipment_lifetime = 10
    sys.simulate()
    tea_int = qs.TEA(sys, discount_rate=r, lifetime=20, simulate_system=False)
    assert tea_dict.annualized_equipment_cost == pytest.approx(tea_int.annualized_equipment_cost)


def test_depreciation_no_effect_without_taxable_income():
    '''With no income tax there is nothing to shield, so depreciation cannot change NPV.'''
    import qsdsan as qs
    from qsdsan.utils import create_example_system
    sys = create_example_system(); sys.simulate()
    npvs = []
    for dep in ('SL', 'MACRS5', 'MACRS7'):
        tea = qs.TEA(sys, discount_rate=0.05, lifetime=20, income_tax=0.,
                     depreciation=dep, simulate_system=False)
        npvs.append(tea.NPV)
    assert npvs[0] == pytest.approx(npvs[1]) == pytest.approx(npvs[2])


def test_depreciation_changes_npv_with_tax_and_profit():
    '''With income tax and a profit, accelerated depreciation raises NPV (defers tax).'''
    import qsdsan as qs
    from qsdsan.utils import create_example_system
    sys = create_example_system(); sys.simulate()
    [s for s in sys.products if s.ID == 'alcohols'][0].price = 5.0   # make it profitable
    sl = qs.TEA(sys, discount_rate=0.05, lifetime=20, income_tax=0.21,
                depreciation='SL', simulate_system=False).NPV
    macrs = qs.TEA(sys, discount_rate=0.05, lifetime=20, income_tax=0.21,
                   depreciation='MACRS7', simulate_system=False).NPV
    assert macrs > sl                       # accelerated schedule defers tax -> higher NPV
    assert macrs != pytest.approx(sl)


def test_depreciation_schedule_must_fit_lifetime():
    '''A MACRS schedule longer than the lifetime is rejected; SL always fits.'''
    import qsdsan as qs
    from qsdsan.utils import create_example_system
    sys = create_example_system(); sys.simulate()
    # MACRS7 is an 8-year schedule -> needs lifetime >= 8
    with pytest.raises(RuntimeError):
        qs.TEA(sys, discount_rate=0.05, lifetime=7, depreciation='MACRS7',
               simulate_system=False).NPV
    assert qs.TEA(sys, discount_rate=0.05, lifetime=8, depreciation='MACRS7',
                  simulate_system=False).NPV < 0
    # straight line spans the whole lifetime, so it fits even when short
    assert qs.TEA(sys, discount_rate=0.05, lifetime=7, depreciation='SL',
                  simulate_system=False).NPV < 0


def test_CEPCI_accessors():
    '''qs.CEPCI is a settable global; CEPCI_by_year is a year->index table.'''
    import qsdsan as qs
    assert isinstance(qs.CEPCI_by_year, dict) and 2023 in qs.CEPCI_by_year
    original = qs.CEPCI
    try:
        qs.CEPCI = qs.CEPCI_by_year[2023]
        assert qs.CEPCI == qs.CEPCI_by_year[2023]
        tea, sys = _example_tea(CEPCI=qs.CEPCI_by_year[2022])
        assert tea.CEPCI == qs.CEPCI_by_year[2022]
    finally:
        qs.CEPCI = original


def test_simpletea_is_deprecated_alias():
    '''SimpleTEA still works but warns; it is a TEA subclass.'''
    import qsdsan as qs
    from qsdsan.utils import create_example_system
    sys = create_example_system(); sys.simulate()
    with pytest.warns(DeprecationWarning):
        tea = qs.SimpleTEA(sys, discount_rate=0.05, lifetime=10, simulate_system=False)
    assert isinstance(tea, qs.TEA)


def test_agile_system_npv_charges_equipment_replacement():
    '''Regression: ``qs.TEA`` on an ``AgileSystem`` reads the equipment lifetime
    from the per-unit ``UnitDesignAndCapital`` objects (attribute
    ``equipment_lifetime``), not a ``lifetime`` attribute. Computing ``NPV``
    previously raised ``AttributeError`` whenever a replacement had to be charged.
    '''
    import qsdsan as qs
    # equipment that lasts the whole project: no replacement charged
    ag_norepl = _example_agile_system('test_tea_agile_norepl', equipment_lifetime=0)
    npv_norepl = qs.TEA(ag_norepl, discount_rate=0.05, lifetime=20,
                        simulate_system=False).NPV
    # an 8-year equipment lifetime forces a mid-project replacement
    ag_repl = _example_agile_system('test_tea_agile_repl', equipment_lifetime=8)
    npv_repl = qs.TEA(ag_repl, discount_rate=0.05, lifetime=20,
                      simulate_system=False).NPV
    # both are pure-cost systems (negative NPV); the replacement makes it more negative
    assert npv_repl < npv_norepl < 0


if __name__ == '__main__':
    for name in __all__:
        globals()[name]()
        print(f'{name}: ok')
