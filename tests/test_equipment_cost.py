#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.

Regression tests for ``Equipment.installed_cost``, which is
``baseline purchase cost * (F_BM + F_D * F_P * F_M - 1)``, the same as what
BioSTEAM calculates for the cost items of a unit.
'''

__all__ = (
    'test_equipment_installed_cost_float_cost',
    'test_equipment_installed_cost_dict_cost',
    'test_equipment_installed_cost_matches_unit',
    'test_equipment_baseline_purchase_cost',
    'test_unit_applies_dict_cost_equipment_factors',
    'test_unit_applies_dict_cost_equipment_lifetime',
    'test_blower_cost_items_have_bare_module_factors',
    )

import pytest


def _equip_cls():
    from qsdsan import Equipment

    class _Equip(Equipment):
        cost = 1000.
        def _design(self): return {}
        def _cost(self): return self.cost

    return _Equip


def test_equipment_installed_cost_float_cost():
    equip = _equip_cls()(ID='eq_float', F_BM=2.3, F_D=1.1, F_P=1.2, F_M=2.)
    equip._baseline_purchase_costs = equip._cost()
    assert equip.installed_cost == pytest.approx(1000 * (2.3 + 1.1*1.2*2. - 1))


def test_equipment_installed_cost_dict_cost():
    from qsdsan import Equipment

    class _Parts(Equipment):
        def _design(self): return {}
        def _cost(self): return {'Vessel': 1000., 'Pump': 500.}

    equip = _Parts(ID='eq_parts', F_BM={'Vessel': 4., 'Pump': 3.}, F_M={'Vessel': 2.})
    equip._baseline_purchase_costs = equip._cost()
    # Vessel: 1000*(4+2-1); Pump: 500*(3+1-1) since its F_M is not given
    assert equip.installed_cost == pytest.approx(1000*5 + 500*3)


def test_equipment_installed_cost_matches_unit():
    '''Installed cost of an equipment should equal what the unit reports for it.'''
    import qsdsan as qs
    from qsdsan import Component, Components, SanUnit
    cmps = Components([Component('Water', search_ID='Water', particle_size='Soluble',
                                 degradability='Undegradable', organic=False)])
    cmps.compile()
    qs.set_thermo(cmps)
    qs.main_flowsheet.set_flowsheet(qs.Flowsheet('test_equipment_cost'))
    Equip = _equip_cls()

    class U(SanUnit):
        _N_ins = 1; _N_outs = 1
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.equipment = [Equip(linked_unit=self, ID='eq', F_BM=2.3, F_D=1.1, F_P=1.2, F_M=2.)]
        def _run(self): self.outs[0].copy_like(self.ins[0])
        def _design(self): self.add_equipment_design()
        def _cost(self): self.add_equipment_cost()

    u = U('U1', ins=qs.WasteStream('in', Water=1000, units='kg/hr'), outs=qs.WasteStream('out'))
    u.simulate()
    assert u.equipment[0].installed_cost == pytest.approx(u.installed_cost)


def test_equipment_baseline_purchase_cost():
    '''``_cost`` can return a float or a dict; the total baseline cost must work for both.'''
    from qsdsan import Equipment
    Flt = _equip_cls()

    class _Parts(Equipment):
        def _design(self): return {}
        def _cost(self): return {'Vessel': 1000., 'Pump': 500.}

    flt = Flt(ID='eq_f')
    assert flt.baseline_purchase_cost == 0. # not calculated yet
    flt._baseline_purchase_costs = flt._cost()
    assert flt.baseline_purchase_cost == pytest.approx(1000.)

    parts = _Parts(ID='eq_p')
    parts._baseline_purchase_costs = parts._cost()
    assert parts.baseline_purchase_cost == pytest.approx(1500.)


def _unit_with_parts_equipment(ID, **equip_kwargs):
    '''A unit with its own cost item and an equipment whose ``_cost`` returns a dict of parts.'''
    import qsdsan as qs
    from qsdsan import Component, Components, SanUnit, Equipment
    cmps = Components([Component('Water', search_ID='Water', particle_size='Soluble',
                                 degradability='Undegradable', organic=False)])
    cmps.compile()
    qs.set_thermo(cmps)
    qs.main_flowsheet.set_flowsheet(qs.Flowsheet(f'test_equipment_cost_{ID}'))

    class _Parts(Equipment):
        def _design(self): return {}
        def _cost(self): return {'Vessel': 1000., 'Pump': 500.}

    class U(SanUnit):
        _N_ins = 1; _N_outs = 1
        _F_BM_default = {'Body': 3.}
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.equipment = [_Parts(linked_unit=self, ID='eq', **equip_kwargs)]
        def _run(self): self.outs[0].copy_like(self.ins[0])
        def _design(self): self.add_equipment_design()
        def _cost(self):
            self.baseline_purchase_costs['Body'] = 2000.
            self.add_equipment_cost()

    return U(ID, ins=qs.WasteStream(f'{ID}_in', Water=1000, units='kg/hr'),
             outs=qs.WasteStream(f'{ID}_out'))


@pytest.mark.parametrize('factors, expected', [
    # Dict factors are looked up per part; Pump has no F_M so it is 1.
    (dict(F_BM={'Vessel': 4., 'Pump': 3.}, F_M={'Vessel': 2.}), 1000*(4+2-1) + 500*(3+1-1)),
    # Single-value factors apply to all parts.
    (dict(F_BM=2.3, F_D=1.1, F_P=1.2, F_M=2.), (1000+500)*(2.3+1.1*1.2*2.-1)),
    ])
def test_unit_applies_dict_cost_equipment_factors(factors, expected):
    '''Factors of an equipment with dict cost must reach the unit's cost items exactly once.'''
    import warnings
    u = _unit_with_parts_equipment(f'U{abs(hash(str(factors)))%10**6}', **factors)
    with warnings.catch_warnings():
        # A missing bare-module factor would warn and silently fall back to 1
        warnings.simplefilter('error', RuntimeWarning)
        u.simulate()
    body = 2000*3.
    assert u.installed_cost == pytest.approx(body + expected)
    assert u.equipment[0].installed_cost == pytest.approx(expected)


@pytest.mark.parametrize('lifetime, expected', [
    ({'Vessel': 20, 'Pump': 5}, {'Vessel': 20, 'Pump': 5}),
    (12, {'Vessel': 12, 'Pump': 12}),
    ])
def test_unit_applies_dict_cost_equipment_lifetime(lifetime, expected):
    '''Lifetime of an equipment with dict cost must be keyed by the unit's cost item names.'''
    u = _unit_with_parts_equipment(f'L{abs(hash(str(lifetime)))%10**6}', lifetime=lifetime)
    u.simulate()
    items = u.baseline_purchase_costs
    for part, years in expected.items():
        name = f'_Parts {u.equipment[0].ID} - {part}'
        assert name in items
        assert u.equipment_lifetime[name] == years


def test_blower_cost_items_have_bare_module_factors():
    '''Every cost item of a ``Blower`` must have a matching entry in its ``F_BM``.'''
    import qsdsan as qs
    from qsdsan import Component, Components, SanUnit
    from qsdsan.equipments import Blower
    cmps = Components([Component('Water', search_ID='Water', particle_size='Soluble',
                                 degradability='Undegradable', organic=False)])
    cmps.compile()
    qs.set_thermo(cmps)
    qs.main_flowsheet.set_flowsheet(qs.Flowsheet('test_blower_F_BM'))

    class U(SanUnit):
        _N_ins = 1; _N_outs = 1
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.equipment = [Blower('blower', linked_unit=self, gas_demand_per_reactor=1000.)]
        def _run(self): self.outs[0].copy_like(self.ins[0])
        def _design(self): self.add_equipment_design()
        def _cost(self): self.add_equipment_cost()

    import warnings
    u = U('U_blower', ins=qs.WasteStream('b_in', Water=1000, units='kg/hr'), outs=qs.WasteStream('b_out'))
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning) # no cost item without a bare-module factor
        u.simulate()
    blower = u.equipment[0]
    assert set(blower.F_BM) == set(blower.baseline_purchase_costs)
    assert u.installed_cost == pytest.approx(blower.installed_cost)
