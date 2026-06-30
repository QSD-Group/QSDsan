#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.

Regression tests for equipment/construction lifetimes reaching the attribute
that TEA and LCA actually read (``equipment_lifetime``). Historically these
lifetimes were written to ``_default_equipment_lifetime`` and then dropped,
because ``SanUnit`` reset ``equipment_lifetime`` to ``{}`` at construction and
nothing synced the two, so equipment replacement costs/impacts were ignored.
'''

__all__ = (
    'test_class_default_equipment_lifetime_applied',
    'test_runtime_equipment_design_lifetime_applied',
    'test_construction_lifetime_applied',
    'test_explicit_lifetime_overrides_class_default',
    'test_equipment_lifetime_is_instance_local',
    'test_real_unit_polishing_filter_lifetime_applied',
    )

import pytest


def _one_component_thermo(flowsheet_ID):
    import qsdsan as qs
    from qsdsan import Component, Components
    cmps = Components([Component('Water', search_ID='Water',
                                 particle_size='Soluble',
                                 degradability='Undegradable', organic=False)])
    cmps.compile()
    qs.set_thermo(cmps)
    qs.main_flowsheet.set_flowsheet(qs.Flowsheet(flowsheet_ID))
    return cmps


def _equip_unit_cls():
    '''A SanUnit subclass that attaches a dynamic Equipment item with a lifetime.'''
    import qsdsan as qs
    from qsdsan import SanUnit, Equipment

    class _Pump(Equipment):
        def _design(self): return {}
        def _cost(self): return 1e6

    class EquipUnit(SanUnit):
        _N_ins = 1; _N_outs = 1; _F_BM_default = {'Body': 1.}
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.equipment = [_Pump(linked_unit=self, ID='pump', lifetime=12, F_BM=1.)]
        def _run(self):
            self.outs[0].copy_like(self.ins[0])
        def _cost(self):
            self.baseline_purchase_costs['Body'] = 1e6
            self.add_equipment_design(); self.add_equipment_cost()

    return EquipUnit


def test_class_default_equipment_lifetime_applied():
    '''A class-level ``_default_equipment_lifetime`` must reach ``equipment_lifetime``.'''
    import qsdsan as qs
    from qsdsan import SanUnit, WasteStream
    _one_component_thermo('test_eqlife_classdefault')

    class Plant(SanUnit):
        _N_ins = 1; _N_outs = 1; _F_BM_default = {'Body': 1.}
        _default_equipment_lifetime = {'Body': 12}
        def _run(self): self.outs[0].copy_like(self.ins[0])
        def _cost(self): self.baseline_purchase_costs['Body'] = 1e6

    u = Plant('u', ins=WasteStream('ws_cd', Water=1000), outs='e')
    assert u.equipment_lifetime == {'Body': 12}
    assert u.lifetime == {'Body': 12}     # the property mirrors equipment_lifetime


def test_runtime_equipment_design_lifetime_applied():
    '''An ``Equipment`` lifetime added during ``_cost`` must reach ``equipment_lifetime``.'''
    import qsdsan as qs
    from qsdsan import WasteStream
    _one_component_thermo('test_eqlife_runtime')
    EquipUnit = _equip_unit_cls()
    u = EquipUnit('u', ins=WasteStream('ws_rt', Water=1000), outs='e')
    qs.System('s', path=(u,)).simulate()
    # the equipment ID is prefixed with the unit ID in the design/cost registers
    assert u.equipment_lifetime.get(u.equipment[0].ID) == 12


def test_construction_lifetime_applied():
    '''A ``Construction`` lifetime added via ``add_construction`` reaches ``equipment_lifetime``.'''
    import qsdsan as qs
    from qsdsan import SanUnit, WasteStream, ImpactItem, Construction
    _one_component_thermo('test_eqlife_constr')
    item = ImpactItem('Steel_cl', functional_unit='kg')

    class Plant(SanUnit):
        _N_ins = 1; _N_outs = 1; _F_BM_default = {'Body': 1.}
        def _run(self): self.outs[0].copy_like(self.ins[0])
        def _cost(self): self.baseline_purchase_costs['Body'] = 1e6

    u = Plant('u', ins=WasteStream('ws_cl', Water=1000), outs='e')
    u.construction = [Construction(linked_unit=u, item=item, quantity=10,
                                   quantity_unit='kg', lifetime=20)]
    u.add_construction(add_unit=False, add_design=False, add_cost=False, add_lifetime=True)
    assert u.equipment_lifetime.get('Steel_cl') == 20


def test_explicit_lifetime_overrides_class_default():
    '''An explicit ``lifetime`` argument takes precedence over the class default.'''
    import qsdsan as qs
    from qsdsan import SanUnit, WasteStream
    _one_component_thermo('test_eqlife_override')

    class Plant(SanUnit):
        _N_ins = 1; _N_outs = 1; _F_BM_default = {'Body': 1.}
        _default_equipment_lifetime = {'Body': 12}
        def _run(self): self.outs[0].copy_like(self.ins[0])
        def _cost(self): self.baseline_purchase_costs['Body'] = 1e6

    u = Plant('u', ins=WasteStream('ws_ov', Water=1000), outs='e', lifetime=30)
    assert u.equipment_lifetime == 30


def test_equipment_lifetime_is_instance_local():
    '''Per-instance lifetimes must not leak through the shared class-level dict.'''
    import qsdsan as qs
    from qsdsan import WasteStream
    EquipUnit = _equip_unit_cls()

    _one_component_thermo('test_eqlife_local_a')
    a = EquipUnit('a', ins=WasteStream('wa', Water=1000), outs='ea')
    qs.System('sa', path=(a,)).simulate()
    a.equipment_lifetime[a.equipment[0].ID] = 99   # mutate one instance

    _one_component_thermo('test_eqlife_local_b')
    b = EquipUnit('b', ins=WasteStream('wb', Water=1000), outs='eb')
    qs.System('sb', path=(b,)).simulate()
    assert b.equipment_lifetime.get(b.equipment[0].ID) == 12   # not contaminated by `a`
    # the class-level default dict was never mutated by either instance
    assert type(b)._default_equipment_lifetime == {}


def test_real_unit_polishing_filter_lifetime_applied():
    '''A shipped unit (``PolishingFilter``) must expose its default equipment lifetimes.'''
    import qsdsan as qs
    qs.main_flowsheet.clear()
    from qsdsan import Component, Components
    from qsdsan.utils import create_example_components
    ref = create_example_components()

    def _cmp(ID, ps, deg, org):
        c = qs.Component(ID, particle_size=ps, degradability=deg, organic=org)
        c.copy_models_from(ref.H2O, names=('V',)); c.default(); return c

    cmps = Components([_cmp('AceticAcid', 'Soluble', 'Readily', True),
                       _cmp('CH4', 'Dissolved gas', 'Undegradable', False),
                       _cmp('CO2', 'Dissolved gas', 'Undegradable', False),
                       _cmp('H2O', 'Soluble', 'Undegradable', False),
                       _cmp('O2', 'Dissolved gas', 'Undegradable', False),
                       _cmp('N2', 'Dissolved gas', 'Undegradable', False),
                       qs.Component('Biomass', phase='s', formula='CH1.8O0.5N0.2',
                                    particle_size='Particulate', degradability='Slowly',
                                    organic=True)])
    cmps.Biomass.copy_models_from(ref.H2O, names=('V',)); cmps.Biomass.default()
    cmps.compile(); qs.set_thermo(cmps)

    split = {'AceticAcid': 0.125, 'CH4': 0.15, 'CO2': 0.15, 'H2O': 0.125,
             'O2': 1.0, 'N2': 1.0, 'Biomass': 0.0}
    pf = qs.sanunits.PolishingFilter(
        'PF', ins=(qs.WasteStream('pf_in', AceticAcid=1, H2O=100, units='kg/hr'),
                   qs.WasteStream('pf_rec'), qs.WasteStream('pf_air', phase='g')),
        outs=('bg', 'eff', 'waste', 'air_out'),
        filter_type='aerobic', biomass_ID='Biomass', split=split, solids=('Biomass',))
    # the default lifetimes declared on the unit must be live, not dropped
    assert pf.equipment_lifetime.get('Pumps') == 15


if __name__ == '__main__':
    for name in __all__:
        globals()[name]()
        print(f'{name}: ok')
