# CODEX adds the followings 164 lines
# -*- coding: utf-8 -*-
'''Tests for the dynamic sand-filtration unit.'''

import numpy as np
import pytest
from numpy.testing import assert_allclose

import qsdsan as qs
from qsdsan import process_models as pc
from qsdsan import unit_operations as su


def _set_thermo():
    cmps = pc.create_asm1_cmps()
    qs.set_thermo(cmps)
    return cmps


def _wastewater(ID, flow=1000., soluble=100., particulate=50.):
    stream = qs.WasteStream(ID)
    stream.set_flow_by_concentration(
        flow, {'S_S': soluble, 'X_S': particulate},
        units=('m3/d', 'mg/L'))
    return stream


def test_constant_splitting_defaults_and_static_mass_balance():
    _set_thermo()
    influent = _wastewater('sand_constant_in')
    unit = su.SandFiltration(
        'sand_constant', ins=influent,
        outs=('sand_constant_eff', 'sand_constant_bw'),
        model='constant_spliting', isdynamic=False,
        energy_consumption=0.04)
    unit.simulate()

    effluent, backwash = unit.outs
    assert unit.model == 'constant_splitting'
    assert len(unit.ins) == 1
    assert len(unit.outs) == 2
    assert unit.backwash_flow_fraction == 0.02
    assert unit.backwash_solids_mass_fraction == 0.8
    assert unit.energy_consumption == 0.04
    assert_allclose(backwash.imass['S_S']/influent.imass['S_S'], 0.02)
    assert_allclose(backwash.imass['X_S']/influent.imass['X_S'], 0.8)
    assert_allclose(effluent.mass+backwash.mass, influent.mass, rtol=1e-12)


def test_constant_splitting_dynamic_algebraic_equations():
    cmps = _set_thermo()
    influent = _wastewater('sand_constant_dynamic_in')
    unit = su.SandFiltration(
        'sand_constant_dynamic', ins=influent,
        outs=('sand_constant_dynamic_eff', 'sand_constant_dynamic_bw'))
    system = qs.System('sand_constant_dynamic_sys', path=(unit,))
    system.set_dynamic_tracker(*unit.outs)
    system.simulate(t_span=(0., 0.1), t_eval=(0., 0.1))

    effluent, backwash = unit.outs
    soluble = cmps.index('S_S')
    particulate = cmps.index('X_S')
    assert system.scope.sol.success
    assert_allclose([effluent.state[-1], backwash.state[-1]], [980., 20.])
    assert_allclose(
        [effluent.state[soluble], backwash.state[soluble]], [100., 100.])
    assert_allclose(
        [effluent.state[particulate], backwash.state[particulate]],
        [50.*0.2/0.98, 50.*0.8/0.02])


def test_empirical_fouling_defaults_curve_and_static_mass_balance():
    _set_thermo()
    influent = _wastewater('sand_empirical_in')
    backwash_in = _wastewater(
        'sand_empirical_bw_in', flow=1000./24.,
        soluble=10., particulate=0.)
    unit = su.SandFiltration(
        'sand_empirical', ins=(influent, backwash_in),
        outs=('sand_empirical_eff', 'sand_empirical_bw_out',
              'sand_empirical_recycle'),
        model='massbalance', isdynamic=False)
    unit.simulate()

    effluent, backwash_out, recycle = unit.outs
    assert unit.model == 'empirical_fouling'
    assert len(unit.ins) == 2
    assert len(unit.outs) == 3
    assert unit.filter_cycle_time == 24.
    assert unit.filter_backwash_time == 1.
    assert unit.timeshift_for_cycles == 0.
    assert unit.fouling_coefficient == 0.03
    assert unit.best_removal_efficiency == 0.75
    assert unit.minimum_removal_efficiency == 0.30
    assert unit.average_backwash_fraction == 1./24.

    e0 = unit.get_removal_efficiency(0.)
    e23 = 1.-0.25*np.exp(0.03*23.)
    assert_allclose(e0, 0.75)
    assert_allclose(unit.get_removal_efficiency(23./24.), e23)
    assert_allclose(
        unit.get_removal_efficiency(23.5/24.), (e23+0.75)/2.)
    assert_allclose(unit.get_removal_efficiency(1.), 0.75)

    assert_allclose(
        effluent.mass+backwash_out.mass+recycle.mass,
        influent.mass+backwash_in.mass, rtol=1e-12, atol=1e-12)
    filtered_XS = effluent.imass['X_S'] + recycle.imass['X_S']
    captured_XS = backwash_out.imass['X_S']-backwash_in.imass['X_S']
    assert_allclose(filtered_XS/influent.imass['X_S'], 0.25)
    assert_allclose(captured_XS/influent.imass['X_S'], 0.75)


def test_empirical_fouling_dynamic_inventory_and_mass_balance():
    _set_thermo()
    influent = _wastewater('sand_empirical_dynamic_in')
    backwash_in = _wastewater(
        'sand_empirical_dynamic_bw_in', flow=1000./24.,
        soluble=10., particulate=0.)
    unit = su.SandFiltration(
        'sand_empirical_dynamic', ins=(influent, backwash_in),
        outs=('sand_empirical_dynamic_eff',
              'sand_empirical_dynamic_bw_out',
              'sand_empirical_dynamic_recycle'),
        model='empirical_fouling')
    system = qs.System('sand_empirical_dynamic_sys', path=(unit,))
    system.set_dynamic_tracker(unit, *unit.outs)
    system.simulate(
        method='BDF', t_span=(0., 2.), t_eval=np.linspace(0., 2., 49),
        rtol=1e-6, atol=1e-8)

    assert system.scope.sol.success
    assert unit.captured_solids_inventory.min() > -1e-8
    assert unit.captured_solids_inventory.max() > 0.
    assert_allclose(
        [stream.state[-1] for stream in unit.outs],
        [1000.*23./24., 1000./24., 1000./24.])

    mass_in = sum(stream.state[-1]*stream.state[:-1]
                  for stream in unit.ins)
    mass_out = sum(stream.state[-1]*stream.state[:-1]
                   for stream in unit.outs)
    assert_allclose(mass_in-mass_out, unit._dstate, rtol=1e-7, atol=1e-7)


# CODEX adds the followings 27 lines
def test_empirical_fouling_closed_backwash_recycle():
    _set_thermo()
    influent = _wastewater('sand_loop_in')
    recycle = qs.WasteStream('sand_loop_recycle')
    unit = su.SandFiltration(
        'sand_loop', ins=(influent, recycle),
        outs=('sand_loop_eff', 'sand_loop_bw_out', recycle),
        model='empirical_fouling')
    system = qs.System(
        'sand_loop_sys', path=(unit,), recycle=recycle)
    system.set_dynamic_tracker(unit, *unit.outs)
    system.simulate(
        method='BDF', t_span=(0., 0.2),
        t_eval=np.linspace(0., 0.2, 5), rtol=1e-6, atol=1e-8)

    assert system.scope.sol.success
    assert_allclose(
        recycle.state[-1],
        influent.state[-1]*unit.average_backwash_fraction)
    assert_allclose(
        sum(stream.state[-1] for stream in unit.ins),
        sum(stream.state[-1] for stream in unit.outs))


def test_sand_filtration_rejects_invalid_options():
    _set_thermo()
    with pytest.raises(ValueError, match='model must be'):
        su.SandFiltration(model='unknown')
    with pytest.raises(ValueError, match='less than filter_cycle_time'):
        su.SandFiltration(
            model='empirical_fouling', filter_cycle_time=1.,
            filter_backwash_time=1.)
    with pytest.raises(ValueError, match='nonzero backwash solids'):
        su.SandFiltration(
            backwash_flow_fraction=0.,
            backwash_solids_mass_fraction=0.8)
    # CODEX adds the followings 3 lines
    with pytest.raises(ValueError, match='cannot exceed best'):
        su.SandFiltration(
            best_removal_efficiency=0.5, minimum_removal_efficiency=0.6)
