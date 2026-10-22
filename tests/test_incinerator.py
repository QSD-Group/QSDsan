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

__all__ = (
    'test_incinerator_mass_balance',
    'test_incinerator_measured_as_conversion',
    'test_incinerator_dynamic_matches_static',
    'test_incinerator_insufficient_oxygen',
    'test_incinerator_energy_balance',
    'test_incinerator_organic_oxygen_demand_equals_cod',
    )


def _make_cmps(adjust_MW=False):
    '''
    Components of the ``Incinerator`` doctest, compiled with or without adjusting
    the MW of components with ``measured_as`` (e.g., S_NH4 as N, S_CH4 as COD)
    to their measured-as unit. Flows of such components are in the measured-as
    unit either way (see the note on ``measured_as`` in :class:`qsdsan.Component`).
    '''
    import qsdsan as qs
    cmps = qs.Components.load_default()
    CO2 = qs.Component.from_chemical('S_CO2', chemical='CO2', particle_size='Soluble',
                                     degradability='Undegradable', organic=False)
    cmps = qs.Components([cmps.S_F, cmps.S_NH4, cmps.X_OHO, cmps.H2O, cmps.S_CH4,
                          cmps.S_O2, cmps.S_N2, cmps.S_H2, cmps.X_Ig_ISS, CO2])
    cmps.default_compile(ignore_inaccurate_molar_weight=True,
                         adjust_MW_to_measured_as=adjust_MW)
    qs.set_thermo(cmps)
    return cmps


def _make_asm2d_cmps():
    '''ASM2d components (the MW is adjusted) with CO2 and ash, as in EXPOsan ``metro``.'''
    import qsdsan as qs
    from qsdsan import process_models as pc
    CO2 = qs.Component.from_chemical('S_CO2', chemical='CO2', particle_size='Dissolved gas',
                                     degradability='Undegradable', organic=False)
    ash = qs.Components.load_default().X_Ig_ISS.copy('ash')
    cmps = qs.Components([*pc.create_asm2d_cmps(), CO2, ash])
    cmps.compile(ignore_inaccurate_molar_weight=True)
    qs.set_thermo(cmps)
    return cmps


def _make_incinerator(sludge, air, fuel, dynamic=False, **kwargs):
    from qsdsan import unit_operations as su
    kwargs.setdefault('ash_component_ID', 'X_Ig_ISS')
    return su.Incinerator('', ins=(sludge, air, fuel), outs=('', ''), isdynamic=dynamic, **kwargs)


def _real_mass(cmps, stream):
    # Flows are in the measured-as unit, i_mass converts them to real mass
    return float((stream.mass * cmps.i_mass).sum())


def test_incinerator_mass_balance():
    '''Real mass is conserved for random influents, with and without MW adjustment.'''
    import numpy as np
    from numpy.testing import assert_allclose
    from qsdsan import WasteStream

    rng = np.random.default_rng(0)
    for adjust_MW in (False, True):
        cmps = _make_cmps(adjust_MW)
        for _ in range(5):
            sludge = WasteStream()
            sludge.set_flow(rng.random(len(cmps)) * 10, 'kg/hr', cmps.IDs)
            air = WasteStream(phase='g', S_O2=5000, S_N2=15000)
            fuel = WasteStream(phase='g', S_CH4=10)
            inc = _make_incinerator(sludge, air, fuel)
            inc.simulate()
            m_in = sum(_real_mass(cmps, i) for i in inc.ins)
            m_out = sum(_real_mass(cmps, i) for i in inc.outs)
            assert_allclose(m_out, m_in, rtol=1e-9)


def test_incinerator_measured_as_conversion():
    '''
    Acetic acid as itself (60.052 kg/hr) and as COD (64 kg/hr), the example in the
    docstring of ``Components.default_compile``, are both 1 kmol/hr of CH3COOH, so
    CH3COOH + 2 O2 -> 2 CO2 + 2 H2O gives the same products with and without
    adjusting the MW to the measured-as unit. (``WasteStream.mol`` of the COD component
    is 1.066 instead of 1 without adjusting, but the combustion must not depend on it.)
    '''
    import qsdsan as qs
    from numpy.testing import assert_allclose
    from qsdsan import WasteStream

    for adjust_MW in (False, True):
        gas = lambda ID, chem: qs.Component.from_chemical(
            ID, chemical=chem, particle_size='Dissolved gas', degradability='Undegradable', organic=False)
        cmps = qs.Components([
            qs.Component('Ac', search_ID='CH3COOH', particle_size='Soluble', degradability='Readily', organic=True),
            qs.Component('Ac_as_COD', search_ID='CH3COOH', measured_as='COD', particle_size='Soluble',
                         degradability='Readily', organic=True),
            qs.Component.from_chemical('H2O', chemical='H2O', particle_size='Soluble',
                                       degradability='Undegradable', organic=False),
            gas('S_O2', 'O2'), gas('S_N2', 'N2'), gas('S_CO2', 'CO2'),
            qs.Components.load_default().X_Ig_ISS.copy('ash'),
            ])
        cmps.default_compile(ignore_inaccurate_molar_weight=True, adjust_MW_to_measured_as=adjust_MW)
        qs.set_thermo(cmps)
        sludge = WasteStream(Ac=60.052, Ac_as_COD=64)
        air = WasteStream(phase='g', S_O2=1000)
        inc = _make_incinerator(sludge, air, WasteStream(phase='g'), ash_component_ID='ash')
        inc.simulate()
        flue_gas, ash = inc.outs
        # 2 kmol/hr of acetic acid
        assert_allclose(flue_gas.imass['S_CO2'], 4 * 44.0095, rtol=1e-4)
        assert_allclose(flue_gas.imass['H2O'], 4 * 18.01528, rtol=1e-4)
        assert_allclose(flue_gas.imass['S_O2'], 1000 - 4 * 31.9988, rtol=1e-4)
        assert ash.F_mass < 1e-9


def test_incinerator_dynamic_matches_static():
    '''
    The algebraic equations in a dynamic system give the same flows as the static
    run (the dynamic equations used to fail for more than one influent, and the
    states of the flue gas and ash were integer arrays that truncated the flows).
    '''
    import numpy as np
    from numpy.testing import assert_allclose
    from qsdsan import System, WasteStream

    _make_cmps()
    def make(dynamic):
        return _make_incinerator(
            WasteStream(S_F=10, S_NH4=20, X_OHO=15, H2O=1000),
            WasteStream(phase='g', S_O2=1500, S_N2=4900),
            WasteStream(phase='g', S_CH4=1000),
            dynamic=dynamic)
    static, dynamic = make(False), make(True)
    System('', path=(static,)).simulate()
    System('', path=(dynamic,)).simulate(t_span=(0, 1), method='BDF')
    IDs = ('S_CO2', 'S_N2', 'H2O', 'S_O2')
    assert_allclose([dynamic.outs[0].imass[i] for i in IDs],
                    [static.outs[0].imass[i] for i in IDs], rtol=1e-9)
    assert_allclose(dynamic.outs[1].imass['X_Ig_ISS'], static.outs[1].imass['X_Ig_ISS'], rtol=1e-9)
    assert dynamic.outs[0].state.dtype == float
    assert dynamic.outs[1].state.dtype == float


def test_incinerator_insufficient_oxygen():
    '''A warning is raised and no negative oxygen is in the flue gas if the air is not enough.'''
    import pytest
    from qsdsan import WasteStream

    _make_cmps()
    inc = _make_incinerator(WasteStream(S_F=10, H2O=100),
                            WasteStream(phase='g', S_O2=10, S_N2=780),
                            WasteStream(phase='g', S_CH4=1000))
    with pytest.warns(UserWarning, match='enough oxygen'):
        inc.simulate()
    assert inc.outs[0].imass['S_O2'] == 0


def test_incinerator_energy_balance():
    '''Heats are from the heating values of the components, or the calorific values if given.'''
    from numpy.testing import assert_allclose
    from qsdsan import WasteStream

    cmps = _make_cmps()
    methane = 100 # kg/hr of real CH4, the flow is in COD
    air_T = 400 # K
    inc = _make_incinerator(
        WasteStream(S_F=10, H2O=1000),
        WasteStream(phase='g', T=air_T, S_O2=1000, S_N2=3000),
        WasteStream(phase='g', S_CH4=methane/cmps.S_CH4.i_mass),
        process_efficiency=0.8)
    inc.simulate()
    LHV, HHV, MW = 802567., 890590., 16.04246 # CH4 [kJ/kmol], [kJ/kmol], [kg/kmol]
    assert_allclose(inc.heat_fuel, methane/MW*LHV, rtol=1e-4)
    # Composites (S_F) are from Thornton's rule, 13.1 MJ/kg O2 [3, 4], and 1 kg COD needs 1 kg O2
    assert_allclose(inc.heat_sludge, 10*13.1e3, rtol=1e-3)
    assert_allclose(inc.heat_combustion, inc.heat_sludge+inc.heat_fuel, rtol=1e-12)
    assert_allclose(inc.heat_air, 4000*1.005*(air_T-298.15), rtol=1e-12)
    heat_in = inc.heat_combustion + inc.heat_air
    assert_allclose(inc.heat_flue_gas, 0.8*heat_in, rtol=1e-12)
    assert_allclose(inc.heat_loss, 0.2*heat_in, rtol=1e-12)
    # 1000 kg/hr of water heated to boiling and evaporated, about 46 kJ/mol
    assert_allclose(inc.heat_evaporation, 1000/18.015*46.3e3, rtol=0.05)

    inc.heating_value_basis = 'HHV'
    inc._design()
    assert_allclose(inc.heat_fuel, methane/MW*HHV, rtol=1e-4)
    inc.heating_value_basis = 'LHV'
    inc.calorific_value_sludge, inc.calorific_value_fuel = 12000, 50000
    inc._design()
    assert_allclose(inc.heat_sludge, 10*cmps.S_F.i_mass*12000, rtol=1e-12) # dry real mass
    assert_allclose(inc.heat_fuel, methane*50000, rtol=1e-12)


def test_incinerator_organic_oxygen_demand_equals_cod():
    '''
    The oxygen consumed by each organic component is its COD, including ions without
    the counter ion in their formula (acetate), whose BioSTEAM combustion reaction
    would need 12.5% less oxygen.
    '''
    import qsdsan as qs
    from numpy.testing import assert_allclose
    from qsdsan import WasteStream

    cmps = _make_asm2d_cmps()
    tested = 0
    for cmp in cmps:
        if not cmp.organic or cmp.i_COD <= 0 or cmp.ID in ('S_O2', 'S_N2', 'S_CO2', 'H2O', 'ash'):
            continue
        sludge = WasteStream()
        sludge.set_flow([1.], 'kg/hr', (cmp.ID,)) # 1 kg/hr as COD
        air = WasteStream(phase='g', S_O2=100)
        inc = _make_incinerator(sludge, air, WasteStream(phase='g'), ash_component_ID='ash')
        inc.simulate()
        assert_allclose(100 - inc.outs[0].imass['S_O2'], cmp.i_COD, rtol=1e-6, err_msg=cmp.ID)
        tested += 1
    assert tested >= 8 # at least the ASM2d organics, including S_A


if __name__ == '__main__':
    test_incinerator_mass_balance()
    test_incinerator_measured_as_conversion()
    test_incinerator_dynamic_matches_static()
    test_incinerator_insufficient_oxygen()
    test_incinerator_energy_balance()
    test_incinerator_organic_oxygen_demand_equals_cod()
