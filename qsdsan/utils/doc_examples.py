#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from .. import (
    Component, Components, SanStream, System, SimpleTEA, Model,
    set_thermo, sanunits as su
    )
from chaospy import distributions as shape

__all__ = ('load_example_cmps', 'load_example_sys', 'load_example_model',)


# %%

# =============================================================================
# Example components
# =============================================================================

def load_example_cmps():
    '''
    Load pre-constructed components for documentation purpose.

    Returns
    -------
    A :class:`qsdsan.CompiledComponents` object with components including
    H2O (alias Water), CO2, N2O, NaCl, H2SO4, CH4 (alias Methane), Methanol, and Ethanol.

    Examples
    --------
    >>> from qsdsan.utils import load_example_cmps
    >>> cmps = load_example_cmps()
    >>> cmps.show()
    CompiledComponents([H2O, CO2, N2O, NaCl, H2SO4, CH4, Methanol, Ethanol])
    '''

    H2O = Component('H2O', search_ID='H2O', particle_size='Soluble',
                     degradability='Undegradable', organic=False)

    CO2 = Component('CO2', search_ID='CO2', phase='g',
                    particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)

    N2O = Component('N2O', search_ID='N2O', phase='g',
                    particle_size='Dissolved gas',
                    degradability='Undegradable', organic=False)

    NaCl = Component('NaCl', search_ID='NaCl', phase='s', particle_size='Soluble',
                     degradability='Undegradable', organic=False)

    H2SO4 = Component('H2SO4', search_ID='H2SO4', phase='s', particle_size='Soluble',
                      degradability='Undegradable', organic=False)

    CH4 = Component('CH4', search_ID='CH4', phase='g', particle_size='Dissolved gas',
                     degradability='Readily', organic=True)

    Methanol = Component('Methanol', search_ID='Methanol', phase='l',
                         particle_size='Soluble',
                         degradability='Readily', organic=True)

    Ethanol = Component('Ethanol', search_ID='Ethanol', phase='l',
                         particle_size='Soluble',
                         degradability='Readily', organic=True)

    cmps = Components((H2O, CO2, N2O, NaCl, H2SO4, CH4, Methanol, Ethanol))
    for cmp in cmps:
        cmp.default()

    cmps.compile()
    cmps.set_synonym('H2O', 'Water')
    cmps.set_synonym('CH4', 'Methane')

    return cmps

# %%

# =============================================================================
# Example systems
# =============================================================================

def load_example_sys(cmps=None):
    '''
    Load a pre-constructed system for documentation purpose.

    Returns
    -------
    A :class:`qsdsan.System` object including a mix tank, a pump, a heat exchanger,
    and various mixer/splitters.

    Parameters
    ----------
    cmps : obj
        If given, will call :func:`qsdsan.set_thermo(cmps)` to set components
        used in system simulation.
        The provided `cmps` must have "Water", "NaCl", "Methanol", and "Ethanol"
        components, which are used in the construction of the system.

    Examples
    --------
    >>> from qsdsan.utils import load_example_cmps, load_example_sys
    >>> cmps = load_example_cmps()
    >>> sys = load_example_sys(cmps)
    >>> sys.path
    (<MixTank: M1>,
     <Pump: P1>,
     <HXutility: H1>,
     <ComponentSplitter: S1>,
     <Mixer: M2>,
     <Splitter: S2>)
    >>> sys.diagram() # doctest: +SKIP
    '''

    if cmps:
        set_thermo(cmps)

    salt_water = SanStream('salt_water', Water=2000, NaCl=50, units='kg/hr')
    methanol = SanStream('methanol', Methanol=20, units='kg/hr')
    ethanol = SanStream('ethanol', Ethanol=10, units='kg/hr')

    M1 = su.MixTank('M1', ins=(salt_water, 'recycled_brine', methanol, ethanol))
    P1 = su.Pump('P1', ins=M1-0)
    H1 = su.HXutility('H1', ins=P1-0, T=350)
    S1 = su.ComponentSplitter('S1', ins=H1-0, split_keys=('Methanol', 'Ethanol'))
    M2 = su.Mixer('M2', ins=(S1-0, S1-1), outs='alcohols')
    S2 = su.Splitter('S2', ins=S1-2, outs=(1-M1, 'waste_brine'), split=0.2)
    sys = System('sys', path=(M1, P1, H1, S1, M2, S2))

    return sys


# %%

# =============================================================================
# Example system model
# =============================================================================

def load_example_model(evaluate=False, N=100, rule='L', seed=554, **sample_kwargs):
    '''
    Load a pre-constructed system model for documentation purpose.

    Parameters
    ----------
    evaluate : bool
        Whether to evaluate the model (i.e., simulate the system and get metrics).
    N : int
        Sample size, will be ignored if `evaluate` is set to False.
    rule : str
        Sampling rule, will be ignored if `evaluate` is set to False.
    seed : int
        Random seed for sample consistency, will be ignored if `evaluate` is set to False.
    sample_kwargs : dict
        Additional keyword arguments that will be passed to :func:`model.sample`

    Returns
    -------
    cmps : obj
        If given, will call :func:`qsdsan.set_thermo(cmps)` to set components
        used in system simulation.
        The provided `cmps` must have "Water", "NaCl", "Methanol", and "Ethanol"
        components, which are used in the construction of the system.

    Examples
    --------
    >>> from qsdsan.utils import load_example_model
    >>> model = load_example_model(N=100, rule='L', seed=554, evaluate=False)
    >>> model.system.path
    (<MixTank: M1>,
     <Pump: P1>,
     <HXutility: H1>,
     <ComponentSplitter: S1>,
     <Mixer: M2>,
     <Splitter: S2>)
    >>> model.parameters
    (<Parameter: [Stream-salt water] Salt flow rate (kg/hr)>,
     <Parameter: [HXutility-H1] Heat exchanger temperature (K)>,
     <Parameter: [Stream-salt water] Salt solution price (USD/kg)>,
     <Parameter: [Mix tank-M1] Mix tank retention time (hr)>,
     <Parameter: [Mix tank-M1] Mix tank mixer power usage (kW/m3)>,
     <Parameter: [Pump-P1] Pump design head (kPa)>)
    >>> model.metrics
    (<Metric: [System] Total heating duty (kJ/yr)>,
     <Metric: [System] Total electricity consumption (kWh/yr)>,
     <Metric: [Simple TEA] Total capital expenditure (USD)>,
     <Metric: [Simple TEA] Net present value (USD)>)
    '''

    cmps = load_example_cmps()
    sys = load_example_sys(cmps)
    M1, P1, H1, S1, M2, S2 = sys.path
    sys.simulate()
    tea = SimpleTEA(sys)
    model = Model(sys)

    # Add parameters
    param = model.parameter
    salt_water = M1.ins[0]
    base = salt_water.imass['NaCl']
    dist = shape.Uniform(lower=45, upper=55)
    @param(name='Salt flow rate', element=salt_water, kind='coupled',
           units='kg/hr', baseline=base, distribution=dist)
    def set_salt_flow(i):
        salt_water.imass['NaCl'] = i

    base = salt_water.price
    dist = shape.Normal(mu=0, sigma=1)
    @param(name='Salt solution price', element=salt_water, kind='cost',
           units='USD/kg', baseline=base, distribution=dist)
    def set_salt_solution_price(i):
        salt_water.price = i

    base = M1.tau
    dist = shape.Triangle(lower=base*0.9, midpoint=base, upper=base*1.1)
    @param(name='Mix tank retention time', element=M1, kind='cost',
           units='hr', baseline=base, distribution=dist)
    def set_mix_tank_tau(i):
        M1.tau = i

    base = M1.kW_per_m3
    dist = shape.Uniform(lower=base*0.8, upper=base*1.2)
    @param(name='Mix tank mixer power usage', element=M1, kind='cost',
           units='kW/m3', baseline=base, distribution=dist)
    def set_mix_tank_mixer_power(i):
        M1.kW_per_m3 = i

    base = P1.dP_design
    dist = shape.Triangle(lower=base*3/4, midpoint=base, upper=base/4*5)
    @param(name='Pump design head', element=P1, kind='design',
           units='kPa', baseline=base, distribution=dist)
    def set_pump_head(i):
        P1.dP_design = i

    base = H1.T
    dist = shape.Uniform(lower=300, upper=400)
    @param(name='Heat exchanger temperature', element=H1, kind='coupled',
           units='K', baseline=base, distribution=dist)
    def set_hx_temperature(i):
        H1.T = i

    # Add metrics
    metric = model.metric
    @metric(name='Total heating duty', units='kJ/yr', element=sys)
    def get_heating_duty():
        return sys.get_heating_duty()

    @metric(name='Total electricity consumption', units='kWh/yr', element=sys)
    def get_electricity_consumption():
        return sys.get_electricity_consumption()

    @metric(name='Total capital expenditure', units='USD', element=tea)
    def get_CAPEX():
        return tea.CAPEX

    @metric(name='Net present value', units='USD', element=tea)
    def get_NPV():
        return tea.NPV

    if evaluate:
        samples = model.sample(N=N, rule=rule, seed=seed, **sample_kwargs)
        model.load_samples(samples)
        model.evaluate()

    return model