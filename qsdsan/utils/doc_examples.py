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

from .. import (
    Component,
    Components,
    get_thermo,
    Model,
    SanStream,
    set_thermo as qs_set_thermo,
    TEA,
    System,
    )

__all__ = (
    'create_example_components',
    'create_example_system',
    'create_example_treatment_systems',
    'create_example_model',
    )


# %%

# =============================================================================
# Example components
# =============================================================================

def create_example_components(set_thermo=True):
    '''
    Load pre-constructed components for documentation purpose.

    Returns
    -------
    A :class:`qsdsan.CompiledComponents` object with components including
    H2O (alias Water), CO2, N2O, NaCl, H2SO4, CH4 (alias Methane), Methanol, and Ethanol.

    Examples
    --------
    >>> from qsdsan.utils import create_example_components
    >>> cmps = create_example_components()
    >>> cmps.show() # doctest: +ELLIPSIS
    CompiledComponents([
        H2O,   CO2, N2O,      NaCl,...
        H2SO4, CH4, Methanol, Ethanol,
    ])
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

    cmps.compile(ignore_inaccurate_molar_weight=True)
    cmps.set_synonym('H2O', 'Water')
    cmps.set_synonym('CH4', 'Methane')
    if set_thermo: qs_set_thermo(cmps)
    return cmps


# %%

# =============================================================================
# Example systems
# =============================================================================

def create_example_system(components=None):
    '''
    Load a pre-constructed system for documentation purpose.

    Returns
    -------
    A :class:`qsdsan.System` object including a mix tank, a pump, a heat exchanger,
    and various mixer/splitters.

    Parameters
    ----------
    components : obj
        If given, will call :func:`qsdsan.set_thermo(components)` to set components
        used in system simulation.
        The provided `components` must have "Water", "NaCl", "Methanol", and "Ethanol"
        components, which are used in the construction of the system.

    Examples
    --------
    >>> from qsdsan.utils import create_example_system
    >>> # Components from `create_example_components` will be loaded if no components are set/given
    >>> sys = create_example_system()
    >>> sys.path
    (<MixTank: M1>,
     <Pump: P1>,
     <HXutility: H1>,
     <ComponentSplitter: S1>,
     <Mixer: M2>,
     <Splitter: S2>)
    >>> sys.diagram() # doctest: +SKIP
    '''
    from .. import unit_operations as su

    if components: qs_set_thermo(components)
    else:
        try:
            thermo = get_thermo()
            if not ('H2O', 'Methanol', 'Ethanol', 'NaCl') in thermo.components.names:
                qs_set_thermo(create_example_components())
        except: qs_set_thermo(create_example_components())

    salt_water = SanStream('salt_water', Water=2000, NaCl=50, units='kg/hr')
    methanol = SanStream('methanol', Methanol=20, units='kg/hr')
    ethanol = SanStream('ethanol', Ethanol=10, units='kg/hr')

    M1 = su.MixTank('M1', ins=(salt_water, 'recycled_brine', methanol, ethanol))
    P1 = su.Pump('P1', dP_design=405300,ins=M1-0)
    H1 = su.HXutility('H1', ins=P1-0, T=350)
    S1 = su.ComponentSplitter('S1', ins=H1-0, split_keys=('Methanol', 'Ethanol'))
    M2 = su.Mixer('M2', ins=(S1-0, S1-1), outs='alcohols')
    S2 = su.Splitter('S2', ins=S1-2, outs=(1-M1, 'waste_brine'), split=0.2)
    sys = System('sys', path=(M1, P1, H1, S1, M2, S2))

    return sys


# %%

# =============================================================================
# Example wastewater treatment systems (shared by the TEA and LCA tutorials)
# =============================================================================

def create_example_treatment_systems(components=None, set_thermo=True):
    '''
    Build the aerobic and anaerobic wastewater treatment systems used in the
    TEA and LCA tutorials, for documentation purpose.

    Two single-unit systems treat the same municipal wastewater
    (4,000 m3/d, COD ~430 mg/L):

    - ``aer_sys``: an aerobic activated-sludge plant (spends electricity on
      aeration, produces waste sludge); and
    - ``ana_sys``: an anaerobic plant (recovers biogas, but needs heating and a
      little alkalinity).

    Both plants share an installed capital and a sludge-disposal cost (declared
    on a common abstract base class), so they are a compact but realistic
    substrate for techno-economic and life cycle analyses. The sizing, energy,
    and cost figures follow Metcalf & Eddy, *Wastewater Engineering* (5th ed.);
    they are order-of-magnitude teaching values, not a design basis.

    Parameters
    ----------
    components : obj
        If given, will call :func:`qsdsan.set_thermo(components)`; otherwise a
        small set of components is created. Must include "H2O", "O2", "CO2",
        "NH3", "CH4", "NaHCO3", "Substrate", and "Biomass".
    set_thermo : bool
        Whether to set the thermo property package to the components used here
        (ignored when `components` is None, in which case it is always set).

    Returns
    -------
    aer_sys : :class:`qsdsan.System`
        The aerobic treatment system (unit ``aer``).
    ana_sys : :class:`qsdsan.System`
        The anaerobic treatment system (unit ``ana``).

    Examples
    --------
    >>> from qsdsan.utils import create_example_treatment_systems
    >>> aer_sys, ana_sys = create_example_treatment_systems()
    >>> aer_sys.simulate()
    >>> ana_sys.simulate()
    >>> ([u.ID for u in aer_sys.units], [u.ID for u in ana_sys.units])
    (['aer'], ['ana'])
    >>> aer_sys.diagram() # doctest: +SKIP
    '''
    from .. import SanUnit, WasteStream, Component, Components, PowerUtility
    from . import get_digestion_rxns

    if components is not None:
        if set_thermo: qs_set_thermo(components)
        cmps = components
    else:
        def make_cmp(ID, formula=None, search_ID=None, phase='l', size='Soluble',
                     deg='Undegradable', org=False):
            return Component(ID, formula=formula, search_ID=search_ID, phase=phase,
                             particle_size=size, degradability=deg, organic=org)
        H2O = make_cmp('H2O', search_ID='H2O')
        O2  = make_cmp('O2',  search_ID='O2',  phase='g', size='Dissolved gas')
        CO2 = make_cmp('CO2', search_ID='CO2', phase='g', size='Dissolved gas')
        NH3 = make_cmp('NH3', search_ID='NH3', phase='g', size='Dissolved gas')
        CH4 = make_cmp('CH4', search_ID='CH4', phase='g', size='Dissolved gas',
                       deg='Readily', org=True)
        NaHCO3 = make_cmp('NaHCO3', search_ID='NaHCO3')
        Substrate = make_cmp('Substrate', formula='C10H19O3N', deg='Readily', org=True)
        Biomass = make_cmp('Biomass', formula='C5H7O2N', phase='s', size='Particulate',
                           deg='Slowly', org=True)
        cmps = Components([H2O, O2, CO2, NH3, CH4, NaHCO3, Substrate, Biomass])
        for c in (NaHCO3, Substrate, Biomass):
            c.copy_models_from(H2O, ('V', 'sigma', 'epsilon', 'kappa', 'Cn', 'mu'))
        cmps.compile(ignore_inaccurate_molar_weight=True)
        qs_set_thermo(cmps)

    Substrate = cmps.Substrate

    class TreatmentPlant(SanUnit, isabstract=True):
        '''Base plant: carries the installed and sludge-disposal costs.'''
        plant_capital = 5e6          # USD installed (M&E Ch. 4)
        sludge_disposal_cost = 0.10  # USD/kg dry solids (M&E Ch. 4)
        _F_BM_default = {'Plant': 1.}
        def _cost(self):
            self.baseline_purchase_costs['Plant'] = self.plant_capital

    class AerobicPlant(TreatmentPlant):
        '''Activated sludge: grows biomass and oxidizes the rest (aeration O2 = COD oxidized).'''
        _N_ins = 1; _N_outs = 2
        X_growth = 0.40; X_oxid = 0.95; O2_per_kWh = 1.2
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.growth_rxns = get_digestion_rxns(self.components, 0., self.X_growth, 'Biomass', 1.)
            self._mixed = WasteStream(f'{self.ID}_mixed')
        def _run(self):
            eff, sludge = self.outs
            m = self._mixed; m.copy_like(self.ins[0])
            self.growth_rxns(m.mol)
            oxidized = m.imass['Substrate']*self.X_oxid
            self._O2_demand = oxidized*self.components.Substrate.i_COD*24
            m.imass['Substrate'] -= oxidized
            sludge.empty(); sludge.phase = 's'; sludge.imass['Biomass'] = m.imass['Biomass']
            m.imass['Biomass'] = 0
            eff.copy_like(m)
        def _cost(self):
            super()._cost()
            self.power_utility(self._O2_demand/self.O2_per_kWh/24)
            self.add_OPEX = {'Sludge disposal': self.outs[1].F_mass*self.sludge_disposal_cost}

    class AnaerobicPlant(TreatmentPlant):
        '''Anaerobic digestion: recovers biogas, but needs heating and alkalinity.'''
        _N_ins = 2; _N_outs = 3
        X_biogas = 0.86; X_growth = 0.05
        CH4_LHV = 50000.; energy_price = 5/1e6
        T_op = 273.15 + 35; HX_eff = 0.80; NaHCO3_dose = 0.10
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.rxns = get_digestion_rxns(self.components, self.X_biogas, self.X_growth, 'Biomass', 1.)
            self._mixed = WasteStream(f'{self.ID}_mixed')
        def _run(self):
            ww, chem = self.ins
            eff, sludge, biogas = self.outs
            chem.empty(); chem.imass['NaHCO3'] = self.NaHCO3_dose*ww.F_vol
            m = self._mixed; m.copy_like(ww)
            self.rxns(m.mol)
            biogas.empty(); biogas.phase = 'g'
            biogas.imass['CH4'] = m.imass['CH4']; biogas.imass['CO2'] = m.imass['CO2']
            biogas.price = (biogas.imass['CH4']*self.CH4_LHV*self.energy_price)/biogas.F_mass \
                           if biogas.F_mass else 0.
            sludge.empty(); sludge.phase = 's'; sludge.imass['Biomass'] = m.imass['Biomass']
            m.imass['CH4'] = m.imass['CO2'] = m.imass['Biomass'] = m.imass['NH3'] = 0
            eff.copy_like(m)
        def _design(self):
            inf = self.ins[0]
            duty = inf.F_mass * inf.Cp * (self.T_op - inf.T)
            if duty > 0:
                self.add_heat_utility(duty, inf.T, T_out=self.T_op,
                                      heat_transfer_efficiency=self.HX_eff)
        def _cost(self):
            super()._cost()
            self.add_OPEX = {'Sludge disposal': self.outs[1].F_mass*self.sludge_disposal_cost}

    def make_influent(ID, COD=430.):
        ww = WasteStream(ID, T=273.15+20); ww.ivol['H2O'] = 4000/24   # 4,000 m3/d at 20 C
        ww.imass['Substrate'] = (COD/1000 * 4000/24)/Substrate.i_COD  # set the COD
        return ww

    aer = AerobicPlant('aer', ins=make_influent('ww_aer'),
                       outs=('aer_effluent', 'aer_sludge'))
    aer_sys = System('aer_sys', path=(aer,))

    NaHCO3_feed = WasteStream('NaHCO3_feed', price=0.90)   # USD/kg (M&E Ch. 10)
    ana = AnaerobicPlant('ana', ins=(make_influent('ww_ana'), NaHCO3_feed),
                         outs=('ana_effluent', 'ana_sludge', 'biogas'))
    ana_sys = System('ana_sys', path=(ana,))

    PowerUtility.price = 0.08   # USD/kWh electricity (M&E Ch. 4)
    return aer_sys, ana_sys


# %%

# =============================================================================
# Example system model
# =============================================================================

def create_example_model(evaluate=False, N=100, rule='L', seed=554, **sample_kwargs):
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
    >>> from qsdsan.utils import create_example_model
    >>> model = create_example_model(N=100, rule='L', seed=554, evaluate=False)
    >>> model.system.path
    (<MixTank: M1>,
     <Pump: P1>,
     <HXutility: H1>,
     <ComponentSplitter: S1>,
     <Mixer: M2>,
     <Splitter: S2>)
    >>> model.parameters # doctest: +SKIP
    (<Parameter: [Stream-salt water] Salt flow rate (kg/hr)>,
     <Parameter: [Stream-salt water] Salt solution price (USD/kg)>,
     <Parameter: [Mix tank-M1] Mix tank retention time (hr)>,
     <Parameter: [Mix tank-M1] Mix tank mixer power usage (kW/m3)>,
     <Parameter: [Pump-P1] Pump design head (kPa)>,
     <Parameter: [HXutility-H1] Heat exchanger temperature (K)>)
    >>> model.metrics # doctest: +SKIP
    (<Metric: [System] Total heating duty (kJ/yr)>,
     <Metric: [System] Total electricity consumption (kWh/yr)>,
     <Metric: [TEA] Total capital expenditure (USD)>,
     <Metric: [TEA] Net present value (USD)>)
    '''
    from chaospy import distributions as shape

    sys = create_example_system()
    M1, P1, H1, S1, M2, S2 = sys.path
    sys.simulate()
    tea = TEA(sys)
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
