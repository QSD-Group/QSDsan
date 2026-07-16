#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from biosteam.units.decorators import cost
from qsdsan import SanUnit, Stream, CEPCI_by_year
from qsdsan.utils import auom
from ._reactor import Reactor
from ..bst import IsothermalCompressor, HXutility, HXprocess

__all__ = ('Hydroprocessing',)

_lb_to_kg = auom('lb').conversion_factor('kg')
_m3perh_to_mmscfd = 1/1177.17 # H2


# %%

# =============================================================================
# Hydroprocessing (general)
# =============================================================================

@cost(basis='Oil lb flowrate', ID='Hydrocracker', units='lb/hr',
      cost=25e6, S=5963, CE=CEPCI_by_year[2007], n=0.75, BM=1)
@cost(basis='Oil lb flowrate', ID='Hydrotreater', units='lb/hr',
      cost=27e6, S=69637, CE=CEPCI_by_year[2007], n=0.68, BM=1)
@cost(basis='PSA H2 lb flowrate', ID='PSA', units='lb/hr',
      cost=1750000, S=5402, CE=CEPCI_by_year[2004], n=0.8, BM=2.47)
class Hydroprocessing(Reactor):
    '''
    General fuel-upgrading reactor (hydrocracking, hydrotreating, or similar):
    influent oil mixed with H2 reacts at elevated temperature and pressure to
    produce a single upgraded effluent stream (mixed gas/oil/aqueous phases,
    blended per `gas_yield`/`oil_yield`/`aq_yield` and their respective
    composition dicts). Co-product includes spent catalyst.

    **On reaction exotherm:** this class does NOT credit the hydroprocessing
    reaction's own exotherm against `inf_hx`'s heating duty by default --
    `inf_hx` costs the full sensible heat needed to bring both the influent
    oil and hydrogen up to `T` (see `inf_T` below for the opt-in mechanism to
    change this). This is a deliberate choice because this class's
    yields/compositions (`gas_yield`/`oil_yield`/`aqueous_composition`/etc.)
    are empirically-fit mass fractions across many named product species, not
    derived from a mass/atom-balanced stoichiometric reaction. Therefore,
    computing the reaction's real heat of reaction from
    each component's tabulated heat of formation (`self.Hnet`) is not
    reliable. If a more accurate influent temperature is known,
    supply it via `inf_T` instead.

    Parameters
    ----------
    ins : Iterable(stream)
        Influent oil, hydrogen, catalyst_in.
    outs : Iterable(stream)
        Effluent (oil + gas + aqueous, blended), catalyst_out.
    T : float
        Reaction temperature, [K].
    P : float
        Reaction pressure, [Pa]; also the hydrogen compressor's target pressure.
    inf_T : float or None
        Influent preheat temperature, [K]; if provided, both the oil and
        hydrogen streams (fresh + recycled) are heated only to `inf_T`
        (not `T`) before entering the reactor, and the remaining temperature
        rise from `inf_T` to `T` is treated as coming from the reaction's own
        heat of reaction (uncosted). Defaults to `None` (no credit; both streams effectively heated
        straight to `T`, matching this class's historical behavior).
    WHSV : float
        Weight hourly space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime : float
        Catalyst lifetime, [hr].
    catalyst_ID : str
        ID of the catalyst.
    hydrogen_rxned_to_inf_oil : float
        Reacted H2 to influent oil mass ratio.
    hydrogen_ratio : float
        Total hydrogen amount = hydrogen_rxned_to_inf_oil * hydrogen_ratio;
        excess hydrogen not recovered by the (optional) PSA leaves with the effluent.
    include_PSA : bool
        Whether to include a pressure swing adsorption (PSA) unit to recover excess H2.
    PSA_efficiency : float
        H2 recovery efficiency of the PSA unit; forced to 0 if `include_PSA` is False.
    gas_yield : float
        Mass ratio of fuel gas to the sum of influent oil and reacted H2.
    oil_yield : float
        Mass ratio of treated oil to the sum of influent oil and reacted H2.
    gas_composition : dict
        Composition of the gas fraction (excluding excess H2), normalized to 100% sum.
    oil_composition : dict
        Composition of the treated-oil fraction, normalized to 100% sum.
    aqueous_composition : dict
        Composition of the aqueous fraction, normalized to 100% sum.
    internal_heat_exchanging : bool
        If True, use the effluent to preheat the influent before the reactor.
    use_decorated_cost : str
        `'Hydrocracker'` or `'Hydrotreater'` to use the corresponding literature
        decorated cost; any other value uses generic `Reactor`/`PressureVessel` costing.

    Examples
    --------
    >>> from qsdsan import Component, Components, set_thermo, Stream
    >>> from qsdsan.unit_operations import Hydroprocessing
    >>> catalyst = Component('Catalyst', phase='s', particle_size='Particulate',
    ...                       degradability='Undegradable', organic=False, formula='Al2O3')
    >>> _ = catalyst.default()
    >>> cmps = Components([
    ...     Component('Octane', search_ID='Octane', particle_size='Soluble',
    ...                degradability='Slowly', organic=True),
    ...     Component('CH4', phase='g', particle_size='Dissolved gas',
    ...                degradability='Slowly', organic=True),
    ...     Component('CO2', phase='g', particle_size='Dissolved gas',
    ...                degradability='Undegradable', organic=False),
    ...     Component('H2', phase='g', particle_size='Dissolved gas',
    ...                degradability='Undegradable', organic=False),
    ...     Component('H2O', particle_size='Soluble', degradability='Undegradable',
    ...                organic=False),
    ...     catalyst,
    ...     ])
    >>> cmps.compile()
    >>> set_thermo(cmps)
    >>> oil = Stream('oil', Octane=100, units='kg/hr')
    >>> H2 = Stream('H2')
    >>> catalyst_in = Stream('catalyst_in', Catalyst=1, units='kg/hr', phase='s')
    >>> U1 = Hydroprocessing('U1', ins=(oil, H2, catalyst_in), outs=('eff', 'catalyst_out'),
    ...                      catalyst_ID='Catalyst',
    ...                      gas_composition={'CO2': 0.5, 'CH4': 0.5},
    ...                      oil_composition={'Octane': 1},
    ...                      aqueous_composition={'H2O': 1})
    >>> U1.simulate()
    >>> round(U1.oil_yield + U1.gas_yield + U1.aq_yield, 6)
    1.0
    >>> sorted(set(type(u).__name__ for u in U1.auxiliary_units))
    ['HXprocess', 'HXutility', 'IsothermalCompressor']

    See Also
    --------
    `saf systems <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/saf>`_

    References
    ----------
    [1] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.;
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.;
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C.
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    '''
    _N_ins = 3
    _N_outs = 2
    _units = {
        'Oil lb flowrate': 'lb/hr',
        'PSA H2 lb flowrate': 'lb/hr',
        }
    _F_BM_default = {
        **Reactor._F_BM_default,
        'Heat exchanger': 3.17,
        'Compressor': 1.1,
        }
    auxiliary_unit_names = ('compressor', 'hx', 'inf_hx', 'inf_hx_H2')

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 include_construction=False,
                 T=451+273.15,
                 P=1034.7*6894.76,
                 inf_T=None,
                 WHSV=0.625,
                 catalyst_lifetime=5*7920,
                 catalyst_ID='HC_catalyst',
                 hydrogen_rxned_to_inf_oil=0.01125,
                 hydrogen_ratio=5.556,
                 include_PSA=False,
                 PSA_efficiency=0.9,
                 gas_yield=0.03880+0.00630, # gas_yield = 1 - oil_yield
                 oil_yield=1-0.03880-0.00630,
                 gas_composition={'CO2': 0.03880, 'CH4': 0.00630},
                 oil_composition={
                     'CYCHEX': 0.03714, 'HEXANE': 0.01111,
                     'HEPTANE': 0.11474, 'OCTANE': 0.08125,
                     'C9H20': 0.09086, 'C10H22': 0.11756,
                     'C11H24': 0.16846, 'C12H26': 0.13198,
                     'C13H28': 0.09302, 'C14H30': 0.04643,
                     'C15H32': 0.03250, 'C16H34': 0.01923,
                     'C17H36': 0.00431, 'C18H38': 0.00099,
                     'C19H40': 0.00497, 'C20H42': 0.00033,
                     },
                 aqueous_composition={'H2O': 1},
                 internal_heat_exchanging=True,
                 use_decorated_cost='Hydrocracker',
                 dynamic_V_wf=False,
                 tau=5, V_wf=0.4,
                 length_to_diameter=2, diameter=None,
                 N=None, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1.5,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 ):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, include_construction=include_construction)
        self.T = T
        self.P = P
        self.inf_T = inf_T
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.catalyst_ID = catalyst_ID
        self.hydrogen_rxned_to_inf_oil = hydrogen_rxned_to_inf_oil
        self.hydrogen_ratio = hydrogen_ratio
        self.include_PSA = include_PSA
        self.PSA_efficiency = PSA_efficiency
        self.gas_yield = gas_yield
        self.oil_yield = oil_yield
        self.gas_composition = gas_composition
        self.oil_composition = oil_composition
        self.aqueous_composition = aqueous_composition
        self.internal_heat_exchanging = internal_heat_exchanging
        IC_in = Stream(f'{ID}_IC_in')
        IC_out = Stream(f'{ID}_IC_out')
        self.compressor = IsothermalCompressor(ID=f'.{ID}_IC', ins=IC_in, outs=IC_out, P=P)
        inf_pre_hx = Stream(f'{ID}_inf_pre_hx')
        eff_pre_hx = Stream(f'{ID}_eff_pre_hx')
        inf_after_hx = Stream(f'{ID}_inf_after_hx')
        eff_after_hx = Stream(f'{ID}_eff_after_hx')
        self.hx = HXprocess(ID=f'.{ID}_hx', ins=(inf_pre_hx, eff_pre_hx), outs=(inf_after_hx, eff_after_hx))
        inf_hx_out = Stream(f'{ID}_inf_hx_out')
        self.inf_hx = HXutility(ID=f'.{ID}_inf_hx', ins=inf_after_hx, outs=inf_hx_out, T=T, rigorous=True)
        inf_hx_H2_out = Stream(f'{ID}_inf_hx_H2_out')
        self.inf_hx_H2 = HXutility(ID=f'.{ID}_inf_hx_H2', ins=Stream(f'{ID}_inf_hx_H2_in'), outs=inf_hx_H2_out, T=T, rigorous=True)
        self.use_decorated_cost = use_decorated_cost
        self.dynamic_V_wf = dynamic_V_wf
        self.tau = tau
        self.V_wf = V_wf
        self.length_to_diameter = length_to_diameter
        self.diameter = diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type

    def _run(self):
        inf_oil, makeup_hydrogen, catalyst_in = self.ins
        eff_oil, catalyst_out = self.outs

        catalyst_in.imass[self.catalyst_ID] = inf_oil.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)

        hydrogen_rxned_to_inf_oil = self.hydrogen_rxned_to_inf_oil
        hydrogen_ratio = self.hydrogen_ratio
        H2_rxned = inf_oil.F_mass * hydrogen_rxned_to_inf_oil
        H2_tot = H2_rxned * hydrogen_ratio
        H2_residual = self._H2_residual = H2_tot - H2_rxned
        H2_recycled = self._H2_recycled = H2_residual * self.PSA_efficiency
        H2_wasted = H2_residual - H2_recycled

        eff_oil.copy_like(inf_oil)
        eff_oil.phase = inf_oil.phase
        eff_oil.empty()
        eff_oil.imass[self.eff_composition.keys()] = self.eff_composition.values()
        eff_oil.F_mass = inf_oil.F_mass*(1 + hydrogen_rxned_to_inf_oil)
        eff_oil.imass['H2'] = H2_wasted
        eff_oil.P = self.P
        eff_oil.T = self.T
        eff_oil.vle(T=eff_oil.T, P=eff_oil.P)

        makeup_hydrogen.imass['H2'] = H2_rxned + H2_wasted
        makeup_hydrogen.phase = 'g'

    def _design(self):
        Design = self.design_results
        Design.clear()
        Design['PSA H2 lb flowrate'] = self._H2_residual / _lb_to_kg
        # Keep `recycled` in kg/hr for the mass balance below; only convert
        # to lb/hr for the Design report entry.
        recycled = self._H2_recycled
        Design['H2 recycled'] = recycled / _lb_to_kg
        Design['Oil lb flowrate'] = self.ins[0].F_mass/_lb_to_kg

        IC = self.compressor
        H2 = self.ins[1]
        IC_ins0, IC_outs0 = IC.ins[0], IC.outs[0]
        IC_ins0.copy_like(H2)
        IC_ins0.F_mass += recycled
        IC_outs0.copy_like(IC_ins0)
        IC_outs0.P = IC.P = self.P
        IC_ins0.phase = IC_outs0.phase = 'g'
        IC.simulate()

        hx = self.hx
        inf_hx = self.inf_hx
        inf_hx_in, inf_hx_out = inf_hx.ins[0], inf_hx.outs[0]
        inf_pre_hx, eff_pre_hx = hx.ins
        inf_after_hx, eff_after_hx = hx.outs
        inf_pre_hx.copy_like(self.ins[0])
        eff_pre_hx.copy_like(self.outs[0])

        if self.internal_heat_exchanging:
            hx.simulate()
        else:
            hx.empty()
            inf_after_hx.copy_like(inf_pre_hx)
            eff_after_hx.copy_like(eff_pre_hx)

        if self.inf_T is not None:
            # Preheat oil and hydrogen (fresh + recycled, already at the
            # compressor's outlet pressure) SEPARATELY, both only to `inf_T`
            # -- the remaining rise from `inf_T` to `T` is treated as coming
            # from the reaction's own exotherm (uncosted); see the class
            # docstring for why `inf_T` is a caller-supplied literature value
            # rather than something computed from thermodynamic data. Two
            # separate exchangers (matching the old, pre-refactor model's
            # topology) rather than one combined mixed-stream exchanger,
            # since mixing oil and H2 before heating produces a real,
            # non-additive lower total duty than heating them apart --
            # confirmed via diagnostics -- and both need to independently
            # participate in whole-plant heat-exchanger-network integration.
            inf_hx_in.copy_like(inf_after_hx)
            inf_hx_out.copy_flow(inf_hx_in)
            inf_hx_out.T = self.inf_T
            inf_hx_out.P = self.P
            inf_hx.simulate_as_auxiliary_exchanger(ins=inf_hx.ins, outs=inf_hx.outs)

            h2_hx = self.inf_hx_H2
            h2_hx_in, h2_hx_out = h2_hx.ins[0], h2_hx.outs[0]
            h2_hx_in.copy_like(IC_outs0)
            h2_hx_out.copy_flow(h2_hx_in)
            # `copy_flow` doesn't copy phase; simulate_as_auxiliary_exchanger's
            # outs-given/duty=None path computes duty from outlet.Hnet directly
            # without re-resolving VLE, so a stale (default liquid) phase on a
            # 100% H2 stream would badly corrupt this exchanger's duty.
            h2_hx_out.phase = h2_hx_in.phase
            h2_hx_out.T = self.inf_T
            h2_hx_out.P = self.P
            h2_hx.simulate_as_auxiliary_exchanger(ins=h2_hx.ins, outs=h2_hx.outs)
        else:
            inf_hx_in.copy_like(inf_after_hx)
            inf_hx_out.copy_flow(inf_hx_in)
            inf_hx_out.T = self.T
            inf_hx_out.P = self.P
            inf_hx.simulate_as_auxiliary_exchanger(ins=inf_hx.ins, outs=inf_hx.outs)
            self.inf_hx_H2.empty() # unused for this path (e.g. saf's own inf_T=None use)

        if self.dynamic_V_wf:
            # Old, pre-refactor model discounted working volume for H2 gas's
            # disproportionate volume vs. liquid oil rather than using a flat
            # constant; confirmed via diagnostics to reproduce the old
            # model's installed cost almost exactly when combined with
            # itemized (non-decorated) reactor costing.
            V_H2 = H2.F_vol/self.hydrogen_ratio*101325/self.P
            V_oil = self.ins[0].F_vol
            self.V_wf = 0.4*V_oil/(V_oil + V_H2)

        Reactor._design(self)

    def _cost(self):
        Cost = self.baseline_purchase_costs
        Cost.clear()

        use_decorated_cost = self.use_decorated_cost
        include_PSA = self.include_PSA
        self._decorated_cost()

        if use_decorated_cost == 'Hydrocracker':
            Cost.pop('Hydrotreater')
        elif use_decorated_cost == 'Hydrotreater':
            Cost.pop('Hydrocracker')
        else:
            Cost.pop('Hydrocracker')
            Cost.pop('Hydrotreater')
            Reactor._cost(self)

        if not include_PSA: Cost.pop('PSA')

    def _normalize_composition(self, dct):
        total = sum(dct.values())
        if total <= 0: raise ValueError(f'Sum of total yields/composition should be positive, not {total}.')
        return {k: v/total for k, v in dct.items()}

    def _normalize_yields(self):
        gas = self._gas_yield
        oil = self._oil_yield
        gas_oil = gas + oil
        aq = 0
        if gas_oil > 1:
            gas /= gas_oil
            oil /= gas_oil
        else:
            aq = 1 - gas_oil
        self._gas_yield = gas
        self._oil_yield = oil
        self._aq_yield = aq

    @property
    def gas_yield(self):
        '''[float] Mass ratio of fuel gas to the sum of influent oil and reacted H2.
        Setting this triggers renormalization: if gas_yield + oil_yield > 1, both are
        rescaled proportionally; otherwise aq_yield absorbs the remainder.'''
        return self._gas_yield
    @gas_yield.setter
    def gas_yield(self, gas):
        self._gas_yield = gas
        if hasattr(self, '_oil_yield'):
            self._normalize_yields()

    @property
    def oil_yield(self):
        '''[float] Mass ratio of treated oil to the sum of influent oil and reacted H2.
        Setting this triggers renormalization: if gas_yield + oil_yield > 1, both are
        rescaled proportionally; otherwise aq_yield absorbs the remainder.'''
        return self._oil_yield
    @oil_yield.setter
    def oil_yield(self, oil):
        self._oil_yield = oil
        if hasattr(self, '_gas_yield'):
            self._normalize_yields()

    @property
    def aq_yield(self):
        '''[float] Mass ratio of aqueous phase to the sum of influent oil and reacted H2.
        Read-only; derived as 1 - gas_yield - oil_yield (or 0 if their sum > 1).'''
        return self._aq_yield

    @property
    def gas_composition(self):
        return self._gas_composition
    @gas_composition.setter
    def gas_composition(self, comp_dct):
        self._gas_composition = self._normalize_composition(comp_dct)

    @property
    def oil_composition(self):
        return self._oil_composition
    @oil_composition.setter
    def oil_composition(self, comp_dct):
        self._oil_composition = self._normalize_composition(comp_dct)

    @property
    def aqueous_composition(self):
        return self._aqueous_composition
    @aqueous_composition.setter
    def aqueous_composition(self, comp_dct):
        self._aqueous_composition = self._normalize_composition(comp_dct)

    @property
    def eff_composition(self):
        '''[dict] Composition of the blended effluent, normalized to 100% sum.'''
        gas_composition = self.gas_composition
        oil_composition = self.oil_composition
        aqueous_composition = self.aqueous_composition
        oil_yield = self.oil_yield
        gas_yield = self.gas_yield
        aq_yield = self.aq_yield
        eff_composition = {k: v*gas_yield for k, v in gas_composition.items()}
        eff_composition.update({k: v*oil_yield for k, v in oil_composition.items()})
        eff_composition.update({k: v*aq_yield for k, v in aqueous_composition.items()})
        return self._normalize_composition(eff_composition)

    @property
    def PSA_efficiency(self):
        '''[float] H2 recovery efficiency of the PSA unit, forced to 0 if `include_PSA` is False.'''
        if self.include_PSA: return self._PSA_efficiency
        return 0
    @PSA_efficiency.setter
    def PSA_efficiency(self, i):
        if i > 1: raise ValueError('PSA_efficiency cannot be larger than 1.')
        self._PSA_efficiency = i
