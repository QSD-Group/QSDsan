#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>

    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from biosteam.units.decorators import cost
from qsdsan import SanUnit, Stream, CEPCI_by_year
from qsdsan.utils import auom
from ._reactor import Reactor
from ..bst import HXutility, HXprocess
from ..dynamic import Pump

__all__ = (
    'CatalyticHydrothermalGasification',
    'HydrothermalLiquefaction',
    'KnockOutDrum',
    )

_lb_to_kg = auom('lb').conversion_factor('kg')
_m3_to_gal = auom('m3').conversion_factor('gallon')
_in_to_m = auom('inch').conversion_factor('m')

# %%

# =============================================================================
# CHG
# =============================================================================

@cost(basis='Treatment capacity', ID='Hydrocyclone', units='lb/h',
      cost=5000000, S=968859,
      CE=CEPCI_by_year[2009], n=0.65, BM=2.1)
class CatalyticHydrothermalGasification(Reactor):
    '''
    CHG serves to reduce the COD content in the aqueous phase and produce fuel
    gas under elevated temperature (350°C) and pressure. The outlet will be
    cooled down and separated by a flash unit.
    
    Parameters
    ----------
    ins : Iterable(stream)
        chg_in, catalyst_in.
    outs : Iterable(stream)
        chg_out, catalyst_out.
    pump_pressure: float
        CHG influent pressure, [Pa].
    heat_temp: float
        CHG influent temperature, [K].
    cool_temp: float
        CHG effluent temperature, [K].
    WHSV: float
        Weight Hourly Space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime: float
        CHG catalyst lifetime, [hr].
    gas_composition: dict
        CHG gas composition.
    gas_C_2_total_C: dict
        CHG gas carbon content to feed carbon content.
    CAPEX_factor: float
        Factor used to adjust CAPEX.
        
    References
    ----------
    [1] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    [2] Davis, R. E.; Grundl, N. J.; Tao, L.; Biddy, M. J.; Tan, E. C.;
        Beckham, G. T.; Humbird, D.; Thompson, D. N.; Roni, M. S. Process
        Design and Economics for the Conversion of Lignocellulosic Biomass
        to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case
        Update; Biochemical Deconstruction and Conversion of Biomass to Fuels
        and Products via Integrated Biorefinery Pathways; NREL/TP--5100-71949,
        1483234; 2018; p NREL/TP--5100-71949, 1483234.
        https://doi.org/10.2172/1483234.
    [3] Elliott, D. C.; Neuenschwander, G. G.; Hart, T. R.; Rotness, L. J.;
        Zacher, A. H.; Santosa, D. M.; Valkenburg, C.; Jones, S. B.;
        Rahardjo, S. A. T. Catalytic Hydrothermal Gasification of Lignin-Rich
        Biorefinery Residues and Algae Final Report. 87.
    '''
    _N_ins = 2
    _N_outs = 2
    
    _F_BM_default = {**Reactor._F_BM_default,
                      'Heat exchanger': 3.17,
                      'Sulfur guard': 2.0}
    _units= {'Treatment capacity': 'lb/h', # hydrocyclone
              'Hydrocyclone weight': 'lb'}
    
    auxiliary_unit_names=('pump','heat_ex_heating','heat_ex_cooling')
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  init_with='Stream',
                  pump_pressure=3089.7*6894.76,
                  heat_temp=350+273.15,
                  cool_temp=60+273.15,
                  WHSV=3.562,
                  catalyst_lifetime=7920, # 1 year [1]
                  gas_composition={'CH4':0.527,
                                   'CO2':0.432,
                                   'C2H6':0.011,
                                   'C3H8':0.030,
                                   'H2':0.0001}, # [1]
                  gas_C_2_total_C=0.5981, # [1]
                  P=None, tau=20/60, void_fraction=0.5, # [2, 3]
                  length_to_diameter=2, diameter=None,
                  N=6, V=None, auxiliary=False,
                  mixing_intensity=None, kW_per_m3=0,
                  wall_thickness_factor=1,
                  vessel_material='Stainless steel 316',
                  vessel_type='Vertical',
                  CAPEX_factor=1):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        
        self.pump_pressure = pump_pressure
        self.heat_temp = heat_temp
        self.cool_temp = cool_temp
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.gas_composition = gas_composition
        self.gas_C_2_total_C = gas_C_2_total_C
        pump_in = Stream(f'{ID}_pump_in')
        pump_out = Stream(f'{ID}_pump_out')
        self.pump = Pump(ID=f'.{ID}_pump', ins=pump_in, outs=pump_out, P=pump_pressure)
        hx_ht_in = Stream(f'{ID}_hx_ht_in')
        hx_ht_out = Stream(f'{ID}_hx_ht_out')
        self.heat_ex_heating = HXutility(ID=f'.{ID}_hx_ht', ins=hx_ht_in, outs=hx_ht_out, T=heat_temp, rigorous=True)
        hx_cl_in = Stream(f'{ID}_hx_cl_in')
        hx_cl_out = Stream(f'{ID}_hx_cl_out')
        self.heat_ex_cooling = HXutility(ID=f'.{ID}_hx_cl', ins=hx_cl_in, outs=hx_cl_out, T=cool_temp, rigorous=True)
        self.P = P
        self.tau = tau
        self.V_wf = void_fraction
        # no headspace, gases produced will be vented, so V_wf = void fraction [2, 3]
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
        self.CAPEX_factor = CAPEX_factor
        
    def _run(self):
        
        chg_in, catalyst_in = self.ins
        chg_out, catalyst_out = self.outs
        
        catalyst_in.imass['CHG_catalyst'] = chg_in.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        # catalysts amount is quite low compared to the main stream, therefore do not consider
        # heating/cooling of catalysts

        chg_out.phase = 'g'

        cmps = self.components
        gas_C_ratio = 0
        for name, ratio in self.gas_composition.items():
            gas_C_ratio += ratio*cmps[name].i_C
            
        gas_mass = chg_in.imass['C']*self.gas_C_2_total_C/gas_C_ratio
        
        for name,ratio in self.gas_composition.items():
            chg_out.imass[name] = gas_mass*ratio
                
        chg_out.imass['H2O'] = chg_in.F_mass - gas_mass
        # all C, N, and P are accounted in H2O here, but will be calculated as properties.

        chg_out.T = self.cool_temp
        chg_out.P = self.pump_pressure

        chg_out.vle(T=chg_out.T, P=chg_out.P)

    @property
    def CHGout_C(self):
        # not include carbon in gas phase
        return self.ins[0].imass['C']*(1 - self.gas_C_2_total_C)
    
    @property
    def CHGout_N(self):
        return self.ins[0].imass['N']
    
    @property
    def CHGout_P(self):
        return self.ins[0].imass['P']
        
    def _design(self):
        Design = self.design_results
        Design['Treatment capacity'] = self.ins[0].F_mass/_lb_to_kg
        
        pump = self.pump
        pump.ins[0].copy_like(self.ins[0])
        pump.simulate()
        
        hx_ht = self.heat_ex_heating
        hx_ht_ins0, hx_ht_outs0 = hx_ht.ins[0], hx_ht.outs[0]
        hx_ht_ins0.copy_like(self.ins[0])
        hx_ht_outs0.copy_like(hx_ht_ins0)
        hx_ht_ins0.T = self.ins[0].T
        hx_ht_outs0.T = hx_ht.T
        hx_ht_ins0.P = hx_ht_outs0.P = pump.P
        
        hx_ht_ins0.vle(T=hx_ht_ins0.T, P=hx_ht_ins0.P)
        hx_ht_outs0.vle(T=hx_ht_outs0.T, P=hx_ht_outs0.P)
        
        hx_ht.simulate_as_auxiliary_exchanger(ins=hx_ht.ins, outs=hx_ht.outs)
            
        hx_cl = self.heat_ex_cooling
        hx_cl_ins0, hx_cl_outs0 = hx_cl.ins[0], hx_cl.outs[0]
        hx_cl_ins0.copy_like(self.outs[0])
        hx_cl_outs0.copy_like(hx_cl_ins0)
        hx_cl_ins0.T = hx_ht.T
        hx_cl_outs0.T = hx_cl.T
        hx_cl_ins0.P = hx_cl_outs0.P = self.outs[0].P

        hx_cl_ins0.vle(T=hx_cl_ins0.T, P=hx_cl_ins0.P)
        hx_cl_outs0.vle(T=hx_cl_outs0.T, P=hx_cl_outs0.P)        
        
        hx_cl.simulate_as_auxiliary_exchanger(ins=hx_cl.ins, outs=hx_cl.outs)

        self.P = self.pump_pressure
        Reactor._design(self)
        Design['Hydrocyclone weight'] = 0.3*Design['Weight']*Design['Number of reactors'] # assume stainless steel
        # based on [1], page 54, the purchase price of hydrocyclone to the purchase price of CHG
        # reactor is around 0.3, therefore, assume the weight of hydrocyclone is 0.3*single CHG weight*number of CHG reactors
        self.construction[0].quantity += Design['Hydrocyclone weight']*_lb_to_kg
    
    def _cost(self):
        Reactor._cost(self)
        purchase_costs = self.baseline_purchase_costs
        current_cost = 0 # cost w/o sulfur guard
        for item in purchase_costs.keys():
            current_cost += purchase_costs[item]
        purchase_costs['Sulfur guard'] = current_cost*0.05
        self._decorated_cost()
        
        purchase_costs = self.baseline_purchase_costs
        for item in purchase_costs.keys():
            purchase_costs[item] *= self.CAPEX_factor
        
# %%

# =============================================================================
# KOdrum
# =============================================================================

class KnockOutDrum(Reactor):
    '''
    Knockout drum is an auxiliary unit for :class:`HydrothermalLiquefaction`,
    used when its cost is calculated using generic pressure vessel algorithms
    (i.e., `HydrothermalLiquefaction`'s `use_decorated_cost` is False).

    Parameters
    ----------
    drum_cost_factor : float
        Cost multiplier applied to the vessel's own baseline purchase cost
        (on top of `vessel_material`'s material factor, which biosteam's
        `PressureVessel` machinery already applies automatically), to match
        the pre-refactor sludge model's own `drum_steel_cost_factor`.
        Defaults to 1.5 -- see [1], page 54: fully scaling the reference
        report's factor from 2000 to 100 tons/day would need ~3, but that is
        too high, so 1.5 is used instead.

    See Also
    --------
    :class:`HydrothermalLiquefaction`

    :class:`Reactor`

    :class:`biosteam.units.design_tools.PressureVessel`

    `saf systems <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/saf>`_

    References
    ----------
    [1] Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels via
        Liquefaction - Hydrothermal Liquefaction Reactor Design: April 5, 2013;
        NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462, 1111191.
        https://doi.org/10.2172/1111191.
    '''
    _N_ins = 3
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', include_construction=False,
                 P=3049.7*6894.76, tau=0, V_wf=0,
                 length_to_diameter=2, diameter=None,
                 N=4, V=None,
                 auxiliary=True,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 drum_cost_factor=1.5,
                 ):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with,
                          include_construction=include_construction)
        self.P = P
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
        self.drum_cost_factor = drum_cost_factor

    def _run(self):
        pass

    def _cost(self):
        # `self.F_M` is NOT usable for this: `vessel_material`'s setter
        # (biosteam.units.design_tools.pressure_vessel.PressureVessel)
        # unconditionally overwrites F_M['Vertical/Horizontal pressure
        # vessel'] with the plain per-material factor as a side effect of
        # being assigned. Applying as a multiplier on baseline_purchase_costs instead.
        Reactor._cost(self)
        purchase_costs = self.baseline_purchase_costs
        purchase_costs[f'{self.vessel_type} pressure vessel'] *= self.drum_cost_factor

# =============================================================================
# HTL
# =============================================================================

@cost(basis='Wet mass flowrate', ID='HTL system', units='lb/hr',
      cost=37486757, S=574476,
      CE=CEPCI_by_year[2011], n=0.77, BM=2.1)
@cost(basis='Wet mass flowrate', ID='Solids filter oil/water separator', units='lb/hr',
      cost=3945523, S=574476,
      CE=CEPCI_by_year[2011], n=0.68, BM=1.9)
@cost(basis='Wet mass flowrate', ID='Hot oil system', units='lb/hr',
      cost=4670532, S=574476,
      CE=CEPCI_by_year[2011], n=0.6, BM=1.4)
class HydrothermalLiquefaction(Reactor):
    '''
    HTL converts feedstock to gas, aqueous, biocrude, and (hydro)char under
    elevated temperature and pressure. Product yields and compositions are
    directly specified via `dw_yields` and the four `*_composition` dicts
    (no feedstock-composition correlation is modeled here — subclass this
    unit to add one, e.g. a sludge-biochemical-composition correlation).

    Parameters
    ----------
    ins : Iterable(stream)
        Feedstock into HTL.
    outs : Iterable(stream)
        Gas, aqueous, biocrude, char.
    T : float
        Temperature of the HTL reaction, [K].
    P : float
        Pressure when the reaction is at temperature, [Pa].
    dw_yields : dict
        Dry weight percentage yields of the four products (gas, aqueous, biocrude, char),
        normalized to 100% sum. Keys must be 'gas', 'aqueous', 'biocrude', and 'char'.
    gas_composition : dict
        Composition of the gaseous products including water, normalized to 100% sum.
    aqueous_composition : dict
        Composition of the aqueous products excluding water, normalized to 100% sum.
        Water not allocated to other products all goes to aqueous.
    biocrude_composition : dict
        Composition of the biocrude products including water, normalized to 100% sum.
    char_composition : dict
        Composition of the char products including water, normalized to 100% sum.
    internal_heat_exchanging : bool
        If True, use the product to preheat the feedstock.
    eff_T : float
        HTL effluent temperature, [K]; if provided, an additional HX controls
        effluent temperature.
    eff_P : float
        HTL effluent pressure, [Pa].
    use_decorated_cost : bool
        If True, use cost scaled per [1]; otherwise use generic `Reactor`
        (`PressureVessel`) costing.
    F_M : dict
        Material factors used to adjust cost (only used when `use_decorated_cost` is False).

    Examples
    --------
    >>> from qsdsan import Components, WasteStream, set_thermo
    >>> from qsdsan.unit_operations import HydrothermalLiquefaction
    >>> cmps = Components.load_default()
    >>> set_thermo(cmps)
    >>> feed = WasteStream('htl_feed', S_F=200, Water=800, units='kg/hr')
    >>> HTL = HydrothermalLiquefaction(
    ...     'HTL', ins=feed, outs=('gas', 'aq', 'crude', 'char'),
    ...     dw_yields={'gas': 0.05, 'aqueous': 0.15, 'biocrude': 0.4, 'char': 0.4},
    ...     gas_composition={'S_CH4': 0.5, 'S_H2': 0.5},
    ...     aqueous_composition={'Water': 1},
    ...     biocrude_composition={'S_F': 0.5, 'Water': 0.5},
    ...     char_composition={'Water': 1},
    ...     T=280+273.15,
    ...     internal_heat_exchanging=False,
    ...     )
    >>> HTL.simulate()
    >>> gas, aq, crude, char = HTL.outs
    >>> round(crude.imass['S_F'], 2)  # kg/hr
    40.0
    >>> sorted(type(u).__name__ for u in HTL.auxiliary_units)
    ['HXprocess', 'HXutility', 'HXutility']

    See Also
    --------
    :class:`KnockOutDrum`

    :class:`Reactor`

    :class:`biosteam.units.design_tools.PressureVessel`

    `saf systems <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/saf>`_

    References
    ----------
    [1] Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels
        via Liquefaction - Hydrothermal Liquefaction Reactor Design:
        April 5, 2013; NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462,
        1111191. https://doi.org/10.2172/1111191.
    '''
    _N_ins = 1
    _N_outs = 4
    _units = {
        'Wet mass flowrate': 'lb/hr',
        'Solid filter and separator weight': 'lb',
        }

    auxiliary_unit_names = ('hx', 'inf_hx', 'eff_hx', 'kodrum')

    _F_BM_default = {
        **Reactor._F_BM_default,
        'Heat exchanger': 3.17,
        }

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', include_construction=False,
                 T=280+273.15,
                 P=101325,
                 dw_yields={
                     'gas': 0,
                     'aqueous': 0,
                     'biocrude': 0,
                     'char': 1,
                     },
                 gas_composition={'HTLgas': 1},
                 aqueous_composition={'HTLaqueous': 1},
                 biocrude_composition={'HTLbiocrude': 1},
                 char_composition={'HTLchar': 1},
                 internal_heat_exchanging=True,
                 eff_T=60+273.15,
                 eff_P=30*6894.76,
                 use_decorated_cost=True,
                 tau=15/60, V_wf=0.45,
                 length_to_diameter=None,
                 diameter=6.875*_in_to_m,
                 N=4, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Horizontal',
                 F_M={
                     'Horizontal pressure vessel': 2.7,
                     'Vertical pressure vessel': 2.7,
                     },
                 ):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with,
                          include_construction=include_construction)
        self.T = T
        self.P = P
        self.dw_yields = dw_yields
        self.gas_composition = gas_composition
        self.aqueous_composition = aqueous_composition
        self.biocrude_composition = biocrude_composition
        self.char_composition = char_composition
        self.internal_heat_exchanging = internal_heat_exchanging
        inf_pre_hx = Stream(f'{ID}_inf_pre_hx')
        eff_pre_hx = Stream(f'{ID}_eff_pre_hx')
        inf_after_hx = Stream(f'{ID}_inf_after_hx')
        eff_after_hx = Stream(f'{ID}_eff_after_hx')
        self.hx = HXprocess(ID=f'.{ID}_hx', ins=(inf_pre_hx, eff_pre_hx), outs=(inf_after_hx, eff_after_hx))
        inf_hx_out = Stream(f'{ID}_inf_hx_out')
        self.inf_hx = HXutility(ID=f'.{ID}_inf_hx', ins=inf_after_hx, outs=inf_hx_out, T=T, rigorous=True)
        self._inf_at_temp = Stream(f'{ID}_inf_at_temp')
        self._eff_at_temp = Stream(f'{ID}_eff_at_temp')
        eff_hx_out = Stream(f'{ID}_eff_hx_out')
        self.eff_T = eff_T
        self.eff_P = eff_P
        self.eff_hx = HXutility(ID=f'.{ID}_eff_hx', ins=eff_after_hx, outs=eff_hx_out, T=eff_T, rigorous=True)
        self.use_decorated_cost = use_decorated_cost
        if not use_decorated_cost:
            self.kodrum = KnockOutDrum(ID=f'.{ID}_KOdrum', include_construction=include_construction)
        else:
            self.kodrum = None
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
        # See KnockOutDrum.__init__'s comment: F_M must be assigned before
        # vessel_material, whose setter mutates F_M in place.
        self.F_M = F_M
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type

    def _run(self):
        feed = self.ins[0]
        gas, aq, crude, char = outs = self.outs
        tot_dw = feed.F_mass - feed.imass['Water']
        comps = (
            self.gas_composition,
            self.aqueous_composition,
            self.biocrude_composition,
            self.char_composition,
            )
        for out, comp in zip(outs, comps):
            out.empty()
            for k, v in comp.items():
                out.imass[k] = v

        dw_yields = self.dw_yields
        gas.F_mass = tot_dw * dw_yields['gas']
        aq.F_mass = tot_dw * dw_yields['aqueous']
        crude.F_mass = tot_dw * dw_yields['biocrude']
        char.F_mass = tot_dw * dw_yields['char']

        aq.imass['Water'] = feed.imass['Water'] - sum(i.imass['Water'] for i in (gas, crude, char))

        for i in outs:
            i.T = self.T
            i.P = self.P

        self._eff_at_temp.mix_from(outs)

        gas.phase = 'g'
        char.phase = 's'
        aq.phase = crude.phase = 'l'

        for attr, val in zip(('T', 'P'), (self.eff_T, self.eff_P)):
            if val:
                for i in self.outs: setattr(i, attr, val)
    
    def _design(self):
        hx = self.hx
        inf_hx = self.inf_hx
        inf_hx_in, inf_hx_out = inf_hx.ins[0], inf_hx.outs[0]
        inf_pre_hx, eff_pre_hx = hx.ins
        inf_after_hx, eff_after_hx = hx.outs
        inf_pre_hx.copy_like(self.ins[0])
        eff_pre_hx.copy_like(self._eff_at_temp)

        if self.internal_heat_exchanging:
            hx.phase0 = hx.phase1 = 'l'
            hx.T_lim1 = self.eff_T
            hx.simulate()
            for i in self.outs:
                i.T = eff_after_hx.T
        else:
            hx.empty()
            inf_after_hx.copy_like(inf_pre_hx)
            eff_after_hx.copy_like(eff_pre_hx)

        inf_hx_in.copy_like(inf_after_hx)
        inf_hx_out.copy_flow(inf_hx_in)
        # `copy_flow` only copies flow rates, not T/P
        inf_hx_out.T = self.T
        # No pressure drop is implied by this exchanger,
        # so carry inf_hx_in.P through (mirroring Hydroprocessing._design's
        # own inf_hx block below, which already does this).
        inf_hx_out.P = inf_hx_in.P
        inf_hx_in.vle(T=inf_hx_in.T, P=inf_hx_in.P)
        inf_hx_out.vle(T=inf_hx_out.T, P=inf_hx_out.P)
        inf_hx.simulate_as_auxiliary_exchanger(ins=inf_hx.ins, outs=inf_hx.outs)

        eff_hx = self.eff_hx
        eff_hx_in, eff_hx_out = eff_hx.ins[0], eff_hx.outs[0]
        eff_hx_in.copy_like(eff_after_hx)
        eff_hx_out.mix_from(self.outs)
        # eff_hx_in/eff_hx_out share composition (this exchanger only cools
        # or heats, it doesn't react), so simulate_as_auxiliary_exchanger's
        # own duty auto-computation (outlet.Hnet - inlet.Hnet) already equals
        # the correct sensible/latent duty.
        eff_hx.simulate_as_auxiliary_exchanger(ins=eff_hx.ins, outs=eff_hx.outs)

        Reactor._design(self)

        Design = self.design_results
        Design['Solid filter and separator weight'] = 0.2*Design['Weight']*Design['Number of reactors']
        if self.include_construction:
            self.construction[0].quantity += Design['Solid filter and separator weight']*_lb_to_kg

        if not self.use_decorated_cost:
            kodrum = self.kodrum
            kodrum.V = self.F_mass_out/_lb_to_kg/1225236*4230/_m3_to_gal
            kodrum.simulate()

    def _cost(self):
        self.baseline_purchase_costs.clear()
        if self.use_decorated_cost:
            Design = self.design_results
            Design.clear()
            ins0 = self.ins[0]
            Design['Wet mass flowrate'] = ins0.F_mass/_lb_to_kg
            self._decorated_cost()
        else:
            # `Reactor._cost` reads `design_results['Total volume']`, populated by
            # `_design`'s `Reactor._design(self)` call -- must not be cleared here.
            Reactor._cost(self)

    def _normalize_composition(self, dct):
        total = sum(dct.values())
        if total <= 0: raise ValueError(f'Sum of total yields/compositions should be positive, not {total}.')
        return {k: v/total for k, v in dct.items()}

    @property
    def dw_yields(self):
        return self._dw_yields
    @dw_yields.setter
    def dw_yields(self, comp_dct):
        self._dw_yields = self._normalize_composition(comp_dct)

    @property
    def gas_composition(self):
        return self._gas_composition
    @gas_composition.setter
    def gas_composition(self, comp_dct):
        self._gas_composition = self._normalize_composition(comp_dct)

    @property
    def aqueous_composition(self):
        return self._aqueous_composition
    @aqueous_composition.setter
    def aqueous_composition(self, comp_dct):
        self._aqueous_composition = self._normalize_composition(comp_dct)

    @property
    def biocrude_composition(self):
        return self._biocrude_composition
    @biocrude_composition.setter
    def biocrude_composition(self, comp_dct):
        self._biocrude_composition = self._normalize_composition(comp_dct)

    @property
    def char_composition(self):
        return self._char_composition
    @char_composition.setter
    def char_composition(self, comp_dct):
        self._char_composition = self._normalize_composition(comp_dct)

    @property
    def biocrude_HHV(self):
        '''[float] Higher heating value of the biocrude, MJ/kg.'''
        crude = self.outs[2]
        return crude.HHV/crude.F_mass/1e3

    @property
    def energy_recovery(self):
        '''[float] Fraction of the feedstock's HHV recovered in the biocrude.'''
        feed = self.ins[0]
        crude = self.outs[2]
        return crude.HHV/feed.HHV