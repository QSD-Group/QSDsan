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
from . import Reactor, IsothermalCompressor, HXutility

__all__ = ('Hydrocracking', 'Hydrotreating')

_lb_to_kg = auom('lb').conversion_factor('kg')
_m3perh_to_mmscfd = 1/1177.17 # H2


# %%

# =============================================================================
# HC
# =============================================================================

class Hydrocracking(Reactor):
    '''
    Biocrude mixed with H2 are hydrotreated at elevated temperature (405°C)
    and pressure to produce upgraded biooil. Co-product includes fuel gas.
    
    Parameters
    ----------
    ins : Iterable(stream)
        heavy_oil, hydrogen, catalyst_in.
    outs : Iterable(stream)
        hc_out, catalyst_out.
    WHSV: float
        Weight Hourly Space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime: float
        CHG catalyst lifetime, [hr].
    hydrogen_P: float
        Hydrogen pressure, [Pa].
    hydrogen_rxned_to_heavy_oil: float
        Reacted H2 to heavy oil mass ratio.
    hydrogen_excess: float
        Actual hydrogen amount = hydrogen_rxned_to_biocrude*hydrogen_excess
    hydrocarbon_ratio: float
        Mass ratio of produced hydrocarbon to the sum of heavy oil and reacted H2.
    HCin_T: float
        HC influent temperature, [K].
    HCrxn_T: float
        HC effluent (after reaction) temperature, [K].
    HC_composition: dict
        HC effluent composition.
        
    References
    ----------
    .. [1] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    '''
    _N_ins = 3
    _N_outs = 2
    
    auxiliary_unit_names=('compressor','heat_exchanger',)
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17,
                     'Compressor': 1.1}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 WHSV=0.625, # wt./hr per wt. catalyst [1]
                 catalyst_lifetime=5*7920, # 5 years [1]
                 hydrogen_P=1039.7*6894.76,
                 hydrogen_rxned_to_heavy_oil=0.01125,
                 hydrogen_excess=5.556,
                 hydrocarbon_ratio=1, # 100 wt% of heavy oil and reacted H2
                 # nearly all input heavy oils and H2 will be converted to
                 # products [1]
                 # spreadsheet HC calculation
                 HCin_T=394+273.15,
                 HCrxn_T=451+273.15,
                 HC_composition={'CO2':0.03880, 'CH4':0.00630,
                                 'CYCHEX':0.03714, 'HEXANE':0.01111,
                                 'HEPTANE':0.11474, 'OCTANE':0.08125,
                                 'C9H20':0.09086, 'C10H22':0.11756,
                                 'C11H24':0.16846, 'C12H26':0.13198,
                                 'C13H28':0.09302, 'C14H30':0.04643,
                                 'C15H32':0.03250, 'C16H34':0.01923,
                                 'C17H36':0.00431, 'C18H38':0.00099,
                                 'C19H40':0.00497, 'C20H42':0.00033},
                 #combine C20H42 and PHYTANE as C20H42
                 # will not be a variable in uncertainty/sensitivity analysis
                 P=None, tau=5, void_fraciton=0.4, # Towler
                 length_to_diameter=2, N=None, V=None, auxiliary=False, mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1.5,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.hydrogen_P = hydrogen_P
        self.hydrogen_rxned_to_heavy_oil = hydrogen_rxned_to_heavy_oil
        self.hydrogen_excess = hydrogen_excess
        self.hydrocarbon_ratio = hydrocarbon_ratio
        self.HCin_T = HCin_T
        self.HCrxn_T = HCrxn_T
        self.HC_composition = HC_composition
        IC_in = Stream(f'{ID}_IC_in')
        IC_out = Stream(f'{ID}_IC_out')
        self.compressor = IsothermalCompressor(ID=f'.{ID}_IC', ins=IC_in,
                                               outs=IC_out, P=None)
        hx_H2_in = Stream(f'{ID}_hx_H2_in')
        hx_H2_out = Stream(f'{ID}_hx_H2_out')
        self.heat_exchanger_H2 = HXutility(ID=f'.{ID}_hx_H2', ins=hx_H2_in, outs=hx_H2_out)
        hx_oil_in = Stream(f'{ID}_hx_oil_in')
        hx_oil_out = Stream(f'{ID}_hx_oil_out')
        self.heat_exchanger_oil = HXutility(ID=f'.{ID}_hx_oil', ins=hx_oil_in, outs=hx_oil_out)
        self.P = P
        self.tau = tau
        self.void_fraciton = void_fraciton
        self.length_to_diameter = length_to_diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        
    def _run(self):
        
        heavy_oil, hydrogen, catalyst_in = self.ins
        hc_out, catalyst_out = self.outs
        
        catalyst_in.imass['HC_catalyst'] = heavy_oil.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        # catalysts amount is quite low compared to the main stream, therefore do not consider
        # heating/cooling of catalysts
        
        hydrogen.imass['H2'] = heavy_oil.F_mass*self.hydrogen_rxned_to_heavy_oil*self.hydrogen_excess
        hydrogen.phase = 'g'

        hydrocarbon_mass = heavy_oil.F_mass*(1 +\
                           self.hydrogen_rxned_to_heavy_oil)*\
                           self.hydrocarbon_ratio

        hc_out.phase = 'g'

        for name, ratio in self.HC_composition.items():
            hc_out.imass[name] = hydrocarbon_mass*ratio
        
        hc_out.imass['H2'] = heavy_oil.F_mass*self.hydrogen_rxned_to_heavy_oil*(self.hydrogen_excess - 1)
        
        hc_out.P = heavy_oil.P
        hc_out.T = self.HCrxn_T
        
        hc_out.vle(T=hc_out.T, P=hc_out.P)
        
        cmps = self.components
        C_in = 0
        total_num = len(list(cmps))
        for num in range(total_num):
            C_in += heavy_oil.imass[str(list(cmps)[num])]*list(cmps)[num].i_C
            
        C_out = self.hydrocarbon_C
        
        if C_out < 0.95*C_in or C_out > 1.05*C_out :
            raise Exception('carbon mass balance is out of +/- 5% for HC')
        # make sure that carbon mass balance is within +/- 5%. Otherwise, an
        # exception will be raised.
        
    @property
    def hydrocarbon_C(self):   
        return sum(self.outs[0].imass[self.HC_composition]*
                   [cmp.i_C for cmp in self.components[self.HC_composition]])

    def _design(self):
        IC = self.compressor
        IC_ins0, IC_outs0 = IC.ins[0], IC.outs[0]
        IC_ins0.copy_like(self.ins[1])
        IC_outs0.copy_like(self.ins[1])
        IC_outs0.P = IC.P = self.hydrogen_P
        IC_ins0.phase = IC_outs0.phase = 'g'
        IC.simulate()
        
        hx_H2 = self.heat_exchanger_H2
        hx_H2_ins0, hx_H2_outs0 = hx_H2.ins[0], hx_H2.outs[0]
        hx_H2_ins0.copy_like(self.ins[1])
        hx_H2_outs0.copy_like(hx_H2_ins0)
        hx_H2_ins0.phase = hx_H2_outs0.phase = 'g'
        hx_H2_outs0.T = self.HCin_T
        hx_H2_ins0.P = hx_H2_outs0.P = IC_outs0.P
        hx_H2.simulate_as_auxiliary_exchanger(ins=hx_H2.ins, outs=hx_H2.outs)

        hx_oil = self.heat_exchanger_oil
        hx_oil_ins0, hx_oil_outs0 = hx_oil.ins[0], hx_oil.outs[0]
        hx_oil_ins0.copy_like(self.ins[0])
        hx_oil_outs0.copy_like(hx_oil_ins0)
        hx_oil_outs0.T = self.HCin_T
        hx_oil_ins0.P = hx_oil_outs0.P = self.ins[0].P
        hx_oil.simulate_as_auxiliary_exchanger(ins=hx_oil.ins, outs=hx_oil.outs)
        
        self.P = min(IC_outs0.P, self.ins[0].P)
        
        V_H2 = self.ins[1].F_vol/self.hydrogen_excess*101325/self.hydrogen_P
        # just account for reacted H2
        V_biocrude = self.ins[0].F_vol
        self.V_wf = self.void_fraciton*V_biocrude/(V_biocrude + V_H2)
        Reactor._design(self)

        
# %%

# =============================================================================
# HT
# =============================================================================

@cost(basis='Hydrogen_PSA', ID='PSA', units='mmscfd',
      cost=1750000, S=10,
      CE=CEPCI_by_year[2004], n=0.8, BM=2.47)
class Hydrotreating(Reactor):
    '''
    Biocrude mixed with H2 are hydrotreated at elevated temperature (405°C)
    and pressure to produce upgraded biooil. Co-product includes fuel gas.
    A pressure swing adsorption (PSA) process can be optionally included
    for H2 recovery.
    
    Parameters
    ----------
    ins : Iterable(stream)
        biocrude, hydrogen, catalyst_in.
    outs : Iterable(stream)
        ht_out, catalyst_out = self.outs.
    WHSV: float
        Weight Hourly Space velocity, [kg feed/hr/kg catalyst].
    catalyst_lifetime: float
        CHG catalyst lifetime, [hr].
    hydrogen_P: float
        Hydrogen pressure, [Pa].
    hydrogen_rxned_to_biocrude: float
        Reacted H2 to biocrude mass ratio.
    hydrogen_excess: float
        Actual hydrogen amount = hydrogen_rxned_to_biocrude*hydrogen_excess
    hydrocarbon_ratio: float
        Mass ratio of produced hydrocarbon to the sum of biocrude and reacted H2.
    HTin_T: float
        HT influent temperature, [K].
    HTrxn_T: float
        HT effluent (after reaction) temperature, [K].
    HT_composition: dict
        HT effluent composition.
    CAPEX_factor: float
        Factor used to adjust CAPEX.
    include_PSA : bool
        Whether to include pressure swing adsorption for H2 recovery.
        
    References
    ----------
    .. [1] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    .. [2] Towler, G.; Sinnott, R. Chapter 14 - Design of Pressure Vessels.
        In Chemical Engineering Design (Second Edition); Towler, G., Sinnott, R.,
        Eds.; Butterworth-Heinemann: Boston, 2013; pp 563–629.
        https://doi.org/10.1016/B978-0-08-096659-5.00014-6.
    '''
    _N_ins = 3
    _N_outs = 2
    auxiliary_unit_names=('compressor','heat_exchanger',)
    
    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17,
                     'Compressor': 1.1}
    
    _units = {'Hydrogen': 'mmscfd',
              'Hydrogen_PSA': 'mmscfd'}
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 WHSV=0.625, # wt./hr per wt. catalyst [1]
                 catalyst_lifetime=2*7920, # 2 years [1]
                 hydrogen_P=1530*6894.76,
                 hydrogen_rxned_to_biocrude=0.046,
                 hydrogen_excess=3,
                 hydrocarbon_ratio=0.875, # 87.5 wt% of biocrude and reacted H2 [1]
                 # spreadsheet HT calculation
                 HTin_T=174+273.15,
                 HTrxn_T=402+273.15, # [1]
                 HT_composition={'CH4':0.02280, 'C2H6':0.02923,
                                 'C3H8':0.01650, 'C4H10':0.00870,
                                 'TWOMBUTAN':0.00408, 'NPENTAN':0.00678,
                                 'TWOMPENTA':0.00408, 'HEXANE':0.00401,
                                 'TWOMHEXAN':0.00408, 'HEPTANE':0.00401,
                                 'CC6METH':0.01020, 'PIPERDIN':0.00408,
                                 'TOLUENE':0.01013, 'THREEMHEPTA':0.01020,
                                 'OCTANE':0.01013, 'ETHCYC6':0.00408,
                                 'ETHYLBEN':0.02040, 'OXYLENE':0.01020,
                                 'C9H20':0.00408, 'PROCYC6':0.00408,
                                 'C3BENZ':0.01020, 'FOURMONAN':0,
                                 'C10H22':0.00240, 'C4BENZ':0.01223,
                                 # C10H22 was originally 0.00203, but it is not
                                 # good for distillation column, the excess amount
                                 # is substracted from HEXANE, HEPTANE, TOLUENE,
                                 # OCTANE, and C9H20, which were originally 0.00408,
                                 # 0.00408, 0.01020, 0.01020, and 0.00408
                                 'C11H24':0.02040, 'C10H12':0.02040,
                                 'C12H26':0.02040, 'OTTFNA':0.01020,
                                 'C6BENZ':0.02040, 'OTTFSN':0.02040,
                                 'C7BENZ':0.02040, 'C8BENZ':0.02040,
                                 'C10H16O4':0.01837, 'C15H32':0.06120,
                                 'C16H34':0.18360, 'C17H36':0.08160, 
                                 'C18H38':0.04080, 'C19H40':0.04080,
                                 'C20H42':0.10200, 'C21H44':0.04080,
                                 'TRICOSANE':0.04080, 'C24H38O4':0.00817,
                                 'C26H42O4':0.01020, 'C30H62':0.00203}, # [1]
                 # spreadsheet HT calculation
                 # will not be a variable in uncertainty/sensitivity analysis
                 P=None, tau=0.5, void_fraciton=0.4, # [2]
                 length_to_diameter=2, N=None, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 CAPEX_factor=1,
                 include_PSA=False,
                 PSA_pre=717.4*6894.76,
                 PSA_efficiency=0.9,):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.WHSV = WHSV
        self.catalyst_lifetime = catalyst_lifetime
        self.hydrogen_P = hydrogen_P
        self.hydrogen_rxned_to_biocrude = hydrogen_rxned_to_biocrude
        self.hydrogen_excess = hydrogen_excess
        self.hydrocarbon_ratio = hydrocarbon_ratio
        self.HTin_T = HTin_T
        self.HTrxn_T = HTrxn_T
        self.HT_composition = HT_composition
        IC_in = Stream(f'{ID}_IC_in')
        IC_out = Stream(f'{ID}_IC_out')
        self.compressor = IsothermalCompressor(ID=f'.{ID}_IC', ins=IC_in,
                                               outs=IC_out, P=None)
        hx_H2_in = Stream(f'{ID}_hx_H2_in')
        hx_H2_out = Stream(f'{ID}_hx_H2_out')
        self.heat_exchanger_H2 = HXutility(ID=f'.{ID}_hx_H2', ins=hx_H2_in, outs=hx_H2_out)
        hx_oil_in = Stream(f'{ID}_hx_oil_in')
        hx_oil_out = Stream(f'{ID}_hx_oil_out')
        self.heat_exchanger_oil = HXutility(ID=f'.{ID}_hx_oil', ins=hx_oil_in, outs=hx_oil_out)
        self.P = P
        self.tau = tau
        self.void_fraciton = void_fraciton
        self.length_to_diameter = length_to_diameter
        self.N = N
        self.V = V
        self.auxiliary = auxiliary
        self.mixing_intensity = mixing_intensity
        self.kW_per_m3 = kW_per_m3
        self.wall_thickness_factor = wall_thickness_factor
        self.vessel_material = vessel_material
        self.vessel_type = vessel_type
        self.CAPEX_factor = CAPEX_factor
        self.include_PSA = include_PSA
        self.PSA_pre = PSA_pre
        self.PSA_efficiency = PSA_efficiency
        
    def _run(self):
        
        biocrude, hydrogen, catalyst_in = self.ins
        ht_out, catalyst_out = self.outs
        
        self.HTL = self.ins[0]._source.ins[0]._source
        
        HT_composition = self.HT_composition
        if self.HTL.biocrude_N == 0:
            remove = HT_composition['PIPERDIN']
            for chemical in self.HT_composition.keys():  
                HT_composition[chemical] /= (1-remove)
            HT_composition['PIPERDIN'] = 0
        
        catalyst_in.imass['HT_catalyst'] = biocrude.F_mass/self.WHSV/self.catalyst_lifetime
        catalyst_in.phase = 's'
        catalyst_out.copy_like(catalyst_in)
        # catalysts amount is quite low compared to the main stream, therefore do not consider
        # heating/cooling of catalysts
        
        hydrogen_excess = self.hydrogen_excess
        H2_rxned =  biocrude.imass['Biocrude']*self.hydrogen_rxned_to_biocrude
        recovered_frac =  (hydrogen_excess - 1)*self.PSA_efficiency*float(self.include_PSA)
        hydrogen.imass['H2'] = H2_rxned*(hydrogen_excess - recovered_frac)
        hydrogen.phase = 'g'

        hydrocarbon_mass = biocrude.imass['Biocrude']*\
                           (1 + self.hydrogen_rxned_to_biocrude)*\
                           self.hydrocarbon_ratio
                           
        ht_out.phase = 'g'
                           
        for name, ratio in self.HT_composition.items():
            ht_out.imass[name] = hydrocarbon_mass*ratio
            
        ht_out.imass['H2'] = H2_rxned*(self.hydrogen_excess - recovered_frac - 1)
        
        ht_out.imass['H2O'] = biocrude.F_mass + hydrogen.F_mass -\
                              hydrocarbon_mass - ht_out.imass['H2']
        # use water to represent HT aqueous phase,
        # C and N can be calculated base on MB closure.
        
        ht_out.P = biocrude.P
        
        ht_out.T = self.HTrxn_T
        
        ht_out.vle(T=ht_out.T, P=ht_out.P)
        
        if self.HTaqueous_C < -0.1*self.HTL.WWTP.sludge_C:
            raise Exception('carbon mass balance is out of +/- 10% for the whole system')
        # allow +/- 10% out of mass balance
        # should be no C in the aqueous phase, the calculation here is just for MB
        
        if self.HTaqueous_N < -0.1*self.HTL.WWTP.sludge_N:
            raise Exception('nitrogen mass balance is out of +/- 10% for the whole system')
        # allow +/- 10% out of mass balance

        # possibility exist that more carbon is in biooil and gas than in
        # biocrude because we use the biooil/gas compositions to calculate
        # carbon. In this case, the C in HT aqueous phase will be negative.
        # It's OK if the mass balance is within +/- 10% of total carbon in 
        # sludge. Otherwise, an exception will be raised.
        
    @property
    def hydrocarbon_C(self):   
        return sum(self.outs[0].imass[self.HT_composition]*
                   [cmp.i_C for cmp in self.components[self.HT_composition]])

    @property
    def hydrocarbon_N(self):
        return sum(self.outs[0].imass[self.HT_composition]*
                   [cmp.i_N for cmp in self.components[self.HT_composition]])

    @property
    def HTaqueous_C(self):
        return self.HTL.biocrude_C - self.hydrocarbon_C
    # should be no C in the aqueous phase, the calculation here is just for MB

    @property
    def HTaqueous_N(self):
        return self.HTL.biocrude_N - self.hydrocarbon_N
    
    @property
    def PSA_efficiency (self):
        return self._PSA_efficiency 
    @PSA_efficiency.setter
    def PSA_efficiency(self, i):
        if i > 1: raise Exception('PSA efficiency cannot be larger than 1.')
        self._PSA_efficiency  = i

    def _design(self):
        
        IC = self.compressor
        IC_ins0, IC_outs0 = IC.ins[0], IC.outs[0]
        IC_ins0.copy_like(self.ins[1])
        IC_outs0.copy_like(self.ins[1])
        IC_outs0.P = IC.P = self.hydrogen_P
        IC_ins0.phase = IC_outs0.phase = 'g'
        IC.simulate()

        hx_H2 = self.heat_exchanger_H2
        hx_H2_ins0, hx_H2_outs0 = hx_H2.ins[0], hx_H2.outs[0]
        hx_H2_ins0.copy_like(self.ins[1])
        hx_H2_outs0.copy_like(hx_H2_ins0)
        hx_H2_ins0.phase = hx_H2_outs0.phase = 'g'
        hx_H2_outs0.T = self.HTin_T
        hx_H2_ins0.P = hx_H2_outs0.P = IC_outs0.P
        hx_H2.simulate_as_auxiliary_exchanger(ins=hx_H2.ins, outs=hx_H2.outs)
        
        hx_oil = self.heat_exchanger_oil
        hx_oil_ins0, hx_oil_outs0 = hx_oil.ins[0], hx_oil.outs[0]
        hx_oil_ins0.copy_like(self.ins[0])
        hx_oil_outs0.copy_like(hx_oil_ins0)
        hx_oil_outs0.T = self.HTin_T
        hx_oil_ins0.P = hx_oil_outs0.P = self.ins[0].P
        hx_oil.simulate_as_auxiliary_exchanger(ins=hx_oil.ins, outs=hx_oil.outs)
        
        self.P = min(IC_outs0.P, self.ins[0].P)
        
        V_H2 = self.ins[1].F_vol/self.hydrogen_excess*101325/self.hydrogen_P
        # just account for reacted H2
        V_biocrude = self.ins[0].F_vol
        self.V_wf = self.void_fraciton*V_biocrude/(V_biocrude + V_H2)
        Reactor._design(self)
        
        Design = self.design_results
        factor = float(self.include_PSA)
        Design['Hydrogen_PSA'] = self.ins[1].F_vol*_m3perh_to_mmscfd*101325/self.PSA_pre*factor
        Design['PSA'] = 0.5*Design['Weight']*Design['Number of reactors']*factor # assume stainless steel
        # based on [1], page 54, the purchase price of PSA to the purchase price of
        # HT reactor (excluding vessels and columns) is around 0.5,
        # therefore, assume the weight of PSA is 0.5*single HT weight*number of HT reactors
        self.construction[0].quantity += Design['PSA']*_lb_to_kg*factor
        
    
    def _cost(self):
        Reactor._cost(self)
        
        purchase_costs = self.baseline_purchase_costs
        CAPEX_factor = self.CAPEX_factor
        if self.include_PSA: self._decorated_cost()
        else: purchase_costs['Hydrogen_PSA'] = 0
        
        for item in purchase_costs.keys():
            purchase_costs[item] *= CAPEX_factor
        
        for aux_unit in self.auxiliary_units:
            purchase_costs = aux_unit.baseline_purchase_costs
            installed_costs = aux_unit.installed_costs
            for item in purchase_costs.keys():
                purchase_costs[item] *= CAPEX_factor
                installed_costs[item] *= CAPEX_factor
