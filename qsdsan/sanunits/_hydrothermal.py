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
from . import Reactor, HXutility, Pump

__all__ = (
    'CatalyticHydrothermalGasification',
    'HydrothermalLiquefaction',
    )

_lb_to_kg = auom('lb').conversion_factor('kg')
_m3_to_gal = auom('m3').conversion_factor('gallon')


# %%

# =============================================================================
# CHG
# =============================================================================

# hydrocyclone
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
    .. [1] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    .. [2] Davis, R. E.; Grundl, N. J.; Tao, L.; Biddy, M. J.; Tan, E. C.;
        Beckham, G. T.; Humbird, D.; Thompson, D. N.; Roni, M. S. Process
        Design and Economics for the Conversion of Lignocellulosic Biomass
        to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case
        Update; Biochemical Deconstruction and Conversion of Biomass to Fuels
        and Products via Integrated Biorefinery Pathways; NREL/TP--5100-71949,
        1483234; 2018; p NREL/TP--5100-71949, 1483234.
        https://doi.org/10.2172/1483234.
    .. [3] Elliott, D. C.; Neuenschwander, G. G.; Hart, T. R.; Rotness, L. J.;
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
                  length_to_diameter=2, N=6, V=None, auxiliary=False,
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
        
        for aux_unit in self.auxiliary_units:
            purchase_costs = aux_unit.baseline_purchase_costs
            installed_costs = aux_unit.installed_costs
            for item in purchase_costs.keys():
                purchase_costs[item] *= self.CAPEX_factor
                installed_costs[item] *= self.CAPEX_factor
                

# %%

# =============================================================================
# KOdrum
# =============================================================================

class KnockOutDrum(Reactor):
    '''
    Konckout drum is an auxiliary unit for :class:`HydrothermalLiquefaction`.
    
    References
    ----------
    .. [1] Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels via
        Liquefaction - Hydrothermal Liquefaction Reactor Design: April 5, 2013;
        NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462, 1111191.
        https://doi.org/10.2172/1111191.
        
    See Also
    --------
    :class:`qsdsan.sanunits.HydrothermalLiquefaction`
    '''
    _N_ins = 3
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream',
                 P=3049.7*6894.76, tau=0, V_wf=0,
                 length_to_diameter=2, N=4, V=None,
                 auxiliary=True,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
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
        pass

# =============================================================================
# HTL (ignore three phase separator for now, ask Yalin)
# =============================================================================

# separator
@cost(basis='Treatment capacity', ID='Solids filter oil/water separator', units='lb/h',
      cost=3945523, S=1219765,
      CE=CEPCI_by_year[2011], n=0.68, BM=1.9)
class HydrothermalLiquefaction(Reactor):
    '''
    HTL converts dewatered sludge to biocrude, aqueous, off-gas, and biochar
    under elevated temperature (350°C) and pressure. The products percentage
    (wt%) can be evaluated using revised MCA model (Li et al., 2017,
    Leow et al., 2018) with known sludge composition (protein%, lipid%,
    and carbohydrate%, all afdw%).
                                                      
    Notice that for HTL we just calculate each phases' total mass (except gas)
    and calculate C, N, and P amount in each phase as properties. We don't
    specify components for oil/char since we want to use MCA model to calculate
    C and N amount and it is not necessary to calculate every possible
    components since they will be treated in HT/AcidEx anyway. We also don't
    specify components for aqueous since we want to calculate aqueous C, N, and
    P based on mass balance closure. But later for CHG, HT, and HC, we specify
    each components (except aqueous phase) for the application of flash,
    distillation column, and CHP units.
    
    Parameters
    ----------
    ins : Iterable(stream)
        dewatered_sludge.
    outs : Iterable(stream)
        biochar, HTLaqueous, biocrude, offgas.
    lipid_2_biocrude: float
        Lipid to biocrude factor.
    protein_2_biocrude: float
        Protein to biocrude factor.
    carbo_2_biocrude: float
        Carbohydrate to biocrude factor.
    protein_2_gas: float
        Protein to gas factor.
    carbo_2_gas: float
        Carbohydrate to gas factor.
    biocrude_C_slope: float
        Biocrude carbon content slope.
    biocrude_C_intercept: float
        Biocrude carbon content intercept.
    biocrude_N_slope: float
        Biocrude nitrogen content slope.
    biocrude_H_slope: float
        Biocrude hydrogen content slope.
    biocrude_H_intercept: float
        Biocrude hydrogen content intercept.
    HTLaqueous_C_slope: float
        HTLaqueous carbon content slope.
    TOC_TC: float   
        HTL TOC/TC.
    biochar_C_slope: float
        Biochar carbon content slope.
    biocrude_moisture_content: float
        Biocrude moisture content.
    biochar_P_recovery_ratio: float
        Biochar phosphorus to total phosphorus ratio.
    gas_composition: dict
        HTL offgas compositions.
    biochar_pre: float
        Biochar pressure, [Pa].
    HTLaqueous_pre: float
        HTL aqueous phase pressure, [Pa].
    biocrude_pre: float
        Biocrude pressure, [Pa].
    offgas_pre: float
        Offgas pressure, [Pa].
    eff_T: float
        HTL effluent temperature, [K].
    CAPEX_factor: float
        Factor used to adjust CAPEX.
        
    References
    ----------
    .. [1] Leow, S.; Witter, J. R.; Vardon, D. R.; Sharma, B. K.;
        Guest, J. S.; Strathmann, T. J. Prediction of Microalgae Hydrothermal
        Liquefaction Products from Feedstock Biochemical Composition.
        Green Chem. 2015, 17 (6), 3584–3599. https://doi.org/10.1039/C5GC00574D.
    .. [2] Li, Y.; Leow, S.; Fedders, A. C.; Sharma, B. K.; Guest, J. S.;
        Strathmann, T. J. Quantitative Multiphase Model for Hydrothermal
        Liquefaction of Algal Biomass. Green Chem. 2017, 19 (4), 1163–1174.
        https://doi.org/10.1039/C6GC03294J.
    .. [3] Li, Y.; Tarpeh, W. A.; Nelson, K. L.; Strathmann, T. J.
        Quantitative Evaluation of an Integrated System for Valorization of
        Wastewater Algae as Bio-Oil, Fuel Gas, and Fertilizer Products.
        Environ. Sci. Technol. 2018, 52 (21), 12717–12727.
        https://doi.org/10.1021/acs.est.8b04035.
    .. [4] Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
        Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
        Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
        Process Design and Economics for the Conversion of Algal Biomass to
        Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
        PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.
    .. [5] Matayeva, A.; Rasmussen, S. R.; Biller, P. Distribution of Nutrients and
        Phosphorus Recovery in Hydrothermal Liquefaction of Waste Streams.
        BiomassBioenergy 2022, 156, 106323.
        https://doi.org/10.1016/j.biombioe.2021.106323.
    .. [6] Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels
        via Liquefaction - Hydrothermal Liquefaction Reactor Design:
        April 5, 2013; NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462,
        1111191. https://doi.org/10.2172/1111191.
    '''
    _N_ins = 1
    _N_outs = 4
    _units= {'Treatment capacity': 'lb/h',
             'Solid filter and separator weight': 'lb'}
    
    auxiliary_unit_names=('heat_exchanger','kodrum')

    _F_BM_default = {**Reactor._F_BM_default,
                     'Heat exchanger': 3.17}

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream',
                 lipid_2_biocrude=0.846, # [1]
                 protein_2_biocrude=0.445, # [1]
                 carbo_2_biocrude=0.205, # [1]
                 protein_2_gas=0.074, # [1]
                 carbo_2_gas=0.418, # [1]
                 biocrude_C_slope=-8.37, # [2]
                 biocrude_C_intercept=68.55, # [2]
                 biocrude_N_slope=0.133, # [2]
                 biocrude_H_slope=-2.61, # [2]
                 biocrude_H_intercept=8.20, # [2]
                 HTLaqueous_C_slope=478, # [2]
                 TOC_TC=0.764, # [3]
                 biochar_C_slope=1.75, # [2]
                 biocrude_moisture_content=0.063, # [4]
                 biochar_P_recovery_ratio=0.86, # [5]
                 gas_composition={'CH4':0.050, 'C2H6':0.032,
                                  'CO2':0.918}, # [4]
                 biochar_pre=3029.7*6894.76, # [4]
                 HTLaqueous_pre=30*6894.76, # [4]
                 biocrude_pre=30*6894.76, # [4]
                 offgas_pre=30*6894.76, # [4]
                 eff_T=60+273.15, # [4]
                 P=None, tau=15/60, V_wf=0.3,
                 length_to_diameter=2, N=4, V=None, auxiliary=False,
                 mixing_intensity=None, kW_per_m3=0,
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical',
                 CAPEX_factor=1):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.lipid_2_biocrude = lipid_2_biocrude
        self.protein_2_biocrude = protein_2_biocrude
        self.carbo_2_biocrude = carbo_2_biocrude
        self.protein_2_gas = protein_2_gas
        self.carbo_2_gas = carbo_2_gas
        self.biocrude_C_slope = biocrude_C_slope
        self.biocrude_C_intercept = biocrude_C_intercept
        self.biocrude_N_slope = biocrude_N_slope
        self.biocrude_H_slope = biocrude_H_slope
        self.biocrude_H_intercept = biocrude_H_intercept
        self.HTLaqueous_C_slope = HTLaqueous_C_slope
        self.TOC_TC = TOC_TC
        self.biochar_C_slope = biochar_C_slope
        self.biocrude_moisture_content = biocrude_moisture_content
        self.biochar_P_recovery_ratio = biochar_P_recovery_ratio
        self.gas_composition = gas_composition
        self.biochar_pre = biochar_pre
        self.HTLaqueous_pre = HTLaqueous_pre
        self.biocrude_pre = biocrude_pre
        self.offgas_pre = offgas_pre
        hx_in = Stream(f'{ID}_hx_in')
        hx_out = Stream(f'{ID}_hx_out')
        self.heat_exchanger = HXutility(ID=f'.{ID}_hx', ins=hx_in, outs=hx_out, T=eff_T, rigorous=True)
        self.kodrum = KnockOutDrum(ID=f'.{ID}_KOdrum')
        self.P = P
        self.tau = tau
        self.V_wf = V_wf
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

    
    def _run(self):
        
        dewatered_sludge = self.ins[0]
        biochar, HTLaqueous, biocrude, offgas = self.outs
        
        self.WWTP = self.ins[0]._source.ins[0]._source.ins[0].\
                         _source.ins[0]._source
        
        dewatered_sludge_afdw = dewatered_sludge.imass['Sludge_lipid'] +\
                                dewatered_sludge.imass['Sludge_protein'] +\
                                dewatered_sludge.imass['Sludge_carbo']
        # just use afdw in revised MCA model, other places use dw
        
        afdw_lipid_ratio = self.WWTP.sludge_afdw_lipid
        afdw_protein_ratio = self.WWTP.sludge_afdw_protein
        afdw_carbo_ratio = self.WWTP.sludge_afdw_carbo

        # the following calculations are based on revised MCA model
        biochar.imass['Biochar'] = 0.377*afdw_carbo_ratio*dewatered_sludge_afdw
        
        HTLaqueous.imass['HTLaqueous'] = (0.481*afdw_protein_ratio +\
                                          0.154*afdw_lipid_ratio)*\
                                          dewatered_sludge_afdw
        # HTLaqueous is TDS in aqueous phase
        # 0.377, 0.481, and 0.154 don't have uncertainties because they are calculated values
         
        gas_mass = (self.protein_2_gas*afdw_protein_ratio + self.carbo_2_gas*afdw_carbo_ratio)*\
                       dewatered_sludge_afdw
                       
        for name, ratio in self.gas_composition.items():
            offgas.imass[name] = gas_mass*ratio
            
        biocrude.imass['Biocrude'] = (self.protein_2_biocrude*afdw_protein_ratio +\
                                      self.lipid_2_biocrude*afdw_lipid_ratio +\
                                      self.carbo_2_biocrude*afdw_carbo_ratio)*\
                                      dewatered_sludge_afdw
        biocrude.imass['H2O'] = biocrude.imass['Biocrude']/(1 -\
                                self.biocrude_moisture_content) -\
                                biocrude.imass['Biocrude']
                                
        HTLaqueous.imass['H2O'] = dewatered_sludge.F_mass - biochar.F_mass -\
                                  biocrude.F_mass - gas_mass - HTLaqueous.imass['HTLaqueous']
        # assume ash (all soluble based on Jones) goes to water
        
        biochar.phase = 's'
        offgas.phase = 'g'
        HTLaqueous.phase = biocrude.phase = 'l'
        
        biochar.P = self.biochar_pre
        HTLaqueous.P = self.HTLaqueous_pre
        biocrude.P = self.biocrude_pre
        offgas.P = self.offgas_pre
        
        for stream in self.outs : stream.T = self.heat_exchanger.T

    @property
    def biocrude_C_ratio(self):
        return (self.WWTP.AOSc*self.biocrude_C_slope + self.biocrude_C_intercept)/100 # [2]
    
    @property
    def biocrude_H_ratio(self):
        return (self.WWTP.AOSc*self.biocrude_H_slope + self.biocrude_H_intercept)/100 # [2]

    @property
    def biocrude_N_ratio(self):
        return self.biocrude_N_slope*self.WWTP.sludge_dw_protein # [2]
    
    @property
    def biocrude_C(self):
        return min(self.outs[2].F_mass*self.biocrude_C_ratio, self.WWTP.sludge_C)


    @property
    def HTLaqueous_C(self):
        return min(self.outs[1].F_vol*1000*self.HTLaqueous_C_slope*\
                   self.WWTP.sludge_dw_protein*100/1000000/self.TOC_TC,
                   self.WWTP.sludge_C - self.biocrude_C)

    @property
    def biocrude_H(self):
        return self.outs[2].F_mass*self.biocrude_H_ratio

    @property
    def biocrude_N(self):
        return min(self.outs[2].F_mass*self.biocrude_N_ratio, self.WWTP.sludge_N)
    
    @property
    def biocrude_HHV(self):
        return 30.74 - 8.52*self.WWTP.AOSc +\
               0.024*self.WWTP.sludge_dw_protein # [2]
               
    @property
    def energy_recovery(self):
        return self.biocrude_HHV*self.outs[2].imass['Biocrude']/\
               (self.WWTP.outs[0].F_mass -\
               self.WWTP.outs[0].imass['H2O'])/self.WWTP.sludge_HHV # [2]
        
    @property
    def offgas_C(self):
        carbon = sum(self.outs[3].imass[self.gas_composition]*
                     [cmp.i_C for cmp in self.components[self.gas_composition]])
        return min(carbon, self.WWTP.sludge_C - self.biocrude_C - self.HTLaqueous_C)
        
    @property
    def biochar_C_ratio(self):
        return min(self.biochar_C_slope*self.WWTP.sludge_dw_carbo, 0.65) # [2]

    @property
    def biochar_C(self):
        return min(self.outs[0].F_mass*self.biochar_C_ratio, self.WWTP.sludge_C -\
                   self.biocrude_C - self.HTLaqueous_C - self.offgas_C)

    @property
    def biochar_P(self):
        return min(self.WWTP.sludge_P*self.biochar_P_recovery_ratio, self.outs[0].F_mass)

    @property
    def HTLaqueous_N(self):
        return self.WWTP.sludge_N - self.biocrude_N
        
    @property
    def HTLaqueous_P(self):
        return self.WWTP.sludge_P*(1 - self.biochar_P_recovery_ratio)

    def _design(self):
        
        Design = self.design_results
        Design['Treatment capacity'] = self.ins[0].F_mass/_lb_to_kg
        
        hx = self.heat_exchanger
        hx_ins0, hx_outs0 = hx.ins[0], hx.outs[0]
        hx_ins0.mix_from((self.outs[1], self.outs[2], self.outs[3]))
        hx_outs0.copy_like(hx_ins0)
        hx_ins0.T = self.ins[0].T # temperature before/after HTL are similar
        hx_outs0.T = hx.T
        hx_ins0.P = hx_outs0.P = self.outs[0].P # cooling before depressurized, heating after pressurized
        # in other words, both heating and cooling are performed under relatively high pressure
        hx_ins0.vle(T=hx_ins0.T, P=hx_ins0.P)
        hx_outs0.vle(T=hx_outs0.T, P=hx_outs0.P)
        hx.simulate_as_auxiliary_exchanger(ins=hx.ins, outs=hx.outs)

        self.P = self.ins[0].P
        Reactor._design(self)
        Design['Solid filter and separator weight'] = 0.2*Design['Weight']*Design['Number of reactors'] # assume stainless steel
        # based on [6], case D design table, the purchase price of solid filter and separator to
        # the purchase price of HTL reactor is around 0.2, therefore, assume the weight of solid filter
        # and separator is 0.2*single HTL weight*number of HTL reactors
        self.construction[0].quantity += Design['Solid filter and separator weight']*_lb_to_kg
        
        self.kodrum.V = self.F_mass_out/_lb_to_kg/1225236*4230/_m3_to_gal
        # in [6], when knockout drum influent is 1225236 lb/hr, single knockout
        # drum volume is 4230 gal
        
        self.kodrum.simulate()
        
    def _cost(self):
        Reactor._cost(self)
        self._decorated_cost()
        
        purchase_costs = self.baseline_purchase_costs
        for item in purchase_costs.keys():
            purchase_costs[item] *= self.CAPEX_factor
        
        for aux_unit in self.auxiliary_units:
            purchase_costs = aux_unit.baseline_purchase_costs
            installed_costs = aux_unit.installed_costs
            for item in purchase_costs.keys():
                purchase_costs[item] *= self.CAPEX_factor
                installed_costs[item] *= self.CAPEX_factor