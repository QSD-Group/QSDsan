#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

    This module contains unit operations for the G2RT non-sewered toilet

    and resource-recovery systems (EXPOsan's ``g2rt`` module).

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''
from warnings import warn
from math import ceil
import numpy as np
import CoolProp.CoolProp as CP
from ..bst import Mixer, IsothermalCompressor
from ... import SanUnit, Construction, WasteStream
from ._excretion import Excretion
from ._non_reactive import Copier
from ...utils import ospath, data_path, load_data, price_ratio

g2rt_su_data_path = ospath.join(data_path, 'sanunit_data/g2rt')

__all__ = ('FWMixer',
           'G2RTBeltSeparation',
           'G2RTControls',
           'G2RTExcretion',
           'G2RThomogenizer',
           'G2RTHousing',
           'G2RTLiquidsTank',
           'G2RTReverseOsmosis',
           'G2RTSolidsSeparation',
           'G2RTSolidsTank',
           'G2RTUltrafiltration',
           'mSCWOConcentratorModule',
           'mSCWOGasModule',
           'mSCWOReactorModule',
           'UFMixer',
           'VolumeReductionCombustor',
           'VolumeReductionFilterPress',
           'VRConcentrator',
           'VRdryingtunnel',
           'VRpasteurization',
           )

#%%
mscwo_gas_module_path = ospath.join(g2rt_su_data_path, '_mscwo_gas_handling_module.csv')
@price_ratio()
class mSCWOGasModule(IsothermalCompressor):
    '''
    Gas handling unit that consisted of compressors and an injection pressure vessel
    
    The following components should be included in system thermo object for simulation:
    O2, N2, CO2  #TODO: update

    The following impact items should be pre-constructed for life cycle assessment:
    Compressor_4kW, Compressor_300kW, StainlessSteel

    Parameters
    ----------
    ins : Iterable(stream)
        Air in
    outs : Iterable(stream)
        Compressed air out
    P : float
        Outlet pressure [Pa].
    eta : float
        Isothermal efficiency.
    vle : bool
        Whether to perform phase equilibrium calculations on
        the outflow. If False, the outlet will be assumed to be the same
        phase as the inlet.
    type: str
        Type of compressor : blower/centrifugal/reciprocating. If None, the type
        will be determined automatically.

    References
    ----------
    [1] YEE et al. Water oxidation non-sewered single unit toilet system. 
    https://patentimages.storage.googleapis.com/57/6a/81/72a168a92be44c/WO2023288331A1.pdf
    '''
    _N_ins = 1
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 P=1.55E7, eta =0.7,vle=False, compressor_type='Reciprocating', 
                 driver=None, material=None, driver_efficiency=None,include_construction= True,
                 user_scale_up=1, extreme= False, capital_scale_func=None, **kwargs):
        IsothermalCompressor.__init__(self, ID=ID, ins=ins, outs=outs,
                                      P = P, eta = eta, vle = vle,
                                      compressor_type = compressor_type, driver= driver,
                                      material = material, driver_efficiency = driver_efficiency
                                      )
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.include_construction = include_construction
        self.user_scale_up = user_scale_up
        self.extreme = extreme
        self.capital_scale_func = capital_scale_func
        data = load_data(path=mscwo_gas_module_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        
        self.ins.P = 101325
        self.outs.P = self.P
    
    def _init_lca(self):
        self.include_construction = True
        self.construction = [
            Construction('compressor_4kW', linked_unit=self, 
                         item='Compressor_4kW', 
                         quantity_unit='ea'),
            Construction('compressor_300kW', linked_unit=self, 
                         item='Compressor_300kW', 
                         quantity_unit='ea'),
            Construction("stainless_steel", linked_unit=self,
                         item = "StainlessSteel", 
                         quantity_unit= "kg"),
            ]
    
    def _run(self):
        air_in, = self.ins
        air_in.P = 101325
        air_out, = self.outs
        air_out.copy_like(air_in)
        air_out.P = self.P
        air_out.T = air_in.T
        if self.vle is True: air_out.vle(T=air_out.T, P=air_out.P)
        self.ideal_power, self.ideal_duty = self._calculate_ideal_power_and_duty()
    
    
    def _design(self):
        self._init_lca()
        air_in, = self.ins
        air_out, = self.outs
        air_in.imass['N2'] = air_out.imass['N2']
        air_in.imass['O2'] = air_out.imass['O2']
        air_out.P = self.P
        air_out.T = air_in.T
        super()._design()
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[2].quantity = self.stainless_steel_weight
        self.add_construction(add_cost = False)
        
    def _cost(self):
        super()._design()
        D = self.design_results
        C = self.baseline_purchase_costs
        C['Compressor'] = self.compressor_cost *(self.user_scale_up ** 0.6)
        C['Valves'] = (self.injection_valve_cost +
                       self.dosing_valve_cost * self.dosing_valve_quantity
                       )
        C['Vessel_cost'] = self.vessel_cost
        C['Misc.parts'] = self.miscellaneous_cost_ratio*(C['Valves'] +
                               C['Vessel_cost'])
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio

        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = ((self.capital_scale_func(total_equipment) if self.capital_scale_func else total_equipment)*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
        if self.extreme:
            self.power_utility(self.dosing_valve_power_demand * self.dosing_valve_daily_operation/24*self.user_scale_up+
                           D['Ideal power']/self.eta/self.extreme_compressor_efficiency
                           ) # kWh/hr
        else:
            self.power_utility(self.dosing_valve_power_demand * self.dosing_valve_daily_operation/24*self.user_scale_up+
                           D['Ideal power']/self.eta/self.compressor_efficiency
                           ) # kWh/hr

    def _calc_maintenance_labor_cost(self): #USD/hr
        if self.extreme:
            maintenance_labor_cost= (self.gas_module_maintenance * self.wages * self.user_scale_up)*0.1
        else:
            maintenance_labor_cost= (self.gas_module_maintenance * self.wages * self.user_scale_up)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        super()._design()
        D = self.design_results
        if self.extreme:
            return (self.dosing_valve_power_demand * self.dosing_valve_daily_operation/24*self.user_scale_up+
                           D['Ideal power']/self.eta/self.extreme_compressor_efficiency)   # kWh/hr
        else:
            return (self.dosing_valve_power_demand * self.dosing_valve_daily_operation/24*self.user_scale_up+
                           D['Ideal power']/self.eta/self.compressor_efficiency) # kWh/hr

#%%
mscwo_reactor_module_path = ospath.join(g2rt_su_data_path, '_mscwo_reactor_module.csv')
@price_ratio()
class mSCWOReactorModule(SanUnit):
    '''
    Reactor unit that performs organic oxidation to CO2 by supercritical water oxidation.
    
    The following components should be included in system thermo object for simulation:
    H2O, OtherSS, N2O, NH3, CO2, O2, sCOD, xCOD, Tissue, N2, O2

    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel

    Parameters
    ----------
    ins : Iterable(stream)
        compressed air (221 bar) , homogenized feces solids
    outs : Iterable(stream)
        gas product mixture, N2O, treated effluent, ash

    References
    ----------
    [1] YEE et al. Water oxidation non-sewered single unit toilet system.
    
    https://patentimages.storage.googleapis.com/57/6a/81/72a168a92be44c/WO2023288331A1.pdf
    '''
    _N_ins = 2
    _N_outs = 4
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', 
                 material_replacement_cost = None, reactor_system_cost = None,
                 user_scale_up=1, extreme=False, capital_scale_func=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.user_scale_up = user_scale_up
        self.extreme = extreme
        self.capital_scale_func = capital_scale_func
        data = load_data(path=mscwo_reactor_module_path)
        for para in data.index:
            if material_replacement_cost is not None and para == "material_replacement_cost":
                self.material_replacement_cost = material_replacement_cost
                continue  # Skip the material_replacement_cost index if it's provided
            if reactor_system_cost is not None and para == "reactor_system_cost":
                self.reactor_system_cost = reactor_system_cost
                continue  # Skip the reactor_system_cost index if it's provided
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])

    def _run(self):
        compressed_air, feces = self.ins
        compressed_air.P = self.operating_pressure * 1e5
        compressed_air_vol = feces.F_vol * (self.reactor_volume-self.feces_batch_volume)/self.feces_batch_volume
        compressed_air.ivol['N2'] = compressed_air_vol * 0.79
        compressed_air.ivol['O2'] = compressed_air_vol* 0.21
        initial_enthalpy_flow = compressed_air.H + feces.H #kJ/hr
        
        gas_product, N2O, liquid_effluent, ash = self.outs
        gas_product.P = N2O.P = liquid_effluent.P = ash.P= 101325 #release into atmosphere
        N2O.phase = 'g'
        ash.phase = 's'
        N2O.imass['N2O'] = (feces.imass['NH3'] + feces.imass['NonNH3']/14*44/2
                            ) * self.N2O_nitrogen_fraction #kg/hr
        
        compressed_air.T = feces.T = self.operating_temperature
        final_enthalpy_flow = (compressed_air.H + 
                               CP.PropsSI('H', 'T', feces.T, 'P', self.operating_pressure * 1e5, 'Water')/1000*
                               feces.imass['H2O']) #kJ/hr
        #use CoolProp package for more accurate enthalpy for supercritical water
        mc = feces.imass['H2O']/ feces.F_mass
        chemical_energy = (feces.F_mass - feces.imass['H2O']) * self.HHV_feces_solids * 1000 #kJ/hr
        required_energy = final_enthalpy_flow - initial_enthalpy_flow #kJ/hr
        if self.extreme:
            recovered_energy = chemical_energy * self.extreme_energy_recovery_efficiency * self.carbon_conversion_efficiency #kJ/hr
            self.power_input = (required_energy/self.extreme_heating_energy_efficiency - recovered_energy)/3600 #kW
        else:
            recovered_energy = chemical_energy * self.energy_recovery_efficiency * self.carbon_conversion_efficiency #kJ/hr
            self.power_input = (required_energy/self.heating_energy_efficiency - recovered_energy)/3600 #kW

        # self.power_input = required_energy/self.heating_energy_efficiency/3600 #kW
        # self.heat_output = recovered_energy #kJ/hr
        
        liquid_effluent.copy_like(feces)
        liquid_effluent.imass['NH3'] = (feces.imass['NH3'] + feces.imass['NonNH3']
                                        ) * self.ammonium_nitrogen_fraction #kg/hr
        liquid_effluent.imass['NonNH3'] = (feces.imass['NH3'] + feces.imass['NonNH3']
                                        ) * (1-self.ammonium_nitrogen_fraction-
                                             self.N2O_nitrogen_fraction) #kg/hr, 
        gas_product.copy_like(compressed_air)
        gas_product.imass['CO2'] = (feces.imass['sCOD'] + 
                                    feces.imass['xCOD'])* self.carbon_conversion_efficiency* self.carbon_COD_ratio/12*44 #kg/hr
        gas_product.imol['O2'] -= gas_product.imol['CO2']
        ash.imass['xCOD'] = liquid_effluent.imass['xCOD'] * (1 - self.carbon_conversion_efficiency)
        ash.imass['P'] = feces.imass['P'] * (1-self.P_loss) #assume P all go to solid phase
        gas_product.imass['P'] = feces.imass['P'] * self.P_loss
        liquid_effluent.imass['K'] = feces.imass['P'] * (1-self.K_loss)
        gas_product.imass['K'] = feces.imass['K'] * self.K_loss
        ash.imass['WoodAsh'] = (feces.F_mass * (1-mc) * self.feces_ash_content) # kg ash /hr
        liquid_effluent.imass['sCOD'] *= (1-self.carbon_conversion_efficiency)
        liquid_effluent.imass['xCOD'] = 0
        liquid_effluent.imass['Tissue'] = 0
        gas_product.P = N2O.P = liquid_effluent.P = ash.P= 101325 #release into atmosphere
        gas_product.T = N2O.T = liquid_effluent.T = ash.T = self.operating_temperature

    def _init_lca(self):
        self.construction = [
            Construction('stainless_steel',linked_unit=self, 
                         item='StainlessSteel', 
                         quantity_unit='kg'),
            ]
        
    def _design(self):
        design=self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.stainless_steel_weight
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Reactor'] = self.reactor_system_cost
        C['Injector'] = self.injector_cost 
        C['Misc.parts'] = self.miscellaneous_cost_ratio*(C['Reactor'] +
                               C['Injector'])
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio
        if self.power_input>=0:
            self.power_utility(self.power_input)
        else:
            self.power_utility.production = -self.power_input # kW
        # self.power_utility(self.power_input)    
        # self.heat_utilities(-self.heat_output)
        
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = ((self.capital_scale_func(total_equipment) if self.capital_scale_func else total_equipment)*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr

    def _calc_maintenance_labor_cost(self): #USD/hr
        if self.extreme:
            maintenance_labor_cost= (self.reactor_maintenance * self.wages *self.user_scale_up)*0.1
        else:
            maintenance_labor_cost= (self.reactor_maintenance * self.wages *self.user_scale_up)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return self.power_input


#%%
mscwo_concentrator_module_path = ospath.join(g2rt_su_data_path, '_mscwo_concentrator_module.csv')
@price_ratio()
class mSCWOConcentratorModule(SanUnit):
    '''
    Concentrator unit that evaporizes water by using heat from mSCWO effluent and external heat.
    
    The following components should be included in system thermo object for simulation:
    H2O, OtherSS, N2O, NH3, CO2, O2, sCOD, xCOD, Tissue

    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, HeatingUnit, ElectricMotor, Pump, Fan, Polyethylene

    Parameters
    ----------
    ins : Iterable(stream)
        RO reject liquid, mscwo effluent liquid, mscwo effluent ash
    outs : Iterable(stream)
        Condensed effluent, recirculation water, fugitive N2O, fugitive CH4, fugitive NH3, water vapor
    excess_water_recirculate_ratio: float, optional, 0 to 1
        fraction of excess water that recirculates to ultrafiltration. Defaults to 1.
    References
    ----------
    [1] YEE et al. Water oxidation non-sewered single unit toilet system.
    
    https://patentimages.storage.googleapis.com/57/6a/81/72a168a92be44c/WO2023288331A1.pdf
    '''
    _N_ins = 3
    _N_outs = 6
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 excess_water_recirculate_ratio = 1,user_scale_up=1, extreme= False,
                 capital_scale_func=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.excess_water_recirculate_ratio = excess_water_recirculate_ratio
        self.user_scale_up = user_scale_up
        self.extreme = extreme
        self.capital_scale_func = capital_scale_func
        data = load_data(path=mscwo_concentrator_module_path)
        
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids)) + ('OtherSS',) #salts oversaturated and become solids
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])
    
    def _run(self):
        ro_waste_in, mscwo_liquid_in, mscwo_ash_in = self.ins
        # waste_out, N2O, CH4 = self.outs
        waste_out, excess_water, N2O, CH4, NH3_gas, water_vapor  = self.outs
        N2O.phase = CH4.phase = NH3_gas.phase = water_vapor.phase = 'g'
        recovered_energy_from_mscwo = (CP.PropsSI('H', 'T', mscwo_liquid_in.T, 'P', mscwo_liquid_in.P, 'Water')/1000*\
        mscwo_liquid_in.imass['H2O']) * self.heat_exchange_efficiency #kJ/hr
        waste_out.mix_from((self.ins[0],self.ins[1],self.ins[2]))
        
        solubles, solids = self.solubles, self.solids
        #Calculate water and solids in the condensed effluent
        mc_in = (ro_waste_in.imass['H2O'] + 
                 mscwo_liquid_in.imass['H2O'] + 
                 mscwo_ash_in.imass['H2O']) / (ro_waste_in.F_mass + 
                                                  mscwo_liquid_in. F_mass+
                                                  mscwo_ash_in.F_mass) # fraction
        mc_out = self.moisture_content_out/100 #convert to fraction
        if mc_in < mc_out*0.999:
            mc_out = mc_in
        TS_in = ro_waste_in.imass[solids].sum() + mscwo_liquid_in.imass[solids].sum() + mscwo_ash_in.imass[solids].sum() # kg TS dry/hr
        
        waste_out.imass['H2O'] = TS_in/(1-mc_out)*mc_out #kg water/hr
        
        #Calculate N2O and CH4 emissions
        CH4.imass['CH4'] = drying_CH4_to_air = \
            self.drying_CH4_emissions * self.carbon_COD_ratio * \
            (ro_waste_in.COD * ro_waste_in.F_vol + 
             mscwo_liquid_in.COD * mscwo_liquid_in.F_vol
             ) / 1000 # kg CH4 /hr
        drying_NH3_to_air = self.drying_NH3_emissions * (ro_waste_in.imass['NH3'] +
                                                         mscwo_liquid_in.imass['NH3']
                                                         ) # kg NH3 /hr
        
        N2O.imass['N2O'] = drying_NH3_to_air * self.NH3_to_N2O /14/2 * 44 # kg N2O /hr
        waste_out.imass['NH3'] = ro_waste_in.imass['NH3'] + mscwo_liquid_in.imass['NH3'] - drying_NH3_to_air

        NH3_gas.imass['NH3'] = drying_NH3_to_air * (1-self.NH3_to_N2O) # kg NH3 /hr
        excess_water.imass['H2O'] = self.excess_water_recirculation * (ro_waste_in.imass['H2O'] + 
                                                                       mscwo_liquid_in.imass['H2O'] - 
                                                                       waste_out.imass['H2O']) * self.excess_water_recirculate_ratio
        excess_water.imass[solubles] = self.excess_water_recirculation * (ro_waste_in.imass[solubles]+
                                                                          mscwo_liquid_in.imass[solubles])* self.excess_water_recirculate_ratio
        water_vapor.imass['H2O'] = (ro_waste_in.imass['H2O'] + 
                                    mscwo_liquid_in.imass['H2O'] - 
                                    waste_out.imass['H2O'] - 
                                    excess_water.imass['H2O']) #kg H2O/hr
        # Store the calculated value of water_vapor.imass['H2O'] for use in the _cost function
        if self.extreme:
            if water_vapor.imass['H2O'] * self.extreme_energy_required_to_evaporize_water >=  (recovered_energy_from_mscwo/3600+ self.fan_power_demand*self.fan_daily_operation/24):
                self.required_energy_input = water_vapor.imass['H2O'] * self.extreme_energy_required_to_evaporize_water - recovered_energy_from_mscwo/3600 #kW
            else:
                self.required_energy_input = self.fan_power_demand*self.fan_daily_operation/24 #minimal fan operation
        else:
            if water_vapor.imass['H2O'] * self.energy_required_to_evaporize_water >=  (recovered_energy_from_mscwo/3600+ self.fan_power_demand*self.fan_daily_operation/24):
                self.required_energy_input = water_vapor.imass['H2O'] * self.energy_required_to_evaporize_water - recovered_energy_from_mscwo/3600 #kW
            else:
                self.required_energy_input = self.fan_power_demand*self.fan_daily_operation/24 #minimal fan operation
        #Calculate COD
        self.drying_CO2_to_air = (self.drying_CO2_emissions * self.carbon_COD_ratio
                             * (ro_waste_in.COD * ro_waste_in.F_vol + 
                              mscwo_liquid_in.COD * mscwo_liquid_in.F_vol
                              ) / 1000) # kg CO2 /hr
        # 44/12/16 are the molecular weights of CO2, C, and CH4, respectively
        waste_out.imass['sCOD'] =  (waste_out.imass['sCOD'] + mscwo_liquid_in.imass['sCOD']) \
            -((self.drying_CO2_to_air/44*12+drying_CH4_to_air/16*12) / self.carbon_COD_ratio + 
              excess_water.imass['sCOD'])
        water_vapor.imass['P'] = self.P_loss * waste_out.imass['P']
        water_vapor.imass['K'] = self.K_loss * waste_out.imass['K']
        waste_out.imass['P'] *= (1-self.P_loss)
        waste_out.imass['K'] *= (1-self.K_loss)
        
    
    def _init_lca(self):
        self.construction = [
            Construction('stainless_steel',linked_unit=self, 
                         item='StainlessSteel', 
                         quantity_unit='kg'),
            Construction('heating_unit',linked_unit=self, 
                         item='HeatingUnit', 
                         quantity_unit='ea'),
            Construction('electric_motor',linked_unit=self, 
                         item='ElectricMotor', 
                         quantity_unit='kg'),
            Construction('pump',linked_unit=self, 
                         item='Pump',
                         quantity_unit='ea'),
            Construction('fan',linked_unit=self, 
                         item='Fan',
                         quantity_unit='kg'),
            Construction('polyethylene',linked_unit=self,
                         item='Polyethylene', 
                         quantity_unit='kg'),
            ]
        
    def _design(self):
        design=self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.stainless_steel_weight
        design['HeatingUnit'] = constr[1].quantity = 1
        design['ElectricMotor'] = constr[2].quantity = self.motor_weight
        design['Pump'] = constr[3].quantity = 1 *(self.user_scale_up ** 0.6)
        design['Fan'] = constr[4].quantity = self.fan_quantity*0.2 #one fan is 0.2 kg
        design['Polyethylene'] = constr[5].quantity = self.polyethylene_weight
        
    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        #C['Stainless steel'] = self.stainless_steel_cost * D['StainlessSteel']
        C['HeatingUnit'] = self.heating_unit_purchase_cost * D['HeatingUnit']
        C['ElectricMotor'] = self.motor_purchase_cost * D['ElectricMotor']
        C['Pump'] = self.pump_purchase_cost * D['Pump'] *(self.user_scale_up ** 0.6)
        #C['Polyethylene'] = self.polyethylene_cost * D['Polyethylene']
        C['Fan'] = self.fan_purchase_cost * self.fan_quantity
        C['Disc'] = self.disc_purchase_cost
        C['PhaseSeparator'] = self.phase_separator_cost
        C['Others'] = (self.enclosure_purchase_cost + 
                       self.open_concentrator_vessel_purchase_cost+
                       self.axle_purchase_cost+
                       self.thermistor_cost                       
                       )
        C["Misc.parts"] = self.miscellaneous_cost_ratio*(C['HeatingUnit'] +
                               C['ElectricMotor']+
                               C['Pump'] +
                               C['Fan'] +
                               C['Disc'] +
                               C['PhaseSeparator'] +
                               C['Others']
                               ) 
        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio
        
        self.power_utility(self.pump_power_demand * self.pump_daily_operation/24+
                           self.required_energy_input) # kW
        
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = ((self.capital_scale_func(total_equipment) if self.capital_scale_func else total_equipment)*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
   
    def _calc_maintenance_labor_cost(self): #USD/hr
        if self.extreme:
            maintenance_labor_cost= (self.concentrator_maintenance * self.wages *self.user_scale_up)*0.1 #USD/yr
        else:
            maintenance_labor_cost= (self.concentrator_maintenance * self.wages *self.user_scale_up) #USD/yr
        return maintenance_labor_cost / (365*24)
   
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.pump_power_demand * self.pump_daily_operation/24 + 
                self.required_energy_input) #kW

#%%
vr_combustor_path = ospath.join(g2rt_su_data_path, '_vr_combustor.csv')
@price_ratio()
class VolumeReductionCombustor(SanUnit):
    '''
    Combustor unit that uses wood pellets biofuels to burn dry feces solid cakes.
    There is no energy or heat recovery in this unit.
    
    The following components should be included in system thermo object for simulation:
    H2O, N, K, P OtherSS, N2O, CH4, NO, SO2, NH3

    The following impact items should be pre-constructed for life cycle assessment:
    Steel, Aluminum

    Parameters
    ----------
    ins : Iterable(WasteStream)
        Dewatered feces solid cakes, wood pellets, air
    outs : Iterable(WasteStream)
        Wood Ash, hot gas, fugitive N2O, fugitive CH4, fugitive NO, fugitive SO2
    if_sludge_service: bool
        If share combustor unit among multiple volume reduction toilets
        (assume 120 users per combustor unit,
        or 20 volume reduction toilets serving a population of 6 users per toilet).
    lifetime : 10 years
        This is used to estimate the amount of wood pellet used for TEA/LCA. Default to 10 years.

    References
    ----------
    [1] YEE et al. VOLUME REDUCTION SOLIDS TREATMENT SYSTEM. 
    https://patents.google.com/patent/WO2023288327A1/en?oq=WO2023288327A1
    '''
    _N_ins = 3
    _N_outs = 7
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_sludge_service = True, lifetime = 10,ppl = None, CH4_emission_factor=None,
                 user_scale_up=1,extreme=False, capital_scale_func=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1, lifetime = lifetime)
        self.if_sludge_service = if_sludge_service
        self.lifetime = lifetime
        self.user_scale_up = user_scale_up
        self.extreme = extreme
        self.capital_scale_func = capital_scale_func
        if ppl is None:
            self.ppl = 6
        else:
            self.ppl = ppl

        data = load_data(path=vr_combustor_path)
        for para in data.index:
            if CH4_emission_factor is not None and para == "fugitive_CH4_emission_factor":
                self.fugitive_CH4_emission_factor = CH4_emission_factor
                continue  # Skip the material_replacement_cost index if it's provided
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    def _run(self):
        solid_cakes, wood_pellets, air = self.ins
        ash, gas, CH4, N2O, NO, SO2, NH3 = self.outs
        cmps = self.components
        ash.copy_like(self.ins[0])
        wood_pellets.imass['WoodPellet'] = solid_cakes.F_mass * self.biofuel_to_solids
        self.wood_pellets_kgphr = wood_pellets.imass['WoodPellet']
        wood_pellets.imass['H2O'] = wood_pellets.imass['WoodPellet'] * self.wood_moisture_content
        gas.phase = N2O.phase = CH4.phase = NO.phase = SO2.phase = NH3.phase = 'g'
        gas.T = self.combustion_temperature
        # Moisture content
        mc = solid_cakes.imass['H2O']/solid_cakes.F_mass if solid_cakes.F_mass!=0 else 0
        
        if mc > 0.351: # allow a small error
            warn(f'Moisture content of the influent is {mc:.1%}, '
                'larger than the maximum allowed level of 35%.')

        # ash.empty()
        
        ash_prcd = (solid_cakes.F_mass * (1-mc) * self.feces_ash_content + 
               wood_pellets.imass['WoodPellet']*(1-self.wood_moisture_content)*self.wood_ash_content) # kg ash /hr
        
        NPKCaMg = ('NH3','NonNH3','N', 'P', 'K','Ca','Mg','OtherSS')
        for cmp in cmps:
            ash.imass[cmp.ID] = 0 if cmp.ID not in NPKCaMg else ash.imass[cmp.ID]
        for element in NPKCaMg:
            if element in ('NH3', 'NonNH3'):
                ash.imass[element] *= (1 - getattr(self, 'combustion_N_loss'))
            elif element in ('Ca','Mg','N'):
                continue
            else:
                gas.imass[element] = getattr(self, f'combustion_{element}_loss')*ash.imass[element]
                ash.imass[element] *= (1 - getattr(self, f'combustion_{element}_loss'))
        if ash_prcd - ash.imass[NPKCaMg].sum() > 0:
            ash.imass['WoodAsh'] = (ash_prcd - ash.imass[NPKCaMg].sum())
        else: ash.imass['WoodAsh'] = 0
        ash.imass['H2O'] = 0.025 * ash.F_mass # kg H2O / hr with 2.5% moisture content
        #CH4 emissions
        CH4.imass['CH4'] = solid_cakes.COD * self.carbon_COD_ratio * solid_cakes.F_vol/1e3* self.fugitive_CH4_emission_factor #kg/hr
        # N2O emissions
        if self.extreme:
            N2O.imass['N2O'] = (solid_cakes.imass['NH3']+solid_cakes.imass['NonNH3']) * self.extreme_N2O_emission_factor/28*44
        else:
            N2O.imass['N2O'] = (solid_cakes.imass['NH3']+solid_cakes.imass['NonNH3']) * self.N2O_emission_factor/28*44
        # NO emissions
        NO.imass['NO'] = (solid_cakes.imass['NH3']+solid_cakes.imass['NonNH3']) * self.NOx_emission_factor/14*30
        #SO2 emissions
        SO2.imass['SO2'] = solid_cakes.F_mass * self.SOx_emission_factor
        #NH3 emissions
        NH3.imass['NH3'] = (solid_cakes.imass['NH3']+solid_cakes.imass['NonNH3']) * self.NH3_emission_factor/14*17
        #gas
        gas.imass['CO2'] = ((solid_cakes.imass['sCOD']+solid_cakes.imass['xCOD']) * self.carbon_COD_ratio * self.combustion_C_loss /12 *44 
                            + cmps['WoodPellet']._i_C * wood_pellets.imass['WoodPellet']
                            #assume all carbon in wood pellets convert to CO2
                            )
        self.CO2_from_solids = (solid_cakes.imass['sCOD']+solid_cakes.imass['xCOD']) * self.carbon_COD_ratio * self.combustion_C_loss /12 *44 #kg/hr
        air.imass['O2'] = gas.imass['CO2']/44*32
        air.imass['N2'] = air.imass['O2']/32/0.21*0.79*28 #assume air content 21% O2 and 79% N2 
        gas.imass['N2'] = ((solid_cakes.imass['NH3']+solid_cakes.imass['NonNH3'])*self.combustion_N_loss - 
                           N2O.imass['N2O']/44*28 - NO.imass['NO']/30*14-NH3.imass['NH3']/17*14) + air.imass['N2']
        
        gas.imass['H2O'] = solid_cakes.imass['H2O'] + wood_pellets.imass['H2O'] - ash.imass['H2O']
        
        
    def _init_lca(self):
        self.construction = [
            Construction('Steel', linked_unit=self, item='Steel', quantity_unit='kg'),
            Construction('Aluminum', linked_unit=self, item='Aluminum', quantity_unit='kg'),
            Construction('WoodPellet', linked_unit=self, item='WoodPellet', quantity_unit='kg'),
            ]
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = self.steel_weight
        design['Aluminum'] = constr[1].quantity = self.aluminum_weight
        design['WoodPellet'] = constr[2].quantity = self.lifetime*365*24* self.wood_pellets_kgphr
        self.add_construction(add_cost = False)
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C["Combustor"] = self.combustor_cost * self.price_ratio #USD
        service_factor = self.ppl/120 if self.if_sludge_service else 1
        total_equipment = 0.
        for equipment, cost in C.items():
            C[equipment] = cost * service_factor
            total_equipment += cost
         
        self.add_OPEX = ((self.capital_scale_func(total_equipment) if self.capital_scale_func else total_equipment)*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_replacement_cost()+
                         self._calc_maintenance_labor_cost()) #USD/hr
        
    def _calc_replacement_cost(self):
        wood_pellet_cost = self.wood_pellets_kgphr * self.wood_pellets_cost #USD/hr
        service_factor = self.ppl/120 if self.if_sludge_service else 1
        return wood_pellet_cost* self.price_ratio * service_factor # USD/hr
            
    def _calc_maintenance_labor_cost(self): #USD/hr
        if self.extreme:
            maintenance_labor_cost= (self.combustor_operation_maintenance * self.wages *self.user_scale_up)*0.1
        else:
            maintenance_labor_cost= (self.combustor_operation_maintenance * self.wages *self.user_scale_up)
        service_factor = self.ppl/120 if self.if_sludge_service else 1
        return maintenance_labor_cost / (365*24) * service_factor
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        '''[float] Power draw, [kW], auto-calculated from `power_utility`.'''
        return self.power_utility.rate

#%%

class UFMixer(Mixer):
    '''
    Mixing ultrafiltration reject and liquid from solids separator before delivering to belt separator
    '''
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 rigorous=False, conserve_phases=False):
        Mixer.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic)
        
        
    # def _run(self):
    #     separated_liquid, uf_retentate = self.ins
    #     waste_out, = self.outs
    #     waste_out.mix_from(self.ins)
    #     # print(f"The flushing water flow is {flushing.imass['H2O']} kg/h.")

# %%
toilet_path = ospath.join(g2rt_su_data_path, '_g2rt_toilet.csv')

class FWMixer(Mixer):
    '''
    Mixing tap water and recycled water from RO to meet the flushing water demand
    '''
    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 rigorous=False, conserve_phases=False,N_user=None, N_toilet=None,
                 N_tot_user = None, if_flushing=True, flushing_water=None):
        Mixer.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic)
        self.N_user = N_user
        self.N_toilet = N_toilet
        self.N_tot_user = N_tot_user
        self.if_flushing = if_flushing
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])
        
        data = load_data(path=toilet_path)    
        for para in data.index:
            if flushing_water is not None and para == "flushing_water":
                self.flushing_water = flushing_water
                continue  # Skip the material_replacement_cost index if it's provided
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
                               
    def _run(self):
        tap_water, ro_permeate = self.ins
        flushing, = self.outs
        solids = self.solids
        solubles = self.solubles
        N_tot_user = self.N_tot_user or self.N_toilet*self.N_user
        flushing.imass['H2O'] = self.flushing_water*N_tot_user
        tap_water.imass['H2O'] = flushing.imass['H2O'] - ro_permeate.imass['H2O']
        for i in solids+solubles:
            flushing.imass[i] = tap_water.imass[i] + ro_permeate.imass[i]

        # flushing.mix_from((tap_water,ro_permeate))
        # print(f"The flushing water flow is {flushing.imass['H2O']} kg/h.")

#%%
class G2RTExcretion(Excretion):
    '''
    G2RT-specific excretion model. Thin subclass of :class:`~.Excretion` that
    represents feces/urine COD as explicit soluble/particulate component mass
    flows (`sCOD`/`xCOD`) instead of a lumped `_COD` concentration override,
    since G2RT's downstream separation units need component-level COD to
    model soluble vs. particulate removal. Uses the same underlying data
    (`_excretion.tsv`) and all other assumptions as :class:`~.Excretion`.

    See Also
    --------
    :class:`~.Excretion`
    '''

    @property
    def power_kW(self):
        '''[float] Power draw, [kW], auto-calculated from `power_utility`.'''
        return self.power_utility.rate

    def _set_COD(self, ur, fec, tot_COD):
        ur.imass['sCOD'] = tot_COD*(1-self.e_fec) # in kg/hr
        fec.imass['xCOD'] = tot_COD*self.e_fec # in kg/hr

#%%
vr_pasteurization_path = ospath.join(g2rt_su_data_path, '_vr_pasteurization.csv')

@price_ratio()
class VRpasteurization(SanUnit):
    '''
    Pasteurizer in volume reduction toilet to remove pathogens in the homogenized 
    solids. Heating is from electricity (Joule heater).
    
    Parameters
    -----------
    ins: Iterable(stream)
        solids: solids produced from the homogenizer in the volume reduction toilet.

    outs: Iterable(stream)
        treated solids: solids treated from pasteurization.
    temp_pasteurization : float
        Pasteurization temperature is 70°C or 343.15 K.
    solids_inlet_temp : float
        Temperature of solids from the inlet.
    heat_loss : float
        Heat loss during pasteurization process is assumed to be 30%
    
    See Also
    --------
    :class:`~.sanunits.SludgePasteurization`
    '''
    
    _N_ins = 1
    _N_outs = 1
    
    # Specific Heat capacity of water
    Cp_w = 4.184 # kJ kg^-1 K^-1
    # Specific Heat capacity of dry matter (sludge)
    Cp_dm = 1.231 # kJ kg^-1 K^-1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 heat_loss=0.3, solids_inlet_temp=273.15+25.6,
                 temp_pasteurization= 273.15+90,user_scale_up=1,  extreme = False,
                 capital_scale_func=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                          F_BM_default=1)
        self.heat_loss = heat_loss
        self.solids_inlet_temp = solids_inlet_temp
        self.temp_pasteurization = temp_pasteurization
        self.user_scale_up = user_scale_up
        self.extreme = extreme
        self.capital_scale_func = capital_scale_func

        data = load_data(path=vr_pasteurization_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    def _run(self):
        solids_in, = self.ins
        solids_out, = self.outs
        solids_out.copy_like(solids_in)
        
    def _init_lca(self):
        self.construction = [Construction("stainless_steel", linked_unit=self,
                                          item = "StainlessSteel", 
                                          quantity_unit= "kg"),
                             Construction("polyethylene", linked_unit=self,
                                          item = "Polyethylene",
                                          quantity_unit= "kg")]
        
    def _cost(self):
        C=self.baseline_purchase_costs
        C["Heating"]= (self.joule_heater_cost+
                       self.thermocouple_port_cost* self.thermocouple_port_quantity+
                       self.temperature_control
                       ) * self.price_ratio
        C["Valves"] = self.valve_cost * self.valve_quantity * self.price_ratio
        C["Tubing"] = self.tubing_cost * self.price_ratio
        C["Misc.parts"] = self.miscellaneous_cost_ratio*(C["Heating"]+
                               C["Valves"]+
                               C["Tubing"])
        solids_in, = self.ins
        solids_in.T = self.solids_inlet_temp + 273.15
        solids_out, = self.outs
        solids_out.T = self.temp_pasteurization + 273.15
        # Overall heat required for pasteurization
        temp_diff = self.temp_pasteurization - self.solids_inlet_temp #K
        Q_d = (solids_in.imass['H2O']*self.Cp_w + 
               (solids_in.F_mass-solids_in.imass['H2O'])*self.Cp_dm)*temp_diff #kJ/hr
        if self.extreme:
            Q_tot = Q_d/(1-self.extreme_heat_loss/100) # kJ/hr, 0% of the total generated is lost
        else:
            Q_tot = Q_d/(1-self.heat_loss/100) # kJ/hr, 30% of the total generated is lost
        self.heating_electricity = Q_tot/3600 #kWh/hr
        # self.add_heat_utility(unit_duty=self.Hnet, T_in=solids_in.T, 
        #                       T_out = solids_out.T, heat_transfer_efficiency=(1-self.heat_loss/100),
        #                       hxn_ok= True)
        self.power_utility(self.heating_electricity) #kWh/hr
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = ((self.capital_scale_func(total_equipment) if self.capital_scale_func else total_equipment)*self.material_replacement_cost/(365*24) + #USD/hr, assume 
                         #replacement cost a fraction of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
        
        # def _calc_replacement_cost(self): #USD/hr, assume 5% of CAPEX per year
        #     replacement_cost = 0.05*(
        #         C["Heating"]+ C["Valves"] + C["Tubing"] + C["Misc.parts"]
        #         ) #USD/yr
        #     return replacement_cost / (365*24)
    
    def _calc_maintenance_labor_cost(self): #USD/hr
        if self.extreme:
            maintenance_labor_cost= (self.pasteurizer_maintenance * self.wages )*0.1 #USD/yr
        else:
            maintenance_labor_cost= (self.pasteurizer_maintenance * self.wages ) #USD/yr
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return self.heating_electricity #kW


#%%
G2RT_homogenizer_path = ospath.join(g2rt_su_data_path, '_g2rt_homogenizer.csv')

@price_ratio()
class G2RThomogenizer(Copier):
    '''
    Homogenizer and buffer tanks in generation II reinveted toilets to  break 
    up solids [1].
    
    .. note:

        This is a non-reactive unit (i.e., the effluent is copied from the influent)

    The following components should be included in system thermo object for simulation:
    H2O, OtherSS.

    The following impact items should be pre-constructed for life cycle assessment:
    Steel.

    Parameters
    ----------
    ins : Iterable(stream)
        Influent stream.
    outs : Iterable(stream)
        Effluent stream, is copied from the influent.
    moisture_content_out : float
        Moisture content of the effluent solids stream.

    References
    ----------
    [1] YEE et al. Buffer tank separation and homogenization system. 
    https://patents.google.com/patent/WO2023288114A1/en?oq=WO2023288114A1
    [2] YEE et al. Volume reduction non-sewered single unit toilet system.
    https://patents.google.com/patent/WO2023288326A1/en?oq=WO2023288326A1
    
    See Also
    ---------
    :class:`~.sanunits.BiogenicRefineryGrinder`
    '''
    _N_ins = 1
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',ppl=None,
                 user_scale_up=1, extreme = False, capital_scale_func=None, **kwargs):
        Copier.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.user_scale_up = user_scale_up
        self.extreme = extreme
        self.capital_scale_func = capital_scale_func
        if ppl is None:
            self.ppl = 6
        else:
            self.ppl = ppl
            
        data = load_data(path=G2RT_homogenizer_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    def _init_lca(self):
        self.construction = [Construction("stainless_steel", linked_unit=self,
                                          item = "StainlessSteel", 
                                          quantity_unit= "kg"),
                             Construction("polyethylene", linked_unit=self,
                                          item = "Polyethylene",
                                          quantity_unit= "kg")]
    
    # def _run(self):
    #     waste_in = self.ins[0]
    #     waste_out = self.outs[0]
    #     waste_out.copy_like(waste_in)

    #     mc_in = waste_in.imass['H2O'] / waste_in.F_mass # fraction
    #     mc_out = self.moisture_content_out/100 #convert to fraction
    #     
    #     if mc_in < mc_out*0.999:
    #         raise RuntimeError(f'Moisture content of the influent stream ({mc_in:.2f}) '
    #                            f'is smaller than the desired moisture content ({mc_out:.2f}).')
    #     TS_in = waste_in.F_mass - waste_in.imass['H2O'] # kg TS dry/hr
    #     waste_out.imass['H2O'] = TS_in/(1-mc_out)*mc_out
    #     waste_out._COD = waste_in.COD * waste_in.F_vol / waste_out.F_vol
        
    def _design(self): #kg
        design = self.design_results
        constr = self.construction    
        self.design_results['StainlessSteel'] = constr[0].quantity = self.stainless_steel_weight
        self.design_results['Polyethylene'] = constr[1].quantity = self.polyethylene_weight
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C["Macerator"] = self.macerator_cost * self.price_ratio
        C["Misc.parts"] = self.miscellaneous_cost_ratio * C["Macerator"]
        
        self.power_utility(self.macerator_power_demand * self.macerator_daily_operation
                           *self.ppl/24) # kWh/hr
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = ((self.capital_scale_func(total_equipment) if self.capital_scale_func else total_equipment)*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost a fraction of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
    
    def _calc_maintenance_labor_cost(self): #USD/hr
        if self.extreme:
            maintenance_labor_cost= (self.homogenizer_maintenance * self.wages)*0.1 #USD/yr
        else:
            maintenance_labor_cost= (self.homogenizer_maintenance * self.wages) #USD/yr
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.macerator_power_demand * self.macerator_daily_operation
                           *self.ppl/24) #kW
    
#%%
vr_filter_press_path = ospath.join(g2rt_su_data_path, '_vr_filter_press.csv')

@price_ratio()
class VolumeReductionFilterPress(SanUnit): 
    
    '''
    A filter press unit for the dewatering of mixed excreta in volume reduction
    generation II reinveted toilet [1]
    
    The following componenet should be included in system thermo object for simulation:
    Water.
    
    The following impact items should be pre-constructed for life cycle assessment:
    Stainless steel, Pump.
    
    Parameters
    ----------
    ins: Iterable(stream)
      Pasteurized solids waste for dewatering treatment
    outs: Iterable(stream)
      Liquids and solids produced from filter press.
    sludge_moisture: float
      Moisture content of the solids cake after filter press [wt% water].
    
    References
    -----------
    [1] YEE et al. VOLUME REDUCTION SOLIDS TREATMENT SYSTEM. 
    https://patents.google.com/patent/WO2023288327A1/en?oq=WO2023288327A1
    
    [2] Bev Express multi-plate sheet filter (ME-10 model)
    https://417cb0.p3cdn1.secureserver.net/wp-content/Documents/Tech%20Sheets/Filtration%20Equipment/Plate%20&%20Frame/Beverage%20Express/ErtelAlsop_Bev_ExPRESS_Tech_Sheet.pdf?dl=1
    
    [3] https://www.earthshields.com/how-much-does-geotextile-fabric-cost/# accessed on yyyy-mm-dd
    
    [4] https://www.suezwaterhandbook.com/processes-and-technologies/liquid-sludge-treatment/filter-press/conventional-recessed-plate-filter-press#:~:text=The%20filter%20press%20consumes%20relatively,solids%20depending%20on%20sludge%20type.
    
    [5] https://multimedia.3m.com/mws/media/2113780O/3m-zeta-plus-vs-filter-press-for-beer-clarification-application-brief.pdf
        
    '''
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', 
                 ins=None, outs=(), thermo=None, 
                 init_with='WasteStream',
                 solids = (),user_scale_up=1,  extreme = False,
                 capital_scale_func=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                                solids=solids,
                                **kwargs)
        self.user_scale_up = user_scale_up
        self.extreme = extreme
        self.capital_scale_func = capital_scale_func
        data = load_data(path=vr_filter_press_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])
        
    def _init_lca(self):
        self.construction = [
            Construction('stainless_steel',linked_unit=self, 
                         item='StainlessSteel', 
                         quantity_unit='kg'),
            Construction('pump', linked_unit=self, 
                         item='Pump', 
                         quantity_unit='ea'),
            ]
    
    def _run(self):
        pasteurized_solids, = self.ins
        supernatant, solid_cakes = self.outs
        solid_cakes.phase = 's'
        solubles, solids = self.solubles, self.solids
        solid_cakes.copy_flow(pasteurized_solids,solids) #all solids go to sludge
        solid_cakes.imass[solids] = pasteurized_solids.imass[solids]*self.TSS_removal/100
        
        mc_in = pasteurized_solids.imass['H2O'] / pasteurized_solids.F_mass # fraction
        mc_out = self.moisture_content_out/100 #convert to fraction
        
        if mc_in < mc_out*0.999:
            mc_out = mc_in
        TS_in = pasteurized_solids.imass[solids].sum() # kg TS dry/hr
        TS_out = solid_cakes.imass[solids].sum()
        #calculate water and solid COD in the solid cakes
        solid_cakes.imass['H2O'] = TS_in/(1-mc_out)*mc_out
        solid_cakes.imass[solubles] = pasteurized_solids.imass[solubles]*\
            (TS_out/(1-mc_out)-TS_out)/(pasteurized_solids.F_mass-TS_in)
        supernatant.mass = pasteurized_solids.mass-solid_cakes.mass
        
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.filterpress_ss_weight
        design['Pump'] = constr[1].quantity = 1 *(self.user_scale_up ** 0.6)
        self.add_construction(add_cost = False)
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C["Filter Press"] = self.filterpress_purchase_cost * self.price_ratio #USD
        if self.extreme:
            self.power_utility(self.extreme_filterpress_energy_per_water*self.outs[0].imass['H2O']
                           /self.extreme_filterpress_energy_efficiency
                           ) # kW
        else:
            self.power_utility(self.filterpress_energy_per_water*self.outs[0].imass['H2O']
                           /self.filterpress_energy_efficiency
                           ) # kW
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = ((self.capital_scale_func(total_equipment) if self.capital_scale_func else total_equipment)*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
            
    def _calc_maintenance_labor_cost(self): #USD/hr
        if self.extreme:
            maintenance_labor_cost= (self.filterpress_maintenance * self.wages * self.user_scale_up)*0.1
        else:
            maintenance_labor_cost= (self.filterpress_maintenance * self.wages * self.user_scale_up)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        if self.extreme:
            return (self.extreme_filterpress_energy_per_water*self.outs[0].imass['H2O']
                           /self.extreme_filterpress_energy_efficiency
                           ) # kW
        else:
            return (self.filterpress_energy_per_water*self.outs[0].imass['H2O']
                           /self.filterpress_energy_efficiency
                           ) # kW

#%%
vr_concentrator_path = ospath.join(g2rt_su_data_path, '_vr_concentrator.csv')

@price_ratio()
class VRConcentrator(SanUnit): 
    '''
    This concentrator unit is used in solids treatmnet in volume reduction 
    generation II reinveted toilet [1]. The heat source is from a heating coil.
#TODO: consider installing a heat exchange to receive heat from pasteurization.
    The following components should be included in system thermo object for simulation:
    H2O, N, CH4, N2O.
    
    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, HeatingUnit, ElectricMotor, Pump, Polyethylene.

    Parameters
    ----------
    ins : Iterable(stream)
        RO reject liquid.
    outs : Iterable(stream)
        Condensed effluent, recirculation water, fugitive N2O, fugitive CH4, fugitive NH3, water vapor
    excess_water_recirculate_ratio: float, optional, 0 to 1
        fraction of excess water that recirculates to ultrafiltration. Defaults to 1.
    Warnings
    --------
    Energy balance is not performed for this unit.

    References
    ----------
    [1] YEE et al. VOLUME REDUCTION SOLIDS TREATMENT SYSTEM. 
       https://patents.google.com/patent/WO2023288327A1/en?oq=WO2023288327A1

    See Also
    --------
    :class:`~.sanunits.BiogenicRefineryHHXdryer`
    :class:`~.sanunits.BiogenicRefineryHHX`
    '''
    _N_ins = 1
    _N_outs = 6
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 excess_water_recirculate_ratio =1 ,user_scale_up=1, extreme = False,
                 capital_scale_func=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.user_scale_up = user_scale_up
        self.excess_water_recirculate_ratio = excess_water_recirculate_ratio
        self.extreme = extreme
        self.capital_scale_func = capital_scale_func
        data = load_data(path=vr_concentrator_path)
        
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids)) + ('OtherSS',) #dissolved salts dry out
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])
    
    def _run(self):
        waste_in, = self.ins
        waste_out, excess_water, N2O, CH4, NH3_gas, water_vapor  = self.outs
        N2O.phase = CH4.phase = NH3_gas.phase = water_vapor.phase = 'g'
        waste_out.copy_like(self.ins[0])
        solubles, solids = self.solubles, self.solids
        #Calculate water and solids in the condensed effluent
        mc_in = waste_in.imass['H2O'] / waste_in.F_mass # fraction
        mc_out = self.moisture_content_out/100 #convert to fraction
        if mc_in < mc_out*0.999:
            mc_out = mc_in
        TS_in = waste_in.imass[solids].sum() # kg TS dry/hr
        waste_out.imass['H2O'] = TS_in/(1-mc_out)*mc_out #kg water/hr
        #Calculate N2O and CH4 emissions
        CH4.imass['CH4'] = drying_CH4_to_air = \
            self.drying_CH4_emissions * self.carbon_COD_ratio * \
            waste_in.COD * waste_in.F_vol / 1000 # kg CH4 /hr
        drying_NH3_to_air = self.drying_NH3_emissions * waste_in.imass['NH3'] # kg NH3 /hr
        N2O.imass['N2O'] = drying_NH3_to_air * self.NH3_to_N2O /14/2*44 # kg N2O /hr
        waste_out.imass['NH3'] = waste_in.imass['NH3'] - drying_NH3_to_air

        NH3_gas.imass['NH3'] = drying_NH3_to_air * (1-self.NH3_to_N2O) # kg NH3 /hr
        excess_water.imass['H2O'] = self.excess_water_recirculation * (waste_in.imass['H2O'] - waste_out.imass['H2O']) * self.excess_water_recirculate_ratio
        excess_water.imass[solubles] = self.excess_water_recirculation * waste_in.imass[solubles] * self.excess_water_recirculate_ratio
        water_vapor.imass['H2O'] = waste_in.imass['H2O'] - waste_out.imass['H2O'] -excess_water.imass['H2O']  #kg H2O/hr
        # Store the calculated value of water_vapor.imass['H2O'] for use in the _cost function
        self.water_vapor_H2O = water_vapor.imass['H2O'] 
        #Calculate COD
        self.drying_CO2_to_air = (self.drying_CO2_emissions * self.carbon_COD_ratio
                             * (waste_in.COD * waste_in.F_vol-
                                excess_water.COD * excess_water.F_vol
                                ) / 1000) # kg CO2 /hr
        # 44/12/16 are the molecular weights of CO2, C, and CH4, respectively
        waste_out.imass['sCOD'] -=  ((self.drying_CO2_to_air/44*12+drying_CH4_to_air/16*12)/ self.carbon_COD_ratio + 
                                     excess_water.imass['sCOD'])
        water_vapor.imass['P'] = self.P_loss * waste_out.imass['P']
        water_vapor.imass['K'] = self.K_loss * waste_out.imass['K']
        waste_out.imass['P'] *= (1-self.P_loss)
        waste_out.imass['K'] *= (1-self.K_loss)
    
    def _init_lca(self):
        self.construction = [
            Construction('stainless_steel',linked_unit=self, 
                         item='StainlessSteel', 
                         quantity_unit='kg'),
            Construction('heating_unit',linked_unit=self, 
                         item='HeatingUnit', 
                         quantity_unit='ea'),
            Construction('electric_motor',linked_unit=self, 
                         item='ElectricMotor', 
                         quantity_unit='kg'),
            Construction('pump',linked_unit=self, 
                         item='Pump',
                         quantity_unit='ea'),
            Construction('fan',linked_unit=self, 
                         item='Fan',
                         quantity_unit='kg'),
            Construction('polyethylene',linked_unit=self,
                         item='Polyethylene', 
                         quantity_unit='kg'),
            ]
        
    def _design(self):
        design=self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.stainless_steel_weight
        design['HeatingUnit'] = constr[1].quantity = 1
        design['ElectricMotor'] = constr[2].quantity = self.motor_weight
        design['Pump'] = constr[3].quantity = 1 *(self.user_scale_up ** 0.6)
        design['Fan'] = constr[4].quantity = self.fan_quantity*0.2 #one fan is 0.2 kg
        design['Polyethylene'] = constr[5].quantity = self.polyethylene_weight
        
    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        #C['Stainless steel'] = self.stainless_steel_cost * D['StainlessSteel']
        C['HeatingUnit'] = self.heating_unit_purchase_cost * D['HeatingUnit']
        C['ElectricMotor'] = self.motor_purchase_cost * D['ElectricMotor']
        C['Pump'] = self.pump_purchase_cost * D['Pump'] *(self.user_scale_up ** 0.6)
        #C['Polyethylene'] = self.polyethylene_cost * D['Polyethylene']
        C['Fan'] = self.fan_purchase_cost * self.fan_quantity
        C['Disc'] = self.disc_purchase_cost
        C['Others'] = (self.enclosure_purchase_cost + 
                       self.open_concentrator_vessel_purchase_cost+
                       self.axle_purchase_cost+
                       self.thermistor_cost                       
                       )
        C["Misc.parts"] = self.miscellaneous_cost_ratio*(C['HeatingUnit'] +
                               C['ElectricMotor']+
                               C['Pump'] +
                               C['Fan'] +
                               C['Disc'] + 
                               C['Others']
                               ) 
        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio
        if self.extreme:
            self.power_utility(self.pump_power_demand * self.pump_daily_operation/24+
                            self.water_vapor_H2O * self.extreme_energy_required_to_evaporize_water
                            ) # kW
        else:
            self.power_utility(self.pump_power_demand * self.pump_daily_operation/24+
                            self.water_vapor_H2O * self.energy_required_to_evaporize_water
                            ) # kW
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = ((self.capital_scale_func(total_equipment) if self.capital_scale_func else total_equipment)*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
   
    def _calc_maintenance_labor_cost(self): #USD/hr
        if self.extreme:
            maintenance_labor_cost= (self.concentrator_maintenance * self.wages*self.user_scale_up)*0.1 #USD/yr
        else:
            maintenance_labor_cost= (self.concentrator_maintenance * self.wages*self.user_scale_up) #USD/yr
        return maintenance_labor_cost / (365*24)
   
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        if self.extreme:
            return (self.pump_power_demand * self.pump_daily_operation/24+
                            self.water_vapor_H2O * self.extreme_energy_required_to_evaporize_water
                            ) # kW
        else:
            return (self.pump_power_demand * self.pump_daily_operation/24+
                            self.water_vapor_H2O * self.energy_required_to_evaporize_water
                            ) # kW
   
    
#%%
g2rt_liquids_tank_path = ospath.join(g2rt_su_data_path, '_g2rt_liquids_tank.csv')
@price_ratio()

class G2RTLiquidsTank(Mixer):
    '''
    Liquids storage unit for generation II reinveted toilets to accumulate enough
    liquid waste before ultrafiltration
    
    This is a non-reactive unit, (i.e., the effluent is copied from the mix of influent).
    
    The following impact items should be pre-constructed for life cycle assessment:
    Pump, Polyethylene, Polycarbonate
    
    References
    ----------
    [1] Shyu et al., The NEWgeneratorTM Non-Sewered Sanitation System:
    Long-Term Field Testing at an Informal Settlement Community in
    EThekwini Municipality, South Africa.
    Journal of Environmental Management 2021, 296, 112921.
    https://doi.org/10.1016/j.jenvman.2021.112921.

    See Also
    --------
    :class:`~.sanunits.Copier`
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', user_scale_up = 1, F_BM_default=None, isdynamic=False,
                 rigorous=False, conserve_phases=False,  extreme = False,):
        Mixer.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic)
        self.user_scale_up = user_scale_up
        self.extreme = extreme
        data = load_data(path=g2rt_liquids_tank_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

    def _init_lca(self):
        self.construction = [
            Construction(item='Polyethylene', linked_unit=self, quantity_unit='kg'),
            Construction(item='Polycarbonate', linked_unit=self, quantity_unit='kg'),
            Construction(item='Pump', linked_unit=self, quantity=1., quantity_unit='ea'),
            ]
    
    def _run(self):
        s_out, = self.outs
        s_out.mix_from(self.ins, vle=self.rigorous,
                       conserve_phases=getattr(self, 'conserve_phases', None))
        V = s_out.vapor_fraction
        if V == 0:
            self._B = 0
        elif V == 1:
            self._B = np.inf
        else:
            self._B = V / (1 - V)
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Polyethylene'] = constr[0].quantity = self.liquids_tank_polyethylene_weight *(self.user_scale_up ** 0.6)
        design['Polycarbonate'] = constr[1].quantity = self.liquids_tank_polycarbonate_weight *(self.user_scale_up ** 0.6)
        design['Pump'] = constr[2].quantity = 1 *(self.user_scale_up ** 0.6)
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tank'] = self.liquids_tank_cost*(self.user_scale_up ** 0.6)
        C['Tubing'] = self.liquids_tank_tubing
        C['Pump'] = self.liquids_tank_feed_pump_cost *(self.user_scale_up ** 0.6)
        C['Misc.parts'] = self.miscellaneous_cost_ratio*(C['Tank']+
                                                         C['Tubing']+
                                                         C['Pump'])
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

        power_demand = self.pump_energy*self.outs[0].F_mass/1000
        self.power_utility(power_demand)  # kW
        
    def _calc_replacement_cost(self):
        pump_replacement_cost = self.liquids_tank_feed_pump_cost/ self.pump_lifetime *(self.user_scale_up ** 0.6)
        return pump_replacement_cost/(365 * 24)* self.price_ratio # USD/hr

    def _calc_maintenance_labor_cost(self):
        if self.extreme:
            liquids_tank_maintenance_labor_cost = (
                self.labor_pump_replacement *self.user_scale_up +
                self.labor_tank_cleaning
                ) * self.wages * 0.1
        else:
            liquids_tank_maintenance_labor_cost = (
                self.labor_pump_replacement *self.user_scale_up +
                self.labor_tank_cleaning
                ) * self.wages
        return liquids_tank_maintenance_labor_cost/(365 * 24) # USD/hr
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return self.pump_energy*self.outs[0].F_mass/1000 #kW

#%%
g2rt_solids_tank_path = ospath.join(g2rt_su_data_path, '_g2rt_solids_tank.csv')
@price_ratio()

class G2RTSolidsTank(Copier):
    '''
    Solids storage unit for generation II reinveted toilets to accumulate enough
    solid waste before homogenizer.
    
    This is a non-reactive unit (i.e., the effluent is copied from the influent).
    
    The following impact items should be pre-constructed for life cycle assessment:
    Polyethylene, Polycarbonate
    
    References
    ----------
    [1] Shyu et al., The NEWgeneratorTM Non-Sewered Sanitation System:
    Long-Term Field Testing at an Informal Settlement Community in
    EThekwini Municipality, South Africa.
    Journal of Environmental Management 2021, 296, 112921.
    https://doi.org/10.1016/j.jenvman.2021.112921.

    See Also
    --------
    :class:`~.sanunits.Copier`
    '''
    def __init__(self, ID='', ins=None, outs=(),  thermo=None, init_with='WasteStream',
                 user_scale_up=1, extreme = False, **kwargs):
        Copier.__init__(self, ID, ins, outs, thermo, init_with)
        self.user_scale_up = user_scale_up
        self.extreme = extreme

        data = load_data(path=g2rt_solids_tank_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [
            Construction(item='Polyethylene', linked_unit=self, quantity_unit='kg'),
            Construction(item='Polycarbonate', linked_unit=self, quantity_unit='kg'),
            ]
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Polyethylene'] = constr[0].quantity = self.solids_tank_polyethylene_weight *(self.user_scale_up ** 0.6)
        design['Polycarbonate'] = constr[1].quantity = self.solids_tank_polycarbonate_weight *(self.user_scale_up ** 0.6)
        self.add_construction(add_cost=False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Tank'] = self.solids_tank_cost *(self.user_scale_up ** 0.6)
        C['Piping'] = self.solids_tank_piping
        C['Misc.parts'] = self.miscellaneous_cost_ratio*(C['Tank']+
                                                         C['Piping'])
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        self.add_OPEX =  self._calc_maintenance_labor_cost()

    def _calc_maintenance_labor_cost(self):
        if self.extreme:
            solids_tank_maintenance_labor_cost = (
                self.labor_tank_cleaning * self.user_scale_up
                ) * self.wages * 0.1
        else:
            solids_tank_maintenance_labor_cost = (
                self.labor_tank_cleaning * self.user_scale_up
                ) * self.wages
        return solids_tank_maintenance_labor_cost/(365 * 24) # USD/hr
        #Macerator in the homogenizer functions as a pump. All the associated cost
        #should refere to sanunits.G2RThomogenizer

    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day

    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day

    @property
    def power_kW(self):
        '''[float] Power draw, [kW], auto-calculated from `power_utility`.'''
        return self.power_utility.rate


#%%
vr_drying_tunnel_path = ospath.join(g2rt_su_data_path, '_vr_dryingtunnel.csv')
@price_ratio()
class VRdryingtunnel(SanUnit):
    '''
    This drying tunnel unit is used to produce solids cakes in volume reduction 
    generation II reinveted toilet [1]. The heat source is from a heating coil.
    #TODO: consider installing a heat exchange to receive heat from pasteurization.
    The following components should be included in system thermo object for simulation:
    H2O, N, CH4, N2O.
    
    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, ConveyorBelt, Fan, Polyethylene.

    Parameters
    ----------
    ins : Iterable(stream)
        condensed effluent from concentrator, filter press cakes.
    outs : Iterable(stream)
        solid cakes, fugitive N2O, fugitive CH4.
    moisture_content_out : float
        Desired moisture content of the Condensed effluent.

    Warnings
    --------
    Energy balance is not performed for this unit.

    References
    ----------
    [1] YEE et al. VOLUME REDUCTION SOLIDS TREATMENT SYSTEM. 
       https://patents.google.com/patent/WO2023288327A1/en?oq=WO2023288327A1

    See Also
    --------
    :class:`~.sanunits.BiogenicRefineryHHXdryer`
    :class:`~.sanunits.BiogenicRefineryHHX`
    '''
    
    _ins_size_is_fixed = False
    _N_outs = 5
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',ppl=None,
                 user_scale_up=1,extreme=False, capital_scale_func=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.user_scale_up = user_scale_up
        self.extreme = extreme
        self.capital_scale_func = capital_scale_func
        data = load_data(path=vr_drying_tunnel_path)
        if ppl is None:
            self.ppl = 6
        else:
            self.ppl = ppl
        
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids)) + ('OtherSS',)
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])

    def _run(self):
        mixture = WasteStream()
        mixture.mix_from(self.ins)
        # condensate, press_cakes = self.ins
        solid_cakes, N2O, CH4, NH3_gas, water_vapor  = self.outs
        solid_cakes.copy_like(mixture)
        solubles, solids = self.solubles, self.solids
        N2O.phase = CH4.phase = NH3_gas.phase = water_vapor.phase
        #Calculate water and solids in the solids cake
        mc_in = solid_cakes.imass['H2O'] / solid_cakes.F_mass # fraction
        mc_out = self.moisture_content_out/100
        if mc_in < mc_out*0.999:
            mc_out = mc_in
        TS_in = solid_cakes.imass[solids].sum()# kg TS dry/hr
        solid_cakes.imass['H2O'] = TS_in/(1-mc_out)*mc_out #kg water/hr
        
        #Calculate N2O and CH4 emissions
        CH4.imass['CH4'] = drying_CH4_to_air = \
            self.drying_CH4_emissions * self.carbon_COD_ratio * \
            solid_cakes.COD * solid_cakes.F_vol / 1000 # kg CH4 /hr
        drying_NH3_to_air = self.drying_NH3_emissions * solid_cakes.imass['NH3'] # kg NH3 /hr
        N2O.imass['N2O'] = drying_NH3_to_air * self.NH3_to_N2O/14/2*44 # kg N2O /hr
        solid_cakes.imass['NH3'] -= drying_NH3_to_air
        NH3_gas.imass['NH3'] = drying_NH3_to_air * (1-self.NH3_to_N2O) # kg NH3 /hr
        water_vapor.imass['H2O'] = (mixture.imass['H2O'] - solid_cakes.imass['H2O']) #kg H2O/hr
        # Store the calculated value of water_vapor.imass['H2O'] for use in the _cost function
        self.water_vapor_H2O = water_vapor.imass['H2O']  #kg H2O/hr
        
        #Calculate COD
        self.drying_CO2_to_air = (self.drying_CO2_emissions * self.carbon_COD_ratio
                             * solid_cakes.COD * solid_cakes.F_vol / 1000) # kg CO2 /hr
        # 44/12/16 are the molecular weights of CO2, C, and CH4, respectively
        solid_cakes.imass['sCOD'] -=  (self.drying_CO2_to_air/44*12+drying_CH4_to_air/16*12) / self.carbon_COD_ratio
    
    def _init_lca(self):
        self.construction = [
            Construction('stainless_steel',linked_unit=self, 
                         item='StainlessSteel', 
                         quantity_unit='kg'),
            Construction('conveyor_belt',linked_unit=self, 
                         item='ConveyorBelt', 
                         quantity_unit='m'),
            Construction('fan',linked_unit=self, 
                         item='Fan',
                         quantity_unit='kg'),
            Construction('polyethylene',linked_unit=self,
                         item='Polyethylene', 
                         quantity_unit='kg'),
            ]
        
    def _design(self):
        design=self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.stainless_steel_weight
        design['ConveyorBelt'] = constr[1].quantity = 0.07 #0.07 m equivalent to 1m based on 0.2m equivalent to 3m width in ecoinvent
        design['Fan'] = constr[2].quantity = 3 # a fan roughly weigh 3 kg
        design['Polyethylene'] = constr[3].quantity = self.polyethylene_weight
        
    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        #C['Stainless steel'] = self.stainless_steel_cost * D['StainlessSteel']
        C['Housing'] = self.housing_cost
        C['AirDuct'] = self.duct_cost + self.ventilator_cost
        C['Conveyor'] = self.conveyor_cost
        C['Misc.parts'] = self.miscellaneous_cost_ratio*(C['Housing'] +
                               C['AirDuct']+
                               C['Conveyor'])
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio
        if self.extreme:
            self.power_utility(self.water_vapor_H2O * self.extreme_energy_required_to_evaporize_water + 
                            self.conveyor_power_demand * self.conveyor_daily_operation*self.ppl/24) # kW
        else:
            self.power_utility(self.water_vapor_H2O * self.energy_required_to_evaporize_water + 
                            self.conveyor_power_demand * self.conveyor_daily_operation*self.ppl/24) # kW
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = ((self.capital_scale_func(total_equipment) if self.capital_scale_func else total_equipment)*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr

    def _calc_maintenance_labor_cost(self): #USD/hr
        if self.extreme:
            maintenance_labor_cost= (self.drying_tunnel_maintenance * self.wages * self.user_scale_up)*0.1
        else:
            maintenance_labor_cost= (self.drying_tunnel_maintenance * self.wages * self.user_scale_up)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        if self.extreme:
            return (self.water_vapor_H2O * self.extreme_energy_required_to_evaporize_water + 
                            self.conveyor_power_demand * self.conveyor_daily_operation*self.ppl/24) # kW
        else:
            return (self.water_vapor_H2O * self.energy_required_to_evaporize_water + 
                            self.conveyor_power_demand * self.conveyor_daily_operation*self.ppl/24) # kW
#%%
g2rt_housing_path = ospath.join(g2rt_su_data_path, '_g2rt_housing.csv')

@price_ratio()
class G2RTHousing(Copier):
    '''
    Housing of the generation II reinveted toilet.
    
    This is a non-reactive unit (i.e., the effluent is copied from the influent).
    
    The following impact items should be pre-constructed for life cycle assessment:
    Aluminum, Steel, ZincCoat, Polyethylene
    
    References
    ----------
    [1] Watabe et al. Advancing the Economic and Environmental Sustainability 
    of the NEWgenerator Nonsewered Sanitation System." ACS Environmental Au 
    3.4 (2023): 209-222.
    https://pubs.acs.org/doi/10.1021/acsenvironau.3c00001
    
    See Also
    --------
    :class:`~.sanunits.NEWgeneratorControls`
    '''
    def __init__(self, ID='', ins=None, outs=(),  thermo=None, init_with='WasteStream',
                 **kwargs):
        data = load_data(path=g2rt_housing_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        Copier.__init__(self, ID, ins, outs, thermo, init_with)
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    def _init_lca(self):
        self.construction = [
            Construction(item='Steel', linked_unit=self, quantity_unit='kg'),
            Construction(item='Aluminum', linked_unit=self, quantity_unit='kg'),
            Construction(item='ZincCoat', linked_unit=self, quantity_unit='m2'),
            Construction(item='Polyethylene', linked_unit=self, quantity_unit='kg'),
            ]
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = self.housing_steel_weight
        design['Aluminum'] = constr[1].quantity = self.housing_aluminum_weight
        design['ZincCoat'] = constr[2].quantity = self.housing_zinc_coat
        design['Polyethylene'] = constr[3].quantity = self.housing_PE_weight
        self.add_construction(add_cost=False)
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C["Housing_frame"] = self.system_housing
        C["Misc.parts"] = self.miscellaneous_cost_ratio* C["Housing_frame"]
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
    
    @property
    def OPEX(self):
        return 0.0 #USD/day

    @property
    def labor_expense(self):
        return 0.0 #USD/day

    @property
    def power_kW(self):
        '''[float] Power draw, [kW], auto-calculated from `power_utility`.'''
        return self.power_utility.rate

#%%
g2rt_controls_path = ospath.join(g2rt_su_data_path, '_g2rt_controls.csv')

@price_ratio()
class G2RTControls(Copier):
    '''
    Electronic control of the generation II reinveted toilet that interface with sensors,
    valves, pumps, and motors.
    
    This is a non-reactive unit (i.e., the effluent is copied from the influent).
    
    The following impact items should be pre-constructed for life cycle assessment:
    ControlUnits, Polycarbonate, Aluminum, ElectricCables, ElectronicsPassive, ElectronicsActive
    
    References
    ----------
    [1] Watabe et al. Advancing the Economic and Environmental Sustainability 
    of the NEWgenerator Nonsewered Sanitation System." ACS Environmental Au 
    3.4 (2023): 209-222.
    https://pubs.acs.org/doi/10.1021/acsenvironau.3c00001
    
    See Also
    --------
    :class:`~.sanunits.NEWgeneratorControls`
    '''
    def __init__(self, ID='', ins=None, outs=(),  thermo=None, init_with='WasteStream',
                 **kwargs):
        data = load_data(path=g2rt_controls_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        Copier.__init__(self, ID, ins, outs, thermo, init_with)

        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    def _init_lca(self):
        self.construction = [
            Construction(item='ControlUnits', linked_unit=self, quantity_unit='kg',
                         lifetime=self.control_system_PLC_lifetime),
            Construction(item='Polycarbonate', linked_unit=self, quantity_unit='kg'),
            Construction(item='Aluminum', linked_unit=self, quantity_unit='kg'),
            Construction(item='ElectricCables', linked_unit=self, quantity_unit='m'),
            Construction(item='ElectronicsPassive', linked_unit=self, quantity_unit='kg'),
            Construction(item='ElectronicsActive', linked_unit=self, quantity_unit='kg'),
            ]
    
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['ControlUnits'] = constr[0].quantity = self.control_control_units_weight
        design['Polycarbonate'] = constr[1].quantity = self.control_polycarbonate_weight
        design['Aluminum'] = constr[2].quantity = self.control_aluminum_weight
        design['ElectricCables'] = constr[3].quantity = self.control_cable_length
        design['ElectronicsPassive'] = constr[4].quantity = self.control_electronics_passive_weight
        design['ElectronicsActive'] = constr[5].quantity = self.control_electronics_active_weight
        self.add_construction(add_cost=False)
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C["Controller"] = self.programmable_logic_controller_cost * self.PLC_quantities
        C["IO_relay_modules"] = self.IO_relay_module_cost * self.relay_module_quantity
        C["Sensors"] = self.sensors_cost
        C["Cables"] = self.cables
        C["Misc.parts"] = self.miscellaneous_cost_ratio*(C["Controller"]+
                                                         C["IO_relay_modules"]+
                                                         C["Sensors"]+
                                                         C["Cables"]
                                                         )
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        power_demand = (
            self.control_system_energy_percycle * self.control_batch_cycle_perday +
            self.control_background_runtime_energy_perday
            )
        power_demand = power_demand / 24  # convert from kWh/d to kW
        self.power_utility(power_demand) # kW
        
        self.add_OPEX =  self._calc_replacement_cost() + self._calc_maintenance_labor_cost()
        
    def _calc_replacement_cost(self):
        control_system_replacement_cost = (
            self.PLC_quantities + self.sensors_cost / self.control_sensor_lifetime
            )
        return control_system_replacement_cost / (365 * 24) * self.price_ratio # USD/hr
    
    def _calc_maintenance_labor_cost(self):
        control_system_maintenance_labor = self.control_labor_replacement_misc_repairs * self.wages
        return control_system_maintenance_labor / (365 * 24) # USD/hr
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.control_system_ORP_energy_percycle * self.control_batch_cycle_perday +
            self.control_background_runtime_energy_perday)/24 #kW

#%%
g2rt_solids_separation_path = ospath.join(g2rt_su_data_path, '_g2rt_solids_separation.csv')

@price_ratio()
class G2RTSolidsSeparation(SanUnit):
    '''
    Solids separation unit in generation II reinveted toilets as a frontend separator [1].
    
    .. note:

    Non-reactive. Moisture content of the effluent solid is adjusted to be 99% [2].

    The following components should be included in system thermo object for simulation:
    H2O

    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, Polyethylene, ElectricMotor, Pump

    Parameters
    ----------
    ins : Iterable(stream)
        Influent stream.
    outs : Iterable(stream)
        liquid stream and solid stream after separation.
    moisture_content_out : float
        Moisture content of the effluent solids stream.

    References
    ----------
    [1] YEE et al. Buffer tank separation and homogenization system. 
    https://patents.google.com/patent/WO2023288114A1/en?oq=WO2023288114A1
    [2] YEE et al. Volume reduction non-sewered single unit toilet system.
    https://patents.google.com/patent/WO2023288326A1/en?oq=WO2023288326A1
    
    See Also
    ---------
    :class:`~.sanunits.BiogenicRefineryGrinder`
    '''
    _N_ins = 2
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', ppl=None,
                 user_scale_up=1, extreme = False, capital_scale_func=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.user_scale_up = user_scale_up
        self.extreme =  extreme
        self.capital_scale_func = capital_scale_func
        if ppl is None:
            self.ppl = 6
        else:
            self.ppl = ppl

        data = load_data(path=g2rt_solids_separation_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])

    def _init_lca(self):
        self.construction = [Construction("stainless_steel", linked_unit=self,
                                          item = "StainlessSteel", 
                                          quantity_unit= "kg"),
                             Construction("polyethylene", linked_unit=self,
                                          item = "Polyethylene",
                                          quantity_unit= "kg"),
                             Construction("electric_motor", linked_unit=self,
                                          item = "ElectricMotor",
                                          quantity_unit= "kg"),
                             Construction("pump", linked_unit=self,
                                          item = "Pump",
                                          quantity_unit= "ea"),
                             ]
    
    def _run(self):
        waste_in,flushing_water = self.ins
        # waste_in.mix_from(self.ins)
        # print(waste_in.imass['H2O']) #TODO:debug
        liquid_stream, solid_stream = self.outs
        solubles, solids = self.solubles, self.solids
        solid_stream.mix_from(self.ins)
        solid_stream.imass[solids] = solid_stream.imass[solids] * self.solids_separator_TSS_removal/100
        mc_in = solid_stream.imass['H2O'] / (waste_in.F_mass + flushing_water.F_mass) # fraction
        mc_out = self.moisture_content_out/100 #convert to fraction
        TS_in = waste_in.imass[solids].sum()+ flushing_water.imass[solids].sum()# kg TS dry/hr
        TS_out = solid_stream.imass[solids].sum()
        if mc_in < mc_out*0.999:
            mc_out = mc_in
        # # solid_stream.copy_flow(waste_in,solids) #all solids go to sludge, remove from waste_in
        # solid_stream.imass[solids] = solid_stream.imass[solids] * self.solids_separator_TSS_removal/100 #add the removed solids back

        # if mc_in < mc_out*0.999:
        #     print(f"Moisture content of the influent stream ({mc_in:.4f})"
        #           f"is smaller than that of the desired effluent stream ({mc_out:.4f})."
        #           "High solids and low flushing event detected, adding more flushing water!")
        #     solid_stream.imass['H2O'] = (waste_in.F_mass + flushing_water.F_mass) * mc_out
        #     # raise RuntimeError(f'Moisture content of the influent stream ({mc_in:.4f}) '
        #     #                    f'is smaller than the desired moisture content ({mc_out:.4f}).')
        
        # TS_in = waste_in.imass[solids].sum()+ flushing_water.imass[solids].sum()# kg TS dry/hr
        # TS_out = solid_stream.imass[solids].sum()
        
        #calculate water in the solid cakes
        solid_stream.imass['H2O'] = TS_out/(1-mc_out)*mc_out

        solid_stream.imass[solubles] = (waste_in.imass[solubles]+flushing_water.imass[solubles])*\
            (TS_out/(1-mc_out)-TS_out)/(waste_in.F_mass+ flushing_water.F_mass-TS_in)
        liquid_stream.mass = flushing_water.mass + waste_in.mass-solid_stream.mass
        
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.stainless_steel_weight
        design['Polyethylene'] = constr[1].quantity = self.polyethylene_weight
        design['ElectricMotor'] = constr[2].quantity = self.actuator_weight *(self.user_scale_up ** 0.6)
        design['Pump'] = constr[1].quantity = 1.*(self.user_scale_up ** 0.6)
        self.add_construction(add_cost = False)
    
    def _cost(self):
        C = self.baseline_purchase_costs
        C["Vacuum tank"] = (self.vacuum_tank_cost + 
                            self.separating_filter_cost +
                            self.solid_outlet_cost + 
                            self.liquid_outlet_cost+
                            self.inlet_chamber_cost)
        C["Vacuum pump"] = self.vacuum_pump_cost *(self.user_scale_up ** 0.6)
        C["Valves"] = self.valve_cost * self.valve_quantity
        C["Actuator"] = self.actuator_cost * self.actuator_quantity *(self.user_scale_up ** 0.6)
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        self.power_utility(self.vacuum_pump_power_demand * 
                           self.vacuum_pump_operation*
                           self.ppl/24) # kW
        
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = ((self.capital_scale_func(total_equipment) if self.capital_scale_func else total_equipment)*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr
    
    def _calc_maintenance_labor_cost(self): #USD/hr
        if self.extreme:
            maintenance_labor_cost= (self.solids_separator_maintenance * self.wages * self.user_scale_up)*0.1
        else:
            maintenance_labor_cost= (self.solids_separator_maintenance * self.wages * self.user_scale_up)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.vacuum_pump_power_demand * 
                           self.vacuum_pump_operation*
                           self.ppl/24) #kW

#%%
g2rt_belt_separation_path = ospath.join(g2rt_su_data_path, '_g2rt_belt_separation.csv')

@price_ratio()
class G2RTBeltSeparation(SanUnit):
    '''
    Belt separation unit in generation II reinveted toilets before the buffer tank [1].
    
    .. note:

    Non-reactive. Moisture content of the effluent solid is adjusted to be 96-99% [2].
    The solids in the liquid stream influent is partially transferred to the solids effluent.

    The following components should be included in system thermo object for simulation:
    H2O, OtherSS.

    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, ConveyorBelt, Polyethylene.

    Parameters
    ----------
    ins : Iterable(stream)
        liquid stream and solid stream from the solid separator, and retentate 
        from the ultrafiltration unit
    outs : Iterable(stream)
        liquid stream and solid stream are copied from the influent.
    moisture_content_out : float
        Moisture content of the effluent solids stream.

    References
    ----------
    [1] YEE et al. Buffer tank separation and homogenization system. 
    https://patents.google.com/patent/WO2023288114A1/en?oq=WO2023288114A1
    [2] YEE et al. Volume reduction non-sewered single unit toilet system.
    https://patents.google.com/patent/WO2023288326A1/en?oq=WO2023288326A1
    
    See Also
    ---------
    :class:`~.sanunits.BiogenicRefineryGrinder`
    '''
    _N_ins = 2
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', ppl=None,
                 user_scale_up = 1, moisture_content_out= None, extreme = False,
                 capital_scale_func=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.user_scale_up = user_scale_up
        self.extreme = extreme
        self.capital_scale_func = capital_scale_func
        if ppl is None:
            self.ppl = 6
        else:
            self.ppl = ppl

        data = load_data(path=g2rt_belt_separation_path)
        for para in data.index:
            if moisture_content_out is not None and para == "moisture_content_out":
                self.moisture_content_out = moisture_content_out
                continue  # Skip the moisture_content_out index if it's provided
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O']) 
        
 
    def _init_lca(self):
        self.construction = [
            Construction('stainless_steel',linked_unit=self, 
                         item='StainlessSteel', 
                         quantity_unit='kg'),
            Construction('conveyor_belt',linked_unit=self, 
                         item='ConveyorBelt', 
                         quantity_unit='m'),
            Construction('fan',linked_unit=self, 
                         item='Fan',
                         quantity_unit='kg'),
            Construction('polyethylene',linked_unit=self,
                         item='Polyethylene', 
                         quantity_unit='kg'),
            ]
    
    def _run(self):
        liquid_in, solid_in = self.ins
        # uf_retentate.show() #TODO:debug
        liquid_out, solid_out = self.outs
        solid_out.phase = 's'
        solubles, solids = self.solubles, self.solids
        solid_out.mix_from(self.ins)
        TS_in = solid_out.imass[solids].sum()
        #solid_stream.copy_flow(solid_stream_in,solids) #all solids in go to solids out
        solid_out.imass[solids] = solid_in.imass[solids] + liquid_in.imass[solids] * self.belt_separator_TSS_removal/100 
        # the removed solids
        TS_out = solid_out.imass[solids].sum() # kg TS dry/hr
        
        # mc_in = solid_stream.imass['H2O'] / (solid_stream_in.F_mass + 
        #                                      liquid_stream_in.F_mass +
        #                                      uf_retentate.F_mass) # fraction
        if self.extreme:
            mc_out = self.extreme_moisture_content_out/100 #convert to fraction; 
        else:
            mc_out = self.moisture_content_out/100 #convert to fraction; 
        mc_in = solid_out.imass['H2O'] / (liquid_in.F_mass + solid_in.F_mass)
        if mc_in < mc_out*0.999:
            mc_out = mc_in
        # if mc_in < mc_out*0.999:
        #     raise RuntimeError(f'Moisture content of the influent stream ({mc_in:.4f}) '
        #                        f'is smaller than the desired moisture content ({mc_out:.4f}).')
        # mix_in = WasteStream(ID='mix_in')
        # mix_in.mix_from((liquid_stream_in, uf_retentate,solid_stream_in))
        
        # TS_in = (solid_stream_in.imass[solids].sum()+ 
        #          liquid_stream_in.imass[solids].sum()+
        #          uf_retentate.imass[solids].sum()
        #          ) # kg TS dry/hr

        #calculate water and solid COD in the solid cakes
        solid_out.imass['H2O'] = TS_out/(1-mc_out)*mc_out

        # solid_stream._COD = (solid_stream_in.COD * solid_stream_in.F_vol +
        #                      (liquid_stream_in.COD * liquid_stream_in.F_vol 
        #                       + uf_retentate.COD * uf_retentate.F_vol) *
        #                      self.belt_separator_TSS_removal/100)/ solid_stream.F_vol

        solid_out.imass[solubles] = (liquid_in.imass[solubles] + 
                                        solid_in.imass[solubles]) *\
            (TS_out/(1-mc_out)-TS_out)/((solid_in.F_mass + liquid_in.F_mass)-TS_in)
        liquid_out.mass = solid_in.mass + liquid_in.mass -solid_out.mass

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.stainless_steel_weight
        design['ConveyorBelt'] = constr[1].quantity = 0.07 #~0.07m equivalent to 0.7m based on 0.3m equivalent to 3m width in ecoinvent
        design['Polyethylene'] = constr[2].quantity = self.polyethylene_weight
        self.add_construction(add_cost=False)
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C["Belt_separator"] = (self.belt_conveyor_cost *(self.user_scale_up ** 0.6)+
                               self.solid_inlet_chamber_cost+
                               self.liquid_inlet_chamber_cost+
                               self.housing_cost+
                               self.squeegee_cost+
                               self.splash_shield_cost
                               )
        C["Misc.parts"] = self.miscellaneous_cost_ratio * C["Belt_separator"]
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        
        power_demand = (
            self.belt_daily_operation * self.belt_power_demand * self.ppl/24
            ) # convert from kWh/d to kW
        self.power_utility(power_demand) # kW
                
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
        self.add_OPEX = ((self.capital_scale_func(total_equipment) if self.capital_scale_func else total_equipment)*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 5% of CAPEX per year
                         self._calc_maintenance_labor_cost()) #USD/hr

    def _calc_maintenance_labor_cost(self): #USD/hr
        if self.extreme:
            maintenance_labor_cost= (self.belt_separator_maintenance * self.wages *self.user_scale_up)*0.1
        else:
            maintenance_labor_cost= (self.belt_separator_maintenance * self.wages *self.user_scale_up)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        return (self.belt_daily_operation * self.belt_power_demand * self.ppl/24) #kW

#%%

ultrafiltration_path = ospath.join(g2rt_su_data_path, '_g2rt_ultrafiltration.csv')

@price_ratio()
class G2RTUltrafiltration(SanUnit):
    '''
    Ultrafiltration in the generation II reinveted toilets is used for removing suspended solids
    with automated backwash.
    
    Modified from ultrafiltration unit in Duke Reclaimer system.

    The following impact items should be pre-constructed for life cycle assessment:
    GFRPlastic, Steel, NylonGlassFilled, CastingBrass

    Parameters
    ----------
    ins : Iterable(stream)
        waste: liquid waste stream to be treated by ultrafiltration unit.
    outs : Iterable(stream)
        treated: treated liquid leaving ultrafiltration unit.
        retentate: concentrated retentate leaving ultrafiltration unit.
    reject_recycle_ratio: float, optional, 0 to 1
        fraction of rejected solid streams that recirculates to belt separation. Defaults to 1.
    ppl: int
        Total number of users for scaling of costs.

    References
    ----------
    [1] Trotochaud et al., Laboratory Demonstration and Preliminary Techno-Economic Analysis of an Onsite
    Wastewater Treatment System Environ. Sci. Technol. 2020, 54, (24), 16147–16155.
    https://dx.doi.org/10.1021/acs.est.0c02755
    
    [2] Duke Center for WaSH-AID Reclaimer design team data and guidance
    https://washaid.pratt.duke.edu/work/water-sanitation/reinvent-toilet-challenge
    
    See Also
    ---------
    :class:`~.sanunits.ReclaimerUltrafiltration`
    
    '''
    _N_ins = 1
    _N_outs = 3

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 reject_recycle_ratio = 1,user_scale_up=1, extreme = False,
                 **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.reject_recycle_ratio = reject_recycle_ratio
        self.user_scale_up = user_scale_up
        self.extreme = extreme
        
        data = load_data(path=ultrafiltration_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components

        self.solids = tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])
        
    def _init_lca(self):
        self.construction = [
            Construction(item='GFRPlastic', linked_unit=self, quantity_unit='kg'),
            Construction(item='Steel', linked_unit=self, quantity_unit='kg'),
            Construction(item='NylonGlassFilled', linked_unit=self, quantity_unit='kg'),
            Construction(item='CastingBrass', linked_unit=self, quantity_unit='kg'),
            ]

    def _run(self):
        waste_in = self.ins[0]
        liquid_stream, recycled_solid_stream, discharged_solid_stream = self.outs
        solid_stream = WasteStream()
        solubles, solids = self.solubles, self.solids
        solid_stream.copy_flow(waste_in,solids) #all solids go to sludge, remove from waste_in
        liquid_stream.copy_flow(waste_in, solubles) #only copy solubles 
        # TS_in = waste_in.imass[solids].sum() # kg TS dry/hr
        # TS_out = solid_stream.imass[solids].sum()
        solid_stream.imass['H2O'] = waste_in.imass['H2O']*(1-self.water_recovery_rate/100)
        solid_stream.imass[solids] = waste_in.imass[solids] * self.TSS_removal/100 # the removed solids
        solid_stream.imass[solubles] = waste_in.imass[solubles]*\
            solid_stream.imass['H2O']/waste_in.imass['H2O']
        solid_stream.imass['sCOD'] += waste_in.imass['sCOD']*self.sCOD_removal/100
        liquid_stream.mass = waste_in.mass-solid_stream.mass
        solid_stream.split_to(recycled_solid_stream, discharged_solid_stream, self.reject_recycle_ratio)
        # print(f"The recycled UF water flow is {solid_stream.imass['H2O']} kg/h.")

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['GFRPlastic'] = constr[0].quantity = self.Plastic_weight
        design['Steel'] = constr[1].quantity = self.Steel_weight
        design['NylonGlassFilled'] = constr[2].quantity = self.Nylon_weight
        design['CastingBrass'] = constr[3].quantity = self.Brass_weight
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs 
        C['Pipes'] = self.one_in_pipe_SCH40 + self.onehalf_in_pipe_SCH40 + self.three_in_pipe_SCH80
        C['Fittings'] = (
            self.one_in_elbow_SCH80 +
            self.one_in_tee_SCH80 +
            self.one_in_SCH80 +
            self.one_onehalf_in_SCH80 +
            self.onehalf_in_SCH80 +
            self.three_in_SCH80_endcap +
            self.one_one_NB_MTA +
            self.one_onehalf_NB_MTA +
            self.foot_valve +
            self.one_onehalf_in_SCH80_threadedtee +
            self.three_in_pipe_clamp +
            self.one_in_pipe_clamp +
            self.onehalf_in_pipe_clamp +
            self.two_way_valve +
            self.UF_brush
            )
        C['UF_unit'] = self.UF_unit
        C['Backwash_parts'] = (self.air_scrubbing_blower_cost + self.back_flush_pump_cost)*(self.user_scale_up ** 0.6)
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio/3 #scaled down from 30 users to 6 users per day

        self.add_OPEX = self._calc_replacement_cost()
        # [W][1 kW/1000 W][hr/d][1 d/ 24 h] = [kW]
        self.back_wash_cycles = ceil(self.ins[0].get_TSS()*self.ins[0].F_vol/1000*24/self.accumulated_TSS_triggering_cleaning)
        if self.extreme:
            power_demand = (self.extreme_ultrafiltration_energy_consumption * self.outs[0].F_mass+
                        (self.backwash_pump_energy_percycle+
                         self.air_scrubbing_energy_percycle)*self.back_wash_cycles/24)  #kW
        else:
            power_demand = (self.ultrafiltration_energy_consumption * self.outs[0].F_mass+
                        (self.backwash_pump_energy_percycle+
                         self.air_scrubbing_energy_percycle)*self.back_wash_cycles/24)  #kW
        self.power_utility(power_demand) #scaled down from 30 users to 6 users per day

    def _calc_replacement_cost(self):
        pipe_replacement_cost = (
            self.one_in_pipe_SCH40 / self.one_in_pipe_SCH40_lifetime +
            self.onehalf_in_pipe_SCH40 / self.onehalf_in_pipe_SCH40_lifetime +
            self.three_in_pipe_SCH80 / self.three_in_pipe_SCH80_lifetime
            )
        
        fittings_replacement_cost = (
            self.one_in_elbow_SCH80 / self.one_in_elbow_SCH80_lifetime +
            self.one_in_tee_SCH80 / self.one_in_tee_SCH80_lifetime +
            self.one_in_SCH80 / self.one_in_SCH80_lifetime +
            self.one_onehalf_in_SCH80 / self.one_onehalf_in_SCH80_lifetime +
            self.onehalf_in_SCH80 / self.onehalf_in_SCH80_lifetime +
            self.three_in_SCH80_endcap / self.three_in_SCH80_endcap_lifetime +
            self.one_one_NB_MTA / self.one_one_NB_MTA_lifetime +
            self.one_onehalf_NB_MTA / self.one_onehalf_NB_MTA_lifetime +
            self.foot_valve / self.foot_valve_lifetime +
            self.one_onehalf_in_SCH80_threadedtee / self.one_onehalf_in_SCH80_threadedtee_lifetime +
            self.three_in_pipe_clamp / self.three_in_pipe_clamp_lifetime +
            self.one_in_pipe_clamp / self.one_in_pipe_clamp_lifetime +
            self.onehalf_in_pipe_clamp / self.onehalf_in_pipe_clamp_lifetime +
            self.two_way_valve / self.two_way_valve_lifetime +
            self.UF_brush / self.UF_brush_lifetime
            )
        if self.extreme:
            uf_replacement_cost = self.UF_unit/ self.extreme_UF_unit_lifetime
            other_replacement = (self.air_scrubbing_blower_cost + self.back_flush_pump_cost)/self.extreme_UF_unit_lifetime*(self.user_scale_up ** 0.6)
        else:
            uf_replacement_cost = self.UF_unit/ self.UF_unit_lifetime
            other_replacement = (self.air_scrubbing_blower_cost + self.back_flush_pump_cost)/self.UF_unit_lifetime*(self.user_scale_up ** 0.6)

        total_replacement_cost = self.price_ratio * (pipe_replacement_cost + 
                                                     fittings_replacement_cost + 
                                                     uf_replacement_cost+ other_replacement)/3 # USD/year, linear scaled to 6 users
        return total_replacement_cost / (365 * 24)  # USD/hr
    
    def _calc_maintenance_labor_cost(self): #USD/hr
        if self.extreme:
            maintenance_labor_cost= (self.ultrafiltration_maintenance * self.wages *self.user_scale_up)*0.1
        else:
            maintenance_labor_cost= (self.ultrafiltration_maintenance * self.wages *self.user_scale_up)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        if self.extreme:
            return (self.extreme_ultrafiltration_energy_consumption * self.outs[0].F_mass+
                        (self.backwash_pump_energy_percycle+
                         self.air_scrubbing_energy_percycle)*self.back_wash_cycles/24)  #kW
        else:
            return (self.ultrafiltration_energy_consumption * self.outs[0].F_mass+
                        (self.backwash_pump_energy_percycle+
                         self.air_scrubbing_energy_percycle)*self.back_wash_cycles/24)  #kW

#%%
reverse_osmosis_path = ospath.join(g2rt_su_data_path, '_g2rt_reverse_osmosis.csv')

@price_ratio()
class G2RTReverseOsmosis(SanUnit):
    '''
    A reverse osmosis unit process to recover water and concentrate liquid stream.
    The model is based on a fraction of water recovered.
    
    Parameters
    ----------
    ins : 
        Inlet fluid to be split.
    outs : 
        * [0] Permeate
        * [1] Brine
    water_recovery : float, optional, 0 to 1
        Water recovered to 0th stream. Defaults to 0.6
    TDS_removal: float, optional, 0 to 1
        rejection rate to total dissolved salts. Defaults to 0.95
    permeate_recycle_ratio: float, optional, 0 to 1
        fraction of permeate that recirculates to front end as flushing water. Defaults to 1.
        
    The following impact items should be pre-constructed for life cycle assessment:
    GFRPlastic, Steel, ReverseOsmosisModule
    '''
    _N_ins = 1
    _N_outs = 3

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 water_recovery=0.6,TDS_removal = 0.95, permeate_recycle_ratio =1,
                 user_scale_up=1,extreme=False, capital_scale_func=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.water_recovery = water_recovery
        self.TDS_removal = TDS_removal
        self.permeate_recycle_ratio = permeate_recycle_ratio
        self.user_scale_up = user_scale_up
        self.extreme = extreme
        self.capital_scale_func = capital_scale_func

        data = load_data(path=reverse_osmosis_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        cmps = self.components
        self.solids = tuple((cmp.ID for cmp in cmps.solids))
        self.solubles = tuple([i.ID for i in cmps if i.ID not in self.solids and i.ID != 'H2O'])

    
    def _run(self):
        waste_in = self.ins[0]
        permeate_recycle, brine, permeate_discharge = self.outs
        solubles, solids = self.solubles, self.solids
        permeate = WasteStream()
        for i in solubles:
            permeate.imass[i] = waste_in.imass[i]

        # permeate.copy_flow(waste_in, solubles) #copy_flow removes stream_impact_item 
        brine.copy_like(waste_in) #copy all the solids and solubles
        if self.extreme:
            permeate.imass['H2O'] = waste_in.imass['H2O'] * self.extreme_water_recovery
        else:
            permeate.imass['H2O'] = waste_in.imass['H2O'] * self.water_recovery
        permeate.imass[solubles] = waste_in.imass[solubles]*(1-self.TDS_removal)
        permeate.imass[solids] = waste_in.imass[solids]*(1-self.TSS_removal)
        permeate.imass['NH3'] = waste_in.imass['NH3']*(1-self.NH3_removal)
        permeate.split_to(permeate_recycle, permeate_discharge, self.permeate_recycle_ratio)
        # permeate.imass['P'] = waste_in.imass['P'] * 0.01 #higher rejection for P
        brine.mass = waste_in.mass - permeate.mass

    def _init_lca(self):
        self.construction = [
            Construction(item='GFRPlastic', linked_unit=self, quantity_unit='kg'),
            Construction(item='Steel', linked_unit=self, quantity_unit='kg'),
            Construction(item='ReverseOsmosisModule', linked_unit=self, quantity_unit='m2'),
            ]
        
    def _design(self):
        design = self.design_results
        constr = self.construction
        design['GFRPlastic'] = constr[0].quantity = self.GFRPlastic_weight
        design['Steel'] = constr[1].quantity = self.Steel_weight
        design['ReverseOsmosisModule'] = constr[2].quantity = self.membrane_area *(self.user_scale_up ** 0.6)
        self.add_construction(add_cost=False)
        
    def _cost(self):
        C = self.baseline_purchase_costs
        C['Pipes'] = self.piping_cost
        C['RO_system'] = self.reverse_osmosis_system_cost
        C['Backwash_parts'] = (self.air_scrubbing_blower_cost + self.back_flush_pump_cost)*(self.user_scale_up ** 0.6)
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
        self.back_wash_cycles = ceil(self.ins[0].COD*self.ins[0].F_vol/1000*24/
                                     self.accumulated_COD_triggering_cleaning/self.membrane_area)
        if self.extreme:
            power_demand = (self.extreme_RO_energy_consumption * (self.outs[0].F_mass+self.outs[2].F_mass)+
                            (self.backwash_pump_energy_percycle+
                             self.air_scrubbing_energy_percycle)*self.back_wash_cycles/24)
        else:
            power_demand = (self.RO_energy_consumption * (self.outs[0].F_mass+self.outs[2].F_mass)+
                            (self.backwash_pump_energy_percycle+
                             self.air_scrubbing_energy_percycle)*self.back_wash_cycles/24)

        self.power_utility(power_demand) #kW
        
        total_equipment = 0.
        for cost in C.values():
           total_equipment += cost
                  
        self.add_OPEX = (((self.capital_scale_func(total_equipment) if self.capital_scale_func else total_equipment)*self.material_replacement_cost/(365*24) + 
                         #USD/hr, assume replacement cost 4% of CAPEX per year
                         self._calc_membrane_replacement_cost() +
                         self._calc_maintenance_labor_cost())) #USD/hr
    def _calc_membrane_replacement_cost(self): #USD/hr
 
        membrane_replacement_cost = self.membrane_cost / self.membrane_life_time *(self.user_scale_up ** 0.6)  #USD/yr
        return membrane_replacement_cost/(365*24) #USD/hr
    
    def _calc_maintenance_labor_cost(self): #USD/hr
        if self.extreme:
            maintenance_labor_cost= (self.reverse_osmosis_maintenance * self.wages *self.user_scale_up)*0.1
        else:
            maintenance_labor_cost= (self.reverse_osmosis_maintenance * self.wages *self.user_scale_up)
        return maintenance_labor_cost / (365*24)
    
    @property
    def OPEX(self):
        return (self.add_OPEX['Additional OPEX']-self._calc_maintenance_labor_cost())*24 #USD/day
    
    @property
    def labor_expense(self):
        return self._calc_maintenance_labor_cost()*24 #USD/day
    
    @property
    def power_kW(self):
        if self.extreme:
            return (self.extreme_RO_energy_consumption * (self.outs[0].F_mass+self.outs[2].F_mass)+
                            (self.backwash_pump_energy_percycle+
                             self.air_scrubbing_energy_percycle)*self.back_wash_cycles/24)
        else:
            return (self.RO_energy_consumption * (self.outs[0].F_mass+self.outs[2].F_mass)+
                            (self.backwash_pump_energy_percycle+
                             self.air_scrubbing_energy_percycle)*self.back_wash_cycles/24)
        