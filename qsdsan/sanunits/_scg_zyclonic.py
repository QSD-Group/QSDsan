#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Hannah Lohman <hlohman94@gmail.com>
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module contains unit operations used in the SCG Zyclonic system
developed by SCG (Siam Cement Group) Chemicals (using Aquonic 1000 L).
https://products.scgchemicals.com/en/products-services/technology-solutions/reinvented-toilet-total-solution

For units without replacement cost calculation,
this is because all of their parts' lifetime equal the system lifetime
(i.e., they do not need to be replaced).

For units without maintenance labor cost calculation, this is because it
is accounted for the entire system and not broken down by the unit
(the labor cost is included in the `SCGZyclonicAquonic1000` unit).

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from .. import SanUnit, Construction
from ..utils import ospath, load_data, data_path, price_ratio


__all__ = (
    'SCGZyclonicAquonic1000',
    'SCGZyclonicControls',
    'SCGZyclonicEqualizerTank',
    'SCGZyclonicPiping',
    'SCGZyclonicSolar',
    'SCGZyclonicTreatedWaterTank',
    )

# sz_data_path = ospath.join(data_path, 'sanunit_data/sz')
sz_data_path = '/Users/saumitrarai/Desktop/Research/EXPOsan-private/exposan/scg_aquonic/data/'


# %%

aquonic_path = ospath.join(sz_data_path, '_sz_aquonic1000.csv')

class SCGZyclonicAquonic1000(SanUnit):
    '''
    Treatment of liquid waste in SCG Aquonic 1000L using 9 chambers:
    1. up-flow filtration;
    2. sedimentation;
    3. anaerobic;
    4. aerobic;
    5. anoxic;
    6. recirculation;
    7. sedimentation;
    8. electrochemical disinfection;
    9. treated water storage

    .. note::

        This is currently a non-reactive unit
        (i.e., the effluent is copied from the influent)
        due to the lack of treatment data.

    The following impact items should be pre-constructed for life cycle assessment:
    PVC, Polyethylene, Polypropylene, SoilMedia, Titanium, Pump, Concrete,
    ElectricCables, Silicone, InjectionMolding.

    Parameters
    ----------
    if_gridtied : Bool
        Set to True if system is using grid electricity and False if system is using photovoltaic electricity

    References
    ----------
    Communication with the design team, communication files include:

        [1] Illinois_Requested_Info_Pack.pdf

        [2] Manual of AQ Tank Installation (16-Apr-2021).pdf

        [3] SCG Follow-up Questions 11-15-2021.docx

        [4] UIUC_info_from_AIT.docx
    '''
    _N_ins = 1  # liquid waste
    _N_outs = 1  # treated water

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_gridtied=True, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.if_gridtied = if_gridtied

        data = load_data(path=aquonic_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _init_lca(self):
        self.construction = [
            Construction(item='PVC', linked_unit=self, quantity_unit='kg'),
            Construction(item='Polyethylene', linked_unit=self, quantity_unit='kg'),
            Construction(item='Polypropylene', linked_unit=self, quantity_unit='kg'),
            Construction(item='SoilMedia', linked_unit=self, quantity_unit='kg'),
            
            Construction(item='Titanium', linked_unit=self, quantity_unit='kg'),
            # Construction(item='ChlorineTablets', linked_unit=self, quantity_unit='kg'),
            
            Construction(item='Pump', linked_unit=self, quantity_unit='ea'),
            Construction(item='Concrete', linked_unit=self, quantity_unit='m3'),
            
            Construction(item='ElectricCables', linked_unit=self, quantity_unit='m'),
            
            Construction(item='Silicone', linked_unit=self, quantity_unit='kg'),
            Construction(item='InjectionMolding', linked_unit=self, quantity_unit='kg'),
            ]

    def _design(self):  # Material quantities for LCA calculations
        design = self.design_results
        constr = self.construction
        design['PVC'] = constr[0].quantity = (
            self.qty_1inch_PVC_pipe * self.mass_1inch_PVC_pipe +
            self.qty_1inch_PVC_valve * self.mass_1inch_PVC_valve +
            self.qty_1inch_PVC_90degree_elbow * self.mass_1inch_PVC_90degree_elbow +
            self.qty_1inch_PVC_T_elbow * self.mass_1inch_PVC_T_elbow +
            self.qty_1inch_PVC_3way_corner_elbow * self.mass_1inch_PVC_3way_corner_elbow
            )
        design['Polyethylene'] = constr[1].quantity = self.qty_PE_housing * self.mass_PE_housing
        design['Polypropylene'] = constr[2].quantity = self.qty_plastic_media
        design['SoilMedia'] = constr[3].quantity = self.qty_soil_media
        
        design['Titanium'] = constr[4].quantity = self.qty_electrode_plates * self.mass_electrode_plates
        # design['ChlorineTablets'] = constr[4].quantity = self.qty_Cl_tablets * self.mass_Cl_tablets # m3 * (kg/m3)
        
        design['Pump'] = constr[5].quantity = self.qty_aerator + self.qty_pump
        design['Concrete'] = constr[6].quantity = self.qty_concrete
        
        design['ElectricCables'] = constr[7].quantity = self.qty_electrical_cable
        
        design['Silicone'] = constr[7].quantity = self.qty_air_hose * self.mass_air_hose  # air hose tubing material
        design['InjectionMolding'] = constr[8].quantity = self.qty_PE_housing * self.mass_PE_housing
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['PipesAndValves'] = (
            self.qty_1inch_PVC_pipe * self.price_1inch_PVC_pipe +
            self.qty_1inch_PVC_valve * self.price_1inch_PVC_valve +
            self.qty_1inch_PVC_90degree_elbow * self.price_1inch_PVC_90degree_elbow +
            self.qty_1inch_PVC_T_elbow * self.price_1inch_PVC_T_elbow +
            self.qty_1inch_PVC_3way_corner_elbow * self.price_1inch_PVC_3way_corner_elbow
            )
        C['AerationAndPumps'] = (
            self.qty_aerator * self.price_aerator + \
            self.qty_air_hose * self.price_air_hose + self.qty_pump * self.price_pump \
            + self.qty_electrical_cable * self.price_electrical_cable
            )
        C['Housing'] = (
            self.qty_PE_housing * self.mass_PE_housing * self.price_PE_housing +
            self.qty_concrete * self.price_concrete
            )
        C['Media'] = (
            self.qty_plastic_media * self.price_plastic_media +
            self.qty_soil_media * self.price_soil_media
            )
        
        C['ECPlates'] = self.qty_electrode_plates * self.price_electrode_plates
        # C['ChlorineTablets'] = self.qty_Cl_tablets * self.mass_Cl_tablets * self.price_Cl_tablets # (m3 * kg/m3 * USD/kg)

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

        if self.if_gridtied:
            power_demand = self.aquonic1000_energy_demand / 24  # kWh/day to # kW
        else:
            power_demand = 0

        self.power_utility(power_demand)  # kW

    def _calc_replacement_cost(self):  # O&M material costs (USD/hour)
        pump_replacement_cost = self.qty_pump * self.price_pump / self.lifetime_pump  # (USD/year)
        soil_media_replacement_cost = self.qty_soil_media * self.price_soil_media / self.lifetime_soil_media  # (USD/year)
        
        ec_plate_replacement_cost = self.qty_electrode_plates * self.price_electrode_plates  # (USD/year)
        # Cl_plate_replacement_cost = self.qty_Cl_tablets * self.mass_Cl_tablets * self.price_Cl_tablets / self.lifetime_Cl_tablets
        
        aquonic1000_replacement_cost = pump_replacement_cost + soil_media_replacement_cost + ec_plate_replacement_cost  # (USD/year)
        # aquonic1000_replacement_cost = pump_replacement_cost + soil_media_replacement_cost + Cl_plate_replacement_cost
        
        aquonic1000_replacement_cost = aquonic1000_replacement_cost / (365 * 24)  # convert from USD/year to USD/hour
        return aquonic1000_replacement_cost

    def _calc_maintenance_labor_cost(self):  # O&M labor costs (USD/hour)
        maintenance_labor_cost = (self.labor_hours_backwash+self.labor_hours_ec_plates) * self.labor_wage # (USD/year)
        # maintenance_labor_cost = (self.labor_hours_backwash+self.labor_hours_Cl_tablets) * self.labor_wage # (USD/year)
        
        maintenance_labor_cost = maintenance_labor_cost / (365 * 24)  # convert from USD/year to USD/hour
        return maintenance_labor_cost


# %%

controls_path = ospath.join(sz_data_path, '_sz_controls.csv')

@price_ratio()
class SCGZyclonicControls(SanUnit):
    '''
    System control box for SCG Zyclonic systems.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    The following impact items should be pre-constructed for life cycle assessment:
    Switches, ElectronicsActive, ElectronicsPassive.

    References
    ----------
    Communication with the design team, communication files include:

        [1] Illinois_Requested_Info_Pack.pdf

        [2] Manual of AQ Tank Installation (16-Apr-2021).pdf

        [3] SCG Follow-up Questions 11-15-2021.docx

        [4] UIUC_info_from_AIT.docx
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        data = load_data(path=controls_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    def _init_lca(self):
        self.construction = [
            Construction(item='Switches', linked_unit=self, quantity_unit='kg'),
            Construction(item='ElectronicsActive', linked_unit=self, quantity_unit='kg'),
            Construction(item='ElectronicsPassive', linked_unit=self, quantity_unit='kg'),
            ]

    def _design(self):  # Material quantities for LCA calculations
        design = self.design_results
        constr = self.construction
        design['Switches'] = constr[0].quantity = (
            self.qty_breaker_1P_16A * self.mass_breaker_1P_16A +
            self.qty_breaker_1P_2A * self.mass_breaker_1P_2A +
            self.qty_breaker_1P_3A * self.mass_breaker_1P_3A +
            self.qty_breaker_2P_25A * self.mass_breaker_2P_25A +
            self.qty_breaker_2P_3A * self.mass_breaker_2P_3A
            )
        design['ElectronicsActive'] = constr[1].quantity = (
            self.qty_inverter * self.mass_inverter +
            self.qty_IO_card * self.mass_IO_card +
            self.qty_level_control_relay * self.mass_level_control_relay +
            self.qty_PLC_controller * self.mass_PLC_controller +
            self.qty_power_supply * self.mass_power_supply +
            self.qty_temp_controller * self.mass_temp_controller +
            self.qty_thermostat * self.mass_thermostat
            )
        design['ElectronicsPassive'] = constr[2].quantity = (
            self.qty_emergency_button * self.mass_emergency_button +
            self.qty_fuse * self.mass_fuse +
            self.qty_fuse_holder * self.mass_fuse_holder +
            self.qty_polot_lamp * self.mass_polot_lamp +
            self.qty_push_button * self.mass_push_button +
            self.qty_push_button_with_lamp * self.mass_push_button_with_lamp +
            self.qty_relay_2_NO * self.mass_relay_2_NO +
            self.qty_relay_4_NO * self.mass_relay_4_NO +
            self.qty_relay_socket * self.mass_relay_socket +
            self.qty_selector_swt_2way * self.mass_selector_swt_2way +
            self.qty_selector_swt_3way * self.mass_selector_swt_3way +
            self.qty_solid_state_relay * self.mass_solid_state_relay
            )
        self.add_construction(add_cost=False)

    def _cost(self):
        self.baseline_purchase_costs['ControlBox'] = (
            self.qty_breaker_1P_16A * self.price_breaker_1P_16A +
            self.qty_breaker_1P_2A * self.price_breaker_1P_2A +
            self.qty_breaker_1P_3A * self.price_breaker_1P_3A +
            self.qty_breaker_2P_25A * self.price_breaker_2P_25A +
            self.qty_breaker_2P_3A * self.price_breaker_2P_3A +
            self.qty_inverter * self.price_inverter +
            self.qty_IO_card * self.price_IO_card +
            self.qty_level_control_relay * self.price_level_control_relay +
            self.qty_PLC_controller * self.price_PLC_controller +
            self.qty_power_supply * self.price_power_supply +
            self.qty_temp_controller * self.price_temp_controller +
            self.qty_thermostat * self.price_thermostat +
            self.qty_emergency_button * self.price_emergency_button +
            self.qty_fuse * self.price_fuse +
            self.qty_fuse_holder * self.price_fuse_holder +
            self.qty_polot_lamp * self.price_polot_lamp +
            self.qty_push_button * self.price_push_button +
            self.qty_push_button_with_lamp * self.price_push_button_with_lamp +
            self.qty_relay_2_NO * self.price_relay_2_NO +
            self.qty_relay_4_NO * self.price_relay_4_NO +
            self.qty_relay_socket * self.price_relay_socket +
            self.qty_selector_swt_2way * self.price_selector_swt_2way +
            self.qty_selector_swt_3way * self.price_selector_swt_3way +
            self.qty_solid_state_relay * self.price_solid_state_relay
            ) * self.price_ratio


# %%

equalizer_tank_path = ospath.join(sz_data_path, '_sz_equalizer_tank.csv')

@price_ratio()
class SCGZyclonicEqualizerTank(SanUnit):
    '''
    Equalizer tank (1200 L volume) for the SCG Zyclonic system
    is used to hold untreated water before treatment steps.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    The following impact items should be pre-constructed for life cycle assessment:
    Polyethylene, InjectionMolding.

    References
    ----------
    Communication with the design team, communication files include:

        [1] Illinois_Requested_Info_Pack.pdf

        [2] Manual of AQ Tank Installation (16-Apr-2021).pdf

        [3] SCG Follow-up Questions 11-15-2021.docx

        [4] UIUC_info_from_AIT.docx
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        data = load_data(path=equalizer_tank_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    def _init_lca(self):
        self.construction = [
            Construction(item='Polyethylene', linked_unit=self, quantity_unit='kg'),
            Construction(item='InjectionMolding', linked_unit=self, quantity_unit='kg'),
            ]

    def _design(self):  # Material quantities for LCA calculations
        design = self.design_results
        constr = self.construction
        design['Polyethylene'] = constr[0].quantity = self.mass_tank
        design['InjectionMolding'] = constr[1].quantity = self.mass_tank
        self.add_construction(add_cost=False)

    def _cost(self):
        self.baseline_purchase_costs['EqualizerTank'] = self.price_tank * self.price_ratio


# %%

piping_path = ospath.join(sz_data_path, '_sz_piping.csv')

@price_ratio()
class SCGZyclonicPiping(SanUnit):
    '''
    System pipeline and transmission pumps connecting all treatment processes
    for the SCG Zyclonic system.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    The following impact items should be pre-constructed for life cycle assessment:
    PVC, Pump.

    References
    ----------
    Communication with the design team, communication files include:

        [1] Illinois_Requested_Info_Pack.pdf

        [2] Manual of AQ Tank Installation (16-Apr-2021).pdf

        [3] SCG Follow-up Questions 11-15-2021.docx

        [4] UIUC_info_from_AIT.docx
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        data = load_data(path=piping_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    def _init_lca(self):
        self.construction = [
            Construction(item='PVC', linked_unit=self, quantity_unit='kg'),
            Construction(item='Pump', linked_unit=self, quantity_unit='ea'),
            ]

    def _design(self):  # Material quantities for LCA calculations
        design = self.design_results
        constr = self.construction
        design['PVC'] = constr[0].quantity = (
            self.qty_1inch_PVC_pipe * self.mass_1inch_PVC_pipe +
            self.qty_4inch_PVC_pipe * self.mass_4inch_PVC_pipe
            )
        design['Pump'] = constr[1].quantity = self.qty_pump
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['Pipeline'] = (
            self.qty_1inch_PVC_pipe * self.price_1inch_PVC_pipe +
            self.qty_4inch_PVC_pipe * self.price_4inch_PVC_pipe
            )
        C['Pump'] = self.qty_pump * self.price_pump * self.price_ratio


    def _calc_replacement_cost(self):  # O&M material costs (USD/hour)
        pump_replacement_cost = self.qty_pump * self.price_pump / self.lifetime_pump  # (USD/year)
        system_pipeline_replacement_cost = pump_replacement_cost / (365 * 24) * self.price_ratio  # convert from USD/year to USD/hour
        return system_pipeline_replacement_cost


# %%

solar_path = ospath.join(sz_data_path, '_sz_solar.csv')

@price_ratio()
class SCGZyclonicSolar(SanUnit):
    '''
    A 5-6 kWh/d photovoltaic (solar) energy system that can provide power
    for the SCG Zyclonic system.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    The following impact items should be pre-constructed for life cycle assessment:
    Solar, Battery.

    References
    ----------
    [1] Based on the Eco-San solar energy system but modified for the SCG Zyclonic system.

    See Also
    --------
    :class:`qsdsan.sanunits.EcoSanSolar`
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        data = load_data(path=solar_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    def _init_lca(self):
        self.construction = [
            Construction(item='Solar', linked_unit=self, quantity_unit='m2'),
            Construction(item='Battery', linked_unit=self, quantity_unit='ea'),
            ]

    def _design(self):  # Material quantities for LCA calculations
        design = self.design_results
        constr = self.construction
        design['Solar'] = constr[0].quantity = self.qty_solar_panel * self.area_solar_panel
        design['Battery'] = constr[1].quantity = self.qty_battery * self.mass_battery
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['SolarPanel'] = self.power_demand_W * self.price_solar_panel
        C['SolarModule'] = self.qty_solar_module * self.price_solar_module
        C['Battery'] = self.qty_battery * self.price_battery
        C['BatteryHolder'] = self.qty_battery_holder * self.price_battery_holder
        C['Inverter'] = self.qty_inverter * self.price_inverter

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

    def _calc_replacement_cost(self):  # O&M material costs (USD/hour)
        solar_panel_replacement_cost = self.power_demand_W * self.price_solar_panel / self.lifetime_solar_panel  # (USD/year)
        battery_replacement_cost = self.qty_battery * self.price_battery / self.lifetime_battery  # (USD/year)
        solar_replacement_cost = (solar_panel_replacement_cost + battery_replacement_cost) / (365 * 24) * self.price_ratio  # convert from USD/year to USD/hour
        return solar_replacement_cost

    def _calc_maintenance_labor_cost(self):
        maintenance_labor_cost = self.labor_hours_pannel_cleaning * self.labor_wage # (USD/year)
        maintenance_labor_cost = maintenance_labor_cost / (365 * 24)  # convert from USD/year to USD/hour
        return maintenance_labor_cost


# %%

treated_water_tank_path = ospath.join(sz_data_path, '_sz_treated_water_tank.csv')

@price_ratio()
class SCGZyclonicTreatedWaterTank(SanUnit):
    '''
    Treated water tank (3000 L volume) for the SCG Zyclonic system is used
    to hold water to be reused as flush water.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    The following impact items should be pre-constructed for life cycle assessment:
    Polyethylene, InjectionMolding.

    References
    ----------
    Communication with the design team, communication files include:

        [1] Illinois_Requested_Info_Pack.pdf

        [2] Manual of AQ Tank Installation (16-Apr-2021).pdf

        [3] SCG Follow-up Questions 11-15-2021.docx

        [4] UIUC_info_from_AIT.docx
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        data = load_data(path=treated_water_tank_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)
            
    def _init_lca(self):
        self.construction = [
            Construction(item='Polyethylene', linked_unit=self, quantity_unit='kg'),
            Construction(item='InjectionMolding', linked_unit=self, quantity_unit='kg'),
            ]

    def _design(self):  # Material quantities for LCA calculations
        design = self.design_results
        constr = self.construction
        design['Polyethylene'] = constr[0].quantity = self.mass_tank
        design['InjectionMolding'] = constr[1].quantity = self.mass_tank
        self.add_construction(add_cost=False)

    def _cost(self):
        self.baseline_purchase_costs['TreatedWaterTank'] = self.price_tank * self.price_ratio