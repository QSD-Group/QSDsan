#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Lewis Rowles <stetsonsc@gmail.com>

    Yalin Li <mailto.yalin.li@gmail.com>

    Hannah Lohman <hlohman94@gmail.com>

    Lane To <lane20@illinois.edu>

This module contains unit operations used in the Biogenic Refinery system
as described in `Rowles et al. <https://doi.org/10.1021/acsenvironau.2c00022>`_

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from warnings import warn
from math import ceil
from .. import SanUnit, Construction
from ._sludge_thickening import SludgeSeparator
from ..utils import ospath, load_data, data_path, price_ratio

__all__ = (
    'BiogenicRefineryCarbonizerBase',
    'BiogenicRefineryControls',
    'BiogenicRefineryGrinder',
    'BiogenicRefineryHHX',
    'BiogenicRefineryHHXdryer',
    'BiogenicRefineryHousing',
    'BiogenicRefineryIonExchange',
    'BiogenicRefineryOHX',
    'BiogenicRefineryPollutionControl',
    'BiogenicRefineryScrewPress',
    'BiogenicRefineryStruvitePrecipitation',
    )

br_su_data_path = ospath.join(data_path, 'sanunit_data/br')


# %%

br_carbonizer_path = ospath.join(br_su_data_path, '_br_carbonizer_base.tsv')

class BiogenicRefineryCarbonizerBase(SanUnit):
    '''
    Carbonizer base in the biogenic refinery is where feedstock is continuously
    fed into the pyrolysis pot.
    The feedstock is exposed to high temperature pyrolysis. This process
    produces biochar and hot gases.

    The carbonizer base is the central location for the combined pyrolysis
    and combustion process. The feedstock is received into the pyrolysis pot
    where it is flash pyrolyzed releasing volatile gases which are then mixed
    with air which is pumped in through the Primary Blower causing combustion
    of the gases. Due to the combustion of gases, the carbonizer base is the
    hottest location of the refinery ranging between 550-900°C. This
    range is closely monitored as pyrolysis usually starts at 350°C and
    all of the volatile gases are released at 550°C. This temperature is
    also important when confirming treatment of the fecal sludge, as it
    serves as our evidence for inactivation of the microbes in the sludge.
    After the gases are combusted, the exhaust travels below a baffle plate
    to encourage the fallout of particulates before it proceeds to the
    pollution control device.

    The following components should be included in system thermo object for simulation:
    H2O, N, K, P OtherSS, N2O.

    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, Steel, ElectricMotor, Electronics.

    Parameters
    ----------
    ins : Iterable(stream)
        Dewatered solids moisture content ≤ 35%.
    outs : Iterable(stream)
        Biochar, hot gas, fugitive N2O.

    References
    ----------
    [1] Rowles et al., Financial viability and environmental sustainability of
    fecal sludge treatment with Omni Processors, ACS Environ. Au, 2022,
    https://doi.org/10.1021/acsenvironau.2c00022
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)

        self.construction = (
            Construction('stainless_steel', linked_unit=self, item='StainlessSteel', quantity_unit='kg'),
            Construction('steel', linked_unit=self, item='Steel', quantity_unit='kg'),
            Construction('electric_motor', linked_unit=self, item='ElectricMotor', quantity_unit='ea'),
            Construction('electronics', linked_unit=self, item='Electronics', quantity_unit='kg'),
            )

        data = load_data(path=br_carbonizer_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    _N_ins = 1
    _N_outs = 3

    def _run(self):
        waste = self.ins[0]
        biochar, gas, N2O = self.outs
        biochar.copy_like(self.ins[0])
        gas.phase = N2O.phase = 'g'
        gas.T = self.pyrolysis_temp

        # Moisture content
        mc = waste.imass['H2O']/waste.F_mass if waste.F_mass!=0 else 0
        if mc > 0.351: # allow a small error
            warn(f'Moisture content of the influent is {mc:.1%}, '
                'larger than the maximum allowed level of 35%.')

        biochar.empty()
        biochar_prcd = waste.F_mass * (1-mc) * self.biochar_production_rate # kg biochar /hr
        biochar.imass['C'] = waste.COD * self.carbon_COD_ratio * waste.F_vol / 1e3 * (1 - self.pyrolysis_C_loss)
        NPK = ('N', 'P', 'K')
        for element in NPK:
            biochar.imass[element] *= 1 - getattr(self, f'pyrolysis_{element}_loss')
        biochar.imass['OtherSS'] = biochar_prcd - biochar.imass[('C', *NPK)].sum()
        biochar.imass['H2O'] = 0.025 * biochar.F_mass # kg H20 / hr with 2.5% moisture content

        # Daily run time for the influent waste
        # One unit is capable of treating 550 kg/d (35% moisture based on 20 hr of run time) or 18 kg dry/hr
        self.daily_run_time = hpd = waste.F_mass*(1-mc) * 24/self.loading_rate # hr/d
        self.uptime_ratio = hpd / 24 # ratio of uptime (all items are per hour)

        # # Energy balance, not really being used
        # thermal_energy_from_feedstock = waste.F_mass * (1-mc) * self.dry_feces_heat_of_combustion # MJ/hr
        # # Calculate thermal energy entering and thermal energy required to dry off water
        # heat_needed_to_dry_0 = waste.F_mass * mc * self.energy_required_to_dry_sludge # MJ/hr
        # self.net_thermal_energy_in = thermal_energy_from_feedstock  - heat_needed_to_dry_0 # MJ/hr

        # N2O emissions
        N_to_gas = waste.imass['N'] * self.pyrolysis_N_loss # kg N / hr
        N2O_from_HCNO = N_to_gas * self.N_to_HCNO * self.HCNO_to_NH3 * self.NH3_to_N2O
        N20_from_NH3 = N_to_gas * self.N_to_NH3 * self.NH3_to_N2O
        N2O_emissions = N2O_from_HCNO + N20_from_NH3 # kg N2O / hr
        N2O.imass['N2O'] = N2O_emissions


    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = \
            self.carbonizer_base_assembly_stainless + \
            self.carbonizer_base_squarebox_stainless + \
            self.carbonizer_base_charbox_stainless
        design['Steel'] = constr[1].quantity = \
            self.carbonizer_base_assembly_steel + \
            self.carbonizer_base_squarebox_steel + \
            self.carbonizer_base_charbox_steel
        design['ElectricMotor'] =  constr[2].quantity = 2.7/5.8 + 6/5.8
        design['Electronics'] = constr[3].quantity = 1
        self.add_construction(add_cost=False)


    def _cost(self):
        D = self.design_results
        C = self.baseline_purchase_costs
        C['Stainless steel'] = self.stainless_steel_cost * D['StainlessSteel']
        C['Steel'] = self.steel_cost * D['Steel']
        C['Electric motors'] = self.char_auger_motor_cost_cb + self.fuel_auger_motor_cost_cb
        C['Misc. parts'] = (
            self.pyrolysis_pot_cost_cb +
            self.primary_air_blower_cost_cb +
            self.thermocouple_cost_cb_pcd +
            self.thermistor_cost_cb_pcd +
            self.forced_air_fan_cost_cb +
            self.airlock_motor_cost_cb +
            self.inducer_fan_cost_cb +
            self.biochar_collection_box_cost_cb +
            self.drive_chain_cost_cb +
            self.chain_guards_cost_cb +
            self.door_cost_cb +
            self.agitator_cost_cb +
            self.combusion_chamber_cost_cb +
            self.sprayer_cost_cb +
            self.vent_cost_cb
            )

        # O&M cost converted to annual basis, labor included,
        # USD/yr only accounts for time running
        tot_hr = self.daily_run_time * 365
        replacement_parts_annual_cost = (
            tot_hr/self.pyrolysis_pot_lifetime_cb * self.pyrolysis_pot_cost_cb +
            tot_hr/self.char_auger_motor_lifetime_cb * self.char_auger_motor_cost_cb +
            tot_hr/self.fuel_auger_motor_lifetime_cb * self.fuel_auger_motor_cost_cb +
            1/self.primary_air_blower_lifetime_cb * self.primary_air_blower_cost_cb +
            1/self.thermocouple_lifetime_cb_2pcd * self.thermocouple_cost_cb_pcd +
            1/self.thermistor_lifetime_cb_2pcd * self.thermistor_cost_cb_pcd +
            1/self.forced_air_fan_lifetime_cb * self.forced_air_fan_cost_cb +
            1/self.airlock_motor_lifetime_cb * self.airlock_motor_cost_cb +
            1/self.inducer_fan_lifetime_cb * self.inducer_fan_cost_cb
            )

        # In min/yr
        annual_maintenance = (
            (self.service_team_greasebearing_cb+self.service_team_tighten_cb+self.service_team_adjustdoor_cb) * 12 +
            (
                self.service_team_replacegasket_cb +
                self.service_team_replacedoor_cb +
                self.service_team_replacechain_cb +
                self.service_team_changefirepot_cb +
                self.service_team_replacecharaugers_cb
            ) * 1/self.frequency_corrective_maintenance +
            tot_hr/self.char_auger_motor_lifetime_cb * self.service_team_replacecharmotor_cb +
            tot_hr/self.fuel_auger_motor_lifetime_cb * self.service_team_replacefuelmotor_cb +
            1/self.primary_air_blower_lifetime_cb * self.service_team_replaceblower_cb +
            1/self.inducer_fan_lifetime_cb * self.service_team_replacefan_cb +
            1/self.frequency_corrective_maintenance * self.service_team_replacepaddleswitch_cb +
            1/self.airlock_motor_lifetime_cb * self.service_team_replaceairlock_cb
            )
        annual_maintenance *= self.service_team_wages / 60

        self.add_OPEX = (replacement_parts_annual_cost + annual_maintenance) / (365 * 24)  # USD/hour

        power_demand = \
            self.carbonizer_biochar_auger_power + \
            self.carbonizer_fuel_auger_power + \
            self.carbonizer_primary_air_blower_power
        self.power_utility(power_demand)  # kW


# %%

br_control_path = ospath.join(br_su_data_path,'_br_controls.tsv')

class BiogenicRefineryControls(SanUnit):
    '''
    Control box (industrial control panel) for the biogenic refinery.
    No process algorithm is included, only design (including) cost algorithms are included.

    This is a non-reactive unit (i.e., the effluent is copied from the influent).

    The following impact items should be pre-constructed for life cycle assessment:
    Electronics, ElectricConnectors, ElectricCables.

    Parameters
    ----------
    ins : Iterable(stream)
        Influent stream.
    outs : Iterable(stream)
        Effluent stream, is copied from the influent stream.

    References
    ----------
    [1] Rowles et al., Financial viability and environmental sustainability of
    fecal sludge treatment with Omni Processors, ACS Environ. Au, 2022,
    https://doi.org/10.1021/acsenvironau.2c00022
    '''

    _N_ins = _N_outs = 1

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        self.construction = (
            Construction('electronics', linked_unit=self, item='Electronics', quantity_unit='kg'),
            Construction('electric_connectors', linked_unit=self, item='ElectricConnectors', quantity_unit='kg'),
            Construction('electric_cables', linked_unit=self, item='ElectricCables', quantity_unit='m'),
            )

        data = load_data(path=br_control_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Electronics'] = constr[0].quantity = 2
        design['ElectricConnectors'] = constr[1].quantity = 0.5
        design['ElectricCables'] = constr[2].quantity = 3
        self.add_construction(add_cost=False)


    def _cost(self):
        self.baseline_purchase_costs['Misc. parts'] = (
            self.icp_controller_board +
            self.icp_variable_frequence_drives +
            self.icp_power_meter +
            self.icp_line_filter +
            self.icp_transformer +
            self.icp_power_meter_transformer +
            self.icp_AC_to_DC +
            self.icp_DC_to_AC +
            self.icp_touch_screen
            )

        # O&M cost converted to annual basis, labor included,
        # USD/yr only accounts for time running
        annual_maintenance = (
            (self.electrician_replacecables_icp+self.electrician_replacewires_icp)*self.certified_electrician_wages +
            self.service_team_replacetouchscreen_icp*self.service_team_wages +
            self.facility_manager_configurevariable_icp*self.facility_manager_wages +
            (self.biomass_controls_replaceboard_icp+self.biomass_controls_codemalfunctioning_icp)*self.biomass_controls_wages
            )

        self.add_OPEX =  annual_maintenance/60/self.frequency_corrective_maintenance/(365*24) # USD/hr

        # kWh/hr
        self.power_utility(self.icp_controller_board_power+self.icp_variable_frequence_drives_power)


# %%

br_grinder_path = ospath.join(br_su_data_path, '_br_grinder.tsv')

@price_ratio()
class BiogenicRefineryGrinder(SanUnit):
    '''
    Grinder in the biogenic refinery is used to break up solids.

    .. note:

        Moisture content of the effluent is adjusted to be 65%, although the grinder
        itself can't change the moisture. This assumption was made based on pilot experiments.

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
    [1] Rowles et al., Financial viability and environmental sustainability of
    fecal sludge treatment with Omni Processors, ACS Environ. Au, 2022,
    https://doi.org/10.1021/acsenvironau.2c00022
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 moisture_content_out=0.65, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.moisture_content_out = moisture_content_out
        self.construction = (
            Construction('steel', linked_unit=self, item='Steel', quantity_unit='kg'),
            )

        data = load_data(path=br_grinder_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 1
    _N_outs = 1

    def _run(self):
        waste_in = self.ins[0]
        waste_out = self.outs[0]
        waste_out.copy_like(waste_in)

        mc_in = waste_in.imass['H2O'] / waste_in.F_mass # fraction
        mc_out = self.moisture_content_out
        if mc_in < mc_out:
            raise RuntimeError(f'Moisture content of the influent stream ({mc_in:.2f}) '
                               f'is smaller than the desired moisture content ({mc_out:.2f}).')
        TS_in = waste_in.F_mass - waste_in.imass['H2O'] # kg TS dry/hr
        waste_out.imass['H2O'] = TS_in/(1-mc_out)*mc_out
        waste_out._COD = waste_in.COD * waste_in.F_vol / waste_out.F_vol

    def _design(self):
        self.design_results['Steel'] = self.construction[0].quantity = self.grinder_steel
        self.add_construction(add_cost=False)


    def _cost(self):
        self.baseline_purchase_costs['Grinder'] = self.grinder * self.price_ratio
        self.power_utility(self.grinder_electricity * self.ins[0].imass['OtherSS']) # kWh/hr


# %%

br_hhx_path = ospath.join(br_su_data_path, '_br_hhx.tsv')

@price_ratio()
class BiogenicRefineryHHX(SanUnit):
    '''
    Hydronic heat exchanger in the biogenic refinery is used to dry
    the feedstock before the refinery is capable of processing the
    material. The heat is exchanged between the exhaust gas and water,
    which is then pumped into radiators connected to a dryer.
    The refinery monitors the temperature of the water to ensure that
    the feedstock is being sufficiently dried before entering the refinery.

    .. note::

        The number of pumps in the design results are floats as the costs
        are scaled based on a pump of different size.
        The exponential scaling method might be considered for a better estimation.


    This class should be used together with :class:`~.sanunits.HHXdryer`.

    The following components should be included in system thermo object for simulation:
    N2O.

    The following impact items should be pre-constructed for life cycle assessment:
    OilHeatExchanger, Pump.

    Parameters
    ----------
    ins : Iterable(stream)
        Hot gas.
    outs : Iterable(stream)
        Hot gas.

    Warnings
    --------
    Energy balance is not performed for this unit.

    References
    ----------
    [1] Rowles et al., Financial viability and environmental sustainability of
    fecal sludge treatment with Omni Processors, ACS Environ. Au, 2022,
    https://doi.org/10.1021/acsenvironau.2c00022

    See Also
    --------
    :class:`~.sanunits.HHXdryer`
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        self.construction = (
            Construction('stainless_steel', linked_unit=self, item='StainlessSteel', quantity_unit='kg'),
            Construction('steel', linked_unit=self, item='Steel', quantity_unit='kg'),
            Construction('hydronic_heat_exchanger', linked_unit=self, item='HydronicHeatExchanger', quantity_unit='ea'),
            Construction('pump', linked_unit=self, item='Pump', quantity_unit='ea'),
            )

        data = load_data(path=br_hhx_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 1
    _N_outs = 1

    def _run(self):
        hot_gas_in = self.ins[0]
        hot_gas_out = self.outs[0]
        hot_gas_out.copy_like(hot_gas_in)
        hot_gas_out.phase = hot_gas_in.phase = 'g'
        hot_gas_out.T = self.hhx_temp

        # # The following codes can be used to calculate the heat that was delivered to the HHX
        # self.heat_output_water = \
        #     self.water_flowrate * \
        #     (self.water_out_temp-self.water_in_temp) * \
        #     self.water_density_kg_m_3 * \
        #     self.water_heat_capacity_k_j_kg_k * (60/1000)  # MJ/hr

        # # The following codes can be used to calculate losses through water pipe
        # self.heat_loss_water_pipe = \
        #     (self.water_out_temp-self.ambient_temp) / \
        #     (self.water_r_pipe_k_k_w+self.water_r_insulation_k_k_w) * 3.6 # MJ/hr


    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.heat_exchanger_hydronic_stainless
        design['Steel'] = constr[1].quantity  = self.heat_exchanger_hydronic_steel
        design['HydronicHeatExchanger'] = constr[2].quantity = 1
        design['Pump'] = constr[3].quantity = 17.2/2.72
        self.add_construction(add_cost=False)


    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        C['Stainless steel'] = self.stainless_steel_cost * D['StainlessSteel']
        C['Steel'] = self.steel_cost * D['Steel']
        C['Misc. parts'] = (
            self.hhx_stack +
            self.hhx_stack_thermocouple +
            self.hhx_oxygen_sensor +
            self.hhx_inducer_fan +
            self.hhx_flow_meter +
            self.hhx_pump +
            self.hhx_water_in_thermistor +
            self.hhx_water_out_thermistor +
            self.hhx_load_tank +
            self.hhx_expansion_tank +
            self.hhx_heat_exchanger +
            self.hhx_values +
            self.hhx_thermal_well +
            self.hhx_hot_water_tank +
            self.hhx_overflow_tank
            )

        ratio = self.price_ratio
        for equipment, cost in C.items():
           C[equipment] = cost * ratio

        # O&M cost converted to annual basis, labor included,
        # USD/yr only accounts for time running
        num = 1 / self.frequency_corrective_maintenance
        annual_maintenance = (
            self.service_team_adjustdoor_hhx*12 +
            num * (self.service_team_replacewaterpump_hhx+self.service_team_purgewaterloop_hhx)
            )

        self.add_OPEX =  annual_maintenance * self.service_team_wages / 60 / (365 * 24) # USD/hr (all items are per hour)

        self.power_utility(self.water_pump_power+self.hhx_inducer_fan_power) # kWh/hr


# %%

br_hhx_dryer_path = ospath.join(br_su_data_path, '_br_hhx_dryer.tsv')

@price_ratio()
class BiogenicRefineryHHXdryer(SanUnit):
    '''
    This dryer unit is used in combination with :class:`~.sanunits.HydronicHeatExchanger`
    in the biogenic refinery.

    The following components should be included in system thermo object for simulation:
    H2O, N, CH4, N2O.

    Parameters
    ----------
    ins : Iterable(stream)
        Dewatered solids, heat.
    outs : Iterable(stream)
        Dried solids, fugitive N2O, fugitive CH4.
    moisture_content_out : float
        Desired moisture content of the effluent.

    Warnings
    --------
    Energy balance is not performed for this unit.

    References
    ----------
    [1] Rowles et al., Financial viability and environmental sustainability of
    fecal sludge treatment with Omni Processors, ACS Environ. Au, 2022,
    https://doi.org/10.1021/acsenvironau.2c00022

    See Also
    --------
    :class:`~.sanunits.HydronicHeatExchanger`
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 moisture_content_out=0.35, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        self.moisture_content_out = moisture_content_out

        data = load_data(path=br_hhx_dryer_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    _N_ins = 2
    _N_outs = 3

    def _run(self):
        waste_in, heat_in = self.ins
        waste_out, N2O, CH4  = self.outs
        waste_out.copy_like(self.ins[0])
        heat_in.phase = N2O.phase = CH4.phase = 'g'

        # Calculate heat needed to dry to the desired moisture content
        mc_in = waste_in.imass['H2O'] / waste_in.F_mass # fraction
        mc_out = self.moisture_content_out
        if mc_in < mc_out*0.999: # allow a small error
            warn(f'Moisture content of the influent stream ({mc_in:.1%}) '
                f'is smaller than the desired moisture content ({mc_out:.1%}).')
        TS_in = waste_in.F_mass - waste_in.imass['H2O'] # kg TS dry/hr
        waste_out.imass['H2O'] = TS_in/(1-mc_out)*mc_out

        # # The following codes can be used to compare if the supplied heat is
        # # sufficient to provide the needed heat for drying
        # water_to_dry = waste_in.imass['H2O'] - waste_out.imass['H2O'] # kg water/hr
        # heat_needed_to_dry_35 = water_to_dry * self.energy_required_to_dry_sludge # MJ/hr
        # heat_supplied = self.dryer_heat_transfer_coeff * self.area_surface * (self.water_out_temp-self.feedstock_temp)
        # if heat_needed_to_dry_35 > heat_supplied:
        #     warn('Heat required exceeds heat supplied.')

        # Emissions
        drying_CO2_to_air = (self.drying_CO2_emissions * self.carbon_COD_ratio
                             * waste_in.COD * waste_in.F_vol / 1000) # kg CO2 /hr

        CH4.imass['CH4'] = drying_CH4_to_air = \
            self.drying_CH4_emissions * self.carbon_COD_ratio * \
            waste_in.COD * waste_in.F_vol / 1000 # kg CH4 /hr

        drying_NH3_to_air = self.drying_NH3_emissions * waste_in.TN * waste_in.F_vol / 1000 # kg NH3 /hr
        N2O.imass['N2O'] = drying_NH3_to_air * self.NH3_to_N2O # kg N2O /hr

        # 44/12/16 are the molecular weights of CO2, C, and CH4, respectively
        waste_out._COD = (
            waste_in.COD*waste_in.F_vol -
            (drying_CO2_to_air/44*12+drying_CH4_to_air/16*12) / self.carbon_COD_ratio
            ) / waste_out.F_vol

        # 17/14 are the molecular weights of NH3/N, respectively
        waste_out.imass['N'] -= drying_NH3_to_air / 17 * 14

        # # The following codes can be used to calculate losses due to convection and radiation
        # self.jacket_heat_loss_conv = self.heat_transfer_coeff * (self.water_air_hx_temp - self.ambient_temp) * self.water_air_hx_area / 1000 / 0.2778 # MJ/hr
        # self.jacket_heat_loss_radiation = self.radiative_emissivity * 5.67e-8 * ((self.water_air_hx_temp1 + 273)**4 - (self.ambient_temp + 273)**4) / 1000 / 0.2778 # MJ/hr
        # self.jacket_heat_loss_sum = self.jacket_heat_loss_conv + self.jacket_heat_loss_radiation # MJ/hr


# %%

br_housing_path = ospath.join(br_su_data_path, '_br_housing.tsv')

@price_ratio()
class BiogenicRefineryHousing(SanUnit):
    '''
    Housing for the biogenic refinery which is composed of the casing around the system,
    containers, and the concrete slab.
    No process algorithm is included, only design (including) cost algorithms are included.

    The following impact items should be pre-constructed for life cycle assessment:
    Steel, StainlessSteelSheet, Concrete.

    Parameters
    ----------
    ins : Iterable(stream)
        Influent stream.
    outs : Iterable(stream)
        Effluent stream, is copied from the influent stream.

    References
    ----------
    [1] Rowles et al., Financial viability and environmental sustainability of
    fecal sludge treatment with Omni Processors, ACS Environ. Au, 2022,
    https://doi.org/10.1021/acsenvironau.2c00022
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 const_wage=15, const_person_days=100, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, F_BM_default=1)
        self.const_wage = const_wage
        self.const_person_days = const_person_days

        self.construction = (
            Construction('steel', linked_unit=self, item='Steel', quantity_unit='kg'),
            Construction('stainless_steel_sheet', item='StainlessSteelSheet', linked_unit=self, quantity_unit='kg'),
            Construction('concrete', item='Concrete', linked_unit=self, quantity_unit='m3'),
            )

        data = load_data(path=br_housing_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    def _design(self):
        design = self.design_results
        constr = self.construction
        design['Steel'] = constr[0].quantity = 2000 + 4000
        design['StainlessSteelSheet'] = constr[1].quantity = 4.88 * 2 * 3 * 4.5
        design['Concrete'] = constr[2].quantity = self.concrete_thickness * self.footprint
        self.add_construction(add_cost=False)


    def _cost(self):
        D = self.design_results
        C = self.baseline_purchase_costs
        C['Containers'] = self.container20ft_cost + self.container40ft_cost
        C['Equip Housing'] = D['StainlessSteelSheet'] / 4.88 * self.stainless_steel_housing
        C['Concrete'] = D['Concrete'] * self.concrete_cost
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        # Labor during initial construction, should not be multiplied by `price_ratio`
        C['Labor'] = self.const_wage * self.const_person_days


# %%

br_ix_path = ospath.join(br_su_data_path, '_br_ion_exchange.tsv')

@price_ratio()
class BiogenicRefineryIonExchange(SanUnit):
    '''
    Ion exchange in the biogenic refinery is used for the recovery of N from the liquid stream.
    Concentrated NH3 is recovered.

    The following components should be included in system thermo object for simulation:
    NH3, Polystyrene, H2SO4.

    The following impact items should be pre-constructed for life cycle assessment:
    PVC, PE.

    Parameters
    ----------
    ins : Iterable (stream)
        Liquid waste, fresh resin, H2SO4.
    outs : Iterable (stream)
        Treated waste, spent resin, concentrated NH3.

    References
    ----------
    [1] Lohman et al., Advancing Sustainable Sanitation and Agriculture
    through Investments in Human-Derived Nutrient Systems.
    Environ. Sci. Technol. 2020, 54, (15), 9217-9227.
    https://dx.doi.org/10.1021/acs.est.0c03764

    [2] Tarpeh et al., Evaluating ion exchange for nitrogen recovery from
    source-separated urine in Nairobi, Kenya. Development Engineering. 2018,
    3, 188–195.
    https://doi.org/10.1016/j.deveng.2018.07.002

    [3] Rowles et al., Financial viability and environmental sustainability of
    fecal sludge treatment with Omni Processors, ACS Environ. Au, 2022,
    https://doi.org/10.1021/acsenvironau.2c00022
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)

        self.construction = (
            Construction('pvc', linked_unit=self, item='PVC', quantity_unit='kg'),
            Construction('pe', linked_unit=self, item='PE', quantity_unit='kg'),
            )

        data = load_data(path=br_ix_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    _N_ins = 3
    _N_outs = 3

    def _run(self):
        waste, resin_in, H2SO4 = self.ins
        treated, resin_out, conc_NH3 = self.outs
        treated.copy_like(waste)

        N_recovered = waste.imass['NH3'] * self.N_rec # kg N / hr
        treated.imass['NH3'] =  waste.imass['NH3'] - N_recovered # kg N / hr

        #!!! Technically, here should be using (NH4)2SO4
        conc_NH3.imass['NH3'] = N_recovered # kg N / hr

        resin_demand_influent = waste.TN / self.resin_lifetime / self.ad_density / 14 # kg resin / m3 treated
        resin_demand_time = resin_demand_influent * waste.F_vol # kg resin / hr
        resin_in.imass['Polystyrene'] = resin_out.imass['Polystyrene'] = resin_demand_time

        acid_demand_influent = waste.TN * self.vol_H2SO4 / self.ad_density / 14 # L acid / L treated
        acid_demand_time = acid_demand_influent * waste.F_vol * 1000 * 1.83 # kg acid / hr
        H2SO4.imass['H2SO4'] = acid_demand_time


    def _design(self):
        design = self.design_results
        constr = self.construction
        N_column = self.N_column

        design['PVC'] = constr[0].quantity = \
            N_column * self.column_length * self.pvc_mass # kg PVC
        tubing_quant = N_column * self.tubing_length * self.tubing_mass # kg PE
        design['PE'] = constr[1].quantity = tubing_quant + self.N_tank*self.tank_mass

        self.add_construction(add_cost=False)


    def _cost(self):
        C = self.baseline_purchase_costs
        N_column = self.N_column
        C['PVC'] = self.cost_PVC_column * self.column_length * N_column
        C['Tubing'] = self.cost_tubing * self.tubing_length*N_column
        C['Tank'] = N_column/3 * self.tank_cost # one tank has three columns

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio


    @property
    def N_column(self):
        '''[int] Number of resin columns.'''
        return ceil(self.ins[0].F_vol*1000*24/self.column_daily_loading_rate)

    @property
    def N_tank(self):
        '''[int] Number of tanks for cost estimation, might be float (instead of int).'''
        return self.N_column/3


# %%

br_ohx_path = ospath.join(br_su_data_path,'_br_ohx.tsv')

@price_ratio()
class BiogenicRefineryOHX(SanUnit):
    '''
    Oil heat exchanger in the biogenic refinery utilizes an organic Rankin cycle.
    This is a combined heat and power system and is
    used to generate additional electricity that the refinery and/or
    facility can use to decrease the units electrical demand on the
    electrical grid. This type of system is required for ISO 31800
    certification as the treatment unity needs to be energy independent
    when processing fecal sludge.

    .. note::

        This unit should be used in conjunction of :class:`~.sanunits.CarbonizerBase`
        as it uses the the heat from that unit for heat-exchanging
        (i.e., itself doesn't generate heat).

        The number of oil heat exchanger and pumps in the design results are floats as the costs
        are scaled based on equipment of different sizes.
        The exponential scaling method might be considered for a better estimation.

    The following components should be included in system thermo object for simulation:
    N2O.

    The following impact items should be pre-constructed for life cycle assessment:
    OilHeatExchanger, Pump.

    Parameters
    ----------
    ins : Iterable(stream)
        Hot gas.
    outs : Iterable(stream)
        Hot gas.

    Warnings
    --------
    Energy balance is not performed for this unit.

    References
    ----------
    [1] Rowles et al., Financial viability and environmental sustainability of
    fecal sludge treatment with Omni Processors, ACS Environ. Au, 2022,
    https://doi.org/10.1021/acsenvironau.2c00022
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        self.construction = (
            Construction('oil_heat_exchanger', linked_unit=self, item='OilHeatExchanger', quantity_unit='ea'),
            Construction('pump', linked_unit=self, item='Pump', quantity_unit='ea'),
            )

        data = load_data(path=br_ohx_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 1
    _N_outs = 1

    def _run(self):
        hot_gas_in = self.ins[0]
        hot_gas_out = self.outs[0]
        hot_gas_out.copy_like(hot_gas_in)
        hot_gas_out.phase = hot_gas_in.phase = 'g'
        hot_gas_out.T = self.ohx_temp

        # # The following codes can be used to calculate the power that was delivered to the ORC
        # self.power_delivery_orc = \
        #     self.oil_flowrate * \
        #     (self.oil_temp_out-self.oil_temp_in) * \
        #     self.oil_density * self.oil_specific_heat * (60/1000)  # MJ/hr

        # # The following codes can be used to calculate losses through pipe
        # self.pipe_heat_loss = \
        #     (self.oil_temp_out-self.amb_temp) / \
        #     (self.oil_r_pipe_k_k_w+self.oil_r_insulation_k_k_w) * 3.6  # MJ/hr


    def _design(self):
        design = self.design_results
        constr = self.construction
        design['OilHeatExchanger'] = constr[0].quantity = 4/200
        design['Pump'] = constr[1].quantity = 2.834/2.27
        self.add_construction(add_cost=False)


    def _cost(self):
        self.baseline_purchase_costs['Oil Heat Exchanger'] = self.orc_cost * self.price_ratio
        self.power_utility(self.oil_pump_power-self.oil_electrical_energy_generated) # kWh/hr


# %%

br_pollution_control_path = ospath.join(br_su_data_path, '_br_pollution_control.tsv')

class BiogenicRefineryPollutionControl(SanUnit):
    '''
    Pollution control device in the biogenic refinery is used for
    the pollution control and pre-heating of the feedstock.

    Due to the inefficiencies of the pyrolysis process, there are typically
    pollutants in the exhaust. In order to treat these pollutants, the
    Biogenic Refinery has a catalyst, similar to a catalytic converter in a
    car, to ensure destruction of the pollutants before they can be released
    into the surrounding environment. The process of destroying the pollutants
    requires the catalyst to maintain temperatures above 315°C, and
    additional energy is released during this process. The temperature of
    the catalysis is closely monitored because the catalyst material
    will start to degrade above 615°C and could cause the feedstock to
    prematurely pyrolyze in the fuel auger.

    The following components should be included in system thermo object for simulation:
    N2O.

    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, Steel, ElectricMotor, CatalyticConverter.

    Parameters
    ----------
    ins : Iterable(stream)
        Hot gas, fugitive N2O (from the carbonizer base).
    outs : Iterable(stream)
        Hot gas, fugitive N2O.

    References
    ----------
    [1] Rowles et al., Financial viability and environmental sustainability of
    fecal sludge treatment with Omni Processors, ACS Environ. Au, 2022,
    https://doi.org/10.1021/acsenvironau.2c00022
    '''


    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)

        self.construction = (
            Construction('stainless_steel', linked_unit=self, item='StainlessSteel', quantity_unit='kg'),
            Construction('steel', linked_unit=self, item='Steel', quantity_unit='kg'),
            Construction('electric_motor', linked_unit=self, item='ElectricMotor', quantity_unit='ea'),
            Construction('catalytic_converter', linked_unit=self, item='CatalyticConverter', quantity_unit='ea'),
            )

        data = load_data(path=br_pollution_control_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 2
    _N_outs = 2

    def _run(self):
        gas_in, N2O_in = self.ins
        gas_out, N2O_out = self.outs
        gas_in.phase = N2O_in.phase = gas_out.phase = N2O_out.phase = 'g'

        self.daily_run_time = self.uptime_ratio * 24 # hr/d

        # Set temperature
        gas_out.T = self.catalyst_temp

        # N2O emissions
        N2O_out.imass['N2O'] = N2O_in.imass['N2O'] * 0.9 # kg N2O / hr

    def _design(self):
        design = self.design_results
        constr = self.construction
        design['StainlessSteel'] = constr[0].quantity = self.pcd_cat_sandwich_stainless
        design['Steel'] = constr[1].quantity = self.pcd_cat_sandwich_steel
        design['ElectricMotor'] =  constr[2].quantity = 5/5.8
        design['CatalyticConverter'] = constr[3].quantity = 1
        self.add_construction(add_cost=False)


    def _cost(self):
        D = self.design_results
        C = self.baseline_purchase_costs
        C['Stainless steel'] = self.stainless_steel_cost * D['StainlessSteel']
        C['Steel'] = self.steel_cost * D['Steel']
        C['Electric motor'] = self.input_auger_motor_pcd
        C['Misc. parts'] = (
            self.o2_sensor_cost_pcd +
            self.thermocouple_cost_cb_pcd +
            self.thermistor_cost_cb_pcd +
            self.input_auger_pcd +
            self.catylist_pcd +
            self.drive_chain_pcd +
            self.catalyst_access_door_pcd +
            self.feedstock_staging_bin_pcd +
            self.bindicator_pcd +
            self.feedstock_staging_assembly_pcd +
            self.temperature_limit_switch_pcd +
            self.airlock_pcd
            )

        # O&M cost converted to annual basis, labor included,
        # USD/yr only accounts for time running
        hpd = self.daily_run_time # hr per day
        replacement_parts_annual_cost = (
            hpd*365/self.o2_sensor_lifetime_pcd * self.o2_sensor_cost_pcd +
            1/self.thermocouple_lifetime_cb_2pcd * 2 * self.thermocouple_cost_cb_pcd +
            1/self.thermistor_lifetime_cb_2pcd * 2 * self.thermistor_cost_cb_pcd
            )

        annual_maintenance = (
            self.service_team_adjustdoor_pcd * 12 +
            1/self.frequency_corrective_maintenance * self.service_team_replacecatalyst_pcd +
            1/self.frequency_corrective_maintenance * self.service_team_replacebrick_pcd +
            hpd*365/self.o2_sensor_lifetime_pcd * self.service_team_replaceo2sensor_pcd
            )
        annual_maintenance *= self.service_team_wages / 60

        self.add_OPEX =  (replacement_parts_annual_cost+annual_maintenance) / (365 * 24) # USD/hr

        power_demand = (self.pcd_auger_power+self.pcd_airlock_power) * hpd / 24
        self.power_utility(power_demand) # kWh/hr, or kW


# %%

br_screw_path = ospath.join(br_su_data_path,'_br_screw_press.tsv')

@price_ratio()
class BiogenicRefineryScrewPress(SludgeSeparator):
    '''
    Screw Press is used for dewatering where sludge, conditioned with cationic
    polymer, is fed into the unit. Sludge is continuously dewatered as it travels
    along the screw.

    The following components should be included in system thermo object for simulation:
    Water, Polyacrylamide.

    The following impact items should be pre-constructed for life cycle assessment:
    Steel.

    Parameters
    ----------
    ins : Iterable(stream)
        Waste for treatment (e.g., wastewater or latrine sludge) and polymer added for treatment.
    outs : Iterable(stream)
        Treated liquids and solids.

    References
    ----------
    [1] Tchobanoglous, G.; Stensel, H. D.; Tsuchihashi, R.; Burton, F.; Abu-Orf,
    M.; Bowden, G.; Pfrang, W. Wastewater Engineering: Treatment and Resource
    Recovery, 5th ed.; Metcalf & Eddy, Inc., AECOM, McGraw-Hill: New York, 2014.
    
    [2] Rowles et al., Financial viability and environmental sustainability of
    fecal sludge treatment with Omni Processors, ACS Environ. Au, 2022,
    https://doi.org/10.1021/acsenvironau.2c00022

    See Also
    --------
    :class:`qsdsan.sanunits.SludgeThickening`
    '''

    def __init__(self, ID='', ins=None, outs=(),thermo=None, init_with='WasteStream',
                 split=None, settled_frac=None, if_N2O_emission=False, **kwargs):
        SludgeSeparator.__init__(self, ID, ins, outs, thermo, init_with,
                                 split, settled_frac, F_BM_default=1)

        self.construction = (
            Construction('steel', linked_unit=self, item='Steel', quantity_unit='kg'))

        data = load_data(path=br_screw_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    _N_ins = 2
    _N_outs = 2

    def _run(self):
        waste, polymer = self.ins
        liq, cake_sol = self.outs
        SludgeSeparator._run(self)

        sol_COD = cake_sol.COD * cake_sol.F_vol if cake_sol.COD else 0
        cake_sol.imass['H2O'] = 0
        cake_sol_mass = cake_sol.F_mass / self.cake_solids_TS
        cake_sol.imass['H2O'] = cake_water = cake_sol_mass - cake_sol.F_mass
        if sol_COD: cake_sol._COD = sol_COD / cake_sol.F_vol

        liq_COD = liq.COD * liq.F_vol if liq.COD else 0
        liq.imass['H2O'] = waste.imass['H2O'] - cake_water
        if liq_COD: liq._COD = liq_COD / liq.F_vol

        waste_TS = waste.F_mass - waste.imass['H2O']
        polymer.imass['Polyacrylamide'] = self.dewatering_polymer_dose * waste_TS / 1000 # kg polymer / hr


    def _design(self):
        self.construction[0].quantity = self.dewatering_screw_press_steel
        self.add_construction(add_cost=False)

    def _cost(self):
        self.baseline_purchase_costs['Screw Press'] = self.dewatering_screw_press_cost*self.price_ratio
        self.power_utility(self.dewatering_energy_demand*self.outs[1].F_mass) # kWh/hr


# %%

br_struvite_path = ospath.join(br_su_data_path, '_br_struvite_precipitation.tsv')

@price_ratio()
class BiogenicRefineryStruvitePrecipitation(SanUnit):
    '''
    Struvite precipitation in the biogenic refinery is used for the recovery of P
    from the liquid stream as solid struvite.

    The following components should be included in system thermo object for simulation:
    P, NH3, K, MagnesiumHydroxide, MagnesiumCarbonate, Struvite, FilterBag.

    The following impact items should be pre-constructed for life cycle assessment:
    StainlessSteel, PVC.


    Parameters
    ----------
    ins : Iterable (stream)
        Liquid waste, Mg(OH)2, MgCO3, filter bag.
    outs : Iterable (stream)
        Treated waste, struvite.
    Mg_molar_split : Iterable(float)
        The molar split between Mg(OH)2 and MgCO3.
        (1, 0) means all Mg is added as Mg(OH)2 and (0,1) means all MgCO3.

    References
    ----------
    [1] Lohman et al., Advancing Sustainable Sanitation and Agriculture
    through Investments in Human-Derived Nutrient Systems.
    Environ. Sci. Technol. 2020, 54, (15), 9217-9227.
    https://dx.doi.org/10.1021/acs.est.0c03764

    [2] Tarpeh et al., Evaluating ion exchange for nitrogen recovery from
    source-separated urine in Nairobi, Kenya. Development Engineering. 2018,
    3, 188–195.
    https://doi.org/10.1016/j.deveng.2018.07.002
    
    [3] Rowles et al., Financial viability and environmental sustainability of
    fecal sludge treatment with Omni Processors, ACS Environ. Au, 2022,
    https://doi.org/10.1021/acsenvironau.2c00022
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 Mg_molar_split=(1,0), **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with,
                         F_BM_default=1)
        self.Mg_molar_split = Mg_molar_split

        self.construction = (
            Construction('stainless_steel', linked_unit=self, item='StainlessSteel', quantity_unit='kg'),
            Construction('pvc', linked_unit=self, item='PVC', quantity_unit='kg'),
            )

        data = load_data(path=br_struvite_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data

        for attr, value in kwargs.items():
            setattr(self, attr, value)


    _N_ins = 4
    _N_outs = 2

    def _run(self):
        waste, magnesium_hydroxide, magnesium_carbonate, bag_filter = self.ins
        treated, struvite = self.outs
        treated.copy_like(waste)
        for ws in self.ins[1:]+self.outs[1:]:
            ws.phase = 's'

        Mg_molar_split = self.Mg_molar_split
        sum_splits = sum(Mg_molar_split)
        Mg_dose, P_rec = self.Mg_dose, self.P_rec
        # Total P recovered [kmol P/hr]
        P_recovered = waste.imol['P'] * P_rec
        treated.imol['P'] =  waste.imol['P'] - P_recovered

        # Total N recovered, 1:1 mol/mol N:P in struvite ((NH4)MgPO4·6(H2O))
        treated.imol['NH3'] =  waste.imol['NH3'] - P_recovered

        # Total K recovered [kg K/hr]
        treated.imol['K'] =  waste.imol['K'] * (1-self.K_rec)

        magnesium_hydroxide.empty()
        magnesium_carbonate.empty()
        magnesium_hydroxide.imol['MagnesiumHydroxide'] = waste.imol['P']*Mg_dose*(Mg_molar_split[0]/sum_splits)
        magnesium_carbonate.imol['MagnesiumCarbonate'] = waste.imol['P']*Mg_dose*(Mg_molar_split[1]/sum_splits)
        struvite.imol['Struvite'] = P_recovered
        bag_filter.imass['FilterBag'] = self.N_tank * self.cycles_per_day / self.filter_reuse / 24 # bags/hr


    def _design(self):
        design = self.design_results
        constr = self.construction
        N_tank = self.N_tank
        design['StainlessSteel'] = constr[0].quantity = N_tank * self.reactor_weight # kg SS
        design['PVC'] = constr[1].quantity = N_tank * self.material_P_pipe * self.pvc_mass # kg PVC
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        N_tank = self.N_tank
        C['Reactor'] = N_tank * self.cost_P_reactor
        C['Stirrer'] = N_tank * self.cost_P_stirrer
        C['PVC'] = N_tank * self.material_P_pipe * self.cost_P_pipe

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = 0.35/4 * self.baseline_purchase_cost / (365 * 24)  # USD/hr (all items are per hour)


    @property
    def N_tank(self):
        '''[int] Number of reactor tanks.'''
        return ceil(self.ins[0].F_vol*1000*24/self.cycles_per_day/self.reactor_volume)