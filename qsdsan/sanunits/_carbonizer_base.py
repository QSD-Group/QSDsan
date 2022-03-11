#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Lewis Rowles <stetsonsc@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>
    Lane To <lane20@illinois.edu>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# %%

from warnings import warn
from .. import SanUnit, Construction
from ..utils import ospath, load_data, data_path, price_ratio

__all__ = ('CarbonizerBase',)

cb_path = ospath.join(data_path, 'sanunit_data/_carbonizer_base.tsv')

@price_ratio(default_price_ratio=1)
class CarbonizerBase(SanUnit):
    '''
    Carbonizer base is where feedstock is continuously fed into the pyrolysis
    pot where it is exposed to high temperature pyrolysis. This process
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
    ins : stream obj
        Dewatered solids moisture content ≤ 35%.
    outs : Iterable(stream obj)
        Biochar, hot gas, fugitive N2O.

    References
    ----------
    Rowles et al., Financial viability and environmental sustainability of fecal sludge
    management with Omni Processors. In Prep.
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

        data = load_data(path=cb_path)
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
        if hpd > 24:
            warn(f'The calcualted run time for `CarbonizerBase` is {hpd:.2f} per day, '
                 'which is unrealistic.')
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
        C['Misc. parts'] = \
            self.pyrolysis_pot_cost_cb + \
            self.primary_air_blower_cost_cb + \
            self.thermocouple_cost_cb_pcd + \
            self.thermistor_cost_cb_pcd + \
            self.forced_air_fan_cost_cb + \
            self.airlock_motor_cost_cb + \
            self.inducer_fan_cost_cb + \
            self.biochar_collection_box_cost_cb + \
            self.drive_chain_cost_cb + \
            self.chain_guards_cost_cb + \
            self.door_cost_cb + \
            self.agitator_cost_cb + \
            self.combusion_chamber_cost_cb + \
            self.sprayer_cost_cb + \
            self.vent_cost_cb

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        # O&M cost converted to annual basis, labor included,
        # USD/yr only accounts for time running
        tot_hr = self.daily_run_time * 365
        replacement_parts_annual_cost = \
            tot_hr/self.pyrolysis_pot_lifetime_cb * self.pyrolysis_pot_cost_cb + \
            tot_hr/self.char_auger_motor_lifetime_cb * self.char_auger_motor_cost_cb + \
            tot_hr/self.fuel_auger_motor_lifetime_cb * self.fuel_auger_motor_cost_cb + \
            1/self.primary_air_blower_lifetime_cb * self.primary_air_blower_cost_cb + \
            1/self.thermocouple_lifetime_cb_2pcd * self.thermocouple_cost_cb_pcd + \
            1/self.thermistor_lifetime_cb_2pcd * self.thermistor_cost_cb_pcd + \
            1/self.forced_air_fan_lifetime_cb * self.forced_air_fan_cost_cb + \
            1/self.airlock_motor_lifetime_cb * self.airlock_motor_cost_cb + \
            1/self.inducer_fan_lifetime_cb * self.inducer_fan_cost_cb
        replacement_parts_annual_cost *= ratio

        annual_maintenance = \
            (self.service_team_greasebearing_cb+self.service_team_tighten_cb+self.service_team_adjustdoor_cb) * 12 + \
            1/self.frequency_corrective_maintenance * (
                self.service_team_replacegasket_cb +
                self.service_team_replacedoor_cb +
                self.service_team_replacechain_cb +
                self.service_team_changefirepot_cb +
                self.service_team_replacecharaugers_cb) + \
            tot_hr/self.char_auger_motor_lifetime_cb * self.service_team_replacecharmotor_cb + \
            tot_hr*self.fuel_auger_motor_lifetime_cb * self.service_team_replacefuelmotor_cb + \
            1/self.primary_air_blower_lifetime_cb * self.service_team_replaceblower_cb + \
            1/self.inducer_fan_lifetime_cb * self.service_team_replacefan_cb + \
            1/self.frequency_corrective_maintenance * self.service_team_replacepaddleswitch_cb + \
            1/self.airlock_motor_lifetime_cb * self.service_team_replaceairlock_cb
        annual_maintenance *= self.service_team_wages / 60

        self.add_OPEX = ( # USD/hr (all items are per hour)
            replacement_parts_annual_cost+annual_maintenance) / (365 * 24)

        power_demand = \
            self.carbonizer_biochar_auger_power + \
            self.carbonizer_fuel_auger_power + \
            self.carbonizer_primary_air_blower_power
        self.power_utility(power_demand)  # kWh