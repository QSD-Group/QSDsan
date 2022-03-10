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

from .. import SanUnit, Construction
from ..utils import ospath, load_data, data_path, price_ratio

__all__ = ('PollutionControlDevice',)

device_path = ospath.join(data_path, 'sanunit_data/_pollution_control_device.tsv')

@price_ratio(default_price_ratio=1)
class PollutionControlDevice(SanUnit):
    '''
    Pollution Control Device’s primary responsibilities include pollution
    control and pre heating of the feedstock.

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
    ins : Iterable(stream obj)
        hot gas, fugitive N2O (from Carbonizer Base)
    outs : Iterable(stream obj)
        hot gas, fugitive N2O.

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
            Construction('catalytic_converter', linked_unit=self, item='CatalyticConverter', quantity_unit='ea'),
            )

        data = load_data(path=device_path)
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
        C['Misc. parts'] = \
            self.o2_sensor_cost_pcd + \
            self.thermocouple_cost_cb_pcd + \
            self.thermistor_cost_cb_pcd + \
            self.input_auger_pcd + \
            self.catylist_pcd + \
            self.drive_chain_pcd + \
            self.catalyst_access_door_pcd + \
            self.feedstock_staging_bin_pcd + \
            self.bindicator_pcd + \
            self.feedstock_staging_assembly_pcd + \
            self.temperature_limit_switch_pcd + \
            self.airlock_pcd

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        # O&M cost converted to annual basis, labor included,
        # USD/yr only accounts for time running
        hpd = self.daily_run_time # hr per day
        replacement_parts_annual_cost = \
            hpd*365/self.o2_sensor_lifetime_pcd * self.o2_sensor_cost_pcd + \
            1/self.thermocouple_lifetime_cb_2pcd * 2 * self.thermocouple_cost_cb_pcd + \
            1/self.thermistor_lifetime_cb_2pcd * 2 * self.thermistor_cost_cb_pcd
        replacement_parts_annual_cost *= ratio

        annual_maintenance = \
            self.service_team_adjustdoor_pcd * 12 + \
            1/self.frequency_corrective_maintenance * self.service_team_replacecatalyst_pcd + \
            1/self.frequency_corrective_maintenance * self.service_team_replacebrick_pcd + \
            hpd*365/self.o2_sensor_lifetime_pcd * self.service_team_replaceo2sensor_pcd
        annual_maintenance *= self.service_team_wages / 60

        self.add_OPEX =  (replacement_parts_annual_cost+annual_maintenance) / (365 * 24) # USD/hr

        power_demand = (self.pcd_auger_power+self.pcd_airlock_power) * hpd / 24
        self.power_utility(power_demand) # kWh/hr