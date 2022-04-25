# !/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Tori Morgan <tvlmorgan@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

from qsdsan import SanUnit, Construction
from ..utils import ospath, load_data, data_path, price_ratio

__all__ = ('SolarReclaimer',)

solar_data_path = ospath.join(data_path, 'sanunit_data/_solar_reclaimer.csv')


@price_ratio(default_price_ratio=1)
class SolarReclaimer(SanUnit):
    '''
    Photovoltaic system for solar power generation.

    The following impact items should be pre-constructed for life cycle assessment:
    Battery, Solar.

    Parameters
    ----------
    ppl: int
        Total number of users for scaling of costs.

    ins: none

    outs: none

    References
    ----------
    [1] Duke Reclaimer team data

    '''

    # Constants
    baseline_ppl = 30  # baseline population served by Reclaimer
    exponent_scale = 0.6  # exponential scaling constant

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', ppl=1, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        self.ppl = ppl

        data = load_data(path=solar_data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def _run(self):
        waste = self.ins[0]
        treated = self.outs[0]
        treated.copy_like(self.ins[0])

    def _design(self):
        design = self.design_results
        design['Solar'] = solar_quant = self.solar_capacity * (self.ppl / self.baseline_ppl)  # linear scale
        design['Battery'] = battery_quant = self.battery_kg * (self.ppl / self.baseline_ppl)  # linear scale

        self.construction = Construction(item='Solar', quantity=solar_quant, quantity_unit='m2')
        self.construction = Construction(item='Battery', quantity=battery_quant, quantity_unit='kg')

        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['Battery System'] = self.battery_storage_cost + self.battery_holder_cost
        C['Solar Cost'] = ((self.solar_cost * self.power_demand_30users) + self.solar_module_system
                           + self.inverter_cost)

        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio * scale

        self.add_OPEX = self._calc_replacement_cost() + self._calc_maintenance_labor_cost()

    def _calc_maintenance_labor_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        labor_cost = (self.wages * self.pannel_cleaning) * scale  # USD/year
        return labor_cost / (365 * 24)  # USD/hr

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        solar_panel_replacement_cost = (self.solar_cost * self.power_demand_30users) / self.solar_lifetime  # USD/year
        battery_replacement_cost = self.battery_storage_cost / self.battery_lifetime  # USD/year
        replacement_cost = (solar_panel_replacement_cost + battery_replacement_cost) / (365 * 24) * self.price_ratio * scale  # USD/hr
        return replacement_cost
