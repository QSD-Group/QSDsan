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

import numpy as np
from qsdsan import SanUnit, Construction
from ..utils import ospath, load_data, data_path, price_ratio


__all__ = ('SystemReclaimer',)

system_data_path = ospath.join(data_path, 'sanunit_data/_system_reclaimer.csv')


@price_ratio(default_price_ratio=1)
class SystemReclaimer(SanUnit):
    '''
    System connection components for the Reclaimer 2.0.

    The following impact items should be pre-constructed for life cycle assessment:
    Steel.

    Parameters
    ----------
    ppl: int
        Total number of users for scaling of costs.
    if_gridtied: bool
        If using grid electricity instead of photovoltaic electricity.

    ins: none

    outs: none

    References
    ----------
    [1] Duke Reclaimer team data
    
    '''

    # Constants
    baseline_ppl = 30  # baseline population served by Reclaimer

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_gridtied=True, ppl=1, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        self.ppl = ppl
        self.if_gridtied = if_gridtied

        self.qty_reclaimers = np.ceil(self.ppl / self.baseline_ppl)  # number of reclaimer units required

        data = load_data(path=system_data_path)
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
        design['Steel'] = steel_quant = self.steel_weight * (self.ppl / self.baseline_ppl)  # linear scale
        self.construction = Construction(item='Steel', quantity=steel_quant, quantity_unit='kg')
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs

        C['System'] = (self.T_nut + self.die_cast_hinge + self.SLS_locks + self.DC_round_key + self.handle_rod +
                       self.eight_mm_bolt + self.button_headed_nut + self.twelve_mm_bolt + self.ten_mm_CSK +
                       self.sixteen_mm_bolt + self.coupling_brass + self.socket + self.onehalf_tank_nipple +
                       self.onehalf_in_coupling_brass + self.onehalf_in_fitting + self.plate + self.pump +
                       self.three_way_valve + self.lofted_tank) * (1 + 0.05 * (self.qty_reclaimers - 1))

        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio

        self.add_OPEX = self._calc_replacement_cost()

        if self.if_gridtied:
            self.power_demand = (self.power_demand_system * 24 / 1000) * self.qty_reclaimers  # [W/day][24 hr operation][1 kW/1000 W] = [kWh/d]
        else:
            self.power_demand = 0

        self.power_utility(self.power_demand)  # kWh/day
        
    def _calc_replacement_cost(self):
        system_capital_cost = (self.T_nut + self.die_cast_hinge + self.SLS_locks + self.DC_round_key + self.handle_rod +
                               self.eight_mm_bolt + self.button_headed_nut + self.twelve_mm_bolt +
                               self.ten_mm_CSK + self.sixteen_mm_bolt + self.coupling_brass + self.socket +
                               self.onehalf_tank_nipple + self.onehalf_in_coupling_brass + self.onehalf_in_fitting +
                               self.plate + self.pump + self.three_way_valve + self.lofted_tank) * (1 + 0.05 * (self.qty_reclaimers - 1)) * self.price_ratio
        system_replacement_cost = system_capital_cost * self.om_capital_ratio
        system_replacement_cost = system_replacement_cost / (365 * 24)  # convert from USD/year to USD/hour
        return system_replacement_cost
