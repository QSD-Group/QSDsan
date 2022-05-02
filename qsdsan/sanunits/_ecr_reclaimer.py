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
from ..utils import ospath, load_data, data_path


__all__ = ('ECR_Reclaimer',)

ECR_data_path = ospath.join(data_path, 'sanunit_data/_ECR_Reclaimer.csv')


class ECR_Reclaimer(SanUnit):
    '''
    Electrochemical treatment with chlorine dosing.

    The following impact items should be pre-constructed for life cycle assessment:
    Titanium.

    Parameters
    ----------
    ppl: int
        Total number of users for scaling of costs.
    if_gridtied: bool
        If using grid electricity instead of photovoltaic electricity.

    ins:
        waste: liquid waste stream to be treated by electrochemical unit

    outs:
        treated: treated liquid leaving electrochemical unit

    References
    ----------
    [1] 2019.06 Technical report for BMGF V3 _ CC 2019.06.13.pdf
    
    '''

    # Constants
    baseline_ppl = 30  # baseline population served by Reclaimer
    exponent_scale = 0.6  # exponential scaling constant

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_gridtied=True, ppl=1, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)
        
        self.ppl = ppl
        self.if_gridtied = if_gridtied

        self.qty_reclaimers = np.ceil(self.ppl / self.baseline_ppl)  # number of reclaimer units required
        self.qty_reclaimers = self.qty_reclaimers.astype(int)  # convert from float to integer

        data = load_data(path=ECR_data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)

    _N_ins = 1
    _N_outs = 1

    def _run(self):
        waste = self.ins[0]
        treated = self.outs[0]
        treated.copy_like(self.ins[0])

    def _design(self):
        design = self.design_results
        design['Titanium'] = electrode_quant = self.Titanium_weight * (self.ppl / self.baseline_ppl)  # linear scale
        self.construction = Construction(item='Titanium', quantity=electrode_quant, quantity_unit='kg')
        self.add_construction(add_cost=False)
 
    def _cost(self):

        C = self.baseline_purchase_costs

        C['EC_brush'] = self.EC_brush
        C['EC_cell'] = self.EC_cell
        
        # Exponentially scale capital cost with number of users
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        for equipment, cost in C.items():
            C[equipment] = cost * scale

        self.add_OPEX = self._calc_replacement_cost()

        if self.if_gridtied:  # scale linearly with number of units
            self.power_demand = (self.power_demand_ecr / 1000) * self.qty_reclaimers  # [W/day][1 kW/1000 W] = [kWh/d]
        else:
            self.power_demand = 0

        self.power_utility(self.power_demand)  # kWh/day
    
    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        ecr_brush_replacement_cost = scale * (self.EC_brush / self.EC_brush_lifetime)
        ecr_cell_replacement_cost = scale * (self.EC_cell / self.EC_cell_lifetime)
        ecr_replacement_cost = ecr_brush_replacement_cost + ecr_cell_replacement_cost
        ecr_replacement_cost = ecr_replacement_cost / (365 * 24)  # convert from USD/year to USD/hour
        return ecr_replacement_cost
