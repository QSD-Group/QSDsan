# !/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Tori Morgan <tvlmorgan@gmail.com>
    Hannah Lohman <hlohman94@gmail.com>
    Lewis Rowles <stetsonsc@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


import numpy as np
from qsdsan import SanUnit, Construction
from ..utils import ospath, load_data, data_path

__all__ = ('UltrafiltrationReclaimer',)

ultrafiltration_data_path = ospath.join(data_path, 'sanunit_data/_ultrafiltration_reclaimer.csv')


class UltrafiltrationReclaimer(SanUnit):
    '''
    Ultrafiltration for removing suspended solids with automated backwash
    to prolong filter

    The following impact items should be pre-constructed for life cycle assessment:
    Plastic, Steel.

    Parameters
    ----------
    ppl: int
        Total number of users for scaling of costs.
    if_gridtied: bool
        If using grid electricity instead of photovoltaic electricity.

    ins:
        waste: liquid waste stream to be treated by ultrafiltration unit

    outs:
        treated: treated liquid leaving ultrafiltration unit
        retentate: concentrated retentate leaving ultrafiltration unit

    References
    ----------
    [1] Duke Reclaimer team data

    '''

    # Constants
    baseline_ppl = 30  # baseline population served by Reclaimer
    exponent_scale = 0.6  # exponential scaling constant

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 if_gridtied=True, ppl=1, **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        self.if_gridtied = if_gridtied
        self.ppl = ppl

        self.qty_reclaimers = np.ceil(self.ppl / self.baseline_ppl)  # number of reclaimer units required
        self.qty_reclaimers = self.qty_reclaimers.astype(int)  # convert from float to integer

        data = load_data(path=ultrafiltration_data_path)
        for para in data.index:
            value = float(data.loc[para]['expected'])
            setattr(self, para, value)
        del data
        
        for attr, value in kwargs.items():
            setattr(self, attr, value)
        
    _N_ins = 1  # waste
    _N_outs = 2  # treated, retentate

    def _run(self):
        waste = self.ins[0]
        treated, retentate = self.outs
        treated.copy_like(self.ins[0])   
        
        self.retentate_prcd = (waste.imass['OtherSS'] * self.TSS_removal)  # mg/L
        retentate.imass['OtherSS'] = self.retentate_prcd
    
    def _design(self):
        design = self.design_results
        design['Plastic'] = plastic_quant = self.Plastic_weight * (self.ppl / self.baseline_ppl)  # linear scale
        design['Steel'] = steel_quant = self.Steel_weight * (self.ppl / self.baseline_ppl)  # linear scale

        self.construction = (Construction(item='Plastic', quantity=plastic_quant, quantity_unit='kg'),
                             Construction(item='Steel', quantity=steel_quant, quantity_unit='kg'))

        self.add_construction(add_cost=False)

    def _cost(self):

        C = self.baseline_purchase_costs
        C['Pipes'] = self.one_in_pipe_SCH40 + self.onehalf_in_pipe_SCH40 + self.three_in_pipe_SCH80
        C['Fittings'] = self.one_in_elbow_SCH80 + self.one_in_tee_SCH80 + self.one_in_SCH80 \
                        + self.one_onehalf_in_SCH80 + self.onehalf_in_SCH80 + self.three_in_SCH80_endcap \
                        + self.one_one_NB_MTA + self.one_onehalf_NB_MTA + self.foot_valve \
                        + self.one_onehalf_in_SCH80_threadedtee + self.three_in_pipe_clamp + self.one_in_pipe_clamp \
                        + self.onehalf_in_pipe_clamp + self.two_way_valve + self.UF_brush
        C['UF_unit'] = self.UF_unit

        # Exponentially scale capital cost with number of users
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale
        for equipment, cost in C.items():
            C[equipment] = cost * scale

        self.add_OPEX = self._calc_replacement_cost()

        if self.if_gridtied:
            if self.qty_reclaimers == 1:
                self.power_demand = self.power_demand_1 / 1000  # [W/day][1 kW/1000 W] = [kW/d]
            elif self.qty_reclaimers == 2:
                self.power_demand = self.power_demand_2 / 1000  # [W/day][1 kW/1000 W] = [kW/d]
            elif self.qty_reclaimers == 3:
                self.power_demand = self.power_demand_3 / 1000  # [W/day][1 kW/1000 W] = [kW/d]
            else:
                self.power_demand = self.power_demand_4 / 1000  # [W/day][1 kW/1000 W] = [kW/d]
        else:
            self.power_demand = 0
            
        self.power_utility(self.power_demand)

    def _calc_replacement_cost(self):
        scale = (self.ppl / self.baseline_ppl) ** self.exponent_scale

        pipe_replacement_cost = (self.one_in_pipe_SCH40 / self.one_in_pipe_SCH40_lifetime +
                                 self.onehalf_in_pipe_SCH40 / self.onehalf_in_pipe_SCH40_lifetime +
                                 self.three_in_pipe_SCH80 / self.three_in_pipe_SCH80_lifetime)

        fittings_replacement_cost = (self.one_in_elbow_SCH80 / self.one_in_elbow_SCH80_lifetime +
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
                                     self.UF_brush / self.UF_brush_lifetime)

        uf_replacement_cost = self.UF_unit / self.UF_unit_lifetime

        total_replacement_cost = scale * (pipe_replacement_cost + fittings_replacement_cost + uf_replacement_cost)  # USD/year
        return total_replacement_cost / (365 * 24)  # USD/hr
