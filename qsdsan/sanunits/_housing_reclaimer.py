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


__all__ = ('HousingReclaimer',)

housing_data_path = ospath.join(data_path,'sanunit_data/_housing_reclaimer.csv')


@price_ratio(default_price_ratio=1)
class HousingReclaimer(SanUnit):
    '''
    Structural housing for the Reclaimer 2.0.

    The following impact items should be pre-constructed for life cycle assessment:
    Steel.

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
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', ppl=1, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo=thermo, init_with=init_with, F_BM_default=1)

        self.ppl = ppl

        self.qty_reclaimers = np.ceil(self.ppl / self.baseline_ppl)  # number of reclaimer units required
        self.qty_toilet = np.ceil(self.ppl / 25)  # assume 25 people per MURT toilet
  
        data = load_data(path=housing_data_path)
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
        design['Steel'] = steel_quant = (self.steel_weight + (self.framework_weight/4) + self.fittings_weight) * \
                                        (self.ppl / self.baseline_ppl)  # linear scale
        self.construction = Construction(item='Steel', quantity=steel_quant, quantity_unit='kg')
        self.add_construction(add_cost=False)

    def _cost(self):
        C = self.baseline_purchase_costs
        C['Housing'] = (self.frame + self.extrusion + self.angle_frame + self.angle + self.door_sheet +
                        self.plate_valve + self.powder) * (1 + 0.1 * (self.qty_reclaimers - 1)) + \
                        self.portable_toilet * self.qty_toilet
        
        ratio = self.price_ratio
        for equipment, cost in C.items():
            C[equipment] = cost * ratio
