#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


# %%

from biosteam._graphics import UnitGraphics
from .. import SanUnit

__all__ = ('LumpedCost',)

class LumpedCost(SanUnit):
    '''
    A unit that does not affect stream flow (e.g., outs will be copied from ins),
    but only for cost calculation purpose.
    
    Parameters
    ----------
    cost_item_name : str
        The name of this unit that will be shown in design, cost, and result dicts.
    CAPEX : float
        Total installed capital cost.
    power : float
        Total electricity usage.
    add_OPEX : float
        Additional operating cost per hour.
    
    '''
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 cost_item_name='Lumped cost',
                 CAPEX=0., power=0., add_OPEX=0., **kwargs):
        try: iter(ins)
        except: ins = (ins,)
        self._N_outs = self._N_ins = len(ins)
        self._graphics = UnitGraphics.box(self._N_ins, self._N_outs)
        SanUnit.__init__(self, ID, ins, outs, thermo)
        self.F_BM = {cost_item_name: 1}
        self.baseline_purchase_costs = {cost_item_name: CAPEX}
        self.power_utility(power)
        self._add_OPEX = add_OPEX
        for attr, val in kwargs.items():
            setattr(self, attr, val)
    
    def _run(self):
        for num, stream in enumerate(self.ins):
            self.outs[num].copy_like(stream)

    def _setup(self):
        pass










