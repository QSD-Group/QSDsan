#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


# %%

from collections.abc import Iterable
from biosteam._graphics import UnitGraphics
from .. import SanUnit

__all__ = (
    'Copier',
    'LumpedCost',
    )


class Copier(SanUnit):
    '''
    A non-reactive unit that simply copy all the influents to the effluents
    (in the same order), i.e., it does not affect stream flow.

    Note that for a unit without the `_run` function and with only one influent,
    the effluent will be automatically copied from the influent, but this
    `Copier` class will work for any number of influents.

    Parameters
    ----------
    ins : Iterable(obj)
        Any number of influents, all will be copied to the effluents.
    ins : Iterable(obj)
        Any number of effluents, all copied from the corresponding influents.
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=1, **kwargs):
        if isinstance(ins, str) or (not isinstance(ins, Iterable)):
            self._N_outs = self._N_ins = 1
        else:
            self._N_outs = self._N_ins = len(ins)
        self._graphics = UnitGraphics.box(self._N_ins, self._N_outs)

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with,
                         F_BM_default=F_BM_default)

        for attr, val in kwargs.items():
            setattr(self, attr, val)


    def _run(self):
        for num, stream in enumerate(self.ins):
            self.outs[num].copy_like(stream)


class LumpedCost(Copier):
    '''
    A non-reactive unit only for cost calculation purpose.

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

    See Also
    --------
    :class:`qsdsan.sanunits.Copier`_

    Examples
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 cost_item_name='Lumped cost',
                 CAPEX=0., power=0., add_OPEX=0., **kwargs):
        Copier.__init__(self, ID, ins, outs, thermo, init_with)
        self.F_BM = {cost_item_name: 1}
        self.CAPEX_dct = {cost_item_name: CAPEX}
        self.power_utility(power)
        self._add_OPEX = add_OPEX
        for attr, val in kwargs.items():
            setattr(self, attr, val)


    def _cost(self):
        self.baseline_purchase_costs.update(self.CAPEX_dct)