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

from ... import SanUnit
from ...utils import auom

__all__ = ('Screening',)


class Screening(SanUnit):
    '''
    A non-reactive unit used to estimate the operating cost of screening.

    Note that only costs from electricity and screened out solids disposal are considered
    (i.e., no equipment cost).

    Parameters
    ----------
    solids_yield : float
        Amount of solids that is screened out, [ft3/hr/MGD].
    compaction : float
        Fraction of the solids that can be compacted
        (i.e., volume after compaction = original volume * (1-compaction)).
    disposal_cost : float
        Cost of compacted solids disposal, [$/ft3].
    power_demand : float
        Power usage for screening, [kW/MGD].

    Examples
    --------
    >>> from qsdsan.utils import create_example_wwt_components
    >>> cmps = create_example_wwt_components()
    >>> from qsdsan import System, WasteStream
    >>> from qsdsan.unit_operations import Screening
    >>> inf = WasteStream('inf')
    >>> inf.set_flow_by_concentration(flow_tot=20,
    ...     concentrations={'Substrate': 300, 'X_inert': 100}, units=('mgd', 'mg/L'))
    >>> SC = Screening('SC', ins=inf, outs='screened')
    >>> sys = System('sys', path=(SC,))
    >>> sys.simulate()
    >>> round(SC.power_utility.rate, 2)  # kW
    14.91
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 solids_yield=2, compaction=0.75,
                 disposal_cost=225/20/27, # converted from $225/20 yd3 container
                 power_demand=1*0.7457): # converted from 1 hp
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self.solids_yield = solids_yield
        self.compaction = compaction
        self.disposal_cost = disposal_cost
        self.power_demand = power_demand


    def _cost(self):
        Q_mgd = self.Q_mgd
        solids = self.solids_yield*Q_mgd*(1-self.compaction)
        self.add_OPEX = {'Solids disposal cost': solids*self.disposal_cost}
        self.power_utility.consumption = Q_mgd*self.power_demand


    @property
    def Q_mgd(self):
        '''
        [float] Influent volumetric flow rate in million gallon per day, [mgd].
        '''
        return auom('m3').convert(self.ins[0].F_vol, 'gallon')*24/1e6