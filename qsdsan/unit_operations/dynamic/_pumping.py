#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

    Joy Zhang <joycheung1994@gmail.com>

    Jianan Feng <jiananf2@illinois.edu>

Part of this module is based on the biosteam package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from ... import SanUnit
from ..bst import Pump

__all__ = ('HydraulicDelay',)


class HydraulicDelay(Pump):
    '''
    A fake unit for implementing hydraulic delay by a first-order reaction
    (i.e., a low-pass filter) with a specified time constant [d].

    See Also
    --------
    `Benchmark Simulation Model No.1 implemented in MATLAB & Simulink <https://www.cs.mcgill.ca/~hv/articles/WWTP/sim_manual.pdf>`
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, t_delay=1e-4, *,
                 init_with='WasteStream', F_BM_default=None, isdynamic=True):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default,
                         isdynamic=isdynamic)
        self.t_delay = t_delay
        self._concs = None

    def set_init_conc(self, **kwargs):
        '''set the initial concentrations [mg/L].'''
        Cs = np.zeros(len(self.components))
        cmpx = self.components.index
        for k, v in kwargs.items(): Cs[cmpx(k)] = v
        self._concs = Cs

    def _init_state(self):
        '''initialize state by specifying or calculating component concentrations
        based on influents. Total flow rate is always initialized as the sum of
        influent wastestream flows.'''
        self._state = self._ins_QC[0]
        self._dstate = self._state * 0
        if self._concs is not None:
            self._state[:-1] = self._concs

    def _run(self):
        s_in, = self.ins
        s_out, = self.outs
        s_out.copy_like(s_in)

    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE

    def _compile_ODE(self):
        T = self.t_delay
        _dstate = self._dstate
        _update_dstate = self._update_dstate
        def dy_dt(t, QC_ins, QC, dQC_ins):
            Q_in = QC_ins[0,-1]
            Q = QC[-1]
            C_in = QC_ins[0,:-1]
            C = QC[:-1]
            if dQC_ins[0,-1] == 0:
                _dstate[-1] = 0
                _dstate[:-1] = (Q_in*C_in - Q*C)/(Q*T)
            else:
                _dstate[-1] = (Q_in - Q)/T
                _dstate[:-1] = Q_in/Q*(C_in - C)/T
            _update_dstate()
        self._ODE = dy_dt

    def _design(self):
        pass

    def _cost(self):
        pass
