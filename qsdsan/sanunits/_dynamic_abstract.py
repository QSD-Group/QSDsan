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
from .. import SanUnit

__all__ = (
    'Sampler',
    )


# %%

class Sampler(SanUnit):
    '''
    A non-reactive (i.e., all the outs at the same as the ins) unit that
    is used in dynamic simulation to record the unit/stream states.
    '''

    _N_ins = 1
    _N_outs = 1

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default,
                         isdynamic=isdynamic)

    def _run(self):
        inf, = self.ins
        out, = self.outs
        out.copy_like(inf)

    @property
    def state(self):
        '''The sampled state, including component concentrations [mg/L] and flow rate [m^3/d].'''
        if self._state is None: return None
        else:
            return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    def _init_state(self):
        self._state = self._ins_QC[0]
        self._dstate = self._state * 0.

    def _update_state(self):
        self._outs[0].state = self._state

    def _update_dstate(self):
        self._outs[0].dstate = self._dstate

    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        def yt(t, QC_ins, dQC_ins):
            _state[:] = QC_ins[0]
            _dstate[:] = dQC_ins[0]
            _update_state()
            _update_dstate()
        self._AE = yt
