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

from biosteam.units import Pump as BSTPump
from .. import SanUnit

__all__ = ('Pump',)


class Pump(SanUnit, BSTPump):
    '''
    Similar to the :class:`biosteam.units.Pump`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`,
    and allows dynamic simulation.

    See Also
    --------
    `biosteam.units.Pump <https://biosteam.readthedocs.io/en/latest/units/Pump.html>`_
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                  P=None, pump_type='Default', material='Cast iron',
                  dP_design=101325, ignore_NPSH=True,
                  init_with='Stream', F_BM_default=None, isdynamic=False):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default,
                         isdynamic=isdynamic)
        self.P = P
        self.pump_type = pump_type
        self.material = material
        self.dP_design = dP_design
        self.ignore_NPSH = ignore_NPSH

    @property
    def state(self):
        '''The state of the Pump, including component concentrations [mg/L] and flow rate [m^3/d].'''
        if self._state is None: return None
        else:
            return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    def _init_state(self):
        self._state = self._ins_QC[0]
        self._dstate = self._state * 0.

    def _update_state(self):
        '''updates conditions of output stream based on conditions of the Pump'''
        self._outs[0].state = self._state

    def _update_dstate(self):
        '''updates rates of change of output stream from rates of change of the Pump'''
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
