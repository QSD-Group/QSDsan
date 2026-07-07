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

Note: `unit_operations/bst/_abstract.py` and `unit_operations/static/_abstract.py`
are unrelated modules that happen to share this filename; each holds a different
set of unit classes. `Mixer`/`Splitter` live here rather than in `bst/_abstract.py`
because they support dynamic simulation.
'''

import numpy as np
from biosteam.units import (
    Mixer as BSTMixer,
    Splitter as BSTSplitter,
    )
from ... import SanUnit

__all__ = (
    'Sampler',
    'HydraulicDelay',
    'Mixer',
    'Splitter',
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


# %%

class HydraulicDelay(SanUnit):
    '''
    A fake unit for implementing hydraulic delay by a first-order reaction
    (i.e., a low-pass filter) with a specified time constant [d].

    See Also
    --------
    `Benchmark Simulation Model No.1 implemented in MATLAB & Simulink <https://www.cs.mcgill.ca/~hv/articles/WWTP/sim_manual.pdf>`
    '''

    _N_ins = 1
    _N_outs = 1

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

    def _update_state(self):
        self._outs[0].state = self._state

    def _update_dstate(self):
        self._outs[0].dstate = self._dstate

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


# %%

class Mixer(SanUnit, BSTMixer):
    '''
    Similar to :class:`biosteam.units.Mixer`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`,
    and allows dynamic simulation.

    See Also
    --------
    `biosteam.units.Mixer <https://biosteam.readthedocs.io/en/latest/units/mixing.html>`_
    '''
    _graphics = BSTMixer._graphics
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 rigorous=False, conserve_phases=False):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic)
        self.rigorous = rigorous
        self.conserve_phases = conserve_phases


    @property
    def state(self):
        '''The state of the Mixer, including component concentrations [mg/L] and flow rate [m^3/d].'''
        if self._state is None: return None
        else:
            return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    def _init_state(self):
        '''initialize state by specifying or calculating component concentrations
        based on influents. Total flow rate is always initialized as the sum of
        influent wastestream flows.'''
        QCs = self._ins_QC
        if QCs.shape[0] <= 1: self._state = QCs[0]
        else:
            Qs = QCs[:,-1]
            Cs = QCs[:,:-1]
            self._state = np.append(Qs @ Cs / Qs.sum(), Qs.sum())
        self._dstate = self._state * 0.

    def _update_state(self):
        '''updates conditions of output stream based on conditions of the Mixer'''
        self._outs[0].state = self._state

    def _update_dstate(self):
        '''updates rates of change of output stream from rates of change of the Mixer'''
        self._outs[0].dstate = self._dstate

    @property
    def AE(self):
        if self._AE is None:
            self._compile_AE()
        return self._AE

    def _compile_AE(self):
        _n_ins = len(self.ins)
        _state = self._state
        _dstate = self._dstate
        _update_state = self._update_state
        _update_dstate = self._update_dstate
        def yt(t, QC_ins, dQC_ins):
            if _n_ins > 1:
                Q_ins = QC_ins[:, -1]
                C_ins = QC_ins[:, :-1]
                dQ_ins = dQC_ins[:, -1]
                dC_ins = dQC_ins[:, :-1]
                Q = Q_ins.sum()
                C = Q_ins @ C_ins / Q
                _state[-1] = Q
                _state[:-1] = C
                Q_dot = dQ_ins.sum()
                C_dot = (dQ_ins @ C_ins + Q_ins @ dC_ins - Q_dot * C)/Q
                _dstate[-1] = Q_dot
                _dstate[:-1] = C_dot
            else:
                _state[:] = QC_ins[0]
                _dstate[:] = dQC_ins[0]
            _update_state()
            _update_dstate()
        self._AE = yt


# %%

class Splitter(SanUnit, BSTSplitter):
    '''
    Similar to :class:`biosteam.units.Splitter`,
    but can be initialized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`,
    and allows dynamic simulation.

    See Also
    --------
    `biosteam.units.Splitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_
    '''
    _graphics = BSTSplitter._graphics
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, split, order=None,
                  init_with='WasteStream', F_BM_default=None, isdynamic=False):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default,
                         isdynamic=isdynamic)
        self._isplit = self.components.isplit(split, order)


    @property
    def state(self):
        '''Component concentrations and total flow rate.'''
        if self._state is None: return None
        else:
            return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    def _init_state(self):
        self._state = self._ins_QC[0]
        self._dstate = self._state * 0.
        s = self.split
        s_flow = s[self.components.index('H2O')]
        self._split_out0_state = np.append(s/s_flow, s_flow)
        self._split_out1_state = np.append((1-s)/(1-s_flow), 1-s_flow)

    def _update_state(self):
        '''updates conditions of output stream based on conditions of the Splitter'''
        arr = self._state
        self._outs[0].state = self._split_out0_state * arr
        self._outs[1].state = self._split_out1_state * arr

    def _update_dstate(self):
        '''updates rates of change of output stream from rates of change of the Splitter'''
        arr = self._dstate
        self._outs[0].dstate = self._split_out0_state * arr
        self._outs[1].dstate = self._split_out1_state * arr

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
