#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

Part of this module is based on the biosteam package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from collections.abc import Iterable
from biosteam.units import Mixer, Splitter, FakeSplitter, ReversedSplitter
from .. import SanUnit

__all__ = (
    'Mixer',
    'Splitter', 'FakeSplitter', 'ReversedSplitter',
    'ComponentSplitter',
    'Sampler',
    )


# %%

class Mixer(SanUnit, Mixer):
    '''
    Similar to :class:`biosteam.units.Mixer`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`,
    and allows dynamic simulation.

    See Also
    --------
    `biosteam.units.Mixer <https://biosteam.readthedocs.io/en/latest/units/mixing.html>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', F_BM_default=None, isdynamic=False,
                 rigorous=False):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with,
                         F_BM_default=F_BM_default, isdynamic=isdynamic)
        self.rigorous = rigorous


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

class Splitter(SanUnit, Splitter):
    '''
    Similar to :class:`biosteam.units.Splitter`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`,
    and allows dynamic simulation.

    See Also
    --------
    `biosteam.units.Splitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, split, order=None,
                  init_with='WasteStream', F_BM_default=None, isdynamic=False):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default,
                         isdynamic=isdynamic)
        self._isplit = self.thermo.chemicals.isplit(split, order)


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


class FakeSplitter(SanUnit, FakeSplitter):
    '''
    Similar to :class:`biosteam.units.FakeSplitter`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.FakeSplitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_
    '''


class ReversedSplitter(SanUnit, ReversedSplitter):
    '''
    Similar to :class:`biosteam.units.ReversedSplitter`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.ReversedSplitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_
    '''


class ComponentSplitter(SanUnit):
    '''
    Split the influent into individual components,
    the last effluent contains all remaining components.

    Parameters
    ----------
    split_keys : iterable
        IDs of components to be splitted to different effluents.
        Element of the item in the iterable can be str or another iterable
        containing component IDs.
        If the item is also iterable, all components whose ID are in the iterable
        will be splitted to the same effluent.
        The split is always 1 for a certain component to an effluent (i.e., complete split).

        .. note::

            Length of the `split_keys()` (which determines size of the outs) \
            cannot be changed after initiation.

    Examples
    --------
    `bwaise systems <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/systems.py>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', split_keys=()):
        if not split_keys:
            raise ValueError('`split_keys` cannot be empty.')

        if isinstance(split_keys, str):
            self._N_outs = 2
        else:
            self._N_outs = len(split_keys) + 1
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

        self._split_keys = split_keys


    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    _graphics = Splitter._graphics


    def _run(self):
        last = self.outs[-1]
        last.mix_from(self.ins)

        splitted = []
        for num, cmps in enumerate(self.split_keys):
            if isinstance(cmps, str):
                cmps = (cmps,)

            elif not isinstance(cmps, Iterable):
                raise ValueError('`split_keys` must be an iterable, '
                                 f'not {type(cmps).__name__}.')

            for cmp in cmps:
                self.outs[num].imass[cmp] = last.imass[cmp]
                last.imass[cmp] = 0
                if cmp in splitted:
                    raise ValueError(f'The component {cmps} appears more than once in `split_keys`.')
                splitted.append(cmp)


    @property
    def split_keys(self):
        '''
        [iterable] IDs of components to be splitted to different effluents.
        Element of the item in the iterable can be str or another iterable
        containing component IDs.
        If the item is also iterable, all components whose ID are in the iterable
        will be splitted to the same effluent.
        The split is always 1 for a certain component to an effluent (i.e., complete split).

        .. note::

            Length of the `split_keys()` (which determines size of the outs) \
                cannot be changed after initiation.
        '''
        return self._split_keys
    @split_keys.setter
    def split_keys(self, i):
        if isinstance(i, str):
            i = (i,)

        if len(i) != len(self.outs):
            raise ValueError('Size of `split_keys` cannot be changed after initiation.')

        self._split_keys = i

class Sampler(SanUnit):
    
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