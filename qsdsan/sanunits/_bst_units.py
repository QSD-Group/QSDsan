#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Joy Cheung <joycheung1994@gmail.com>

Part of this module is based on the biosteam package:
https://github.com/BioSTEAMDevelopmentGroup/biosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
import biosteam as bst
from .. import SanUnit, WasteStream


__all__ = (
    'Mixer',
    'Splitter', 'FakeSplitter', 'ReversedSplitter',
    'Pump',
    'Tank', 'StorageTank', 'MixTank',
    'HXutility', 'HXprocess',
    )

# =============================================================================
# Units subclassed from biosteam
# =============================================================================

class Mixer(SanUnit, bst.units.Mixer):
    '''
    Similar to :class:`biosteam.units.Mixer`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`,
    and allows dynamic simulation.

    See Also
    --------
    `biosteam.units.Mixer <https://biosteam.readthedocs.io/en/latest/units/mixing.html>`_
    '''

    @property
    def state(self):
        '''The state of the Mixer, including component concentrations [mg/L] and flow rate [m^3/d].'''
        if self._state is None: return None
        else:
            return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    @state.setter
    def state(self, QCs):
        QCs = np.asarray(QCs)
        if QCs.shape != (len(self.components)+1, ):
            raise ValueError(f'state must be a 1D array of length {len(self.components) + 1},'
                              'indicating component concentrations [mg/L] and total flow rate [m^3/d]')
        self._state = QCs

    def _init_state(self, state=None):
        '''initialize state by specifiying or calculating component concentrations
        based on influents. Total flow rate is always initialized as the sum of
        influent wastestream flows.'''
        mixed = WasteStream()
        mixed.mix_from(self.ins)
        Q = mixed.get_total_flow('m3/d')
        self._state = np.append(state or mixed.Conc, Q)

    def _state_locator(self, arr):
        '''derives conditions of output stream from conditions of the Mixer'''
        dct = {}
        dct[self.outs[0].ID] = arr
        dct[self.ID] = arr
        return dct

    def _dstate_locator(self, arr):
        '''derives rates of change of output stream from rates of change of the Mixer'''
        return self._state_locator(arr)

    def _load_state(self):
        '''returns a dictionary of values of state variables within the CSTR and in the output stream.'''
        if self._state is None: self._init_state()
        return self._state_locator(self._state)

    @property
    def _ODE(self):
        _n_ins = len(self.ins)
        _n_state = len(self.components)+1
        def dy_dt(QC_ins, QC, dQC_ins):
            if _n_ins > 1:
                QC_ins = QC_ins.reshape((_n_ins, _n_state))
                Q_ins = QC_ins[:, -1]
                C_ins = QC_ins[:, :-1]
                dQC_ins = dQC_ins.reshape((_n_ins, _n_state))
                dQ_ins = dQC_ins[:, -1]
                dC_ins = dQC_ins[:, :-1]
                Q = Q_ins.sum()
                C = Q_ins @ C_ins / Q
                Q_dot = dQ_ins.sum()
                C_dot = (dQ_ins @ C_ins + Q_ins @ dC_ins - Q_dot * C)/Q
                return np.append(C_dot, Q_dot)
            else:
                return dQC_ins
        return dy_dt


class Splitter(SanUnit, bst.units.Splitter):
    '''
    Similar to :class:`biosteam.units.Splitter`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`,
    and allows dynamic simulation.

    See Also
    --------
    `biosteam.units.Splitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, split, order=None,
                  init_with='Stream', F_BM_default=None):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default)
        self._isplit = self.thermo.chemicals.isplit(split, order)


    @property
    def state(self):
        '''Component concentrations in each layer and total flow rate.'''
        if self._state is None: return None
        else:
            return dict(zip(list(self.components.IDs) + ['Q'], self._state))

    @state.setter
    def state(self, QCs):
        QCs = np.asarray(QCs)
        if QCs.shape != (len(self.components)+1, ):
            raise ValueError(f'state must be a 1D array of length {len(self.components) + 1},'
                              'indicating component concentrations [mg/L] and total flow rate [m^3/d]')
        self._state = QCs

    def _init_state(self):
        Cs = self.ins[0].Conc
        Q = self.ins[0].get_total_flow('m3/d')
        self._state = np.append(Cs, Q)

    def _state_locator(self, arr):
        '''derives conditions of output stream from conditions of the Splitter'''
        dct = {}
        Q1 = arr[-1] * self.split   # assuming split is a single value for all components
        Q2 = arr[-1] - Q1
        Cs = arr[:-1]
        dct[self.outs[0].ID] = np.append(Cs[0], Q1)
        dct[self.outs[1].ID] = np.append(Cs[-1], Q2)
        dct[self.ID] = arr
        return dct

    def _dstate_locator(self, arr):
        '''derives rates of change of output streams from rates of change of the Splitter'''
        return self._state_locator(arr)

    def _load_state(self):
        '''returns a dictionary of values of state variables within the clarifer and in the output streams.'''
        if self._state is None: self._init_state()
        return self._state_locator(self._state)

    @property
    def _ODE(self):
        def dy_dt(QC_ins, QC, dQC_ins):
            return dQC_ins
        return dy_dt


class FakeSplitter(SanUnit, bst.units.FakeSplitter):
    '''
    Similar to :class:`biosteam.units.FakeSplitter`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.FakeSplitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_
    '''


class ReversedSplitter(SanUnit, bst.units.ReversedSplitter):
    '''
    Similar to :class:`biosteam.units.ReversedSplitter`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.ReversedSplitter <https://biosteam.readthedocs.io/en/latest/units/splitting.html>`_
    '''


class Pump(SanUnit, bst.units.Pump):
    '''
    Similar to the :class:`biosteam.units.Pump`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`,
    and allows dynamic simulation.

    See Also
    --------
    `biosteam.units.Pump <https://biosteam.readthedocs.io/en/latest/units/Pump.html>`_
    '''
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                  P=None, pump_type='Default', material='Cast iron',
                  dP_design=405300, ignore_NPSH=True,
                  init_with='Stream', F_BM_default=None):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default)
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

    @state.setter
    def state(self, QCs):
        QCs = np.asarray(QCs)
        if QCs.shape != (len(self.components)+1, ):
            raise ValueError(f'state must be a 1D array of length {len(self.components) + 1},'
                              'indicating component concentrations [mg/L] and total flow rate [m^3/d]')
        self._state = QCs

    def _init_state(self):
        Cs = self.ins[0].Conc
        Q = self.ins[0].get_total_flow('m3/d')
        self._state = np.append(Cs, Q)

    def _state_locator(self, arr):
        '''derives conditions of output stream from conditions of the Mixer'''
        dct = {}
        dct[self.outs[0].ID] = arr
        dct[self.ID] = arr
        return dct

    def _dstate_locator(self, arr):
        '''derives rates of change of output stream from rates of change of the Mixer'''
        return self._state_locator(arr)

    def _load_state(self):
        '''returns a dictionary of values of state variables within the CSTR and in the output stream.'''
        if self._state is None: self._init_state()
        return self._state_locator(self._state)

    @property
    def _ODE(self):
        def dy_dt(QC_ins, QC, dQC_ins):
            return dQC_ins
        return dy_dt


class Tank(SanUnit, bst.units.Tank, isabstract=True):
    '''
    Similar to the :class:`biosteam.units.Tank`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.Tank <https://biosteam.readthedocs.io/en/latest/units/Tank.html>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                  vessel_type=None, tau=None, V_wf=None,
                  vessel_material=None, kW_per_m3=0.,
                  init_with='Stream', F_BM_default=None):

        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default)

        self.vessel_type = vessel_type or self._default_vessel_type
        self.tau = tau or self._default_tau
        self.V_wf = V_wf or self._default_V_wf
        self.vessel_material = vessel_material or self._default_vessel_material
        self.kW_per_m3 = kW_per_m3 or self._default_kW_per_m3


class StorageTank(Tank, bst.units.StorageTank):
    '''
    Similar to the :class:`biosteam.units.StorageTank`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.StorageTank <https://biosteam.readthedocs.io/en/latest/units/Tank.html>`_
    '''


class MixTank(Tank, bst.units.MixTank):
    '''
    Similar to the :class:`biosteam.units.MixTank`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.MixTank <https://biosteam.readthedocs.io/en/latest/units/Tank.html>`_
    '''


class HXutility(SanUnit, bst.units.HXutility):
    '''
    Similar to :class:`biosteam.units.HXutility`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.HXutility <https://biosteam.readthedocs.io/en/latest/units/heat_exchange.html>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', F_BM_default=None,
                 *, T=None, V=None, rigorous=False, U=None, H=None,
                 heat_exchanger_type="Floating head",
                 material="Carbon steel/carbon steel",
                 N_shells=2, ft=None, heat_only=None, cool_only=None):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default)
        self.T = T #: [float] Temperature of outlet stream (K).
        self.V = V #: [float] Vapor fraction of outlet stream.
        self.H = H #: [float] Enthalpy of outlet stream.

        #: [bool] If true, calculate vapor liquid equilibrium
        self.rigorous = rigorous

        #: [float] Enforced overall heat transfer coefficent (kW/m^2/K)
        self.U = U

        #: Number of shells for LMTD correction factor method.
        self.N_shells = N_shells

        #: User imposed correction factor.
        self.ft = ft

        #: [bool] If True, heat exchanger can only heat.
        self.heat_only = heat_only

        #: [bool] If True, heat exchanger can only cool.
        self.cool_only = cool_only

        self.material = material
        self.heat_exchanger_type = heat_exchanger_type


class HXprocess(SanUnit, bst.units.HXprocess):
    '''
    Similar to :class:`biosteam.units.HXprocess`,
    but can be initilized with :class:`qsdsan.SanStream` and :class:`qsdsan.WasteStream`.

    See Also
    --------
    `biosteam.units.HXprocess <https://biosteam.readthedocs.io/en/latest/units/heat_exchange.html>`_
    '''

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='Stream', F_BM_default=None,
                 *, U=None, dT=5., T_lim0=None, T_lim1=None,
                 material="Carbon steel/carbon steel",
                 heat_exchanger_type="Floating head",
                 N_shells=2, ft=None,
                 phase0=None,
                 phase1=None,
                 H_lim0=None,
                 H_lim1=None):
        SanUnit.__init__(self, ID, ins, outs, thermo,
                         init_with=init_with, F_BM_default=F_BM_default)

        #: [float] Enforced overall heat transfer coefficent (kW/m^2/K)
        self.U = U

        #: [float] Total heat transfered.
        self.Q = None

        #: Number of shells for LMTD correction factor method.
        self.N_shells = N_shells

        #: User imposed correction factor.
        self.ft = ft

        #: [float] Pinch temperature difference.
        self.dT = dT

        #: [float] Temperature limit of outlet stream at index 0.
        self.T_lim0 = T_lim0

        #: [float] Temperature limit of outlet stream at index 1.
        self.T_lim1 = T_lim1

        #: [float] Temperature limit of outlet stream at index 0.
        self.H_lim0 = H_lim0

        #: [float] Temperature limit of outlet stream at index 1.
        self.H_lim1 = H_lim1

        #: Enforced phase of outlet at index 0
        self.phase0 = phase0

        #: Enforced phase of outlet at index 1
        self.phase1 = phase1

        self.material = material
        self.heat_exchanger_type = heat_exchanger_type