# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Cheung <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('FlatBottomCircularClarifier', )

from .. import SanUnit
from math import exp
# from sympy import symbols, lambdify, Matrix
# from scipy.integrate import solve_ivp
# from warnings import warn
import numpy as np
import pandas as pd


def _settling_flux(X, v_max, v_max_practical, X_min, rh, rp):
    X_star = max(X-X_min, 0)
    v = min(v_max_practical, v_max*(exp(-rh*X_star) - exp(-rp*X_star)))
    return X*max(v, 0)

class FlatBottomCircularClarifier(SanUnit):
    """
    A flat-bottom circular clarifier with a simple 1-dimensional
    N-layer settling model.

    Parameters
    ----------
    ID : str
        ID for the clarifier. The default is ''.
    ins : :class:`WasteStream`
        Influent to the clarifier. Expected number of influent is 1.
    outs : :class:`WasteStream`
        Treated effluent and sludge.
    sludge_flow_rate : float, optional
        Designed sludge flowrate (WAS + RAS), in [m^3/d]. The default is 2000.
    surface_area : float, optional
        Surface area of the clarifier, in [m^2]. The default is 1500.
    height : float, optional
        Height of the clarifier, in [m]. The default is 4.
    N_layer : int, optional
        The number of layers to model settling. The default is 10.
    feed_layer : int, optional
        The feed layer counting from top to bottom. The default is 4.
    X_threshold : float, optional
        Threshold suspended solid cocentration, in [g/m^3]. The default is 3000.
    v_max : float, optional
        Maximum theoretical (i.e. Vesilind) settling velocity, in [m/d]. The
        default is 474.
    v_max_practical : float, optional
        Maximum practical settling velocity, in [m/d]. The default is 250.
    rh : float, optional
        Hindered zone settling parameter in the double-exponential settling velocity
        function, in [m^3/g]. The default is 5.76e-4.
    rp : float, optional
        Flocculant zone settling parameter in the double-exponential settling velocity
        function, in [m^3/g]. The default is 2.86e-3.
    fns : float, optional
        Non-settleable fraction of the suspended solids, dimensionless. Must be within
        [0, 1]. The default is 2.28e-3.

    References
    ----------
    .. [1] Takács, I.; Patry, G. G.; Nolasco, D. A Dynamic Model of the Clarification
        -Thickening Process. Water Res. 1991, 25 (10), 1263–1271.
        https://doi.org/10.1016/0043-1354(91)90066-Y.

    """

    _N_ins = 1
    _N_outs = 2

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', sludge_flow_rate=2000,
                 surface_area=1500, height=4, N_layer=10, feed_layer=4,
                 X_threshold=3000, v_max=474, v_max_practical=250,
                 rh=5.76e-4, rp=2.86e-3, fns=2.28e-3,
                 isdynamic=True, **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, isdynamic=isdynamic)
        self._Qs = sludge_flow_rate
        self._V = surface_area * height
        self._A = surface_area
        self._hj = height/N_layer
        self._N_layer = N_layer
        self._feed_layer = feed_layer
        self._v_max = v_max
        self._v_max_p = v_max_practical
        self._X_t = X_threshold
        self._rh = rh
        self._rp = rp
        self._fns = fns
        # self._cache_state = cache_state
        for attr, value in kwargs.items():
            setattr(self, attr, value)
        # self._state = None

    @property
    def state(self):
        '''Component concentrations in each layer and total flow rate.'''
        if self._state is None: return None
        else:
            Cs = pd.DataFrame(self._state[:-1].reshape((self._N_layer, len(self.components))),
                              columns=self.components.IDs,
                              index={'layer': range(1, self._N_layer + 1)})
            Q = self._state[-1]
            return {'Concentrations [mg/L]': Cs,
                    'Total flow rate [m3/d]': Q}

    @state.setter
    def state(self, QCs):
        QCs = np.asarray(QCs)
        if QCs.shape != (self._N_layer * len(self.components) + 1,):
            raise ValueError(f'state must be a 1d-array of shape {(self._N_layer * len(self.components) + 1,)}')
        self._state = QCs

    def _init_state(self, state=None):
        if state: Cs = state.flatten()
        else: Cs = np.repeat(self.ins[0].Conc, self._N_layer, axis=0)
        Q = self.ins[0].get_total_flow('m3/d')
        self._state = np.append(Cs, Q)

    def _state_locator(self, arr):
        '''derives conditions of output stream from conditions within the clarifier'''
        dct = {}
        Q_s = self._Qs
        Q = arr[-1]
        Q_e = Q - Q_s
        Cs = arr[:-1].reshape((self._N_layer, len(self.components)))
        dct[self.outs[0].ID] = np.append(Cs[0], Q_e)
        dct[self.outs[1].ID] = np.append(Cs[-1], Q_s)
        dct[self.ID] = arr
        return dct

    def _dstate_locator(self, arr):
        '''derives rates of change of output streams from rates of change within the clarifier'''
        dct = {}
        dQ = arr[-1]
        dCs = arr[:-1].reshape((self._N_layer, len(self.components)))
        dct[self.outs[0].ID] = np.append(dCs[0], dQ)
        dct[self.outs[1].ID] = np.append(dCs[-1], 0)
        dct[self.ID] = dCs
        return dct

    def _load_state(self):
        '''returns a dictionary of values of state variables within the clarifer and in the output streams.'''
        if self._state is None: self._init_state()
        return self._state_locator(self._state)

    def _run(self):
        '''only to converge volumetric flows.'''
        inf, = self.ins
        Q_in = inf.get_total_flow('m3/d')
        eff, sludge = self.outs
        Q_s = self._Qs
        Q_e = Q_in - Q_s
        inf.split_to(eff, sludge, split=Q_e/Q_in)


    @property
    def _ODE(self):
        n = self._N_layer
        jf = self._feed_layer - 1
        if jf not in range(self._N_layer):
            raise ValueError(f'feed layer {self._feed_layer} is out of range.'
                              f'must be an integer between 1 and {self._N_layer}.')
        x = self.components.x
        imass = self.components.i_mass
        Q_s = self._Qs
        fns = self._fns
        vmax = self._v_max
        vmaxp = self._v_max_p
        rh = self._rh
        rp = self._rp
        X_t = self._X_t
        A = self._A
        hj = self._hj

        def dy_dt(QC_ins, QC, dQC_ins):
            Q_in = QC_ins[-1]
            Q_e = max(Q_in - Q_s, 0)
            C_in = QC_ins[:-1]
            X_in = sum(C_in*imass*x)
            X_composition = C_in*x/X_in
            Z_in = C_in*(1-x)
            X_min = X_in * fns
            C = QC[:-1].reshape((n,len(x)))
            X = np.dot(C, imass*x)      # (n,) array, TSS for each layer
            Z = C*(1-x)                 # (n, m) array, solubles for each layer
            Q_jout = np.array([Q_e if j < jf else Q_in if j == jf else Q_s for j in range(n)])
            #*********particulates***********
            flow_out = X*Q_jout
            flow_in = np.array([Q_e*X[j+1] if j < jf else Q_in*X_in if j == jf else Q_s*X[j-1] for j in range(n)])
            VX = [_settling_flux(xj, vmax, vmaxp, X_min, rh, rp) for xj in X]
            J = [VX[j] if X[j+1] <= X_t and j < jf else min(VX[j], VX[j+1]) for j in range(n-1)]
            settle_out = np.array(J + [0])
            settle_in = np.array([0] + J)
            X_dot = ((flow_in - flow_out)/A + settle_in - settle_out)/hj        # (n,)
            C_x_dot = np.array([X_composition*x_dot for x_dot in X_dot])        # (n, m), 0 wherever not particulate
            #********non-particulates***********
            flow_out = np.array([zj*qjout for zj, qjout in zip(Z, Q_jout)])
            flow_in = np.array([Q_e*Z[j+1] if j < jf else Q_in*Z_in if j == jf else Q_s*Z[j-1] for j in range(n)])
            C_nx_dot = (flow_in - flow_out)/A/hj                                # (n, m), 0 wherever particulate

            C_dot = (C_x_dot + C_nx_dot).flatten()
            Q_dot = dQC_ins[-1]

            return np.append(C_dot, Q_dot)

        return dy_dt


    def _design(self):
        pass