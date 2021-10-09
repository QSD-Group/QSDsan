# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from .. import SanUnit
import numpy as np
# import pandas as pd

__all__ = ('FlatBottomCircularClarifier', 
           'IdealClarifier',)

def _settling_flux(X, v_max, v_max_practical, X_min, rh, rp):
    X_star = max(X-X_min, 0)
    v = min(v_max_practical, v_max*(np.exp(-rh*X_star) - np.exp(-rp*X_star)))
    return X*max(v, 0)

class FlatBottomCircularClarifier(SanUnit):
    """
    A flat-bottom circular clarifier with a simple 1-dimensional
    N-layer settling model. [1]_

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
        self.feed_layer = feed_layer
        self._v_max = v_max
        self._v_max_p = v_max_practical
        self._X_t = X_threshold
        self._rh = rh
        self._rp = rp
        self._fns = fns
        self._ODE = None
        self._solids = None
        self._solubles = None
        self._X_comp = np.zeros(len(self.components))
        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def reset_cache(self):
        '''Reset cached states.'''
        self._state = None
        for s in self.outs:
            s.empty()

    @property
    def sludge_flow_rate(self):
        '''[float] The designed sludge flow rate (wasted + recycled) in m3/d.'''
        return self._Qs

    @sludge_flow_rate.setter
    def sludge_flow_rate(self, Qs):
        self._Qs = Qs

    @property
    def V_settle(self):
        '''[float] Total volume modeled for settling in m^3, calculated based on surface area and height.'''
        return self._V

    @property
    def A_settle(self):
        '''[float] The surface area for settling in m^2, i.e., the area of the clarifer's flat bottom.'''
        return self._A

    @A_settle.setter
    def A_settle(self, A):
        self._A = A

    @property
    def h_layer(self):
        '''[float] The height of each layer in the settling model, in m.'''
        return self._hj

    @h_layer.setter
    def h_layer(self, h):
        self._hj = h

    @property
    def N_layer(self):
        '''[int] The number of layers into which the clarifier is divided in the settling model.'''
        return self._N_layer

    @N_layer.setter
    def N_layer(self, N):
        self._N_layer = int(N)

    @property
    def feed_layer(self):
        '''[int] The feed layer counting from top to bottom.'''
        return self._feed_layer

    @feed_layer.setter
    def feed_layer(self, jf):
        jf = int(jf)
        if jf > self._N_layer or jf < 1:
            raise ValueError(f'feed layer {self._feed_layer} is out of range.'
                             f'must be an integer between 1 and {self._N_layer}.')
        self._feed_layer = jf

    @property
    def v_max(self):
        '''[float] Maximum theoretical (i.e. Vesilind) settling velocity, in m/d'''
        return self._v_max

    @v_max.setter
    def v_max(self, vm):
        if vm < self._v_max_p:
            raise ValueError('v_max must be greater than or equal to v_max_p.')
        self._v_max = vm

    @property
    def v_max_p(self):
        '''[float] Maximum practical settling velocity, in m/d'''
        return self._v_max_p

    @v_max_p.setter
    def v_max_p(self, vmp):
        if vmp > self._v_max or vmp <= 0:
            raise ValueError('v_max_p must be within (0, v_max].')
        self._v_max_p = vmp

    @property
    def X_t(self):
        '''[float] Threshold suspended solid cocentration, in g/m^3.'''
        return self._X_t

    @X_t.setter
    def X_t(self, xt):
        if xt < 0: raise ValueError('X_t must be positive.')
        self._X_t = xt

    @property
    def rh(self):
        '''[float] Hindered zone settling parameter in the double-exponential settling velocity function, in m^3/g.'''
        return self._rh

    @rh.setter
    def rh(self, rh):
        if rh > self._rp: raise ValueError('rh must be less than or equal to rp.')
        self._rh = rh

    @property
    def rp(self):
        '''[float] Flocculant zone settling parameter in the double-exponential settling velocity function, in m^3/g.'''
        return self._rp

    @rp.setter
    def rp(self, rp):
        if rp < self._rh: raise ValueError('rp must be greater than or equal to rp.')
        self._rp = rp

    @property
    def fns(self):
        '''[float] Non-settleable fraction of the suspended solids'''
        return self._fns

    @fns.setter
    def fns(self, fns):
        if fns < 0 or fns > 1: raise ValueError('fns must be within [0,1].')
        self._fns = fns

    def set_init_solubles(self, **kwargs):
        '''set the initial concentrations [mg/L] of solubles in the clarifier.'''
        Cs = np.zeros(len(self.components))
        cmpx = self.components.index
        x = self.components.x
        for k, v in kwargs.items(): Cs[cmpx(k)] = v
        # self._solubles = np.tile(Cs*(1-x), self._N_layer)
        self._solubles = Cs*(1-x)

    def set_init_TSS(self, arr):
        '''set the initial TSS [mg/L] in each layer of the clarifier.'''
        if len(arr) != self._N_layer:
            raise ValueError(f'expects an iterable of length {self._N_layer}, not {len(arr)}')
        self._solids = np.asarray(arr, dtype=float)

    def _init_state(self):
        n = self._N_layer
        x = self.components.x
        imass = self.components.i_mass
        QCs = self._state_tracer()[0]
        Q = QCs[-1]
        Z = self._solubles if self._solubles is not None \
            else QCs[:-1]*(1-x)
            # else np.tile(QCs[:-1]*(1-x), n)
        TSS_in = sum(QCs[:-1] * x * imass)
        TSS = self._solids if self._solids is not None \
            else np.array([TSS_in*f for f in 20**np.linspace(-1,1,n)])
        ZQs = np.append(Z, Q)
        self._state = np.append(ZQs, TSS)
        if TSS_in != 0: self._X_comp = QCs[:-1] * x / TSS_in

    def _state_locator(self, arr):
        x = self.components.x
        n = self._N_layer
        dct = {}
        Q = arr[-(1+n)]
        Q_e = Q - self._Qs
        Z = arr[:len(x)]
        X_composition = self._X_comp # (m, ), mg COD/ mg TSS
        X_e = arr[-n] * X_composition
        X_s = arr[-1] * X_composition
        dct[self.ID] = arr
        dct[self.outs[0].ID] = np.append(Z+X_e, Q_e)
        dct[self.outs[1].ID] = np.append(Z+X_s, self._Qs)        
        return dct
    
    def _dstate_locator(self, arr):
        x = self.components.x
        n = self._N_layer
        dct = {}
        dQ = arr[-(1+n)]
        dZ = arr[:len(x)]
        X_composition = self._X_comp # (m, ), mg COD/ mg TSS
        dX_e = arr[-n] * X_composition
        dX_s = arr[-1] * X_composition
        dct[self.ID] = arr
        dct[self.outs[0].ID] = np.append(dZ+dX_e, dQ)
        dct[self.outs[1].ID] = np.append(dZ+dX_s, 0)
        return dct

    def _load_state(self):
        '''returns a dictionary of values of state variables within the clarifer and in the output streams.'''
        if self._state is None: self._init_state()
        return {self.ID: self._state}

    def _run(self):
        '''only to converge volumetric flows.'''
        inf, = self.ins
        Q_in = inf.get_total_flow('m3/d')
        eff, sludge = self.outs
        Q_s = self._Qs
        Q_e = max(Q_in - Q_s, 0)
        inf.split_to(eff, sludge, split=Q_e/Q_in)


    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE
    
    def _compile_ODE(self):
        n = self._N_layer
        jf = self._feed_layer - 1
        x = self.components.x
        m = len(x)
        imass = self.components.i_mass
        fns = self._fns
        vmax = self._v_max
        vmaxp = self._v_max_p
        rh = self._rh
        rp = self._rp
        X_t = self._X_t
        A = self._A
        hj = self._hj
        Q_s = self._Qs
        
        def dy_dt(t, QC_ins, QC, dQC_ins):
            dQC = np.zeros_like(QC)
            dQC[-(n+1)] = dQC_ins[-1]
            Q_in = QC_ins[-1]
            Q_e = Q_in - Q_s
            C_in = QC_ins[:-1]
            Z_in = C_in*(1-x)
            X_in = sum(C_in*imass*x)           # influent TSS
            if X_in != 0: self._X_comp = C_in * x / X_in     # g COD/g TSS for solids in influent
            X_min = X_in * fns
            X = QC[-n:]                        # (n, ), TSS for each layer
            Z = QC[:m] * (1-x)
            #***********TSS*************
            Q_jout = np.array([Q_e if j < jf else Q_in if j == jf else Q_s for j in range(n)])
            flow_out = X*Q_jout
            flow_in = np.array([Q_e*X[j+1] if j < jf else Q_in*X_in if j == jf else Q_s*X[j-1] for j in range(n)])
            VX = [_settling_flux(xj, vmax, vmaxp, X_min, rh, rp) for xj in X]
            J = [VX[j] if X[j+1] <= X_t and j < jf else min(VX[j], VX[j+1]) for j in range(n-1)]
            settle_out = np.array(J + [0])
            settle_in = np.array([0] + J)
            dQC[-n:] = ((flow_in - flow_out)/A + settle_in - settle_out)/hj        # (n,)
            #*********solubles**********
            dQC[:m] = Q_in*(Z_in - Z)/A/(hj*n)
            return dQC
        
        self._ODE = dy_dt      
    
    def _define_outs(self):
        dct_y = self._state_locator(self._state)
        for out in self.outs:
            Q = dct_y[out.ID][-1]
            Cs = dict(zip(self.components.IDs, dct_y[out.ID][:-1]))
            Cs.pop('H2O', None)
            out.set_flow_by_concentration(Q, Cs, units=('m3/d', 'mg/L'))

    def _design(self):
        pass


class IdealClarifier(SanUnit):
    
    _N_ins = 1
    _N_outs = 2

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 sludge_flow_rate=2000, solids_removal_efficiency=.995,
                 sludge_MLSS=None, isdynamic=False, init_with='WasteStream', 
                 F_BM_default=None, **kwargs):
        
        SanUnit.__init__(self, ID, ins, outs, thermo, isdynamic=isdynamic,
                         init_with=init_with, F_BM_default=F_BM_default)
        self.sludge_flow_rate = sludge_flow_rate
        self.solids_removal_efficiency = solids_removal_efficiency
        self.sludge_MLSS = sludge_MLSS
        
    @property
    def sludge_flow_rate(self):
        '''[float] The designed sludge flow rate (wasted + recycled) in m3/d.'''
        return self._Qs

    @sludge_flow_rate.setter
    def sludge_flow_rate(self, Qs):
        if Qs is not None: self._Qs = Qs
        elif self.ins[0].isempty(): self._Qs = None
        else: self._Qs = self._calc_Qs()

    @property
    def solids_removal_efficiency(self):
        return self._e_rmv
    
    @solids_removal_efficiency.setter
    def solids_removal_efficiency(self, f):
        if f is not None:
            if f > 1 or f < 0:
                raise ValueError(f'solids removal efficiency must be within [0, 1], not {f}')
            self._e_rmv = f
        elif self.ins[0].isempty(): self._e_rmv = None
        else: self._e_rmv = self._calc_ermv()
    
    @property
    def sludge_MLSS(self):
        return self._MLSS
    
    @sludge_MLSS.setter
    def sludge_MLSS(self, MLSS):
        if MLSS is not None: self._MLSS = MLSS
        elif self.ins[0].isempty(): self._MLSS = None
        else: self._MLSS = self._calc_SS()[1]
    
    def _calc_Qs(self, TSS_in=None, Q_in=None):
        if Q_in is None: Q_in = self.ins[0].get_total_flow('m3/d')
        if TSS_in is None: TSS_in = self.ins[0].get_TSS()
        return Q_in*TSS_in*self._e_rmv/(self._MLSS-TSS_in)
    
    def _calc_ermv(self, TSS_in=None, Q_in=None):
        if Q_in is None: Q_in = self.ins[0].get_total_flow('m3/d')
        if TSS_in is None: TSS_in = self.ins[0].get_TSS()
        return self._Qs*(self._MLSS-TSS_in)/TSS_in/(Q_in-self._Qs)        
    
    def _calc_SS(self, SS_in=None, Q_in=None):
        if Q_in is None: Q_in = self.ins[0].get_total_flow('m3/d')
        if SS_in is None: SS_in = self.ins[0].get_TSS()
        SS_e = (1-self._e_rmv)*SS_in
        Qs = self._Qs
        Qe = Q_in - Qs
        return SS_e, (Q_in*SS_in - Qe*SS_e)/Qs
    
    def _run(self):
        inf, = self.ins
        eff, sludge = self.outs
        cmps = self.components
        Q_in = inf.get_total_flow('m3/d')
        TSS_in = (inf.conc*cmps.x*cmps.i_mass).sum()
        params = (Qs, e_rmv, MLSS) = self._Qs, self._e_rmv, self._MLSS
        if sum([i is None for i in params]) > 1: 
            raise RuntimeError('must specify two of the following parameters: '
                               'sludge_flow_rate, solids_removal_efficiency, sludge_MLSS')
        if Qs is None:
            Qs = self._calc_Qs(TSS_in, Q_in)
            Xs = MLSS / TSS_in * inf.conc * cmps.x
            Xe = (1-e_rmv) * inf.conc * cmps.x
        elif e_rmv is None:
            e_rmv = self._calc_ermv(TSS_in, Q_in)
            Xs = MLSS / TSS_in * inf.conc * cmps.x
            Xe = (1-e_rmv) * inf.conc * cmps.x
        else:
            Xe, Xs = self._calc_SS(inf.conc * cmps.x, Q_in)
        Zs = Ze = inf.conc * (1-cmps.x)
        Ce = dict(zip(cmps.IDs, Ze+Xe))
        Cs = dict(zip(cmps.IDs, Zs+Xs))
        Ce.pop('H2O', None)
        Cs.pop('H2O', None)
        eff.set_flow_by_concentration(Q_in-Qs, Ce, units=('m3/d', 'mg/L'))
        sludge.set_flow_by_concentration(Qs, Cs, units=('m3/d', 'mg/L'))
    
    def _design(self):
        pass
        