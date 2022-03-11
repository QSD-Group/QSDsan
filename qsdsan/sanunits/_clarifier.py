# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>    

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from numpy import maximum as npmax, minimum as npmin, exp as npexp
from .. import SanUnit, WasteStream
import numpy as np

__all__ = ('FlatBottomCircularClarifier',
           'IdealClarifier',)


def _settling_flux(X, v_max, v_max_practical, X_min, rh, rp, n0):
    X_star = npmax(X-X_min, n0)
    v = npmin(v_max_practical, v_max*(npexp(-rh*X_star) - npexp(-rp*X_star)))
    return X*npmax(v, n0)

# from math import exp
# def _settling_flux(X, v_max, v_max_practical, X_min, rh, rp, n0):
#     X_star = max(X-X_min, 0)
#     v = min(v_max_practical, v_max*(exp(-rh*X_star) - exp(-rp*X_star)))
#     return X*max(v, 0)


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
    underflow : float, optional
        Designed recycling sludge flowrate (RAS), in [m^3/d]. The default is 2000.
    wastage : float, optional
        Designed wasted sludge flowrate (WAS), in [m^3/d]. The default is 385.
    surface_area : float, optional
        Surface area of the clarifier, in [m^2]. The default is 1500.
    height : float, optional
        Height of the clarifier, in [m]. The default is 4.
    N_layer : int, optional
        The number of layers to model settling. The default is 10.
    feed_layer : int, optional
        The feed layer counting from top to bottom. The default is 4.
    X_threshold : float, optional
        Threshold suspended solid concentration, in [g/m^3]. The default is 3000.
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
    _N_outs = 3

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', underflow=2000, wastage=385,
                 surface_area=1500, height=4, N_layer=10, feed_layer=4,
                 X_threshold=3000, v_max=474, v_max_practical=250,
                 rh=5.76e-4, rp=2.86e-3, fns=2.28e-3, isdynamic=True, **kwargs):

        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, isdynamic=isdynamic)
        self._Qras = underflow
        self._Qwas = wastage
        self._sludge = WasteStream()
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
        self._solids = None
        self._solubles = None
        self._X_comp = np.zeros(len(self.components))
        self._dX_comp = self._X_comp.copy()
        header = self._state_header
        self._state_header = list(header) + [f'TSS{i+1} [mg/L]' for i in range(N_layer)]
        for attr, value in kwargs.items():
            setattr(self, attr, value)


    @property
    def underflow(self):
        '''[float] The designed recycling sludge flow rate in m3/d.'''
        return self._Qras

    @underflow.setter
    def underflow(self, ras):
        self._Qras = ras

    @property
    def wastage(self):
        '''[float] The designed wasted sludge flow rate in m3/d.'''
        return self._Qwas

    @wastage.setter
    def wastage(self, was):
        self._Qwas = was

    @property
    def V_settle(self):
        '''[float] Total volume modeled for settling in m^3, calculated based on surface area and height.'''
        return self._V

    @property
    def A_settle(self):
        '''[float] The surface area for settling in m^2, i.e., the area of the clarifier's flat bottom.'''
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
        '''[float] Threshold suspended solid concentration, in g/m^3.'''
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
        self._solubles = Cs*(1-x)

    def set_init_sludge_solids(self, **kwargs):
        '''set the initial concentrations [mg/L] of solids in the underflow sludge.'''
        cmpx = self.components.index
        x = self.components.x
        Xs = np.zeros_like(x)
        for k, v in kwargs.items(): Xs[cmpx(k)] = v
        Xs *= x
        tss = sum(Xs * self.components.i_mass)
        if tss != 0: self._X_comp = Xs / tss

    def set_init_TSS(self, arr):
        '''set the initial TSS [mg/L] in each layer of the clarifier.'''
        if len(arr) != self._N_layer:
            raise ValueError(f'expects an iterable of length {self._N_layer}, not {len(arr)}')
        self._solids = np.asarray(arr, dtype=float)

    def _init_state(self):
        n = self._N_layer
        x = self.components.x
        imass = self.components.i_mass
        QCs = self._ins_QC[0]
        Q = QCs[-1]
        Z = self._solubles if self._solubles is not None \
            else QCs[:-1]*(1-x)
        TSS_in = sum(QCs[:-1] * x * imass)
        TSS = self._solids if self._solids is not None \
            else TSS_in*(20**np.linspace(-1,1,n))
        ZQs = np.append(Z, Q)
        self._state = np.append(ZQs, TSS)
        self._dstate = self._state * 0.
        if TSS_in != 0: self._X_comp = QCs[:-1] * x / TSS_in

    def _update_state(self):
        arr = self._state
        x = self.components.x
        n = self._N_layer
        Q_e = arr[-(1+n)] - self._Qras - self._Qwas
        Z = arr[:len(x)]
        inf, = self.ins
        imass = self.components.i_mass
        C_in = inf.state[:-1]
        X_composition = self._X_comp = C_in*x/sum(C_in*imass*x)
        X_e = arr[-n] * X_composition
        C_s = Z + arr[-1] * X_composition
        eff, ras, was = self._outs
        if eff.isproduct() and eff.state is None: 
            eff.state = np.append(Z+X_e, Q_e)
        else:
            eff.state[:-1] = Z+X_e   # not sure if this works for a setter
            eff.state[-1] = Q_e
        #!!! might need to enable dynamic sludge volume flows
        if ras.isproduct() and ras.state is None: 
            ras.state = np.append(C_s, self._Qras)
        else: 
            ras.state[:-1] = C_s
            ras.state[-1] = self._Qras
        if was.isproduct() and was.state is None: 
            was.state = np.append(C_s, self._Qwas)
        else: 
            was.state[:-1] = C_s
            was.state[-1] = self._Qwas

    def _update_dstate(self):
        arr = self._dstate
        x = self.components.x
        n = self._N_layer
        dQ = arr[-(1+n)]
        dZ = arr[:len(x)]
        TSS_e, TSS_s = self._state[-n], self._state[-1]
        X_composition = self._X_comp # (m, ), mg COD/ mg TSS
        dX_composition = self._dX_comp
        dC_e = dZ + arr[-n] * X_composition + dX_composition * TSS_e
        dC_s = dZ + arr[-1] * X_composition + dX_composition * TSS_s
        eff, ras, was = self._outs
        if eff.isproduct() and eff.dstate is None: 
            eff.dstate = np.append(dC_e, dQ)
        else:
            eff.dstate[:-1] = dC_e # not sure if this works for a setter
            eff.dstate[-1] = dQ
        #!!! might need to enable dynamic sludge volume flows
        if ras.isproduct() and ras.dstate is None: 
            ras.dstate = np.append(dC_s, 0.)
        else: 
            ras.dstate[:-1] = dC_s
        if was.isproduct() and was.dstate is None: 
            was.dstate = np.append(dC_s, 0.)
        else: 
            was.dstate[:-1] = dC_s

    def _run(self):
        '''only to converge volumetric flows.'''
        inf, = self.ins
        sludge = self._sludge
        Q_in = inf.get_total_flow('m3/d')
        eff, ras, was = self.outs
        Q_ras = self._Qras
        Q_was = self._Qwas
        s_e = 1 - (Q_ras+Q_was)/Q_in
        inf.split_to(eff, sludge, s_e)
        sludge.split_to(ras, was, Q_ras/(Q_ras+Q_was))

    def get_retained_mass(self, biomass_IDs):
        cmps = self.components
        tss = self._state[-self._N_layer:].mean()
        mass = cmps.i_mass * self._X_comp * tss
        return self._V * mass[cmps.indices(biomass_IDs)].sum()

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
        Q_s = self._Qras + self._Qwas

        dQC = self._dstate
        dX_comp = self._dX_comp
        _update_dstate = self._update_dstate

        nzeros = np.zeros(n)
        Q_jout = nzeros.copy()
        X_rolled = nzeros.copy()
        X_min_arr = nzeros.copy()
        settle_out = nzeros.copy()
        settle_in = nzeros.copy()

        # Make these constants into arrays so it'll be faster in `dy_dt`
        vmax_arr = np.full_like(nzeros, self._v_max)
        vmaxp_arr = np.full_like(nzeros, self._v_max_p)
        rh_arr = np.full_like(nzeros, self._rh)
        rp_arr = np.full_like(nzeros, self._rp)
        func_vx = lambda x_arr, xmin_arr : _settling_flux(x_arr, vmax_arr, vmaxp_arr, xmin_arr, rh_arr, rp_arr, nzeros)
        
        A, hj, V = self._A, self._hj, self._V
        A_arr = np.full_like(nzeros, A)
        hj_arr = np.full_like(nzeros, hj)
        J = np.zeros(n-1)
        X_t_arr = np.full(jf, self._X_t)
        Q_in_arr = np.zeros(m)
        V_arr = np.full(m, V)

        def dy_dt(t, QC_ins, QC, dQC_ins):
            dQC[-(n+1)] = dQC_ins[0,-1]
            Q_in = QC_ins[0,-1]
            Q_e = Q_in - Q_s
            C_in = QC_ins[0,:-1]
            dC_in = dQC_ins[0,:-1]
            Z_in = C_in*(1-x)
            X_in = sum(C_in*imass*x)           # influent TSS
            dX_in = sum(dC_in*imass*x)
            X_min_arr[:] = X_in * fns
            X = QC[-n:]                        # (n, ), TSS for each layer
            Z = QC[:m] * (1-x)
            #***********TSS*************
            Q_jout[:jf] = Q_e
            Q_jout[jf] = Q_in
            Q_jout[jf+1:] = Q_s
            flow_out = X * Q_jout
            X_rolled[:jf] = X[1: jf+1]
            X_rolled[jf] = X_in
            X_rolled[jf+1:] = X[jf: -1]
            flow_in = X_rolled * Q_jout
            VX = func_vx(X, X_min_arr)
            J[:] = npmin(VX[:-1], VX[1:])
            condition = (X_rolled[:jf]<X_t_arr)
            J[:jf][condition] = VX[:jf][condition]
            settle_out[:-1] = J
            settle_in[1:] = J
            dQC[-n:] = ((flow_in - flow_out)/A_arr + settle_in - settle_out)/hj_arr       # (n,)
            #*********solubles**********
            Q_in_arr[:] = Q_in
            dQC[:m] = Q_in_arr*(Z_in - Z)/V_arr
            # instrumental variables
            dX_comp[:] = (dC_in * X_in - dX_in * C_in) * x / X_in**2
            _update_dstate()

        self._ODE = dy_dt

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