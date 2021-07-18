# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Cheung <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.


Reference
---------
[1] Takács, I.; Patry, G. G.; Nolasco, D. A Dynamic Model of the 
Clarification-Thickening Process. Water Res. 1991, 25 (10), 1263–1271. 
https://doi.org/10.1016/0043-1354(91)90066-Y.

'''

__all__ = ('FlatBottomCircularClarifier', )

from .. import SanUnit
from math import exp
from sympy import symbols, lambdify, Matrix
from scipy.integrate import solve_ivp
import numpy as np

def settling_flux(X, v_max, v_max_practical, X_min, rh, rp):
    X_star = max(X-X_min, 0)
    v = min(v_max_practical, v_max*(exp(-rh*X_star) - exp(-rp*X_star)))
    return X*max(v, 0)

class FlatBottomCircularClarifier(SanUnit):
    
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), sludge_flow_rate=2000, 
                 surface_area=1500, height=4, N_layer=10, feed_layer=4, 
                 X_threshold=3000, v_max=474, v_max_practical=250, 
                 rh=5.76e-4, rp=2.86e-3, fns=2.28e-3, **kwargs):
        SanUnit.__init__(self, ID, ins, outs)
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
        for attr, value in kwargs.items():
            setattr(self, attr, value)
        self._init_X = None   #TSS
        self._init_Z = None   #TDS
        
    def _run(self, t_bound=10, steady_state=False, cache_state=True):
        inf = self.ins[0]
        X_in = inf.get_TSS()
        Q_in = inf.get_total_flow('m3/d')
        eff, sludge = self.outs
        eff.copy_like(inf)
        sludge.copy_like(inf)
        Q_s = self._Qs
        Q_e = Q_in - Q_s
        jf = self._feed_layer - 1
        if jf not in range(self._N_layer): 
            raise ValueError(f'feed layer {self._feed_layer} is out of range.'
                             f'must be an integer between 1 and {self._N_layer}.')
        if self._init_X: X_0 = self._init_X
        else: X_0 = X_in  
        X_min = X_in * self._fns
        n = self._N_layer
        def dX_dt(t, X):
            flow_out = X*Q_e*[j <= jf for j in range(n)] + X*Q_s*[j >= jf for j in range(n)]
            flow_in = np.array([flow_out[j+1] if j < jf else flow_out[j-1] for j in range(n)])
            flow_in[jf] = Q_in * X_in            
            VX = [settling_flux(x, self._v_max, self._v_max_p, X_min, self._rh, self._rp) for x in X]
            J = [VX[j] if X[j+1] <= X_t and j < jf else min(VX[j], VX[j+1]) for j in range(n-1)]
            settle_out = np.array(J + [0])
            settle_in = np.array([0] + J)
            dXdt = ((flow_in - flow_out)/self._A + settle_in - settle_out)/self._hj
            return dXdt
        def limit(t, X):
            dXdt = dX_dt(0, X)
            if np.allclose(dXdt, np.zeros(len(X)), atol=1e-2): return 0
            else: return 1
        limit.terminal = True
        
        sol = solve_ivp(dX_dt, (0, t_bound), np.ones(n)*X_0, method='BDF', events=limit)
        X_out = sol.y.transpose()[-1]
        if steady_state:
            while len(sol.t_events) == 0:
                sol = solve_ivp(dX_dt, (0, t_bound), X_out, method='BDF', events=limit)
                X_out = sol.y.transpose()[-1]
        else:
            if len(sol.t_events) == 0: 
                warn(f'{self.ID} did not reach steady state in this run.')
        if cache_state: self._init_X = X_out
        
        cmps = self.components
        TSS_ratio_eff_in = X_out[1]*Q_e/(X_in*Q_in)
        eff_mass_flow = (TSS_ratio_eff_in*cmps.x + Q_e/Q_in*(1-cmps.x))*inf.mass
        eff.set_flow(eff_mass_flow, 'kg/hr', cmps.IDs)
        sludge.separate_out(eff)
    
    def _design(self):
        pass