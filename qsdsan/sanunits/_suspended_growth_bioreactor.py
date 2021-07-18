# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Cheung <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from .. import SanUnit, Construction, WasteStream, Process, Processes
from sympy import symbols, lambdify, Matrix
from scipy.integrate import solve_ivp
from warnings import warn
import numpy as np

__all__ = ('CSTR', 
           'SBRs', 
           # 'PFR',
           )

class CSTR(SanUnit):

    _N_ins = 2
    _N_outs = 1
    
    def __init__(self, ID='', ins=None, outs=(), V_max=1000, aeration=2.0, DO_ID='S_O2',
                 suspended_growth_model=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs)
        self._V_max = V_max
        self._aeration = aeration
        self._DO_ID = DO_ID
        self._model = suspended_growth_model
        for attr, value in kwargs.items():
            setattr(self, attr, value)
        self._init_C = None
        
    def _run(self, t_bound=10, steady_state=False, cache_state=True):
        isa = isinstance        
        if self._model:
            mixed = WasteStream()
            mixed.mix_from(self.ins)
            treated = self.outs[0]
            treated.copy_like(mixed)
            Q = mixed.get_total_flow('m3/d')
            tau = self._V_max / Q
            C_in = mixed.mass / mixed.F_vol * 1e3    # concentrations in g/m3
            if self._init_C: C_0 = self._init_C
            else: C_0 = C_in            
            C = list(symbols(self.components.IDs))
            
            if isa(self._aeration, Process):
                processes = Processes(self._model.tuple)
                processes.append(self._aeration)
                processes.compile()
            else:
                processes = self._model                
                
            mass_balance_terms = list(zip(C_in, C, processes.production_rates.rate_of_production))
            C_dot_eqs = [(cin-c)/tau + r for cin, c, r in mass_balance_terms]
            
            if isa(self._aeration, (float, int)):
                i = self.components.index(self._DO_ID)
                C_0[i] = C_in[i] = self._aeration
                C_dot_eqs[i] = 0
            
            def dC_dt(t, y):
                C_dot = lambdify(C, C_dot_eqs)
                return C_dot(*y)
            def limit(t, y):
                dCdt = np.array(dC_dt(0, y))
                if np.allclose(dCdt, np.zeros(len(y)), atol=1e-3): return 0
                else: return 1
            limit.terminal = True            
            J = Matrix(dC_dt(None, C)).jacobian(C)
            def J_func(t, y):
                J_func = lambdify(C, J)
                return J_func(*y)            
            
            sol = solve_ivp(dC_dt, (0, t_bound), C_0, method='BDF', jac=J_func, events=limit)
            C_out = sol.y.transpose()[-1]
            if steady_state:
                while len(sol.t_events) == 0:
                    sol = solve_ivp(dC_dt, (0, t_bound), C_out, method='BDF', jac=J_func, events=limit)
                    C_out = sol.y.transpose()[-1]
            else:
                if len(sol.t_events) == 0: 
                    warn(f'{self.ID} did not reach steady state in this run.')
            if cache_state: self._init_C = C_out
            treated.set_flow(C_out*treated.F_vol, 'g/hr', self.components.IDs)
            
        else:
            raise RuntimeError(f'{self.ID} was initiated without a suspended growth model.')
    
    def _design(self):
        pass
    

class SBRs(SanUnit):
    
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), **kwargs):
        SanUnit.__init__(self, ID, ins, outs)
    
    def _run(self, steady_state=True):
        pass
    
    def _design(self):
        pass


# class PFR(SanUnit):
    
#     _N_ins = 1
#     _N_outs = 2
    
#     def __init__(self, ID='', ins=None, outs=(), **kwargs):
#         SanUnit.__init__(self, ID, ins, outs)
    
#     def _run(self, steady_state=True):
#         pass
    
#     def _design(self):
#         pass
    