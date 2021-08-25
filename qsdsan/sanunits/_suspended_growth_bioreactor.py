# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Cheung <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from .. import SanUnit, WasteStream, Process, Processes
from ._clarifier import _settling_flux
from sympy import symbols, lambdify, Matrix
from scipy.integrate import solve_ivp
from warnings import warn
from math import floor, ceil
import numpy as np
import pandas as pd

__all__ = ('CSTR', 
           'SBR', 
           # 'PFR',
           )

def _add_aeration_to_growth_model(aer, model):
    if isinstance(aer, Process):
        processes = Processes(model.tuple)
        processes.append(aer)
        processes.compile()
        processes.set_parameters(**aer.parameters, **model.parameters)
    else:
        processes = model
        processes.compile()
    return processes

#%%
class CSTR(SanUnit):
    '''
    A single continuous stirred tank reactor.

    Parameters
    ----------
    ID : str
        ID for the reactor.
    ins : :class:`WasteStream`
        Influents to the reactor. Can be an array of up to 3 WasteStream objects by
        default, typically wastewater to be treated, recycled effluent, recycled 
        activated sludge.
    outs : :class:`WasteStream`
        Treated effluent.
    V_max : float
        Designed volume, in [m^3]. The default is 1000.
    aeration : float or :class:`Process`, optional
        Aeration setting. Either specify a targeted dissolved oxygen concentration
        in [mg O2/L] or provide a :class:`Process` object to represent aeration,
        or None for no aeration. The default is 2.0.
    DO_ID : str, optional
        The :class:`Component` ID for dissolved oxygen, only relevant when the 
        reactor is aerated. The default is 'S_O2'.
    suspended_growth_model : :class:`Processes`, optional
        The suspended growth biokinetic model. The default is None.

    '''

    _N_ins = 3
    _N_outs = 1
    _ins_size_is_fixed = False
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', 
                 V_max=1000, aeration=2.0, DO_ID='S_O2',
                 suspended_growth_model=None, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)
        self._V_max = V_max
        self._aeration = aeration
        self._DO_ID = DO_ID
        self._model = suspended_growth_model
        for attr, value in kwargs.items():
            setattr(self, attr, value)
        self._init_C = None
        
    def _run(self, t_bound=10, steady_state=False, cache_state=True):
        isa = isinstance        
        if self._model is not None:
            mixed = WasteStream()
            mixed.mix_from(self.ins)
            treated = self.outs[0]
            treated.copy_like(mixed)
            Q = mixed.get_total_flow('m3/d')
            tau = self._V_max / Q
            C_in = mixed.mass / mixed.F_vol * 1e3    # concentrations in g/m3
            if self._init_C is not None: C_0 = self._init_C
            else: C_0 = C_in            
            C = list(symbols(self.components.IDs))
            
            processes = _add_aeration_to_growth_model(self._aeration, self._model)                               
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
            C_out = np.array([max(y, 0) for y in sol.y.transpose()[-1]])
            if steady_state:
                while len(sol.t_events) == 0:
                    sol = solve_ivp(dC_dt, (0, t_bound), C_out, method='BDF', jac=J_func, events=limit)
                    C_out = np.array([max(y, 0) for y in sol.y.transpose()[-1]])
            else:
                if len(sol.t_events) == 0: 
                    warn(f'{self.ID} did not reach steady state in this run.')
            if cache_state: self._init_C = C_out
            treated.set_flow(C_out*treated.F_vol, 'g/hr', self.components.IDs)
        else:
            raise RuntimeError(f'{self.ID} was initiated without a suspended growth model.')
    
    def _design(self):
        pass
    
    @property
    def init_state(self):
        return dict(zip(self.components.IDs, self._init_C))
    
    @init_state.setter
    def init_state(self, C):
        if len(C) == len(self.components.IDs):
            self._init_C = np.asarray(C)
        else: 
            raise ValueError(f'Must be a 1D array of length {len(self.components.IDs)}')

class SBR(SanUnit):
    '''
    Sequential batch reactors operated in parallel. The number of reactors is 
    determined by operation cycle and influent flowrate. 

    Parameters
    ----------
    ID : str
        ID for the reactors. The default is ''.
    ins : :class:`WasteStream`
        Influent to the reactor. Expected number of influent is 1.
    outs : :class:`WasteStream`
        Treated effluent and wasted sludge.
    surface_area : float, optional
        Surface area of the reactor bottom, in [m^2]. The reactor is assumed 
        to be cylinder. The default is 1500.
    height : float, optional
        Height of the reactor, in [m]. The default is 4.
    operation_cycle : iterable of float, optional
        Operation cycle of the SBR, time for each stage specified in [h]. There
        are 7 stages: 1 - fill, 2 - fill, 3 - mix, 4 - mix, 5 - settle, 6 - decant, 
        7 - desludge. The first 4 stages are modeled as a biological reactor. 
        The 5th stage is modeled as a 1D N-layer settler. The last 2 stages are 
        assumed inactive. The default is (0.5, 1.5, 2.0, 0, 1.0, 0.5, 0.1).
    aeration : iterable of float and/or :class:`Process`, optional
        Aeration settings for the first 4 stages of the operation cycle. Either 
        specify a targeted dissolved oxygen concentration in [mg O2/L] or provide 
        a :class:`Process` object to represent aeration, or None for no aeration. 
        The default is (None, None, None, 2.0).
    DO_ID : str, optional
        The :class:`Component` ID for dissolved oxygen, only relevant when the 
        reactor is aerated. The default is 'S_O2'.
    suspended_growth_model : :class:`Processes`, optional
        The suspended growth biokinetic model. The default is None.
    N_layer : int, optional
        The number of layers to model settling. The default is 10.
    pumped_flow : float, optional
        Designed effluent flowrate, in [m^3/d]. The default is None.
    underflow : float, optional
        Designed wasted activated sludge flowrate, in [m^3/d]. The default is None.
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
        [0, 1], The default is 2.28e-3.

    References
    ----------
    .. [1] Takács, I.; Patry, G. G.; Nolasco, D. A Dynamic Model of the Clarification
        -Thickening Process. Water Res. 1991, 25 (10), 1263–1271. 
        https://doi.org/10.1016/0043-1354(91)90066-Y.

    '''
        
    _N_ins = 1
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), surface_area=1500, height=4, 
                 operation_cycle=(0.5, 1.5, 2.0, 0, 1.0, 0.5, 0.1), 
                 aeration=(None, None, None, 2.0), DO_ID='S_O2',
                 suspended_growth_model=None, N_layer=10, 
                 pumped_flow=None, underflow=None,
                 X_threshold=3000, v_max=474, v_max_practical=250, 
                 rh=5.76e-4, rp=2.86e-3, fns=2.28e-3, **kwargs):
        SanUnit.__init__(self, ID, ins, outs)
        self._V = surface_area * height
        self._A = surface_area
        self._h = height
        self._operation_cycle = operation_cycle
        self._aeration = aeration
        self._DO_ID = DO_ID
        self._model = suspended_growth_model
        self._N_layer = N_layer
        self._Q_e = pumped_flow
        self._Q_WAS = underflow
        self._X_t = X_threshold
        self._v_max = v_max
        self._v_max_p = v_max_practical
        self._rh = rh
        self._rp = rp
        self._fns = fns
        for attr, value in kwargs.items():
            setattr(self, attr, value)
        self._init_Vas = None
        self._init_Cas = None
        self._dynamic_composition = None
        
    @property
    def operation_cycle(self):
        return dict(zip(('fill_1', 'fill_2', 'mix_1', 'mix_2', 'settle', 'decant', 'desludge'),
                        self._operation_cycle))
    @property
    def total_cycle_time(self):
        return sum(self._operation_cycle)
    
    @property
    def aeration(self):
        return dict(zip(('fill_1', 'fill_2', 'mix_1', 'mix_2'),
                        self._aeration[:4]))
    
    @property
    def C_t(self):
        if self._dynamic_composition:
            return pd.DataFrame(self._dynamic_composition, 
                                columns = ['Time[d]'] + list(self.components.IDs))
        else: return None
    
    def _run(self, cache_state=True):
        isa = isinstance        
        if self._model is not None:
            inf = self.ins[0]
            Q_in = inf.get_total_flow('m3/d')
            eff, sludge = self.outs
            eff.copy_like(inf)
            sludge.copy_like(inf)
            C_in = inf.mass / inf.F_vol * 1e3    # concentrations in g/m3
            cmps = self.components
            C = list(symbols(cmps.IDs))
            if self._init_Vas is not None:
                V_0 = self._init_Vas
                C_0 = self._init_Cas
            else: 
                V_0 = 0
                C_0 = C_in
            n = self._N_layer
            if self._aeration.count(None) == len(self._aeration): 
                Vmax = self._V
                hj = self._h/n
            else: 
                Vmax = self._V*0.75
                hj = self._h*0.75/n
            
            # fill and mix/aerate stages
            T_fill = (Vmax - V_0)/Q_in # maximum total fill time in day
            T = [t/24 for t in self._operation_cycle]  # operation cycle in day
            if T_fill <= T[0]: 
                schedule = [T_fill, T[0]-T_fill] + T[1:4]
                aer = [self._aeration[0], self._aeration[0]] + list(self._aeration[1:4])
                fill = [True] + [False]*4
                V_total = Vmax
            elif T_fill <= T[0]+T[1]:
                schedule = [T[0], T_fill-T[0], T[0]+T[1]-T_fill] + T[2:4]
                aer = list(self._aeration[:2]) + [self._aeration[1]] + list(self._aeration[2:4])
                fill = [True]*2 + [False]*3
                V_total = Vmax
            else:
                schedule = T[:4]
                aer = list(self._aeration[:4])
                fill = [True]*2 + [False]*2
                V_total = Q_in*(T[0]+T[1])+V_0
                hj = V_total/self._V*self._h/n
            
            for i in range(1, len(schedule)):
                if fill[-i] == fill[-i-1] and aer[-i] == aer[-i-1]:
                    schedule[-i-1] += schedule[-i]
                    schedule[-i] = 0
            
            t_arr = np.array([])
            y_mat = np.ndarray([])
            for i in range(len(schedule)):
                if schedule[i] > 0:
                    dC_dt, J_func = self._compile_dC_dt(V_0, Q_in, C_in, C, fill[i], aer[i])
                    if isa(aer[i], (float, int)): C_0[cmps.index(self._DO_ID)] = aer[i]
                    sol = solve_ivp(dC_dt, (0, schedule[i]), C_0, method='BDF', jac=J_func)
                    C_0 = sol.y.transpose()[-1]
                    V_0 += Q_in * schedule[i] * fill[i]
                    t_arr = np.concatenate((t_arr, sol.t + t_arr[-1]))
                    y_mat = np.hstack((y_mat, sol.y))
            self._dynamic_composition = np.vstack((t_arr, y_mat)).transpose()
            
            # settle, decant, desludge
            eff.set_flow(C_0*eff.F_vol, 'g/hr', self.components.IDs)
            X_0 = eff.get_TSS()
            X_min = X_0 * self._fns
            T_settle = T[4]
            def dX_dt(t, X):         
                VX = [_settling_flux(x, self._v_max, self._v_max_p, X_min, self._rh, self._rp) for x in X]
                J = [VX[j] if X[j+1] <= self._X_t else min(VX[j], VX[j+1]) for j in range(n-1)]
                settle_out = np.array(J + [0])
                settle_in = np.array([0] + J)
                dXdt = (settle_in - settle_out)/hj
                return dXdt
            sol = solve_ivp(dX_dt, (0, T_settle), np.ones(n)*X_0)
            X = sol.y.transpose()[-1]
            
            V_eff = min(T[5]*self._Q_e, V_total*(n-1)/n)
            n_eff = V_eff/V_total
            w_eff = np.array([1]*floor(n_eff)+[n_eff-floor(n_eff)])
            X_eff = np.average(X[:ceil(n_eff)], weights=w_eff)
            eff_mass_flow = (X_eff/X_0*cmps.x + (1-cmps.x))*C_0*V_eff/T[5]
            eff.set_flow(eff_mass_flow, 'g/d', cmps.IDs)
            
            V_was = min(T[6]*self._Q_WAS, V_total-V_eff)
            X_as = (V_total*X_0 - V_eff*X_eff) / (V_total-V_eff)
            C_as = (X_as/X_0*cmps.x + (1-cmps.x))*C_0
            was_mass_flow = C_as*V_was/T[6]
            sludge.set_flow(was_mass_flow, 'g/d', cmps.IDs)
            
            if cache_state: 
                self._init_Vas = V_total - V_eff - V_was
                self._init_Cas = C_as
        else: 
            raise RuntimeError(f'{self.ID} was initiated without a suspended growth model.')

    def _design(self):
        pass
    
    def _compile_dC_dt(self, V0, Qin, Cin, C, fill, aer):
        isa = isinstance
        processes = _add_aeration_to_growth_model(aer, self._model)
        if fill:
            t = symbols('t')
            mass_balance_terms = list(zip(Cin, C, processes.production_rates.rate_of_production))
            C_dot_eqs = [(cin-c)/(t+V0/Qin) + r for cin, c, r in mass_balance_terms]
            if isa(aer, (float, int)): C_dot_eqs[self.components.index(self._DO_ID)] = 0
            def dC_dt(t, y):
                C_dot = lambdify([t]+C, C_dot_eqs)
                return C_dot(t, *y)
            J = Matrix(dC_dt(t, C)).jacobian(C)
        else:
            C_dot_eqs = processes.production_rates.rate_of_production
            if isa(aer, (float, int)): C_dot_eqs[self.components.index(self._DO_ID)] = 0
            def dC_dt(t, y):
                C_dot = lambdify(C, C_dot_eqs)
                return C_dot(*y)
            J = Matrix(dC_dt(None, C)).jacobian(C)
        def J_func(t, y):
            J_func = lambdify(C, J)
            return J_func(*y)
        return (dC_dt, J_func)

# class PFR(SanUnit):
    
#     _N_ins = 1
#     _N_outs = 2
    
#     def __init__(self, ID='', ins=None, outs=(), **kwargs):
#         SanUnit.__init__(self, ID, ins, outs)
    
#     def _run(self, steady_state=True):
#         pass
    
#     def _design(self):
#         pass
    