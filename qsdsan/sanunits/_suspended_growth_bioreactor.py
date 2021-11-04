# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from .. import SanUnit, WasteStream, Process, Processes, CompiledProcesses
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
    
    # cache_state : bool, optional
    #     Whether to store the states of stream composition in the tank from
    #     most recent run. The default is True.

    _N_ins = 3
    _N_outs = 1
    _ins_size_is_fixed = False

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 V_max=1000, aeration=2.0, DO_ID='S_O2', suspended_growth_model=None,
                 isdynamic=True, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, isdynamic=isdynamic)
        self._V_max = V_max
        self._aeration = aeration
        self._DO_ID = DO_ID
        self._model = suspended_growth_model
        self._concs = None
        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def reset_cache(self):
        '''Reset cached states.'''
        self._state = None
        for s in self.outs:
            s.empty()

    @property
    def V_max(self):
        '''[float] The designed maximum liquid volume, not accounting for increased volume due to aeration, in m^3.'''
        return self._V_max

    @V_max.setter
    def V_max(self, Vm):
        self._V_max = Vm

    @property
    def aeration(self):
        '''[:class:`Process` or float or NoneType] Aeration model.'''
        return self._aeration

    @aeration.setter
    def aeration(self, ae):
        if ae is None or isinstance(ae, Process): self._aeration = ae
        elif isinstance(ae, (float, int)):
            if ae < 0:
                raise ValueError('targeted dissolved oxygen concentration for aeration must be non-negative.')
            else:
                if ae > 14:
                    warn(f'targeted dissolved oxygen concentration for {self.ID} might exceed the saturated level.')
                self._aeration = ae
        else:
            raise TypeError(f'aeration must be one of the following types: float, '
                            f'int, Process, NoneType. Not {type(ae)}')

    @property
    def suspended_growth_model(self):
        '''[:class:`CompiledProcesses` or NoneType] Suspended growth model.'''
        return self._model

    @suspended_growth_model.setter
    def suspended_growth_model(self, model):
        if isinstance(model, CompiledProcesses) or model is None: self._model = model
        else: raise TypeError(f'suspended_growth_model must be one of the following '
                              f'types: CompiledProesses, NoneType. Not {type(model)}')

    @property
    def DO_ID(self):
        '''[str] The `Component` ID for dissolved oxygen used in the suspended growth model and the aeration model.'''
        return self._DO_ID

    @DO_ID.setter
    def DO_ID(self, doid):
        if doid not in self.components.IDs:
            raise ValueError(f'DO_ID must be in the set of `CompiledComponents` used to set thermo, '
                             f'i.e., one of {self.components.IDs}.')
        self._DO_ID = doid

    @property
    def state(self):
        '''The state of the CSTR, including component concentrations [mg/L] and flow rate [m^3/d].'''
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

    def set_init_conc(self, **kwargs):
        '''set the initial concentrations [mg/L] of the CSTR.'''
        Cs = np.zeros(len(self.components))
        cmpx = self.components.index
        for k, v in kwargs.items(): Cs[cmpx(k)] = v
        self._concs = Cs

    def _init_state(self, state=None):
        '''initialize state by specifiying or calculating component concentrations
        based on influents. Total flow rate is always initialized as the sum of
        influent wastestream flows.'''
        mixed = WasteStream()
        mixed.mix_from(self.ins)
        Q = mixed.get_total_flow('m3/d')
        if state is not None: Cs = state
        elif self._concs is not None: Cs = self._concs
        else: Cs = mixed.conc
        self._state = np.append(Cs, Q)

    def _state_locator(self, arr):
        '''derives conditions of output stream from conditions within the CSTR'''
        dct = {}
        dct[self.outs[0].ID] = arr
        dct[self.ID] = arr
        return dct

    def _dstate_locator(self, arr):
        '''derives rates of change of output stream from rates of change within the CSTR'''
        return self._state_locator(arr)

    def _load_state(self):
        '''returns a dictionary of values of state variables within the CSTR and in the output stream.'''
        if self._state is None: self._init_state()
        return {self.ID: self._state}

    def _run(self):
        '''Only to converge volumetric flows.'''
        mixed = WasteStream()
        mixed.mix_from(self.ins)
        treated, = self.outs
        treated.copy_like(mixed)

    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE

    def _compile_ODE(self):
        isa = isinstance
        V = self._V_max
        C = list(symbols(self.components.IDs))
        if self._model is None:
            warn(f'{self.ID} was initiated without a suspended growth model, '
                 f'and thus run as a non-reactive unit')
            r = lambda *args: np.zeros(len(C))
        else:
            processes = _add_aeration_to_growth_model(self._aeration, self._model)
            r_eqs = list(processes.production_rates.rate_of_production)
            r = lambdify(C, r_eqs)

        _n_ins = len(self.ins)
        _n_state = len(C) + 1

        if isa(self._aeration, (float, int)):
            i = self.components.index(self._DO_ID)
            fixed_DO = self._aeration
            def dy_dt(t, QC_ins, QC, dQC_ins):
                if _n_ins > 1:
                    QC_ins = QC_ins.reshape((_n_ins, _n_state))
                    Q_ins = QC_ins[:, -1]
                    C_ins = QC_ins[:, :-1]
                    flow_in = Q_ins @ C_ins / V
                    Q_e = Q_ins.sum()
                    dQC_ins = dQC_ins.reshape((_n_ins, _n_state))
                    Q_dot = dQC_ins[:, -1].sum()
                else:
                    Q_ins = QC_ins[-1]
                    C_ins = QC_ins[:-1]
                    flow_in = Q_ins * C_ins / V
                    Q_e = Q_ins
                    Q_dot = dQC_ins[-1]
                Cs = QC[:-1]
                C[i] = fixed_DO
                flow_out = Q_e * Cs / V
                react = np.asarray(r(*Cs))
                C_dot = flow_in - flow_out + react
                C_dot[i] = 0.0
                return np.append(C_dot, Q_dot)
        else:
            def dy_dt(t, QC_ins, QC, dQC_ins):
                if _n_ins > 1:
                    QC_ins = QC_ins.reshape((_n_ins, _n_state))
                    Q_ins = QC_ins[:, -1]
                    C_ins = QC_ins[:, :-1]
                    flow_in = Q_ins @ C_ins / V
                    Q_e = Q_ins.sum()
                    dQC_ins = dQC_ins.reshape((_n_ins, _n_state))
                    Q_dot = dQC_ins[:, -1].sum()
                else:
                    Q_e = Q_ins = QC_ins[-1]
                    C_ins = QC_ins[:-1]
                    flow_in = Q_ins * C_ins / V
                    Q_dot = dQC_ins[-1]
                C = QC[:-1]
                flow_out = Q_e * C / V
                react = np.asarray(r(*C))
                C_dot = flow_in - flow_out + react
                return np.append(C_dot, Q_dot)

        self._ODE = dy_dt

    def _define_outs(self):
        dct_y = self._state_locator(self._state)
        out, = self.outs
        Q = dct_y[out.ID][-1]
        Cs = dict(zip(self.components.IDs, dct_y[out.ID][:-1]))
        Cs.pop('H2O', None)
        out.set_flow_by_concentration(Q, Cs, units=('m3/d', 'mg/L'))

    def _design(self):
        pass


class SBR(SanUnit):
    '''
    Sequential batch reactors operated in parallel. The number of reactors is
    determined by operation cycle and influent flowrate. [1]_

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
        [0, 1]. The default is 2.28e-3.
    cache_state : bool, optional
        Whether to store volume and composition of retained sludge in the tank from
        most recent run. The default is True.

    References
    ----------
    .. [1] Takács, I.; Patry, G. G.; Nolasco, D. A Dynamic Model of the Clarification
        -Thickening Process. Water Res. 1991, 25 (10), 1263–1271.
        https://doi.org/10.1016/0043-1354(91)90066-Y.

    '''

    _N_ins = 1
    _N_outs = 2

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 surface_area=1500, height=4,
                 operation_cycle=(0.5, 1.5, 2.0, 0, 1.0, 0.5, 0.1),
                 aeration=(None, None, None, 2.0), DO_ID='S_O2',
                 suspended_growth_model=None, N_layer=10,
                 pumped_flow=None, underflow=None,
                 X_threshold=3000, v_max=474, v_max_practical=250,
                 rh=5.76e-4, rp=2.86e-3, fns=2.28e-3,
                 cache_state=True, **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with)

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
        self._cache_state = cache_state

        for attr, value in kwargs.items():
            setattr(self, attr, value)
        self._init_Vas = None
        self._init_Cas = None
        self._dynamic_composition = None


    def reset_cache(self):
        '''Reset cached states.'''
        self._init_Vas = self._init_Cas = None
        self._state = None
        for s in self.outs:
            s.empty()

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
        if self._model is None:
            raise RuntimeError(f'{self.ID} was initiated without a suspended growth model.')
        else:
            isa = isinstance
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

            # ********fill and mix/aerate stages***********
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

            # *********settle, decant, desludge**********
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

            if self._cache_state:
                self._init_Vas = V_total - V_eff - V_was
                self._init_Cas = C_as


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