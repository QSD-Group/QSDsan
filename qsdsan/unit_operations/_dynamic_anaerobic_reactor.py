#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
from .. import SanUnit, WasteStream
from ..process_models import T_correction_factor
from ..unit_operations import CSTR
from ..utils import auom, ExogenousDynamicVariable as EDV
from scipy.optimize import newton

__all__ = ('AnaerobicCSTR',)


# %%

class AnaerobicCSTR(CSTR):

    '''
    An anaerobic continuous stirred tank reactor with biogas in headspace. [1]_, [2]_

    Parameters
    ----------
    ins : :class:`WasteStream`
        Influent to the reactor.
    outs : Iterable
        Biogas and treated effluent(s).
    V_liq : float, optional
        Liquid-phase volume [m^3]. The default is 3400.
    V_gas : float, optional
        Headspace volume [m^3]. The default is 300.
    model : :class:`Processes`, optional
        The kinetic model, typically ADM1-like. The default is None.
    T : float, optional
        Operation temperature [K]. The default is 308.15.
    headspace_P : float, optional
        Headspace pressure, if fixed [bar]. The default is 1.013.
    external_P : float, optional
        External pressure, typically atmospheric pressure [bar]. The default is 1.013.
    pipe_resistance : float, optional
        Biogas extraction pipe resistance [m3/d/bar]. The default is 5.0e4.
    fixed_headspace_P : bool, optional
        Whether to assume fixed headspace pressure. The default is False.
    retain_cmps : Iterable[str], optional
        IDs of the components that are assumed to be retained in the reactor, ideally.
        The default is ().
    fraction_retain : float, optional
        The assumed fraction of ideal retention of select components. The default is 0.95.

    References
    ----------
    .. [1] Batstone, D. J.; Keller, J.; Angelidaki, I.; Kalyuzhnyi, S. V;
        Pavlostathis, S. G.; Rozzi, A.; Sanders, W. T. M.; Siegrist, H.;
        Vavilin, V. A. The IWA Anaerobic Digestion Model No 1 (ADM1).
        Water Sci. Technol. 2002, 45 (10), 65–73.
    .. [2] Rosen, C.; Jeppsson, U. Aspects on ADM1 Implementation within
        the BSM2 Framework; Lund, 2006.
    '''

    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False
    _R = 8.3145e-2 # Universal gas constant, [bar/M/K]
    algebraic_h2 = False

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                 init_with='WasteStream', V_liq=3400, V_gas=300, model=None,
                 T=308.15, headspace_P=1.013, external_P=1.013,
                 pipe_resistance=5.0e4, fixed_headspace_P=False,
                 retain_cmps=(), fraction_retain=0.95,
                 isdynamic=True, exogenous_vars=(),
                 pH_ctrl=None, **kwargs):
        if len(exogenous_vars) == 0:
            exogenous_vars = (EDV('T', function=lambda t: T), )
        super().__init__(ID=ID, ins=ins, outs=outs, thermo=thermo,
                         init_with=init_with, V_max=V_liq, aeration=None,
                         DO_ID=None, suspended_growth_model=None,
                         isdynamic=isdynamic, exogenous_vars=exogenous_vars, **kwargs)
        self.V_gas = V_gas
        self.T = T
        # self._S_gas = None
        self._q_gas = 0
        self._n_gas = None
        self._gas_cmp_idx = None
        self._state_keys = None
        self._S_vapor = None
        self.model = model
        self._biogas = WasteStream(phase='g')
        self.headspace_P = headspace_P
        self.external_P = external_P
        self.pipe_resistance = pipe_resistance
        self.fixed_headspace_P = fixed_headspace_P
        self._f_retain = np.array([fraction_retain if cmp.ID in retain_cmps \
                                   else 0 for cmp in self.components])
        self.pH_ctrl = pH_ctrl
        self._mixed = WasteStream()
        self._tempstate = {}

    def ideal_gas_law(self, p=None, S=None):
        '''Calculates partial pressure [bar] given concentration [M] at
        operation temperature or vice versa according to the ideal gas law .'''
        # p in bar, S in M
        if p: return p/self._R/self.T
        elif S: return S*self._R*self.T

    def p_vapor(self, convert_to_bar=True):
        '''Calculates the saturated vapor pressure at operation temperature.'''
        p = self.components.H2O.Psat(self.T)
        if convert_to_bar:
            return p*auom('Pa').conversion_factor('bar')
        else: return p

    @property
    def DO_ID(self):
        '''Not applicable.'''
        return None
    @DO_ID.setter
    def DO_ID(self, doid):
        '''Does nothing.'''
        pass

    @property
    def aeration(self):
        '''Not applicable'''
        return None
    @aeration.setter
    def aeration(self, ae):
        '''Does nothing.'''
        pass

    V_liq = property(CSTR.V_max.fget)
    @V_liq.setter
    def V_liq(self, V):
        '''[float] The liquid-phase volume, in m^3.'''
        CSTR.V_max.fset(self, V)

    model = property(CSTR.suspended_growth_model.fget)
    @model.setter
    def model(self, model):
        '''[:class:`CompiledProcesses` or NoneType] Anaerobic digestion model.'''
        CSTR.suspended_growth_model.fset(self, model)
        if model is not None:
            #!!! how to make unit conversion generalizable to all models?
            self._S_vapor = self.ideal_gas_law(p=self.p_vapor())
            self._n_gas = ng = len(model._biogas_IDs)
            self._state_keys = keys = list(self.components.IDs) \
                + [ID+'_gas' for ID in self.model._biogas_IDs] \
                + ['Q']
            self._gas_cmp_idx = self.components.indices(self.model._biogas_IDs)
            units = ['kg/m3']*len(self.components) + ['M']*ng + ['m3/d']
            self._state_header = [f'{name} [{unit}]' for name, unit in zip(keys, units)]

    split = property(CSTR.split.fget)
    @split.setter
    def split(self, split):
        if split is None: self._split = split
        else:
            if len(split) != len(self._outs)-1:
                raise ValueError('split and outs must have the same size')
            self._split = np.array(split)/sum(split)

    @property
    def headspace_P(self):
        '''Headspace total pressure [bar].'''
        return self._P_gas
    @headspace_P.setter
    def headspace_P(self, P):
        self._P_gas = P

    @property
    def external_P(self):
        '''External (atmospheric) pressure [bar].'''
        return self._P_atm
    @external_P.setter
    def external_P(self, P):
        self._P_atm = P

    @property
    def pipe_resistance(self):
        '''Gas pipe resistance coefficient [m3/d/bar].'''
        return self._k_p
    @pipe_resistance.setter
    def pipe_resistance(self, k):
        self._k_p = k

    @property
    def fixed_headspace_P(self):
        '''Headspace total pressure [bar].'''
        return self._fixed_P_gas
    @fixed_headspace_P.setter
    def fixed_headspace_P(self, b):
        self._fixed_P_gas = bool(b)

    def set_retention_efficacy(self, i):
        if i < 0 or i > 1:
            raise ValueError('retention efficacy must be within [0,1]')
        self._f_retain = (self._f_retain > 0) * i

    @property
    def state(self):
        '''The state of the anaerobic CSTR, including component concentrations [kg/m3],
        biogas concentrations in the headspace [M biogas], and liquid flow rate [m^3/d].'''
        if self._state is None: return None
        else:
            return dict(zip(self._state_keys, self._state))

    @state.setter
    def state(self, arr):
        arr = np.asarray(arr)
        n_state = len(self._state_keys)
        if arr.shape != (n_state, ):
            raise ValueError(f'state must be a 1D array of length {n_state}')
        self._state = arr

    def _run(self):
        '''Only to converge volumetric flows.'''
        mixed = self._mixed # avoid creating multiple new streams
        mixed.mix_from(self.ins)
        # mixed.mix_from(self.ins, energy_balance=False)
        if self.split is None:
            gas, liquid = self.outs
            liquid.copy_like(mixed)
            liquid.T = self.T
            # self._rQ = liquid.F_vol / mixed.F_vol
        else:
            gas = self.outs[0]
            liquids = self._outs[1:]
            Q = mixed.F_vol # m3/hr
            for liquid, spl in zip(liquids, self.split):
                liquid.copy_like(mixed)
                liquid.set_total_flow(Q*spl, 'm3/hr')
                liquid.T = self.T
            # self._rQ = sum([ws.F_vol for ws in liquids]) / mixed.F_vol
        gas.copy_like(self._biogas)
        gas.T = self.T
        if self._fixed_P_gas:
            gas.P = self.headspace_P * auom('bar').conversion_factor('Pa')

    def _init_state(self):
        mixed = self._mixed
        Q = mixed.get_total_flow('m3/d')
        #!!! how to make unit conversion generalizable to all models?
        if self._concs is not None: Cs = self._concs * 1e-3 # mg/L to kg/m3
        else: Cs = mixed.conc * 1e-3 # mg/L to kg/m3
        Gs = [0]*self._n_gas  # initial gas phase concentrations [M]
        Gs[0] = 0.041*0.01
        Gs[1] = 0.041*0.57
        Gs[2] = 0.041*0.4
        self._state = np.append(Cs, Gs + [Q]).astype('float64')
        self._dstate = self._state * 0.

    def _update_state(self):
        y = self._state
        y[-1] = sum(ws.state[-1] for ws in self.ins)
        y[y<1e-16] = 0.
        f_rtn = self._f_retain
        i_mass = self.components.i_mass
        chem_MW = self.components.chem_MW
        n_cmps = len(self.components)
        Cs = y[:n_cmps]*(1-f_rtn)*1e3 # kg/m3 to mg/L
        pH = self.pH_ctrl or self._tempstate.pop('pH', 7)
        if self.split is None:
            gas, liquid = self._outs
            if liquid.state is None:
                liquid.state = np.append(Cs, y[-1])
            else:
                liquid.state[:n_cmps] = Cs
                liquid.state[-1] = y[-1]
            liquid._pH = pH
        else:
            gas = self._outs[0]
            liquids = self._outs[1:]
            for liquid, spl in zip(liquids, self.split):
                if liquid.state is None:
                    liquid.state = np.append(Cs, y[-1]*spl)
                else:
                    liquid.state[:n_cmps] = Cs
                    liquid.state[-1] = y[-1]*spl
                liquid._pH = pH
        if gas.state is None:
            gas.state = np.zeros(n_cmps+1)
        gas.state[self._gas_cmp_idx] = y[n_cmps:(n_cmps + self._n_gas)]
        gas.state[self.components.index('H2O')] = self._S_vapor
        gas.state[-1] = self._q_gas
        gas.state[:n_cmps] = gas.state[:n_cmps] * chem_MW / i_mass * 1e3 # i.e., M biogas to mg (measured_unit) / L

    def _update_dstate(self):
        self._tempstate = self.model.rate_function._params['root'].data.copy()
        dy = self._dstate
        f_rtn = self._f_retain
        n_cmps = len(self.components)
        dCs = dy[:n_cmps]*(1-f_rtn)*1e3
        if self.split is None:
            gas, liquid = self._outs
            if liquid.dstate is None:
                liquid.dstate = np.append(dCs, dy[-1])
            else:
                liquid.dstate[:n_cmps] = dCs
                liquid.dstate[-1] = dy[-1]
        else:
            gas = self._outs[0]
            liquids = self._outs[1:]
            for liquid, spl in zip(liquids, self.split):
                if liquid.dstate is None:
                    liquid.dstate = np.append(dCs, dy[-1]*spl)
                else:
                    liquid.dstate[:n_cmps] = dCs
                    liquid.dstate[-1] = dy[-1]*spl
        if gas.dstate is None:
            # contains no info on dstate
            gas.dstate = np.zeros(n_cmps+1)


    def f_q_gas_fixed_P_headspace(self, rhoTs, S_gas, T):
        cmps = self.components
        gas_mass2mol_conversion = (cmps.i_mass / cmps.chem_MW)[self._gas_cmp_idx]
        self._q_gas = self._R*T/(self._P_gas-self.p_vapor(convert_to_bar=True))\
                                *self.V_liq*sum(rhoTs*gas_mass2mol_conversion)
        return self._q_gas

    def f_q_gas_var_P_headspace(self, rhoTs, S_gas, T):
        p_gas = S_gas * self._R * T
        self._P_gas = P = sum(p_gas) + self.p_vapor(convert_to_bar=True)
        self._q_gas = max(0, self._k_p * (P - self._P_atm))
        return self._q_gas

    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE(self.algebraic_h2, self.pH_ctrl)
        return self._ODE

    def _compile_ODE(self, algebraic_h2=True, pH_ctrl=None):
        if self._model is None:
            CSTR._compile_ODE(self)
        else:
            cmps = self.components
            f_rtn = self._f_retain
            _state = self._state
            _dstate = self._dstate
            _update_dstate = self._update_dstate
            h = None
            if pH_ctrl:
                _params = self.model.rate_function.params
                h = 10**(-pH_ctrl)
                _f_rhos = lambda state_arr: self.model.flex_rhos(state_arr, _params, h=h)
            else:
                _f_rhos = self.model.rate_function
            _f_param = self.model.params_eval
            _M_stoichio = self.model.stoichio_eval
            n_cmps = len(cmps)
            n_gas = self._n_gas
            V_liq = self.V_liq
            V_gas = self.V_gas
            gas_mass2mol_conversion = (cmps.i_mass / cmps.chem_MW)[self._gas_cmp_idx]
            hasexo = bool(len(self._exovars))
            f_exovars = self.eval_exo_dynamic_vars
            if self._fixed_P_gas:
                f_qgas = self.f_q_gas_fixed_P_headspace
            else:
                f_qgas = self.f_q_gas_var_P_headspace
            if self.model._dyn_params:
                def M_stoichio(state_arr):
                    _f_param(state_arr)
                    return self.model.stoichio_eval().T
            else:
                _M_stoichio = self.model.stoichio_eval().T
                M_stoichio = lambda state_arr: _M_stoichio

            h2_idx = cmps.index('S_h2')
            if algebraic_h2:
                params = self.model.rate_function.params
                if self.model._dyn_params:
                    def h2_stoichio(state_arr):
                        return M_stoichio(state_arr)[h2_idx]
                else:
                    _h2_stoichio = _M_stoichio[h2_idx]
                    h2_stoichio = lambda state_arr: _h2_stoichio
                unit_conversion = cmps.i_mass / cmps.chem_MW
                solve_pH = self.model.solve_pH
                dydt_Sh2_AD = self.model.dydt_Sh2_AD
                grad_dydt_Sh2_AD = self.model.grad_dydt_Sh2_AD
                def solve_h2(QC, S_in, T, h=h):
                    if h == None:
                        Ka = params['Ka_base'] * T_correction_factor(params['T_base'], T, params['Ka_dH'])
                        h = solve_pH(QC, Ka, unit_conversion)
                    # S_h2_0 = QC[h2_idx]
                    S_h2_0 = 2.8309E-07
                    S_h2_in = S_in[h2_idx]
                    S_h2 = newton(
                        dydt_Sh2_AD, S_h2_0, grad_dydt_Sh2_AD,
                        args=(QC, h, params, h2_stoichio, V_liq, S_h2_in),
                                  )
                    return S_h2
                def update_h2_dstate(dstate):
                    dstate[h2_idx] = 0.
            else:
                solve_h2 = lambda QC, S_in, T: QC[h2_idx]
                def update_h2_dstate(dstate):
                    pass
            def dy_dt(t, QC_ins, QC, dQC_ins):
                # QC[QC < 0] = 0.
                Q_ins = QC_ins[:, -1]
                S_ins = QC_ins[:, :-1] * 1e-3  # mg/L to kg/m3
                Q = sum(Q_ins)
                S_in = Q_ins @ S_ins / Q
                if hasexo:
                    exo_vars = f_exovars(t)
                    QC = np.append(QC, exo_vars)
                    T = exo_vars[0]
                else: T = self.T
                QC[h2_idx] = _state[h2_idx] = solve_h2(QC, S_in, T)
                rhos =_f_rhos(QC)
                S_liq = QC[:n_cmps]
                S_gas = QC[n_cmps: (n_cmps+n_gas)]
                _dstate[:n_cmps] = (Q_ins @ S_ins - Q*S_liq*(1-f_rtn))/V_liq \
                    + np.dot(M_stoichio(QC), rhos)
                q_gas = f_qgas(rhos[-3:], S_gas, T)
                _dstate[n_cmps: (n_cmps+n_gas)] = - q_gas*S_gas/V_gas \
                    + rhos[-3:] * V_liq/V_gas * gas_mass2mol_conversion
                # _dstate[-1] = dQC_ins[0,-1]
                _dstate[-1] = 0.
                update_h2_dstate(_dstate)
                _update_dstate()

            self._ODE = dy_dt

    def get_retained_mass(self, biomass_IDs):
        cmps = self.components
        mass = cmps.i_mass * self._state[:len(cmps)] * 1e3 # kg/m3 to mg/L
        return self._V_max * mass[cmps.indices(biomass_IDs)].sum()

    def _design(self):
        inf = self.ins[0]
        T_in = inf.T
        T_target = self.T
        if T_target > T_in:
            unit_duty = inf.F_mass * inf.Cp * (T_target - T_in) #kJ/hr
            self.add_heat_utility(unit_duty, T_in, T_out=T_target,
                                  heat_transfer_efficiency=0.8)
