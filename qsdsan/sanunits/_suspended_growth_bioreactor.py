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
from ..sanunits import dydt_cstr
from sympy import symbols, lambdify, Matrix
from scipy.integrate import solve_ivp
from warnings import warn
from math import floor, ceil
import numpy as np
import pandas as pd
from numba import njit

__all__ = ('CSTR',
           'BatchExperiment',
           'SBR',
           'PFR',
           )

def _add_aeration_to_growth_model(aer, model):
    if isinstance(aer, Process):
        processes = Processes(model.tuple)
        processes.append(aer)
        processes.compile()
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
    split : iterable of float
        Volumetric splits of effluent flows if there are more than one effluent.
        The default is None.
    V_max : float
        Designed volume, in [m^3]. The default is 1000.
        
    # W_tank : float
    #     The design width of the tank, in [m]. The default is 6.4 m (21 ft). [1, Yalin's adaptation of code]
    # D_tank : float
    #     The design depth of the tank in [m]. The default is 3.65 m (12 ft). [1, Yalin's adaptation of code]
    # freeboard : float
    #     Freeboard added to the depth of the reactor tank, [m]. The default is 0.61 m (2 ft). [1, Yalin's adaptation of code]
        
    aeration : float or :class:`Process`, optional
        Aeration setting. Either specify a targeted dissolved oxygen concentration
        in [mg O2/L] or provide a :class:`Process` object to represent aeration,
        or None for no aeration. The default is 2.0.
    DO_ID : str, optional
        The :class:`Component` ID for dissolved oxygen, only relevant when the
        reactor is aerated. The default is 'S_O2'.
    suspended_growth_model : :class:`Processes`, optional
        The suspended growth biokinetic model. The default is None.
    exogenous_var : iterable[:class:`ExogenousDynamicVariable`], optional
        Any exogenous dynamic variables that affect the process mass balance,
        e.g., temperature, sunlight irradiance. Must be independent of state 
        variables of the suspended_growth_model (if has one).
    
    References:
        
     [1] Shoener, B. D.; Zhong, C.; Greiner, A. D.; Khunjar, W. O.; Hong, P.-Y.; Guest, J. S.
         Design of Anaerobic Membrane Bioreactors for the Valorization
         of Dilute Organic Carbon Waste Streams.
         Energy Environ. Sci. 2016, 9 (3), 1102–1112.
         https://doi.org/10.1039/C5EE03715H.
    
    '''
    _N_ins = 3
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = False

    def __init__(self, ID='', ins=None, outs=(), split=None, thermo=None,
                 init_with='WasteStream', V_max=1000, W_tank = 6.4, D_tank = 3.65,
                 freeboard = 0.61, t_wall = None, t_slab = None, aeration=2.0, 
                 DO_ID='S_O2', suspended_growth_model=None, isdynamic=True, exogenous_vars=(), **kwargs):
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, isdynamic=isdynamic,
                         exogenous_vars=exogenous_vars, **kwargs)
        self._V_max = V_max
        self._aeration = aeration
        self._DO_ID = DO_ID
        self._model = suspended_growth_model
        self._concs = None
        self._mixed = WasteStream()
        self.split = split

        # # Design parameters 
        # self._W_tank = W_tank
        # self._D_tank = D_tank
        # self._freeboard = freeboard
        # self._t_wall = t_wall
        # self._t_slab = t_slab
    
        # for attr, value in kwargs.items():
        #     setattr(self, attr, value)
    
    @property
    def V_max(self):
        '''[float] The designed maximum liquid volume, not accounting for increased volume due to aeration, in m^3.'''
        return self._V_max

    @V_max.setter
    def V_max(self, Vm):
        self._V_max = Vm
        
    # @property
    # def W_tank(self):
    #     '''[float] The design width of the tank, in m.'''
    #     return self._W_tank

    # @W_tank.setter
    # def W_tank(self, W_tank):
    #     self._W_tank = W_tank
        
    # @property
    # def D_tank(self):
    #     '''[float] The design depth of the tank, in m.'''
    #     return self._D_tank

    # @D_tank.setter
    # def D_tank(self, D_tank):
    #     self._D_tank = D_tank
        
    # @property
    # def freeboard(self):
    #     '''[float] Freeboard added to the depth of the reactor tank, [m].'''
    #     return self._freeboard
    
    # @freeboard.setter
    # def freeboard(self, i):
    #     self._freeboard = i
        
    # @property
    # def t_wall(self):
    #     '''
    #     [float] Thickness of the wall concrete, [m].
    #     default to be minimum of 1 ft with 1 in added for every ft of depth over 12 ft.
    #     '''
    #     D_tank = self.D_tank*39.37 # m to inches 
    #     return self._t_wall or (1 + max(D_tank - 12, 0)/12)*0.3048 # from feet to m
    
    # @t_wall.setter
    # def t_wall(self, i):
    #     self._t_wall = i

    # @property
    # def t_slab(self):
    #     '''
    #     [float] Concrete slab thickness, [m],
    #     default to be 2 in thicker than the wall thickness.
    #     '''
    #     return self._t_slab or (self.t_wall + 2/12)*0.3048 # from feet to m
    
    # @t_slab.setter
    # def t_slab(self, i):
    #     self._t_slab = i
     
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
    def split(self):
        '''[numpy.1darray or NoneType] The volumetric split of outs.'''
        return self._split

    @split.setter
    def split(self, split):
        if split is None: self._split = split
        else:
            if len(split) != len(self._outs):
                raise ValueError('split and outs must have the same size')
            self._split = np.array(split)/sum(split)

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
        self._concs = self.components.kwarray(kwargs)

    def _init_state(self):
        mixed = self._mixed
        Q = mixed.get_total_flow('m3/d')
        if self._concs is not None: Cs = self._concs
        else: Cs = mixed.conc
        self._state = np.append(Cs, Q).astype('float64')
        self._dstate = self._state * 0.

    def _update_state(self):
        arr = self._state
        arr[-1] = sum(ws.state[-1] for ws in self.ins)
        if self.split is None: self._outs[0].state = arr
        else:
            for ws, spl in zip(self._outs, self.split):
                y = arr.copy()
                y[-1] *= spl
                ws.state = y

    def _update_dstate(self):
        arr = self._dstate
        if self.split is None: self._outs[0].dstate = arr
        else:
            for ws, spl in zip(self._outs, self.split):
                y = arr.copy()
                y[-1] *= spl
                ws.dstate = y

    def _run(self):
        '''Only to converge volumetric flows.'''
        mixed = self._mixed # avoid creating multiple new streams
        mixed.mix_from(self.ins)
        Q = mixed.F_vol # m3/hr
        if self.split is None: self.outs[0].copy_like(mixed)
        else:
            for ws, spl in zip(self._outs, self.split):
                ws.copy_like(mixed)
                ws.set_total_flow(Q*spl, 'm3/hr')

    def get_retained_mass(self, biomass_IDs):
        cmps = self.components
        mass = cmps.i_mass * self._state[:len(cmps)]
        return self._V_max * mass[cmps.indices(biomass_IDs)].sum()

    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE

    def _compile_ODE(self):
        isa = isinstance
        C = list(symbols(self.components.IDs))
        m = len(C)
        if self._model is None:
            warn(f'{self.ID} was initialized without a suspended growth model, '
                 f'and thus run as a non-reactive unit')
            r = lambda state_arr: np.zeros(m)
            
        else:
            processes = _add_aeration_to_growth_model(self._aeration, self._model)
            r = processes.production_rates_eval

        _dstate = self._dstate
        _update_dstate = self._update_dstate
        V = self._V_max
        hasexo = bool(len(self._exovars))
        f_exovars = self.eval_exo_dynamic_vars
        
        if isa(self._aeration, (float, int)):
            i = self.components.index(self._DO_ID)
            fixed_DO = self._aeration
            def dy_dt(t, QC_ins, QC, dQC_ins):
                QC[i] = fixed_DO
                dydt_cstr(QC_ins, QC, V, _dstate)
                if hasexo: QC = np.append(QC, f_exovars(t))
                _dstate[:-1] += r(QC)
                _dstate[i] = 0
                _update_dstate()
        else:
            def dy_dt(t, QC_ins, QC, dQC_ins):
                dydt_cstr(QC_ins, QC, V, _dstate)
                if hasexo: QC = np.append(QC, f_exovars(t))
                _dstate[:-1] += r(QC)
                _update_dstate()

        self._ODE = dy_dt
        
    # _units = {
    #     'Tank volume': 'm3',
    #     'Tank width': 'm',
    #     'Tank depth': 'm',
    #     'Tank length': 'm',
    #     'Volume of concrete wall': 'm3',
    #     'Volume of concrete slab': 'm3' 
    #     }

    def _design(self):
        pass
        # self._mixed.mix_from(self.ins)
        # # mixed = self._mixed
        # D = self.design_results
        
        # D['Tank volume'] = V = self.V_max
        # D['Tank width'] = W = self.W_tank
        # D['Tank depth'] = depth = self.D_tank
        # D['Tank length'] = L = V/(W*depth)
        
        # t_wall, t_slab = self.t_wall, self.t_slab
        # t = t_wall + t_slab
        # D_tot = depth + self.freeboard 
        
        # # get volume of wall concrete
        # VWC = 2*((L + 2*t_wall)*t_wall*D_tot) + 2*(W*t_wall*D_tot)
        
        # # get volume of slab concrete
        # VSC = (L + 2*t_wall)*(W + 2*t_wall)*t
        
        # D['Volume of concrete wall'] = VWC
        # D['Volume of concrete slab'] = VSC
            
    def _cost(self):
        pass
        # self._mixed.mix_from(self.ins)
       
        # D = self.design_results
        # C = self.baseline_purchase_costs
       
        # # Construction of concrete and stainless steel walls
        # C['Wall concrete'] = D['Volume of concrete wall']*self.wall_concrete_unit_cost
        # C['Slab concrete'] = D['Volume of concrete slab']*self.slab_concrete_unit_cost


#%%
class BatchExperiment(SanUnit):
    
    _N_ins = 0
    _N_outs = 0
    # _ins_size_is_fixed = True
    # _outs_size_is_fixed = True

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream', 
                 model=None, isdynamic=True, exogenous_vars=(), **kwargs):
        '''
        A batch reactor in experimental settings.

        Parameters
        ----------
        model : :class:`CompiledProcesses`, optional
            Process model that describes the dynamics of state variables. 
            The `state` of the batch reactor is entirely determined by the 
            stoichiometry and rate function in this model.
        '''
        SanUnit.__init__(self, ID, None, (), thermo, init_with, isdynamic=isdynamic,
                         exogenous_vars=exogenous_vars)
        self._model = model
        self._concs = None
        for attr, value in kwargs.items():
            setattr(self, attr, value)
    
    model = property(CSTR.suspended_growth_model.fget, CSTR.suspended_growth_model.fset)
    
    @property
    def state(self):
        '''The state of the BatchExperiment, i.e., component concentrations.'''
        if self._state is None: return None
        else:
            return dict(zip(self.components.IDs, self._state[:-1]))

    @state.setter
    def state(self, concs):
        concs = np.asarray(concs)
        if concs.shape != (len(self.components), ):
            raise ValueError(f'state must be a 1D array of length {len(self.components)},'
                              'indicating component concentrations.')
        self._state[:-1] = concs
    
    def set_init_conc(self, **kwargs):
        '''set the initial concentrations of the BatchExperiment.'''
        Cs = np.zeros(len(self.components))
        cmpx = self.components.index
        for k, v in kwargs.items(): Cs[cmpx(k)] = v
        self._concs = Cs
    
    def _init_state(self):
        if self._concs is None: 
            raise RuntimeError('must `set_init_conc` before starting simulation.')
        self._state = np.append(self._concs, 0)
        self._dstate = self._state * 0.

    def _update_state(self):
        pass

    def _update_dstate(self):
        pass

    def _run(self):
        pass

    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE

    #!!! how to considered sealed vs. open batch reactor (e.g., gaseous products/reactants)
    def _compile_ODE(self):
        if self._model is None:
            warn(f'{self.ID} was initialized without a kinetic model, '
                 f'and thus run as a non-reactive unit')
            rs = np.zeros(len(self.components)+1) 
            r = lambda state_arr: rs
        else:
            r = self.model.production_rates_eval

        _dstate = self._dstate
        hasexo = bool(len(self._exovars))
        f_exovars = self.eval_exo_dynamic_vars
        
        def dy_dt(t, QC_ins, QC, dQC_ins):
            if hasexo: QC = np.append(QC, f_exovars(t))
            _dstate[:-1] = r(QC)
            
        self._ODE = dy_dt
        
    #TODO: add functions for convenient model calibration
    
#%% NOT READY
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
            raise RuntimeError(f'{self.ID} was initialized without a suspended growth model.')
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

#%%
class PFR(SanUnit):

    _N_ins = 1
    _N_outs = 1
    _ins_size_is_fixed = False
    _outs_size_is_fixed = True

    def __init__(self, ID='', ins=None, outs=(), thermo=None, init_with='WasteStream',
                 N_tanks_in_series=5, V_tanks=[1500, 1500, 3000, 3000, 3000], 
                 influent_fractions=[[1.0, 0,0,0,0]], internal_recycles=[(4,0,35000),],
                 DO_setpoints=[], kLa=[0, 0, 120, 120, 60], DO_sat=8.0,
                 DO_ID='S_O2', suspended_growth_model=None, 
                 isdynamic=True, exogenous_vars=(), **kwargs):
        if exogenous_vars:
            warn(f'currently exogenous dynamic variables are not supported in process simulation for {self.ID}')
        SanUnit.__init__(self, ID, ins, outs, thermo, init_with, isdynamic=isdynamic,
                         exogenous_vars=exogenous_vars, **kwargs)
        self.N_tanks_in_series = N_tanks_in_series
        self.V_tanks = V_tanks
        self.influent_fractions = influent_fractions
        self.internal_recycles = internal_recycles
        self.DO_setpoints = DO_setpoints
        self.kLa = kLa
        self.DO_sat = DO_sat
        self.DO_ID = DO_ID
        self.suspended_growth_model = suspended_growth_model
        self._concs = None
        self._Qs = self.V_tanks * 0

    @property
    def V_tanks(self):
        '''[iterable[float]] Volumes of CSTRs in series [m3]'''
        return self._Vs
    @V_tanks.setter
    def V_tanks(self, Vs):
        if not iter(Vs): 
            raise TypeError(f'V_tanks must be an iterable, not {type(Vs).__name__}')
        elif len(Vs) != self.N_tanks_in_series:
            raise RuntimeError(f'cannot set the volumes of {self.N_tanks_in_series} tanks'
                               f'in series with {len(Vs)} value(s).')
        else: self._Vs = np.asarray(Vs)

    @property
    def influent_fractions(self):
        '''[iterable[float]] Fractions of influents fed into different zones [unitless]'''
        return self._f_ins
    @influent_fractions.setter
    def influent_fractions(self, fs):
        fs = np.asarray(fs)
        if len(fs.shape) == 1:
            if len(self.ins) != 1:
                raise RuntimeError(f'influent fractions must be a 2d array, with {len(self.ins)} row(s)'
                                   f'and {self.N_tanks_in_series} columns')
            if fs.shape[0] != self.N_tanks_in_series:
                raise RuntimeError('cannot set the fractions of influent fed into '
                                   f'{self.N_tanks_in_series} tanks in series with {len(fs)} value(s).')
            fs = fs.reshape(1, len(fs))
        elif len(fs.shape) > 2: 
            raise RuntimeError(f'influent fractions must be a 2d array, with {len(self.ins)} row(s)'
                               f'and {self.N_tanks_in_series} columns')
        if (fs < 0).any():
            raise ValueError('influent fractions must not have negative value(s)')
        for row in fs:
            rowsum = row.sum()
            if rowsum != 1: row /= rowsum
        self._f_ins = fs
    
    @property
    def internal_recycles(self):
        '''[list[3-tuple[int, int, float]]] A list of three-element tuples (i, j, Q) indicating internal recycling
         streams from zone i to zone j with a flowrate of Q [m3/d], respectively. 
         Indices i,j start from 0 not 1.'''
        return self._rcy
    @internal_recycles.setter
    def internal_recycles(self, rcy):
        isa = isinstance
        if isa(rcy, tuple):
            if len(rcy) != 3: 
                raise TypeError('internal recycles must be indicated by a list of 3-tuples')
            rcy = [rcy]
        elif isa(rcy, list):
            for row in rcy:
                if len(row) != 3:
                    raise TypeError('internal recycles must be indicated by a list of 3-tuples')
        else:
            raise TypeError('internal recycles must be indicated by a list of 3-tuples')
        _rcy = []
        for row in rcy:
            i, j, q = row
            _rcy.append((int(i), int(j), q))
        self._rcy = _rcy

    @property
    def DO_setpoints(self):
        '''[iterable[float]] Dissolve oxygen setpoints of CSTRs in series [mg-O2/L]. 
        0 is treated as no active aeration. DO setpoints take priority over kLa values.'''
        return self._DOs
    @DO_setpoints.setter
    def DO_setpoints(self, DOs):
        if not iter(DOs): 
            raise TypeError(f'DO setpoints must be an iterable, not {type(DOs).__name__}')
        elif 0 < len(DOs) != self.N_tanks_in_series:
            raise RuntimeError(f'cannot set the DO setpoints of {self.N_tanks_in_series} tanks'
                               f'in series with {len(DOs)} value(s).')
        else: self._DOs = np.asarray(DOs)

    @property
    def kLa(self):
        '''[iterable[float]] Aeration kLa values of CSTRs in series [d^(-1)]. If DO
        setpoints are specified, kLa values would be ignored in process simulation.'''
        return self._kLa
    @kLa.setter
    def kLa(self, ks):
        if not iter(ks): 
            raise TypeError(f'V_tanks must be an iterable, not {type(ks).__name__}')
        elif len(ks) != self.N_tanks_in_series:
            raise RuntimeError(f'cannot set kLa of {self.N_tanks_in_series} tanks'
                               f'in series with {len(ks)} value(s).')
        else: self._kLa = np.asarray(ks)

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

                
    def _run(self):
        out, = self.outs
        out.mix_from(self.ins)
    
    @property
    def state(self):
        '''The state of the PFR, including component concentrations [mg/L] and flow rate [m^3/d] for each zone.'''
        if self._state is None: return None
        else:
            N = self.N_tanks_in_series
            y = self._state.copy()
            y = y.reshape((N, int(len(y)/N)))
            y = pd.DataFrame(y, columns=self._state_header)
            return y

    def set_init_conc(self, concentrations=None, i_zone=None, **kwargs):
        '''set the initial concentrations [mg/L] of specific zones.'''
        isa = isinstance
        cmps = self.components
        N = self.N_tanks_in_series
        if self._concs is None: self._concs = np.zeros((N, len(cmps)))
        if concentrations is None:
            concs = cmps.kwarray(kwargs)
            if i_zone is None:
                self._concs[:] = concs
            else:
                self._concs[i_zone] = concs
        elif isa(concentrations, pd.DataFrame):
            concentrations.index = range(N)
            dct = concentrations.to_dict('index')
            for i, concs in dct:
                self._concs[i] = cmps.kwarray(concs)
        elif isa(concentrations, np.ndarray):
            if concentrations.shape != self._concs.shape:
                raise RuntimeError(f'cannot set the concentrations of {len(cmps)} '
                                   f'components across {N} with a {concentrations.shape} array')
            self._concs = concentrations
        else:
            raise TypeError('specify initial concentrations with pandas.DataFrame, numpy.ndarray'
                            'or "<component ID> = value" kwargs')

    def _init_state(self):
        self._state_header = ny = list(self.components.IDs) + ['Q']
        self._Qs_idx = list(len(ny) * np.arange(1, 1+self.N_tanks_in_series) - 1)
        out, = self.outs
        y = np.zeros((self.N_tanks_in_series, len(ny)))
        if self._concs is not None:
            y[:,:-1] = self._concs
        else: 
            y[:,:-1] = out.conc
        y[:,-1] = out.F_vol*24
        self._state = y.flatten()
        self._dstate = self._state * 0.

    def _update_state(self):
        out, = self.outs
        n_cmp = len(self.components)
        self._state[self._Qs_idx] = self._Qs
        out.state[:-1] = self._state[-(n_cmp+1):-1]
        out.state[-1] = sum(ws.state[-1] for ws in self.ins)

    def _update_dstate(self):
        out, = self.outs
        n_cmp = len(self.components)
        out.dstate[:-1] = self._dstate[-(n_cmp+1):-1]
        out.dstate[-1] = sum(ws.dstate[-1] for ws in self.ins)

    @property
    def ODE(self):
        if self._ODE is None:
            self._compile_ODE()
        return self._ODE

    def _compile_ODE(self):
        _Qs = self._Qs
        _dstate = self._dstate
        _update_dstate = self._update_dstate
        N = self.N_tanks_in_series
        _1_ov_V = np.diag(1/self.V_tanks)
        f_in = self.influent_fractions
        DO = self.DO_setpoints
        kLa = self.kLa
        rcy = self.internal_recycles
        DO_idx = self.components.index(self.DO_ID)
        DOsat = self.DO_sat
        Q_internal = np.zeros(N)
        for i_from, i_to, qr in rcy:
            Q_internal[i_to: i_from] += qr

        if self._model is None:
            warn(f'{self.ID} was initialized without a suspended growth model, '
                 f'and thus run as a non-reactive unit')
            Rs = lambda Cs: 0.
        else:
            f_rho = self._model.rate_function
            M_stoi = self._model.stoichio_eval()
            # @njit
            def Rs(Cs):
                rhos = np.vstack([f_rho(c) for c in Cs])
                rxn = rhos @ M_stoi
                return rxn
        
        if any(DO):
            aerated_zones = (DO > 0)
            aerated_DO = DO[aerated_zones]
            # @njit
            def dy_dt(t, QC_ins, QC, dQC_ins):
                y = QC.reshape((N, int(len(QC)/N)))
                Cs = y[:,:-1]
                Cs[aerated_zones, DO_idx] = aerated_DO
                _Qs[:] = Qs = np.dot(QC_ins[:,-1], f_in).cumsum() + Q_internal
                M_ins = f_in.T @ np.diag(QC_ins[:,-1]) @ QC_ins[:,:-1]  # N * n_cmps
                M_outs = np.diag(Qs) @ Cs
                M_ins[1:] += M_outs[:-1]
                for i_from, i_to, qr in rcy:
                    M_rcy = Cs[i_from] * qr
                    M_ins[i_to] += M_rcy
                    if i_from + 1 < N: M_ins[i_from+1] -= M_rcy
                rxn = Rs(Cs)
                dy = np.zeros_like(y)
                dy[:,:-1] = _1_ov_V @ (M_ins - M_outs) + rxn
                dy[aerated_zones, DO_idx] = 0.
                _dstate[:] = dy.flatten()
                _update_dstate()

        else:
            if not any(kLa): kLa = np.zeros(N)
            # @njit
            def dy_dt(t, QC_ins, QC, dQC_ins):
                # breakpoint()
                y = QC.reshape((N, int(len(QC)/N)))
                Cs = y[:,:-1]    
                do = Cs[:, DO_idx]
                aer = kLa*(DOsat-do)
                # aer[aer < 0] = 0.
                _Qs[:] = Qs = np.dot(QC_ins[:,-1], f_in).cumsum() + Q_internal
                M_ins = f_in.T @ np.diag(QC_ins[:,-1]) @ QC_ins[:,:-1]  # N * n_cmps
                M_outs = np.diag(Qs) @ Cs
                M_ins[1:] += M_outs[:-1]
                for i_from, i_to, qr in rcy:
                    M_rcy = Cs[i_from] * qr
                    M_ins[i_to] += M_rcy
                    if i_from + 1 < N: M_ins[i_from+1] -= M_rcy
                rxn = Rs(Cs)
                dy = np.zeros_like(y)
                dy[:,:-1] = _1_ov_V @ (M_ins - M_outs) + rxn
                dy[:,DO_idx] += aer 
                _dstate[:] = dy.flatten()
                _update_dstate()
                
        self._ODE = dy_dt

    def _design(self):
        pass