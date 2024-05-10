# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

Part of this module is based on the Thermosteam package:
https://github.com/BioSTEAMDevelopmentGroup/thermosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from warnings import warn
from . import Component, Components
from .utils import load_data, get_stoichiometric_coeff
from thermosteam.utils import chemicals_user, read_only
from thermosteam import settings
from sympy import symbols, Matrix, simplify, lambdify
from sympy.parsing.sympy_parser import parse_expr
import numpy as np
import pandas as pd

__all__ = ('DynamicParameter', 'Kinetics', 'MultiKinetics',
           'Process', 'Processes', 'CompiledProcesses', )

_load_components = settings.get_default_chemicals

class UndefinedProcess(AttributeError):
    '''AttributeError regarding undefined Component objects.'''
    def __init__(self, ID):
        super().__init__(repr(ID))

#%%
class DynamicParameter:
    """
    Create a :class:`DynamicParameter` object which defines a parameter
    in a :class:`Process`'s stoichiometry that dynamically depends on the state
    variables.

    Parameters
    ----------
    symbol : str
        The symbol of the dynamic parameter in the stoichiometric coefficient.
    function : Callable or float or int, optional
        The function that returns the value of the dynamic parameter when given
        the array of state variables and additional parameters as positional
        arguments. The default is None.
    params : dict, optional
        Additional parameters input to the function for the calculation of the
        dynamic parameters. The default is {}.

    Examples
    --------
    None.

    """
    def __init__(self, symbol, function=None, params={}):
        self.symbol = symbol
        self.function = function
        self._params = params

    @property
    def symbol(self):
        '''[sympy.Symbol] Symbol of the dynamic parameter in stoichiometric coefficients.'''
        return self._symbol
    @symbol.setter
    def symbol(self, s):
        self._symbol = symbols(s)

    @property
    def function(self):
        '''[Callable] Function that evaluates the dynamic parameter.'''
        return self._function

    @function.setter
    def function(self, f):
        if callable(f):
            nargs = f.__code__.co_argcount
            if nargs > 2:
                raise ValueError(f'function for the {self.__repr__()} must take '
                                 f'at most 2 positional arguments: an array '
                                 f'of state variables and a dictionary of parameters, '
                                 f'not {nargs} positional arguments.')
            elif nargs == 0:
                self._function = lambda state_arr, params: f()
            elif nargs == 1:
                self._function = lambda state_arr, params: f(state_arr)
            else:
                self._function = f
        elif isinstance(f, (float, int)):
            self._function = lambda state_arr, params: f
        elif f is None:
            self._function = lambda state_arr, params: None
        else:
            raise TypeError(f'function must be a callable or a float or int, if not None, '
                            f'not {type(f)}')

    @property
    def params(self):
        '''[dict] Additional parameters input to the function.'''
        return self._params

    def set_params(self, **params):
        '''
        Set values for the parameters needed for the evaluation of the dynamic parameter.

        Parameters
        ----------
        **params : float or int
            Parameters and values.
        '''
        self._params.update(params)

    def __call__(self, state_arr):
        '''Return the evaluated dynamic parameter when given an array of state variables.'''
        return self.function(state_arr, self.params)

    def __repr__(self):
        return f'<DynamicParameter: {str(self.symbol)}>'

#%%
class Kinetics(DynamicParameter):
    """
    Create a :class:`Kinetics` object which defines a :class:`Process` object's
    kinetic rate. It is a subclass of :class:`DynamicParameter`.

    Parameters
    ----------
    process : :class:`Process`
        The process that features this kinetics.
    function : Callable or float or int, optional
        The function that evaluates the kinetic rate of reaction when given
        the array of state variables and additional parameters as positional
        arguments. The default is None.
    params : dict, optional
        Additional parameters input to the function for the calculation of the
        kinetic rate of reaction. The default is {}.

    Examples
    --------
    None.

    """

    def __init__(self, process, function=None, params={}):
        super().__init__(symbol=f'rho_{process.ID}', function=function,
                         params=params)
        self.process = process

    # @property
    # def process(self):
    #     '''[:class:`Process`] The process whose reaction rate is defined by this Kinetics object.'''
    #     return self._process
    # @process.setter
    # def process(self, pc):
    #     if not isinstance(pc, Process):
    #         raise TypeError(f'processes must be of type `Process`, not {type(pc)}')
    #     self._process = pc

    def copy(self, new_process=None):
        '''Return a copy.'''
        pc = new_process or self.process
        copy = object.__new__(Kinetics, process=pc,
                              function=self.function,
                              params=self.params)
        for attr, v in self.__dict__.items():
            if attr not in ('_process', '_function', '_params'):
                setattr(copy, attr, v)
        return copy

    def __repr__(self):
        return f'<Kinetics: {self.process.ID}>'

class MultiKinetics:
    """
    Create a :class:`MultiKinetics` object which defines a :class:`Processes` object's
    kinetic rates.

    Parameters
    ----------
    processes : :class:`Processes`
        The process that features this kinetics.
    function : Callable or Iterable[float or int], optional
        The function that returns the array of kinetic rates of the processes when
        given the array of state variables and additional parameters as positional
        arguments. The default is None.
    params : dict, optional
        Additional parameters input to the function for the calculation of the
        kinetic rates. The default is {}.

    Examples
    --------
    None.

    """
    def __init__(self, processes, function=None, params={}):
        self.processes = processes
        self.function = function
        self._params = params

    # @property
    # def processes(self):
    #     '''[:class:`Processes`] The process whose reaction rate is defined
    #     by this Kinetics object.'''
    #     return self._processes
    # @processes.setter
    # def processes(self, pc):
    #     if not isinstance(pc, (Processes, CompiledProcesses)):
    #         raise TypeError(f'processes must be of type `Processes` or '
    #                         f'`CompiledProcesses`, not {type(pc)}')
    #     self._processes = pc

    @property
    def function(self):
        '''[Callable] Function that evaluates the kinetic rates.'''
        return self._function

    @function.setter
    def function(self, f):
        if callable(f):
            nargs = f.__code__.co_argcount
            if nargs > 2:
                raise ValueError(f'function for the {self.__repr__()} must take '
                                 f'at most 2 positional arguments: an array '
                                 f'of state variables and a dictionary of parameters, '
                                 f'not {nargs} positional arguments.')
            elif nargs == 0:
                self._function = lambda state_arr, params: f()
            elif nargs == 1:
                self._function = lambda state_arr, params: f(state_arr)
            else:
                self._function = f
        elif f is None:
            self._function = self._collect_kinetics()
        else:
            try: f = np.array(f, dtype=float)
            except TypeError:
                raise TypeError(f'function must be a callable or an array '
                                f'of length {len(self.processes)}, if not None, '
                                f'not {type(f)}')
            if f.shape == (len(self.processes),):
                self._function = lambda state_arr, params: f
            else:
                raise ValueError(f'function, if provided as an array, '
                                f'must have equal length as {self.processes}, '
                                f'not of shape {f.shape}')

    def _collect_kinetics(self):
        pcs = self.processes
        if any([pc.rate_function == None for pc in pcs]): return None
        rho_arr = np.empty(len(pcs))
        def f(state_arr, params):
            rho_arr[:] = [p.rate_function.function(state_arr,
                                                   {**p.rate_function.params, **params})
                          for p in pcs]
            return rho_arr
        return f

    @property
    def params(self):
        '''[dict] Additional parameters input to the function.'''
        return self._params

    def set_params(self, **params):
        '''
        Set values for the parameters needed for the evaluation of the kinetic rates.

        Parameters
        ----------
        **params : float or int
            Parameters and values.
        '''
        self._params.update(params)

    def __call__(self, state_arr):
        '''Return the evaluated array of kinetic rates when given
        an array of state variables.'''
        #!!! consider allow function to return multiple outputs including rho values,
        # updated parameter value, and/or intermediate variables that are useful but
        # don't need integration
        return self.function(state_arr, self.params)

    def __repr__(self):
        return f'<MultiKinetics({[p.ID for p in self.processes]})>'

#%%
@chemicals_user
class Process:
    """
    Create a :class:`Process` object which defines a stoichiometric process and its kinetics.
    A :class:`Process` object is capable of reacting the component flow rates of a :class:`WasteStream`
    object.

    Parameters
    ----------
    ID : str
        A unique identification.
    reaction : dict, str, or numpy.ndarray
        A dictionary of stoichiometric coefficients with component IDs as
        keys, or a numeric array of stoichiometric coefficients, or a string
        of a stoichiometric equation written as:
        i1 R1 + ... + in Rn -> j1 P1 + ... + jm Pm.
        Stoichiometric coefficients can be symbolic or numerical.
        Unknown stoichiometric coefficients to solve for should be expressed as "?".
    ref_component : str
        ID of the reference :class:`Component` object of the process rate.
    rate_equation : str, optional
        The kinetic rate equation of the process. The default is None.
    components=None : class:`CompiledComponents`, optional
        Components corresponding to each entry in the stoichiometry array,
        defaults to all components set in the system (i.e., through :func:`set_thermo`).
    conserved_for : tuple[str], optional
        Materials subject to conservation rules, must be an 'i\_' attribute of
        the components. The default is ("COD", "N", "P", "charge").
    parameters : Iterable[str], optional
        Symbolic parameters in stoichiometry coefficients (both constant and dynamic)
        and/or rate equation. The default is None.

    Examples
    --------
    To create a :class:`Process` object, basic information including stoichiometry,
    and reference component must be specified. Unknown stoichiometric coefficients
    can be solved automatically based on conservation of materials.

    >>> import qsdsan as qs
    >>> cmps = qs.processes.create_asm1_cmps()
    >>> pc1 = qs.Process(ID='aerobic_hetero_growth',
    ...                  reaction='[1/Y_H]S_S + [1-1/Y_H]S_O + [?]S_NH + [?]S_ALK -> X_BH',
    ...                  ref_component='X_BH',
    ...                  rate_equation='mu_H*S_S/(K_S+S_S)*S_O/(K_O_H+S_O)*S_NH/(K_NH+S_NH)*X_BH',
    ...                  conserved_for=('COD', 'N', 'charge'),
    ...                  parameters=('Y_H', 'mu_H', 'K_S', 'K_O_H', 'K_NH'))
    >>> pc1.show()
    Process: aerobic_hetero_growth
    [stoichiometry]      S_S: -1/Y_H
                         X_BH: 1.00
                         S_O: -1.0 + 1/Y_H
                         S_NH: -0.0860
                         S_ALK: -0.0737
    [reference]          X_BH
    [rate equation]      S_NH*S_O*S_S*X_BH*mu_H/((K_N...
    [parameters]         Y_H: Y_H
                         mu_H: mu_H
                         K_S: K_S
                         K_O_H: K_O_H
                         K_NH: K_NH
    [dynamic parameters]

    If all stoichiometric coefficients and relevant parameters are specified, one can also
    check if materials are conserved. No returns when all materials are conserved.

    >>> pc2 = qs.Process(ID='hydrolysis',
    ...                  reaction={'S_S':1, 'X_S':-1},
    ...                  ref_component='S_S',
    ...                  rate_equation='k_h*X_S/(K_X*X_BH+X_S)*(S_O/(K_O_H+S_O)+eta_h*K_O_H/(K_O_H+S_O)*S_NO/(K_NO+S_NO))*X_BH',
    ...                  conserved_for=(),
    ...                  parameters=('k_h', 'K_X', 'K_O_H', 'eta_h', 'K_NO'))
    >>> pc2.show()
    Process: hydrolysis
    [stoichiometry]      S_S: 1
                         X_S: -1
    [reference]          S_S
    [rate equation]      X_BH*X_S*k_h*(K_O_H*S_NO*eta...
    [parameters]         k_h: k_h
                         K_X: K_X
                         K_O_H: K_O_H
                         eta_h: eta_h
                         K_NO: K_NO
    [dynamic parameters]

    >>> pc2.conserved_for = ('COD', 'N')
    >>> pc2.check_conservation()

    Raise warning when materials are not strictly conserved based on the given
    stoichiometric coefficients with defined parameters. An error will be raised
    instead for stoichiometric coefficients that are purely numerical.

    >>> pc3 = qs.Process(ID='decay_hetero',
    ...                  reaction='X_BH -> [f_P]X_P + [1-f_P]X_S + [0.05-0.05*f_P]X_ND',
    ...                  ref_component='X_BH',
    ...                  rate_equation='b_H*X_BH',
    ...                  conserved_for=('COD', 'N'),
    ...                  parameters=('f_P', 'b_H'))
    >>> pc3.check_conservation() # doctest: +SKIP
    UserWarning: The following materials aren't strictly conserved by the stoichiometric coefficients.
    A positive value means the material is created, a negative value means the material is destroyed:
      N: 0.01*f_P - 0.036

    Parameter values can be set. Changes will be reflected in the stoichiometry
    and/or the rate equation:

    >>> pc2.set_parameters(eta_h = 0.8)
    >>> pc2.parameters
    {'k_h': k_h, 'K_X': K_X, 'K_O_H': K_O_H, 'eta_h': 0.8, 'K_NO': K_NO}

    >>> str(pc2.rate_equation)
    'X_BH*X_S*k_h*(0.8*K_O_H*S_NO/((K_NO + S_NO)*(K_O_H + S_O)) + S_O/(K_O_H + S_O))/(K_X*X_BH + X_S)'

    An reversed process can be easily created:

    >>> pc2_reversed = pc2.copy('hydrolysis_reversed')
    >>> pc2_reversed.reverse()
    >>> pc2_reversed.show()
    Process: hydrolysis_reversed
    [stoichiometry]      S_S: -1
                         X_S: 1
    [reference]          S_S
    [rate equation]      -X_BH*X_S*k_h*(K_O_H*S_NO*et...
    [parameters]         k_h: k_h
                         K_X: K_X
                         K_O_H: K_O_H
                         eta_h: 0.8
                         K_NO: K_NO
    [dynamic parameters]

    See Also
    --------
    `numpy.isclose <https://numpy.org/doc/stable/reference/generated/numpy.isclose.html>`_ for ``rtol`` and ``atol`` settings.
    """
    def __init__(self, ID, reaction, ref_component, rate_equation=None,
                 components=None, conserved_for=('COD', 'N', 'P', 'charge'),
                 parameters=()):
        self._ID = ID
        self._reaction = reaction
        self.components = self._load_chemicals(components)
        self._ref_component = ref_component
        self.conserved_for = conserved_for
        self._parameters = {p: symbols(p) for p in parameters}
        self._stoichiometry = get_stoichiometric_coeff(
            reaction, self._ref_component, self._components, self._conserved_for,
            self.parameters)
        self._parse_rate_eq(rate_equation)
        self._dyn_params={}
        self.rate_function=None

    def get_conversion_factors(self, as_matrix=False):
        '''
        Return conversion factors (i.e., the 'i\_' attributes of the components)
        as a numpy.ndarray or a SymPy Matrix.
        '''
        conserved_for = self._conserved_for
        if conserved_for:
            cmps = self._components
            getfield = getattr
            arr = np.array([getfield(cmps, f'i_{x}') for x in conserved_for])
            if as_matrix: return Matrix(arr.tolist())
            return arr
        else: return None

    def check_conservation(self, rtol=1e-5, atol=1e-8):
        '''
        Check conservation of materials (based on the rules set in the `conserved_for` attribute)
        given purely numerical stoichiometric coefficients or stoichiometric coefficients
        with defined parameters.

        No return indicates that all rules are satisfied.

        Parameters
        ----------
        rtol : float
            Relative tolerance. Only applicable to purely numerical coefficients.
        atol : float
            Absolute tolerance. Only applicable to purely numerical coefficients.

        See Also
        --------
        `numpy.isclose <https://numpy.org/doc/stable/reference/generated/numpy.isclose.html>`_ for ``rtol`` and ``atol`` settings.
        '''
        isa = isinstance
        v = self._stoichiometry
        ic = self.get_conversion_factors()
        if isa(v, np.ndarray):
            if ic is None: # no `conserved_for` attribute
                warn('No available `conserved_for` attributes, cannot check conservation.')
                return
            ic_dot_v = ic @ v
            ic_dot_v = np.array([ic_dot_v]) if ic_dot_v.size==1 else ic_dot_v
            conserved_arr = np.isclose(ic_dot_v, np.zeros(ic_dot_v.shape), rtol=rtol, atol=atol)
            conserved_arr = np.array([conserved_arr]) if conserved_arr.size==1 else conserved_arr
            if not conserved_arr.all():
                materials = self._conserved_for
                unconserved = [(materials[i], ic_dot_v[i]) for i, conserved in enumerate(conserved_arr) if not conserved]
                raise RuntimeError("The following materials are unconserved by the "
                                   "stoichiometric coefficients. A positive value "
                                   "means the material is created, a negative value "
                                   "means the material is destroyed:\n "
                                   + "\n ".join([f"{material}: {value:.2f}" for material, value in unconserved]))
        elif isa(v, list):
            if ic is None:
                warn('No available `conserved_for` attributes, cannot check conservation.')
                return
            ic_dot_v = list(simplify(ic * Matrix(v)))
            materials = self._conserved_for
            unconserved = [(materials[i], prod) for i, prod in enumerate(ic_dot_v) if prod != 0]
            if len(unconserved) > 0:
                warn("The following materials aren't strictly conserved by the "
                     "stoichiometric coefficients. A positive value "
                     "means the material is created, a negative value "
                     "means the material is destroyed:\n "
                     + "\n ".join([f"{material}: {value}" for material, value in unconserved]))
        else:
            raise RuntimeError("Can only check conservations with purely numerical "
                               "stoichiometric coefficients or coefficients with defined parameters.")

    def reverse(self):
        '''Reverse the process as to flip the signs of stoichiometric coefficients of all components.'''
        if isinstance(self._stoichiometry, np.ndarray):
            self._stoichiometry = -self._stoichiometry
        else:
            self._stoichiometry = [-v for v in self._stoichiometry]
        if self._rate_equation:
            self._rate_equation = -self._rate_equation
        self._rate_function = None

    @property
    def ID(self):
        '''[str] A unique identification.'''
        return self._ID
    @ID.setter
    def ID(self, ID):
        self._ID = ID

    @property
    def reaction(self):
        '''
        [dict, str, or numpy.ndarray]

        A dictionary of stoichiometric coefficients with component IDs as
        keys, or a numeric array of stoichiometric coefficients, or a string
        of a stoichiometric equation written as:
        i1 R1 + ... + in Rn -> j1 P1 + ... + jm Pm.
        Stoichiometric coefficients can be symbolic or numerical.
        Unknown stoichiometric coefficients to solve for should be expressed as "?".
        '''
        return self._reaction
    @reaction.setter
    def reaction(self, rxn):
        self._reaction = rxn
        if rxn:
            self._stoichiometry = get_stoichiometric_coeff(
                rxn, self._ref_component, self._components, self._conserved_for, self._parameters)

    @property
    def ref_component(self):
        '''
        [str] ID of the reference component.

        .. note::

            When a new value is assigned, all stoichiometric coefficient will be
            normalized so that the new stoichiometric coefficient of the new reference
            component is 1 or -1. The rate equation will also be updated automatically.
        '''
        return getattr(self._components, self._ref_component).ID
    @ref_component.setter
    def ref_component(self, ref_cmp):
        if ref_cmp:
            self._ref_component = str(ref_cmp) # in case a component obj is passed
            self._normalize_rate_eq(ref_cmp)
            self._normalize_stoichiometry(ref_cmp)

    @property
    def components(self):
        '''
        [:class:`CompiledComponents`]

        Components corresponding to each entry in the stoichiometry array,
        defaults to all components set in the system (i.e., through :func:`set_thermo`).
        '''
        return self._components
    @components.setter
    def components(self, cmps):
        self._components = _load_components(cmps)

    @property
    def conserved_for(self):
        '''
        [tuple] Materials subject to conservation rules, must have corresponding
        'i\_' attributes for the components.
        '''
        return self._conserved_for
    @conserved_for.setter
    def conserved_for(self, materials):
        get = getattr
        for mat in materials:
            try: get(Component, f'i_{mat}')
            except: raise ValueError(f'Components do not have i_{mat} attribute.')
        self._conserved_for = materials

    @property
    def parameters(self):
        '''[dict] Symbolic parameters in stoichiometric coefficients (both
        constant and dynamic) and/or rate equation.'''
        return self._parameters

    def append_parameters(self, *new_pars):
        '''Append new symbolic parameters'''
        for p in new_pars:
            self._parameters[p] = symbols(p)

    def set_parameters(self, **parameters):
        '''
        Set values for symbolic stoichiometric and/or kinetic parameters.

        Parameters
        ----------
        **parameters : float or int
            Parameters and values.
        '''
        self._parameters.update(parameters)

    def dynamic_parameter(self, function=None, symbol=None, params={}):
        '''Add a function for the evaluation of a dynamic parameter in the
        stoichiometry of the process. Creates a :class:`DynamicParameter` object.'''
        if symbol:
            if symbol not in self.parameters.keys():
                warn(f'new symbolic parameter {symbol} added.')
                self.append_parameters(symbol)
            if not function:
                return lambda f: self.dynamic_parameter(f, symbol, params=params)
            dp = DynamicParameter(symbol, function, params)
            self._dyn_params[symbol] = dp
            return dp

    @property
    def stoichiometry(self):
        '''[dict] Non-zero stoichiometric coefficients.'''
        allcmps = dict(zip(self._components.IDs, self._stoichiometry))
        active_cmps = {k:v for k,v in allcmps.items() if v != 0}
        isa = isinstance
        if isa(self._stoichiometry, list):
            active_cmps = {k:v.subs(self._parameters) \
                           if not isa(v, (float, int)) else v \
                               for k,v in active_cmps.items()}
        return active_cmps

    @property
    def rate_equation(self):
        '''
        [SymPy expression] Kinetic rate equation of the process. Also the rate in
        which the reference component is reacted or produced in the process. Kinetic
        parameters in the equation are replaced with their assigned values.
        '''
        if self._rate_equation:
            return self._rate_equation.subs(self._parameters)

    def kinetics(self, function=None, parameters={}):
        '''Add a function for the evaluation of the kinetic rate of the process.
        Creates a :class:`Kinetics` object.'''
        if not function:
            return lambda f: self.kinetics(f, parameters=parameters)
        else:
            self._rate_function = Kinetics(self, function, parameters)

    @property
    def rate_function(self):
        '''[:class:`Kinetics`] The function to evaluate the kinetic rate of the process.'''
        if self._rate_function is None and self._rate_equation:
            self._rate_eq2func()
        return self._rate_function
    @rate_function.setter
    def rate_function(self, k):
        if k is None:
            self._rate_function = None
        elif isinstance(k, Kinetics):
            if k.process is self:
                self._rate_function = k
            else:
                warn(f'attempted to use {k.__repr__()} for Process: {self.ID}. '
                     f'A copy was created instead.')
                self._rate_function = k.copy(self)
        elif callable(k):
            self.kinetics(function=k)
        else:
            raise TypeError(f'rate_function must be a `Kinetics` object, or '
                            f'a function that takes exactly 1 input argument '
                            f'(i.e., an array of state variables), or None, '
                            f'not {type(k)}')

    def _parse_rate_eq(self, eq):
        if eq is not None:
            state_symbols = {c: symbols(c) for c in self._components.IDs}
            params = self._parameters
            self._rate_equation = parse_expr(str(eq), {**state_symbols, **params})
        else: self._rate_equation = None

    def _rate_eq2func(self):
        var_kw = self._components.IDs
        var = list(symbols(var_kw)) + [symbols(p) for p in self._parameters.keys()]
        lamb = lambdify(var, self._rate_equation, 'numpy')
        def f(state_arr, params={}):
            states = dict(zip(var_kw, state_arr))
            return lamb(**states, **params)

        self.kinetics(function=f, parameters=self.parameters)

    def _normalize_stoichiometry(self, new_ref):
        isa = isinstance
        new_ref = str(new_ref)
        stoich = self._stoichiometry
        factor = abs(stoich[self._components.index(new_ref)])
        if isa(self._stoichiometry, np.ndarray):
            stoich /= factor
        elif isa(stoich, list):
            self._stoichiometry = [v/factor for v in stoich]

    def _normalize_rate_eq(self, new_ref):
        if self._rate_equation:
            factor = abs(self._stoichiometry[self._components.index(str(new_ref))])
            self._rate_equation *= factor

    def show(self):
        info = f"Process: {self.ID}"
        header = '\n[stoichiometry]'.ljust(22)
        section = []
        for cmp, stoichio in tuple(zip(self._components.IDs, self._stoichiometry)):
            if stoichio != 0:
                if isinstance(stoichio, (int, float)): line = f"{cmp}: {stoichio:.3g}"
                else: line = f"{cmp}: {stoichio.evalf(n=3)}"
                section.append(line)
        info += header + ("\n" + 21*" ").join(section)
        info += '\n[reference]'.ljust(22)+ f"{self.ref_component}"
        line = '\n[rate equation]'.ljust(22) + f"{self._rate_equation}"
        if len(line) > 50: line = line[:50] + '...'
        info += line
        header = '\n[parameters]'.ljust(22)
        section = []
        for k,v in self.parameters.items():
            if isinstance(v, (int, float)): line = f"{k}: {v:.3g}"
            else: line = f"{k}: {v}"
            section.append(line)
        info += header + ("\n" + 21*" ").join(section)
        header = '\n[dynamic parameters]'.ljust(22)
        section = []
        for v in self._dyn_params.values():
            line = v.__repr__()
            section.append(line)
        info += header + ("\n" + 21*" ").join(section)
        print(info)

    _ipython_display_ = show

    def copy(self, new_ID=''):
        '''
        Return a copy of the process with the same attributes (other than `ID`).

        Parameters
        ----------
        new_ID : str
            ID of the new process, will be set to the original process'ID with
            a "_copy" suffix if not provided.
        '''
        cls = self.__class__
        new = cls.__new__(cls)
        new_ID = new_ID or self.ID+'_copy'
        new.__init__(
            new_ID, reaction=self.reaction, ref_component=self.ref_component,
            # rate_equation=self._rate_equation, rate_function=self._rate_function,
            rate_equation=self._rate_equation,
            components=self._components, conserved_for=self.conserved_for,
            # parameters=self.parameters.keys(), dynamic_params=self._dyn_params.keys()
            parameters=self.parameters.keys()
            )
        new._parameters.update(self.parameters)
        new._dyn_params.update(self._dyn_params)
        rho = self._rate_function
        if rho is not None:
            new.rate_function(rho.function, rho.params)
        return new
    __copy__ = copy


#%%
setattr = object.__setattr__

@chemicals_user
class Processes:
    """
    Create a :class:`Processes` object that contains :class:`Process` objects as attributes.

    Parameters
    ----------
    processes : Iterable[:class:`Process`]

    Examples
    --------
    >>> import qsdsan as qs
    >>> cmps = qs.processes.create_asm1_cmps()
    >>> pc1 = qs.Process(ID='aerobic_hetero_growth',
    ...                  reaction='[1/Y_H]S_S + [1-1/Y_H]S_O + [?]S_NH + [?]S_ALK -> X_BH',
    ...                  ref_component='X_BH',
    ...                  rate_equation='mu_H*S_S/(K_S+S_S)*S_O/(K_O_H+S_O)*S_NH/(K_NH+S_NH)*X_BH',
    ...                  conserved_for=('COD', 'N', 'charge'),
    ...                  parameters=('Y_H', 'mu_H', 'K_S', 'K_O_H', 'K_NH'))
    >>> pc2 = qs.Process(ID='hydrolysis',
    ...                  reaction={'S_S':1, 'X_S':-1},
    ...                  ref_component='S_S',
    ...                  rate_equation='k_h*X_S/(K_X*X_BH+X_S)*(S_O/(K_O_H+S_O)+eta_h*K_O_H/(K_O_H+S_O)*S_NO/(K_NO+S_NO))*X_BH',
    ...                  conserved_for=(),
    ...                  parameters=('k_h', 'K_X', 'K_O_H', 'eta_h', 'K_NO'))
    >>> pcs = qs.Processes([pc1, pc2])
    >>> pcs.show()
    Processes([aerobic_hetero_growth, hydrolysis])

    >>> pcs.hydrolysis.show()
    Process: hydrolysis
    [stoichiometry]      S_S: 1
                         X_S: -1
    [reference]          S_S
    [rate equation]      X_BH*X_S*k_h*(K_O_H*S_NO*eta...
    [parameters]         k_h: k_h
                         K_X: K_X
                         K_O_H: K_O_H
                         eta_h: eta_h
                         K_NO: K_NO
    [dynamic parameters]

    >>> pc3 = qs.Process(ID='decay_hetero',
    ...                  reaction='X_BH -> [f_P]X_P + [1-f_P]X_S + [?]X_ND',
    ...                  ref_component='X_BH',
    ...                  rate_equation='b_H*X_BH',
    ...                  conserved_for=('COD', 'N'),
    ...                  parameters=('f_P', 'b_H'))
    >>> pcs.append(pc3)
    >>> pcs.compile()
    >>> pcs.show()
    CompiledProcesses([aerobic_hetero_growth, hydrolysis, decay_hetero])

    Once the processes are compiled, corresponding attributes become accessible:

    >>> pcs.parameters
    {'Y_H': Y_H,
      'mu_H': mu_H,
      'K_S': K_S,
      'K_O_H': K_O_H,
      'K_NH': K_NH,
      'k_h': k_h,
      'K_X': K_X,
      'eta_h': eta_h,
      'K_NO': K_NO,
      'f_P': f_P,
      'b_H': b_H}

    >>> pcs.production_rates.rate_of_production.loc['X_S']
    -1.0*X_BH*X_S*k_h*(K_O_H*S_NO*eta_h/((K_NO + S_NO)*(K_O_H + S_O)) + S_O/(K_O_H + S_O))/(K_X*X_BH + X_S) + X_BH*b_H*(1 - f_P)
    """

    def __new__(cls, processes):

        self = super().__new__(cls)
        setfield = setattr
        hasfield = hasattr
        duplicates = []
        for i in processes:
            if hasfield(self, i.ID): duplicates.append(i.ID)
            else: setfield(self, i.ID, i)
        if set(duplicates):
            raise ValueError(f'Processes with duplicate IDs were found: {set(duplicates)}')
        return self

    # def __getnewargs__(self):
    #     return(tuple(self),)

    def __setattr__(self, ID, process):
        raise TypeError("can't set attribute; use <Processes>.append or <Processes>.extend instead")

    def __setitem__(self, ID, process):
        raise TypeError("can't set attribute; use <Processes>.append or <Processes>.extend instead")

    def __getitem__(self, key):
        """
        Return a :class:`Process` object or a list of :class:`Process` objects.

        Parameters
        ----------
        key : Iterable[str] or str
              Process IDs.

        """
        dct = self.__dict__
        try:
            if isinstance(key, str):
                return dct[key]
            else:
                return [dct[i] for i in key]
        except KeyError:
            raise KeyError(f"undefined process {key}")

    def copy(self):
        """Return a copy."""
        copy = object.__new__(Processes)
        for proc in self: setattr(copy, proc.ID, proc)
        return copy

    def append(self, process):
        """Append a :class:`Process`."""
        if not isinstance(process, Process):
            raise TypeError("only 'Process' objects can be appended, "
                           f"not '{type(process).__name__}'")
        ID = process.ID
        if ID in self.__dict__:
            raise ValueError(f"{ID} already defined in processes")
        setattr(self, ID, process)

    def extend(self, processes):
        """Extend with more :class:`Process` objects."""
        if isinstance(processes, Processes):
            self.__dict__.update(processes.__dict__)
        else:
            for process in processes: self.append(process)

    def insert(self, index, process):
        """Insert a :class:`Process` object at a given index."""
        if not isinstance(process, Process):
            raise TypeError("only 'Process' objects can be inserted, "
                           f"not '{type(process).__name__}'")
        ID = process.ID
        if ID in self.__dict__:
            raise ValueError(f"{ID} already defined in processes")
        processes = list(self.__dict__.values())
        processes.insert(index, process)
        setattr(self, '__dict__', {p.ID: p for p in processes})

    def subgroup(self, IDs):
        """
        Create a new subgroup of processes.

        Parameters
        ----------
        IDs : Iterable[str]
              Process IDs.

        """
        return Process([getattr(self, i) for i in IDs])

    def compile(self, skip_checks=False, to_class=None):
        '''Cast as a :class:`CompiledProcesses` object unless otherwise specified.'''
        processes = tuple(self)
        setattr(self, '__class__', CompiledProcesses)
        try:
            self._compile(processes, skip_checks)
            if to_class is not None:
                setattr(self, '__class__', to_class)
        except Exception as error:
            setattr(self, '__class__', Processes)
            setattr(self, '__dict__', {i.ID: i for i in processes})
            raise error

    # kwarray = array = index = indices = must_compile

    def show(self):
        print(self)

    _ipython_display_ = show

    def __len__(self):
        return len(self.__dict__)

    def __contains__(self, process):
        if isinstance(process, str):
            return process in self.__dict__
        elif isinstance(process, Process):
            return process in self.__dict__.values()
        else: # pragma: no cover
            return False

    def __iter__(self):
        yield from self.__dict__.values()

    def __repr__(self):
        return f"{type(self).__name__}([{', '.join(self.__dict__)}])"

    _default_data = None

    @classmethod
    def load_from_file(cls, path='', components=None,
                       conserved_for=('COD', 'N', 'P', 'charge'), parameters=(),
                       use_default_data=False, store_data=False, compile=True,
                       **compile_kwargs):
        """
        Create :class:`CompiledProcesses` object from a table of process IDs, stoichiometric
        coefficients, and rate equations stored in a .tsv, .csv, or Excel file.

        Parameters
        ----------
        path : str, optional
            File path.
        components : :class:`CompiledComponents`, optional
            Components corresponding to the columns in the stoichiometry matrix,
            to all components set in the system (i.e., through :func:`set_thermo`).
        conserved_for : tuple[str], optional
            Materials subject to conservation rules, must have corresponding 'i\_'
            attributes for the components. Applied to all processes.
            The default is ('COD', 'N', 'P', 'charge').
        parameters : Iterable[str], optional
            Symbolic parameters in stoichiometry coefficients and/or rate equation.
        use_default_data : bool, optional
            Whether to use default data. The default is False.
        store_data : bool, optional
            Whether to store the file as default data. The default is False.
        compile : bool, optional
            Whether to compile processes. The default is True.

        .. note::

            [1] First column of the table should be process IDs, followed by stoichiometric
            coefficient matrix with corresponding component IDs as column names, and rate
            equations as the last column.

            [2] Entries of stoichiometric coefficients can be symbolic or numerical.
            Blank cells are considered zero.

            [3] Unknown stoichiometric coefficients to solve for using conservation
            rules should be uniformly written as '?'.

            [4] For each process, the first component with stoichiometric coefficient
            of -1 or 1 is considered the reference component. If none of the components
            has -1 or 1 stoichiometric coefficient, the first component with non-zero
            coefficient is considered the reference.
        """
        if use_default_data and cls._default_data is not None:
            data = cls._default_data
        else:
            data = load_data(path=path, index_col=None, na_values=0)

        cmps = _load_components(components)

        cmp_IDs = [i for i in data.columns if i in cmps.IDs]
        data.dropna(how='all', subset=cmp_IDs, inplace=True)
        new = cls(())
        for i, proc in data.iterrows():
            ID = proc[0]
            stoichio = proc[cmp_IDs]
            if data.columns[-1] in cmp_IDs: rate_eq = None
            else:
                if pd.isna(proc[-1]): rate_eq = None
                else: rate_eq = proc[-1]
            stoichio = stoichio[-pd.isna(stoichio)].to_dict()
            ref = None
            for k,v in stoichio.items():
                try:
                    v = stoichio[k] = float(v)
                    if ref is None and v in (-1, 1): ref = k
                except: continue
            if ref is None: ref = stoichio.keys()[0]
            process = Process(ID, stoichio,
                              ref_component=ref,
                              rate_equation=rate_eq,
                              components=cmps,
                              conserved_for=conserved_for,
                              parameters=parameters)
            new.append(process)

        if store_data:
            cls._default_data = data

        if compile: new.compile(**compile_kwargs)
        return new


#%%
@read_only(methods=('append', 'extend'))
class CompiledProcesses(Processes):
    """
    Create a :class:`CompiledProcesses` object that contains :class:`Process` objects as attributes.

    Parameters
    ----------
    processes : Iterable[:class:`Process`]

    """

    _cache = {}

    def __new__(cls, processes):
        cache = cls._cache
        processes = tuple(processes)
        if processes in cache:
            self = cache[processes]
        else:
            self = object.__new__(cls)
            setfield = setattr
            for i in processes:
                setfield(self, i.ID, i)
            self._compile()
            cache[processes] = self
        return self

    # def __dir__(self):
    #     pass

    def compile(self, skip_checks=False):
        """Do nothing, :class:`CompiledProcesses` objects are already compiled."""
        pass

    def _compile(self, processes, skip_checks=False):
        dct = self.__dict__
        isa = isinstance
        tuple_ = tuple
        # processes = tuple_(dct.values())
        # missing_rate_proc = []
        # for process in processes:
        #     if process.rate_equation is None:
        #         missing_rate_proc.append(process.ID)
        # if not skip_checks and len(missing_rate_proc) > 0:
        #     raise RuntimeError(f"The following processes are missing rate equations:"
        #                        f"{missing_rate_proc}")
        IDs = tuple_([i.ID for i in processes])
        size = len(IDs)
        index = tuple_(range(size))
        dct['tuple'] = processes
        dct['size'] = size
        dct['IDs'] = IDs
        dct['_index'] = index = dict(zip(IDs, index))
        if len(set([i._components for i in processes])) > 1:
            cmps = Components([cmp for i in processes for cmp in i._components])
        else:
            cmps = processes[0]._components
        dct['_components'] = _load_components(cmps)
        M_stch = []
        params = {}
        dyn_params = {}
        rate_eqs = tuple_([i._rate_equation for i in processes])
        all_numeric = True
        for i in processes:
            stch = [0]*cmps.size
            params.update(i._parameters)
            dyn_params.update(i._dyn_params)
            if all_numeric and isa(i._stoichiometry, (list, tuple)): all_numeric = False
            for cmp, coeff in i.stoichiometry.items():
                stch[cmps._index[cmp]] = coeff
            M_stch.append(stch)
        dct['_parameters'] = params
        dct['_dyn_params'] = dyn_params
        for i in processes:
            i._parameters = dct['_parameters']
            i._dyn_params = dct['_dyn_params']
        if all_numeric: M_stch = np.asarray(M_stch)
        dct['_stoichiometry'] = M_stch
        dct['_stoichio_lambdified'] = None
        dct['_rate_equations'] = rate_eqs
        if all(rate_eqs):
            dct['_production_rates'] = list(Matrix(M_stch).T * Matrix(rate_eqs))
        else:
            dct['_production_rates'] = None
        dct['_rate_function'] = None

    @property
    def parameters(self):
        '''[dict] All symbolic stoichiometric and kinetic parameters.'''
        return self._parameters

    def append_parameters(self, *new_pars):
        '''Append new symbolic parameters'''
        for p in new_pars:
            self._parameters[p] = symbols(p)

    def set_parameters(self, **parameters):
        '''Set values to stoichiometric and/or kinetic parameters.'''
        self._parameters.update(parameters)
        if self._stoichio_lambdified is not None:
            self.__dict__['_stoichio_lambdified'] = None

    def dynamic_parameter(self, function=None, symbol=None, params={}):
        '''Add a function for the evaluation of a dynamic parameter in the
        stoichiometry of the process. Creates a :class:`DynamicParameter` object.'''
        if symbol:
            if symbol not in self.parameters.keys():
                warn(f'new symbolic parameter {symbol} added.')
                self.append_parameters(symbol)
            if not function:
                return lambda f: self.dynamic_parameter(f, symbol, params=params)
            dp = DynamicParameter(symbol, function, params)
            self.__dict__['_dyn_params'][symbol] = dp
            return dp

    def params_eval(self, state_arr):
        '''Evaluate the dynamic parameters in the stoichiometry given an array of state variables.'''
        dct = self._parameters
        for k, p in self._dyn_params.items():
            dct[k] = p(state_arr)

    @property
    def stoichiometry(self):
        '''[pandas.DataFrame] Stoichiometric coefficients.'''
        stoichio = self._stoichiometry
        isa = isinstance
        v_params = {k:v for k,v in self._parameters.items() if isa(v, (float, int))}
        if isa(stoichio, list) and len(v_params) > 0:
            stoichio_vals = []
            for row in stoichio:
                # stoichio_vals.append([v.subs(v_params) if not isa(v, (float, int)) else v for v in row])
                stoichio_vals.append([v.evalf(subs=v_params) if not isa(v, (float, int)) else v for v in row])
            try: 
                stoichio_vals = np.asarray(stoichio_vals, dtype=float)
                #!!! round-off error
                # stoichio_vals[np.abs(stoichio_vals) < 2.22044604925e-16] = 0.0
            except TypeError: pass
            return pd.DataFrame(stoichio_vals, index=self.IDs, columns=self._components.IDs)
        else: return pd.DataFrame(stoichio, index=self.IDs, columns=self._components.IDs)

    def _lambdify_stoichio(self):
        dct = self._dyn_params
        dct_vals = self._parameters
        if dct:
            static_params = {k:v for k,v in dct_vals.items() if k not in dct}
            stoichio = []
            isa = isinstance
            for row in self._stoichiometry:
                stoichio.append([v.evalf(subs=static_params) if not isa(v, (float, int)) else v for v in row])
            sbs = [symbols(i) for i in dct.keys()]
            lamb = lambdify(sbs, stoichio, 'numpy')
            arr = np.empty((self.size, len(self._components)))
            def f():
                v = [dct_vals[k] for k in dct.keys()]
                arr[:,:] = lamb(*v)
                return arr
            self.__dict__['_stoichio_lambdified'] = f
        else:
            try:
                stoichio_arr = self.stoichiometry.to_numpy(dtype=float)
            except TypeError:
                isa = isinstance
                undefined = [k for k, v in dct_vals if not isa(v, (float, int))]
                raise TypeError(f'Undefined static parameters: {undefined}')
            self.__dict__['_stoichio_lambdified'] = lambda : stoichio_arr

    def stoichio_eval(self):
        '''Return the stoichiometric coefficients.'''
        if self._stoichio_lambdified is None: self._lambdify_stoichio()
        return self._stoichio_lambdified()

    @property
    def rate_equations(self):
        '''[pandas.DataFrame] Rate equations.'''
        if all(self._rate_equations):
            rate_eqs = [eq.subs(self._parameters) for eq in self._rate_equations]
            return pd.DataFrame(rate_eqs, index=self.IDs, columns=('rate_equation',))

    @property
    def rate_function(self):
        '''[:class:`MultiKinetics`] The function to evaluate the kinetic rates of the processes.'''
        if self._rate_function is None:
            self._collect_rate_func()
        return self._rate_function

    def set_rate_function(self, k):
        dct = self.__dict__
        if k is None:
            dct['_rate_function'] = None
        elif isinstance(k, MultiKinetics):
            dct['_rate_function'] = k
        elif callable(k):
            dct['_rate_function'] = MultiKinetics(self, function=k)
        else:
            raise TypeError(f'rate_function must be a `MultiKinetics` object, or '
                            f'a function that takes exactly 1 input argument '
                            f'(i.e., an array of state variables), or None, '
                            f'not {type(k)}')

    def _collect_rate_func(self):
        self.__dict__['_rate_function'] = MultiKinetics(self)

    # def rate_eval(self, state_arr):
    #     '''Return the kinetic rates given an array of state variables.'''
    #     return self.rate_function(state_arr)

    @property
    def production_rates(self):
        '''[pandas.DataFrame] The rates of production of the components.'''
        if self._production_rates is None:
            return None
        else:
            rates = [r.subs(self._parameters) for r in self._production_rates]
            return pd.DataFrame(rates, index=self._components.IDs, columns=('rate_of_production',))

    def production_rates_eval(self, state_arr):
        '''Return the rates of production or consumption of the components.'''
        self.params_eval(state_arr)
        M_stoichio = self.stoichio_eval()
        rho_arr = self.rate_function(state_arr)
        return np.dot(M_stoichio.T, rho_arr)

    def subgroup(self, IDs):
        '''Create a new subgroup of :class:`CompiledProcesses` objects.'''
        processes = self[IDs]
        new = Processes(processes)
        new.compile()
        return new

    def index(self, ID):
        '''Return index of specified process.'''
        try: return self._index[ID]
        except KeyError:
            raise UndefinedProcess(ID)

    def indices(self, IDs):
        '''Return indices of multiple processes.'''
        try:
            dct = self._index
            return [dct[i] for i in IDs]
        except KeyError as key_error:
            raise UndefinedProcess(key_error.args[0])

    def __contains__(self, process):
        if isinstance(process, str):
            return process in self.__dict__
        elif isinstance(process, Process):
            return process in self.tuple
        else: # pragma: no cover
            return False

    def __len__(self):
        return self.size

    def __iter__(self):
        return iter(self.tuple)

    def copy(self):
        '''Return a copy.'''
        copy = Processes(self.tuple)
        copy.compile()
        return copy

    def __repr__(self):
        return f"{type(self).__name__}([{', '.join(self.IDs)}])"