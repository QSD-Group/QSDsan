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
from sympy import symbols, Matrix
from sympy.parsing.sympy_parser import parse_expr
import numpy as np
import pandas as pd
    
__all__ = ('Process', 'Processes', 'CompiledProcesses', )

class UndefinedProcess(AttributeError):
    '''AttributeError regarding undefined Component objects.'''
    def __init__(self, ID):
        super().__init__(repr(ID))
        
_load_components = settings.get_default_chemicals
#%%
@chemicals_user        
class Process():
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
        Symbolic parameters in stoichiometry coefficients and/or rate equation.
        The default is None.

    Examples
    --------
    None.

    """
    def __init__(self, ID, reaction, ref_component, rate_equation=None, components=None, 
                 conserved_for=('COD', 'N', 'P', 'charge'), parameters=None):
        self._ID = ID
        self._reaction = reaction
        self._ref_component = ref_component # might be better to set the components first, but this should be fine considering legacy
        self.components = self._load_chemicals(components)
        self.conserved_for = conserved_for
        self._parameters = {p: symbols(p) for p in parameters}
        self._stoichiometry = get_stoichiometric_coeff(
            reaction, self._ref_component, self._components, self._conserved_for, self._parameters)
        self._parse_rate_eq(rate_equation)
                
    def get_conversion_factors(self, as_matrix=False):        
        '''
        return conversion factors (i.e., the 'i\_' attributes of the components) 
        as a numpy.ndarray or a SymPy Matrix.
        '''
        conserved_for = self._conserved_for
        if conserved_for:
            cmps = self._components
            arr = getattr(cmps, 'i_'+conserved_for[0])
            for c in conserved_for[1:]:
                arr = np.vstack((arr, getattr(cmps, 'i_'+c)))
            if as_matrix: return Matrix(arr.tolist())
            return arr
        else: return None
                    
    def check_conservation(self, rtol=1e-5, atol=1e-8):
        '''
        Check conservation of materials (based on the rules set in the `conserved_for` attribute)
        given numerical stoichiometric coefficients.
        
        No return indicates that all rules are satisfied.

        Parameters
        ----------
        rtol : float
            Relative tolerance.
        atol : float
            Absolute tolerance.

        Examples
        --------
        >>> import qsdsan as qs
        >>> cmps = qs.Components.load_default()
        >>> qs.set_thermo(cmps)

        Prompt warning when no conservation rule is set.

        >>> aer = qs.Process('aer', reaction={'S_O2':1}, ref_component='S_O2',
        ...                          rate_equation='KLa*(DOsat - S_O2)',
        ...                          conserved_for=(), parameters=('KLa', 'DOsat'))

        >>> aer.check_conservation() # doctest: +SKIP
        
        Raise error when materials are not conserved.
        
        >>> den1 = qs.Process('den1', reaction={'S_NO3':-1, 'S_NO2':0.5}, ref_component='S_NO3',
        ...                   rate_equation='k1*S_NO3',
        ...                   conserved_for=('N',), parameters=('k1',))
        >>> den1.check_conservation() # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        RuntimeError: The following materials are unconserved by the stoichiometric coefficients.

        No returns when all materials are conserved.
        
        >>> den1.reaction = {'S_NO3':-1, 'S_NO2':1}
        >>> den1.check_conservation()

        See Also
        --------
        `numpy.isclose <https://numpy.org/doc/stable/reference/generated/numpy.isclose.html>`_ for ``rtol`` and ``atol`` settings.
        '''
        isa = isinstance
        if isa(self._stoichiometry, np.ndarray):
            ic = self.get_conversion_factors()
            if ic is None: # no `conserved_for` attribute
                warn('No available `conserved_for` attributes, cannot check conservation.')
                return
            v = self._stoichiometry
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
        else: 
            raise RuntimeError("Can only check conservations with numerical "
                               "stoichiometric coefficients.")
    
    def reverse(self):
        '''reverse the process as to flip the signs of stoichiometric coefficients of all components.'''
        if isinstance(self._stoichiometry, np.ndarray):
            self._stoichiometry = -self._stoichiometry
        else:
            self._stoichiometry = [-v for v in self._stoichiometry]
        self._rate_equation = -self._rate_equation
        
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
            self._normalize_stoichiometry(ref_cmp)
            self._normalize_rate_eq(ref_cmp)

    @property
    def components(self):
        '''
        [class:`CompiledComponents`]
        
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
        '''[dict] Symbolic parameters in stoichiometric coefficients and/or rate equation.'''
        return self._parameters
    
    def append_parameters(self, *new_pars):
        '''append new symbolic parameters'''
        for p in new_pars:
            self._parameters[p] = symbols(p)
    
    def set_parameters(self, **parameters):
        '''
        Set values for symbolic stoichiometric and/or kinetic parameters.

        Parameters
        ----------
        **parameters : float, int, or str
            Parameters and values.
        '''
        self._parameters.update(parameters)
    
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
        which the reference component is reacted or produced in the process.
        '''
        return self._rate_equation.subs(self._parameters)
    
    def _parse_rate_eq(self, eq):
        cmpconc_symbols = {c: symbols(c) for c in self._components.IDs}
        params = self._parameters
        self._rate_equation = parse_expr(str(eq), {**cmpconc_symbols, **params})
    
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
        factor = self._stoichiometry[self._components.index(str(new_ref))]
        self._rate_equation *= factor
    
    def show(self):
        info = f"Process: {self.ID}"
        header = '\n[stoichiometry] '
        section = []
        for cmp, stoichio in tuple(zip(self._components.IDs, self._stoichiometry)):
            if stoichio != 0:
                if isinstance(stoichio, (int, float)): line = f"{cmp}: {stoichio:.3g}"
                else: line = f"{cmp}: {stoichio.evalf(n=3)}"
                section.append(line)
        info += header + ("\n" + 16*" ").join(section)
        info += '\n[reference]     ' + f"{self.ref_component}"
        line = '\n[rate equation] ' + f"{self._rate_equation}"
        if len(line) > 47: line = line[:47] + '...'
        info += line
        header = '\n[parameters]    '
        section = []
        for k,v in self.parameters.items():
            if isinstance(v, (int, float)): line = f"{k}: {v:.3g}"
            else: line = f"{k}: {v}"
            section.append(line)
        info += header + ("\n" + 16*" ").join(section)
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

        Examples
        --------
        >>> import qsdsan as qs
        >>> cmps = qs.Components.load_default()
        >>> qs.set_thermo(cmps)
        >>> aer = qs.Process('aer', reaction={'S_O2':1}, ref_component='S_O2',
        ...                          rate_equation='KLa*(DOsat - S_O2)',
        ...                          conserved_for=(), parameters=('KLa', 'DOsat'))
        >>> aer.show()
        Process: aer
        [stoichiometry] S_O2: 1
        [reference]     S_O2
        [rate equation] KLa*(DOsat - S_O2)
        [parameters]    KLa: KLa
                        DOsat: DOsat
        >>> aer_cp = aer.copy()
        >>> aer_cp.show()
        Process: aer_copy
        [stoichiometry] S_O2: 1
        [reference]     S_O2
        [rate equation] KLa*(DOsat - S_O2)
        [parameters]    KLa: KLa
                        DOsat: DOsat
        '''
        cls = self.__class__
        new = cls.__new__(cls)
        new_ID = new_ID or self.ID+'_copy'
        new.__init__(
            new_ID, reaction=self.reaction, ref_component=self.ref_component,
            rate_equation=self.rate_equation, components=self._components,
            conserved_for=self.conserved_for, parameters=self.parameters
            )
        return new
    __copy__ = copy


#%%
setattr = object.__setattr__
class Processes():
    """
    Create a :class:`Processes` object that contains :class:`Process` objects as attributes.

    Parameters
    ----------
    processes : Iterable[:class:`Process`]

    Examples
    --------
    None.

    """
    
    def __new__(cls, processes):

        self = super().__new__(cls)
        #!!! add function to detect duplicated processes
        setfield = setattr
        for i in processes:
            setfield(self, i.ID, i)
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
    
    def subgroup(self, IDs):
        """
        Create a new subgroup of processes.
        
        Parameters
        ----------
        IDs : Iterable[str]
              Process IDs.
              
        """
        return Process([getattr(self, i) for i in IDs])
    
    def compile(self, skip_checks=False):
        '''Cast as a :class:`CompiledProcesses` object.'''
        processes = tuple(self)
        setattr(self, '__class__', CompiledProcesses)
        try: self._compile(processes, skip_checks)
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
                       conserved_for=('COD', 'N', 'P', 'charge'), parameters=None,
                       use_default_data=False, store_data=False, compile=True):
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
            attributes for the components. The default is ('COD', 'N', 'P', 'charge').
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
        
        cmp_IDs = data.columns[1:-1]
        data.dropna(how='all', subset=cmp_IDs, inplace=True)
        new = cls(())
        for i, proc in data.iterrows():
            ID = proc[0]
            stoichio = proc[1:-1]
            if pd.isna(proc[-1]): rate_eq = None
            else: rate_eq = proc[-1]
            stoichio = stoichio[-pd.isna(stoichio)].to_dict()
            for k,v in stoichio.items():
                try: stoichio[k] = float(v)
                except: continue
            ref = [k for k,v in stoichio.items() if v in (-1, 1)]
            if len(ref) == 0: ref = list(stoichio.keys())[0]                
            else: ref = ref[0]            
            process = Process(ID, stoichio, 
                              ref_component=ref, 
                              rate_equation=rate_eq,
                              components=components,
                              conserved_for=conserved_for,
                              parameters=parameters)
            new.append(process)

        if store_data:
            cls._default_data = data
        
        if compile: new.compile()
        return new
        
            
#%%
@read_only(methods=('append', 'extend', '__setitem__'))
class CompiledProcesses(Processes):
    """
    Create a :class:`CompiledProcesses` object that contains :class:`Process` objects as attributes.

    Parameters
    ----------
    processes : Iterable[:class:`Process`]

    Examples
    --------
    None.

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
        missing_rate_proc = []
        for process in processes:
            if process.rate_equation == None: 
                missing_rate_proc.append(process.ID)
        if not skip_checks and len(missing_rate_proc) > 0: 
            raise RuntimeError(f"The following processes are missing rate equations:"
                               f"{missing_rate_proc}")        
        IDs = tuple_([i.ID for i in processes])
        size = len(IDs)
        index = tuple_(range(size))
        dct['tuple'] = processes
        dct['size'] = size
        dct['IDs'] = IDs
        dct['_index'] = index = dict(zip(IDs, index))
        cmps = Components([cmp for i in processes for cmp in i._components])
        dct['_components'] = _load_components(cmps)
        M_stch = []
        params = {}
        rate_eqs = tuple_([i._rate_equation for i in processes])
        all_numeric = True
        for i in processes:
            stch = [0]*cmps.size
            params.update(i._parameters)
            if all_numeric and isa(i._stoichiometry, (list, tuple)): all_numeric = False
            for cmp, coeff in i.stoichiometry.items():
                stch[cmps._index[cmp]] = coeff
            M_stch.append(stch)
        dct['_parameters'] = params
        if all_numeric: M_stch = np.asarray(M_stch)
        dct['_stoichiometry'] = M_stch
        dct['_rate_equations'] = rate_eqs
        dct['_production_rates'] = list(Matrix(M_stch).T * Matrix(rate_eqs))
        
    @property
    def parameters(self):
        '''[dict] All symbolic stoichiometric and kinetic parameters.'''
        return self._parameters
    
    @property
    def stoichiometry(self):
        '''[pandas.DataFrame] Stoichiometric coefficients.'''
        stoichio = self._stoichiometry
        isa = isinstance
        v_params = {k:v for k,v in self._parameters.items() if isa(v, (float, int))}
        if isa(stoichio, list) and len(v_params) > 0:
            stoichio_vals = []
            for row in stoichio:
                stoichio_vals.append([v.subs(v_params) if not isa(v, (float, int)) else v for v in row])
            return pd.DataFrame(stoichio_vals, index=self.IDs, columns=self._components.IDs)
        else: return pd.DataFrame(stoichio, index=self.IDs, columns=self._components.IDs)
    
    @property
    def rate_equations(self):
        '''[pandas.DataFrame] Rate equations.'''
        rate_eqs = [eq.subs(self._parameters) for eq in self._rate_equations]
        return pd.DataFrame(rate_eqs, index=self.IDs, columns=('rate_equation',))
    
    # @property
    # def rate_equations_function(self):
    #     '''
    #     [function] A function object of the rate equations. When evaluated, 
    #     returns a list of process rates.
    #     '''
    #     return self._rate_equations_function
    
    @property
    def production_rates(self):
        '''[pandas.DataFrame] The rates of production of the components.'''
        rates = [r.subs(self._parameters) for r in self._production_rates]
        return pd.DataFrame(rates, index=self._components.IDs, columns=('rate_of_production',))
    
    # @property
    # def production_rates_function(self):
    #     '''
    #     [function] A function object of the components' rates of production. 
    #     When evaluated, returns a list of production rates.
    #     '''
    #     args = list(symbols(self._components.IDs)) + [v for k,v in self._parameters.items() if k == str(v)]
    #     return [lambdify(args, r) for r in self._production_rates]
    
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
    
    def set_parameters(self, **parameters):
        '''Set values to stoichiometric and/or kinetic parameters.'''
        self._parameters.update(parameters)
    
    def __repr__(self):
        return f"{type(self).__name__}([{', '.join(self.IDs)}])"