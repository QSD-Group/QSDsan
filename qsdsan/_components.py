#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

Part of this module is based on the Thermosteam package:
https://github.com/BioSTEAMDevelopmentGroup/thermosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np
import pandas as pd
import thermosteam as tmo
from . import _component, Chemical, Chemicals, CompiledChemicals, Component
from .utils import load_data

__all__ = ('Components', 'CompiledComponents')

setattr = object.__setattr__

utils = tmo.utils
_component_properties = _component._component_properties
_num_component_properties = _component._num_component_properties
_key_component_properties = _component._key_component_properties
# _TMH = tmo.base.thermo_model_handle.ThermoModelHandle
_PH = tmo.base.phase_handle.PhaseHandle
DomainError = tmo.exceptions.DomainError


# %%

class UndefinedComponent(AttributeError):
    '''AttributeError regarding undefined :class:`Component` objects.'''
    def __init__(self, ID):
        super().__init__(repr(ID))

def must_compile(*args, **kwargs): # pragma: no cover
    raise TypeError('Method valid only for CompiledChemicals, '
                    'run <Components>.compile() to compile first.')


# =============================================================================
# Define the Components class
# =============================================================================

class Components(Chemicals):
    '''
    A subclass of :class:`thermosteam.Chemicals`, contains :class:`Component`
    objects as attributes.

    Examples
    --------
    `Component <https://qsdsan.readthedocs.io/en/latest/tutorials/Component.html>`_

    See Also
    --------
    `thermosteam.Chemicals <https://thermosteam.readthedocs.io/en/latest/Chemicals.html>`_

    '''

    def __new__(cls, components, cache=False):
        self = super(Chemicals, cls).__new__(cls)
        isa = isinstance
        setfield = setattr
        IDs = set()
        # CASs = set()
        for i in components:
            if isa(i, Component):
                ID = i.ID
                if ID in IDs:
                    raise ValueError(f'More than one `Component` has the ID {ID}.')
                # CAS = i.CAS
                # if CAS in CASs: continue
                # CASs.add(CAS)
                setfield(self, i.ID, i)
            elif isa(i, Chemical):
                raise TypeError(f'{i} is a `thermosteam.Chemical` object, '
                                'use `Component.from_chemical` to define a `Component` object.')
            else:
                raise TypeError(f'Only `Component` objects can be included, not a `{type(i).__name__}` object.')

        return self


    def __setattr__(self, ID, component):
        raise TypeError('Cannot set attribute; use `<Components>.append/extend` instead.')


    def __setitem__(self, ID, component):
        raise TypeError('Cannot set item; use `<Components>.append/extend` instead.')


    def __getitem__(self, key):
        '''Return a :class:`Component` object or a list of :class:`Component` objects.'''
        dct = self.__dict__
        try:
            if isinstance(key, str):
                return dct[key]
            else:
                return [dct[i] for i in key]
        except KeyError as key_error:
            raise UndefinedComponent(key_error.args[0])


    def __contains__(self, component):
        if isinstance(component, str):
            return component in self.__dict__
        elif isinstance(component, Component):
            return component in self.__dict__.values()
        else: # pragma: no cover
            return False


    def __repr__(self):
        return f"Components([{', '.join(self.__dict__)}])"


    def copy(self):
        '''Return a copy.'''
        copy = object.__new__(Components)
        for cmp in self: setattr(copy, cmp.ID, cmp)
        return copy


    def append(self, component):
        '''Append a Component'''
        if not isinstance(component, Component):
            if isinstance(component, Chemical):
                raise TypeError(f'{component} is a `Chemical` object, '
                                'use `Component.from_chemical` to define a `Component` object.')
            else:
                raise TypeError("only `Component` objects can be appended, "
                               f"not `{type(component).__name__}` object.")
        ID = component.ID
        if ID in self.__dict__:
            raise ValueError(f"{ID} already defined in this `Components` object.")
        setattr(self, ID, component)


    def extend(self, components):
        '''Extend with more :class:`Component` objects.'''
        if isinstance(components, Components):
            self.__dict__.update(components.__dict__)
        else:
            for component in components: self.append(component)


    def compile(self, skip_checks=False):
        '''Cast as a :class:`CompiledComponents` object.'''
        components = tuple(self)
        setattr(self, '__class__', CompiledComponents)
        try: self._compile(components, skip_checks)
        except Exception as error:
            setattr(self, '__class__', Components)
            setattr(self, '__dict__', {i.ID: i for i in components})
            raise error

    kwarray = array = index = indices = must_compile

    _default_data = None


    @classmethod
    def load_from_file(cls, path_or_df, index_col=None,
                       use_default_data=False, store_data=False):
        '''
        Create and return a :class:`Components` objects based on properties
        defined in a datasheet.

        Parameters
        ----------
        path_or_df : str or :class:`pandas.DataFrame`
            File path, the file should end with ".cvs", ".xls", or "xlsx".
        index_col : None or int
            Index column of the :class:`pandas.DataFrame`.
        use_default_data : bool
            Whether to use the default components.
        store_data : bool
            Whether to store this as the default components.

        Returns
        -------
        A :class:`Components` object that contains all created Component objects.


        .. note::

            The :class:`Components` object needs to be compiled before it is used in simulation.

        '''
        if use_default_data and cls._default_data is not None:
            data = cls._default_data
        elif isinstance(path_or_df, str):
            data = load_data(path_or_df, index_col=index_col)
        else:
            data = path_or_df

        new = cls(())

        for i, cmp in data.iterrows():
            if pd.isna(cmp.measured_as):
                cmp.measured_as = None
            try:
                component = Component(ID = cmp.ID,
                                      search_ID = str(cmp.CAS),
                                      measured_as = cmp.measured_as)
            except LookupError:
                try:
                    component = Component(ID = cmp.ID,
                                          search_ID = 'PubChem='+str(int(cmp.PubChem)),
                                          measured_as = cmp.measured_as)
                except:
                    if not pd.isna(cmp.formula):
                        component = Component(ID = cmp.ID,
                                              formula = cmp.formula,
                                              measured_as = cmp.measured_as)
                    else:
                        component = Component(ID = cmp.ID,
                                              measured_as = cmp.measured_as)
            for j in _component_properties:
                field = '_' + j
                if pd.isna(cmp[j]): setattr(component, j, None)
                else: setattr(component, field, cmp[j])
            new.append(component)

        if store_data:
            cls._default_data = data
        return new


    @classmethod
    def load_default(cls, use_default_data=True, store_data=True, default_compile=True):
        '''
        Create and return a :class:`Components` or :class:`CompiledComponents`
        object containing all default :class:`Component` objects.

        Parameters
        ----------
        use_default_data : bool, optional
            Whether to use default cache data. The default is True.
        store_data : bool, optional
            Whether to store the default data as cache. The default is True.
        default_compile : bool, optional
            Whether to compile the default :class:`Components`. The default is True.

        Returns
        -------
        A :class:`Components` or :class:`CompiledComponents` object with
        default :class:`Component` objects.


        .. note::

            [1] Component-specific properties are defined in ./data/component.cvs.

            [2] When `default_compile` is True, all essential chemical-specific properties
            (except molar volume model and normal boiling temperature) that are missing will
            be defaulted to those of water.

            [3] When `default_compile` is True, missing molar volume models will be defaulted
            according to particle sizes: particulate or colloidal -> 1.2e-5 m3/mol,
            soluble -> copy from urea, dissolved gas -> copy from CO2.

            [4] When `default_compile` is True, missing normal boiling temoerature will be
            defaulted according to particle sizes: particulate or colloidal -> copy from NaCl,
            soluble -> copy from urea, dissolved gas -> copy from CO2.


        '''
        import os
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/_components.tsv')
        del os
        new = cls.load_from_file(path, index_col=None, use_default_data=True, store_data=True)

        H2O = Component.from_chemical('H2O', Chemical('H2O'),
                                      i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                                      f_Vmass_Totmass=0, description="Water",
                                      particle_size='Soluble',
                                      degradability='Undegradable', organic=False)
        new.append(H2O)

        if default_compile:
            isa = isinstance
            for i in new:
                i.default()

                if i.particle_size == 'Soluble':
                    ref_chem = Chemical('urea')
                elif i.particle_size == 'Dissolved gas':
                    ref_chem = Chemical('CO2')
                else:
                    ref_chem = Chemical('NaCl')

                i.Tb = ref_chem.Tb if not i.Tb else i.Tb

                COPY_V = False # if don't have V model, set those to default
                if isa(i.V, _PH):
                    if not i.V.l.valid_methods(298.15): COPY_V = True
                else:
                    if not i.V.valid_methods(298.15): COPY_V = True

                if COPY_V:
                    if i.particle_size in ('Soluble', 'Dissolved gas'):
                        i.copy_models_from(ref_chem, names=('V',))
                    else:
                        try: i.V.add_model(1.2e-5)
                        except AttributeError:
                            i.V.l.add_model(1.2e-5)    # m^3/mol
                            i.V.s.add_model(1.2e-5)

                i.copy_models_from(H2O)

                try:
                    i.Hvap(i.Tb)
                except RuntimeError: # Hvap model of H2O cannot be extrapolated to Tb
                    i.copy_models_from(ref_chem, names=('Hvap',))

            new.compile()
        return new


    @classmethod
    def from_chemicals(cls, chemicals, **data):
        '''
        Return a new :class:`Components` from a :class:`thermosteam.Chemicals`
        or :class:`thermosteam.CompiledChemicals` object.

        Parameters
        ----------
        chemicals: thermosteam.Chemicals
            The :class:`thermosteam.Chemicals` object as the basis
            for the new :class:`~.Components` object.
            :class:`Component` objects will have the same ID as the corresponding
            :class:`thermosteam.Chemical` object in the :class:`thermosteam.Chemicals`
            object.
        data : dict
            A nested dict with keys being the new components and values being the inner dict,
            keys and values of the inner dict are the attribute names and values, respectively.

        Examples
        --------
        >>> import qsdsan as qs
        >>> chems = qs.Chemicals((qs.Chemical('Water'), qs.Chemical('Ethanol')))
        >>> data = {'Water': {'particle_size': 'Soluble',
        ...                   'degradability': 'Undegradable',
        ...                   'organic': False},
        ...         'Ethanol': {'particle_size': 'Soluble',
        ...                     'degradability': 'Readily',
        ...                     'organic': False}}
        >>> cmps = qs.Components.from_chemicals(chems, **data)
        >>> cmps
        Components([Water, Ethanol])
        '''
        cmps = cls.__new__(cls, ())
        for i in chemicals.__iter__():
            val_dct = data.get(i.ID)
            cmp = Component.from_chemical(i.ID, i)
            if val_dct:
                for k, v in val_dct.items():
                    setattr(cmp, k, v)
            cmps.append(cmp)

        return cmps


# %%

# =============================================================================
# Define the CompiledComponents class
# =============================================================================

chemical_data_array = tmo._chemicals.chemical_data_array

def component_data_array(components, attr):
    data = chemical_data_array(components, attr)
    return data


class CompiledComponents(CompiledChemicals):
    '''
    A subclass of :class:`thermosteam.CompiledChemicals`, contains `Component` objects as attributes.

    Examples
    --------
    `Component <https://qsdsan.readthedocs.io/en/latest/tutorials/Component.html>`_

    See Also
    --------
    `thermosteam.CompiledChemicals <https://thermosteam.readthedocs.io/en/latest/Chemicals.html>`_

    '''

    _cache = {}

    def __new__(cls, components, cache=None):
        isa = isinstance
        components = tuple([cmp if isa(cmp, Component) else Component(cmp, cache)
                           for cmp in components])
        cache = cls._cache
        if components in cache:
            self = cache[components]
        else:
            self = object.__new__(cls)
            setfield = setattr
            for cmp in components:
                setfield(self, cmp.ID, cmp)
            self._compile(components)
            cache[components] = self
        return self

    def __reduce__(self):
        return CompiledComponents, (self.tuple, )


    def __contains__(self, component):
        if isinstance(component, str):
            return component in self.__dict__
        elif isinstance(component, Component):
            return component in self.tuple
        else: # pragma: no cover
            return False


    def __repr__(self):
        return f"CompiledComponents([{', '.join(self.IDs)}])"


    def refresh_constants(self):
        '''
        Refresh constant arrays of :class:`Components` objects,
        including all chemical and component-specific properties.
        '''
        super().refresh_constants()
        dct = self.__dict__
        components = self.tuple
        for i in _num_component_properties:
            dct[i] = component_data_array(components, i)

    def compile(self, skip_checks=False):
        '''Skip, :class:`CompiledComponents` have already been compiled.'''
        pass


    def _compile(self, components, skip_checks=False):
        dct = self.__dict__
        tuple_ = tuple # this speeds up the code
        components = tuple_(dct.values())
        CompiledChemicals._compile(self, components, skip_checks)
        for component in components:
            missing_properties = component.get_missing_properties(_key_component_properties)
            if not missing_properties: continue
            missing = utils.repr_listed_values(missing_properties)
            raise RuntimeError(f'{component} is missing key component-related properties ({missing}).')

        for i in _num_component_properties:
            dct[i] = component_data_array(components, i)

        dct['s'] = np.asarray([1 if cmp.particle_size == 'Soluble' else 0 for cmp in components])
        dct['c'] = np.asarray([1 if cmp.particle_size == 'Colloidal' else 0 for cmp in components])
        dct['x'] = np.asarray([1 if cmp.particle_size == 'Particulate' else 0 for cmp in components])
        dct['b'] = np.asarray([1 if cmp.degradability != 'Undegradable' else 0 for cmp in components])
        dct['rb'] = np.asarray([1 if cmp.degradability == 'Readily' else 0 for cmp in components])
        dct['org'] = np.asarray([int(cmp.organic) for cmp in components])


    def subgroup(self, IDs):
        '''Create a new subgroup of :class:`Component` objects.'''
        components = self[IDs]
        new = Components(components)
        new.compile()
        for i in new.IDs:
            for j in self.get_synonyms(i):
                try: new.set_synonym(i, j)
                except: pass
        return new


    def index(self, ID):
        '''Return index of specified component.'''
        try: return self._index[ID]
        except KeyError:
            raise UndefinedComponent(ID)


    def indices(self, IDs):
        '''Return indices of multiple components.'''
        try:
            dct = self._index
            return [dct[i] for i in IDs]
        except KeyError as key_error:
            raise UndefinedComponent(key_error.args[0])


    def copy(self):
        '''Return a copy.'''
        copy = Components(self)
        copy.compile()
        return copy