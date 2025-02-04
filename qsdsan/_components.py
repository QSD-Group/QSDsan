#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>
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
from .utils import add_V_from_rho, load_data

__all__ = ('Components', 'CompiledComponents')

setattr = object.__setattr__

utils = tmo.utils
_component_properties = _component._component_properties
_num_component_properties = _component._num_component_properties
_key_component_properties = _component._key_component_properties
# _TMH = tmo.base.thermo_model_handle.ThermoModelHandle
_PH = tmo.base.phase_handle.PhaseHandle


# %%

class UndefinedComponent(tmo.exceptions.UndefinedChemical):
    '''AttributeError regarding undefined :class:`Component` objects.'''


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
    `Component <https://qsdsan.readthedocs.io/en/latest/tutorials/2_Component.html>`_

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


    def compile(self, skip_checks=False, ignore_inaccurate_molar_weight=False, 
                adjust_MW_to_measured_as=False):
        '''Cast as a :class:`CompiledComponents` object.'''
        components = tuple(self)
        tmo._chemicals.prepare(components, skip_checks)
        setattr(self, '__class__', CompiledComponents)
        
        try: self._compile(components, ignore_inaccurate_molar_weight, adjust_MW_to_measured_as)
        except Exception as error:
            setattr(self, '__class__', Components)
            setattr(self, '__dict__', {i.ID: i for i in components})
            raise error

    kwarray = array = index = indices = must_compile

    _default_data = None


    def default_compile(self, lock_state_at='l',
                        soluble_ref='Urea', gas_ref='CO2', particulate_ref='Stearin', 
                        ignore_inaccurate_molar_weight=False,
                        adjust_MW_to_measured_as=False):
        '''
        Auto-fill of the missing properties of the components and compile,
        boiling point (Tb) and molar volume (V) will be copied from the reference component,
        the remaining missing properties will be copied from those of water.

        Parameters
        ----------
        lock_state_at : str
            Lock the state of components at a certain phase,
            can be 'g', 'l', 's', or left as empty to avoid locking state.
            Components that have already been locked will not be affected.
        soluble_ref : obj or str
            Reference component (or chemical ID) for those with `particle_size` == 'Soluble'.
        gas_ref : obj or str
            Reference component (or chemical ID) for those with `particle_size` == 'Dissolved gas'.
        particulate_ref : obj or str
            Reference component (or chemical ID) for those with `particle_size` == 'Particulate'.
        ignore_inaccurate_molar_weight : bool
            Default is False, need to be manually set to True if having components
            with `measured_as` attributes. This is to alert the users that
            calculations for attributes using molecular weight will be inacurate,
            unless all components have sensible MWs and `adjust_MW_to_measured_as`
            is set to True.
        adjust_MW_to_measured_as : bool
            Default is False. Manually set it to True to adjust the MW data of
            components with `measured_as` attributes and chemical formulas. This 
            is to enable correct calculations of component molar flows when possible.
            For components without a sensible MW, their MWs will remain 1 by default.
            This is pertinent for calculations of molar flows and any thermodynamic
            property of a stream.

        Examples
        --------
        >>> from qsdsan import Component, Components, Stream, set_thermo
        >>> X = Component('X', phase='s', measured_as='COD', i_COD=0, description='Biomass',
        ...               organic=True, particle_size='Particulate', degradability='Readily')
        >>> X_inert = Component('X_inert', phase='s', description='Inert biomass', i_COD=0,
        ...                     organic=True, particle_size='Particulate', degradability='Undegradable')
        >>> Substrate = Component('Substrate', phase='s', measured_as='COD', i_mass=18.3/300,
        ...                       organic=True, particle_size='Particulate', degradability='Readily')
        >>> cmps = Components([X, X_inert, Substrate])
        >>> # As none of the components above has a chemical formula, i.e., no sensible MW, 
        >>> # simply set `ignore_inaccurate_molar_weight` to True to bypass error.
        >>> cmps.default_compile(ignore_inaccurate_molar_weight=True)
        >>> cmps
        CompiledComponents([X, X_inert, Substrate])
        
        >>> Ac = Component('Ac', search_ID='CH3COOH', particle_size='Soluble', 
        ...                degradability='Readily', organic=True)
        >>> Ac_as_COD = Component('Ac_as_COD', search_ID='CH3COOH', measured_as='COD', 
        ...                       particle_size='Soluble', degradability='Readily', organic=True)
        >>> Acs = Components([Ac, Ac_as_COD])
        >>> Acs.default_compile(ignore_inaccurate_molar_weight=True, 
        ...                     adjust_MW_to_measured_as=False)
        >>> set_thermo(Acs)
        >>> # Create a WasteStream object with 1 kmol/hr of each acetic acid component, 
        >>> # knowing 1 mol acetate is equivalent to 2 mol of O2 demand
        >>> s1 = Stream('s1', Ac=60.052, Ac_as_COD=2 * 32, units='kg/hr')
        >>> s1.mass
        sparse([60.052, 64.   ])
        
        >>> # However, the calculated molar flow is inaccurate because the MW for Ac_as_COD
        >>> # is inconsistent with its `measured_as`. This also affects the
        >>> # calculation of other thermodynamic properties.
        >>> s1.mol
        sparse([1.   , 1.066])
        >>> s1.vol
        sparse([0.05 , 0.054])
        
        >>> # To fix the molar flow calculation, simply set `adjust_MW_to_measured_as` to True when compile.
        >>> Acs_MW_adjusted = Components([Ac, Ac_as_COD])
        >>> Acs_MW_adjusted.default_compile(adjust_MW_to_measured_as=True)
        >>> set_thermo(Acs_MW_adjusted)
        >>> s2 = Stream('s2', Ac=60.052, Ac_as_COD=2 * 32, units='kg/hr')
        >>> s2.mol
        sparse([1., 1.])
        >>> s2.vol
        sparse([0.05, 0.05])
        
        '''
        isa = isinstance
        get = getattr
        if isa(soluble_ref, str):
            sol = Chemical(soluble_ref)
            if soluble_ref.lower() == 'urea': sol.at_state('l')
        if isa(gas_ref, str):
            gas = Chemical(gas_ref)
            if gas_ref.lower() == 'co2': gas.at_state('g')
        if isa(particulate_ref, str):
            par = Chemical(particulate_ref)
            if particulate_ref.lower() == 'stearin':
                # 0.8559 at 90 °C https://pubchem.ncbi.nlm.nih.gov/compound/Tristearin#section=Density
                # avg ~0.9 http://www.dgfett.de/material/physikalische_eigenschaften.pdf
                add_V_from_rho(par, 0.9, 'g/ml', 's')
                par.at_state('s')

        TPkwargs = dict(T=298.15, P=101325)
        def get_constant_V_model(ref_cmp, phase=''):
            if ref.locked_state: return ref.V(**TPkwargs)
            else: return ref.V(phase, **TPkwargs)

        water = Chemical('Water')
        for cmp in self:
            particle_size = cmp.particle_size
            ref = sol if particle_size=='Soluble' \
                else gas if particle_size=='Dissolved gas' else par
            if lock_state_at:
                try: cmp.at_state(lock_state_at)
                except TypeError: pass
            cmp.Tb = cmp.Tb or ref.Tb

            # If don't have model for molar volume, set those to default
            COPY_V = False
            cmp_V = cmp.V if cmp.locked_state else cmp.V.l
            try: cmp_V(**TPkwargs)
            except: COPY_V = True
            if COPY_V:
                locked_state = cmp.locked_state
                if locked_state:
                    V_const = get_constant_V_model(ref, locked_state)*(cmp.chem_MW/ref.MW)
                    cmp.V.add_model(V_const)
                else:
                    for phase in ('g', 'l', 's'): # iterate through phases
                        backup_ref = gas if phase=='g' else sol if phase=='l' else par
                        try: V_const = get_constant_V_model(ref, phase)
                        except: V_const = get_constant_V_model(backup_ref, phase)
                        V_const *= (cmp.chem_MW/ref.MW)
                        get(cmp.V, phase).add_model(V_const)

            if not cmp.Hvap.valid_methods():
                try:
                    ref.Hvap(cmp.Tb)
                    cmp.copy_models_from(ref, names=('Hvap', ))
                except RuntimeError: # Hvap model cannot be extrapolated to Tb
                    cmp.copy_models_from(water, names=('Hvap', ))

            # Copy all remaining properties from water
            cmp.copy_models_from(water)
        for cmp in self: cmp.default()
        self.compile(ignore_inaccurate_molar_weight=ignore_inaccurate_molar_weight,
                     adjust_MW_to_measured_as=adjust_MW_to_measured_as)


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
            Whether to use the cached default components.
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
                try:
                    if pd.isna(cmp[j]): setattr(component, j, None)
                    else: setattr(component, field, cmp[j])
                except KeyError:
                    continue
            new.append(component)

        if store_data:
            cls._default_data = data
        return new


    @classmethod
    def load_default(cls, use_default_data=True, store_data=False, default_compile=True):
        '''
        Create and return a :class:`Components` or :class:`CompiledComponents`
        object containing all default :class:`Component` objects based on
        `Reiger et al. <https://iwaponline.com/ebooks/book/630/Guidelines-for-Using-Activated-Sludge-Models>`_

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
            according to particle sizes: particulate or colloidal -> copy from NaCl,
            soluble -> copy from urea, dissolved gas -> copy from CO2.

            [4] When `default_compile` is True, missing normal boiling temperature will be
            defaulted according to particle sizes: particulate or colloidal -> copy from NaCl,
            soluble -> copy from urea, dissolved gas -> copy from CO2.


        See Also
        --------
        :func:`~.Components.default_compile`

        References
        ----------
        [1] Rieger, L.; Gillot, S.; Langergraber, G.; Ohtsuki, T.; Shaw, A.; Tak´cs, I.; Winkler, S.
        Guidelines for Using Activated Sludge Models; IWA Publishing, 2012.
        https://doi.org/10.2166/9781780401164.
        '''
        import os
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data/_components.tsv')
        del os
        new = cls.load_from_file(path, index_col=None,
                                 use_default_data=use_default_data, store_data=store_data)

        H2O = Component.from_chemical('H2O', Chemical('H2O'),
                                      i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                                      f_Vmass_Totmass=0, description="Water",
                                      particle_size='Soluble',
                                      degradability='Undegradable', organic=False)
        new.append(H2O)

        if default_compile:
            new.default_compile(lock_state_at='', particulate_ref='NaCl', ignore_inaccurate_molar_weight=True)
            new.compile(ignore_inaccurate_molar_weight=True)
            # Add aliases
            new.set_alias('H2O', 'Water')
            # Pre-define groups
            new.define_group('substrates',
                             ('S_CH3OH', 'S_Ac', 'S_Prop', 'S_F', 'C_B_Subst', 'X_B_Subst'))
            new.define_group('biomass',
                             ('X_OHO', 'X_AOO', 'X_NOO', 'X_AMO', 'X_PAO',
                              'X_MEOLO', 'X_FO', 'X_ACO', 'X_HMO', 'X_PRO'))
            new.define_group('S_VFA', ('S_Ac', 'S_Prop'))
            new.define_group('X_Stor', ('X_OHO_PHA', 'X_GAO_PHA', 'X_PAO_PHA',
                                        'X_GAO_Gly', 'X_PAO_Gly'),)
            new.define_group('X_ANO', ('X_AOO', 'X_NOO'))
            new.define_group('X_Bio', ('X_OHO', 'X_AOO', 'X_NOO', 'X_AMO', 'X_PAO',
                                       'X_MEOLO', 'X_ACO', 'X_HMO', 'X_PRO', 'X_FO'))
            new.define_group('S_NOx', ('S_NO2', 'S_NO3'))
            new.define_group('X_PAO_PP', ('X_PAO_PP_Lo', 'X_PAO_PP_Hi'))
            new.define_group('TKN', [i.ID for i in new if i.ID not in ('S_N2','S_NO2','S_NO3')])
        return new


    @staticmethod
    def append_combustion_components(components, alt_IDs={},
                                     try_default_compile=True,
                                     **default_compile_kwargs):
        '''
        Return a new :class:`~.Components` object with the given components
        and those needed for combustion reactions (complete oxidation with O2),
        namely O2, CO2 (for C), H2O (for H), N2 (for N), P4O10 (for P), and SO2 (for S).

        If the combustion components are already in the given collection,
        they will NOT be overwritten.

        Parameters
        ----------
        components : Iterable(obj)
            The original components to be appended.
        alt_IDs : dict
            Alternative IDs for the combustion components to be added as aliases,
            e.g., if "S_O2" is used instead of "O2", then pass {'O2': 'S_O2'}.
        default_compile : bool
            Whether to try default compile when some components
            are missing key properties for compiling.
        default_compile_kwargs : dict
            Keyword arguments to pass to `default_compile` if needed.

        See Also
        --------
        :func:`default_compile`

        Examples
        --------
        >>> from qsdsan import Components
        >>> cmps = Components.load_default()
        >>> cmps
        CompiledComponents([S_H2, S_CH4, S_CH3OH, S_Ac, S_Prop, S_F, S_U_Inf, S_U_E, C_B_Subst, C_B_BAP, C_B_UAP, C_U_Inf, X_B_Subst, X_OHO_PHA, X_GAO_PHA, X_PAO_PHA, X_GAO_Gly, X_PAO_Gly, X_OHO, X_AOO, X_NOO, X_AMO, X_PAO, X_MEOLO, X_FO, X_ACO, X_HMO, X_PRO, X_U_Inf, X_U_OHO_E, X_U_PAO_E, X_Ig_ISS, X_MgCO3, X_CaCO3, X_MAP, X_HAP, X_HDP, X_FePO4, X_AlPO4, X_AlOH, X_FeOH, X_PAO_PP_Lo, X_PAO_PP_Hi, S_NH4, S_NO2, S_NO3, S_PO4, S_K, S_Ca, S_Mg, S_CO3, S_N2, S_O2, S_CAT, S_AN, H2O])
        >>> CH4 = cmps.S_CH4.copy('CH4', phase='g')
        >>> cmps = Components.append_combustion_components([*cmps, CH4], alt_IDs=dict(O2='S_O2'))
        >>> cmps
        CompiledComponents([S_H2, S_CH4, S_CH3OH, S_Ac, S_Prop, S_F, S_U_Inf, S_U_E, C_B_Subst, C_B_BAP, C_B_UAP, C_U_Inf, X_B_Subst, X_OHO_PHA, X_GAO_PHA, X_PAO_PHA, X_GAO_Gly, X_PAO_Gly, X_OHO, X_AOO, X_NOO, X_AMO, X_PAO, X_MEOLO, X_FO, X_ACO, X_HMO, X_PRO, X_U_Inf, X_U_OHO_E, X_U_PAO_E, X_Ig_ISS, X_MgCO3, X_CaCO3, X_MAP, X_HAP, X_HDP, X_FePO4, X_AlPO4, X_AlOH, X_FeOH, X_PAO_PP_Lo, X_PAO_PP_Hi, S_NH4, S_NO2, S_NO3, S_PO4, S_K, S_Ca, S_Mg, S_CO3, S_N2, S_O2, S_CAT, S_AN, H2O, CH4, CO2, N2, P4O10, SO2])
        >>> cmps.O2 is cmps.S_O2
        True
        '''
        cmps = components if isinstance(components, (Components, CompiledComponents)) \
            else Components(components)
        comb_cmps = ['O2', 'CO2', 'H2O', 'N2', 'P4O10', 'SO2']
        aliases = dict(H2O='Water')
        aliases.update(alt_IDs)
        get = getattr
        for k, v in alt_IDs.items():
            try:
                get(cmps, v)
                aliases[k] = comb_cmps.pop(comb_cmps.index(k))
            except AttributeError:
                pass
        for ID in comb_cmps:
            try: get(cmps, ID)
            except AttributeError:
                phase = 'g' if ID in ('O2', 'CO2', 'N2', 'SO2') else 's' if ID=='P4O10' else ''
                ps = 'Dissolved gas' if phase == 'g' else 'Particulate' if phase=='s' else 'Soluble'
                cmp = Component(ID, phase=phase, organic=False, particle_size=ps,
                                degradability='Undegradable')
                if isinstance(cmps, CompiledComponents): cmps = Components(cmps)
                cmps.append(cmp)
        add_V_from_rho(cmps.P4O10, rho=2.39, rho_unit='g/mL') # http://www.chemspider.com/Chemical-Structure.14128.html
        try:
            for cmp in cmps: cmp.default()
            cmps.compile(ignore_inaccurate_molar_weight=True)
        except RuntimeError: # cannot compile due to missing properties
            cmps.default_compile(ignore_inaccurate_molar_weight=True, **default_compile_kwargs)
        for k, v in aliases.items():
            cmps.set_alias(k, v)
        return cmps


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
        for i in chemicals:
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
    `Component <https://qsdsan.readthedocs.io/en/latest/tutorials/2_Component.html>`_

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

    def compile(self, skip_checks=False, ignore_inaccurate_molar_weight=False):
        '''Do nothing, :class:`CompiledComponents` have already been compiled.'''
        pass


    def _compile(self, components, ignore_inaccurate_molar_weight, adjust_MW_to_measured_as):
        dct = self.__dict__
        tuple_ = tuple # this speeds up the code
        components = tuple_(dct.values())
        CompiledChemicals._compile(self, components)
        for component in components:
            missing_properties = component.get_missing_properties(_key_component_properties)
            if component.measured_as:
                inaccurate = True
                if component.formula and adjust_MW_to_measured_as:
                    component._MW = component.chem_MW / component.i_mass
                    inaccurate = False
                if (not ignore_inaccurate_molar_weight) and inaccurate: 
                    raise RuntimeError(f'{component} does not have a sensible molar weight. Set ignore_inaccurate_molar_weight=True to bypass this error.')
            if not missing_properties: continue
            missing = utils.repr_listed_values(missing_properties)
            raise RuntimeError(f'{component} is missing key component-related properties ({missing}).')

        if adjust_MW_to_measured_as:
            dct['MW'] = component_data_array(components, 'MW')

        for i in _num_component_properties:
            dct[i] = component_data_array(components, i)

        dct['g'] = np.asarray([1 if cmp.particle_size == 'Dissolved gas' else 0 for cmp in components])
        s = dct['s'] = np.asarray([1 if cmp.particle_size == 'Soluble' else 0 for cmp in components])
        c = dct['c'] = np.asarray([1 if cmp.particle_size == 'Colloidal' else 0 for cmp in components])
        dct['x'] = np.asarray([1 if cmp.particle_size == 'Particulate' else 0 for cmp in components])
        b = dct['b'] = np.asarray([1 if cmp.degradability != 'Undegradable' else 0 for cmp in components])
        dct['rb'] = np.asarray([1 if cmp.degradability == 'Readily' else 0 for cmp in components])
        org = dct['org'] = np.asarray([int(cmp.organic) for cmp in components])
        inorg = dct['inorg'] = np.ones_like(org) - org
        ID_arr = dct['_ID_arr'] = np.asarray([i.ID for i in components])
        dct['chem_MW'] = component_data_array(components, 'chem_MW')

        # Inorganic degradable non-gas, incorrect
        inorg_b = inorg * b * (s+c)
        if inorg_b.sum() > 0:
            bad_IDs = ID_arr[np.where(inorg_b==1)[0]]
            raise ValueError(f'Components {bad_IDs} are inorganic, degradable, and not gas, '
                             'which is not correct.')


    def subgroup(self, IDs):
        '''Create a new subgroup of :class:`Component` objects.'''
        components = self[IDs]
        new = Components(components)
        new.compile(ignore_inaccurate_molar_weight=True)
        for i in new.IDs:
            for j in self.get_aliases(i):
                try: new.set_alias(i, j)
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
        copy.compile(ignore_inaccurate_molar_weight=True)
        return copy


    def get_IDs_from_array(self, array):
        '''
        Get the IDs of a group of components based on the 1/0 or True/False array.

        Parameters
        ----------
        array : Iterable(1/0)
            1D collection of 1/0 or True/False with the same length
            as the IDs.

        Examples
        --------
        >>> from qsdsan import Components
        >>> cmps = Components.load_default()
        >>> cmps.get_IDs_from_array(cmps.g)
        ('S_H2', 'S_CH4', 'S_N2', 'S_O2')
        '''
        return tuple(self._ID_arr[np.asarray(array).astype(bool)])


    def get_array_from_IDs(self, IDs):
        '''
        Generate a ``numpy`` array in the same shape as ``CompiledComponents.IDs``,
        where the values would be 1 for components whose IDs are in the given ID iterable
        and 0 for components not in the given ID iterable.

        Parameters
        ----------
        IDs : Iterable(str)
            IDs of select components within this ``qsdsan.CompiledComponents``.

        Examples
        --------
        >>> from qsdsan import Components
        >>> cmps = Components.load_default()
        >>> IDs = ('S_H2', 'S_CH4', 'S_N2', 'S_O2')
        >>> cmps.get_array_from_IDs(IDs)
        array([1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0])
        '''
        arr = np.zeros_like(self._ID_arr, dtype='int')
        arr[self.indices(IDs)] = 1
        return arr


    def remove_alias(self, component, alias):
        '''
        Remove the alias of a component.

        Parameters
        ----------
        component : str or obj
            The component (or the ID of which) whose alias will be removed.
        alias : str
            The alias of the component to be removed.

        Examples
        --------
        >>> from qsdsan.utils import create_example_components
        >>> cmps = create_example_components()
        >>> cmps.H2O is cmps.Water
        True
        >>> cmps.remove_alias(cmps.H2O, 'Water')
        >>> cmps.Water # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        AttributeError: 'CompiledComponents' object has no attribute 'Water'
        '''
        if isinstance(component, str): component = getattr(self, component)
        if alias not in component.aliases:
            raise ValueError(f'The component "{component.ID}" does not have the alias "{alias}".')
        self._index.pop(alias)
        self.__dict__.pop(alias)
    remove_synonym = remove_alias


    @property
    def aliases(self):
        '''All of the aliases of the components.'''
        return set(sum([list(cmp.aliases) for cmp in self], []))

    @property
    def names(self):
        '''All of the names and aliases of the components.'''
        return set(self.IDs).union(self.aliases)

    @property
    def gases(self):
        '''[list] Gas components.'''
        return self[self.get_IDs_from_array(self.g)]

    @property
    def solids(self):
        '''[list] Solids (particulate) components.'''
        return self[self.get_IDs_from_array(self.x)]

    @property
    def inorganics(self):
        '''[list] Inorganic components.'''
        return self[self.get_IDs_from_array(self.inorg)]

    @property
    def inorganic_solids(self):
        '''[list] Inorganic solids (particulate & inorganic, all undegradable) components.'''
        return self[self.get_IDs_from_array(self.x*self.inorg)]

    @property
    def organic_solids(self):
        '''[list] Organic solids (particulate & organic) components.'''
        return self[self.get_IDs_from_array(self.x*self.org)]

    @property
    def substrates(self):
        '''
        [list] Substrate components.
        '''
        try: return self.__dict__['substrates']
        except:
            raise ValueError('The `substrates` group is not set, '
                             'use `define_group` to define the `substrates` group.')
    @substrates.setter
    def substrates(self, i):
        raise RuntimeError('Use `define_group` to define the `substrates` group.')

    @property
    def biomass(self):
        '''
        [list] Biomass components.
        '''
        try: return self.__dict__['biomass']
        except:
            raise ValueError('The `biomass` group is not set, '
                             'use `define_group` to define the `biomass` group.')
    @biomass.setter
    def biomass(self, i):
        raise RuntimeError('Use `define_group` to define the `biomass` group.')