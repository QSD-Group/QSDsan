#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Oct 19 21:12:51 2020

@author: yalinli_cabbi
"""

import thermosteam as tmo
from thermosteam import Chemical, Chemicals, CompiledChemicals
from . import _component

__all__ = ('Components', 'CompiledComponents')

utils = tmo.utils
Component = _component.Component
_key_component_properties = _component._key_component_properties
setattr = object.__setattr__


# %%

class UndefinedComponent(AttributeError):
    '''AttributeError regarding undefined Component objects'''
    def __init__(self, ID): super(Chemical, self).__init__(repr(ID))

# =============================================================================
# Define the Components class
# =============================================================================

class Components(Chemicals):
    '''
    A subclass of the Chemicals object in the thermosteam package, contains Component objects as attributes
    '''
    
    def __new__(cls, components, cache=False):
        self = super(Chemicals, cls).__new__(cls)
        for component in components:
            if isinstance(component, Component):
                setattr(self, component.ID, component)
            else:
                if isinstance(component, Chemical):
                    raise TypeError(f'{component} is a Chemical object, use <Component>.from_chemical to define a Component object')
                raise TypeError(f'Only Component objects can be included, not a {type(component).__name__} object')
        return self
    
    def __setattr__(self, ID, component):
        raise TypeError("Cannot set attribute; use <Components>.append instead")
    
    def __setitem__(self, ID, component):
        raise TypeError("Cannot set item; use <Components>.append instead")
    
    def __getitem__(self, key):
        '''Return a Component object or a list of Component objects'''
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
    
    def append(self, component):
        '''Append a Component'''
        if not isinstance(component, Component):
            if isinstance(component, Chemical):
                raise TypeError(f'{component} is a Chemical object, use <Component>.from_chemical to define a Component object')
            else:
                raise TypeError("only 'Component' objects can be appended, "
                               f"not '{type(component).__name__}'")
        ID = component.ID
        if ID in self.__dict__:
            raise ValueError(f"{ID} already defined in Components")
        setattr(self, ID, component)
    
    def extend(self, components):
        '''Extend with more Component objects'''
        if isinstance(components, Components):
            self.__dict__.update(components.__dict__)
        else:
            for component in components: self.append(component)
    
    def compile(self):
        '''Cast as a CompiledComponents object'''
        CompiledComponents._compile(self)
        setattr(self, '__class__', CompiledComponents)
        


# %%

# =============================================================================
# Define the CompiledComponents class
# =============================================================================

chemical_data_array = tmo._chemicals.chemical_data_array

def component_data_array(components, attr):
    data = chemical_data_array(components, attr)
    return data

class CompiledComponents(CompiledChemicals):
    ''' A subclass of the CompiledChemicals object in the thermosteam package, contains Component objects as attributes'''
    
    def refresh_constants(self):
        '''Refresh constant arrays of Components, including all Chemical and Component-specific properties'''
        super().refresh_constants()
        dct = self.__dict__
        components = self.tuple
        dct['charge'] = chemical_data_array([i.charge for i in components])
        # TODO: add C, N, contents, etc.

    def _compile(self):
        dct = self.__dict__
        tuple_ = tuple # to differentiate it from dct['tuple']?
        components = tuple_(dct.values())
        CompiledChemicals._compile(self)
        
        for component in components:
            missing_properties = component.get_missing_properties(_key_component_properties)
            if not missing_properties: continue
            if 'i_charge' in missing_properties:
                if getattr(component, 'i_charge') == 0:
                    continue
            if 'organic' in missing_properties:
                if str(getattr(component, 'organic')) in ('True', 'False'):
                    continue
            missing = utils.repr_listed_values(missing_properties)
            raise RuntimeError(f'{component} is missing key Component-related properties ({missing})')
        
        # TODO: add other ones
        dct['i_charge'] = component_data_array(components, 'i_charge')

    
    def subgroup(self, IDs):
        '''Create a new subgroup of Component objects'''
        components = self[IDs]
        new = Components(components)
        new.compile()
        for i in new.IDs:
            for j in self.get_synonyms(i):
                try: new.set_synonym(i, j)
                except: pass
        return new
    
    def index(self, ID):
        '''Return index of specified Component'''
        try: return self._index[ID]
        except KeyError:
            raise UndefinedComponent(ID)

    def indices(self, IDs):
        try:
            dct = self._index
            return [dct[i] for i in IDs]
        except KeyError as key_error:
            raise UndefinedComponent(key_error.args[0])
    
    def __contains__(self, component):
        if isinstance(component, str):
            return component in self.__dict__
        elif isinstance(component, Component):
            return component in self.tuple
        else: # pragma: no cover
            return False
    
    
    
    
    
    
    
    