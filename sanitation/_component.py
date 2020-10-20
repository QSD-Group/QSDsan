#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Oct  8 08:33:29 2020

TODO:
    Need to define measure unit


@author: yalinli_cabbi
"""

import thermosteam as tmo

__all__ = ('Component',)

_chemical_fields = tmo._chemical._chemical_fields
_checked_properties = tmo._chemical._checked_properties
display_asfunctor = tmo._chemical.display_asfunctor
chemical_units_of_measure = tmo._chemical.chemical_units_of_measure
copy_maybe = tmo.utils.copy_maybe


# %%

# =============================================================================
# Representation
# =============================================================================

def component_identity(component, pretty=False):
    typeheader = f"{type(component).__name__}:"
    full_ID = f"{typeheader} {component.ID} (phase_ref={repr(component.phase_ref)})"
    phase = component.locked_state
    state = ' at ' + f"phase={repr(phase)}" if phase else ""
    return full_ID + state


# Component-related properties
_component_properties = ('i_charge', 'description', 'particle_size',
                         'degradability', 'organic',)

# Fields that cannot be left as None
_key_component_properties = ('i_charge', 'particle_size', 'degradability', 'organic',)

_checked_properties = (*_checked_properties, *_key_component_properties)

AbsoluteUnitsOfMeasure = tmo.units_of_measure.AbsoluteUnitsOfMeasure
component_units_of_measure = {
    'i_charge': AbsoluteUnitsOfMeasure('mol/mol')
    }


# %%

allowed_values = {
    'particle_size': ('Dissolved gas', 'Soluble', 'Colloidal', 'Particulate'),
    'degradability': ('Biological', 'Chemical', 'Undegradable'),
    'organic': (True, False)
    }

def check_property(name, value):
    if name in ('i_charge', ):
        try: float(value)
        except: raise TypeError(f'{name} must be a number, not a {type(value).__name__}')
    elif name in allowed_values.values():
        assert value in allowed_values[name], \
            f'{name} must be in {allowed_values[name]}'

# =============================================================================
# Define the Component class
# =============================================================================

# Component should not exist in the chemical database, otherwise should just use
# that chemical
class Component(tmo.Chemical):
    '''A subclass of the Chemical object in the thermosteam package with additional attributes and methods for waste treatment'''

    def __new__(cls, ID, formula=None, search_ID=None, phase='l', 
                i_charge=None, description=None, particle_size=None, 
                degradability=None, organic=None, **chemical_properties):
        if search_ID:
            self = super().__new__(cls, ID=ID, search_ID=search_ID,
                                   search_db=True, **chemical_properties)
        else:
            self = super().__new__(cls, ID=ID, search_db=False, **chemical_properties)
            if formula:
                self._formula = formula
        self._ID = ID
        tmo._chemical.lock_phase(self, phase)
        
        self._i_charge = i_charge
        self._description = description
        self._particle_size = particle_size
        self._degradability = degradability
        self._organic = organic
        return self


    @property
    def i_charge(self):
        '''Charge content of the Component, [mol +/mol], negative values indicate anions'''
        return self._i_charge
    @i_charge.setter
    def i_charge(self, i_charge):
        check_property('i_charge', i_charge)
        self._i_charge = float(i_charge)

    @property
    def description(self):
        '''[str] Description of the Component'''
        return self._description
    @description.setter
    def description(self, description):
        self._description = description

    @property
    def particle_size(self):
        '''[str] Dissolved gas, Soluble, Colloidal, or Particulate'''
        return self._particle_size
    @particle_size.setter
    def particle_size(self, particle_size):
        check_property('particle_size', particle_size)
        self._particle_size = particle_size

    @property
    def degradability(self):
        '''[str] Biological (counted in BOD), Chemical (counted in COD), or Undegradable'''
        return self._degradability
    @degradability.setter
    def degradability(self, degradability):
        check_property('degradability', degradability)
        self._degradability = degradability

    @property
    def organic(self):
        '''[bool] True (organic) or False (inorganic)'''
        return self._organic
    @organic.setter
    def organic(self, organic):
        check_property('organic', organic)
        self._organic = organic


    def show(self, chemical_info=False):
        '''Show Component properties'''
        info = component_identity(self, pretty=True)
        if chemical_info:
            for header, fields in _chemical_fields.items():
                section = []
                for field in fields:
                    value = getattr(self, field)
                    field = field.lstrip('_')
                    if value is None:
                        line = f"{field}: None"
                    if callable(value):
                        line = f"{display_asfunctor(value, name=field, var=field, show_var=False)}"
                    else:
                        if isinstance(value, (int, float)):
                            line = f"{field}: {value:.5g}"
                            units = chemical_units_of_measure.get(field, "")
                            if units: line += f' {units}'
                        else:
                            value = str(value)
                            line = f"{field}: {value}"
                            if len(line) > 27: line = line[:27] + '...'
                    section.append(line)
                if section:
                    info += header + ("\n" + 9*" ").join(section)
        info += '\nComponent-specific properties:\n'
        header = '[Others] '
        section = []
        for field in _component_properties:
            value = getattr(self, field)
            field = field.lstrip('_')
            if value is None:
                line = f"{field}: None"
            elif str(value) in ('True', 'False'):
                line = f"{field}: {value}" 
            else:
                if isinstance(value, (int, float)):
                    line = f"{field}: {value:.5g}"
                    units = component_units_of_measure.get(field, '')
                    if units: line += f' {units}'
                else:
                    value = str(value)
                    line = f"{field}: {value}"
                    if len(line) > 40: line = line[:40] + '...'
            section.append(line)
        info += header + ("\n" + 9*" ").join(section)
        print(info)
        
    _ipython_display_ = show

    def get_missing_properties(self, properties=None):
        return [i for i in (properties or _checked_properties) if not getattr(self, i)]

    @classmethod
    def from_chemical(cls, ID, chemical, i_charge=None, description=None, 
                      particle_size=None, degradability=None, organic=None,
                      **data):
        '''Make a new Component from Chemical'''
        new = cls.__new__(cls, ID=ID)
        for field in chemical.__slots__: 
            value = getattr(chemical, field)
            setattr(new, field, copy_maybe(value))
        new._ID = ID
        new._locked_state = new._locked_state
        new._init_energies(new.Cn, new.Hvap, new.Psat, new.Hfus, new.Tm,
                           new.Tb, new.eos, new.eos_1atm, new.phase_ref)
        new._label_handles()
        
        # TODO: add other properties
        new.i_charge = i_charge
        new.description = description
        new.particle_size = particle_size
        new.degradability = degradability
        new.organic = organic
        for i,j in data.items(): setattr(new, i , j)
        return new


















