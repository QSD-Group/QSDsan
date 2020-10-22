#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Oct  8 08:33:29 2020
Modified on Wed Oct 21 2020

TODO:
    Need to define measure unit


@author: yalinli_cabbi, joy_c
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


# Will stored as an array when compiled
_num_component_properties = ('i_C', 'i_N', 'i_P', 'i_K', 'i_mass', 'i_charge',
                             'f_BOD5_COD', 'f_uBOD_COD', 'f_Vmass_Totmass', )

# Fields that cannot be left as None
_key_component_properties = (*_num_component_properties,
                             'particle_size', 'degradability', 'organic', 'measure_unit')

# All Component-related properties
_component_properties = (*_key_component_properties,
                         'description', )

_checked_properties = (*_checked_properties, *_key_component_properties)

AbsoluteUnitsOfMeasure = tmo.units_of_measure.AbsoluteUnitsOfMeasure
component_units_of_measure = {
    'i_C': AbsoluteUnitsOfMeasure('g C'), 
    'i_N': AbsoluteUnitsOfMeasure('g N'), 
    'i_P': AbsoluteUnitsOfMeasure('g P'), 
    'i_K': AbsoluteUnitsOfMeasure('g K'), 
    'i_mass': AbsoluteUnitsOfMeasure('g'), 
    'i_charge': AbsoluteUnitsOfMeasure('mol'),
    # 'f_BOD5_COD': AbsoluteUnitsOfMeasure(''), 
    # 'f_uBOD_COD': AbsoluteUnitsOfMeasure(''), 
    # 'f_Vmass_Totmass': AbsoluteUnitsOfMeasure('')
    }


# %%

allowed_values = {
    'particle_size': ('Dissolved gas', 'Soluble', 'Colloidal', 'Particulate'),
    'degradability': ('Biological', 'Chemical', 'Undegradable'),
    'organic': (True, False)
    }

def check_property(name, value):
    if name in ('i_C', 'i_N', 'i_P', 'i_K', 'i_mass', 'i_charge',):
        try: float(value)
        except: raise TypeError(f'{name} must be a number, not a {type(value).__name__}')
    elif name in ('f_BOD5_COD', 'f_uBOD_COD', 'f_Vmass_Totmass',):
        try: float(value)
        except: raise TypeError(f'{name} must be a number, not a {type(value).__name__}')        
        if value>1 or value<0:
            raise ValueError(f'{name} must be within [0,1].')
    elif name in allowed_values.keys():
        assert value in allowed_values[name], \
            f'{name} must be in {allowed_values[name]}'

# =============================================================================
# Define the Component class
# =============================================================================

class Component(tmo.Chemical):
    '''A subclass of the Chemical object in the thermosteam package with additional attributes and methods for waste treatment'''

    def __new__(cls, ID, formula=None, search_ID=None, phase='l', measure_unit='g COD', 
                i_C=None, i_N=None, i_P=None, i_K=None, i_mass=None, i_charge=None, 
                f_BOD5_COD=None, f_uBOD_COD=None, f_Vmass_Totmass=None,
                description=None, particle_size=None, 
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
        self._i_C = i_C
        self._i_N = i_N
        self._i_P = i_P
        self._i_K = i_K
        self._i_mass = i_mass
        self._i_charge = i_charge
        
        self._f_BOD5_COD = f_BOD5_COD
        self._f_uBOD_COD = f_uBOD_COD
        self._f_Vmass_Totmass = f_Vmass_Totmass
        
        
        self._particle_size = particle_size
        self._degradability = degradability
        self._organic = organic
        
        # self._measure_unit = AbsoluteUnitsOfMeasure(measure_unit)
        self._measure_unit = measure_unit
        self._description = description
        return self

    @property
    def i_C(self):
        '''Carbon content of the Component, [g C/(measure unit)].'''
        return self._i_C
    @i_C.setter
    def i_C(self, i):
        check_property('i_C', i)
        self._i_C = float(i)

    @property
    def i_N(self):
        '''Nitrogen content of the Component, [g N/(measure unit)].'''
        return self._i_N
    @i_N.setter
    def i_N(self, i):
        check_property('i_N', i)
        self._i_N = float(i)

    @property
    def i_P(self):
        '''Phosphorus content of the Component, [g P/(measure unit)].'''
        return self._i_P
    @i_P.setter
    def i_P(self, i):
        check_property('i_P', i)
        self._i_P = float(i)

    @property
    def i_K(self):
        '''Potassium content of the Component, [g K/(measure unit)].'''
        return self._i_K
    @i_K.setter
    def i_K(self, i):
        check_property('i_K', i)
        self._i_K = float(i)

    @property
    def i_mass(self):
        '''Mass content of the Component, [g component/(measure unit)].'''
        return self._i_mass
    @i_mass.setter
    def i_mass(self, i):
        check_property('i_mass', i)
        self._i_mass = float(i)
        
    @property
    def i_charge(self):
        '''Charge content of the Component, [mol +/(measure unit)], negative values indicate anions'''
        return self._i_charge
    @i_charge.setter
    def i_charge(self, i):
        check_property('i_charge', i)
        self._i_charge = float(i)

    @property
    def f_BOD5_COD(self):
        '''BOD5 fraction in COD of the Component, unitless.'''
        return self._f_BOD5_COD
    @f_BOD5_COD.setter
    def f_BOD5_COD(self, f):
        check_property('f_BOD5_COD', f)
        self._f_BOD5_COD = float(f)

    @property
    def f_uBOD_COD(self):
        '''ultimate BOD fraction in COD of the Component, unitless.'''
        return self._f_uBOD_COD
    @f_uBOD_COD.setter
    def f_uBOD_COD(self, f):
        check_property('f_uBOD_COD', f)
        if f < self.f_BOD5_COD:
            raise ValueError('f_uBOD_COD cannot be less than f_BOD5_COD')
        self._f_uBOD_COD = float(f)

    @property
    def f_Vmass_Totmass(self):
        '''Volatile fraction of the mass of the Component, unitless.'''
        return self._f_Vmass_Totmass
    @f_Vmass_Totmass.setter
    def f_Vmass_Totmass(self, f):
        check_property('f_Vmass_Totmass', f)
        self._f_Vmass_Totmass = float(f)
        
    @property
    def description(self):
        '''[str] Description of the Component'''
        return self._description
    @description.setter
    def description(self, description):
        self._description = description

    @property
    def measure_unit(self):
        '''[str] The measuring unit of the Component'''
        return self._measure_unit
    @measure_unit.setter
    def measure_unit(self, measure_unit):
        self._measure_unit = measure_unit
        # self._measure_unit = AbsoluteUnitsOfMeasure(measure_unit)
        
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
        info = ''
        if chemical_info:
            super().show()
        else:
            info = component_identity(self, pretty=True)
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
                    if units: 
                        if field in ('i_C', 'i_N', 'i_P', 'i_K', 'i_mass', 'i_charge',):
                            line += f' {units}/{self._measure_unit}'
                        else: line += f' {units}'
                else:
                    value = str(value)
                    line = f"{field}: {value}"
                    if len(line) > 40: line = line[:40] + '...'
            section.append(line)
        info += header + ("\n" + 9*" ").join(section)
        print(info)
        
    _ipython_display_ = show

    def get_missing_properties(self, properties=None):
        missing = []
        for i in (properties or _checked_properties):
            if getattr(self, i) == 0:
                continue
            elif str(getattr(self, i)) in ('True', 'False'):
                continue
            elif not getattr(self, i):
                missing.append(i)
        return missing

    @classmethod
    def from_chemical(cls, ID, chemical, measure_unit='g', 
                      i_C=None, i_N=None, i_P=None, i_K=None, i_mass=None, 
                      i_charge=None, f_BOD5_COD=None, f_uBOD_COD=None, 
                      f_Vmass_Totmass=None, description=None, 
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
        
        new.i_C = i_C
        new.i_N = i_N
        new.i_P = i_P
        new.i_K = i_K
        new.i_mass = i_mass
        new.i_charge = i_charge
        new.f_BOD5_COD = f_BOD5_COD
        new.f_uBOD_COD = f_uBOD_COD
        new.f_Vmass_Totmass = f_Vmass_Totmass
        new.description = description
        new.particle_size = particle_size
        new.degradability = degradability
        new.organic = organic
        new.measure_unit = measure_unit
        
        for i,j in data.items(): setattr(new, i , j)
        return new


















