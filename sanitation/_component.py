#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Sanitation Explorer: Sustainable design of non-sewered sanitation technologies
Copyright (C) 2020, Sanitation Explorer Development Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    Joy Cheung

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-for-WaSH/sanitation/blob/master/LICENSE.txt
for license details.
'''


import thermosteam as tmo
from chemicals.elements import mass_fractions as get_mass_frac

__all__ = ('Component',)

_chemical_fields = tmo._chemical._chemical_fields
_checked_properties = tmo._chemical._checked_properties
display_asfunctor = tmo._chemical.display_asfunctor
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
_num_component_properties = ('i_C', 'i_N', 'i_P', 'i_K', 'i_Mg', 'i_Ca',
                             'i_mass', 'i_charge', 'i_COD',
                             'f_BOD5_COD', 'f_uBOD_COD', 'f_Vmass_Totmass', )

# Fields that cannot be left as None
_key_component_properties = ('particle_size', 'degradability', 'organic',
                             *_num_component_properties)

# All Component-related properties
_component_properties = ('measured_as', 'description', 
                         *_key_component_properties)
                         
_component_slots = (*tmo.Chemical.__slots__,
                    *tuple('_'+i for i in _component_properties))

_checked_properties = (*_checked_properties, *_key_component_properties)

AbsoluteUnitsOfMeasure = tmo.units_of_measure.AbsoluteUnitsOfMeasure
component_units_of_measure = {
    'i_C': AbsoluteUnitsOfMeasure('g'), 
    'i_N': AbsoluteUnitsOfMeasure('g'), 
    'i_P': AbsoluteUnitsOfMeasure('g'), 
    'i_K': AbsoluteUnitsOfMeasure('g'), 
    'i_Mg': AbsoluteUnitsOfMeasure('g'), 
    'i_Ca': AbsoluteUnitsOfMeasure('g'), 
    'i_mass': AbsoluteUnitsOfMeasure('g'), 
    'i_charge': AbsoluteUnitsOfMeasure('mol'),
    'i_COD': AbsoluteUnitsOfMeasure('g'),
    }


# %%

#!!! What should the gas/solid-phase Component value (e.g., CH4 in biogas)
# for particle_size and degradability? Can we put None or NA?
allowed_values = {
    'particle_size': ('Dissolved gas', 'Soluble', 'Colloidal', 'Particulate'),
    'degradability': ('Readily', 'Slowly', 'Undegradable'),
    'organic': (True, False),
    # 'measured_as': (None, 'COD', 'C', 'N', 'P'),
    }

def check_return_property(name, value):
    if name.startswith('i_') or name.startswith('f_'):
        if not value: return None
        try: return float(value)
        except: raise TypeError(f'{name} must be a number, not a {type(value).__name__}.')        
        if name.startswith('f_') and (value>1 or value<0):
            raise ValueError(f'{name} must be within [0,1].')
    elif name in allowed_values.keys():
        assert value in allowed_values[name], \
            f'{name} must be in {allowed_values[name]}.'
        return value

# =============================================================================
# Define the Component class
# =============================================================================

class Component(tmo.Chemical):
    '''A subclass of the Chemical object in the thermosteam package with additional attributes and methods for waste treatment'''         

    __slots__ = _component_slots

    def __new__(cls, ID='', search_ID=None, formula=None, phase=None, measured_as=None, 
                i_C=None, i_N=None, i_P=None, i_K=None, i_Mg=None, i_Ca=None,
                i_mass=None, i_charge=None, i_COD=None, 
                f_BOD5_COD=None, f_uBOD_COD=None, f_Vmass_Totmass=None,
                description=None, particle_size=None,
                degradability=None, organic=None, **chemical_properties):
        
        if search_ID:
            self = super().__new__(cls, ID=ID, search_ID=search_ID,
                                   search_db=True, **chemical_properties)
        else:
            self = super().__new__(cls, ID=ID, search_db=False, **chemical_properties)

        self._ID = ID
        if formula:
            self._formula = None
            self.formula = formula
        if phase: tmo._chemical.lock_phase(self, phase)
        self.measured_as = measured_as
        self.i_mass = i_mass
        self.i_C = i_C
        self.i_N = i_N
        self.i_P = i_P
        self.i_K = i_K
        self.i_Mg = i_Mg
        self.i_Ca = i_Ca
        self.i_charge = i_charge        
        self.f_BOD5_COD = f_BOD5_COD
        self.f_uBOD_COD = f_uBOD_COD
        self.f_Vmass_Totmass = f_Vmass_Totmass
        
        #!!! Need to check if this works
        #!!! This is necessary for correct molar flow calculation.
        if measured_as:
            if measured_as == 'COD': self._MW = tmo.Chemical('O').MW
            else:
                try: self._MW = tmo.Chemical(measured_as).MW
                except LookupError: raise LookupError(f"measure unit {measured_as} not recognized")

        self._particle_size = particle_size
        self._degradability = degradability
        self._organic = organic
        self.description = description
        if not self.MW and not self.formula: self.MW = 1.
        
        self.i_COD = i_COD
        
        return self

    #!!! use i_mass
    def _atom_frac_setter(self, atom=None, frac=None):
        if self.formula:
            if frac:
                raise AttributeError('This Component has formula, '
                                     f'i_{atom} is calculated based on formula, '
                                     'cannot be set.')
            else:
                if atom in self.atoms.keys():
                    # if self.measured_as:
                    #     if atom == self.measured_as:
                    #         return 1.0
                    #     elif atom =='O' and self.measured_as == 'COD':
                    #         return 1.0  #!!! this is not true.
                    #     else:
                    #         return 0.
                    try: return get_mass_frac(self.atoms)[atom] * self.i_mass
                    except: return None
                return 0. # does not have this atom
        else:
            return check_return_property(f'i_{atom}', frac)

    #!!! i_{} does not have to be within [0,1].
    @property
    def i_C(self):
        '''
        [float] Carbon content of the Component, [g C/g measure unit].

        Notes
        -------
        Will be calculated based on formula and measured_as if given.
        '''
        return self._i_C or 0.
    @i_C.setter
    def i_C(self, i):
        self._i_C = self._atom_frac_setter('C', i)


    @property
    def i_N(self):
        '''
        [float] Nitrogen content of the Component, [g N/g measure unit].

        Notes
        -------
        [1] If the Component is measured as N, then i_N is 1.
        [2] Will be calculated based on formula and measured_as if given.
        '''
        return self._i_N or 0.
    @i_N.setter
    def i_N(self, i):
        self._i_N = self._atom_frac_setter('N', i)

    @property
    def i_P(self):
        '''
        [float] Phosphorus content of the Component, [g P/g measure unit].

        Notes
        -------
        [1] If the Component is measured as P, then i_P is 1.
        [2] Will be calculated based on formula and measured_as if given.
        '''
        return self._i_P or 0.
    @i_P.setter
    def i_P(self, i):
        self._i_P = self._atom_frac_setter('P', i)

    @property
    def i_K(self):
        '''
        [float] Potassium content of the Component, [g K/g measure unit].

        Notes
        -------
        Will be calculated based on formula and measured_as if given.
        '''
        return self._i_K or 0.
    @i_K.setter
    def i_K(self, i):
        self._i_K = self._atom_frac_setter('K', i)

    @property
    def i_Mg(self):
        '''
        [float] Magnesium content of the Component, [g Mg/g measure unit].

        Notes
        -------
        Will be calculated based on formula and measured_as if given.
        '''
        return self._i_Mg or 0.
    @i_Mg.setter
    def i_Mg(self, i):
        self._i_Mg = self._atom_frac_setter('Mg', i)
        
    @property
    def i_Ca(self):
        '''
        [float] Calcium content of the Component, [g Ca/g measure unit].

        Notes
        -------
        Will be calculated based on formula and measured_as if given.
        '''
        return self._i_Ca or 0.
    @i_Ca.setter
    def i_Ca(self, i):
        self._i_Ca = self._atom_frac_setter('Ca', i)

    @property
    def i_mass(self):
        '''[float] Mass content of the Component, [g Component/g measure unit].'''
        if self.measured_as:
            if self.measured_as in self.atoms.keys():
                return 1/get_mass_frac(self.atoms)[self.measured_as]
            elif not self._i_mass:                 
                if self.measured_as == 'COD' and self.combustion:
                    chem_MW = tmo.Chemical(self.ID, search_ID=self.CAS).MW
                    i = -chem_MW/self.combustion['O2'] * tmo.Chemical('O2').MW
                    return check_return_property('i_mass', i)
                else: raise ValueError(f"Must specify i_mass for Component {self.ID} "
                                       f"measured as {self.measured_as}.")
        return self._i_mass or 1.
    @i_mass.setter
    def i_mass(self, i):
        if i and self.measured_as and self.formula:
            if self.measured_as != 'COD':
                raise AttributeError(f'This Component is measured as {self.measured_as} '
                                     f'and has formula, i_mass is calculated, '
                                     'cannot be set.')
        self._i_mass = check_return_property('i_mass', i)
        
    @property
    def i_charge(self):
        '''
        [float] Charge content of the Component, [mol +/g measure unit].

        Notes
        -------
        Positive values indicate cations and negative values indicate anions.
        '''
        return self._i_charge or 0.
    @i_charge.setter
    def i_charge(self, i):
        self._i_charge = check_return_property('i_charge', i)

    @property
    def f_BOD5_COD(self):
        '''
        BOD5 fraction in COD of the Component, unitless.

        Notes
        -------
        [1] Must be within [0,1].
        [2] Must be less than or equal to f_uBOD_COD.
        '''
        return self._f_BOD5_COD or 0.
    @f_BOD5_COD.setter
    def f_BOD5_COD(self, f):
        self._f_BOD5_COD = check_return_property('f_BOD5_COD', f)

    @property
    def f_uBOD_COD(self):
        '''
        [fload] Ultimate BOD fraction in COD of the Component, unitless.

        Notes
        -------
        [1] Must be within [0,1].
        [2] Must be greater than or equal to f_BOD5_COD.
        '''
        return self._f_uBOD_COD or 0.
    @f_uBOD_COD.setter
    def f_uBOD_COD(self, f):
        frac = f or 0.
        if frac < self.f_BOD5_COD:
            raise ValueError('f_uBOD_COD cannot be less than f_BOD5_COD.')
        self._f_uBOD_COD = check_return_property('f_uBOD_COD', frac)

    @property
    def f_Vmass_Totmass(self):
        '''
        [fload] Volatile fraction of the mass of the Component, unitless.

        Notes
        -------
        Must be within [0,1].
        '''
        return self._f_Vmass_Totmass or 0.
    @f_Vmass_Totmass.setter
    def f_Vmass_Totmass(self, f):
        self._f_Vmass_Totmass = check_return_property('f_Vmass_Totmass', f)
        
    @property
    def description(self):
        '''[str] Description of the Component.'''
        return self._description
    @description.setter
    def description(self, description):
        self._description = description

    @property
    def measured_as(self):
        '''
        [str] The unit as which the Component is measured.

        Notes
        -------
        Can be left as blank or chosen from 'COD', 'C', 'N', 'P' 
        or a chemical element of the Component.
        '''
        return self._measured_as
    @measured_as.setter
    def measured_as(self, measured_as):
        self._measured_as = measured_as
        
    @property
    def particle_size(self):
        '''
        [str] Size of the Component based on the type.

        Notes
        -------
        Must be chosen from 'Dissolved gas', 'Soluble', 'Colloidal', or 'Particulate'.
        '''
        return self._particle_size
    @particle_size.setter
    def particle_size(self, particle_size):
        self._particle_size = check_return_property('particle_size', particle_size)

    @property
    def degradability(self):
        '''
        [str] Degradability of the Component.

        Notes
        -------
        Must be chosen from 'Readily', 'Slowly', or 'Undegradable'.
        '''
        return self._degradability
    @degradability.setter
    def degradability(self, degradability):
        self._degradability = check_return_property('degradability', degradability)

    @property
    def organic(self):
        '''[bool] True (organic) or False (inorganic)'''
        return self._organic
    @organic.setter
    def organic(self, organic):
        self._organic = bool(check_return_property('organic', organic))

    @property
    def i_COD(self):
        '''[float] COD content, calculated based on measured_as, organic and formula.'''
        return self._i_COD

    @i_COD.setter
    def i_COD(self, i):
        if i: self._i_COD = check_return_property('i_COD', i)
        else: 
            if self.organic: 
                if self.measured_as == 'COD': self._i_COD = 1.
                elif self.measured_as == None: 
                    self._i_COD = -self.combustion['O2'] * tmo.Chemical('O2').MW/self.MW
                else:
                    chem_MW = tmo.Chemical(self.ID, search_ID=self.CAS).MW
                    self._i_COD = -self.combustion['O2'] * tmo.Chemical('O2').MW/chem_MW * self.i_mass 
            else: self._i_COD = 0.
       

    def show(self, chemical_info=False):
        '''
        Show Component properties.

        Parameters
        ----------
        chemical_info : [bool]
            Whether to show properties associated with the corresponding
            Chemical object of the Component. The default is False.
        '''
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
                        if field.startswith('i_'):
                            if self._measured_as: denom = self._measured_as
                            else: denom = ''
                            if field == 'i_charge': line += f' {units} +/g {denom}'
                            else: line += f' {units} {field[2:]}/g {denom}'
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

    def copy(self, ID, **data):
        new = self.__class__.__new__(cls=self.__class__, ID=ID)

        for field in self.__slots__:
            if field == '_CAS': continue
            value = getattr(self, field)
            setattr(new, field, copy_maybe(value))
        new._ID = ID
        new._locked_state = self._locked_state
        new._init_energies(new.Cn, new.Hvap, new.Psat, new.Hfus, new.Sfus, new.Tm,
                           new.Tb, new.eos, new.eos_1atm, new.phase_ref)
        new._label_handles()
        for i,j in data.items(): setattr(new, i , j)
        return new
    __copy__ = copy

    @classmethod
    def from_chemical(cls, ID, chemical, phase=None, measured_as=None, 
                      i_C=None, i_N=None, i_P=None, i_K=None, i_Mg=None, i_Ca=None,
                      i_mass=None, i_charge=None, i_COD=None,
                      f_BOD5_COD=None, f_uBOD_COD=None, f_Vmass_Totmass=None, 
                      description=None, particle_size=None, degradability=None, 
                      organic=None, **data):
        '''Return a new Component from Chemical'''
        new = cls.__new__(cls, ID=ID, phase=phase)
        for field in chemical.__slots__:
            value = getattr(chemical, field)
            setattr(new, field, copy_maybe(value))
        new._ID = ID
        if phase: new._locked_state = phase
        new._init_energies(new.Cn, new.Hvap, new.Psat, new.Hfus, new.Sfus, new.Tm,
                           new.Tb, new.eos, new.eos_1atm, new.phase_ref)
        new._label_handles()
        new.measured_as = measured_as        
        new.i_mass = i_mass
        new.i_C = i_C
        new.i_N = i_N
        new.i_P = i_P
        new.i_K = i_K
        new.i_Mg = i_Mg
        new.i_Ca = i_Ca
        new.i_charge = i_charge
        new.f_BOD5_COD = f_BOD5_COD
        new.f_uBOD_COD = f_uBOD_COD
        new.f_Vmass_Totmass = f_Vmass_Totmass
        new.description = description
        new.particle_size = particle_size
        new.degradability = degradability
        new.organic = organic
        new.i_COD = i_COD
        for i,j in data.items():
            if i == 'formula':
                new._formula = j
            else: setattr(new, i , j)
        return new


















