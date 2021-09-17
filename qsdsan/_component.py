#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Cheung <joycheung1994@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

Part of this module is based on the Thermosteam package:
https://github.com/BioSTEAMDevelopmentGroup/thermosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


import thermosteam as tmo
from chemicals.elements import (
    mass_fractions as get_mass_frac,
    molecular_weight,
    charge_from_formula
    )
from . import Chemical
from .utils import auom, copy_attr, cod_test_stoichiometry, electron_acceptor_cod

__all__ = ('Component',)

_chemical_fields = tmo._chemical._chemical_fields
_checked_properties = tmo._chemical._checked_properties
lock_phase = tmo._chemical.lock_phase
display_asfunctor = tmo._chemical.display_asfunctor
copy_maybe = tmo.utils.copy_maybe
get_component_data = tmo._chemical.get_chemical_data


# %%
# =============================================================================
# Filling missing properties
# =============================================================================

def unpickle_component(cmp_data):
    cmp = object.__new__(Component)
    setfield = setattr
    for field, value in cmp_data.items():
        setfield(cmp, field, value)
    return cmp

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
                             'i_mass', 'i_charge', 'i_COD', 'i_NOD',
                             'f_BOD5_COD', 'f_uBOD_COD', 'f_Vmass_Totmass', )

# Fields that cannot be left as None
_key_component_properties = ('particle_size', 'degradability', 'organic',
                             *_num_component_properties)

# All Component-related properties
_component_properties = ('measured_as', 'description',
                         *_key_component_properties)

_component_slots = (*Chemical.__slots__,
                    *tuple('_'+i for i in _component_properties))

_checked_properties = (*_checked_properties, *_key_component_properties)

component_units_of_measure = {
    'i_C': auom('g'),
    'i_N': auom('g'),
    'i_P': auom('g'),
    'i_K': auom('g'),
    'i_Mg': auom('g'),
    'i_Ca': auom('g'),
    'i_mass': auom('g'),
    'i_charge': auom('mol'),
    'i_COD': auom('g'),
    'i_NOD': auom('g'),
    }


# %%

#!!! What should the gas/solid-phase Component value (e.g., CH4 in biogas)
# for particle_size and degradability? Can we put None or NA?
allowed_values = {
    'particle_size': ('Dissolved gas', 'Soluble', 'Colloidal', 'Particulate'),
    'degradability': ('Readily', 'Slowly', 'Undegradable'),
    'organic': (True, False),
    }

def check_return_property(name, value):
    if name.startswith(('i_', 'f_')):
        try: return float(value)
        except:
            if not value: return None
            raise TypeError(f'{name} must be a number, not a {type(value).__name__}.')
        if name.startswith('f_') and (value>1 or value<0):
            raise ValueError(f'{name} must be within [0,1].')
    elif name in allowed_values.keys():
        assert value in allowed_values[name], \
            f'{name} must be in {allowed_values[name]}.'
        return value

# =============================================================================
# Define the Component class
# =============================================================================

class Component(Chemical):
    '''
    A subclass of :class:`thermosteam.Chemical` with additional attributes
    and methods for waste treatment.

    Parameters
    ----------
    ID : str
        ID for the component, must be unique.
    search_ID : str
        ID that will be passed to :class:`thermosteam.Chemical` to search the database.
    formula : str
        Formula for the component, formula from the database will be used
        if the component is constructed from the database and it has a formula in the database.
    phase : str
        If provided, this component will be assumed to only exist in the given phase.
    i_C : float
        Carbon content of the component, [g C/g measure unit].
    i_N : float
        Nitrogen content of the component, [g N/g measure unit].
    i_P : float
        Phosphorus content of the component, [g P/g measure unit].
    i_K : float
        Potassium content of the component, [g K/g measure unit].
    i_Mg : float
        Magnesium content of the component, [g Mg/g measure unit].
    i_Ca : float
        Calcium content of the component, [g Ca/g measure unit].
    i_mass : float
        Mass content of the component, [g Component/g measure unit].
    i_charge : float
        Charge content of the component, [mol +/g measure unit].
        Positive values indicate cations and negative values indicate anions.
    i_COD : float
        COD content, calculated based on `measured_as`, `organic`, and `formula` if not given.
    i_NOD : float
        Nitrogenous oxygen demand, calculated based on `measured_as`,
        `degradability`, and `formula` if not given.
    f_BOD5_COD : float
        BOD5 fraction in COD of the component, unitless.
    f_uBOD_COD : float
        Ultimate BOD fraction in COD of the component, unitless.
    f_Vmass_Totmass : float
        Volatile fraction of the mass of the component, unitless.
    description : str
        Description of the component.
    measured_as : str
        The unit as which the component is measured.
        Can be left as blank or chosen from 'COD', or a constituent
        element of the component.
    particle_size : str
        Size of the component based on the type.
        Must be chosen from 'Dissolved gas', 'Soluble', 'Colloidal', or 'Particulate'.
    degradability : str
        Degradability of the Component.
        Must be chosen from 'Readily', 'Slowly', or 'Undegradable'.
    organic : bool
        True (organic) or False (inorganic).
    chemical_properties : kwargs
        Will be passed to :class:`thermosteam.Chemical`.

    .. note::

        [1] Element ratios like `i_C`, `i_N`, `i_P`, `i_K`, `i_Mg`, and `i_Ca` will
        be calculated based on `formula` and `measured_as` if given; and the ratio
        will be 1 if the component is measured as this element.

        [2] For fractions including `f_BOD5_COD`, `f_uBOD_COD`, and `f_Vmass_Totmass`,
        their values must be within [0, 1].

        [3] If no formula or MW is provided, then MW of this component is assumed to
        1 to avoid ZeroDivisionError exception in calculation.


    Examples
    --------
    `Component <https://qsdsan.readthedocs.io/en/latest/tutorials/Component.html>`_


    See Also
    --------
    `thermosteam.Chemical <https://thermosteam.readthedocs.io/en/latest/Chemical.html>`_
    '''

    __slots__ = _component_slots

    # ID must be provided
    def __new__(cls, ID, search_ID=None, formula=None, phase=None, measured_as=None,
                i_C=None, i_N=None, i_P=None, i_K=None, i_Mg=None, i_Ca=None,
                i_mass=None, i_charge=None, i_COD=None, i_NOD=None,
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
        if phase: lock_phase(self, phase)

        self._measured_as = measured_as
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
        self._particle_size = particle_size
        self._degradability = degradability
        self._organic = organic
        self.description = description
        if not self.MW and not self.formula: self.MW = 1.
        self.i_COD = i_COD
        self.i_NOD = i_NOD
        return self


    def __reduce__(self):
        return unpickle_component, (get_component_data(self),)


    def _atom_frac_setter(self, atom=None, frac=None):
        if self.formula:
            if frac:
                raise AttributeError('This component has formula, '
                                     f'i_{atom} is calculated based on formula, '
                                     'cannot be set.')
            else:
                if atom in self.atoms.keys():
                    try: return get_mass_frac(self.atoms)[atom] * self.i_mass
                    except: return None
                return 0. # does not have this atom
        else:
            return check_return_property(f'i_{atom}', frac)


    @property
    def i_C(self):
        '''[float] Carbon content of the component, [g C/g measure unit].'''
        return self._i_C or 0.
    @i_C.setter
    def i_C(self, i):
        self._i_C = self._atom_frac_setter('C', i)

    @property
    def i_N(self):
        '''[float] Nitrogen content of the component, [g N/g measure unit].'''
        return self._i_N or 0.
    @i_N.setter
    def i_N(self, i):
        self._i_N = self._atom_frac_setter('N', i)

    @property
    def i_P(self):
        '''[float] Phosphorus content of the component, [g P/g measure unit].'''
        return self._i_P or 0.
    @i_P.setter
    def i_P(self, i):
        self._i_P = self._atom_frac_setter('P', i)

    @property
    def i_K(self):
        '''[float] Potassium content of the component, [g K/g measure unit].'''
        return self._i_K or 0.
    @i_K.setter
    def i_K(self, i):
        self._i_K = self._atom_frac_setter('K', i)

    @property
    def i_Mg(self):
        '''[float] Magnesium content of the component, [g Mg/g measure unit].'''
        return self._i_Mg or 0.
    @i_Mg.setter
    def i_Mg(self, i):
        self._i_Mg = self._atom_frac_setter('Mg', i)

    @property
    def i_Ca(self):
        '''[float] Calcium content of the component, [g Ca/g measure unit].'''
        return self._i_Ca or 0.
    @i_Ca.setter
    def i_Ca(self, i):
        self._i_Ca = self._atom_frac_setter('Ca', i)

    @property
    def i_mass(self):
        '''[float] Mass content of the component, [g Component/g measure unit].'''
        if self._i_mass is None: return 1
        else: return self._i_mass
    @i_mass.setter
    def i_mass(self, i):
        if self.atoms:
            if i: raise AttributeError(f'Component {self.ID} has formula, i_mass '
                                       f'is calculated, cannot be set.')
            else:
                if self.measured_as in self.atoms:
                    i = 1/get_mass_frac(self.atoms)[self.measured_as]
                elif self.measured_as == 'COD':
                    chem_MW = molecular_weight(self.atoms)
                    chem_charge = charge_from_formula(self.formula)
                    Cr2O7 = - cod_test_stoichiometry(self.atoms, chem_charge)['Cr2O7-2']
                    cod = Cr2O7 * 1.5 * molecular_weight({'O':2})
                    i = chem_MW/cod
                elif self.measured_as:
                    raise AttributeError(f'Must specify i_mass for component {self.ID} '
                                         f'measured as {self.measured_as}.')
        if self.measured_as == None:
            if i and i != 1:
                raise AttributeError(f'Component {self.ID} is measured as itself, '
                                     f'i_mass cannot be set to values other than 1.')
            i = 1
        self._i_mass = check_return_property('i_mass', i)

    #!!! Need to enable calculation from formula and water chemistry equilibria
    @property
    def i_charge(self):
        '''
        [float] Charge content of the component, [mol +/g measure unit].
        Positive values indicate cations and negative values indicate anions.
        '''
        return self._i_charge or 0.
    @i_charge.setter
    def i_charge(self, i):
        self._i_charge = check_return_property('i_charge', i)
        if not self._i_charge:
            if self.formula:
                charge = charge_from_formula(self.formula)
                chem_MW = molecular_weight(self.atoms)
                i = charge/chem_MW * self.i_mass
                self._i_charge = check_return_property('i_charge', i)
            else: self._i_charge = 0.

    @property
    def f_BOD5_COD(self):
        '''
        BOD5 fraction in COD of the component, unitless.
        Must be within [0, 1] and must be less than or equal to `f_uBOD_COD`.
        '''
        return self._f_BOD5_COD or 0.
    @f_BOD5_COD.setter
    def f_BOD5_COD(self, f):
        self._f_BOD5_COD = check_return_property('f_BOD5_COD', f)

    @property
    def f_uBOD_COD(self):
        '''
        [float] Ultimate BOD fraction in COD of the component, unitless.
        Must be within [0, 1] and must be larger than or equal to `f_BOD5_COD`.
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
        [float] Volatile fraction of the mass of the component, unitless.
        Must be within [0, 1].
        '''
        return self._f_Vmass_Totmass or 0.
    @f_Vmass_Totmass.setter
    def f_Vmass_Totmass(self, f):
        self._f_Vmass_Totmass = check_return_property('f_Vmass_Totmass', f)

    @property
    def description(self):
        '''[str] Description of the component.'''
        return self._description
    @description.setter
    def description(self, description):
        self._description = description

    @property
    def measured_as(self):
        '''
        [str] The unit as which the component is measured.
        Can be left as blank or chosen from 'COD', or a constituent
        element of the component.
        '''
        return self._measured_as
    @measured_as.setter
    def measured_as(self, measured_as):
        '''
        When measured_as is set to a different value, all i_{} values will
        be automatically updated.
        '''
        if measured_as:
            if measured_as == 'COD':
                self._MW = molecular_weight({'O':2})
            elif measured_as in self.atoms or 'i_'+measured_as in _num_component_properties:
                self._MW = molecular_weight({measured_as:1})
            else:
                raise AttributeError(f"Component {self.ID} must be measured as "
                                     f"either COD or one of its constituent atoms, "
                                     f"if not as itself.")

        if self._measured_as != measured_as:
            self._convert_i_attr(measured_as)

        self._measured_as = measured_as

    def _convert_i_attr(self, new):
        if new == None:
            denom = self._i_mass
        elif new == 'COD':
            denom = self._i_COD
        elif new in self.atoms or 'i_'+new in _num_component_properties:
            try: denom = getattr(self, '_i_'+new)
            except AttributeError:
                denom = get_mass_frac(self.atoms)[new] * self._i_mass
        else:
            raise AttributeError(f"Component {self.ID} must be measured as "
                                 f"either COD or one of its constituent atoms, "
                                 f"if not as itself.")

        if denom == 0:
            raise ValueError(f'{self.ID} cannot be measured as {new}')

        for field in _num_component_properties:
            if field.startswith('i_'):
                new_i = getattr(self, '_'+field)/denom
                setattr(self, '_'+field, new_i)

    @property
    def particle_size(self):
        '''
        [str] Size of the component based on the type.
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
        Must be chosen from 'Readily', 'Slowly', or 'Undegradable'.
        '''
        return self._degradability
    @degradability.setter
    def degradability(self, degradability):
        self._degradability = check_return_property('degradability', degradability)

    @property
    def organic(self):
        '''[bool] True (organic) or False (inorganic).'''
        return self._organic
    @organic.setter
    def organic(self, organic):
        self._organic = bool(check_return_property('organic', organic))


    @property
    def i_COD(self):
        '''[float] COD content, calculated based on `measured_as`, `organic`, and `formula` if not given.'''
        return self._i_COD or 0.
    @i_COD.setter
    def i_COD(self, i):
        if i != None: self._i_COD = check_return_property('i_COD', i)
        else:
            if self.organic or self.formula in ('H2', 'O2', 'N2', 'NO2-', 'NO3-'):
                if self.measured_as == 'COD': self._i_COD = 1.
                elif not self.atoms:
                    raise AttributeError(f"Must specify `i_COD` for organic component {self.ID}, "
                                         f"which is not measured as COD and has no formula.")
                else:
                    chem_MW = molecular_weight(self.atoms)
                    chem_charge = charge_from_formula(self.formula)
                    if self.formula in ('O2', 'N2', 'NO2-', 'NO3-'):
                        cod = electron_acceptor_cod(self.atoms, chem_charge) * molecular_weight({'O':2})
                    else:
                        Cr2O7 = - cod_test_stoichiometry(self.atoms, chem_charge)['Cr2O7-2']
                        cod = Cr2O7 * 1.5 * molecular_weight({'O':2})
                    self._i_COD = check_return_property('i_COD', cod/chem_MW * self.i_mass)
            else: self._i_COD = 0.

    @property
    def i_NOD(self):
        '''
        [float] Nitrogenous oxygen demand, calculated based on `measured_as`,
        `degradability`, and `formula` if not given.'''
        return self._i_NOD or 0.
    @i_NOD.setter
    def i_NOD(self, i):
        if i == None:
            if self.degradability in ('Readily', 'Slowly') or self.formula in ('H3N', 'NH4', 'NH3', 'NH4+'):
                i = self.i_N * molecular_weight({'O':4}) / molecular_weight({'N':1})
            elif self.formula in ('NO2-', 'HNO2'):
                i = self.i_N * molecular_weight({'O':1}) / molecular_weight({'N':1})
            else:
                i = 0.
        self._i_NOD = check_return_property('i_NOD', i)


    def __str__(self):
        return self._ID


    def __repr__(self):
        return f"Component('{self}')"


    def show(self, chemical_info=False):
        '''
        Show component properties.

        Parameters
        ----------
        chemical_info : bool
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
        '''Return a list of all missing thermodynamic properties.'''
        missing = []
        for i in (properties or _checked_properties):
            if getattr(self, i) == 0:
                continue
            elif str(getattr(self, i)) in ('True', 'False'):
                continue
            elif not getattr(self, i):
                missing.append(i)
        return missing


    def copy(self, new_ID, **data):
        '''
        Return a new :class:`Component` object with the same settings with
        alternative data set by kwargs.
        '''
        new = self.__class__.__new__(cls=self.__class__, ID=new_ID)
        new = copy_attr(new, self, skip=('_ID', '_CAS', '_N_solutes', '_locked_state'))
        new._ID = new_ID

        phase = data.get('phase') or self._locked_state
        new._locked_state = phase

        new._init_energies(new.Cn, new.Hvap, new.Psat, new.Hfus, new.Sfus,
                           new.Tm, new.Tb, new.eos, new.phase_ref)
        new._label_handles()

        for i,j in data.items():
            if i == 'phase':
                continue
            setattr(new, i , j)
        return new

    __copy__ = copy


    @classmethod
    def from_chemical(cls, ID, chemical=None, phase=None, measured_as=None,
                      i_C=None, i_N=None, i_P=None, i_K=None, i_Mg=None, i_Ca=None,
                      i_mass=None, i_charge=None, i_COD=None, i_NOD=None,
                      f_BOD5_COD=None, f_uBOD_COD=None, f_Vmass_Totmass=None,
                      description=None, particle_size=None, degradability=None,
                      organic=None, **data):
        '''Return a new :class:`Component` from a :class:`thermosteam.Chemical` object.'''
        new = cls.__new__(cls, ID=ID, phase=phase)

        if chemical is None:
            chemical = ID

        if isinstance(chemical, str):
            chemical = Chemical(chemical)

        for field in chemical.__slots__:
            value = getattr(chemical, field, None)
            setattr(new, field, copy_maybe(value))
        new._ID = ID

        if phase: new._locked_state = phase

        new._init_energies(new.Cn, new.Hvap, new.Psat, new.Hfus, new.Sfus,
                           new.Tm, new.Tb, new.eos, new.phase_ref)
        new._label_handles()
        new._measured_as = measured_as
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
        new._particle_size = particle_size
        new._degradability = degradability
        new._organic = organic
        new.i_COD = i_COD
        new.i_NOD = i_NOD
        for i,j in data.items():
            if i == 'formula':
                new._formula = j
            else: setattr(new, i , j)
        return new