#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 08:33:29 2020

TODO:
    Need to define measure unit


@author: yalinli_cabbi
"""

import os, sys
import pandas as pd
import thermosteam as tmo

__all__ = ('Component',)


# %%

# =============================================================================
# Representation
# =============================================================================

# Component-related properties
_component_fields = ('ID', 'formula', 'phase', 'i_C', 'i_N', 'i_P', 'i_K',
                     'i_mass', 'i_charge', 'f_BOD5_COD', 'f_BODu_COD', 'f_Vmass_Totmass',
                     'description', 'particle_size', 'degradability', 'organic',
                     'measure_unit')

#!!! Feel like these units should be for the stream (will know COD, etc.)
AbsoluteUnitsOfMeasure = tmo.units_of_measure.AbsoluteUnitsOfMeasure
component_units_of_measure = {
    # 'i_C': AbsoluteUnitsOfMeasure('kg C/(measure unit)'),
    # 'i_N': AbsoluteUnitsOfMeasure('kg N/(measure unit)'),
    # 'i_P': AbsoluteUnitsOfMeasure('kg P/(measure unit)'),
    # 'i_K': AbsoluteUnitsOfMeasure('kg K/(measure unit)'),
    # #!!! This is just the mass, might not need this
    # 'i_mass': AbsoluteUnitsOfMeasure('kg component/(measure unit)'),
    # 'i_charge': AbsoluteUnitsOfMeasure('mol +/(measure unit)')
    'i_C': AbsoluteUnitsOfMeasure('kg C/kg'),
    'i_N': AbsoluteUnitsOfMeasure('kg N/kg'),
    'i_P': AbsoluteUnitsOfMeasure('kg P/kg'),
    'i_K': AbsoluteUnitsOfMeasure('kg K/kg'),
    'i_charge': AbsoluteUnitsOfMeasure('mol/kg')
    }


# %%

# =============================================================================
# Define the Component class
# =============================================================================

# Component should not exist in the chemical database, otherwise should just use
# that chemical
class Component(tmo.Chemical):
    '''
    A subclass of the Chemical object in the thermosteam package with additional attributes and methods for waste treatment
    '''

    #!!! Some of the properties should be calculable
    def __new__(cls, ID, formula=None, phase='l', i_C=0, i_N=0, i_P=0, i_K=0, i_mass=0,
                i_charge=0, f_BOD5_COD=0, f_BODu_COD=0, f_Vmass_Totmass=0,
                description=None, particle_size='Soluble', degradability='Undegradable',
                organic=False, measure_unit=None):

        self = super().__new__(cls, ID=ID, search_db=False)
        self._ID = ID
        if formula:
            self._formula = formula
        tmo._chemical.lock_phase(self, phase)
        self._i_C = i_C
        self._i_N = i_N
        self._i_P = i_P
        self._i_K = i_K
        self._i_mass = i_mass
        self._i_charge = i_charge
        self._f_BOD5_COD = f_BOD5_COD
        self._f_BODu_COD = f_BODu_COD
        self._f_Vmass_Totmass = f_Vmass_Totmass
        self._description = description
        assert particle_size in ('Dissolved gas', 'Soluble', 'Colloidal', 'Particulate'), \
            "particle_size must be 'Dissolved gas', 'Soluble', 'Colloidal', or 'Particulate'"
        self._particle_size = particle_size
        assert degradability in ('Biological', 'Chemical', 'Undegradable'), \
            "degradability must be 'Biological', 'Chemical', or 'Undegradable'"
        self._degradability = degradability
        assert organic.__class__ is bool, "organic must be 'True' or 'False'"
        self._organic = organic
        self._measure_unit = measure_unit
        return self

    @property
    def i_C(self):
        '''Carbon content of the component, [g C/(measure unit)]'''
        return self._i_C
    @i_C.setter
    def i_C(self, i_C):
        self._i_C = float(i_C)
        
    @property
    def i_N(self):
        '''Nitrogen content of the component, [g N/(measure unit)]'''
        return self._i_N
    @i_N.setter
    def i_N(self, i_N):
        self._i_N = float(i_N)

    @property
    def i_P(self):
        '''Phosphorus content of the component, [g P/(measure unit)]'''
        return self._i_P
    @i_P.setter
    def i_P(self, i_P):
        self._i_P = float(i_P)
    
    @property
    def i_K(self):
        '''Potassium content of the component, [g K/(measure unit)]'''
        return self._i_K
    @i_K.setter
    def i_K(self, i_K):
        self._i_K = float(i_K)

    @property
    def i_mass(self):
        '''Mass content of the component, [g component/(measure unit)]'''
        return self._i_mass
    @i_mass.setter
    def i_mass(self, i_mass):
        self._i_mass = float(i_mass)

    @property
    def i_charge(self):
        '''Charge content of the component, [mole+/(measure unit)]'''
        return self._i_charge
    @i_charge.setter
    def i_charge(self, i_charge):
        self._i_charge = float(i_charge)

    @property
    def f_BOD5_COD(self):
        '''Five-day biochemcial oxygen demand to chemical oxygen demand ratio [g BOD5/g COD]'''
        return self._f_BOD5_COD
    @f_BOD5_COD.setter
    def f_BOD5_COD(self, f_BOD5_COD):
        self._f_BOD5_COD = float(f_BOD5_COD)

    @property
    def f_BODu_COD(self):
        '''Ultimate biochemcial oxygen demand to chemical oxygen demand ratio [g BODu/g COD]'''
        return self._f_BODu_COD
    @f_BODu_COD.setter
    def f_BODu_COD(self, f_BODu_COD):
        self._f_BODu_COD = float(f_BODu_COD)

    @property
    def f_Vmass_Totmass(self):
        '''Volatile mass to total mass ratio'''
        return self._f_Vmass_Totmass
    @f_Vmass_Totmass.setter
    def f_Vmass_Totmass(self, f_Vmass_Totmass):
        self._f_Vmass_Totmass = float(f_Vmass_Totmass)

    @property
    def description(self):
        '''[str] Description of the component'''
        return self._description
    @description.setter
    def description(self, description):
        self._description = description

    @property
    def particle_size(self):
        '''[str] Dissolved gas, soluble, colloidal, or particulate'''
        return self._particle_size
    @particle_size.setter
    def particle_size(self, particle_size):
        self._particle_size = particle_size

    @property
    def degradability(self):
        '''[str] Biologicallly (counted in BOD), chemically (counted in COD), or undegradable'''
        return self._degradability
    @degradability.setter
    def degradability(self, degradability):
        self._degradability = degradability

    @property
    def organic(self):
        '''[bool] True (organic) or False (inorganic)'''
        return self._organic
    @organic.setter
    def organic(self, organic):
        self._organic = organic

    @property
    def measure_unit(self):
        '''[bool] True (organic) or False (inorganic)'''
        return self._measure_unit
    @measure_unit.setter
    def measure_unit(self, measure_unit):
        self._measure_unit = measure_unit

    def show(self, chemical_specifications=True):
        '''Print all properties'''
        if chemical_specifications:
            super().show()
        info = str()
        header = '[Others] '
        section = []
        for field in _component_fields:
            if field in ('ID', 'formula', 'phase'): continue
            value = getattr(self, field)
            field = field.lstrip('_')
            if value is None:
                line = f"{field}: None"
            else:
                if isinstance(value, (int, float)):
                    line = f"{field}: {value:.5g}"
                    # units = chemical_units_of_measure.get(field, "")
                    # if units: line += f' {units}'
                else:
                    value = str(value)
                    line = f"{field}: {value}"
                    if len(line) > 40: line = line[:40] + '...'
            section.append(line)
        info += header + ("\n" + 9*" ").join(section)
        print(info)
        
    _ipython_display_ = show


    @classmethod
    def from_dict(cls, ID, property_dict):
        '''Create and return a Component object from dictionary'''
        self = Component.__new__(cls, ID)
        for i, j in property_dict.items():
            setattr(self, i, j)
        return self

    @classmethod
    def load_default_components(cls):
        '''Create and return a Chemicals object that contains all default components,
            note that the Chemicals object needs to be compiled before using in simulation'''
        path = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), 'default_components.xlsx')
        # path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'default_components.xlsx')
        
        data = pd.read_excel(path, sheet_name='components')
        components = tmo.Chemicals(())
        for i in range(data.shape[0]):
            component = Component.__new__(cls, ID=data['ID'][i])
            for j in _component_fields:
                field = '_' + j
                if pd.isna(data[j][i]): continue
                setattr(component, field, data[j][i])
            components.append(component)
        return components




















