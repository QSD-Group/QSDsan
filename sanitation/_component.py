#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 08:33:29 2020

TODO:
    Need to define measure unit


@author: yalinli_cabbi
"""
import os
import sys
import pandas as pd
import thermosteam as tmo

__all__ = ('Component',)


# %%

# =============================================================================
# Representation
# =============================================================================

# Component-related properties
_component_fields = ('charge', 'description', 'particle_size', 'degradability',
                     'organic',)

AbsoluteUnitsOfMeasure = tmo.units_of_measure.AbsoluteUnitsOfMeasure
component_units_of_measure = {
    'charge': AbsoluteUnitsOfMeasure('mol/mol')
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

    def __new__(cls, ID, formula=None, search_ID=None, phase='l', 
                charge=0, description=None, particle_size='Soluble', 
                degradability='Undegradable', organic=False,
                **chemical_properties):
        if search_ID:
            self = super().__new__(cls, ID=ID, search_ID=search_ID,
                                   search_db=True, **chemical_properties)
        else:
            self = super().__new__(cls, ID=ID, search_db=False, **chemical_properties)
            if formula:
                self._formula = formula
        self._ID = ID
        tmo._chemical.lock_phase(self, phase)
        self._charge = charge
        self._description = description
        assert particle_size in ('Dissolved gas', 'Soluble', 'Colloidal', 'Particulate'), \
            "particle_size must be 'Dissolved gas', 'Soluble', 'Colloidal', or 'Particulate'"
        self._particle_size = particle_size
        assert degradability in ('Biological', 'Chemical', 'Undegradable'), \
            "degradability must be 'Biological', 'Chemical', or 'Undegradable'"
        self._degradability = degradability
        assert organic.__class__ is bool, "organic must be 'True' or 'False'"
        self._organic = organic
        return self


    @property
    def charge(self):
        '''Charge content of the Component, [mol +/mol]'''
        return self._charge
    @charge.setter
    def charge(self, charge):
        self._charge = float(charge)

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
        self._particle_size = particle_size

    @property
    def degradability(self):
        '''[str] Biological (counted in BOD), Chemical (counted in COD), or Undegradable'''
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

    # TODO: in the future only show Component-related info, not chemical ones
    # If it's a Chemcial, then it should be defined as a Chemical rather than Component
    def show(self, chemical_specifications=True):
        '''Print all properties'''
        if chemical_specifications:
            super().show()
        info = str()
        header = '[Others] '
        section = []
        for field in _component_fields:
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


    @classmethod
    def load_components_from_excel(cls, path=None):
        '''
        Create Component objects based on properties defined in an Excel spreadsheet,
        return a Chemicals object that contains all created Component objects,
        note that the Chemicals object needs to be compiled before using in simulation'''
        if not path:
            path = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), 'components.xlsx')
            # path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'components.xlsx')
        
        data = pd.read_excel(path, sheet_name='components')
        components = tmo.Chemicals(())
        for i in range(data.shape[0]):
            component = Component.__new__(cls, ID=data['ID'][i])
            for j in {'ID', 'formula', 'phase', 'charge', 'description',
                      'particle_size', 'organic', 'degradability'}:
                field = '_' + j
                if pd.isna(data[j][i]): continue
                setattr(component, field, data[j][i])
            components.append(component)
        return components

















