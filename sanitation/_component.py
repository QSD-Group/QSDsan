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

_chemical_fields = tmo._chemical._chemical_fields
display_asfunctor = tmo._chemical.display_asfunctor
chemical_units_of_measure = tmo._chemical.chemical_units_of_measure

__all__ = ('Component',)


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
_component_fields = ('i_charge', 'description', 'particle_size', 'degradability',
                     'organic',)

AbsoluteUnitsOfMeasure = tmo.units_of_measure.AbsoluteUnitsOfMeasure
component_units_of_measure = {
    'i_charge': AbsoluteUnitsOfMeasure('mol/mol')
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
                i_charge=None, description=None, particle_size=None, 
                degradability=None, organic=None,
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
        self._i_charge = i_charge
        self._description = description
        assert particle_size in ('Dissolved gas', 'Soluble', 'Colloidal', 'Particulate'), \
            "particle_size must be 'Dissolved gas', 'Soluble', 'Colloidal', or 'Particulate'"
        self._particle_size = particle_size
        assert degradability in ('Biological', 'Chemical', 'Undegradable'), \
            "degradability must be 'Biological', 'Chemical', or 'Undegradable'"
        self._degradability = degradability
        assert organic.__class__ is bool, "organic must be True or False"
        self._organic = organic
        return self


    @property
    def i_charge(self):
        '''Charge content of the Component, [mol +/mol], negative values indicate anions'''
        return self._i_charge
    @i_charge.setter
    def i_charge(self, i_charge):
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


    def show(self, chemical_specifications=False):
        '''Print all properties'''
        info = component_identity(self, pretty=True)
        if chemical_specifications:
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
            for j in _component_fields:
                field = '_' + j
                if pd.isna(data[j][i]): continue
                setattr(component, field, data[j][i])
            components.append(component)
        return components

















