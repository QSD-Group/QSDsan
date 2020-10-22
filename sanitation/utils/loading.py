#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 08:33:29 2020

@author: yalinli_cabbi
"""

import pandas as pd
import thermosteam as tmo
import sanitation
from sanitation._component import _component_properties

__all__ = ('load_components_from_excel', 'load_default_components')


# %%

def load_components_from_excel(path='default'):    
    '''
    Create Component objects based on properties defined in an Excel spreadsheet,
    return a Components object that contains all created Component objects,
    note that the Components object needs to be compiled before using in simulation

    Parameters
    ----------
    path : [str], optional
        Path of the Excel file, when set to 'default' (the default value),
        will return default components.

    Returns
    -------
    A Components object that contains all created Component objects.
        
    Notes
    -------
    The Components object needs to be compiled before using in simulation

    '''

    import os
    # import sys
    if path == 'default':
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'default_components.xlsx')
        # path = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), 'default_components.xlsx')
    data = pd.read_excel(path, sheet_name='components')
    components = sanitation.Components(())
    for i in range(data.shape[0]):
        component = sanitation.Component(ID=data['ID'][i])
        for j in _component_properties:
            field = '_' + j
            if pd.isna(data[j][i]): continue
            setattr(component, field, data[j][i])
        components.append(component)
    del data, os
    # del sys
    return components

def load_default_components(compile_components=True):
    
    '''
    Create default Component objects.

    Parameters
    ----------
    compile_components : [bool], optional
        Whether to compile the default. The default is True.

    Returns
    -------
    A Components or CompiledComponents object with default Component objects.

    Notes
    -------
    [1] Component-specific properties are defined in the default_component.xlsx Excel file.
    [2] When compile_components is True, Chemical-specific properties are defaulted to those of water.

    '''
    assert compile_components.__class__ is bool, \
        f'compile_components must be True or False, not {compile_components}'
    components = load_components_from_excel()
    H2O = sanitation.Component.from_chemical('H2O', tmo.Chemical('H2O'),
                                  i_C=0, i_N=0, i_P=0, i_K=0, i_mass=1,
                                  i_charge=0, f_BOD5_COD=0, f_uBOD_COD=0,
                                  f_Vmass_Totmass=0,
                                  particle_size='Soluble',
                                  degradability='Undegradable', organic=False)
    components.append(H2O)
    if compile_components:    
        for i in components:
            i.default()
            i.copy_models_from(components.H2O, ['sigma', 'epsilon', 'kappa', 'V', 'Cn', 'mu'])
        components.compile()
    return components








