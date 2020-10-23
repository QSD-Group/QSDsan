#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 08:33:29 2020

@author: yalinli_cabbi
"""
import os
import pandas as pd
import sanitation
from sanitation._component import _component_properties

__all__ = ('load_components',)


# %%

def load_components(path=None):
    '''
    Create Component objects based on properties defined in a .csv file,
    return a Components object that contains all created Component objects,
    note that the Components object needs to be compiled before using in simulation'''
    # import sys
    if not path:
        path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'default_components.csv')
    data = pd.read_csv(path)
    components = sanitation.Components(())
    for i, cmp in data.iterrows():
        if pd.isna(cmp.measured_as):
            cmp.measured_as = None
        try:
            component = sanitation.Component(ID = cmp.ID, 
                                             search_ID = str(cmp.CAS), 
                                             measured_as = cmp.measured_as)
        except LookupError:
            try:
                component = sanitation.Component(ID = cmp.ID, 
                                                 search_ID = 'PubChem='+str(int(cmp.PubChem)), 
                                                 measured_as = cmp.measured_as) 
            except:
                if not pd.isna(cmp.formula):
                    component = sanitation.Component(ID = cmp.ID, 
                                                     formula = cmp.formula, 
                                                     measured_as = cmp.measured_as)            
                else:
                    component = sanitation.Component(ID = cmp.ID, 
                                                     measured_as = cmp.measured_as)
        for j in _component_properties:
            field = '_' + j
            if pd.isna(cmp[j]): continue
            setattr(component, field, cmp[j])
        components.append(component)
    
    del data
    return components
    # del sys, os
