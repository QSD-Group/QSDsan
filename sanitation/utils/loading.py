#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 08:33:29 2020

@author: yalinli_cabbi
"""

import os
import sys
import pandas as pd
import sanitation

__all__ = ('load_components_from_excel',)
_component_properties = sanitation._component._component_properties


# %%

def load_components_from_excel(path):
    '''
    Create Component objects based on properties defined in an Excel spreadsheet,
    return a Components object that contains all created Component objects,
    note that the Components object needs to be compiled before using in simulation'''
    # if not path:
    #     path = os.path.join(os.path.abspath(os.path.dirname(sys.argv[0])), 'default_components.xlsx')
    #     # path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'default_components.xlsx')
    
    data = pd.read_excel(path, sheet_name='components')
    components = sanitation.Components(())
    for i in range(data.shape[0]):
        component = sanitation.Component(ID=data['ID'][i])
        for j in _component_properties:
            field = '_' + j
            if pd.isna(data[j][i]): continue
            setattr(component, field, data[j][i])
        components.append(component)
    return components











