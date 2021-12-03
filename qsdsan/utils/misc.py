#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

Part of the code is based on modules in BioSTEAM,
https://github.com/BioSTEAMDevelopmentGroup/biosteam


This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


__all__ = ('default_component_dict', )


def default_component_dict(dct={}, *, cmps, gas=None, soluble=None, solid=None):
    '''
    Fill out a dict with default values of components based on its property
    (values already in the dict will NOT be updated).
    
    Parameters
    ----------
    dct : dict
        Dict to be updated with the default values.
    g : value
        Default value for gas components (components in `cmps.gases`).
    s : value
        Default value for soluble components (components NOT in `cmps.gases` nor `cmps.solids`).
    x : value
        Default value for solid components (components in `cmps.solids`).
    '''
    gases = cmps.gases
    solids = cmps.solids
    for i in cmps:
        ID = i.ID
        if ID not in dct:
            if ID in gases:
                dct[ID] = gas
            elif ID in solids:
                dct[ID] = solid
            else:
                dct[ID] = soluble
    return dct