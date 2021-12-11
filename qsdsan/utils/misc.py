#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('register_with_prefix', )


def register_with_prefix(obj, prefix, ID):
    '''
    Register the object with a prefix (and a "_" between the prefix and the ID).

    Parameters
    ----------
    obj : obj
        The object to be registered, must has the `registry` attribute.
    prefix : str
        Prefix of the ID.
    ID : str
        The original ID.
    '''
    registry = obj.registry
    if ID == '' or None:
        data = registry.data
        ID = obj._take_ticket()
        full_ID = prefix+'_'+ID if prefix else ID
        while full_ID in data:
            ID = obj._take_ticket()
            full_ID = prefix+'_'+ID if prefix else ID
        registry.register(full_ID, obj)
    else:
        full_ID = prefix+'_'+ID if prefix else ID
        registry.register_safely(full_ID, obj)