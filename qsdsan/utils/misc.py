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

__all__ = ('copy_attr', 'register_with_prefix', )


def copy_attr(new, original, skip=(), same=()):
    '''
    Set the attributes of a new object based on an original one:

        - If one attribute is in `skip`, it will not be copied to the new object.
        - If one attribute is in `same`, the attribute of the new object will be \
        the same as the original object.
        - For remaining attributes, if it has :func:`copy`, then the attribute \
        of the new object will be set as the copy of the original one; otherwise, \
        it will be the same as the original one.

    Parameters
    ----------
    new : obj
        The new object.
    origin : obj
        The original object.
    skip : Iterable
        Attributes that will not be copied.
    same : Iterable
        Attributes that will be the same for the original one and the copy.
    '''

    for slot in original.__slots__:
        if slot in skip:
            continue
        else:
            value = getattr(original, slot)
            if slot in same:
                setattr(new, slot, value)
                return new
            else:
                if hasattr(value, 'copy'):
                    new_value = value.copy()
                else:
                    new_value = value
            setattr(new, slot, new_value)
    return new


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