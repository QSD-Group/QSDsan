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

from datetime import timedelta
from biosteam.utils import TicToc

__all__ = (
    'copy_attr',
    'register_with_prefix',
    'time_printer',
    'ords',
    )


def copy_attr(new, original, skip=(), same=(), slots=None):
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
    slots : Iterable[str]
        All fields of the original object, will be set to `original.__slots__` if not provided.
    '''
    slots = slots or original.__slots__
    for slot in slots:
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


# Allow functions to print execution time with a `print_time` kwarg
def time_printer(func):
    '''
    Allow functions to print execution time with a `print_time` kwarg.

    Examples
    --------
    >>> from qsdsan.utils import time_printer
    >>> @time_printer
    ... def foo(a=1, **kwargs):
    ...     return a
    >>> # This will print run time
    >>> print(foo(a=5))
    function `foo`
    Total time: 0:00:00.
    5
    >>> # This will NOT print run time
    >>> print(foo(a=5, print_time=False))
    5
    '''
    def inner(*args, **kwargs):
        print_time = kwargs.get('print_time')
        if print_time is not False:
            timer = TicToc()
            timer.tic()
        output = func(*args, **kwargs)
        if print_time is not False:
            time = str(timedelta(seconds=round(timer.elapsed_time)))
            name = str(func).split(' ')[1]
            print(f'function `{name}`')
            print(f'Total time: {time}.')
        return output
    inner.__doc__ = func.__doc__
    return inner


# Return the sum of unicode of a string, more for fun
def ords(string):
    string = str(string)
    added = sum(ord(i) for i in string)
    return added