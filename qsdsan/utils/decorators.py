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

# %%

import biosteam as bst
from datetime import timedelta

__all__ = ('time_printer',
           # 'NonNegativeFloat', 'NonNgeativeInt', 'Fraction',
           )


# %%

# =============================================================================
# Allow functions to print execution time with a `print_time` kwargs
# =============================================================================

TicToc = bst.utils.TicToc

def time_printer(func):
    '''
    Allow functions to print execution time with a `print_time` kwargs.

    Examples
    --------
    >>> from qsdsan.utils.decorators import time_printer
    >>> @time_printer
    ... def foo(a=1, print_time=False):
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


# %%

# =============================================================================
# Checkers
# =============================================================================

'''
Checkers are decorators to check if the given value satisfy a set of constrains,
a great reference can be found online. [1]_

.. note::
    Using these decorators are much slower than just adding the checks when defining
    them, but probably won't matter since it's still < 1 microseconds.

References
----------
.. [1] https://www.programiz.com/python-programming/decorator

'''

def NonNegativeFloat(obj):
    def inner_setter(outer_setter):
        data = outer_setter()
        attr = '_' + outer_setter.__name__
        if data < 0:
            raise ValueError(f'Value for {attr[1:]} cannot be negative.')
        setattr(obj, attr, float(data))
    return inner_setter


def NonNgeativeInt(obj):
    def inner_setter(outer_setter):
        data = outer_setter()
        attr = '_' + outer_setter.__name__
        if data < 0:
            raise ValueError(f'Value for {attr[1:]} cannot be negative.')
        setattr(obj, attr, int(data))
    return inner_setter


def Fraction(obj):
    def inner_setter(outer_setter):
        data = outer_setter()
        attr = '_' + outer_setter.__name__
        if not 0 <=data <=1:
            raise ValueError(f'Value for {attr[1:]} must be between [0, 1].')
        setattr(obj, attr, float(data))
    return inner_setter