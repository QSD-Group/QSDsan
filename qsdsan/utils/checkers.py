#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''


__all__ = ('NonNegativeFloat', 'NonNgeativeInt', 'Fraction')


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




