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


from collections import defaultdict
from weakref import WeakKeyDictionary

__all__ = ('Validator', 'NonNegativeFloat', 'NonNgeativeInt', 'Fraction',
           'BareModule')


# %%

class Validator:
    '''
    Descriptors can be used to make resusable property logics,
    a great reference can be found online. [1]_
    
    .. note::
        
        Using Descriptors to check value properties is slightly slower than using
        the decorators functions in qsdsan.utils.checkers.
    
    
    References
    ----------
    .. [1] https://nbviewer.jupyter.org/urls/gist.github.com/ChrisBeaumont/5758381/raw/descriptor_writeup.ipynb
     
    '''
    
    __slots__ = ['default', 'name', 'instance', 'owner']
    
    def __init__(self, default=None, name=None):
        self.default = default
        self.name = name
        self.data = WeakKeyDictionary()
        
    def __get__(self, instance, owner):
        return self.data.get(instance, self.default)


class NonNegativeFloat(Validator):
    '''For non-negative floats.'''
    
    def __set__(self, instance, value):
        if value < 0:
            if self.name:
                raise ValueError(f'Value for {self.name} cannot be negative.')
            else:
                raise ValueError('Value cannot be negative.')
        self.data[instance] = float(value)


class NonNgeativeInt(Validator):
    '''For non-negative integers.'''
    
    def __set__(self, instance, value):
        if value < 0:
            if self.name:
                raise ValueError(f'Value for {self.name} cannot be negative.')
            else:
                raise ValueError('Value cannot be negative.')
        self.data[instance] = int(value)


class Fraction(Validator):
    '''For values in [0,1].'''
    
    def __set__(self, instance, value):
        if not 0 <= value <= 1:
            if self.name:
                raise ValueError(f'Value for {self.name} must be in [0,1].')
            else:
                raise ValueError('Value must be in [0,1].')
        self.data[instance] = float(value)
        


# # Not in use
# class BareModule:
#     '''For unit BM (bare module factors).'''
    
#     def __get__(self, obj, objtype=None):
#         # Copy the _BM dict from the class attribute to prevent overwriting,
#         # need to in both __get__ and __set__ as not sure which one will be called first
#         if '_BM' not in obj.__dict__:
#             _BM = defaultdict(lambda:1)
#             _BM.update(obj._BM)
#             obj._BM = _BM
        
#         # If don't want defaults
#         # if '_BM' not in obj.__dict__:
#         #     obj._BM = obj._BM.copy()

#         return obj._BM
    
#     def __set__(self, obj, BM:dict):
#         # Copy the _BM dict from the class attribute to prevent overwriting,
#         # need to in both __get__ and __set__ as not sure which one will be called first
#         if '_BM' not in obj.__dict__:
#             _BM = defaultdict(lambda:1)
#             _BM.update(obj._BM)
#             obj._BM = _BM
        
#         # If don't want defaults
#         # if '_BM' not in obj.__dict__:
#         #     obj._BM = obj._BM.copy()
        
#         # if not isinstance(BM, dict):
#         #     raise TypeError('`BM` must be a dict, not {type(i).__name__}.')
        
#         # for i in self.BM.keys():
#         #     if i not in self.purhcase_costs.keys():
#         #         raise ValueError(f'The item "{i}" does not exist in `purchase_costs`.')

#         obj._BM.update(BM)



        
    