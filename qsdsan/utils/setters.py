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


# %%

__all__ = ('AttrSetter', 'AttrFuncSetter', 'DictAttrSetter')

setattr = setattr
getattr = getattr
isinstance = isinstance

class AttrSetter:
    __slots__ = ('obj', 'attrs')
    def __init__(self, obj, attrs):
        self.obj = obj
        if isinstance(attrs, str):
            attrs = (attrs,)
        self.attrs = attrs
        
    def __call__(self, value):
        for attr in self.attrs:
            setattr(self.obj, attr, value)

class AttrFuncSetter:
    __slots__ = ('obj', 'attrs', 'funcs')
    def __init__(self, obj, attrs, funcs):
        self.obj = obj
        if isinstance(attrs, str):
            attrs = (attrs,)
        if callable(funcs):
            funcs = (funcs,)
        self.attrs = attrs
        self.funcs = funcs
        
    def __call__(self, value):
        attrs = self.attrs
        funcs = self.funcs
        obj = self.obj
        
        if len(funcs) == 1:
            func = funcs[0]
            for attr in attrs:
                setattr(obj, attr, func(value))
        elif len(funcs) == len(attrs):
            for num, func in enumerate(funcs):
                setattr(obj, attrs[num], func(value))
        else:
            raise ValueError('Number of functions does not match number of attributes.')

class DictAttrSetter:
    __slots__ = ('obj', 'dict_attr', 'keys')
    def __init__(self, obj, dict_attr, keys):
        self.dict_attr = getattr(obj, dict_attr)
        if isinstance(keys, str):
            keys = (keys,)
        self.keys = keys

    def __call__(self, value):
        for key in self.keys:
            self.dict_attr[key] = value

