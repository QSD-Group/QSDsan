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


__all__ = (
    'AttrGetter', 'FuncGetter',
    'AttrSetter', 'AttrFuncSetter', 'DictAttrSetter',
    'copy_samples',
    )


# %%

# =============================================================================
# Getter Functions
# =============================================================================

class AttrGetter:
    __slots__ = ('obj', 'attr', 'hook', 'hook_param')
    def __init__(self, obj, attr, hook=lambda i: i, hook_param=()):
        self.obj = obj
        self.attr = attr
        self.hook = hook
        self.hook_param = hook_param

    def __call__(self):
        return self.hook(getattr(self.obj, self.attr), *self.hook_param)

# # The below one needs updating
# class AttrGetter:
#     __slots__ = ('obj', 'attrs')
#     def __init__(self, obj, attrs):
#         self.obj = obj
#         if isinstance(attrs, str):
#             attrs = (attrs,)
#         self.attrs = attrs

#     def __call__(self):
#         if len(self.attrs) == 1:
#             print(getattr(self.obj, self.attrs[0]))
#             return getattr(self.obj, self.attrs[0])
#         return [getattr(self.obj, attr) for attr in self.attrs]


class FuncGetter:
    __slots__ = ('func', 'params')
    def __init__(self, func, params):
        self.func = func
        self.params = params

    def __call__(self):
        return self.func(*self.params)


# %%

# =============================================================================
# Setter Functions
# =============================================================================

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


# %%

# =============================================================================
# Sample Manipulation
# =============================================================================

def copy_samples(original, new, exclude=()):
    '''
    Copy samples of the shared parameters in the original model to the new model.
    Parameters in `exclude` will be excluded (i.e., not copied).
    '''
    col0 = original.table.columns.get_level_values(1)[:len(original.parameters)]
    col1 = new.table.columns.get_level_values(1)[:len(new.parameters)]
    shared = col0.intersection(col1)
    shared = shared.difference([i.name_with_units for i in exclude])
    idx0 = original.table.columns.get_locs([slice(None), shared])
    idx1 = new.table.columns.get_locs([slice(None), shared])
    new.table[new.table.columns[idx1]] = new._samples[:, idx1] \
        = original._samples[:, idx0]