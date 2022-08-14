#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import pandas as pd

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
    def __init__(self, func, params=()):
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

def copy_samples(original, new, exclude=(), only_same_baseline=False):
    '''
    Copy samples of the shared parameters in the original model to the new model.
    Parameters in `exclude` will be excluded (i.e., will not be copied).

    Parameters
    ----------
    original : obj
        Original model where the samples will be copied from.
    new : obj
        New model whose samples of the shared parameters
        will be copied from the original model.
    exclude : tuple(obj)
        Parameters that will be excluded from copying.
    only_same_baseline : bool
        If True, will only copy parameters with the same
        name, units, and baseline values.

    Examples
    --------
    Create two models with shared and different parameters

    >>> from qsdsan.utils import load_example_model, copy_samples
    >>> original = load_example_model()
    >>> # You might get some replacing warnings as we are making two exact systems,
    >>> # it's OK
    >>> new = load_example_model()
    >>> # Sort the parameter orders
    >>> for model in (original, new):
    ...    model.parameters = sorted(model.parameters, key=lambda p: p.name)


    Let's make the 1st/2nd parameters of the original model the same as
    the 0th/1st parameters of the new model.

    >>> original.parameters = original.parameters[:3]
    >>> original.parameters # doctest: +SKIP
    (<Parameter: [HXutility-H1] Heat exchanger temperature (K)>,
     <Parameter: [Mix tank-M1] Mix tank mixer power usage (kW/m3)>,
     <Parameter: [Mix tank-M1] Mix tank retention time (hr)>)
    >>> new.parameters = new.parameters[1:4]
    >>> new.parameters # doctest: +SKIP
    (<Parameter: [Mix tank-M1] Mix tank mixer power usage (kW/m3)>,
     <Parameter: [Mix tank-M1] Mix tank retention time (hr)>,
     <Parameter: [Pump-P1] Pump design head (kPa)>)

    Before copying, samples of the shared parameter are not the same

    >>> original_samples = original.sample(N=100, rule='L')
    >>> original.load_samples(original_samples)
    >>> new_samples = new.sample(N=100, rule='L')
    >>> new.load_samples(new_samples)
    >>> (original.table.values[:, 1:3]==new.table.values[:, 0:2]).all()
    False

    Let's copy the samples, but exclude the "Mix tank retention time"

    >>> copy_samples(original, new, exclude=(original.parameters[2],))
    >>> # After copying, all of the samples of the "Mix tank mixer power usage"
    >>> # parameter are the same
    >>> (original.table.values[:, 1]==new.table.values[:, 0]).all()
    True
    >>> # But the samples are not the same for the excluded parameter
    >>> # "Mix tank retention time" parameter
    >>> (original.table.values[:, 2]==new.table.values[:, 1]).all()
    False

    When `only_same_baseline` is True, only values for samples with the same
    name and unit will not be copied if their baseline values are different.

    >>> # Let's change the baseline values for the "Mix tank retention time" parameter
    >>> # to be different for the two models
    >>> new.parameters[1].baseline = original.parameters[2].baseline + 10
    >>> copy_samples(original, new, exclude=(), only_same_baseline=True)
    >>> # Samples are not copied even if they have the same names due to the unequal baselines
    >>> (original.table.values[:, 2]==new.table.values[:, 1]).any()
    False
    '''
    try: iter(exclude)
    except: exclude = (exclude,)
    col0 = original.table.columns.get_level_values(1)[:len(original.parameters)]
    col1 = new.table.columns.get_level_values(1)[:len(new.parameters)]
    shared = col0.intersection(col1)
    shared = shared.difference([i.name_with_units for i in exclude])
    if only_same_baseline:
        original_b = {p.name_with_units: p.baseline for p in original.parameters
                      if p.name_with_units in shared}
        original_b = pd.DataFrame.from_dict(original_b, columns=['val'], orient='index').sort_index()
        new_b = {p.name_with_units: p.baseline for p in new.parameters
                 if p.name_with_units in shared}
        new_b = pd.DataFrame.from_dict(new_b, columns=['val'], orient='index').sort_index()
        shared = original_b[original_b.val==new_b.val].index
    idx0 = original.table.columns.get_locs([slice(None), shared])
    idx1 = new.table.columns.get_locs([slice(None), shared])
    new.table[new.table.columns[idx1]] = new._samples[:, idx1] \
        = original._samples[:, idx0]