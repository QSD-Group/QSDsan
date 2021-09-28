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

__all__ = ('AttrGetter', 'FuncGetter')



class AttrGetter:
    __slots__ = ('obj', 'attr', 'hook', 'hook_param')
    def __init__(self, obj, attr, hook=lambda i: i, hook_param=None):
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