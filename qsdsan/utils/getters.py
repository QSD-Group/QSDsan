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

__all__ = ('FuncGetter',)


class FuncGetter:
    __slots__ = ('func', 'params')
    def __init__(self, func, params):
        self.func = func
        self.params = params

    def __call__(self):
        return self.func(*self.params)