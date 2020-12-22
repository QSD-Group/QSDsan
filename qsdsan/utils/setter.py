#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''


# %%

from biosteam.evaluation.evaluation_tools.parameter import Setter


__all__ = ('DictAttrSetter',)

class DictAttrSetter(Setter):
    __slots__ = ('obj', 'dict_attr', 'key')
    def __init__(self, obj, dict_attr, key):
        self.dict_attr = getattr(obj, dict_attr)
        self.key = key

    def __call__(self, value):
        self.dict_attr[self.key] = value