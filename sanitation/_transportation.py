#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Sanitation Explorer: Sustainable design of non-sewered sanitation technologies
Copyright (C) 2020, Sanitation Explorer Development Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-for-WaSH/sanitation/blob/master/LICENSE.txt
for license details.

'''


# %%

from . import ImpactIndicator, ImpactItem
from ._units_of_measure import auom
from .utils.formatting import format_number as f_num

indicators = ImpactIndicator._indicators

__all__ = ('Transportation',)


class Transportation:
    '''Transportation cost and environmental impacts.'''
    '''TO BE ADDED'''