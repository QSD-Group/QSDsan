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

__all__ = ('format_number', )


def format_number(number):
    number = float(number)
    if number > 1e6 or number < 1e-4:
        return str(f'{number:.2E}')
    elif number == int(number):
        return str(int(number))
    else:
        return str(round(number, 4))