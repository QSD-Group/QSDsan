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

__all__ = ('format_number', )


def format_number(number):
    number = float(number)
    if abs(number) == 0 or str(number) == 'nan':
        return '0'
    elif abs(number) > 1e6 or abs(number) < 1e-4:
        return str(f'{number:.2E}')
    elif 1e-3 < abs(number) < 0.01:
        return str(round(number, 4))
    elif number == int(number):
        return str(int(number))
    else:
        return str(round(number, 2))