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

import biosteam as bst
from datetime import timedelta

__all__ = ('time_printer',)


# %%

TicToc = bst.utils.TicToc

def time_printer(func):
    def inner(*args, **kwargs):
        try: print_time = kwargs['print_time']
        except: print_time = False
        if print_time:
            timer = TicToc()
            timer.tic()
        output = func(*args, **kwargs)
        if print_time:
            time = str(timedelta(seconds=round(timer.elapsed_time)))
            print(f'\nTotal time: {time}.')
        return output
    return inner

@time_printer
def foo(a=1, print_time=False):
    return a