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

import biosteam as bst
from datetime import timedelta

__all__ = ('time_printer',)


# %%

# =============================================================================
# Allow functions to print execution time with a `print_time` kwargs
# =============================================================================

TicToc = bst.utils.TicToc

def time_printer(func):
    '''
    Allow functions to print execution time with a `print_time` kwargs.
    
    Examples
    --------
    >>> from qsdsan.utils.decorators import time_printer
    >>> @time_printer
    ... def foo(a=1, print_time=False):
    ...    return a
    >>> # This will not print run time
    >>> print(foo(a=5))
    5
    >>> # This will print run time
    >>> print(foo(a=5, print_time=True))
    Total time: 0:00:00.
    5
    '''
    
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
    inner.__doc__ = func.__doc__
    return inner

