#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>

    Joy Zhang <joycheung1994@gmail.com>

    Lewis Rowles <stetsonsc@gmail.com>

    Hannah Lohman <hlohman94@gmail.com>

    Tori Morgan <vlmorgan@illinois.edu>

    Shion Watabe <swatabe2@illinois.edu>

    Lane To <lane20@illinois.edu>

    Smiti Mittal <smitimittal@gmail.com>

    Anna Kogler <akogler@stanford.edu>

    Jianan Feng <jiananf2@illinois.edu>

    Saumitra Rai <raisaumitra9@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''
# %%
# dydt_cstr must be defined before subpackage imports so subpackages can
# retrieve it via `from .. import dydt_cstr` during their own initialization.
from numba import njit
@njit(cache=True)
def dydt_cstr(QC_ins, QC, V, _dstate):
    Q_ins = QC_ins[:, -1]
    C_ins = QC_ins[:, :-1]
    _dstate[-1] = 0
    _dstate[:-1] = (Q_ins @ C_ins - sum(Q_ins)*QC[:-1])/V

#%%
from .bst import *
from .static import *
from .dynamic import *

from . import bst, static, dynamic

__all__ = (
    *bst.__all__,
    *static.__all__,
    *dynamic.__all__,
    'bst',
    'static',
    'dynamic',
)
