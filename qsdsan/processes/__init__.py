#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from numpy import arange, cumprod, exp

def ion_speciation(h_ion, *Kas):
    n = len(Kas)
    out = h_ion ** arange(n, -1, -1) * cumprod([1.0, *Kas])
    return out/sum(out)

substr_inhibit = Monod = lambda S, K: S/(S+K)

non_compet_inhibit = lambda S, K: K/(K+S)

grad_non_compet_inhibit = lambda S, K: -K/(K+S)**2

grad_substr_inhibit = lambda S, K: K/(K+S)**2

def mass2mol_conversion(cmps):
    '''conversion factor from kg[measured_as]/m3 to mol[component]/L'''
    return cmps.i_mass / cmps.chem_MW

R = 8.3145e-2 # Universal gas constant, [bar/M/K]

def T_correction_factor(T1, T2, delta_H):
    """
    Returns temperature correction factor for equilibrium constants based on
    the Van't Holf equation.

    Parameters
    ----------
    T1 : float
        Base temperature, in K.
    T2 : float
        Actual temperature, in K.
    delta_H : float
        Heat of reaction, in J/mol.
    """
    if T1 == T2: return 1
    return exp(delta_H/(R*100) * (1/T1 - 1/T2))  # R converted to SI

class TempState:
    def __init__(self):
        self.data = {}

#%%
from ._aeration import *
from ._asm1 import *
from ._asm2d import *
from ._adm1 import *
from ._adm1_p_extension import *
# from ._madm1 import *
from ._decay import *
from ._kinetic_reaction import *
from ._pm2 import *

from . import (
    _aeration,
    _asm1,
    _asm2d,
    _adm1,
    _adm1_p_extension,
    # _madm1,
    _decay,
    _kinetic_reaction,
    _pm2
    )

__all__ = (
    *_aeration.__all__,
    *_asm1.__all__,
    *_asm2d.__all__,
    *_adm1.__all__,
    *_adm1_p_extension.__all__,
    # *_madm1.__all__,
    *_decay.__all__,
    *_kinetic_reaction.__all__,
    *_pm2.__all__,
    )
