#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>
    
    Joy Zhang <joycheung1994@gmail.com>

    Ga-Yeong Kim <gayeong1225@gmail.com>

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

def Mbamba_rhos(S_Ca, S_Mg, co3, nh4, po4, hpo4, # mass concentrations
                X_CaCO3, X_struv, X_newb, X_ACP, X_MgCO3, 
                k_mmp, Ksp):        # mass-based solubility product
    rhos = [0]*5
    if S_Ca > 0 and co3 > 0:
        SI = (S_Ca * co3 / Ksp[0])**(1/2)
        if SI > 1: rhos[0] = X_CaCO3 * (SI-1)**2
    if S_Mg > 0 and nh4 > 0 and po4 > 0:
        SI = (S_Mg * nh4 * po4 / Ksp[1])**(1/3)
        if SI > 1: rhos[1] = X_struv * (SI-1)**3
    if S_Mg > 0 and hpo4 > 0:
        SI = (S_Mg * hpo4 / Ksp[2])**(1/2)
        if SI > 1: rhos[2] =  X_newb * (SI-1)**2
    if S_Ca > 0 and po4 > 0:
        SI = (S_Ca**3 * po4**2 / Ksp[3])**(1/5)
        if SI > 1: rhos[3] = X_ACP * (SI-1)**2
    if S_Mg > 0 and co3 > 0:
        SI = (S_Mg * co3 / Ksp[4])**(1/2)
        if SI > 1: rhos[4] = X_MgCO3 * (SI-1)**2
    rhos[:] *= k_mmp
    return rhos     # [mass conc]/d

def Musvoto_rhos(Ca, Mg, co3, nh4, po4, hpo4, k_mmp, Ksp): # molar concentrations & solubility product
    rhos = [0]*5
    if Ca > 0 and co3 > 0:
        SI = (Ca * co3)**(1/2) - Ksp[0]**(1/2)
        if SI > 0: rhos[0] = SI**2
    if Mg > 0 and nh4 > 0 and po4 > 0:
        SI = (Mg * nh4 * po4)**(1/3) - Ksp[1]**(1/3)
        if SI > 0: rhos[1] = SI**3
    if Mg > 0 and hpo4 > 0:
        SI = (Mg * hpo4)**(1/2) - Ksp[2]**(1/2)
        if SI > 0: rhos[2] = SI**2
    if Ca > 0 and po4 > 0:
        SI = (Ca**3 * po4**2)**(1/5) - Ksp[3]**(1/5)
        if SI > 0: rhos[3] = SI**2
    if Mg > 0 and co3 > 0:
        SI = (Mg * co3)**(1/2) - Ksp[4]**(1/2)
        if SI > 0: rhos[4] = SI**2
    rhos[:] *= k_mmp
    return rhos   # [molar conc]/d

#%%
from ._aeration import *
from ._aerobic_digestion_addon import *
from ._asm1 import *
from ._asm2d import *
from ._adm1 import *
from ._adm1_p_extension import *
# from ._madm1 import *
from ._decay import *
from ._kinetic_reaction import *
from ._pm2 import *
from ._pm2asm2d import *

from . import (
    _aeration,
    _aerobic_digestion_addon,
    _asm1,
    _asm2d,
    _adm1,
    _adm1_p_extension,
    # _madm1,
    _decay,
    _kinetic_reaction,
    _pm2,
    _pm2asm2d,
    )

__all__ = (
    *_aeration.__all__,
    *_aerobic_digestion_addon.__all__,
    *_asm1.__all__,
    *_asm2d.__all__,
    *_adm1.__all__,
    *_adm1_p_extension.__all__,
    # *_madm1.__all__,
    *_decay.__all__,
    *_kinetic_reaction.__all__,
    *_pm2.__all__,
    *_pm2asm2d.__all__,
    )