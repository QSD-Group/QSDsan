# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

# from thermosteam.utils import chemicals_user
from chemicals.elements import molecular_weight as get_mw
from qsdsan import Component, Components, Process, Processes, CompiledProcesses
import numpy as np, qsdsan.processes as pc, qsdsan as qs
# from qsdsan.utils import ospath, data_path
# from scipy.optimize import brenth
# from warnings import warn

# from qsdsan.processes import create_adm1_cmps


__all__ = ('create_madm1_cmps',
           )

#%%
# C_mw = get_mw({'C':1})
# N_mw = get_mw({'N':1})
# P_mw = get_mw({'P':1})
# S_mw = get_mw({'S':1})
Fe_mw = get_mw({'Fe':1})
O_mw = get_mw({'O':1})

def create_madm1_cmps(set_thermo=True):
    
    # Components from the original ADM1
    # *********************************
    _cmps = pc.create_adm1_cmps(False)
    S_aa = _cmps.S_aa
    X_pr = _cmps.X_pr
    S_aa.i_C = X_pr.i_C = 0.36890
    S_aa.i_N = X_pr.i_N = 0.11065
    S_aa.i_mass = X_pr.i_mass = 1/1.35566
    
    S_fa = _cmps.S_fa
    S_fa.formula = 'C25H52O3'
    S_bu = _cmps.S_bu
    S_bu.formula = 'C4H8O2'
    S_pro = _cmps.S_pro
    S_pro.formula = 'C3H6O2'
    S_ac = _cmps.S_ac
    S_ac.formula = 'C2H4O2'
    
    S_I = _cmps.S_I
    X_I = _cmps.X_I
    S_I.i_C = X_I.i_C = 0.36178
    S_I.i_N = X_I.i_N = 0.06003
    S_I.i_P = X_I.i_P = 0.00649
    S_I.i_mass = X_I.i_mass = 1/1.54100

    X_ch = _cmps.X_ch
    X_ch.formula = 'C24H48O24'
    # _cmps.X_li.formula = 'C64H119O7.5P'
    X_li = X_pr.copy('X_li')
    X_li.i_C = 0.26311
    X_li.i_N = 0.
    X_li.i_P = 0.01067
    X_li.i_mass = 1/2.81254
    
    adm1_biomass = (_cmps.X_su, _cmps.X_aa, _cmps.X_fa, _cmps.X_c4, _cmps.X_pro, _cmps.X_ac, _cmps.X_h2)
    for bio in adm1_biomass:
        # bio.formula = 'C5H7O2NP0.113'
        bio.i_C = 0.36612
        bio.i_N = 0.08615
        bio.i_P = 0.02154
        bio.i_mass = 1/1.39300
    
    # P related components from ASM2d
    # *******************************
    asm_cmps = pc.create_asm2d_cmps(False)
    X_PHA = asm_cmps.X_PHA
    X_PHA.formula = '(C2H4O)n'
    # X_PHA.i_C = 0.3
    # X_PHA.i_mass = 0.55
    
    X_PAO = _cmps.X_su.copy('X_PAO')
    X_PAO.description = 'Phosphorus-accumulating organism biomass'
    
    # Additional components for P, S, Fe extensions
    # *********************************************
    S_IP = asm_cmps.S_PO4.copy('S_IP')
    S_K = Component.from_chemical('S_K', chemical='K+',
                                  description='Potassium ion',
                                  measured_as='K',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
    S_Mg = Component.from_chemical('S_Mg', chemical='Mg2+',
                                  description='Magnesium ion',
                                  measured_as='Mg',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
    S_SO4 = Component.from_chemical('S_SO4', chemical='SO4-2',
                                  description='Sulfate',
                                  measured_as='S',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
    S_IS = Component.from_chemical('S_IS', chemical='H2S',
                                  description='Hydrogen sulfide',
                                  measured_as='COD',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
    
    X_hSRB = X_PAO.copy('X_hSRB')
    X_hSRB.description = 'sulfate-reducing biomass, utilizing H2'
    X_aSRB = X_PAO.copy('X_aSRB')
    X_aSRB.description = 'sulfate-reducing biomass, utilizing acetate'
    X_pSRB = X_PAO.copy('X_pSRB')
    X_pSRB.description = 'sulfate-reducing biomass, utilizing propionate'
    X_c4SRB = X_PAO.copy('X_c4SRB')
    X_c4SRB.description = 'sulfate-reducing biomass, utilizing butyrate and valerate'
    
    S_S0 = Component.from_chemical('S_S0', chemical='S',
                                  description='Sulfur',
                                  measured_as='COD',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
    S_Fe3 = Component.from_chemical('S_Fe3', chemical='Fe3+',
                                  description='Iron (III)',
                                  measured_as='Fe',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
    S_Fe2 = Component.from_chemical('S_Fe2', chemical='Fe2+',
                                  description='Iron (II)',
                                  measured_as='Fe',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
    S_Fe2.i_COD = 0.5*O_mw/Fe_mw
    S_Fe2.measured_as = 'COD'
    
    cmps_madm1 = Components([_cmps.S_su, S_aa, S_fa, _cmps.S_va, S_bu, 
                             S_pro, S_ac, _cmps.S_h2, _cmps.S_ch4, 
                             _cmps.S_IC, _cmps.S_IN, S_IP, S_I, 
                             X_ch, X_pr, X_li, *adm1_biomass, X_I,
                             X_PHA, asm_cmps.X_PP, X_PAO, S_K, S_Mg, 
                             S_SO4, S_IS, X_hSRB, X_aSRB, X_pSRB, X_c4SRB,
                             S_S0, S_Fe3, S_Fe2,
                             _cmps.H2O])
    cmps_madm1.default_compile()
    
    if set_thermo: qs.set_thermo(cmps_madm1)
    return cmps_madm1

#%%
