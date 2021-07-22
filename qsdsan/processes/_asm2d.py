# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Cheung <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

import thermosteam as tmo    
# import os
# os.chdir("C:/Users/joy_c/Dropbox/PhD/Research/QSD/codes_developing/QSDsan")
from qsdsan import Components, Process, Processes
from ..utils import data_path

__all__ = ('cmps_asm2d', 'asm2d')


data_path += 'process_data/_asm2d.tsv'
# data_path = "qsdsan/data/process_data/_asm2d.tsv"

############# Components with default notation #############
cmps = Components.load_default()

S_A = cmps.S_Ac.copy('S_A')
S_ALK = cmps.S_CO3.copy('S_ALK')      # measured as g C
S_F = cmps.S_F.copy('S_F')
S_I = cmps.S_U_E.copy('S_I')
S_N2 = cmps.S_N2.copy('S_N2')
S_NH4 = cmps.S_NH4.copy('S_NH4')
S_NO3 = cmps.S_NO3.copy('S_NO3')
S_O2 = cmps.S_O2.copy('S_O2')
S_PO4 = cmps.S_PO4.copy('S_PO4')

X_AUT = cmps.X_AOO.copy('X_AUT')
X_H = cmps.X_OHO.copy('X_H')
X_I = cmps.X_U_OHO_E.copy('X_I')
X_MeOH = cmps.X_FeOH.copy('X_MeOH')
X_MeP = cmps.X_FePO4.copy('X_MeP')
X_PAO = cmps.X_PAO.copy('X_PAO')
X_PHA = cmps.X_PAO_PHA.copy('X_PHA')
X_PP = cmps.X_PAO_PP_Lo.copy('X_PP')
X_S = cmps.X_B_Subst.copy('X_S')

S_I.i_N = 0.01
S_F.i_N = 0.03
X_I.i_N = 0.02
X_S.i_N = 0.04
X_H.i_N = X_PAO.i_N = X_AUT.i_N = 0.07

S_I.i_P = 0.00
S_F.i_P = 0.01
X_I.i_P = 0.01
X_S.i_P = 0.01
X_H.i_P = X_PAO.i_P = X_AUT.i_P = 0.02

cmps_asm2d = Components([S_O2, S_F, S_A, S_NH4, S_NO3, S_PO4, S_I, S_ALK, S_N2,
                         X_I, X_S, X_H, X_PAO, X_PP, X_PHA, X_AUT, X_MeOH, X_MeP,
                         cmps.H2O])

cmps_asm2d.compile()
tmo.settings.set_thermo(cmps_asm2d)

############ Processes in ASM2d #################
params = ('f_SI', 'Y_H', 'f_XI', 'Y_PAO', 'Y_PO4', 'Y_PHA', 'Y_A', 
          'K_h', 'eta_NO3', 'eta_fe', 'K_O2', 'K_NO3', 'K_X',
          'mu_H', 'q_fe', 'eta_NO3_deni', 'b_H', 'K_F', 'K_fe', 'K_A',
          'K_NH4', 'K_P', 'K_ALK', 'q_PHA', 'q_PP', 'mu_PAO', 'b_PAO',
          'b_PP', 'b_PHA', 'K_PS', 'K_PP', 'K_MAX', 'K_IPP', 'K_PHA',
          'mu_AUT', 'b_AUT', 'K_O2_AUT', 'K_NH4_AUT', 'K_ALK_2',
          'k_PRE', 'k_RED')

asm2d = Processes.load_from_file(data_path,
                                 conserved_for=('COD', 'N', 'P', 'charge'),
                                 parameters=params,
                                 compile=False)

p12 = Process('anox_storage_PP',
              'S_PO4 + [Y_PHA]X_PHA + [?]S_NO3 -> X_PP + [?]S_N2 + [?]S_NH4 + [?]S_ALK',
              ref_component='X_PP',
              rate_equation='q_PP * S_O2/(K_O2+S_O2) * S_PO4/(K_PS+S_PO4) * S_ALK/(K_ALK+S_ALK) * (X_PHA/X_PAO)/(K_PHA+X_PHA/X_PAO) * (K_MAX-X_PP/X_PAO)/(K_IPP+K_MAX-X_PP/X_PAO) * X_PAO * eta_NO3 * K_O2/S_O2 * S_NO3/(K_NO3+S_NO3)',
              parameters=('Y_PHA', 'q_PP', 'K_O2', 'K_PS', 'K_ALK', 'K_PHA', 'eta_NO3', 'K_IPP', 'K_NO3'),
              conserved_for=('COD', 'N', 'P', 'NOD', 'charge'))

p14 = Process('PAO_anox_growth',
             '[1/Y_PAO]X_PHA + [?]S_NO3 + [?]S_PO4 -> X_PAO + [?]S_N2 + [?]S_NH4  + [?]S_ALK',
             ref_component='X_PAO',
             rate_equation='mu_PAO * S_O2/(K_O2 + S_O2) * S_NH4/(K_NH4 + S_NH4) * S_PO4/(K_P + S_PO4) * S_ALK/(K_ALK + S_ALK) * (X_PHA/X_PAO)/(K_PHA + X_PHA/X_PAO) * X_PAO * eta_NO3 * K_O2/S_O2 * S_NO3/(K_NO3 + S_NO3)',
             parameters=('Y_PAO', 'mu_PAO', 'K_O2', 'K_NH4', 'K_P', 'K_ALK', 'K_PHA', 'eta_NO3', 'K_NO3'),
             conserved_for=('COD', 'N', 'P', 'NOD', 'charge'))

asm2d.extend([p12, p14])
asm2d.compile()

# ASM2d typical values at 20 degree C
asm2d.set_parameters(
    f_SI = 0,                    # production of soluble inerts in hydrolysis = 0.0 gCOD/gCOD
    Y_H = 0.625,                 # heterotrophic yield = 0.625 gCOD/gCOD
    Y_PAO = 0.625,               # PAO yield = 0.625 gCOD/gCOD
    f_XI = 0.1,                  # fraction of inert COD generated in biomass lysis = 0.1 gCOD/gCOD
    Y_PO4 = 0.4,                 # PP requirement (PO4 release) per PHA stored = 0.4 gP/gCOD
    Y_PHA = 0.2,                 # PHA requirement for PP storage = 0.2 gCOD/gP
    Y_A = 0.24,                  # autotrophic yield = 0.24 gCOD/gN
    K_h = 3,                     # hydrolysis rate constant = 3.0 d^(-1)
    eta_NO3 = 0.6,               # reduction factor for anoxic activity and anoxic hydrolysis = 0.6
    eta_fe = 0.4,                # anaerobic hydrolysis reduction factor = 0.4
    K_O2 = 0.2,                  # O2 half saturation coefficient = 0.2 mgO2/L
    K_NO3 = 0.5,                 # nitrate half saturation coefficient = 0.5 mgN/L
    K_X = 0.1,                   # slowly biodegradable substrate half saturation coefficient for hydrolysis = 0.1 gCOD/gCOD
    mu_H = 6,                    # heterotrophic maximum specific growth rate = 6.0 d^(-1)
    q_fe = 3,                    # fermentation maximum rate = 3.0 d^(-1)
    eta_NO3_deni = 0.8,          # denitrification reduction factor = 0.8
    b_H = 0.4,                   # lysis and decay rate constant = 0.4 d^(-1)
    K_F = 4,                     # fermentable substrate half saturation coefficient for heterotrophic growth = 4.0 mgCOD/L
    K_fe = 4,                    # fermentable substrate half saturation coefficient for fermentation = 4.0 mgCOD/L
    K_A = 4,                     # VFA half saturation coefficient = 4.0 mgCOD/L
    K_NH4 = 0.05,                # ammonium (nutrient) half saturation coefficient = 0.05 mgN/L
    K_P = 0.01,                  # phosphorus (nutrient) half saturation coefficient = 0.01 mgP/L
    K_ALK = 0.1*12,              # alkalinity half saturation coefficient = 0.1 mol(HCO3-)/m^3 = 1.2 gC/m^3
    q_PHA = 3,                   # rate constant for storage of PHA = 3.0 d^(-1)
    q_PP = 1.5,                  # rate constant for storage of PP = 1.5 d^(-1)
    mu_PAO = 1,                  # PAO maximum specific growth rate = 1.0 d^(-1)
    b_PAO = 0.2,                 # PAO lysis rate = 0.2 d^(-1)
    b_PP = 0.2,                  # PP lysis rate = 0.2 d^(-1)
    b_PHA = 0.2,                 # PHA lysis rate = 0.2 d^(-1)
    K_PS = 0.2,                  # phosphorus half saturation coefficient for storage of PP = 0.2 mgP/L
    K_PP = 0.01,                 # PP half saturation coefficient for storage of PHA = 0.01 gP/gCOD (?). gCOD/gCOD in GPS-X
    K_MAX = 0.34,                # maximum ratio of X_PP/X_PAO = 0.34 gX_PP/gX_PAO
    K_IPP = 0.02,                # inhibition coefficient for PP storage = 0.02 gP/gCOD
    K_PHA = 0.01,                # PHA half saturation coefficient = 0.01 gCOD/gCOD
    mu_AUT = 1,                  # autotrophic maximum specific growth rate = 1.0 d^(-1)
    b_AUT = 0.15,                # autotrophic decay rate = 0.15 d^(-1)
    K_O2_AUT = 0.5,              # O2 half saturation coefficient for autotrophic growth = 0.5 mgO2/L
    K_NH4_AUT = 1,               # ammonium (substrate) half saturation coefficient for autotrophic growth = 1.0 mgN/L
    K_ALK_2 = 0.5*12,            # alkalinity half saturation coefficient for autotrophic growth and phosphates redissolution = 0.5 mol(HCO3-)/m^3 = 6.0 gC/m^3
    k_PRE = 1,                   # phosphorus precipitation with MeOH rate constant = 1.0 m^3/g/d
    k_RED = 0.6                  # redissoluation of phosphates rate constant = 0.6 d^(-1)
    )

# ASM2d typical values at 10 degree C
# asm2d.set_parameters(
#     f_SI = 0,                    # production of soluble inerts in hydrolysis = 0.0 gCOD/gCOD
#     Y_H = 0.625,                 # heterotrophic yield = 0.625 gCOD/gCOD
#     Y_PAO = 0.625,               # PAO yield = 0.625 gCOD/gCOD
#     f_XI = 0.1,                  # fraction of inert COD generated in biomass lysis = 0.1 gCOD/gCOD
#     Y_PO4 = 0.4,                 # PP requirement (PO4 release) per PHA stored = 0.4 gP/gCOD
#     Y_PHA = 0.2,                 # PHA requirement for PP storage = 0.2 gCOD/gP
#     Y_A = 0.24,                  # autotrophic yield = 0.24 gCOD/gN
#     K_h = 2,                     # hydrolysis rate constant = 2.0 d^(-1)
#     eta_NO3 = 0.6,               # reduction factor for anoxic activity and anoxic hydrolysis = 0.6
#     eta_fe = 0.4,                # anaerobic hydrolysis reduction factor = 0.4
#     K_O2 = 0.2,                  # O2 half saturation coefficient = 0.2 mgO2/L
#     K_NO3 = 0.5,                 # nitrate half saturation coefficient = 0.5 mgN/L
#     K_X = 0.1,                   # slowly biodegradable substrate half saturation coefficient for hydrolysis = 0.1 gCOD/gCOD
#     mu_H = 6,                    # heterotrophic maximum specific growth rate = 6.0 d^(-1)
#     q_fe = 1.5,                  # fermentation maximum rate = 1.5 d^(-1)
#     eta_NO3_deni = 0.8,          # denitrification reduction factor = 0.8
#     b_H = 0.2,                   # lysis and decay rate constant = 0.2 d^(-1)
#     K_F = 4,                     # fermentable substrate half saturation coefficient for heterotrophic growth = 4.0 mgCOD/L
#     K_fe = 4,                    # fermentable substrate half saturation coefficient for fermentation = 4.0 mgCOD/L
#     K_A = 4,                     # VFA half saturation coefficient = 4.0 mgCOD/L
#     K_NH4 = 0.05,                # ammonium (nutrient) half saturation coefficient = 0.05 mgN/L
#     K_P = 0.01,                  # phosphorus (nutrient) half saturation coefficient = 0.01 mgP/L
#     K_ALK = 0.1*12,              # alkalinity half saturation coefficient = 0.1 mol(HCO3-)/m^3 = 1.2 gC/m^3
#     q_PHA = 2,                   # rate constant for storage of PHA = 2.0 d^(-1)
#     q_PP = 1.0,                  # rate constant for storage of PP = 1.0 d^(-1)
#     mu_PAO = 0.67,               # PAO maximum specific growth rate = 0.67 d^(-1)
#     b_PAO = 0.1,                 # PAO lysis rate = 0.1 d^(-1)
#     b_PP = 0.1,                  # PP lysis rate = 0.1 d^(-1)
#     b_PHA = 0.1,                 # PHA lysis rate = 0.1 d^(-1)
#     K_PS = 0.2,                  # phosphorus half saturation coefficient for storage of PP = 0.2 mgP/L
#     K_PP = 0.01,                 # PP half saturation coefficient for storage of PHA = 0.01 gP/gCOD (?). gCOD/gCOD in GPS-X
#     K_MAX = 0.34,                # maximum ratio of X_PP/X_PAO = 0.34 gX_PP/gX_PAO
#     K_IPP = 0.02,                # inhibition coefficient for PP storage = 0.02 gP/gCOD
#     K_PHA = 0.01,                # PHA half saturation coefficient = 0.01 gCOD/gCOD
#     mu_AUT = 0.35,               # autotrophic maximum specific growth rate = 0.35 d^(-1)
#     b_AUT = 0.05,                # autotrophic decay rate = 0.05 d^(-1)
#     K_O2_AUT = 0.5,              # O2 half saturation coefficient for autotrophic growth = 0.5 mgO2/L
#     K_NH4_AUT = 1,               # ammonium (substrate) half saturation coefficient for autotrophic growth = 1.0 mgN/L
#     K_ALK_2 = 0.5*12,            # alkalinity half saturation coefficient for autotrophic growth and phosphates redissolution = 0.5 mol(HCO3-)/m^3 = 6.0 gC/m^3
#     k_PRE = 1,                   # phosphorus precipitation with MeOH rate constant = 1.0 m^3/g/d
#     k_RED = 0.6                  # redissoluation of phosphates rate constant = 0.6 d^(-1)
#     )