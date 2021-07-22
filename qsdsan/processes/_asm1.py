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
from qsdsan import Components, Processes#, Process
from ..utils import data_path

__all__ = ('cmps_asm1', 'asm1')


data_path += 'process_data/_asm1.tsv'
# data_path = "qsdsan/data/process_data/_asm1.tsv"

############# Components with default notation #############
cmps = Components.load_default()

S_I = cmps.S_U_Inf.copy('S_I')
S_I.description = 'Soluble inert organic matter'

X_I = cmps.X_U_Inf.copy('X_I')
X_I.description = 'Particulate inert organic matter'

S_S = cmps.S_F.copy('S_S')
S_S.description = 'Readily biodegradable substrate'
S_S.i_N = 0

X_S = cmps.X_B_Subst.copy('X_S')
X_S.description = 'Slowly biodegradable substrate'
X_S.i_N = 0

X_BH = cmps.X_OHO.copy('X_BH')
X_BH.description = 'Active heterotrophic biomass'

X_BA = cmps.X_AOO.copy('X_BA')
X_BA.description = 'Active autotrophic biomass'
X_BH.i_N = X_BA.i_N = 0.086     # i_XB

X_P = cmps.X_U_OHO_E.copy('X_P')
X_P.description = 'Particulate products arising from biomass decay'
X_P.i_N = 0.06                  # i_XP

S_O = cmps.S_O2.copy('S_O')

S_NO = cmps.S_NO3.copy('S_NO')
S_NO.description = 'Nitrate and nitrite nitrogen'

S_NH = cmps.S_NH4.copy('S_NH')

S_ND = cmps.S_F.copy('S_ND')
S_ND.description = 'Soluble biodegradable organic nitrogen'
S_ND.measured_as = 'N'
S_ND.i_COD = S_ND.i_C = S_ND.i_P = 0

X_ND = cmps.X_B_Subst.copy('X_ND')
X_ND.description = 'Particulate biodegradable organic nitrogen'
X_ND.measured_as = 'N'
X_ND.i_COD = X_ND.i_C = X_ND.i_P = 0

S_ALK = cmps.S_CO3.copy('S_ALK')      # measured as g C

# add S_N2 to close mass balance
# add water for the creation of WasteStream objects
cmps_asm1 = Components([S_I, S_S, X_I, X_S, X_BH, X_BA, X_P, 
                        S_O, S_NO, S_NH, S_ND, X_ND, S_ALK, 
                        cmps.S_N2, cmps.H2O])         

cmps_asm1.compile()
tmo.settings.set_thermo(cmps_asm1)

############ Processes in ASM1 #################
params = ('Y_H', 'Y_A', 'f_P', 
          'mu_H', 'K_S', 'K_O_H', 'K_NO', 'b_H', 
          'mu_A', 'K_NH', 'K_O_A', 'b_A',
          'eta_g', 'k_a', 'k_h', 'K_X', 'eta_h')

asm1 = Processes.load_from_file(data_path,
                                conserved_for=('COD', 'charge', 'N'),
                                parameters=params,
                                compile=True)

# ASM1 typical values at 20 degree C
asm1.set_parameters(
    Y_A = 0.24,                  # autotrophic yield = 0.24 gCOD/gN
    Y_H = 0.67,                  # heterotrophic yield = 0.67 gCOD/gCOD
    f_P = 0.08,                  # fraction of biomass yielding particulate products = 0.08, unitless 
    mu_H = 6,                    # heterotrophic maximum specific growth rate = 6.0 d^(-1)
    K_S = 20,                    # readily biodegradable substrate half saturation coefficient = 20.0 gCOD/m3
    K_O_H = 0.2,                 # O2 half saturation coefficient = 0.2 gO2/m3
    K_NO = 0.5,                  # nitrate half saturation coefficient = 0.5 gN/m3
    b_H = 0.62,                  # heterotrophic biomass decay rate constant = 0.62 d^(-1)
    eta_g = 0.8,                 # reduction factor for anoxic growth of heterotrophs = 0.8, unitless
    eta_h = 0.4,                 # anoxic hydrolysis reduction factor = 0.4, unitless
    k_h = 3.0,                   # hydrolysis rate constant = 3.0 d^(-1)
    K_X = 0.03,                  # slowly biodegradable substrate half saturation coefficient for hydrolysis = 0.03 gCOD/gCOD
    mu_A = 0.8,                  # autotrophic maximum specific growth rate = 0.8 d^(-1)
    K_NH = 1.0,                  # ammonium (nutrient) half saturation coefficient = 1.0 gN/m3
    K_O_A = 0.4,                 # O2 half saturation coefficient for autotrophic growth = 0.4 gO2/m3
    k_a = 0.08                   # ammonification rate constant = 0.08 d^(-1)/(gCOD/m3)
    )

# ASM1 typical values at 10 degree C
# asm1.set_parameters(
#     Y_A = 0.24,                  # autotrophic yield = 0.24 gCOD/gN
#     Y_H = 0.67,                  # heterotrophic yield = 0.67 gCOD/gCOD
#     f_P = 0.08,                  # fraction of biomass yielding particulate products = 0.08, unitless 
#     mu_H = 3,                    # heterotrophic maximum specific growth rate = 3.0 d^(-1)
#     K_S = 20,                    # readily biodegradable substrate half saturation coefficient = 20.0 gCOD/m3
#     K_O_H = 0.2,                 # O2 half saturation coefficient = 0.2 gO2/m3
#     K_NO = 0.5,                  # nitrate half saturation coefficient = 0.5 gN/m3
#     b_H = 0.20,                  # heterotrophic biomass decay rate constant = 0.20 d^(-1)
#     eta_g = 0.8,                 # reduction factor for anoxic growth of heterotrophs = 0.8, unitless
#     eta_h = 0.4,                 # anoxic hydrolysis reduction factor = 0.4, unitless
#     k_h = 1.0,                   # hydrolysis rate constant = 1.0 d^(-1)
#     K_X = 0.01,                  # slowly biodegradable substrate half saturation coefficient for hydrolysis = 0.01 gCOD/gCOD
#     mu_A = 0.3,                  # autotrophic maximum specific growth rate = 0.3 d^(-1)
#     K_NH = 1.0,                  # ammonium (nutrient) half saturation coefficient = 1.0 gN/m3
#     K_O_A = 0.4,                 # O2 half saturation coefficient for autotrophic growth = 0.4 gO2/m3
#     k_a = 0.04                   # ammonification rate constant = 0.04 d^(-1)/(gCOD/m3)
#     )