# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import os
from thermosteam.utils import chemicals_user
from thermosteam import settings
from qsdsan import Components, Process, Processes, _pk
from ..utils import data_path, save_pickle, load_pickle

__all__ = ('load_asm2d_cmps', 'ASM2d')

_path = data_path + 'process_data/_asm2d.tsv'
_path_cmps = os.path.join(data_path, '_asm2d_cmps.pckl')
_load_components = settings.get_default_chemicals

############# Components with default notation #############
def _create_asm2d_cmps(pickle=False):
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

    X_I.i_mass = 0.75
    X_S.i_mass = 0.75
    X_H.i_mass = X_PAO.i_mass = X_AUT.i_mass = 0.9

    cmps_asm2d = Components([S_O2, S_F, S_A, S_NH4, S_NO3, S_PO4, S_I, S_ALK, S_N2,
                             X_I, X_S, X_H, X_PAO, X_PP, X_PHA, X_AUT, X_MeOH, X_MeP,
                             cmps.H2O])
    cmps_asm2d.compile()

    if pickle:
        save_pickle(cmps_asm2d, _path_cmps)
    return cmps_asm2d


# _pickle_asm2d_cmps()

def load_asm2d_cmps():
    if _pk:
        return load_pickle(_path_cmps)
    else:
        return _create_asm2d_cmps(pickle=False)


############ Processes in ASM2d #################
# params = ('f_SI', 'Y_H', 'f_XI_H', 'Y_PAO', 'Y_PO4', 'Y_PHA', 'f_XI_PAO', 'Y_A', 'f_XI_AUT',
#           'K_h', 'eta_NO3', 'eta_fe', 'K_O2', 'K_NO3', 'K_X',
#           'mu_H', 'q_fe', 'eta_NO3_H', 'b_H', 'K_O2_H', 'K_F', 'K_fe', 'K_A_H',
#           'K_NO3_H', 'K_NH4_H', 'K_P_H', 'K_ALK_H',
#           'q_PHA', 'q_PP', 'mu_PAO', 'eta_NO3_PAO', 'b_PAO', 'b_PP', 'b_PHA',
#           'K_O2_PAO', 'K_NO3_PAO', 'K_A_PAO', 'K_NH4_PAO', 'K_PS',' K_P_PAO',
#           'K_ALK_PAO', 'K_PP', 'K_MAX', 'K_IPP', 'K_PHA',
#           'mu_AUT', 'b_AUT', 'K_O2_AUT', 'K_NH4_AUT', 'K_ALK_AUT', 'K_P_AUT',
#           'k_PRE', 'k_RED', 'K_ALK_PRE')

# asm2d = Processes.load_from_file(data_path,
#                                  conserved_for=('COD', 'N', 'P', 'charge'),
#                                  parameters=params,
#                                  compile=False)

# p12 = Process('anox_storage_PP',
#               'S_PO4 + [Y_PHA]X_PHA + [?]S_NO3 -> X_PP + [?]S_N2 + [?]S_NH4 + [?]S_ALK',
#               ref_component='X_PP',
#               rate_equation='q_PP * S_O2/(K_O2_PAO+S_O2) * S_PO4/(K_PS+S_PO4) * S_ALK/(K_ALK_PAO+S_ALK) * (X_PHA/X_PAO)/(K_PHA+X_PHA/X_PAO) * (K_MAX-X_PP/X_PAO)/(K_IPP+K_MAX-X_PP/X_PAO) * X_PAO * eta_NO3_PAO * K_O2_PAO/S_O2 * S_NO3/(K_NO3_PAO+S_NO3)',
#               parameters=('Y_PHA', 'q_PP', 'K_O2_PAO', 'K_PS', 'K_ALK_PAO', 'K_PHA', 'eta_NO3_PAO', 'K_IPP', 'K_NO3_PAO'),
#               conserved_for=('COD', 'N', 'P', 'NOD', 'charge'))

# p14 = Process('PAO_anox_growth',
#               '[1/Y_PAO]X_PHA + [?]S_NO3 + [?]S_PO4 -> X_PAO + [?]S_N2 + [?]S_NH4  + [?]S_ALK',
#               ref_component='X_PAO',
#               rate_equation='mu_PAO * S_O2/(K_O2_PAO + S_O2) * S_NH4/(K_NH4_PAO + S_NH4) * S_PO4/(K_P_PAO + S_PO4) * S_ALK/(K_ALK_PAO + S_ALK) * (X_PHA/X_PAO)/(K_PHA + X_PHA/X_PAO) * X_PAO * eta_NO3_PAO * K_O2_PAO/S_O2 * S_NO3/(K_NO3_PAO + S_NO3)',
#               parameters=('Y_PAO', 'mu_PAO', 'K_O2_PAO', 'K_NH4_PAO', 'K_P_PAO', 'K_ALK_PAO', 'K_PHA', 'eta_NO3_PAO', 'K_NO3_PAO'),
#               conserved_for=('COD', 'N', 'P', 'NOD', 'charge'))

# asm2d.extend([p12, p14])
# asm2d.compile()

# # ASM2d typical values at 20 degree C
# asm2d.set_parameters(
#     f_SI = 0,                    # production of soluble inerts in hydrolysis = 0.0 gCOD/gCOD
#     Y_H = 0.625,                 # heterotrophic yield = 0.625 gCOD/gCOD
#     f_XI_H=0.1,                  # fraction of inert COD generated in heterotrophic biomass lysis = 0.1 gCOD/gCOD
#     Y_PAO = 0.625,               # PAO yield = 0.625 gCOD/gCOD
#     Y_PO4 = 0.4,                 # PP requirement (PO4 release) per PHA stored = 0.4 gP/gCOD
#     Y_PHA = 0.2,                 # PHA requirement for PP storage = 0.2 gCOD/gP
#     f_XI_PAO=0.1,                # fraction of inert COD generated in PAO biomass lysis = 0.1 gCOD/gCOD
#     Y_A = 0.24,                  # autotrophic yield = 0.24 gCOD/gN
#     f_XI_AUT=0.1,                # fraction of inert COD generated in autotrophic biomass lysis = 0.1 gCOD/gCOD
#     K_h = 3,                     # hydrolysis rate constant = 3.0 d^(-1)
#     eta_NO3 = 0.6,               # reduction factor for anoxic hydrolysis = 0.6
#     eta_fe = 0.4,                # anaerobic hydrolysis reduction factor = 0.4
#     K_O2 = 0.2,                  # O2 half saturation coefficient of hydrolysis = 0.2 mgO2/L
#     K_NO3 = 0.5,                 # nitrate half saturation coefficient of hydrolysis = 0.5 mgN/L
#     K_X = 0.1,                   # slowly biodegradable substrate half saturation coefficient for hydrolysis = 0.1 gCOD/gCOD
#     mu_H = 6,                    # heterotrophic maximum specific growth rate = 6.0 d^(-1)
#     q_fe = 3,                    # fermentation maximum rate = 3.0 d^(-1)
#     eta_NO3_H = 0.8,             # denitrification reduction factor for heterotrophic growth = 0.8
#     b_H = 0.4,                   # lysis and decay rate constant = 0.4 d^(-1)
#     K_O2_H=0.2,                  # O2 half saturation coefficient of heterotrophs = 0.2 mgO2/L
#     K_F = 4,                     # fermentable substrate half saturation coefficient for heterotrophic growth = 4.0 mgCOD/L
#     K_fe = 4,                    # fermentable substrate half saturation coefficient for fermentation = 4.0 mgCOD/L
#     K_A_H = 4,                   # VFA half saturation coefficient for heterotrophs = 4.0 mgCOD/L
#     K_NO3_H=0.5,                 # nitrate half saturation coefficient = 0.5 mgN/L
#     K_NH4_H = 0.05,              # ammonium (nutrient) half saturation coefficient for heterotrophs = 0.05 mgN/L
#     K_P_H = 0.01,                # phosphorus (nutrient) half saturation coefficient for heterotrophs = 0.01 mgP/L
#     K_ALK_H = 0.1*12,            # alkalinity half saturation coefficient for heterotrophs = 0.1 mol(HCO3-)/m^3 = 1.2 gC/m^3
#     q_PHA = 3,                   # rate constant for storage of PHA = 3.0 d^(-1)
#     q_PP = 1.5,                  # rate constant for storage of PP = 1.5 d^(-1)
#     mu_PAO = 1,                  # PAO maximum specific growth rate = 1.0 d^(-1)
#     eta_NO3_PAO=0.6,             # denitrification reduction factor for PAO growth = 0.8
#     b_PAO = 0.2,                 # PAO lysis rate = 0.2 d^(-1)
#     b_PP = 0.2,                  # PP lysis rate = 0.2 d^(-1)
#     b_PHA = 0.2,                 # PHA lysis rate = 0.2 d^(-1)
#     K_O2_PAO=0.2,                # O2 half saturation coefficient for PAOs = 0.2 mgO2/L
#     K_NO3_PAO=0.5,               # nitrate half saturation coefficient for PAOs = 0.5 mgN/L
#     K_A_PAO=4.0,                 # VFA half saturation coefficient for PAOs = 4.0 mgCOD/L
#     K_NH4_PAO=0.05,              # ammonium (nutrient) half saturation coefficient for PAOs = 0.05 mgN/L
#     K_PS = 0.2,                  # phosphorus half saturation coefficient for storage of PP = 0.2 mgP/L
#     K_P_PAO=0.01,                # phosphorus (nutrient) half saturation coefficient for heterotrophs = 0.01 mgP/L
#     K_ALK_PAO=0.1*12,            # alkalinity half saturation coefficient for PAOs = 0.1 mol(HCO3-)/m^3 = 1.2 gC/m^3
#     K_PP = 0.01,                 # PP half saturation coefficient for storage of PHA = 0.01 gP/gCOD (?). gCOD/gCOD in GPS-X
#     K_MAX = 0.34,                # maximum ratio of X_PP/X_PAO = 0.34 gX_PP/gX_PAO
#     K_IPP = 0.02,                # inhibition coefficient for PP storage = 0.02 gP/gCOD
#     K_PHA = 0.01,                # PHA half saturation coefficient = 0.01 gCOD/gCOD
#     mu_AUT = 1,                  # autotrophic maximum specific growth rate = 1.0 d^(-1)
#     b_AUT = 0.15,                # autotrophic decay rate = 0.15 d^(-1)
#     K_O2_AUT = 0.5,              # O2 half saturation coefficient for autotrophic growth = 0.5 mgO2/L
#     K_NH4_AUT = 1,               # ammonium (substrate) half saturation coefficient for autotrophic growth = 1.0 mgN/L
#     K_ALK_AUT=0.5*12,            # alkalinity half saturation coefficient for autotrophic growth = 0.5 mol(HCO3-)/m^3 = 6.0 gC/m^3
#     K_P_AUT=0.01,                # phosphorus (nutrient) half saturation coefficient for autotrophic growth = 0.01 mgP/L
#     k_PRE = 1,                   # phosphorus precipitation with MeOH rate constant = 1.0 m^3/g/d
#     k_RED = 0.6,                 # redissoluation of phosphates rate constant = 0.6 d^(-1)
#     K_ALK_PRE=0.5*12             # alkalinity half saturation coefficient for phosphate precipitation = 0.5 mol(HCO3-)/m^3 = 6.0 gC/m^3
#     )

# ASM2d typical values at 10 degree C
# asm2d.set_parameters(
#     f_SI = 0,                    # production of soluble inerts in hydrolysis = 0.0 gCOD/gCOD
#     Y_H = 0.625,                 # heterotrophic yield = 0.625 gCOD/gCOD
#     f_XI_H=0.1,                  # fraction of inert COD generated in heterotrophic biomass lysis = 0.1 gCOD/gCOD
#     Y_PAO = 0.625,               # PAO yield = 0.625 gCOD/gCOD
#     Y_PO4 = 0.4,                 # PP requirement (PO4 release) per PHA stored = 0.4 gP/gCOD
#     Y_PHA = 0.2,                 # PHA requirement for PP storage = 0.2 gCOD/gP
#     f_XI_PAO=0.1,                # fraction of inert COD generated in PAO biomass lysis = 0.1 gCOD/gCOD
#     Y_A = 0.24,                  # autotrophic yield = 0.24 gCOD/gN
#     f_XI_AUT=0.1,                # fraction of inert COD generated in autotrophic biomass lysis = 0.1 gCOD/gCOD
#     K_h = 2,                     # hydrolysis rate constant = 2.0 d^(-1)
#     eta_NO3 = 0.6,               # reduction factor for anoxic hydrolysis = 0.6
#     eta_fe = 0.4,                # anaerobic hydrolysis reduction factor = 0.4
#     K_O2 = 0.2,                  # O2 half saturation coefficient of hydrolysis = 0.2 mgO2/L
#     K_NO3 = 0.5,                 # nitrate half saturation coefficient of hydrolysis = 0.5 mgN/L
#     K_X = 0.1,                   # slowly biodegradable substrate half saturation coefficient for hydrolysis = 0.1 gCOD/gCOD
#     mu_H = 3,                    # heterotrophic maximum specific growth rate = 3.0 d^(-1)
#     q_fe = 1.5,                  # fermentation maximum rate = 1.5 d^(-1)
#     eta_NO3_H = 0.8,             # denitrification reduction factor for heterotrophic growth = 0.8
#     b_H = 0.4,                   # lysis and decay rate constant = 0.4 d^(-1)
#     K_O2_H=0.2,                  # O2 half saturation coefficient of heterotrophs = 0.2 mgO2/L
#     K_F = 4,                     # fermentable substrate half saturation coefficient for heterotrophic growth = 4.0 mgCOD/L
#     K_fe = 4,                    # fermentable substrate half saturation coefficient for fermentation = 4.0 mgCOD/L
#     K_A_H = 4,                   # VFA half saturation coefficient for heterotrophs = 4.0 mgCOD/L
#     K_NO3_H=0.5,                 # nitrate half saturation coefficient = 0.5 mgN/L
#     K_NH4_H = 0.05,              # ammonium (nutrient) half saturation coefficient for heterotrophs = 0.05 mgN/L
#     K_P_H = 0.01,                # phosphorus (nutrient) half saturation coefficient for heterotrophs = 0.01 mgP/L
#     K_ALK_H = 0.1*12,            # alkalinity half saturation coefficient for heterotrophs = 0.1 mol(HCO3-)/m^3 = 1.2 gC/m^3
#     q_PHA = 2,                   # rate constant for storage of PHA = 2.0 d^(-1)
#     q_PP = 1.0,                  # rate constant for storage of PP = 1.0 d^(-1)
#     mu_PAO = 0.67,               # PAO maximum specific growth rate = 0.67 d^(-1)
#     eta_NO3_PAO=0.6,             # denitrification reduction factor for PAO growth = 0.8
#     b_PAO = 0.1,                 # PAO lysis rate = 0.1 d^(-1)
#     b_PP = 0.1,                  # PP lysis rate = 0.1 d^(-1)
#     b_PHA = 0.1,                 # PHA lysis rate = 0.1 d^(-1)
#     K_O2_PAO=0.2,                # O2 half saturation coefficient for PAOs = 0.2 mgO2/L
#     K_NO3_PAO=0.5,               # nitrate half saturation coefficient for PAOs = 0.5 mgN/L
#     K_A_PAO=4.0,                 # VFA half saturation coefficient for PAOs = 4.0 mgCOD/L
#     K_NH4_PAO=0.05,              # ammonium (nutrient) half saturation coefficient for PAOs = 0.05 mgN/L
#     K_PS = 0.2,                  # phosphorus half saturation coefficient for storage of PP = 0.2 mgP/L
#     K_P_PAO=0.01,                # phosphorus (nutrient) half saturation coefficient for heterotrophs = 0.01 mgP/L
#     K_ALK_PAO=0.1*12,            # alkalinity half saturation coefficient for PAOs = 0.1 mol(HCO3-)/m^3 = 1.2 gC/m^3
#     K_PP = 0.01,                 # PP half saturation coefficient for storage of PHA = 0.01 gP/gCOD (?). gCOD/gCOD in GPS-X
#     K_MAX = 0.34,                # maximum ratio of X_PP/X_PAO = 0.34 gX_PP/gX_PAO
#     K_IPP = 0.02,                # inhibition coefficient for PP storage = 0.02 gP/gCOD
#     K_PHA = 0.01,                # PHA half saturation coefficient = 0.01 gCOD/gCOD
#     mu_AUT = 0.35,               # autotrophic maximum specific growth rate = 0.35 d^(-1)
#     b_AUT = 0.05,                # autotrophic decay rate = 0.05 d^(-1)
#     K_O2_AUT = 0.5,              # O2 half saturation coefficient for autotrophic growth = 0.5 mgO2/L
#     K_NH4_AUT = 1,               # ammonium (substrate) half saturation coefficient for autotrophic growth = 1.0 mgN/L
#     K_ALK_AUT=0.5*12,            # alkalinity half saturation coefficient for autotrophic growth = 0.5 mol(HCO3-)/m^3 = 6.0 gC/m^3
#     K_P_AUT=0.01,                # phosphorus (nutrient) half saturation coefficient for autotrophic growth = 0.01 mgP/L
#     k_PRE = 1,                   # phosphorus precipitation with MeOH rate constant = 1.0 m^3/g/d
#     k_RED = 0.6,                 # redissoluation of phosphates rate constant = 0.6 d^(-1)
#     K_ALK_PRE=0.5*12             # alkalinity half saturation coefficient for phosphate precipitation = 0.5 mol(HCO3-)/m^3 = 6.0 gC/m^3
#     )

@chemicals_user
class ASM2d(Processes):
    '''
    Activated Sludge Model No. 2d in orginal notation. [1]_, [2]_

    Parameters
    ----------
    components: class:`CompiledComponents`, optional
        Components corresponding to each entry in the stoichiometry array,
        defaults to thermosteam.settings.chemicals.
    iN_SI : float, optional
        Nitrogen content of inert soluble COD, in [g N/g COD]. The default is 0.01.
    iN_SF : float, optional
        Nitrogen content of fermentable substrate, in [g N/g COD]. The default is 0.03.
    iN_XI : float, optional
        Nitrogen content of inert particulate COD, in [g N/g COD]. The default is 0.02.
    iN_XS : float, optional
        Nitrogen content of slowly biodegradable substrate, in [g N/g COD]. The default
        is 0.04.
    iN_BM : float, optional
        Nitrogen content of biomass, in [g N/g COD]. The default is 0.07.
    iP_SI : float, optional
        Phosphorus content of inert soluble COD, in [g P/g COD]. The default is 0.0.
    iP_SF : float, optional
        Phosphorus content of fermentable substrate, in [g P/g COD]. The default
        is 0.01.
    iP_XI : float, optional
        Phosphorus content of inert particulate COD, in [g P/g COD]. The default
        is 0.01.
    iP_XS : float, optional
        Phosphorus content of slowly biodegradable substrate, in [g P/g COD]. The
        default is 0.01.
    iP_BM : float, optional
        Phosphorus content of biomass, in [g P/g COD]. The default is 0.02.
    iTSS_XI : float, optional
        TSS to COD ratio for inert particulate COD, in [g TSS/g COD]. The default
        is 0.75.
    iTSS_XS : float, optional
        TSS to COD ratio for slowly biodegradable substrate, in [g TSS/g COD]. The
        default is 0.75.
    iTSS_BM : float, optional
        TSS to COD ratio for biomass, in [g TSS/g COD]. The default is 0.9.
    f_SI : float, optional
        Production of soluble inerts in hydrolysis, in [g COD/g COD]. The default
        is 0.0.
    Y_H : float, optional
        Heterotrophic yield coefficient, in[g COD/g COD]. The default is 0.625.
    f_XI_H : float, optional
        Fraction of inert COD generated in heterotrophic biomass lysis,
        in [g COD/g COD]. The default is 0.1.
    Y_PAO : float, optional
        PAO yield coefficient, in[g COD/g COD]. The default is 0.625.
    Y_PO4 : float, optional
        PP requirement (PO4 release) per PHA stored, in [g P/g COD].
        The default is 0.4.
    Y_PHA : float, optional
        PHA requirement for PP storage, in [g COD/g P]. The default is 0.2.
    f_XI_PAO : float, optional
        Fraction of inert COD generated in PAO lysis, in [g COD/g COD].
        The default is 0.1.
    Y_A : float, optional
        Autotrophic yield, in [g COD/g N]. The default is 0.24.
    f_XI_AUT : float, optional
        Fraction of inert COD generated in autotrophic biomass lysis,
        in [g COD/g COD]. The default is 0.1.
    K_h : float, optional
        Hydrolysis rate constant, in [d^(-1)]. The default is 3.0.
    eta_NO3 : float, optional
        Reduction factor for anoxic hydrolysis, dimensionless. The default is 0.6.
    eta_fe : float, optional
        Anaerobic hydrolysis reduction factor, dimensionless. The default is 0.4.
    K_O2 : float, optional
        O2 half saturation coefficient for hydrolysis, in [g O2/m^3]. The default is 0.2.
    K_NO3 : float, optional
        Nitrate half saturation coefficient for hydrolysis, in [g N/m^3].
        The default is 0.5.
    K_X : float, optional
        Slowly biodegradable substrate half saturation coefficient for hydrolysis,
        in [g COD/g COD]. The default is 0.1.
    mu_H : float, optional
        Heterotrophic maximum specific growth rate, in [d^(-1)]. The default is 6.0.
    q_fe : float, optional
        Fermentation maximum rate, in [d^(-1)]. The default is 3.0.
    eta_NO3_H : float, optional
        Reduction factor for anoxic heterotrophic growth, dimensionless.
        The default is 0.8.
    b_H : float, optional
        Lysis and decay rate constant, in [d^(-1)]. The default is 0.4.
    K_O2_H : float, optional
        O2 half saturation coefficient for heterotrophs, in [g O2/m^3].
        The default is 0.2.
    K_F : float, optional
        Fermentable substrate half saturation coefficient for heterotrophic growth,
        in [g COD/m^3]. The default is 4.0.
    K_fe : float, optional
        Fermentable substrate half saturation coefficient for fermentation,
        in [g COD/m^3]. The default is 4.0.
    K_A_H : float, optional
        VFA half saturation coefficient for heterotrophs, in [g COD/m^3].
        The default is 4.0.
    K_NO3_H : float, optional
        Nitrate half saturation coefficient for heterotrophs, in [g N/m^3].
        The default is 0.5.
    K_NH4_H : float, optional
        Ammonium (nutrient) half saturation coefficient for heterotrophs, in [g N/m^3].
        The default is 0.05.
    K_P_H : float, optional
        Phosphorus (nutrient) half saturation coefficient for heterotrophs,
        in [g P/m^3]. The default is 0.01.
    K_ALK_H : float, optional
        Alkalinity half saturation coefficient for heterotrophs, in [mol(HCO3-)/m^3].
        The default is 0.1.
    q_PHA : float, optional
        Rate constant for storage of PHA, in [d^(-1)]. The default is 3.0.
    q_PP : float, optional
        Rate constant for storage of PP, in [d^(-1)]. The default is 1.5.
    mu_PAO : float, optional
        PAO maximum specific growth rate, in [d^(-1)]. The default is 1.0.
    eta_NO3_PAO : float, optional
        Reduction factor for anoxic growth of PAO, dimensionless. The default is 0.6.
    b_PAO : float, optional
        PAO lysis rate constant, in [d^(-1)]. The default is 0.2.
    b_PP : float, optional
        PP lysis rate constant, in [d^(-1)]. The default is 0.2.
    b_PHA : float, optional
        PHA lysis rate constant, in [d^(-1)]. The default is 0.2.
    K_O2_PAO : float, optional
        O2 half saturation coefficient for PAO, in [g O2/m^3]. The default is 0.2.
    K_NO3_PAO : float, optional
        Nitrate half saturation coefficient for PAO, in [g N/m^3]. The default is 0.5.
    K_A_PAO : float, optional
        VFA half saturation coefficient for PAO, in [g COD/m^3]. The default is 4.0.
    K_NH4_PAO : float, optional
        Ammonium (nutrient) half saturation coefficient for PAO, in [g N/m^3].
        The default is 0.05.
    K_PS : float, optional
        Phosphorus half saturation coefficient for storage of PP, in [g P/m^3].
        The default is 0.2.
    K_P_PAO : float, optional
        Phosphorus (nutrient) half saturation coefficient for PAO, in [g P/m^3].
        The default is 0.01.
    K_ALK_PAO : float, optional
        Alkalinity half saturation coefficient for PAO, in [mol(HCO3-)/m^3].
        The default is 0.1.
    K_PP : float, optional
        Poly-phosphate half saturation coefficient for storage of PHA, in [g P/g COD].
        The default is 0.01.
    K_MAX : float, optional
        Maximum ratio of poly-phosphate to PAO, in [g P/g COD]. The default is 0.34.
    K_IPP : float, optional
        Inhibition coefficient for poly-phosphate storage, in [g P/g COD].
        The default is 0.02.
    K_PHA : float, optional
        PHA half saturation coefficient, in [g COD/g COD]. The default is 0.01.
    mu_AUT : float, optional
        Autotrophic maximum specific growth rate, in [d^(-1)]. The default is 1.0.
    b_AUT : float, optional
        Autotrophic decay rate, in [d^(-1)]. The default is 0.15.
    K_O2_AUT : float, optional
        O2 half saturation coefficient for autotrophs, in [g O2/m^3].
        The default is 0.5.
    K_NH4_AUT : float, optional
        Ammonium (nutrient) half saturation coefficient for autotrophs, in [g N/m^3].
        The default is 1.0.
    K_ALK_AUT : float, optional
        Alkalinity half saturation coefficient for autotrophs, in [mol(HCO3-)/m^3].
        The default is 0.5.
    K_P_AUT : float, optional
        Phosphorus (nutrient) half saturation coefficient for autotrophs, in [g P/m^3].
        The default is 0.01.
    k_PRE : float, optional
        Rate constant for phosphorus precipitation with Fe(OH)3, in [m^3/g Fe(OH)3/d].
        The default is 1.0.
    k_RED : float, optional
        Rate constant for redissoluation of FePO4, in [d^(-1)]. The default is 0.6.
    K_ALK_PRE : float, optional
        Alkalinity half saturation coefficient for P precipitation, in [mol(HCO3-)/m^3].
        The default is 0.5.
    path : str, optional
        Alternative file path for the Gujer matrix. The default is None.


    References
    ----------
    .. [1] Henze, M.; Gujer, W.; Mino, T.; Loosdrecht, M. van. Activated Sludge
        Models: ASM1, ASM2, ASM2d and ASM3; IWA task group on mathematical modelling
        for design and operation of biological wastewater treatment, Ed.; IWA
        Publishing: London, 2000.
    .. [2] Rieger, L.; Gillot, S.; Langergraber, G.; Ohtsuki, T.; Shaw, A.; Takács,
        I.; Winkler, S. Guidelines for Using Activated Sludge Models; IWA Publishing:
        London, New York, 2012; Vol. 11.
        https://doi.org/10.2166/9781780401164.
    '''
    _params = ('f_SI', 'Y_H', 'f_XI_H', 'Y_PAO', 'Y_PO4', 'Y_PHA', 'f_XI_PAO', 'Y_A', 'f_XI_AUT',
               'K_h', 'eta_NO3', 'eta_fe', 'K_O2', 'K_NO3', 'K_X',
               'mu_H', 'q_fe', 'eta_NO3_H', 'b_H', 'K_O2_H', 'K_F', 'K_fe', 'K_A_H',
               'K_NO3_H', 'K_NH4_H', 'K_P_H', 'K_ALK_H',
               'q_PHA', 'q_PP', 'mu_PAO', 'eta_NO3_PAO', 'b_PAO', 'b_PP', 'b_PHA',
               'K_O2_PAO', 'K_NO3_PAO', 'K_A_PAO', 'K_NH4_PAO', 'K_PS',' K_P_PAO',
               'K_ALK_PAO', 'K_PP', 'K_MAX', 'K_IPP', 'K_PHA',
               'mu_AUT', 'b_AUT', 'K_O2_AUT', 'K_NH4_AUT', 'K_ALK_AUT', 'K_P_AUT',
               'k_PRE', 'k_RED', 'K_ALK_PRE')



    def __new__(cls, components=None,
                iN_SI=0.01, iN_SF=0.03, iN_XI=0.02, iN_XS=0.04, iN_BM=0.07,
                iP_SI=0.0, iP_SF=0.01, iP_XI=0.01, iP_XS=0.01, iP_BM=0.02,
                iTSS_XI=0.75, iTSS_XS=0.75, iTSS_BM=0.9,
                f_SI=0.0, Y_H=0.625, f_XI_H=0.1,
                Y_PAO=0.625, Y_PO4=0.4, Y_PHA=0.2, f_XI_PAO=0.1,
                Y_A=0.24, f_XI_AUT=0.1,
                K_h=3.0, eta_NO3=0.6, eta_fe=0.4, K_O2=0.2, K_NO3=0.5, K_X=0.1,
                mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4, K_O2_H=0.2, K_F=4.0,
                K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_P_H=0.01, K_ALK_H=0.1,
                q_PHA=3.0, q_PP=1.5, mu_PAO=1.0, eta_NO3_PAO=0.6, b_PAO=0.2, b_PP=0.2,
                b_PHA=0.2, K_O2_PAO=0.2, K_NO3_PAO=0.5, K_A_PAO=4.0, K_NH4_PAO=0.05,
                K_PS=0.2, K_P_PAO=0.01, K_ALK_PAO=0.1,
                K_PP=0.01, K_MAX=0.34, K_IPP=0.02, K_PHA=0.01,
                mu_AUT=1.0, b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,
                k_PRE=1.0, k_RED=0.6, K_ALK_PRE=0.5,
                path=None, **kwargs):

        if not path: path = _path
        
        cmps = _load_components(components)
        cmps.S_I.i_N = iN_SI
        cmps.S_F.i_N = iN_SF
        cmps.X_I.i_N = iN_XI
        cmps.X_S.i_N = iN_XS
        cmps.X_H.i_N = cmps.X_PAO.i_N = cmps.X_AUT.i_N = iN_BM
        cmps.S_I.i_P = iP_SI
        cmps.S_F.i_P = iP_SF
        cmps.X_I.i_P = iP_XI
        cmps.X_S.i_P = iP_XS
        cmps.X_H.i_P = cmps.X_PAO.i_P = cmps.X_AUT.i_P = iP_BM
        cmps.X_I.i_mass = iTSS_XI
        cmps.X_S.i_mass = iTSS_XS
        cmps.X_H.i_mass = cmps.X_PAO.i_mass = cmps.X_AUT.i_mass = iTSS_BM

        self = Processes.load_from_file(path,
                                        components=cmps,
                                        conserved_for=('COD', 'N', 'P', 'charge'),
                                        parameters=cls._params,
                                        compile=False)

        if path == None:
            _p12 = Process('anox_storage_PP',
                           'S_PO4 + [Y_PHA]X_PHA + [?]S_NO3 -> X_PP + [?]S_N2 + [?]S_NH4 + [?]S_ALK',
                           components=cmps,
                           ref_component='X_PP',
                           rate_equation='q_PP * S_O2/(K_O2_PAO+S_O2) * S_PO4/(K_PS+S_PO4) * S_ALK/(K_ALK_PAO+S_ALK) * (X_PHA/X_PAO)/(K_PHA+X_PHA/X_PAO) * (K_MAX-X_PP/X_PAO)/(K_IPP+K_MAX-X_PP/X_PAO) * X_PAO * eta_NO3_PAO * K_O2_PAO/S_O2 * S_NO3/(K_NO3_PAO+S_NO3)',
                           parameters=('Y_PHA', 'q_PP', 'K_O2_PAO', 'K_PS', 'K_ALK_PAO', 'K_PHA', 'eta_NO3_PAO', 'K_IPP', 'K_NO3_PAO'),
                           conserved_for=('COD', 'N', 'P', 'NOD', 'charge'))

            _p14 = Process('PAO_anox_growth',
                           '[1/Y_PAO]X_PHA + [?]S_NO3 + [?]S_PO4 -> X_PAO + [?]S_N2 + [?]S_NH4  + [?]S_ALK',
                           components=cmps,
                           ref_component='X_PAO',
                           rate_equation='mu_PAO * S_O2/(K_O2_PAO + S_O2) * S_NH4/(K_NH4_PAO + S_NH4) * S_PO4/(K_P_PAO + S_PO4) * S_ALK/(K_ALK_PAO + S_ALK) * (X_PHA/X_PAO)/(K_PHA + X_PHA/X_PAO) * X_PAO * eta_NO3_PAO * K_O2_PAO/S_O2 * S_NO3/(K_NO3_PAO + S_NO3)',
                           parameters=('Y_PAO', 'mu_PAO', 'K_O2_PAO', 'K_NH4_PAO', 'K_P_PAO', 'K_ALK_PAO', 'K_PHA', 'eta_NO3_PAO', 'K_NO3_PAO'),
                           conserved_for=('COD', 'N', 'P', 'NOD', 'charge'))
            self.extend([_p12, _p14])

        self.compile()
        self.set_parameters(f_SI=f_SI, Y_H=Y_H, f_XI_H=f_XI_H, Y_PAO=Y_PAO, Y_PO4=Y_PO4,
                            Y_PHA=Y_PHA, f_XI_PAO=f_XI_PAO, Y_A=Y_A, f_XI_AUT=f_XI_AUT,
                            K_h=K_h, eta_NO3=eta_NO3, eta_fe=eta_fe, K_O2=K_O2,
                            K_NO3=K_NO3, K_X=K_X,
                            mu_H=mu_H, q_fe=q_fe, eta_NO3_H=eta_NO3_H, b_H=b_H,
                            K_O2_H=K_O2_H, K_F=K_F, K_fe=K_fe, K_A_H=K_A_H, K_NO3_H=K_NO3_H,
                            K_NH4_H=K_NH4_H, K_P_H=K_P_H, K_ALK_H=K_ALK_H*12,
                            q_PHA=q_PHA, q_PP=q_PP, mu_PAO=mu_PAO, eta_NO3_PAO=eta_NO3_PAO,
                            b_PAO=b_PAO, b_PP=b_PP, b_PHA=b_PHA, K_O2_PAO=K_O2_PAO,
                            K_NO3_PAO=K_NO3_PAO, K_A_PAO=K_A_PAO, K_NH4_PAO=K_NH4_PAO,
                            K_PS=K_PS, K_P_PAO=K_P_PAO, K_ALK_PAO=K_ALK_PAO*12, K_PP=K_PP,
                            K_MAX=K_MAX, K_IPP=K_IPP, K_PHA=K_PHA,
                            mu_AUT=mu_AUT, b_AUT=b_AUT, K_O2_AUT=K_O2_AUT, K_NH4_AUT=K_NH4_AUT,
                            K_ALK_AUT=K_ALK_AUT*12, K_P_AUT=K_P_AUT,
                            k_PRE=k_PRE, k_RED=k_RED, K_ALK_PRE=K_ALK_PRE*12,
                            **kwargs)
        return self