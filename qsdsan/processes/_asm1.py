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
from qsdsan import Components, Processes, _pk
from ..utils import data_path, save_pickle, load_pickle

__all__ = ('load_asm1_cmps', 'ASM1')

_path = data_path + 'process_data/_asm1.tsv'
_path_cmps = os.path.join(data_path, '_asm1_cmps.pckl')
_load_components = settings.get_default_chemicals

############# Components with default notation #############
def _create_asm1_cmps(pickle=False):
    cmps = Components.load_default()

    S_I = cmps.S_U_Inf.copy('S_I')
    S_I.description = 'Soluble inert organic matter'

    X_I = cmps.X_U_Inf.copy('X_I')
    X_I.description = 'Particulate inert organic matter'


    S_S = cmps.S_F.copy('S_S')
    S_S.description = 'Readily biodegradable substrate'
    S_S.i_N = S_I.i_N = 0

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
    X_P.i_N = X_I.i_N = 0.06                  # i_XP

    # X_I.i_mass = X_S.i_mass = X_BH.i_mass = X_BA.i_mass = X_P.i_mass = .75    # fr_COD_SS

    S_O = cmps.S_O2.copy('S_O')

    S_NO = cmps.S_NO3.copy('S_NO')
    S_NO.description = 'Nitrate and nitrite nitrogen'

    S_NH = cmps.S_NH4.copy('S_NH')

    S_ND = cmps.S_F.copy('S_ND')
    S_ND.description = 'Soluble biodegradable organic nitrogen'
    S_ND.measured_as = 'N'
    S_ND.i_COD = S_ND.i_C = S_ND.i_P = S_ND.i_mass = 0

    X_ND = cmps.X_B_Subst.copy('X_ND')
    X_ND.description = 'Particulate biodegradable organic nitrogen'
    X_ND.measured_as = 'N'
    X_ND.i_COD = X_ND.i_C = X_ND.i_P = X_ND.i_mass = 0

    S_ALK = cmps.S_CO3.copy('S_ALK')      # measured as g C
    S_ALK.description = 'Alkalinity, assumed to be HCO3-'

    # add S_N2 to close mass balance
    # add water for the creation of WasteStream objects
    cmps_asm1 = Components([S_I, S_S, X_I, X_S, X_BH, X_BA, X_P,
                            S_O, S_NO, S_NH, S_ND, X_ND, S_ALK,
                            cmps.S_N2, cmps.H2O])
    cmps_asm1.compile()

    if pickle:
        save_pickle(cmps_asm1, _path_cmps)
    return cmps_asm1


#_create_asm1_cmps(True)

def load_asm1_cmps():
    if _pk:
        return load_pickle(_path_cmps)
    else:
        return _create_asm1_cmps(pickle=False)



############ Processes in ASM1 #################
# params = ('Y_H', 'Y_A', 'f_P',
#           'mu_H', 'K_S', 'K_O_H', 'K_NO', 'b_H',
#           'mu_A', 'K_NH', 'K_O_A', 'b_A',
#           'eta_g', 'k_a', 'k_h', 'K_X', 'eta_h')

# asm1 = Processes.load_from_file(data_path,
#                                 conserved_for=('COD', 'charge', 'N'),
#                                 parameters=params,
#                                 compile=True)

# ASM1 typical values at 20 degree C
# asm1.set_parameters(
#     Y_A = 0.24,                  # autotrophic yield = 0.24 gCOD/gN
#     Y_H = 0.67,                  # heterotrophic yield = 0.67 gCOD/gCOD
#     f_P = 0.08,                  # fraction of biomass yielding particulate products = 0.08, unitless
#     mu_H = 6,                    # heterotrophic maximum specific growth rate = 6.0 d^(-1)
#     K_S = 20,                    # readily biodegradable substrate half saturation coefficient = 20.0 gCOD/m3
#     K_O_H = 0.2,                 # O2 half saturation coefficient = 0.2 gO2/m3
#     K_NO = 0.5,                  # nitrate half saturation coefficient = 0.5 gN/m3
#     b_H = 0.62,                  # heterotrophic biomass decay rate constant = 0.62 d^(-1)
#     eta_g = 0.8,                 # reduction factor for anoxic growth of heterotrophs = 0.8, unitless
#     eta_h = 0.4,                 # anoxic hydrolysis reduction factor = 0.4, unitless
#     k_h = 3.0,                   # hydrolysis rate constant = 3.0 d^(-1)
#     K_X = 0.03,                  # slowly biodegradable substrate half saturation coefficient for hydrolysis = 0.03 gCOD/gCOD
#     mu_A = 0.8,                  # autotrophic maximum specific growth rate = 0.8 d^(-1)
#     K_NH = 1.0,                  # ammonium (nutrient) half saturation coefficient = 1.0 gN/m3
#     K_O_A = 0.4,                 # O2 half saturation coefficient for autotrophic growth = 0.4 gO2/m3
#     b_A = 0.05,                  # !!! BSM1 value
#     k_a = 0.08                   # ammonification rate constant = 0.08 d^(-1)/(gCOD/m3)
#     )

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
#     b_A = 0.05,                  # !!! BSM1 value
#     k_a = 0.04                   # ammonification rate constant = 0.04 d^(-1)/(gCOD/m3)
#     )

@chemicals_user
class ASM1(Processes):
    '''
    Activated Sludge Model No. 1 in original notation. [1]_, [2]_

    Parameters
    ----------
    components: class:`CompiledComponents`, optional
        Components corresponding to each entry in the stoichiometry array,
        defaults to thermosteam.settings.chemicals.
    Y_A : float, optional
        Autotrophic yield, in [g COD/g N]. The default is 0.24.
    Y_H : float, optional
        Heterotrophic yield, in [g COD/g COD]. The default is 0.67.
    f_P : float, optional
        Fraction of biomass yielding particulate products. The default is 0.08.
    i_XB : float, optional
        Nitrogen content of biomass, in [g N/g COD]. The default is 0.08.
    i_XP : float, optional
        Nitrogen content of particulate products arising from biomass decay,
        in [g N/g COD]. The default is 0.06.
    mu_H : float, optional
        Heterotrophic maximum specific growth rate, in [d^(-1)]. The default
        is 4.0.
    K_S : float, optional
        Readily biodegradable substrate half saturation coefficient, in [g COD/m^3].
        The default is 10.0.
    K_O_H : float, optional
        Oxygen half saturation coefficient for heterotrophic growth, in [g O2/m^3].
        The default is 0.2.
    K_NO : float, optional
        Nitrate half saturation coefficient, in [g N/m^3]. The default is 0.5.
    b_H : float, optional
        Heterotrophic biomass decay rate constant, in [d^(-1)]. The default is 0.3.
    eta_g : float, optional
        Reduction factor for anoxic growth of heterotrophs. The default is 0.8.
    eta_h : float, optional
        Reduction factor for anoxic hydrolysis. The default is 0.8.
    k_h : float, optional
        Hydrolysis rate constant, in [d^(-1)]. The default is 3.0.
    K_X : float, optional
        Slowly biodegradable substrate half saturation coefficient for hydrolysis,
        in [g COD/g COD]. The default is 0.1.
    mu_A : float, optional
        Autotrophic maximum specific growth rate, in [d^(-1)]. The default is 0.5.
    K_NH : float, optional
        Ammonium (nutrient) half saturation coefficient, in [g N/m^3]. The default
        is 1.0.
    b_A : float, optional
        Autotrophic biomass decay rate constant, in [d^(-1)]. The default is 0.05.
    K_O_A : float, optional
        Oxygen half saturation coefficient for autotrophic growth, in [g O2/m^3].
        The default is 0.4.
    k_a : float, optional
        Ammonification rate constant, in [m^3/g COD/d]. The default is 0.05.
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

    _params = ('Y_H', 'Y_A', 'f_P',
              'mu_H', 'K_S', 'K_O_H', 'K_NO', 'b_H',
              'mu_A', 'K_NH', 'K_O_A', 'b_A',
              'eta_g', 'k_a', 'k_h', 'K_X', 'eta_h')

    def __new__(cls, components=None, Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
                mu_H=4.0, K_S=10.0, K_O_H=0.2, K_NO=0.5, b_H=0.3, eta_g=0.8, eta_h=0.8,
                k_h=3.0, K_X=0.1, mu_A=0.5, K_NH=1.0, b_A=0.05, K_O_A=0.4, k_a=0.05,
                fr_SS_COD=0.75, path=None, **kwargs):
        if not path: path = _path
        
        cmps = _load_components(components)
        cmps.X_BH.i_N = cmps.X_BA.i_N = i_XB
        cmps.X_P.i_N = cmps.X_I.i_N = i_XP
        cmps.X_I.i_mass = cmps.X_S.i_mass = cmps.X_P.i_mass = cmps.X_BH.i_mass = cmps.X_BA.i_mass = fr_SS_COD
        cmps.refresh_constants()
        
        self = Processes.load_from_file(path,
                                        conserved_for=('COD', 'charge', 'N'),
                                        parameters=cls._params,
                                        components=cmps,
                                        compile=True)

        self.set_parameters(Y_A=Y_A, Y_H=Y_H, f_P=f_P, mu_H=mu_H, K_S=K_S, K_O_H=K_O_H,
                            K_NO=K_NO, b_H=b_H, eta_g=eta_g, eta_h=eta_h, k_h=k_h,
                            K_X=K_X, mu_A=mu_A, K_NH=K_NH, b_A=b_A, K_O_A=K_O_A, k_a=k_a,
                            **kwargs)
        return self