# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''
import numpy as np
from thermosteam.utils import chemicals_user
from thermosteam import settings
from qsdsan import Component, Components, Processes, CompiledProcesses
from ..utils import ospath, data_path, load_data
from . import Monod, ion_speciation
from scipy.optimize import brenth
# from math import log10


__all__ = ('create_asm2d_cmps', 'ASM2d',
           'create_masm2d_cmps', 'mASM2d')

_path = ospath.join(data_path, 'process_data/_asm2d.tsv')
_load_components = settings.get_default_chemicals

############# Components with default notation #############
def create_asm2d_cmps(set_thermo=True):
    cmps = Components.load_default()

    S_A = cmps.S_Ac.copy('S_A')
    # S_A.i_charge = -1/64
    S_ALK = cmps.S_CO3.copy('S_ALK')      # measured as g C
    S_F = cmps.S_F.copy('S_F')
    S_I = cmps.S_U_E.copy('S_I')
    S_NH4 = cmps.S_NH4.copy('S_NH4')
    # S_NH4.i_charge = 1/14

    X_AUT = cmps.X_AOO.copy('X_AUT')
    X_H = cmps.X_OHO.copy('X_H')
    X_I = cmps.X_U_OHO_E.copy('X_I')
    X_MeOH = cmps.X_FeOH.copy('X_MeOH')
    X_MeP = cmps.X_FePO4.copy('X_MeP')
    X_PAO = cmps.X_PAO.copy('X_PAO')
    X_PHA = cmps.X_PAO_PHA.copy('X_PHA')
    X_PP = cmps.X_PAO_PP_Lo.copy('X_PP')
    X_PP.i_charge = -1/31
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

    cmps_asm2d = Components([cmps.S_O2, cmps.S_N2, S_NH4, cmps.S_NO3, cmps.S_PO4,
                             S_F, S_A, S_I, S_ALK, X_I, X_S, X_H, X_PAO, X_PP,
                             X_PHA, X_AUT, X_MeOH, X_MeP, cmps.H2O])

    cmps_asm2d.compile()
    if set_thermo: settings.set_thermo(cmps_asm2d)

    return cmps_asm2d

def create_masm2d_cmps(set_thermo=True):
    c2d = create_asm2d_cmps(False)
    ion_kwargs = dict(particle_size='Soluble',
                      degradability='Undegradable',
                      organic=False)
    mineral_kwargs = dict(particle_size='Particulate',
                         degradability='Undegradable',
                         organic=False)
    S_K = Component.from_chemical('S_K', chemical='K+', measured_as='K',
                                  description='Potassium', **ion_kwargs)
    
    S_Mg = Component.from_chemical('S_Mg', chemical='Mg+2', measured_as='Mg',
                                   description='Magnesium', **ion_kwargs)

    S_IC = c2d.S_ALK.copy('S_IC')
    c2d.S_PO4.formula = 'HPO4-2'
    c2d.S_PO4.measured_as = 'P'    
    
    c2d.S_F.i_C = c2d.X_S.i_C = 0.31843
    c2d.S_F.i_N = c2d.X_S.i_N = 0.03352
    c2d.S_F.i_P = c2d.X_S.i_P = 5.59e-3
    
    c2d.S_I.i_C = c2d.X_I.i_C = 0.36178
    c2d.S_I.i_N = c2d.X_I.i_N = 0.06003
    c2d.S_I.i_P = c2d.X_I.i_P = 6.49e-3
    c2d.S_I.i_K = c2d.X_I.i_K = 0.0
    c2d.S_F.i_mass = c2d.X_S.i_mass = c2d.S_I.i_mass = c2d.X_I.i_mass = 0.75
    c2d.S_F.f_Vmass_Totmass = c2d.X_S.f_Vmass_Totmass = c2d.S_I.f_Vmass_Totmass = c2d.X_I.f_Vmass_Totmass = 0.85
    
    c2d.X_H.i_C = c2d.X_AUT.i_C = c2d.X_PAO.i_C = 0.36612
    c2d.X_H.i_N = c2d.X_AUT.i_N = c2d.X_PAO.i_N = 0.08615
    c2d.X_H.i_P = c2d.X_AUT.i_P = c2d.X_PAO.i_P = 0.02154
    c2d.X_H.i_K = c2d.X_AUT.i_K = c2d.X_PAO.i_K = 0.0
    c2d.X_H.i_Mg = c2d.X_AUT.i_Mg = c2d.X_PAO.i_Mg = 0.0
    c2d.X_H.i_Ca = c2d.X_AUT.i_Ca = c2d.X_PAO.i_Ca = 0.0
    c2d.X_H.i_mass = c2d.X_AUT.i_mass = c2d.X_PAO.i_mass = 0.90
    c2d.X_H.f_Vmass_Totmass = c2d.X_AUT.f_Vmass_Totmass = c2d.X_PAO.f_Vmass_Totmass = 0.85
    
    c2d.X_PHA.i_C = 0.3
    c2d.X_PHA.i_mass = 0.55
    c2d.X_PHA.f_Vmass_Totmass = 0.92727
    c2d.X_PP.i_charge = 0
    
    S_Ca = Component.from_chemical('S_Ca', chemical='Ca+2', measured_as='Ca',
                                   description='Calcium', **ion_kwargs)
    
    X_CaCO3 = Component.from_chemical('X_CaCO3', chemical='CaCO3', 
                                   description='Calcite', **mineral_kwargs)
    X_struv = Component.from_chemical('X_struv', chemical='MgNH4PO4(H2O)6', 
                                   description='Struvite', **mineral_kwargs)
    X_newb = Component.from_chemical('X_newb', chemical='MgHPO4(H2O)3', 
                                   description='Newberyite', **mineral_kwargs)
    X_ACP = Component.from_chemical('X_ACP', chemical='Ca3P2O8', 
                                   description='Amorphous calcium phosphate', 
                                   **mineral_kwargs)
    X_MgCO3 = Component.from_chemical('X_MgCO3', chemical='MgCO3', 
                                   description='Magnesite', **mineral_kwargs)
    X_AlOH = Component.from_chemical('X_AlOH', chemical='Al(OH)3', 
                                   description='Aluminum hydroxide', **mineral_kwargs)
    X_AlPO4 = Component.from_chemical('X_AlPO4', chemical='AlPO4', 
                                   description='Aluminum phosphate', **mineral_kwargs)    
    X_FeOH = c2d.X_MeOH.copy('X_FeOH')
    X_FePO4 = c2d.X_MeP.copy('X_FePO4')
    
    S_Na = Component.from_chemical('S_Na', chemical='Na+', description='Sodium', **ion_kwargs)
    S_Cl = Component.from_chemical('S_Cl', chemical='Cl-', description='Chloride', **ion_kwargs)
    
    H2O = c2d.H2O
    
    for cmp in (c2d.S_F, c2d.X_S, c2d.S_I, c2d.X_I, c2d.X_H, c2d.X_AUT, c2d.X_PHA):
        cmp.i_NOD = None
    c2d.refresh_constants()
    c2d = [*c2d]
    solubles = c2d[:8]  # replace S_ALK with S_IC
    particulates = c2d[9:-3] # exclude X_MeOH, X_MeP, H2O
    
    cmps = Components([*solubles, S_IC, S_K, S_Mg, *particulates, 
                       S_Ca, X_CaCO3, X_struv, X_newb, X_ACP, X_MgCO3, 
                       X_AlOH, X_AlPO4, X_FeOH, X_FePO4, S_Na, S_Cl, H2O])
    cmps.default_compile()
    if set_thermo: settings.set_thermo(cmps)

    return cmps


#%%
_rhos = np.zeros(21)
def rhos_asm2d(state_arr, params):
    S_O2, S_N2, S_NH4, S_NO3, S_PO4, S_F, S_A, S_I, S_ALK, \
        X_I, X_S, X_H, X_PAO, X_PP, X_PHA, X_AUT, X_MeOH, X_MeP = state_arr[:18]
      
    _rhos[:19] = 0.
    if X_H > 0:
        K_h = params['K_h']
        K_O2 = params['K_O2']
        K_X = params['K_X']
        K_NO3 = params['K_NO3']
        eta_NO3 = params['eta_NO3']
        eta_fe = params['eta_fe']
        _rhos[:3] = K_h*(X_S/X_H)/(K_X+X_S/X_H)*X_H
        _rhos[0] *= S_O2/(K_O2+S_O2)
        _rhos[1] *= eta_NO3*K_O2/(K_O2+S_O2)*S_NO3/(K_NO3+S_NO3)
        _rhos[2] *= eta_fe*K_O2/(K_O2+S_O2)*K_NO3/(K_NO3+S_NO3)
        
        mu_H = params['mu_H']
        K_O2_H = params['K_O2_H']
        K_F = params['K_F']
        K_A_H = params['K_A_H']
        K_NH4_H = params['K_NH4_H']
        K_P_H = params['K_P_H']
        K_ALK_H = params['K_ALK_H']
        K_NO3_H = params['K_NO3_H']
        eta_NO3_H = params['eta_NO3_H']                
        _rhos[3:7] = mu_H*S_NH4/(K_NH4_H+S_NH4)*S_PO4/(K_P_H+S_PO4)*S_ALK/(K_ALK_H+S_ALK)*X_H
        _rhos[[3,5]] *= S_F/(K_F+S_F)*S_F/(S_F+S_A)
        _rhos[[4,6]] *= S_A/(K_A_H+S_A)*S_A/(S_F+S_A)
        _rhos[3:5] *= S_O2/(K_O2_H+S_O2)
        _rhos[5:7] *= eta_NO3_H*K_O2_H/(K_O2_H+S_O2)*S_NO3/(K_NO3_H+S_NO3)
    
        q_fe = params['q_fe']
        K_fe = params['K_fe']
        b_H = params['b_H']
        _rhos[7] = q_fe*K_O2_H/(K_O2_H+S_O2)*K_NO3_H/(K_NO3_H+S_NO3)*S_F/(K_fe+S_F)*S_ALK/(K_ALK_H+S_ALK)*X_H
        _rhos[8] = b_H*X_H
    
    K_ALK_PAO = params['K_ALK_PAO']
    if X_PAO > 0:
        q_PHA = params['q_PHA']
        K_A_PAO = params['K_A_PAO']
        K_PP = params['K_PP']
        _rhos[9] = q_PHA*S_A/(K_A_PAO+S_A)*S_ALK/(K_ALK_PAO+S_ALK)*(X_PP/X_PAO)/(K_PP+X_PP/X_PAO)*X_PAO
    
        q_PP = params['q_PP']
        K_O2_PAO = params['K_O2_PAO']
        K_PS = params['K_PS']
        K_PHA = params['K_PHA']
        K_MAX = params['K_MAX']
        K_IPP = params['K_IPP']
        eta_NO3_PAO = params['eta_NO3_PAO']
        K_NO3_PAO = params['K_NO3_PAO']
        mu_PAO = params['mu_PAO']
        K_NH4_PAO = params['K_NH4_PAO']
        K_P_PAO = params['K_P_PAO']
        _rhos[10:12] = q_PP*S_PO4/(K_PS+S_PO4)*S_ALK/(K_ALK_PAO+S_ALK)*(X_PHA/X_PAO)/(K_PHA+X_PHA/X_PAO)*(K_MAX-X_PP/X_PAO)/(K_IPP+K_MAX-X_PP/X_PAO)*X_PAO
        _rhos[12:14] = mu_PAO*S_NH4/(K_NH4_PAO+S_NH4)*S_PO4/(K_P_PAO+S_PO4)*S_ALK/(K_ALK_PAO+S_ALK)*(X_PHA/X_PAO)/(K_PHA+X_PHA/X_PAO)*X_PAO
        _rhos[[10,12]] *= S_O2/(K_O2_PAO+S_O2)
        _rhos[[11,13]] *= eta_NO3_PAO*K_O2_PAO/(K_O2_PAO+S_O2)*S_NO3/(K_NO3_PAO+S_NO3)

        b_PAO = params['b_PAO']
        _rhos[14] = b_PAO*X_PAO
    
    b_PP = params['b_PP']
    b_PHA = params['b_PHA']
    _rhos[15] = b_PP*X_PP
    _rhos[16] = b_PHA*X_PHA
    _rhos[14:17] *= S_ALK/(K_ALK_PAO+S_ALK)

    if X_AUT > 0:
        mu_AUT = params['mu_AUT']
        K_O2_AUT = params['K_O2_AUT']
        K_NH4_AUT = params['K_NH4_AUT']
        K_P_AUT = params['K_P_AUT']
        K_ALK_AUT = params['K_ALK_AUT']
        b_AUT = params['b_AUT']
        _rhos[17] = mu_AUT*S_O2/(K_O2_AUT+S_O2)*S_NH4/(K_NH4_AUT+S_NH4)*S_PO4/(K_P_AUT+S_PO4)*S_ALK/(K_ALK_AUT+S_ALK)*X_AUT
        _rhos[18] = b_AUT*X_AUT
    
    k_PRE = params['k_PRE']
    k_RED = params['k_RED']
    K_ALK_PRE = params['K_ALK_PRE']
    _rhos[19] = k_PRE*S_PO4*X_MeOH
    _rhos[20] = k_RED*X_MeP*S_ALK/(K_ALK_PRE+S_ALK)
    
    return _rhos

@chemicals_user
class ASM2d(CompiledProcesses):
    '''
    Activated Sludge Model No. 2d in original notation.

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

    Examples
    --------
    >>> from qsdsan import processes as pc
    >>> cmps = pc.create_asm2d_cmps()
    >>> asm2d = pc.ASM2d()
    >>> asm2d.show()
    ASM2d([aero_hydrolysis, anox_hydrolysis, anae_hydrolysis, hetero_growth_S_F, hetero_growth_S_A, denitri_S_F, denitri_S_A, ferment, hetero_lysis, PAO_storage_PHA, aero_storage_PP, anox_storage_PP, PAO_aero_growth_PHA, PAO_anox_growth, PAO_lysis, PP_lysis, PHA_lysis, auto_aero_growth, auto_lysis, precipitation, redissolution])

    References
    ----------
    [1] Henze, M.; Gujer, W.; Mino, T.; Loosdrecht, M. van. Activated Sludge
    Models: ASM1, ASM2, ASM2d and ASM3; IWA task group on mathematical modelling
    for design and operation of biological wastewater treatment, Ed.; IWA
    Publishing: London, 2000.
    
    [2] Rieger, L.; Gillot, S.; Langergraber, G.; Ohtsuki, T.; Shaw, A.; Takács,
    I.; Winkler, S. Guidelines for Using Activated Sludge Models; IWA Publishing:
    London, New York, 2012; Vol. 11.
    https://doi.org/10.2166/9781780401164.
    '''
    _params = ('f_SI', 'Y_H', 'f_XI_H', 'Y_PAO', 'Y_PO4', 'Y_PHA', 'f_XI_PAO', 'Y_A', 'f_XI_AUT',
               'K_h', 'eta_NO3', 'eta_fe', 'K_O2', 'K_NO3', 'K_X',
               'mu_H', 'q_fe', 'eta_NO3_H', 'b_H', 'K_O2_H', 'K_F', 'K_fe', 'K_A_H',
               'K_NO3_H', 'K_NH4_H', 'K_P_H', 'K_ALK_H',
               'q_PHA', 'q_PP', 'mu_PAO', 'eta_NO3_PAO', 'b_PAO', 'b_PP', 'b_PHA',
               'K_O2_PAO', 'K_NO3_PAO', 'K_A_PAO', 'K_NH4_PAO', 'K_PS','K_P_PAO',
               'K_ALK_PAO', 'K_PP', 'K_MAX', 'K_IPP', 'K_PHA',
               'mu_AUT', 'b_AUT', 'K_O2_AUT', 'K_NH4_AUT', 'K_ALK_AUT', 'K_P_AUT',
               'k_PRE', 'k_RED', 'K_ALK_PRE', 'COD_deN')


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

        self.compile(to_class=cls)
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
                            COD_deN=cmps.S_N2.i_COD-cmps.S_NO3.i_COD,
                            **kwargs)
        self.set_rate_function(rhos_asm2d)
        self.rate_function._params = self.parameters

        return self
    
#%%
_mpath = ospath.join(data_path, 'process_data/_masm2d.tsv')
_mmp = ospath.join(data_path, 'process_data/_mmp.tsv')

def acid_base_rxn(h_ion, ionic_states, Ka):
    K, Mg, Ca, Na, Cl, NOx, NH, IC, IP, Ac = ionic_states  # in M
    Kw, Knh, Kc1, Kc2, Kp1, Kp2, Kp3, Kac = Ka
    oh_ion = Kw/h_ion
    nh4 = NH * h_ion/(Knh + h_ion)
    ac = Ac * Kac/(Kac + h_ion)
    co2, hco3, co3 = ion_speciation(h_ion, Kc1, Kc2) * IC
    h3po4, h2po4, hpo4, po4 = ion_speciation(h_ion, Kp1, Kp2, Kp3) * IP
    return K + 2*Mg + 2*Ca + Na + h_ion + nh4 - Cl - NOx - oh_ion - ac - hco3 - 2*co3 - h2po4 - 2*hpo4 - 3*po4

def solve_pH(state_arr, Ka, unit_conversion):
    cmps_in_M = state_arr[:31] * unit_conversion *1e-3
    # S_K, S_Mg, S_Ca, S_Na, S_Cl, S_NO3, S_NH4, S_IC, S_PO4, S_A
    ions = cmps_in_M[[9, 10, 18, 28, 29, 3, 2, 8, 4, 6]]
    h = brenth(acid_base_rxn, 1e-14, 1.0,
               args=(ions, Ka),
               xtol=1e-12, maxiter=100)
    return h

# rhos = np.zeros(19+7+2) # 19 biological processes, 7 precipitation/dissociation, 2 gas stripping
rhos = np.zeros(19+7) # 19 biological processes, 7 precipitation/dissociation
def _rhos_masm2d(state_arr, params, acceptor_dependent_decay=True, h=None):
    if 'ks' not in params:
        k_h, mu_H, mu_PAO, mu_AUT, \
        q_fe, q_PHA, q_PP, \
        b_H, b_PAO, b_PP, b_PHA, b_AUT, \
        eta_NO3, eta_fe, eta_NO3_H, eta_NO3_PAO, \
        eta_NO3_Hl, eta_NO3_PAOl, eta_NO3_PPl, eta_NO3_PHAl, eta_NO3_AUTl, \
        K_O2, K_O2_H, K_O2_PAO, K_O2_AUT, \
        K_NO3, K_NO3_H, K_NO3_PAO, K_NO3_AUT, \
        K_X, K_F, K_fe, K_A_H, K_A_PAO, \
        K_NH4_H, K_NH4_PAO, K_NH4_AUT, \
        K_P_H, K_P_PAO, K_P_AUT, K_P_S, \
        K_PP, K_MAX, K_IPP, K_PHA, \
            = list(params.values())[:45]
        
        cmps = params['cmps']
        params['mass2mol'] = cmps.i_mass / cmps.chem_MW
        
        params['ks'] = ks = np.zeros(19)
        # rate constants
        ks[:3] = k_h
        ks[3:7] = mu_H
        ks[7:19] = (q_fe, b_H, q_PHA, q_PP, q_PP, mu_PAO, mu_PAO, b_PAO, b_PP, b_PHA, mu_AUT, b_AUT)
        # rate reduction factors
        ks[1] *= eta_NO3
        ks[2] *= eta_fe
        ks[5:7] *= eta_NO3_H
        ks[[11,13]] *= eta_NO3_PAO
        
        # half saturation / inhibition factors
        params['Ks_o2'] = np.array([K_O2, K_O2_H, K_O2_H, K_O2_PAO, K_O2_PAO, K_O2_AUT])
        params['Ks_no3'] = np.array([K_NO3, K_NO3_H, K_NO3_H, K_NO3_PAO, K_NO3_PAO, K_NO3_AUT])
        params['Ks_nh4'] = np.array([K_NH4_H, K_NH4_PAO, K_NH4_AUT])
        params['Ks_po4'] = np.array([K_P_H, K_P_PAO, K_P_AUT])
        params['eta_decay'] = np.array([eta_NO3_Hl, eta_NO3_PAOl, eta_NO3_PPl, eta_NO3_PHAl, eta_NO3_AUTl])
        
    ks = params['ks']
    Ks_o2 = params['Ks_o2']
    Ks_no3 = params['Ks_no3']
    Ks_nh4 = params['Ks_nh4']
    Ks_po4 = params['Ks_po4']
    eta_decay = params['eta_decay']
    
    Kx = params['K_X']
    Kf = params['K_F'] 
    Kfe = params['K_fe']
    Ka_H = params['K_A_H']
    Ka_PAO = params['K_A_PAO']
    Kp_stor = params['K_P_S']
    Kpp = params['K_PP']
    Kmax = params['K_MAX']
    Kipp = params['K_IPP']
    Kpha = params['K_PHA']
    
    S_O2, S_N2, S_NH4, S_NO3, S_PO4, S_F, S_A, S_I, S_IC, S_K, S_Mg, \
        X_I, X_S, X_H, X_PAO, X_PP, X_PHA, X_AUT, \
            S_Ca, X_CaCO3, X_struv, X_newb, X_ACP, X_MgCO3, \
                X_AlOH, X_AlPO4, X_FeOH, X_FePO4, S_Na, S_Cl \
                    = state_arr[:30]
    
    ############# biological processes ###############
    nutrients = Monod(S_NH4, Ks_nh4) * Monod(S_PO4, Ks_po4)

    rhos[:19] = ks
    rhos[:9] *= X_H 
    rhos[3:7] *= nutrients[0]
    rhos[9:15] *= X_PAO
    rhos[12:14] *= nutrients[1]
    rhos[15] *= X_PP
    rhos[16] *= X_PHA
    rhos[17:19] *= X_AUT
    rhos[17] *= nutrients[2]
    
    aero = Monod(S_O2, Ks_o2)
    anox = Monod(S_NO3, Ks_no3)
    rhos[[0,3,4,10,12,17]] *= aero                          # aerobic
    rhos[[1,5,6,11,13]] *= (1-aero[:5]) * anox[:5]          # anoxic
    rhos[[2,7]] *= (1-aero[:2]) * (1-anox[:2])              # anaerobic/fermentation
    
    if X_H > 0: rhos[:3] *= Monod(X_S/X_H, Kx)
    if S_F+S_A == 0: 
        rhos[3:7] = 0
    else:
        rhos[[3,5]] *= Monod(S_F, Kf) * S_F/(S_F+S_A)
        rhos[[4,6]] *= Monod(S_A, Ka_H) * S_A/(S_F+S_A)
    
    rhos[7] *= Monod(S_F, Kfe)
    if X_PAO > 0:
        pha = Monod(X_PHA/X_PAO, Kpha)
        rhos[9] *= Monod(S_A, Ka_PAO) * Monod(X_PP/X_PAO, Kpp)
        rhos[[10,11]] *= pha * Monod(Kmax-X_PP/X_PAO, Kipp) * Monod(S_PO4, Kp_stor)
        rhos[[12,13]] *= pha
    
    if acceptor_dependent_decay:
        rhos[8] *= (aero[1] + eta_decay[0]*(1-aero[1])*anox[1])
        rhos[14:17] *= (aero[3] +eta_decay[1:4]*(1-aero[3])*anox[3])
        rhos[18] *= (aero[5] + eta_decay[4]*(1-aero[5])*anox[5])
    
    ######### pH ############
    mass2mol = params['mass2mol']
    Ka = params['Ka']
    Kw, Knh, Kc1, Kc2, Kp1, Kp2, Kp3, Kac = Ka
    if h == None: h = solve_pH(state_arr, Ka, mass2mol)
    nh4 = state_arr[2] * h/(Knh + h)
    co2, hco3, co3 = state_arr[8] * ion_speciation(h, Kc1, Kc2)
    h3po4, h2po4, hpo4, po4 = state_arr[4] * ion_speciation(h, Kp1, Kp2, Kp3)
    
    ########## precipitation-dissolution #############
    k_mmp = params['k_mmp']
    Ksp = params['Ksp']
    # K_dis = params['K_dis']
    K_AlOH = params['K_AlOH']
    K_FeOH = params['K_FeOH']
    # f_dis = Monod(state_arr[19:24], K_dis[:5])
    # if X_CaCO3 > 0: rhos[19] = (S_Ca * co3 - Ksp[0]) * f_dis[0]
    # else: rhos[19] = S_Ca * co3
    # if X_struv > 0: rhos[20] = (S_Mg * nh4 * po4 - Ksp[1]) * f_dis[1]
    # else: rhos[20] = S_Mg * nh4 * po4
    # if X_newb > 0: rhos[21] = (S_Mg * hpo4 - Ksp[2]) * f_dis[2]
    # else: rhos[21] = S_Mg * hpo4
    # if X_ACP > 0: rhos[22] = (S_Ca**3 * po4**2 - Ksp[3]) * f_dis[3]
    # else: rhos[22] = S_Ca**3 * po4**2
    # if X_MgCO3 > 0: rhos[23] = (S_Mg * co3 - Ksp[4]) * f_dis[4]
    # else: rhos[23] = S_Mg * co3

    rhos[19:26] = 0.
    SI = (S_Ca * co3 / Ksp[0])**(1/2)
    if SI > 1: rhos[19] = X_CaCO3 * (SI-1)**2

    SI = (S_Mg * nh4 * po4 / Ksp[1])**(1/3)
    if SI > 1: rhos[20] = X_struv * (SI-1)**3

    SI = (S_Mg * hpo4 / Ksp[2])**(1/2)
    if SI > 1: rhos[21] =  X_newb * (SI-1)**2
    
    SI = (S_Ca**3 * po4**2 / Ksp[3])**(1/5)
    if SI > 1: rhos[22] = X_ACP * (SI-1)**2
    
    SI = (S_Mg * co3 / Ksp[4])**(1/2)
    if SI > 1: rhos[23] = X_MgCO3 * (SI-1)**2
    
    rhos[24] = X_AlOH * po4 * Monod(X_AlOH, K_AlOH)
    rhos[25] = X_FeOH * po4 * Monod(X_FeOH, K_FeOH)
    rhos[19:26] *= k_mmp
    
    return rhos

#%%
@chemicals_user
class mASM2d(CompiledProcesses):
    '''
    Modified ASM2d. Compatible with `ADM1p` for plant-wide simulations.
    Includes an algebraic pH solver and precipitation/dissolution of common minerals.

    Parameters
    ----------
    components : :class:`CompiledComponents`, optional
        Can be created with the `create_masm2d_cmps` function.
    path : str, optional
        File path for an alternative Petersen Matrix. The default is None.
    electron_acceptor_dependent_decay : bool, optional
        Whether biomass decay kinetics is dependent on concentrations of 
        electron acceptors. The default is True.
    pH_ctrl : float or None, optional
        Whether to fix pH at a specific value or solve for pH (`None`). The default is 7.0.
    k_h : float, optional
        Hydrolysis rate constant, in [d^(-1)]. The default is 3.0.
    eta_NO3_Hl : float, optional
        Anoxic reduction factor for endogenous respiration of heterotrophs, 
        unitless. The default is 0.5.
    eta_NO3_PAOl : float, optional
        Anoxic reduction factor for lysis of PAOs, unitless. The default is 0.33.
    eta_NO3_PPl : float, optional
        Anoxic reduction factor for lysis of PP, unitless. The default is 0.33.
    eta_NO3_PHAl : float, optional
        Anoxic reduction factor for lysis of PHA, unitless. The default is 0.33.
    eta_NO3_AUTl : float, optional
        Anoxic reduction factor for decay of autotrophs, unitless. The default is 0.33.
    K_NO3_AUT : float, optional
        Half saturation coefficient of NOx- for autotrophs [mg-N/L]. The default is 0.5.
    K_P_S : float, optional
        Half saturation coefficient of ortho-P for PP storage [mg-P/L]. The default is 0.2.
    k_mmp : iterable[float], optional
        Rate constants for multi-mineral precipitation/dissolution 
        [mg-precipitate/L/(unit of solubility product)/d]. Follows the exact order
        of `mASM2d._precipitates`. The default is (5.0, 300, 0.05, 150, 50, 1.0, 1.0).
    pKsp : iterable[float], optional
        Solubility of minerals, in order of `mASM2d._precipitates`. 
        The default is (6.45, 13.16, 5.8, 23, 7, 21, 26).
    K_dis : iterable[float], optional
        Saturation coefficient for the switching function of mineral dissolution
        [mg-precipitate/L]. The default is (1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0).
    K_AlOH : float, optional
        Half saturation coefficient of aluminum hydroxide for AlPO4 precipitation
        [mg-Al(OH)3/L]. The default is 0.001.
    K_FeOH : float, optional
        Half saturation coefficient of ferric hydroxide for FePO4 precipitation
        [mg-Fe(OH)3/L]. The default is 0.001.
    pKa : iterable[float], optional
        Equilibrium coefficient values of acid-base pairs, unitless, 
        following the order of `mASM2d._acid_base_pairs`. 
        The default is (14, 9.25, 6.37, 10.32, 2.12, 7.21, 12.32, 4.76).
    
    
    See Also
    --------
    :class:`qsdsan.processes.ASM2d`
    :class:`qsdsan.processes.ADM1p`

    Examples
    --------
    >>> import qsdsan.processes as pc
    >>> cmps = pc.create_masm2d_cmps()
    >>> asm = pc.mASM2d()
    >>> asm.show()
    mASM2d([aero_hydrolysis, anox_hydrolysis, anae_hydrolysis, hetero_growth_S_F, hetero_growth_S_A, denitri_S_F, denitri_S_A, ferment, hetero_lysis, storage_PHA, aero_storage_PP, anox_storage_PP, PAO_aero_growth_PHA, PAO_anox_growth, PAO_lysis, PP_lysis, PHA_lysis, auto_aero_growth, auto_lysis, CaCO3_precipitation_dissolution, struvite_precipitation_dissolution, newberyite_precipitation_dissolution, ACP_precipitation_dissolution, MgCO3_precipitation_dissolution, AlPO4_precipitation_dissolution, FePO4_precipitation_dissolution])

    >>> # Calculate process rate given state variable values and fixed pH.
    >>> import numpy as np
    >>> state_arr = np.ones(len(cmps))
    >>> rhos = asm.rate_function(state_arr)  # reaction rate for each process
    >>> for i,j in zip(asm.IDs, rhos): 
    ...     print(f'{i}{(40-len(i))*" "}{j:.3g}')
    aero_hydrolysis                         2.27
    anox_hydrolysis                         0.182
    anae_hydrolysis                         0.0606
    hetero_growth_S_F                       0.471
    hetero_growth_S_A                       0.471
    denitri_S_F                             0.0503
    denitri_S_A                             0.0503
    ferment                                 0.0333
    hetero_lysis                            0.356
    storage_PHA                             0.594
    aero_storage_PP                         1.06
    anox_storage_PP                         0.0851
    PAO_aero_growth_PHA                     0.778
    PAO_anox_growth                         0.0622
    PAO_lysis                               0.174
    PP_lysis                                0.174
    PHA_lysis                               0.174
    auto_aero_growth                        0.33
    auto_lysis                              0.111
    CaCO3_precipitation_dissolution         0
    struvite_precipitation_dissolution      0
    newberyite_precipitation_dissolution    0
    ACP_precipitation_dissolution           0
    MgCO3_precipitation_dissolution         0
    AlPO4_precipitation_dissolution         1.82e-11
    FePO4_precipitation_dissolution         1.82e-11

    >>> # Estimate pH given state variable values.
    >>> Ka = asm.rate_function.params['Ka']
    >>> unit_conversion = asm.rate_function.params['mass2mol']
    >>> h_ion = asm.solve_pH(state_arr, Ka, unit_conversion)
    >>> pH = -np.log10(h_ion)
    >>> print(f'{pH:.2f}')
    8.40

    References
    ----------
    [1] Henze, M., Gujer, W., Mino, T., & van Loosdrecht, M. (2000). 
    Activated Sludge Models: ASM1, ASM2, ASM2d and ASM3. In IWA task group 
    on mathematical modelling for design and operation of biological 
    wastewater treatment (Ed.), Scientific and Technical Report No. 9. 
    IWA Publishing.
    
    [2] Solon, K., Flores-Alsina, X., Kazadi Mbamba, C., Ikumi, D., Volcke, 
    E. I. P., Vaneeckhaute, C., Ekama, G., Vanrolleghem, P. A., Batstone, 
    D. J., Gernaey, K. V., & Jeppsson, U. (2017). Plant-wide modelling 
    of phosphorus transformations in wastewater treatment systems: 
    Impacts of control and operational strategies. Water Research, 113, 
    97–110. https://doi.org/10.1016/j.watres.2017.02.007
    '''
    _stoichio_params = ('f_SI', 'Y_H', 'Y_PAO', 'Y_PO4', 'Y_PHA', 'Y_A', 
                        'f_XI_H', 'f_XI_PAO', 'f_XI_AUT', 'COD_deN', 'K_XPP', 'Mg_XPP')
    _kinetic_params = ('k_h', 'mu_H', 'mu_PAO', 'mu_AUT', 
                       'q_fe', 'q_PHA', 'q_PP', 
                       'b_H', 'b_PAO', 'b_PP', 'b_PHA', 'b_AUT', 
                       'eta_NO3', 'eta_fe', 'eta_NO3_H', 'eta_NO3_PAO', 
                       'eta_NO3_Hl', 'eta_NO3_PAOl', 'eta_NO3_PPl', 'eta_NO3_PHAl', 'eta_NO3_AUTl',
                       'K_O2', 'K_O2_H', 'K_O2_PAO', 'K_O2_AUT', 
                       'K_NO3', 'K_NO3_H', 'K_NO3_PAO', 'K_NO3_AUT', 
                       'K_X', 'K_F', 'K_fe', 'K_A_H', 'K_A_PAO', 
                       'K_NH4_H', 'K_NH4_PAO', 'K_NH4_AUT', 
                       'K_P_H', 'K_P_PAO', 'K_P_AUT', 'K_P_S', 
                       'K_PP', 'K_MAX', 'K_IPP', 'K_PHA',
                       'k_mmp', 'Ksp', 'K_dis', 'K_AlOH', 'K_FeOH', 
                       # 'kLa_min', 'f_kLa', 'K_Henry', 
                       'Ka', 'cmps')
    
    _acid_base_pairs = (('H+', 'OH-'), ('NH4+', 'NH3'), 
                        ('CO2', 'HCO3-'), ('HCO3-', 'CO3-2'), 
                        ('H3PO4', 'H2PO4-'), ('H2PO4-', 'HPO4-2'), ('HPO4-2', 'PO4-3'),
                        ('HAc', 'Ac-'),)
    _precipitates = ('X_CaCO3', 'X_struv', 'X_newb', 'X_ACP', 'X_MgCO3', 'X_AlPO4', 'X_FePO4')
    
    gas_IDs = ['S_N2', 'S_IC']
    kLa_min = [3.0, 3.0]
    K_Henry = [6.5e-4, 3.5e-2]  # 20 degree C
    D_gas = [1.88e-9, 1.92e-9]  # diffusivity
    p_gas_atm = [0.78, 3.947e-4]# partial pressure in air
    
    def __new__(cls, components=None, path=None, 
                electron_acceptor_dependent_decay=True, pH_ctrl=7.0, 
                f_SI=0.0, Y_H=0.625, Y_PAO=0.625, Y_PO4=0.4, Y_PHA=0.2, Y_A=0.24, 
                f_XI_H=0.1, f_XI_PAO=0.1, f_XI_AUT=0.1,
                k_h=3.0, mu_H=6.0, mu_PAO=1.0, mu_AUT=1.0, 
                q_fe=3.0, q_PHA=3.0, q_PP=1.5, 
                b_H=0.4, b_PAO=0.2, b_PP=0.2, b_PHA=0.2, b_AUT=0.15, 
                eta_NO3=0.6, eta_fe=0.4, eta_NO3_H=0.8, eta_NO3_PAO=0.6, 
                eta_NO3_Hl=0.5, eta_NO3_PAOl=0.33, eta_NO3_PPl=0.33, eta_NO3_PHAl=0.33, eta_NO3_AUTl=0.33,
                K_O2=0.2, K_O2_H=0.2, K_O2_PAO=0.2, K_O2_AUT=0.5, 
                K_NO3=0.5, K_NO3_H=0.5, K_NO3_PAO=0.5, K_NO3_AUT=0.5, 
                K_X=0.1, K_F=4.0, K_fe=4.0, K_A_H=4.0, K_A_PAO=4.0, 
                K_NH4_H=0.05, K_NH4_PAO=0.05, K_NH4_AUT=1.0, 
                K_P_H=0.01, K_P_PAO=0.01, K_P_AUT=0.01, K_P_S=0.2, 
                K_PP=0.01, K_MAX=0.34, K_IPP=0.02, K_PHA=0.01,
                # k_mmp=(5.0, 300, 0.05, 150, 50, 1.0, 1.0),
                # pKsp=(6.45, 13.16, 5.8, 23, 7, 21, 26),
                # k_mmp=(0.024, 120, 0.024, 72, 0.024, 0.024, 0.024),  # Flores-Alsina 2016
                # pKsp=(8.3, 13.6, 18.175, 28.92, 7.46, 18.2, 37.76),  # Flores-Alsina 2016
                k_mmp=(8.4, 240, 1.0, 72, 1.0, 1.0e-5, 1.0e-5),              # MATLAB
                pKsp=(8.45, 13.5, 5.7, 29.1, 7.4, 18.2, 26.4),               # MINTEQ (except newberyite), 20 C    
                K_dis=(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
                K_AlOH=0.001, K_FeOH=0.001, 
                pKa=(14, 9.25, 6.37, 10.32, 2.12, 7.21, 12.32, 4.76),
                **kwargs):       
        

        if not path: path = _mpath

        cmps = _load_components(components)

        self = Processes.load_from_file(path,
                                        components=cmps,
                                        conserved_for=('COD', 'C', 'N', 'P'),
                                        parameters=cls._stoichio_params,
                                        compile=False)

        mmp = Processes.load_from_file(_mmp, components=cmps, 
                                       conserved_for=(), compile=False)
        mmp_stoichio = {}
        df = load_data(_mmp)
        for i, j in df.iterrows():
            j.dropna(inplace=True)
            key = j.index[j == 1][0]
            j = j.to_dict()
            j.pop(key)
            mmp_stoichio[key] = j
        mol_to_mass = cmps.chem_MW / cmps.i_mass
        Ksp_mass = np.array([10**(-p) for p in pKsp])     # mass in mg/L or g/m3
        i = 0
        for pd, xid in zip(mmp, cls._precipitates):
            for k,v in mmp_stoichio[xid].items():
                m2m = mol_to_mass[cmps.index(k)] * 1e3
                Ksp_mass[i] *= m2m**abs(v)
            i += 1
            pd._stoichiometry *= mol_to_mass
            pd.ref_component = xid
            
        self.extend(mmp)
        
        # for gas in cls._gas_stripping:
        #     new_p = Process('%s_stripping' % gas.lstrip('S_'),
        #                     reaction={gas:-1},
        #                     ref_component=gas,
        #                     conserved_for=(),)
        #     self.append(new_p)
        self.compile(to_class=cls)
        
        dct = self.__dict__
        dct.update(kwargs)
        dct['K_Henry'] = [K*mol_to_mass[cmps.index(i)]*1000 for K, i in zip(cls.K_Henry, cls.gas_IDs)]
        dct['mmp_stoichio'] = mmp_stoichio
        stoichio_vals = (f_SI, Y_H, Y_PAO, Y_PO4, Y_PHA, Y_A, 
                         f_XI_H, f_XI_PAO, f_XI_AUT, cmps.S_N2.i_COD-cmps.S_NO3.i_COD,
                         cmps.X_PP.i_K, cmps.X_PP.i_Mg)
        dct['_parameters'] = dict(zip(cls._stoichio_params, stoichio_vals))
        dct['_edecay'] = bool(electron_acceptor_dependent_decay)
        dct['pH_ctrl'] = pH_ctrl
        if pH_ctrl: h = 10**(-pH_ctrl)
        else: h = None
        rhos_masm2d = lambda state_arr, params: _rhos_masm2d(state_arr, params, electron_acceptor_dependent_decay, h)
        self.set_rate_function(rhos_masm2d)
        Ka = np.array([10**(-p) for p in pKa])
        # f_kLa = np.array(cls.D_gas)/cls.D_O2
        kinetic_vals = (k_h, mu_H, mu_PAO, mu_AUT, 
                        q_fe, q_PHA, q_PP, 
                        b_H, b_PAO, b_PP, b_PHA, b_AUT, 
                        eta_NO3, eta_fe, eta_NO3_H, eta_NO3_PAO, 
                        eta_NO3_Hl, eta_NO3_PAOl, eta_NO3_PPl, eta_NO3_PHAl, eta_NO3_AUTl,
                        K_O2, K_O2_H, K_O2_PAO, K_O2_AUT, 
                        K_NO3, K_NO3_H, K_NO3_PAO, K_NO3_AUT, 
                        K_X, K_F, K_fe, K_A_H, K_A_PAO, 
                        K_NH4_H, K_NH4_PAO, K_NH4_AUT, 
                        K_P_H, K_P_PAO, K_P_AUT, K_P_S, 
                        K_PP, K_MAX, K_IPP, K_PHA,
                        np.array(k_mmp), Ksp_mass, 
                        np.array(K_dis), K_AlOH, K_FeOH, 
                        # kLa_min, f_kLa, K_Henry, 
                        Ka, cmps,
                        )
        self.rate_function._params = dict(zip(cls._kinetic_params, kinetic_vals))
        dct['solve_pH'] = solve_pH
        return self
    
    @property
    def electron_acceptor_dependent_decay(self):
        '''[bool] Whether the decay rate is dependent on electron acceptor (O2, NO3-) concentrations'''
        return self._edecay
    @electron_acceptor_dependent_decay.setter
    def electron_acceptor_dependent_decay(self, dependent):
        self._edecay = edecay = bool(dependent)
        rhos_masm2d = lambda state_arr, params: _rhos_masm2d(state_arr, params, edecay)
        self.set_rate_function(rhos_masm2d)
    
    def set_parameters(self, **parameters):
        stoichio = self._parameters
        kinetic = self.rate_function._params
        if set(parameters.keys()).intersection(set(kinetic.keys())):
            for key in ('ks', 'Ks_o2', 'Ks_no3', 'Ks_nh4', 'Ks_po4', 'eta_decay'):
                kinetic.pop(key, None)
        for k,v in parameters.items():
            if k in self._kinetic_params: kinetic[k] = v
            else: stoichio[k] = v
        if self._stoichio_lambdified is not None:
            self.__dict__['_stoichio_lambdified'] = None
    
    def set_pKsps(self, ps):
        cmps = self.components
        mol_to_mass = cmps.chem_MW / cmps.i_mass
        idxer = cmps.index
        stoichio = self.mmp_stoichio
        Ksp_mass = []    # mass in mg/L or g/m3
        for xid, p in zip(self._precipitates, ps):
            K = 10**(-p)
            for cmp, v in stoichio[xid]:
                m2m = mol_to_mass[idxer(cmp)] * 1e3
                K *= m2m**abs(v)
            Ksp_mass.append(K)
        self.rate_function._params['Ksp'] = np.array(Ksp_mass)
