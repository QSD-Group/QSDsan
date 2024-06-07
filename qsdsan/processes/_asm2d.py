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
from qsdsan import Component, Components, Process, Processes, CompiledProcesses
from ..utils import ospath, data_path

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

# create_asm2d_cmps()

def create_masm2d_cmps(set_thermo=True):
    c2d = create_asm2d_cmps(False)
    S_K = Component.from_chemical('S_K', chemical='K',
                                  measured_as='K',
                                  description='Potassium',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
    
    S_Mg = Component.from_chemical('S_Mg', chemical='Mg',
                                   measured_as='Mg',
                                   description='Magnesium',
                                   particle_size='Soluble',
                                   degradability='Undegradable',
                                   organic=False)

    S_IC = c2d.S_ALK.copy('S_IC')
    
    c2d.S_F.i_C = c2d.X_S.i_C = 0.31843
    c2d.S_F.i_N = c2d.X_S.i_N = 0.03352
    c2d.S_F.i_P = c2d.X_S.i_P = 5.59e-3
    
    c2d.S_I.i_C = c2d.X_I.i_C = 0.36178
    c2d.S_I.i_N = c2d.X_I.i_N = 0.06003
    c2d.S_I.i_P = c2d.X_I.i_P = 6.49e-3
    c2d.S_F.i_mass = c2d.X_S.i_mass = c2d.S_I.i_mass = c2d.X_I.i_mass = 0.75
    c2d.S_F.f_Vmass_Totmass = c2d.X_S.f_Vmass_Totmass = c2d.S_I.f_Vmass_Totmass = c2d.X_I.f_Vmass_Totmass = 0.85
    
    c2d.X_H.i_C = c2d.X_AUT.i_C = c2d.X_PAO.i_C = 0.36612
    c2d.X_H.i_N = c2d.X_AUT.i_N = c2d.X_PAO.i_N = 0.08615
    c2d.X_H.i_P = c2d.X_AUT.i_P = c2d.X_PAO.i_P = 0.02154
    c2d.X_H.i_mass = c2d.X_AUT.i_mass = c2d.X_PAO.i_mass = 0.90
    c2d.X_H.f_Vmass_Totmass = c2d.X_AUT.f_Vmass_Totmass = c2d.X_PAO.f_Vmass_Totmass = 0.85
    
    c2d.X_PHA.i_C = 0.3
    c2d.X_PHA.i_mass = 0.55
    c2d.X_PHA.f_Vmass_Totmass = 0.92727
    c2d.X_PP.i_charge = 0
    
    for cmp in (c2d.S_F, c2d.X_S, c2d.S_I, c2d.X_I, c2d.X_H, c2d.X_AUT, c2d.X_PHA):
        cmp.i_NOD = None
    c2d.refresh_constants()
    c2d = [*c2d]
    solubles = c2d[:8]  # replace S_ALK with S_IC
    others = c2d[9:]
    
    cmps = Components([*solubles, S_IC, S_K, S_Mg, *others])
    cmps.default_compile()
    if set_thermo: settings.set_thermo(cmps)

    return cmps


#%%
@chemicals_user
class ASM2d(CompiledProcesses):
    '''
    Activated Sludge Model No. 2d in original notation. [1]_, [2]_

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
    >>> from qsdsan import processes as pc, set_thermo
    >>> cmps = pc.create_asm2d_cmps()
    >>> asm2d = pc.ASM2d()
    >>> asm2d.show()
    ASM2d([aero_hydrolysis, anox_hydrolysis, anae_hydrolysis, hetero_growth_S_F, hetero_growth_S_A, denitri_S_F, denitri_S_A, ferment, hetero_lysis, PAO_storage_PHA, aero_storage_PP, PAO_aero_growth_PHA, PAO_lysis, PP_lysis, PHA_lysis, auto_aero_growth, auto_lysis, precipitation, redissolution, anox_storage_PP, PAO_anox_growth])

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
               'K_O2_PAO', 'K_NO3_PAO', 'K_A_PAO', 'K_NH4_PAO', 'K_PS','K_P_PAO',
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

        if path == _path:
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
                            **kwargs)
        return self
    
#%%
_mpath = ospath.join(data_path, 'process_data/_masm2d.tsv')

Monod = lambda S, K: S/(S+K)

rhos = np.zeros(19) # 19 biological processes, no precipitation/dissociation or gas stripping yet
def _rhos_masm2d(state_arr, params, acceptor_dependent_decay=True):
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
            = params.values()
        
        params['ks'] = ks = rhos * 0
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
    
    S_O2, S_NH4, S_NO3, S_PO4, S_F, S_A, \
        X_S, X_H, X_PAO, X_PP, X_PHA, X_AUT \
            = state_arr[[0,2,3,4,5,6,12,13,14,15,16,17]]
    
    nutrients = Monod(S_NH4, Ks_nh4) * Monod(S_PO4, Ks_po4)

    rhos[:] = ks
    rhos[:9] *= X_H * nutrients[0]
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
    if S_F > 0: rhos[[3,5]] *= Monod(S_F, Kf) * S_F/(S_F+S_A)
    else: rhos[[3,5]] = 0.        
    if S_A > 0: rhos[[4,6]] *= Monod(S_A, Ka_H) * S_A/(S_F+S_A)
    else: rhos[[4,6]] = 0.
    
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
        
    return rhos

@chemicals_user
class mASM2d(CompiledProcesses):
    '''
    Modified ASM2d. [1]_, [2]_ Compatible with `ADM1p` for plant-wide simulations.

    Parameters
    ----------

    Examples
    --------

    References
    ----------
    .. [1] Henze, M., Gujer, W., Mino, T., & van Loosdrecht, M. (2000). 
        Activated Sludge Models: ASM1, ASM2, ASM2d and ASM3. In IWA task group 
        on mathematical modelling for design and operation of biological 
        wastewater treatment (Ed.), Scientific and Technical Report No. 9. 
        IWA Publishing.
    .. [2] Solon, K., Flores-Alsina, X., Kazadi Mbamba, C., Ikumi, D., Volcke, 
        E. I. P., Vaneeckhaute, C., Ekama, G., Vanrolleghem, P. A., Batstone, 
        D. J., Gernaey, K. V., & Jeppsson, U. (2017). Plant-wide modelling 
        of phosphorus transformations in wastewater treatment systems: 
        Impacts of control and operational strategies. Water Research, 113, 
        97–110. https://doi.org/10.1016/j.watres.2017.02.007
    
    '''
    _stoichio_params = ('f_SI', 'Y_H', 'Y_PAO', 'Y_PO4', 'Y_PHA', 'Y_A', 
                        'f_XI_H', 'f_XI_PAO', 'f_XI_AUT','K_XPP', 'Mg_XPP')
    _kinetic_params = ('k_h', 'mu_H', 'mu_PAO', 'mu_AUT', 
                       'q_fe', 'q_PHA', 'q_PP', 
                       'b_H', 'b_PAO', 'b_PP', 'b_PHA', 'b_AUT', 
                       # k_PRE, k_RED, 
                       'eta_NO3', 'eta_fe', 'eta_NO3_H', 'eta_NO3_PAO', 
                       'eta_NO3_Hl', 'eta_NO3_PAOl', 'eta_NO3_PPl', 'eta_NO3_PHAl', 'eta_NO3_AUTl',
                       'K_O2', 'K_O2_H', 'K_O2_PAO', 'K_O2_AUT', 
                       'K_NO3', 'K_NO3_H', 'K_NO3_PAO', 'K_NO3_AUT', 
                       'K_X', 'K_F', 'K_fe', 'K_A_H', 'K_A_PAO', 
                       'K_NH4_H', 'K_NH4_PAO', 'K_NH4_AUT', 
                       'K_P_H', 'K_P_PAO', 'K_P_AUT', 'K_P_S', 
                       'K_PP', 'K_MAX', 'K_IPP', 'K_PHA',
                       # 'K_ALK_PRE'
                       )
    
    decay_dependon_electron_acceptor = True

    def __new__(cls, components=None, path=None, 
                f_SI=0.0, Y_H=0.625, Y_PAO=0.625, Y_PO4=0.4, Y_PHA=0.2, Y_A=0.24, 
                f_XI_H=0.1, f_XI_PAO=0.1, f_XI_AUT=0.1,
                k_h=3.0, mu_H=6.0, mu_PAO=1.0, mu_AUT=1.0, 
                q_fe=3.0, q_PHA=3.0, q_PP=1.5, 
                b_H=0.4, b_PAO=0.2, b_PP=0.2, b_PHA=0.2, b_AUT=0.15, 
                # k_PRE=1.0, k_RED=0.6, 
                eta_NO3=0.6, eta_fe=0.4, eta_NO3_H=0.8, eta_NO3_PAO=0.6, 
                eta_NO3_Hl=0.5, eta_NO3_PAOl=0.33, eta_NO3_PPl=0.33, eta_NO3_PHAl=0.33, eta_NO3_AUTl=0.33,
                K_O2=0.2, K_O2_H=0.2, K_O2_PAO=0.2, K_O2_AUT=0.5, 
                K_NO3=0.5, K_NO3_H=0.5, K_NO3_PAO=0.5, K_NO3_AUT=0.5, 
                K_X=0.1, K_F=4.0, K_fe=4.0, K_A_H=4.0, K_A_PAO=4.0, 
                K_NH4_H=0.05, K_NH4_PAO=0.05, K_NH4_AUT=1.0, 
                K_P_H=0.01, K_P_PAO=0.01, K_P_AUT=0.01, K_P_S=0.2, 
                K_PP=0.01, K_MAX=0.34, K_IPP=0.02, K_PHA=0.01,
                # K_ALK_PRE=0.5,
                #!!! kLa and/or solubility values for gas stripping
                #!!! precipitation kinetics
                **kwargs):

        if not path: path = _mpath

        cmps = _load_components(components)

        self = Processes.load_from_file(path,
                                        components=cmps,
                                        conserved_for=('COD', 'C', 'N', 'P',),
                                        parameters=cls._stoichio_params,
                                        compile=False)

        if path == _mpath:
            _p12 = Process('anox_storage_PP',
                           'S_PO4 + [K_XPP]S_K + [Mg_XPP]S_Mg +[Y_PHA]X_PHA + [?]S_NO3 -> X_PP + [?]S_N2 + [?]S_NH4 + [?]S_IC',
                           components=cmps,
                           ref_component='X_PP',
                           conserved_for=('C', 'N', 'NOD', 'COD'))

            _p14 = Process('PAO_anox_growth',
                           '[1/Y_PAO]X_PHA + [?]S_NO3 + [?]S_PO4 -> X_PAO + [?]S_N2 + [?]S_NH4  + [?]S_IC',
                           components=cmps,
                           ref_component='X_PAO',
                           conserved_for=('C', 'N', 'P', 'NOD', 'COD'))
            
            self.insert(11, _p12)
            self.insert(13, _p14)

        #!!! add gas stripping

        self.compile(to_class=cls)
        
        dct = self.__dict__
        dct.update(kwargs)
        stoichio_vals = (f_SI, Y_H, Y_PAO, Y_PO4, Y_PHA, Y_A, 
                         f_XI_H, f_XI_PAO, f_XI_AUT,
                         cmps.X_PP.i_K, cmps.X_PP.i_Mg)
        dct['_parameters'] = dict(zip(cls._stoichio_params, stoichio_vals))
        rhos_masm2d = lambda state_arr, params: _rhos_masm2d(state_arr, params, cls.decay_dependon_electron_acceptor)
        self.set_rate_function(rhos_masm2d)
        kinetic_vals = (k_h, mu_H, mu_PAO, mu_AUT, 
                        q_fe, q_PHA, q_PP, 
                        b_H, b_PAO, b_PP, b_PHA, b_AUT, 
                        # k_PRE, k_RED, 
                        eta_NO3, eta_fe, eta_NO3_H, eta_NO3_PAO, 
                        eta_NO3_Hl, eta_NO3_PAOl, eta_NO3_PPl, eta_NO3_PHAl, eta_NO3_AUTl,
                        K_O2, K_O2_H, K_O2_PAO, K_O2_AUT, 
                        K_NO3, K_NO3_H, K_NO3_PAO, K_NO3_AUT, 
                        K_X, K_F, K_fe, K_A_H, K_A_PAO, 
                        K_NH4_H, K_NH4_PAO, K_NH4_AUT, 
                        K_P_H, K_P_PAO, K_P_AUT, K_P_S, 
                        K_PP, K_MAX, K_IPP, K_PHA,
                        # K_ALK_PRE
                        )
        self.rate_function._params = dict(zip(cls._kinetic_params, kinetic_vals))
        
        return self
    
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