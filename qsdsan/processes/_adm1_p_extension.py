# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>
    Saumitra Rai <raisaumitra9@gmail.com>
    

Part of this module is based on the Thermosteam package:
https://github.com/BioSTEAMDevelopmentGroup/thermosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from thermosteam.utils import chemicals_user
from thermosteam import settings
from chemicals.elements import molecular_weight as get_mw
from qsdsan import Component, Components, Process, Processes, CompiledProcesses
from qsdsan.processes import (
    create_adm1_cmps, 
    create_asm2d_cmps,
    T_correction_factor,
    non_compet_inhibit,
    grad_non_compet_inhibit,
    substr_inhibit,
    grad_substr_inhibit,
    mass2mol_conversion,
    Hill_inhibit,
    ADM1,
    TempState
    )
from qsdsan.utils import ospath, data_path
from scipy.optimize import brenth
from warnings import warn
import numpy as np

__all__ = ('create_adm1_p_extension_cmps', 
           'ADM1_p_extension',
           'rhos_adm1_p_extension')

_path = ospath.join(data_path, 'process_data/_adm1_p_extension.tsv')
_load_components = settings.get_default_chemicals

#%%
# =============================================================================
# ADM1 (with P extension) -specific components
# =============================================================================

C_mw = get_mw({'C':1})
N_mw = get_mw({'N':1})
P_mw = get_mw({'P':1})

def create_adm1_p_extension_cmps(set_thermo=True):
    c1 = create_adm1_cmps(False)
    c2d = create_asm2d_cmps(False)
    
    S_IP = c2d.S_PO4.copy('S_IP')
    
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
    
    c = [*c1]
    Ss = c[:11]
    Xs = c[13:-3]   # X_c is excluded
    others = c[-3:]

    cmps_adm1_p_extension = Components([*Ss, S_IP, c1.S_I,
                                        *Xs, c2d.X_PHA, c2d.X_PP, c2d.X_PAO, 
                                        S_K, S_Mg, c2d.X_MeOH, c2d.X_MeP,
                                        *others])
    cmps_adm1_p_extension.default_compile()
    if set_thermo: settings.set_thermo(cmps_adm1_p_extension)
    return cmps_adm1_p_extension

# create_adm1_p_extension_cmps()


#%%
# =============================================================================
# kinetic rate functions
# =============================================================================

R = 8.3145e-2 # Universal gas constant, [bar/M/K]


def acid_base_rxn(h_ion, weak_acids_tot, Kas):
    # S_cat, S_K, S_Mg, S_an, S_IN, S_IP, S_IC, S_ac, S_pro, S_bu, S_va = weak_acids_tot  # in M
    S_cat, S_K, S_Mg, S_an, S_IN, S_IP = weak_acids_tot[:6]
    # Kw, Ka_nh, Ka_h2po4, Ka_co2, Ka_ac, Ka_pr, Ka_bu, Ka_va = Kas
    Kw = Kas[0]
    oh_ion = Kw/h_ion
    nh3, hpo4, hco3, ac, pro, bu, va = Kas[1:] * weak_acids_tot[4:] / (Kas[1:] + h_ion)
    return S_cat + S_K + 2*S_Mg + h_ion + (S_IN - nh3) - S_an - oh_ion - hco3 - ac - pro - bu - va - 2*hpo4 - (S_IP - hpo4)

# The function 'fprime_abr' is not used in the code
def fprime_abr(h_ion, weak_acids_tot, Kas):
    # S_cat, S_K, S_Mg, S_an, S_IN, S_IP = weak_acids_tot[:6]
    Kw = Kas[0]
    doh_ion = - Kw / h_ion ** 2
    dnh3, dhpo4, dhco3, dac, dpro, dbu, dva = - Kas[1:] * weak_acids_tot[4:] / (Kas[1:] + h_ion)**2
    return 1 + (-dnh3) - doh_ion - dhco3 - dac - dpro - dbu - dva - dhpo4


rhos = np.zeros(28) # 28 kinetic processes (25 as defined in modified ADM1 + 3 for gases)
Cs = np.empty(25) # 25 processes as defined in modified ADM1

def solve_pH(state_arr, Ka, unit_conversion):
    cmps_in_M = state_arr[:34] * unit_conversion
    # S_cat, S_K, S_Mg, S_an, S_IN, S_IP, S_IC, S_ac, S_pro, S_bu, S_va
    # Kw, Ka_nh, Ka_h2po4, Ka_co2, Ka_ac, Ka_pr, Ka_bu, Ka_va = Kas
    weak_acids = cmps_in_M[[31, 27, 28, 32, 10, 11, 9, 6, 5, 4, 3]]
    h = brenth(acid_base_rxn, 1e-14, 1.0,
            args=(weak_acids, Ka),
            xtol=1e-12, maxiter=100)
    nh3 = Ka[1] * weak_acids[4] / (Ka[1] + h)
    co2 = weak_acids[6] - Ka[3] * weak_acids[6] / (Ka[3] + h)
    return h, nh3, co2

rhos_adm1_p_extension = lambda state_arr, params: _rhos_adm1p(state_arr, params, h=None)

def _rhos_adm1p(state_arr, params, h=None):
    ks = params['rate_constants']
    Ks = params['half_sat_coeffs']
    
    cmps = params['components']
    pH_ULs = params['pH_ULs']
    pH_LLs = params['pH_LLs']
    KS_IN = params['KS_IN']
    KS_IP = params['KS_IP']
    KI_nh3 = params['KI_nh3']
    KIs_h2 = params['KIs_h2']
    KHb = params['K_H_base']
    Kab = params['Ka_base']
    KH_dH = params['K_H_dH']
    Ka_dH = params['Ka_dH']
    kLa = params['kLa']
    T_base = params['T_base']
    root = params['root']
    
    if 'unit_conv' not in params: params['unit_conv'] = mass2mol_conversion(cmps)
    unit_conversion = params['unit_conv']
    
    # state_arr_cmps stated just for readability of code 
    # {0: 'S_su',
    #  1: 'S_aa',
    #  2: 'S_fa',
    #  3: 'S_va',
    #  4: 'S_bu',
    #  5: 'S_pro',
    #  6: 'S_ac',
    #  7: 'S_h2',
    #  8: 'S_ch4',
    #  9: 'S_IC',
    #  10: 'S_IN',
    #  11: 'S_IP',
    #  12: 'S_I',
    #  13: 'X_ch',
    #  14: 'X_pr',
    #  15: 'X_li',
    #  16: 'X_su',
    #  17: 'X_aa',
    #  18: 'X_fa',
    #  19: 'X_c4',
    #  20: 'X_pro',
    #  21: 'X_ac',
    #  22: 'X_h2',
    #  23: 'X_I',
    #  24: 'X_PHA',
    #  25: 'X_PP',
    #  26: 'X_PAO',
    #  27: 'S_K',
    #  28: 'S_Mg',
    #  29: 'X_MeOH',
    #  30: 'X_MeP',
    #  31: 'S_cat',
    #  32: 'S_an',
    #  33: 'H2O'}

    Cs[:7] = state_arr[13:20]
    Cs[7:11] = state_arr[19:23]
    Cs[11:18] = state_arr[16:23]
    Cs[18:23] = X_PAO = state_arr[26]
    Cs[23] = X_PP = state_arr[25]
    Cs[24] = state_arr[24]
    
    substrates = state_arr[:8]
    
    S_va, S_bu, S_h2, S_IN, S_IP = state_arr[[3,4,7,10,11]]
    
    T_op = state_arr[-1]
    if T_op == T_base:
        Ka = Kab
        KH = KHb / unit_conversion[7:10]
    else:
        T_temp = params.pop('T_op', None)
        if T_op != T_temp:
            params['Ka'] = Kab * T_correction_factor(T_base, T_op, Ka_dH)
            params['KH'] = KHb * T_correction_factor(T_base, T_op, KH_dH) / unit_conversion[7:10]        
        params['T_op'] = T_op
        Ka = params['Ka']
        KH = params['KH']

    rhos[:-3] = ks * Cs
    rhos[3:11] *= substr_inhibit(substrates, Ks[:8])
    if S_va > 0: rhos[6] *= 1/(1+S_bu/S_va)
    if S_bu > 0: rhos[7] *= 1/(1+S_va/S_bu)
    
    vfas = state_arr[3:7]
    
    if X_PAO > 0:
        K_A, K_PP = Ks[-2:]
        rhos[18:22] *= substr_inhibit(vfas, K_A) * substr_inhibit(X_PP/X_PAO, K_PP)
    
    if sum(vfas) > 0:
        rhos[18:22] *= vfas/sum(vfas)
    
    biogas_S = state_arr[7:10].copy()
    biogas_p = R * T_op * state_arr[34:37]

    if h is None:
        h, nh3, co2 = solve_pH(state_arr, Ka, unit_conversion)
    else:
        nh3 = Ka[1] * S_IN * unit_conversion[10] / (Ka[1] + h)
        S_IC = state_arr[9] * unit_conversion[9]
        co2 = S_IC - Ka[3] * S_IC / (Ka[3] + h)
    biogas_S[-1] = co2 / unit_conversion[9]
    
    Iph = Hill_inhibit(h, pH_ULs, pH_LLs)
    Iin = substr_inhibit(S_IN, KS_IN)
    Iip = substr_inhibit(S_IP, KS_IP)
    Ih2 = non_compet_inhibit(S_h2, KIs_h2)
    Inh3 = non_compet_inhibit(nh3, KI_nh3)
    root.data = {
        'pH':-np.log10(h), 
        'Iph':Iph, 
        'Ih2':Ih2, 
        'Iin':Iin, 
        'Inh3':Inh3,
        }
    rhos[3:11] *= Iph * Iin * Iip
    rhos[5:9] *= Ih2
    rhos[9] *= Inh3
    rhos[-3:] = kLa * (biogas_S - KH * biogas_p)
    
    # print(rhos)
    return rhos

def dydt_Sh2_AD(S_h2, state_arr, h, params, f_stoichio, V_liq, S_h2_in):
    state_arr[7] = S_h2
    Q = state_arr[37]
    rxn = _rhos_adm1p(state_arr, params, h=h)
    stoichio = f_stoichio(state_arr)  # should return the stoichiometric coefficients of S_h2 for all processes
    return Q/V_liq*(S_h2_in - S_h2) + np.dot(rxn, stoichio)

grad_rhos = np.zeros(5)
X_bio = np.zeros(5)
def grad_dydt_Sh2_AD(S_h2, state_arr, h, params, f_stoichio, V_liq, S_h2_in):
    state_arr[7] = S_h2
    ks = params['rate_constants'][[5,6,7,8,10]]
    Ks = params['half_sat_coeffs'][2:6]
    K_h2 = params['half_sat_coeffs'][7]
    pH_ULs = params['pH_ULs']
    pH_LLs = params['pH_LLs']
    KS_IN = params['KS_IN']
    KIs_h2 = params['KIs_h2']
    kLa = params['kLa']
    
    X_bio[:] = state_arr[[18,19,19,20,22]]
    substrates = state_arr[2:6]
    S_va, S_bu, S_IN = state_arr[[3,4,10]]
    Iph = Hill_inhibit(h, pH_ULs, pH_LLs)[[2,3,4,5,7]]
    Iin = substr_inhibit(S_IN, KS_IN)
    grad_Ih2 = grad_non_compet_inhibit(S_h2, KIs_h2)

    grad_rhos[:] = ks * X_bio * Iph * Iin
    grad_rhos[:-1] *= substr_inhibit(substrates, Ks) * grad_Ih2
    if S_va > 0: grad_rhos[1] *= 1/(1+S_bu/S_va)
    if S_bu > 0: grad_rhos[2] *= 1/(1+S_va/S_bu)
    
    grad_rhos[-1] *= grad_substr_inhibit(S_h2, K_h2)
    stoichio = f_stoichio(state_arr)

    Q = state_arr[37]
    return -Q/V_liq + np.dot(grad_rhos, stoichio[[5,6,7,8,10]]) + kLa*stoichio[-3]

#%%
# =============================================================================
# ADM1_p_extension class
# =============================================================================

@chemicals_user
class ADM1_p_extension(ADM1):
    """
    Anaerobic Digestion Model No.1 with P extension. [1]_, [2]_, [3]_

    Parameters
    ----------
    components : class:`CompiledComponents`, optional
        Components corresponding to each entry in the stoichiometry array,
        defaults to thermosteam.settings.chemicals.
    path : str, optional
        Alternative file path for the Petersen matrix. The default is None.
    f_sI_xb : float, optional
        fraction of soluble inerts from biomass. The default is 0.
    f_ch_xb : float, optional
        fraction of carbohydrates from biomass. The default is 0.275. 
    f_pr_xb : float, optional
        fraction of proteins from biomass. The default is 0.275. 
    f_li_xb : float, optional
        fraction of lipids from biomass. The default is 0.35. 
    f_xI_xb : float, optional
        fraction of particulate inerts from biomass. The default is 0.1. 
    f_ac_PHA : float, optional
        Yield of acetate on PHA [kg COD/kg COD]. The default is 0.4.
    f_bu_PHA : float, optional
        Yield of butyrate on PHA [kg COD/kg COD]. The default is 0.1.
    f_pro_PHA : float, optional
        Yield of propionate on PHA [kg COD/kg COD]. The default is 0.4.
    f_va_PHA : float, optional
        Yield of valerate on PHA [kg COD/kg COD]. The default is 0.1.
    Y_PO4 : float, optional
        Yield of biomass on phosphate [kmol P/kg COD]. The default is 0.013. 
    K_A : float, optional
        VFAs half saturation coefficient for PHA storage [kg COD/m3]. The default is 0.004.
    K_PP : float, optional
        Half saturation coefficient for polyphosphate [kmol PP/kg PAO COD].
        The default is 0.00032.
    q_PHA : float, optional
        Rate constant for storage of PHA [d^(-1)]. The default is 3.
    b_PAO : float, optional
        Lysis rate of PAOs [d^(-1)]. The default is 0.2.
    b_PP : float, optional
        Lysis rate of polyphosphates [d^(-1)]. The default is 0.2.
    b_PHA : float, optional
        Lysis rate of PHAs [d^(-1)]. The default is 0.2.
    KS_IP : float, optional
        P limitation for inorganic phosphorous [kmol P/m3]. 
        The default is 2e-5. 
    pKa_base : iterable[float], optional
        pKa (equilibrium coefficient) values of acid-base pairs at the base 
        temperature, unitless, following the order of `ADM1_p_extension._acid_base_pairs`.
        The default is [14, 9.25, 7.20, 6.35, 4.76, 4.88, 4.82, 4.86].
    Ka_dH : iterable[float], optional
        Heat of reaction of each acid-base pair at base temperature [J/mol], 
        following the order of `ADM1_p_extension._acid_base_pairs`. The default is 
        [55900, 51965, 3600, 7646, 0, 0, 0, 0].

    See Also
    --------
    :class:`qsdsan.processes.ADM1`
    
    Examples
    --------
    >>> from qsdsan import processes as pc
    >>> cmps = pc.create_adm1_p_extension_cmps()
    >>> adm1_p = pc.ADM1_p_extension()
    >>> adm1_p.show()
    ADM1_p_extension([hydrolysis_carbs, hydrolysis_proteins, hydrolysis_lipids, uptake_sugars, uptake_amino_acids, uptake_LCFA, uptake_valerate, uptake_butyrate, uptake_propionate, uptake_acetate, uptake_h2, decay_Xsu, decay_Xaa, decay_Xfa, decay_Xc4, decay_Xpro, decay_Xac, decay_Xh2, storage_Sva_in_XPHA, storage_Sbu_in_XPHA, storage_Spro_in_XPHA, storage_Sac_in_XPHA, lysis_XPAO, lysis_XPP, lysis_XPHA, h2_transfer, ch4_transfer, IC_transfer])
    
    References
    ----------
    .. [1] Batstone, D. J.; Keller, J.; Angelidaki, I.; Kalyuzhnyi, S. V; 
        Pavlostathis, S. G.; Rozzi, A.; Sanders, W. T. M.; Siegrist, H.; 
        Vavilin, V. A. The IWA Anaerobic Digestion Model No 1 (ADM1). 
        Water Sci. Technol. 2002, 45 (10), 65–73.
    .. [2] Rosen, C.; Jeppsson, U. Aspects on ADM1 Implementation within 
        the BSM2 Framework; Lund, 2006.
    .. [3] Flores-Alsina, X.; Solon, K.; Kazadi Mbamba, C.; Tait, S.; 
        Gernaey, K. V.; Jeppsson, U.; Batstone, D. J. 
        Modelling phosphorus (P), sulfur (S) and iron (FE) interactions for 
        dynamic simulations of anaerobic digestion processes. Water Research. 2016,
        95, 370–382.
    """

    _stoichio_params = (*ADM1._stoichio_params[5:],
                        'f_sI_xb', 'f_ch_xb', 'f_pr_xb', 'f_li_xb', 'f_xI_xb', 
                        'f_ac_PHA', 'f_bu_PHA', 'f_pro_PHA', 'f_va_PHA', 
                        'Y_PO4', 'K_XPP', 'Mg_XPP')
    
    _kinetic_params = (*ADM1._kinetic_params, 'KS_IP', )
    
    _acid_base_pairs = (('H+', 'OH-'), ('NH4+', 'NH3'), ('H2PO4-', 'HPO4-2'), 
                        ('CO2', 'HCO3-'), ('HAc', 'Ac-'), ('HPr', 'Pr-'),
                        ('HBu', 'Bu-'), ('HVa', 'Va-'))
    
    _biogas_IDs = ('S_h2', 'S_ch4', 'S_IC')

    def __new__(cls, components=None, path=None, N_xc=2.686e-3, N_I=4.286e-3, N_aa=7e-3,
                f_sI_xb=0, f_ch_xb=0.275, f_pr_xb=0.275, f_li_xb=0.350,
                f_fa_li=0.95, f_bu_su=0.13, f_pro_su=0.27, f_ac_su=0.41,
                f_va_aa=0.23, f_bu_aa=0.26, f_pro_aa=0.05, f_ac_aa=0.4,
                f_ac_fa=0.7, f_pro_va=0.54, f_ac_va=0.31, f_ac_bu=0.8, f_ac_pro=0.57,
                f_ac_PHA=0.4, f_bu_PHA=0.1, f_pro_PHA=0.4,
                Y_su=0.1, Y_aa=0.08, Y_fa=0.06, Y_c4=0.06, Y_pro=0.04, Y_ac=0.05, Y_h2=0.06, Y_PO4=0.013,
                q_dis=0.5, q_ch_hyd=10, q_pr_hyd=10, q_li_hyd=10,
                k_su=30, k_aa=50, k_fa=6, k_c4=20, k_pro=13, k_ac=8, k_h2=35,
                K_su=0.5, K_aa=0.3, K_fa=0.4, K_c4=0.2, K_pro=0.1, K_ac=0.15, K_h2=7e-6, 
                K_A=4e-3, K_PP=32e-5,
                b_su=0.02, b_aa=0.02, b_fa=0.02, b_c4=0.02, b_pro=0.02, b_ac=0.02, b_h2=0.02,
                q_PHA=3, b_PAO=0.2, b_PP=0.2, b_PHA=0.2, 
                KI_h2_fa=5e-6, KI_h2_c4=1e-5, KI_h2_pro=3.5e-6, KI_nh3=1.8e-3, KS_IN=1e-4, KS_IP=2e-5, 
                pH_limits_aa=(4,5.5), pH_limits_ac=(6,7), pH_limits_h2=(5,6),
                T_base=298.15, pKa_base=[14, 9.25, 7.20, 6.35, 4.76, 4.88, 4.82, 4.86],
                Ka_dH=[55900, 51965, 3600, 7646, 0, 0, 0, 0],
                kLa=200, K_H_base=[7.8e-4, 1.4e-3, 3.5e-2],
                K_H_dH=[-4180, -14240, -19410],
                **kwargs):
        
        cmps = _load_components(components)
        # Sure that some things are missing here! (Saumitra)
        # cmps.X_c.i_N = N_xc * N_mw
        cmps.X_I.i_N = cmps.S_I.i_N = N_I * N_mw
        cmps.S_aa.i_N = cmps.X_pr.i_N = N_aa * N_mw

        if not path: path = _path
        self = Processes.load_from_file(path,
                                        components=cmps,
                                        conserved_for=('C', 'N', 'P'),
                                        parameters=cls._stoichio_params,
                                        compile=False)
            
        gas_transfer = []
        for i in cls._biogas_IDs:
            new_p = Process('%s_transfer' % i.lstrip('S_'),
                            reaction={i:-1},
                            ref_component=i,
                            conserved_for=(),
                            parameters=())
            gas_transfer.append(new_p)
        self.extend(gas_transfer)
        self.compile(to_class=cls)

        stoichio_vals = (f_fa_li, f_bu_su, f_pro_su, f_ac_su, 1-f_bu_su-f_pro_su-f_ac_su,
                         f_va_aa, f_bu_aa, f_pro_aa, f_ac_aa, 1-f_va_aa-f_bu_aa-f_pro_aa-f_ac_aa,
                         f_ac_fa, 1-f_ac_fa, f_pro_va, f_ac_va, 1-f_pro_va-f_ac_va,
                         f_ac_bu, 1-f_ac_bu, f_ac_pro, 1-f_ac_pro,
                         Y_su, Y_aa, Y_fa, Y_c4, Y_pro, Y_ac, Y_h2, 
                         # new parameters
                         f_sI_xb, f_ch_xb, f_pr_xb, f_li_xb, 1-f_sI_xb-f_ch_xb-f_pr_xb-f_li_xb,
                         f_ac_PHA, f_bu_PHA, f_pro_PHA, 1-f_ac_PHA-f_bu_PHA-f_pro_PHA,
                         Y_PO4*P_mw, cmps.X_PP.i_K, cmps.X_PP.i_Mg)
        pH_LLs = np.array([pH_limits_aa[0]]*6 + [pH_limits_ac[0], pH_limits_h2[0]])
        pH_ULs = np.array([pH_limits_aa[1]]*6 + [pH_limits_ac[1], pH_limits_h2[1]])
        
        ks = np.array((q_ch_hyd, q_pr_hyd, q_li_hyd,
                       k_su, k_aa, k_fa, k_c4, k_c4, k_pro, k_ac, k_h2,
                       b_su, b_aa, b_fa, b_c4, b_pro, b_ac, b_h2,
                       q_PHA, q_PHA, q_PHA, q_PHA, b_PAO, b_PP, b_PHA))
        
        Ks = np.array((K_su, K_aa, K_fa, K_c4, K_c4, K_pro, K_ac, K_h2, 
                       #!!! new
                       K_A, K_PP))
        
        KIs_h2 = np.array((KI_h2_fa, KI_h2_c4, KI_h2_c4, KI_h2_pro))
        K_H_base = np.array(K_H_base)
        K_H_dH = np.array(K_H_dH)
        Ka_base = np.array([10**(-pKa) for pKa in pKa_base])
        Ka_dH = np.array(Ka_dH)
        root = TempState()
        dct = self.__dict__
        dct.update(kwargs)

        self.set_rate_function(rhos_adm1_p_extension)
        dct['_parameters'] = dict(zip(cls._stoichio_params, stoichio_vals))
        self.rate_function._params = dict(zip(cls._kinetic_params,
                                              [ks, Ks, pH_ULs, pH_LLs, KS_IN*N_mw, 
                                               KI_nh3, KIs_h2, Ka_base, Ka_dH,
                                               K_H_base, K_H_dH, kLa,
                                               T_base, self._components, root, 
                                               #!!! new parameter
                                               KS_IP*P_mw]))
        dct['solve_pH'] = solve_pH
        dct['dydt_Sh2_AD'] = dydt_Sh2_AD
        dct['grad_dydt_Sh2_AD'] = grad_dydt_Sh2_AD
        return self

    def set_half_sat_K(self, K, process):
        '''Set the substrate half saturation coefficient [kg/m3] for a process given its ID.'''
        i = self._find_index(process)
        if i < 11:
            self.rate_function._params['half_sat_coeffs'][i-3] = K
        else:
            ValueError('To set "K_A", specify process = -2; to set "K_PP", specify process = -1,'
                       f'not {process}')

    def set_pH_inhibit_bounds(self, process, lower=None, upper=None):
        '''Set the upper and/or lower limit(s) of pH inhibition [unitless] for a process given its ID.'''
        i = self._find_index(process) - 3
        dct = self.rate_function._params
        if lower is None: lower = dct['pH_LLs'][i]
        else: dct['pH_LLs'][i] = lower
        if upper is None: upper = dct['pH_ULs'][i]
        else: dct['pH_ULs'][i] = upper
        if lower >= upper:
            raise ValueError(f'lower limit for pH inhibition of {process} must '
                             f'be lower than the upper limit, not {[lower, upper]}')

    def set_h2_inhibit_K(self, KI, process):
        '''Set the H2 inhibition coefficient [kg/m3] for a process given its ID.'''
        i = self._find_index(process)
        self.rate_function._params['KIs_h2'][i-5] = KI

        
    def set_KS_IP(self, K):
        '''Set inhibition coefficient for inorganic phosphorous as a secondary
        substrate [M phosphorous].'''
        self.rate_function._params['KS_IP'] = K * P_mw


    def check_stoichiometric_parameters(self):
        '''Check whether product COD fractions sum up to 1 for each process.'''
        stoichio = self.parameters
        subst = ('xb', 'su', 'aa', 'fa', 'va', 'bu', 'pro', 'PHA')
        for s in subst:
            f_tot = sum([stoichio[k] for k in self._stoichio_params[:-7] \
                         if k.endswith(s)])
            if f_tot != 1:
                raise ValueError(f"the sum of 'f_()_{s}' values must equal 1")