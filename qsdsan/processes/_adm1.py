# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

Part of this module is based on the Thermosteam package:
https://github.com/BioSTEAMDevelopmentGroup/thermosteam

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from thermosteam.utils import chemicals_user
from thermosteam import settings
# from chemicals.elements import molecular_weight as get_mw
from qsdsan import Component, Components, Process, Processes, CompiledProcesses
import numpy as np
from qsdsan.utils import ospath, data_path
from scipy.optimize import brenth
from warnings import warn
from . import (
    non_compet_inhibit, grad_non_compet_inhibit, 
    substr_inhibit, grad_substr_inhibit,
    mass2mol_conversion, 
    T_correction_factor, R,
    TempState
    ) 

__all__ = ('create_adm1_cmps', 'ADM1',
           'pH_inhibit', 'Hill_inhibit', 
           'rhos_adm1', )

_path = ospath.join(data_path, 'process_data/_adm1.tsv')
_load_components = settings.get_default_chemicals

#%%
# =============================================================================
# ADM1-specific components
# =============================================================================

# C_mw = get_mw({'C':1})
# N_mw = get_mw({'N':1})
C_mw = 12
N_mw = 14

def create_adm1_cmps(set_thermo=True, adjust_MW_to_measured_as=True):
    cmps_all = Components.load_default()

    # varies
    X_c = cmps_all.X_OHO.copy('X_c')
    X_c.description = 'Composite'
    X_c.i_C = 0.02786 * C_mw
    X_c.i_N = 0.0376
    X_ch = Component.from_chemical('X_ch', chemical='glycogen', Tc=1011.4, # glucan
                                    description='Carbohydrates',
                                    measured_as='COD',
                                    particle_size='Particulate',
                                    degradability='Slowly',
                                    organic=True)
    X_ch.copy_models_from(X_c, names=('mu',))
    # X_ch = cmps_all.X_B_Subst.copy('X_ch')
    # X_ch.i_N = 0
    # X_ch.i_C = 0.0313 * C_mw

    # varies
    X_pr = cmps_all.X_B_Subst.copy('X_pr')
    X_pr.i_N = 0.007 * N_mw
    X_pr.i_C = 0.03 * C_mw

    X_li = Component.from_chemical('X_li', chemical='tripalmitin',
                                    description='Lipids',
                                    measured_as='COD',
                                    particle_size='Particulate',
                                    degradability='Slowly',
                                    organic=True)

    # both varies
    X_I = cmps_all.X_U_Inf.copy('X_I')
    S_I = cmps_all.S_U_Inf.copy('S_I')
    X_I.i_C = S_I.i_C = 0.03 * C_mw
    X_I.i_N = S_I.i_N = 0.06

    S_su = Component.from_chemical('S_su', chemical='glucose',
                                    description='Monosaccharides',
                                    measured_as='COD',
                                    particle_size='Soluble',
                                    degradability='Readily',
                                    organic=True)
    # S_su = cmps_all.S_F.copy('S_su')
    # S_su.i_N = 0
    # S_su.i_C = 0.0313 * 12

    # varies
    S_aa = cmps_all.S_F.copy('S_aa')
    S_aa.i_N = 0.007 * N_mw
    S_aa.i_P = 0
    S_aa.i_C = 0.03 * C_mw

    S_fa = Component.from_chemical('S_fa', chemical='palmitate',
                                    description='Total long-chain fatty acids',
                                    measured_as='COD',
                                    particle_size='Soluble',
                                    degradability='Readily',
                                    organic=True)

    S_va = Component.from_chemical('S_va', chemical='valerate',
                                    description='Total valerate',
                                    measured_as='COD',
                                    particle_size='Soluble',
                                    degradability='Readily',
                                    organic=True)

    S_bu = Component.from_chemical('S_bu', chemical='butyrate',
                                    description='Total butyrate',
                                    measured_as='COD',
                                    particle_size='Soluble',
                                    degradability='Readily',
                                    organic=True)

    S_pro = cmps_all.S_Prop.copy('S_pro')
    S_ac = cmps_all.S_Ac.copy('S_ac')
    S_h2 = cmps_all.S_H2.copy('S_h2')
    S_ch4 = cmps_all.S_CH4.copy('S_ch4')

    S_IC = Component.from_chemical('S_IC', chemical='CO2',
                                    measured_as='C',
                                    description='Inorganic carbon',
                                    particle_size='Dissolved gas',
                                    degradability='Undegradable',
                                    organic=False)
    S_IC.copy_models_from(S_ch4, ('Cn',))

    S_IN = Component.from_chemical('S_IN', chemical='NH3',
                                    measured_as='N',
                                    description='Inorganic nitrogen',
                                    particle_size='Soluble',
                                    degradability='Undegradable',
                                    organic=False)

    X_su = cmps_all.X_FO.copy('X_su')
    X_su.description = 'Biomass uptaking sugars'

    X_aa = cmps_all.X_FO.copy('X_aa')
    X_aa.description = 'Biomass uptaking amino acids'

    X_fa = cmps_all.X_FO.copy('X_fa')
    X_fa.description = 'Biomass uptaking long chain fatty acids'

    X_c4 = cmps_all.X_FO.copy('X_c4')
    X_c4.description = 'Biomass uptaking c4 fatty acids'

    X_pro = cmps_all.X_PRO.copy('X_pro')
    X_ac = cmps_all.X_ACO.copy('X_ac')
    X_h2 = cmps_all.X_HMO.copy('X_h2')

    for bio in (X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2):
        # bio.formula = 'C5H7O2N'
        bio.i_C = 0.0313 * C_mw
        bio.i_N = 0.08

    S_cat = cmps_all.S_CAT.copy('S_cat')
    S_an = cmps_all.S_AN.copy('S_an')
    S_cat.i_mass = S_an.i_mass = 1

    cmps_adm1 = Components([S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2,
                            S_ch4, S_IC, S_IN, S_I, X_c, X_ch, X_pr, X_li,
                            X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I,
                            S_cat, S_an, cmps_all.H2O])
    cmps_adm1.default_compile(ignore_inaccurate_molar_weight=True,
                              adjust_MW_to_measured_as=adjust_MW_to_measured_as)
    if set_thermo: settings.set_thermo(cmps_adm1)
    return cmps_adm1

# create_adm1_cmps()


#%%
# =============================================================================
# kinetic rate functions
# =============================================================================

def acid_base_rxn(h_ion, weak_acids_tot, Kas):
    # h, nh4, hco3, ac, pr, bu, va = mols
    # S_cat, S_an, S_IN, S_IC, S_ac, S_pro, S_bu, S_va = weak_acids_tot  # in M
    S_cat, S_an, S_IN = weak_acids_tot[:3]
    # Kw, Ka_nh, Ka_co2, Ka_ac, Ka_pr, Ka_bu, Ka_va = Kas
    Kw = Kas[0]
    oh_ion = Kw/h_ion
    nh3, hco3, ac, pro, bu, va = Kas[1:] * weak_acids_tot[2:] / (Kas[1:] + h_ion)
    return S_cat + h_ion + (S_IN - nh3) - S_an - oh_ion - hco3 - ac - pro - bu - va

def fprime_abr(h_ion, weak_acids_tot, Kas):
    S_cat, S_an, S_IN = weak_acids_tot[:3]
    Kw = Kas[0]
    doh_ion = - Kw / h_ion ** 2
    dnh3, dhco3, dac, dpro, dbu, dva = - Kas[1:] * weak_acids_tot[2:] / (Kas[1:] + h_ion)**2
    return 1 + (-dnh3) - doh_ion - dhco3 - dhco3 - dac - dpro - dbu - dva

def pH_inhibit(pH, ul, ll, lower_only=True):
    if lower_only:
        low_by = np.minimum(pH-ul, 0)
        return np.exp(-3 * (low_by/(ul-ll))**2)
    else:
        return (1+2*10**(0.5*(ll-ul)))/(1+10**(pH-ul)+10**(ll-pH))

def Hill_inhibit(H_ion, ul, ll):
    n = 3/(ul-ll)
    K = 10**(-(ul+ll)/2)
    return 1/(1+(H_ion/K) ** n)

rhos = np.zeros(22) # 22 kinetic processes
Cs = np.empty(19)

def solve_pH(state_arr, Ka, unit_conversion):
    cmps_in_M = state_arr[:27] * unit_conversion
    weak_acids = cmps_in_M[[24, 25, 10, 9, 6, 5, 4, 3]]
    h = brenth(acid_base_rxn, 1e-14, 1.0,
               args=(weak_acids, Ka),
               xtol=1e-12, maxiter=100)
    return h
rhos_adm1 = lambda state_arr, params: _rhos_adm1(state_arr, params, h=None)

def _rhos_adm1(state_arr, params, h=None):
    ks = params['rate_constants']
    Ks = params['half_sat_coeffs']
    cmps = params['components']
    pH_ULs = params['pH_ULs']
    pH_LLs = params['pH_LLs']
    KS_IN = params['KS_IN']
    KI_nh3 = params['KI_nh3']
    KIs_h2 = params['KIs_h2']
    KHb = params['K_H_base']
    Kab = params['Ka_base']
    KH_dH = params['K_H_dH']
    Ka_dH = params['Ka_dH']
    kLa = params['kLa']
    T_base = params['T_base']
    root = params['root']
    if 'unit_conv' in params:
        unit_conversion = params['unit_conv']
    else:
        unit_conversion = params['unit_conv'] = mass2mol_conversion(cmps)

    Cs[:8] = state_arr[12:20]
    Cs[8:12] = state_arr[19:23]
    Cs[12:] = state_arr[16:23]

    substrates = state_arr[:8]
    S_va, S_bu, S_h2, S_IN = state_arr[[3,4,7,10]]

    T_op = state_arr[-1]
    # Ka, KH = T_corrected_params(T_op, params)
    if T_op == T_base:
        Ka = Kab
        KH = KHb / unit_conversion[7:10]
    else:
        T_temp = params.pop('T_op', None)
        if T_op == T_temp:
            params['T_op'] = T_op
            Ka = params['Ka']
            KH = params['KH']
        else:
            params['T_op'] = T_op
            Ka = params['Ka'] = Kab * T_correction_factor(T_base, T_op, Ka_dH)
            KH = params['KH'] = KHb * T_correction_factor(T_base, T_op, KH_dH) / unit_conversion[7:10]

    biogas_S = state_arr[7:10].copy()
    biogas_p = R * T_op * state_arr[27:30]

    rhos[:-3] = ks * Cs
    Monod = substr_inhibit(substrates, Ks)
    rhos[4:12] *= Monod
    if S_va > 0: rhos[7] *= 1/(1+S_bu/S_va)
    if S_bu > 0: rhos[8] *= 1/(1+S_va/S_bu)

    if h is None: h = solve_pH(state_arr, Ka, unit_conversion)
    nh3 = S_IN * unit_conversion[10] * Ka[1] / (Ka[1] + h)
    co2 = state_arr[9] * h / (Ka[2] + h)
    biogas_S[-1] = co2
    
    Iph = Hill_inhibit(h, pH_ULs, pH_LLs)
    Iin = substr_inhibit(S_IN, KS_IN)
    Ih2 = non_compet_inhibit(S_h2, KIs_h2)
    Inh3 = non_compet_inhibit(nh3, KI_nh3)
    rhos[4:12] *= Iph * Iin
    rhos[6:10] *= Ih2
    rhos[10] *= Inh3
    root.data = {
        'pH':-np.log10(h), 
        'Iph':Iph, 
        'Ih2':Ih2, 
        'Iin':Iin, 
        'Inh3':Inh3,
        'Monod':Monod,
        'rhos':rhos[4:12].copy()
        }
    rhos[-3:] = kLa * (biogas_S - KH * biogas_p)
    return rhos

def dydt_Sh2_AD(S_h2, state_arr, h, params, f_stoichio, V_liq, S_h2_in):
    state_arr[7] = S_h2
    Q = state_arr[30]
    rxn = _rhos_adm1(state_arr, params, h=h)
    stoichio = f_stoichio(state_arr)  # should return the stoichiometric coefficients of S_h2 for all processes
    return Q/V_liq*(S_h2_in - S_h2) + np.dot(rxn, stoichio)

grad_rhos = np.zeros(5)
X_bio = np.zeros(5)
def grad_dydt_Sh2_AD(S_h2, state_arr, h, params, f_stoichio, V_liq, S_h2_in):
    state_arr[7] = S_h2
    ks = params['rate_constants'][[6,7,8,9,11]]
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

    Q = state_arr[30]
    return -Q/V_liq + np.dot(grad_rhos, stoichio[[6,7,8,9,11]]) + kLa*stoichio[-3]
    

#%%
# =============================================================================
# ADM1 class
# =============================================================================

@chemicals_user
class ADM1(CompiledProcesses):
    """
    Anaerobic Digestion Model No.1. [1]_, [2]_

    Parameters
    ----------
    components : class:`CompiledComponents`, optional
        Components corresponding to each entry in the stoichiometry array,
        defaults to thermosteam.settings.chemicals.
    path : str, optional
        Alternative file path for the Petersen matrix. The default is None.
    N_xc : float, optional
        Nitrogen content of composite materials [kmol N/kg COD]. The default is 2.686e-3.
    N_I : float, optional
        Nitrogen content of inert organics [kmol N/kg COD]. The default is 4.286e-3.
    N_aa : float, optional
        Nitrogen content of amino acids [kmol N/kg COD]. The default is 7e-3.
    f_ch_xc : float, optional
        Fraction of carbohydrates from composite disintegration [kg COD/kg COD]. 
        The default is 0.2.
    f_pr_xc : float, optional
        Fraction of proteins from composite disintegration [kg COD/kg COD].
        The default is 0.2.
    f_li_xc : float, optional
        Fraction of lipids from composite disintegration [kg COD/kg COD]. 
        The default is 0.3.
    f_xI_xc : float, optional
        Fraction of inert particulates from composite disintegration 
        [kg COD/kg COD]. The default is 0.2.
    f_fa_li : float, optional
        Fraction of long chain fatty acids (LCFAs) from hydrolysis of lipids
        [kg COD/kg COD]. The default is 0.95.
    f_bu_su : float, optional
        Fraction of butyrate from sugars [kg COD/kg COD]. The default is 0.13.
    f_pro_su : float, optional
        Fraction of propionate from sugars [kg COD/kg COD]. The default is 0.27.
    f_ac_su : float, optional
        Fraction of acetate from sugars [kg COD/kg COD]. The default is 0.41.
    f_va_aa : float, optional
        Fraction of valerate from amino acids [kg COD/kg COD]. The default is 0.23.
    f_bu_aa : float, optional
        Fraction of butyrate from amino acids [kg COD/kg COD]. The default is 0.26.
    f_pro_aa : float, optional
        Fraction of propionate from amino acids [kg COD/kg COD]. The default is 0.05.
    f_ac_aa : float, optional
        Fraction of acetate from amino acids [kg COD/kg COD]. The default is 0.4.
    f_ac_fa : float, optional
        Fraction of acetate from LCFAs [kg COD/kg COD]. The default is 0.7.
    f_pro_va : float, optional
        Fraction of propionate from LCFAs [kg COD/kg COD]. The default is 0.54.
    f_ac_va : float, optional
        Fraction of acetate from valerate [kg COD/kg COD]. The default is 0.31.
    f_ac_bu : float, optional
        Fraction of acetate from butyrate [kg COD/kg COD]. The default is 0.8.
    f_ac_pro : float, optional
        Fraction of acetate from propionate [kg COD/kg COD]. The default is 0.57.
    Y_su : float, optional
        Biomass yield of sugar uptake [kg COD/kg COD]. The default is 0.1.
    Y_aa : float, optional
        Biomass yield of amino acid uptake [kg COD/kg COD]. The default is 0.08.
    Y_fa : float, optional
        Biomass yield of LCFA uptake [kg COD/kg COD]. The default is 0.06.
    Y_c4 : float, optional
        Biomass yield of butyrate or valerate uptake [kg COD/kg COD]. 
        The default is 0.06.
    Y_pro : float, optional
        Biomass yield of propionate uptake [kg COD/kg COD]. The default is 0.04.
    Y_ac : float, optional
        Biomass yield of acetate uptake [kg COD/kg COD]. The default is 0.05.
    Y_h2 : float, optional
        Biomass yield of H2 uptake [kg COD/kg COD]. The default is 0.06.
    q_dis : float, optional
        Composites disintegration rate constant [d^(-1)]. The default is 0.5.
    q_ch_hyd : float, optional
        Carbohydrate hydrolysis rate constant [d^(-1)]. The default is 10.
    q_pr_hyd : float, optional
        Protein hydrolysis rate constant [d^(-1)]. The default is 10.
    q_li_hyd : float, optional
        Lipid hydrolysis rate constant [d^(-1)]. The default is 10.
    k_su : float, optional
        Sugar uptake rate constant [d^(-1)]. The default is 30.
    k_aa : float, optional
        Amino acid uptake rate constant [d^(-1)]. The default is 50.
    k_fa : float, optional
        LCFA uptake rate constant [d^(-1)]. The default is 6.
    k_c4 : float, optional
        Butyrate or valerate uptake rate constant [d^(-1)]. The default is 20.
    k_pro : float, optional
        Propionate uptake rate constant [d^(-1)]. The default is 13.
    k_ac : float, optional
        Acetate uptake rate constant [d^(-1)]. The default is 8.
    k_h2 : float, optional
        H2 uptake rate constant [d^(-1)]. The default is 35.
    K_su : float, optional
        Half saturation coefficient of sugar uptake [kg COD/m3]. 
        The default is 0.5.
    K_aa : float, optional
        Half saturation coefficient of amino acid uptake [kg COD/m3]. 
        The default is 0.3.
    K_fa : float, optional
        Half saturation coefficient of LCFA uptake [kg COD/m3]. 
        The default is 0.4.
    K_c4 : float, optional
        Half saturation coefficient of butyrate or valerate uptake [kg COD/m3]. 
        The default is 0.2.
    K_pro : float, optional
        Half saturation coefficient of propionate uptake [kg COD/m3]. 
        The default is 0.1.
    K_ac : float, optional
        Half saturation coefficient of acetate uptake [kg COD/m3]. 
        The default is 0.15.
    K_h2 : float, optional
        Half saturation coefficient of H2 uptake [kg COD/m3]. 
        The default is 7e-6.
    b_su : float, optional
        Decay rate constant of sugar-uptaking biomass [d^(-1)]. 
        The default is 0.02.
    b_aa : float, optional
        Decay rate constant of amino-acid-uptaking biomass [d^(-1)]. 
        The default is 0.02.
    b_fa : float, optional
        Decay rate constant of LCFA-uptaking biomass [d^(-1)].
        The default is 0.02.
    b_c4 : float, optional
        Decay rate constant of valerate- or butyrate-uptaking biomass [d^(-1)]. 
        The default is 0.02.
    b_pro : float, optional
        Decay rate constant of propionate-uptaking biomass [d^(-1)]. 
        The default is 0.02.
    b_ac : float, optional
        Decay rate constant of acetate-uptaking biomass [d^(-1)]. 
        The default is 0.02.
    b_h2 : float, optional
        Decay rate constant of H2-uptaking biomass [d^(-1)]. The default is 0.02.
    KI_h2_fa : float, optional
        H2 inhibition coefficient for LCFA uptake [kg COD/m3]. The default is 5e-6.
    KI_h2_c4 : float, optional
        H2 inhibition coefficient for butyrate or valerate uptake [kg COD/m3]. 
        The default is 1e-5.
    KI_h2_pro : float, optional
        H2 inhibition coefficient for propionate uptake [kg COD/m3]. 
        The default is 3.5e-6.
    KI_nh3 : float, optional
        Free ammonia inhibition coefficient for acetate uptake [M]. 
        The default is 1.8e-3.
    KS_IN : float, optional
        Inorganic nitrogen (nutrient) inhibition coefficient for soluble 
        substrate uptake [M]. The default is 1e-4.
    pH_limits_aa : 2-tuple, optional
        Lower and upper limits of pH inhibition for acidogens and acetogens, 
        unitless. The default is (4,5.5).
    pH_limits_ac : 2-tuple, optional
        Lower and upper limits of pH inhibition for aceticlastic methanogens, 
        unitless. The default is (6,7).
    pH_limits_h2 : 2-tuple, optional
        Lower and upper limits of pH inhibition for H2-utilizing methanogens, 
        unitless. The default is (5,6).
    T_base : float, optional
        Base temperature for kinetic parameters [K]. The default is 298.15.
    pKa_base : iterable[float], optional
        pKa (equilibrium coefficient) values of acid-base pairs at the base 
        temperature, unitless, following the order of `ADM1._acid_base_pairs`.
        The default is [14, 9.25, 6.35, 4.76, 4.88, 4.82, 4.86].
    Ka_dH : iterable[float], optional
        Heat of reaction of each acid-base pair at base temperature [J/mol], 
        following the order of `ADM1._acid_base_pairs`. The default is 
        [55900, 51965, 7646, 0, 0, 0, 0].
    kLa : float, optional
        Liquid-gas mass transfer coefficient [d^(-1)]. The default is 200.
    K_H_base : iterable[float], optional
        Henry's Law coefficients of biogas species at the base temperature 
        [M dissolved in liquid/bar]. Follows the order of `ADM1._biogas_IDs`. 
        The default is [7.8e-4, 1.4e-3, 3.5e-2].
    K_H_dH : iterable[float], optional
        Heat of reaction of liquid-gas transfer of biogas species [J/mol]. 
        Follows the order of `ADM1._biogas_IDs`. The default is 
        [-4180, -14240, -19410].
    
    Examples
    --------
    >>> from qsdsan import processes as pc
    >>> cmps = pc.create_adm1_cmps()
    >>> adm1 = pc.ADM1()
    >>> adm1.show()
    ADM1([disintegration, hydrolysis_carbs, hydrolysis_proteins, hydrolysis_lipids, uptake_sugars, uptake_amino_acids, uptake_LCFA, uptake_valerate, uptake_butyrate, uptake_propionate, uptake_acetate, uptake_h2, decay_Xsu, decay_Xaa, decay_Xfa, decay_Xc4, decay_Xpro, decay_Xac, decay_Xh2, h2_transfer, ch4_transfer, IC_transfer])
    
    References
    ----------
    .. [1] Batstone, D. J.; Keller, J.; Angelidaki, I.; Kalyuzhnyi, S. V; 
        Pavlostathis, S. G.; Rozzi, A.; Sanders, W. T. M.; Siegrist, H.; 
        Vavilin, V. A. The IWA Anaerobic Digestion Model No 1 (ADM1). 
        Water Sci. Technol. 2002, 45 (10), 65–73.
    .. [2] Rosen, C.; Jeppsson, U. Aspects on ADM1 Implementation within 
        the BSM2 Framework; Lund, 2006.
    """

    _stoichio_params = ('f_ch_xc', 'f_pr_xc', 'f_li_xc', 'f_xI_xc', 'f_sI_xc',
                        'f_fa_li', 'f_bu_su', 'f_pro_su', 'f_ac_su', 'f_h2_su',
                        'f_va_aa', 'f_bu_aa', 'f_pro_aa', 'f_ac_aa', 'f_h2_aa',
                        'f_ac_fa', 'f_h2_fa', 'f_pro_va', 'f_ac_va', 'f_h2_va',
                        'f_ac_bu', 'f_h2_bu', 'f_ac_pro', 'f_h2_pro',
                        'Y_su', 'Y_aa', 'Y_fa', 'Y_c4', 'Y_pro', 'Y_ac', 'Y_h2')
    _kinetic_params = ('rate_constants', 'half_sat_coeffs', 'pH_ULs', 'pH_LLs',
                       'KS_IN', 'KI_nh3', 'KIs_h2',
                       'Ka_base', 'Ka_dH', 'K_H_base', 'K_H_dH', 'kLa',
                       'T_base', 'components', 'root')
    _acid_base_pairs = (('H+', 'OH-'), ('NH4+', 'NH3'), ('CO2', 'HCO3-'),
                        ('HAc', 'Ac-'), ('HPr', 'Pr-'),
                        ('HBu', 'Bu-'), ('HVa', 'Va-'))
    _biogas_IDs = ('S_h2', 'S_ch4', 'S_IC')
    _biomass_IDs = ('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2')

    def __new__(cls, components=None, path=None, N_xc=2.686e-3, N_I=4.286e-3, N_aa=7e-3,
                f_ch_xc=0.2, f_pr_xc=0.2, f_li_xc=0.3, f_xI_xc=0.2,
                f_fa_li=0.95, f_bu_su=0.13, f_pro_su=0.27, f_ac_su=0.41,
                f_va_aa=0.23, f_bu_aa=0.26, f_pro_aa=0.05, f_ac_aa=0.4,
                f_ac_fa=0.7, f_pro_va=0.54, f_ac_va=0.31, f_ac_bu=0.8, f_ac_pro=0.57,
                Y_su=0.1, Y_aa=0.08, Y_fa=0.06, Y_c4=0.06, Y_pro=0.04, Y_ac=0.05, Y_h2=0.06,
                q_dis=0.5, q_ch_hyd=10, q_pr_hyd=10, q_li_hyd=10,
                k_su=30, k_aa=50, k_fa=6, k_c4=20, k_pro=13, k_ac=8, k_h2=35,
                K_su=0.5, K_aa=0.3, K_fa=0.4, K_c4=0.2, K_pro=0.1, K_ac=0.15, K_h2=7e-6,
                b_su=0.02, b_aa=0.02, b_fa=0.02, b_c4=0.02, b_pro=0.02, b_ac=0.02, b_h2=0.02,
                KI_h2_fa=5e-6, KI_h2_c4=1e-5, KI_h2_pro=3.5e-6, KI_nh3=1.8e-3, KS_IN=1e-4,
                pH_limits_aa=(4,5.5), pH_limits_ac=(6,7), pH_limits_h2=(5,6),
                T_base=298.15, pKa_base=[14, 9.25, 6.35, 4.76, 4.88, 4.82, 4.86],
                Ka_dH=[55900, 51965, 7646, 0, 0, 0, 0],
                kLa=200, K_H_base=[7.8e-4, 1.4e-3, 3.5e-2],
                K_H_dH=[-4180, -14240, -19410],
                **kwargs):
        
        cmps = _load_components(components)
        cmps.X_c.i_N = round(N_xc * N_mw, 4)
        cmps.X_I.i_N = cmps.S_I.i_N = round(N_I * N_mw, 4)
        cmps.S_aa.i_N = cmps.X_pr.i_N = round(N_aa * N_mw, 4)
        # cmps.X_c.i_N = N_xc * N_mw
        # cmps.X_I.i_N = cmps.S_I.i_N = N_I * N_mw
        # cmps.S_aa.i_N = cmps.X_pr.i_N = N_aa * N_mw

        if not path: path = _path
        self = Processes.load_from_file(path,
                                        components=cmps,
                                        conserved_for=('C', 'N'),
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

        stoichio_vals = (f_ch_xc, f_pr_xc, f_li_xc, f_xI_xc, 1.0-f_ch_xc-f_pr_xc-f_li_xc-f_xI_xc,
                         f_fa_li, f_bu_su, f_pro_su, f_ac_su, 1.0-f_bu_su-f_pro_su-f_ac_su,
                         f_va_aa, f_bu_aa, f_pro_aa, f_ac_aa, 1.0-f_va_aa-f_bu_aa-f_pro_aa-f_ac_aa,
                         f_ac_fa, 1.0-f_ac_fa, f_pro_va, f_ac_va, 1.0-f_pro_va-f_ac_va,
                         f_ac_bu, 1.0-f_ac_bu, f_ac_pro, 1.0-f_ac_pro,
                         Y_su, Y_aa, Y_fa, Y_c4, Y_pro, Y_ac, Y_h2)
        pH_LLs = np.array([pH_limits_aa[0]]*6 + [pH_limits_ac[0], pH_limits_h2[0]])
        pH_ULs = np.array([pH_limits_aa[1]]*6 + [pH_limits_ac[1], pH_limits_h2[1]])
        ks = np.array((q_dis, q_ch_hyd, q_pr_hyd, q_li_hyd,
                       k_su, k_aa, k_fa, k_c4, k_c4, k_pro, k_ac, k_h2,
                       b_su, b_aa, b_fa, b_c4, b_pro, b_ac, b_h2))
        Ks = np.array((K_su, K_aa, K_fa, K_c4, K_c4, K_pro, K_ac, K_h2))
        KIs_h2 = np.array((KI_h2_fa, KI_h2_c4, KI_h2_c4, KI_h2_pro))
        K_H_base = np.array(K_H_base)
        K_H_dH = np.array(K_H_dH)
        Ka_base = np.array([10**(-pKa) for pKa in pKa_base])
        Ka_dH = np.array(Ka_dH)
        root = TempState()
        dct = self.__dict__
        dct.update(kwargs)

        self.set_rate_function(rhos_adm1)
        dct['_parameters'] = dict(zip(cls._stoichio_params, stoichio_vals))
        self.rate_function._params = dict(zip(cls._kinetic_params,
                                              [ks, Ks, pH_ULs, pH_LLs, KS_IN*N_mw,
                                               KI_nh3, KIs_h2, Ka_base, Ka_dH,
                                               K_H_base, K_H_dH, kLa,
                                               T_base, self._components, root]))

        dct['flex_rhos'] = _rhos_adm1
        dct['solve_pH'] = solve_pH
        dct['dydt_Sh2_AD'] = dydt_Sh2_AD
        dct['grad_dydt_Sh2_AD'] = grad_dydt_Sh2_AD
        return self

    def set_pKas(self, pKas):
        '''Set the pKa values of the acid-base reactions at the base temperature.'''
        if len(pKas) != 7:
            raise ValueError(f'pKas must be an array of 7 elements, one for each '
                             f'acid-base pair, not {len(pKas)} elements.')
        dct = self.rate_function._params
        dct['Ka_base'] = np.array([10**(-pKa) for pKa in pKas])

    def _find_index(self, process):
        isa = isinstance
        if isa(process, int): return process
        elif isa(process, str): return self.index(process)
        elif isa(process, Process): return self.index(process.ID)
        else: raise TypeError(f'must input an int or str or :class:`Process`, '
                              f'not {type(process)}')

    def set_rate_constant(self, k, process):
        '''Set the reaction rate constant [d^(-1)] for a process given its ID.'''
        i = self._find_index(process)
        self.rate_function._params['rate_constants'][i] = k

    def set_half_sat_K(self, K, process):
        '''Set the substrate half saturation coefficient [kg/m3] for a process given its ID.'''
        i = self._find_index(process)
        self.rate_function._params['half_sat_coeffs'][i-4] = K

    def set_pH_inhibit_bounds(self, process, lower=None, upper=None):
        '''Set the upper and/or lower limit(s) of pH inhibition [unitless] for a process given its ID.'''
        i = self._find_index(process)
        dct = self.rate_function._params
        if lower is None: lower = dct['pH_LLs'][i-4]
        else: dct['pH_LLs'][i-4] = lower
        if upper is None: upper = dct['pH_ULs'][i-4]
        else: dct['pH_ULs'][i-4] = upper
        if lower >= upper:
            raise ValueError(f'lower limit for pH inhibition of {process} must '
                             f'be lower than the upper limit, not {[lower, upper]}')

    def set_h2_inhibit_K(self, KI, process):
        '''Set the H2 inhibition coefficient [kg/m3] for a process given its ID.'''
        i = self._find_index(process)
        self.rate_function._params['KIs_h2'][i-6] = KI

    def set_KS_IN(self, K):
        '''Set inhibition coefficient for inorganic nitrogen as a secondary
        substrate [M nitrogen].'''
        self.rate_function._params['KS_IN'] = K * N_mw

    def set_KI_nh3(self, K):
        '''Set inhibition coefficient for free ammonia [M].'''
        self.rate_function._params['KI_nh3'] = K

    def set_parameters(self, **parameters):
        '''Set values to stoichiometric parameters in `ADM1._stoichio_params`.'''
        non_stoichio = {}
        for k,v in parameters.items():
            if k in self._stoichio_params:
                if v >= 0 : self._parameters[k] = v
                else: raise ValueError(f'{k} must >= 0, not {v}')
            else: non_stoichio[k] = v
        if non_stoichio:
            warn(f'ignored value setting for non-stoichiometric parameters {non_stoichio}')
        self.check_stoichiometric_parameters()
        if self._stoichio_lambdified is not None:
            self.__dict__['_stoichio_lambdified'] = None

    def check_stoichiometric_parameters(self):
        '''Check whether product COD fractions sum up to 1 for each process.'''
        stoichio = self.parameters
        subst = ('xc', 'su', 'aa', 'fa', 'va', 'bu', 'pro')
        for s in subst:
            f_tot = sum([stoichio[k] for k in self._stoichio_params[:-7] \
                         if k.endswith(s)])
            if f_tot != 1:
                raise ValueError(f"the sum of 'f_()_{s}' values must equal 1")