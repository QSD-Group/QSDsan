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
import numpy as np
from qsdsan.utils import ospath, data_path
from scipy.optimize import brenth
from warnings import warn

__all__ = ('create_adm1_p_extension_cmps', 'ADM1_p_extension',
           'non_compet_inhibit', 'substr_inhibit',
           'T_correction_factor', 
           'pH_inhibit', 'Hill_inhibit', 
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
    cmps_all = Components.load_default()

    # varies
    # X_c = cmps_all.X_OHO.copy('X_c')
    # X_c.description = 'Composite'
    # X_c.i_C = 0.02786 * C_mw
    # X_c.i_N = 0.0376
    
    X_ch = Component.from_chemical('X_ch', chemical='glycogen', Tc=1011.4, # glucan
                                    description='Carbohydrates',
                                    measured_as='COD',
                                    particle_size='Particulate',
                                    degradability='Slowly',
                                    organic=True)
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
    
    S_IP = cmps_all.S_PO4.copy('S_IP')
  
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
    
    X_PHA = cmps_all.X_PAO_PHA.copy('X_PHA')
    
    X_PP = cmps_all.X_PAO_PP_Lo.copy('X_PP')
    
    X_PAO = cmps_all.X_PAO.copy('X_PAO')
    
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
    
    # X_MeOH and X_MeP are added at a later iteration of p_extension ADM1. 
    # They do not participate in any ADM1 process, and would be directly
    # mapped to ASM2d components at the ASM-ADM-ASM interface. 
    X_MeOH = cmps_all.X_FeOH.copy('X_MeOH')
    X_MeP = cmps_all.X_FePO4.copy('X_MeP')

    for bio in (X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2):
        # bio.formula = 'C5H7O2N'
        bio.i_C = 0.0313 * C_mw
        bio.i_N = 0.08

    S_cat = cmps_all.S_CAT.copy('S_cat')
    S_an = cmps_all.S_AN.copy('S_an')
    S_cat.i_mass = S_an.i_mass = 1

    cmps_adm1_p_extension = Components([S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2,
                            S_ch4, S_IC, S_IN, S_IP, S_I, X_ch, X_pr, X_li,
                            X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I,
                            X_PHA, X_PP, X_PAO, S_K, S_Mg, X_MeOH, X_MeP, 
                            S_cat, S_an, cmps_all.H2O])
    cmps_adm1_p_extension.default_compile()
    if set_thermo: settings.set_thermo(cmps_adm1_p_extension)
    return cmps_adm1_p_extension

# create_adm1_p_extension_cmps()


#%%
# =============================================================================
# kinetic rate functions
# =============================================================================

R = 8.3145e-2 # Universal gas constant, [bar/M/K]

def non_compet_inhibit(Si, Ki):
    return Ki/(Ki+Si)

def substr_inhibit(Si, Ki):
    return Si/(Ki+Si)

def mass2mol_conversion(cmps):
    '''conversion factor from kg[measured_as]/m3 to mol[component]/L'''
    return cmps.i_mass / cmps.chem_MW

# def T_correction_factor(T1, T2, theta):
#     return np.exp(theta * (T2-T1))

def T_correction_factor(T1, T2, delta_H):
    return np.exp(delta_H/(R*100) * (1/T1 - 1/T2))  # R converted to SI

# def calc_Kas(pKas, T_base, T_op, theta):
#     pKas = np.asarray(pKas)
#     return 10**(-pKas) * T_correction_factor(T_base, T_op, theta)

def acid_base_rxn(h_ion, weak_acids_tot, Kas):
    # h, nh4, hco3, ac, pr, bu, va = mols
    # S_cat, S_K, S_Mg, S_an, S_IN, S_IP, S_IC, S_ac, S_pro, S_bu, S_va = weak_acids_tot  # in M
    S_cat, S_K, S_Mg, S_an, S_IN, S_IP = weak_acids_tot[:6]
    # Kw, Ka_nh, Ka_h3po4, Ka_co2, Ka_ac, Ka_pr, Ka_bu, Ka_va = Kas
    Kw = Kas[0]
    oh_ion = Kw/h_ion
    nh3, h2po4, hco3, ac, pro, bu, va = Kas[1:] * weak_acids_tot[4:] / (Kas[1:] + h_ion)
    return S_cat + S_K + 2*S_Mg + h_ion + (S_IN - nh3) - S_an - oh_ion - hco3 - ac - pro - bu - va - (3*S_IP + h2po4)

# The function 'fprime_abr' is not used in the code
def fprime_abr(h_ion, weak_acids_tot, Kas):
    S_cat, S_K, S_Mg, S_an, S_IN, S_IP = weak_acids_tot[:6]
    Kw = Kas[0]
    doh_ion = - Kw / h_ion ** 2
    dnh3, dh2po4, dhco3, dac, dpro, dbu, dva = - Kas[1:] * weak_acids_tot[4:] / (Kas[1:] + h_ion)**2
    return 1 + (-dnh3) - dh2po4 - doh_ion - dhco3 - dhco3 - dac - dpro - dbu - dva

def pH_inhibit(pH, ul, ll, lower_only=True):
    if lower_only:
        # if pH >= ul: return 1
        # else: return exp(-3 * ((pH-ul)/(ul-ll))**2)
        low_by = np.minimum(pH-ul, 0)
        return np.exp(-3 * (low_by/(ul-ll))**2)
    else:
        return (1+2*10**(0.5*(ll-ul)))/(1+10**(pH-ul)+10**(ll-pH))

def Hill_inhibit(H_ion, ul, ll):
    n = 3/(ul-ll)
    K = 10**(-(ul+ll)/2)
    return 1/(1+(H_ion/K) ** n)

rhos = np.zeros(28) # 28 kinetic processes (25 as defined in modified ADM1 + 3 for gases)
Cs = np.empty(25) # 25 processes as defined in modified ADM1

def rhos_adm1_p_extension(state_arr, params):
    ks = params['rate_constants']
    Ks = params['half_sat_coeffs']
    cmps = params['components']
    # n = len(cmps)
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
    
    # state_arr_cmps stated just for readability of code 
    # state_arr_cmps = [S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2, S_ch4, S_IC, S_IN,
    #             S_IP, S_I, X_ch, X_pr, X_li, X_su, X_aa, X_fa, X_c4, X_pro, X_ac, 
    #             X_h2, X_I, X_PHA, X_PP, X_PAO, S_K, S_Mg, X_MeOH, X_MeP]
    
    # Cs_ids = cmps.indices(['X_ch', 'X_pr', 'X_li', 'X_su', 'X_aa',
    #                        'X_fa', 'X_c4', 'X_c4', 'X_pro', 'X_ac', 'X_h2',
    #                        'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2',
    #                        'X_PAO', 'X_PAO', 'X_PAO', 'X_PAO', 'X_PAO', 'X_PP', 'X_PHA'])
    # Cs = state_arr[Cs_ids]
    Cs[:7] = state_arr[13:20]
    Cs[7:11] = state_arr[19:23]
    Cs[11:18] = state_arr[16:23]
    Cs[18:23] = state_arr[26]
    Cs[23] = state_arr[25]
    Cs[24] = state_arr[24]
    
    # substrates_ids = cmps.indices(['S_su', 'S_aa', 'S_fa', 'S_va',
    #                                'S_bu', 'S_pro', 'S_ac', 'S_h2'])
    # substrates = state_arr[substrates_ids]
    substrates = state_arr[:8]
    
    # S_va, S_bu, S_h2, S_IN = state_arr[cmps.indices(['S_va', 'S_bu', 'S_h2', 'S_IN'])]
    # S_va, S_bu, S_h2, S_ch4, S_IC, S_IN = state_arr[[3,4,7,8,9,10]]
    S_va, S_bu, S_h2, S_IN, S_IP = state_arr[[3,4,7,10,11]]
    unit_conversion = mass2mol_conversion(cmps)
    # Should I let X_MeOH, and X_MeP unit convert? (ask Joy)
    cmps_in_M = state_arr[:34] * unit_conversion 
    # weak acids (ADM1) = [S_ca, S_an, S_IN, S_IC, S_ac, S_pro, S_bu, S_va]
    # weak acids (modified_ADM1) = [S_ca, S_K, S_Mg, S_an, S_IN, S_IP, S_IC, S_ac, S_pro, S_bu, S_va]
    weak_acids = cmps_in_M[[31, 27, 28, 32, 10, 11, 9, 6, 5, 4, 3]]

    T_op = state_arr[-1]
    biogas_S = state_arr[7:10].copy()
    biogas_p = R * T_op * state_arr[34:37] 
    Kas = Kab * T_correction_factor(T_base, T_op, Ka_dH)
    KH = KHb * T_correction_factor(T_base, T_op, KH_dH) / unit_conversion[7:10]

    rhos[:-3] = ks * Cs
    rhos[3:11] *= substr_inhibit(substrates, Ks)
    if S_va > 0: rhos[6] *= 1/(1+S_bu/S_va)
    if S_bu > 0: rhos[7] *= 1/(1+S_va/S_bu)
    
    # substrates_ids = cmps.indices(['S_va', 'S_bu', 'S_pro', 'S_ac'])
    # substrates_modified = state_arr[substrates_ids]
    substrates_modified = state_arr[3:7]
    # Ka = K_a 
    Ka = Ks[-2]
    rhos[18:22] *= substr_inhibit(substrates_modified, Ka)
    
    # PP_PAO = X_PP/X_PAO
    PP_PAO = state_arr[25]/state_arr[26]
    # Kpp = K_pp
    Kpp = Ks[-1]
    rhos[18:22] *= substr_inhibit(PP_PAO, Kpp)
    
    # Multiplication by {Sva, Sbu, Spro, Sac}/(Sva + Sbu + Spro + Sac)
    transformation_array = state_arr[3:7]/sum(state_arr[3:7])
    rhos[18:22] *= transformation_array
    
    h = brenth(acid_base_rxn, 1e-14, 1.0,
               args=(weak_acids, Kas),
               xtol=1e-12, maxiter=100)

    nh3 = Kas[1] * weak_acids[4] / (Kas[1] + h)
    co2 = weak_acids[6] - Kas[2] * weak_acids[6] / (Kas[2] + h)
    biogas_S[-1] = co2 / unit_conversion[9]
    
    Iph = Hill_inhibit(h, pH_ULs, pH_LLs)
    Ih2 = non_compet_inhibit(S_h2, KIs_h2)
    root.data = [-np.log10(h), Iph, Ih2]
    rhos[3:11] *= Iph * substr_inhibit(S_IN, KS_IN) * substr_inhibit(S_IP, KS_IP)
    rhos[5:9] *= Ih2
    # rhos[4:12] *= Hill_inhibit(h, pH_ULs, pH_LLs) * substr_inhibit(S_IN, KS_IN)
    # rhos[6:10] *= non_compet_inhibit(S_h2, KIs_h2)
    rhos[9] *= non_compet_inhibit(nh3, KI_nh3)
    rhos[-3:] = kLa * (biogas_S - KH * biogas_p)
    # print(rhos)
    return rhos
#%%
# =============================================================================
# ADM1_p_extension class
# =============================================================================
class TempState:
    def __init__(self):
        self.data = []
    
    # def append(self, value):
    #     self.data += [value]

@chemicals_user
class ADM1_p_extension(CompiledProcesses):
    """
    Anaerobic Digestion Model No.1. [1]_, [2]_, [3]_

    Parameters
    ----------
    components : class:`CompiledComponents`, optional
        Components corresponding to each entry in the stoichiometry array,
        defaults to thermosteam.settings.chemicals.
    path : str, optional
        Alternative file path for the Petersen matrix. The default is None.
    N_I : float, optional
        Nitrogen content of inert organics [kmol N/kg COD]. The default is 4.286e-3.
    N_aa : float, optional
        Nitrogen content of amino acids [kmol N/kg COD]. The default is 7e-3.
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
    f_ac_PHA : float, optional
        Yield of acetate on PHA [kg COD/kg COD]. The default is 0.4.
    f_bu_PHA : float, optional
        Yield of butyrate on PHA [kg COD/kg COD]. The default is 0.1.
    f_pro_PHA : float, optional
        Yield of propionate on PHA [kg COD/kg COD]. The default is 0.4.
    f_va_PHA : float, optional
        Yield of valerate on PHA [kg COD/kg COD]. The default is 0.1.
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
    Y_po4 : float, optional
        Yield of biomass on phosphate [kmol P/kg COD]. The default is 0.013. 
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
    K_a : float, optional
        Saturation coefficient for acetate [kg COD/m3].
        The default is 0.004.
    K_pp : float, optional
        Saturation coefficient for polyphosphate [kg COD/m3].
        The default is 0.00032.
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
    q_PHA : float, optional
        Rate constant for storage of PHA [d^(-1)]. The default is 3.
    b_PAO : float, optional
        Lysis rate of PAOs [d^(-1)]. The default is 0.2.
    b_PP : float, optional
        Lysis rate of polyphosphates [d^(-1)]. The default is 0.2.
    b_PHA : float, optional
        Lysis rate of PHAs [d^(-1)]. The default is 0.2.
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
    KS_IP : float, optional
        P limitation for inorganic phosphorous [kmol P/m3]. 
        The default is 2e-5. 
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
        The default is [14, 9.25, 2.12, 6.35, 4.76, 4.88, 4.82, 4.86].
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

    _stoichio_params = ('f_sI_xb', 'f_ch_xb', 'f_pr_xb', 'f_li_xb', 'f_xI_xb', 
                        'f_fa_li', 'f_bu_su', 'f_pro_su', 'f_ac_su', 'f_h2_su',
                        'f_va_aa', 'f_bu_aa', 'f_pro_aa', 'f_ac_aa', 'f_h2_aa',
                        'f_ac_fa', 'f_h2_fa', 'f_pro_va', 'f_ac_va', 'f_h2_va',
                        'f_ac_bu', 'f_h2_bu', 'f_ac_pro', 'f_h2_pro',
                        'f_ac_PHA', 'f_bu_PHA', 'f_pro_PHA', 'f_va_PHA', 
                        'Y_su', 'Y_aa', 'Y_fa', 'Y_c4', 'Y_pro', 'Y_ac', 'Y_h2', 'Y_po4')
    _kinetic_params = ('rate_constants', 'half_sat_coeffs', 'pH_ULs', 'pH_LLs',
                       'KS_IN', 'KS_IP', 'KI_nh3', 'KIs_h2',
                       'Ka_base', 'Ka_dH', 'K_H_base', 'K_H_dH', 'kLa',
                       'T_base', 'components', 'root')
    _acid_base_pairs = (('H+', 'OH-'), ('NH4+', 'NH3'), ('H3PO4', 'H2PO4-'), 
                        ('CO2', 'HCO3-'), ('HAc', 'Ac-'), ('HPr', 'Pr-'),
                        ('HBu', 'Bu-'), ('HVa', 'Va-'))
    
    _biogas_IDs = ('S_h2', 'S_ch4', 'S_IC')

    def __new__(cls, components=None, path=None, N_xc=2.686e-3, N_I=4.286e-3, N_aa=7e-3,
                f_sI_xb=0, f_ch_xb=0.275, f_pr_xb=0.275, f_li_xb=0.350,
                f_fa_li=0.95, f_bu_su=0.13, f_pro_su=0.27, f_ac_su=0.41,
                f_va_aa=0.23, f_bu_aa=0.26, f_pro_aa=0.05, f_ac_aa=0.4,
                f_ac_fa=0.7, f_pro_va=0.54, f_ac_va=0.31, f_ac_bu=0.8, f_ac_pro=0.57,
                f_ac_PHA=0.4, f_bu_PHA=0.1, f_pro_PHA=0.4,
                Y_su=0.1, Y_aa=0.08, Y_fa=0.06, Y_c4=0.06, Y_pro=0.04, Y_ac=0.05, Y_h2=0.06, Y_po4=0.013,
                q_dis=0.5, q_ch_hyd=10, q_pr_hyd=10, q_li_hyd=10,
                k_su=30, k_aa=50, k_fa=6, k_c4=20, k_pro=13, k_ac=8, k_h2=35,
                K_su=0.5, K_aa=0.3, K_fa=0.4, K_c4=0.2, K_pro=0.1, K_ac=0.15, K_h2=7e-6, K_a=4e-3, K_pp=32e-5,
                b_su=0.02, b_aa=0.02, b_fa=0.02, b_c4=0.02, b_pro=0.02, b_ac=0.02, b_h2=0.02,
                q_PHA=3, b_PAO=0.2, b_PP=0.2, b_PHA=0.2, 
                KI_h2_fa=5e-6, KI_h2_c4=1e-5, KI_h2_pro=3.5e-6, KI_nh3=1.8e-3, KS_IN=1e-4, KS_IP=2e-5, 
                pH_limits_aa=(4,5.5), pH_limits_ac=(6,7), pH_limits_h2=(5,6),
                T_base=298.15, pKa_base=[14, 9.25, 2.12, 6.35, 4.76, 4.88, 4.82, 4.86],
                Ka_dH=[55900, 51965, 7646, 0, 0, 0, 0],
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
        
        # Need to add 5 processes separately for stoichiometry. See ASM2d for reference. 
        if path == _path:
            _p19 = Process('storage_Sva_in_XPHA',
                           'S_va + [Y_po4]X_PP -> X_PHA + [?]S_K + [?]S_Mg',
                           components=cmps,
                           ref_component='X_PHA',
                           rate_equation='q_PHA * S_va/(K_a+S_va) * (X_PP/X_PAO)/(K_PP+(X_PP/X_PAO)) * X_PAO * S_va/(S_va+S_bu+S_pro+S_ac)',
                           parameters=('Y_po4', 'q_PHA', 'K_a', 'K_PP'),
                           conserved_for=('K', 'Mg'))
            
            _p20 = Process('storage_Sbu_in_XPHA',
                           'S_bu + [Y_po4]X_PP -> X_PHA + [?]S_K + [?]S_Mg',
                           components=cmps,
                           ref_component='X_PHA',
                           rate_equation='q_PHA * S_bu/(K_a+S_bu) * (X_PP/X_PAO)/(K_PP+(X_PP/X_PAO)) * X_PAO * S_bu/(S_va+S_bu+S_pro+S_ac)',
                           parameters=('Y_po4', 'q_PHA', 'K_a', 'K_PP'),
                           conserved_for=('K', 'Mg'))
            
            _p21 = Process('storage_Spro_in_XPHA',
                           'S_pro + [Y_po4]X_PP -> X_PHA + [?]S_K + [?]S_Mg',
                           components=cmps,
                           ref_component='X_PHA',
                           rate_equation='q_PHA * S_pro/(K_a+S_pro) * (X_PP/X_PAO)/(K_PP+(X_PP/X_PAO)) * X_PAO * S_pro/(S_va+S_bu+S_pro+S_ac)',
                           parameters=('Y_po4', 'q_PHA', 'K_a', 'K_PP'),
                           conserved_for=('K', 'Mg'))
            
            _p22 = Process('storage_Sac_in_XPHA',
                           'S_ac + [Y_po4]X_PP -> X_PHA + [?]S_K + [?]S_Mg',
                           components=cmps,
                           ref_component='X_PHA',
                           rate_equation='q_PHA * S_ac/(K_a+S_ac) * (X_PP/X_PAO)/(K_PP+(X_PP/X_PAO)) * X_PAO * S_ac/(S_va+S_bu+S_pro+S_ac)',
                           parameters=('Y_po4', 'q_PHA', 'K_a', 'K_PP'),
                           conserved_for=('K', 'Mg'))
            
            _p24 = Process('lysis_XPP',
                           'X_PP -> [?]S_K + [?]S_Mg',
                           components=cmps,
                           ref_component='X_PP',
                           rate_equation='b_PP*X_PP',
                           parameters=('b_PP'),
                           conserved_for=('K', 'Mg'))
            
            self.insert(18, _p19)
            self.insert(19, _p20)
            self.insert(20, _p21)
            self.insert(21, _p22)
            self.insert(23, _p24)
            
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

        stoichio_vals = (f_sI_xb, f_ch_xb, f_pr_xb, f_li_xb, 1-f_sI_xb-f_ch_xb-f_pr_xb-f_li_xb,
                         f_fa_li, f_bu_su, f_pro_su, f_ac_su, 1-f_bu_su-f_pro_su-f_ac_su,
                         f_va_aa, f_bu_aa, f_pro_aa, f_ac_aa, 1-f_va_aa-f_bu_aa-f_pro_aa-f_ac_aa,
                         f_ac_fa, 1-f_ac_fa, f_pro_va, f_ac_va, 1-f_pro_va-f_ac_va,
                         f_ac_bu, 1-f_ac_bu, f_ac_pro, 1-f_ac_pro,
                         f_ac_PHA, f_bu_PHA, f_pro_PHA, 1-f_ac_PHA-f_bu_PHA-f_pro_PHA,
                         Y_su, Y_aa, Y_fa, Y_c4, Y_pro, Y_ac, Y_h2, Y_po4)
        pH_LLs = np.array([pH_limits_aa[0]]*6 + [pH_limits_ac[0], pH_limits_h2[0]])
        pH_ULs = np.array([pH_limits_aa[1]]*6 + [pH_limits_ac[1], pH_limits_h2[1]])
        
        ks = np.array((q_ch_hyd, q_pr_hyd, q_li_hyd,
                       k_su, k_aa, k_fa, k_c4, k_c4, k_pro, k_ac, k_h2,
                       b_su, b_aa, b_fa, b_c4, b_pro, b_ac, b_h2,
                       q_PHA, q_PHA, q_PHA, q_PHA, b_PAO, b_PP, b_PHA))
        
        Ks = np.array((K_su, K_aa, K_fa, K_c4, K_c4, K_pro, K_ac, K_h2, K_a, K_pp))
        KIs_h2 = np.array((KI_h2_fa, KI_h2_c4, KI_h2_c4, KI_h2_pro))
        K_H_base = np.array(K_H_base)
        K_H_dH = np.array(K_H_dH)
        Ka_base = np.array([10**(-pKa) for pKa in pKa_base])
        Ka_dH = np.array(Ka_dH)
        root = TempState()
        # root.data = 10**(-7.4655)
        dct = self.__dict__
        dct.update(kwargs)

        self.set_rate_function(rhos_adm1_p_extension)
        dct['_parameters'] = dict(zip(cls._stoichio_params, stoichio_vals))
        self.rate_function._params = dict(zip(cls._kinetic_params,
                                              [ks, Ks, pH_ULs, pH_LLs, KS_IN*N_mw, KS_IP*P_mw, 
                                               KI_nh3, KIs_h2, Ka_base, Ka_dH,
                                               K_H_base, K_H_dH, kLa,
                                               T_base, self._components, root]))

        return self

    def set_pKas(self, pKas):
        '''Set the pKa values of the acid-base reactions at the base temperature.'''
        if len(pKas) != 8:
            raise ValueError(f'pKas must be an array of 8 elements, one for each '
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
        
    def set_KS_IP(self, K):
        '''Set inhibition coefficient for inorganic phosphorous as a secondary
        substrate [M phosphorous].'''
        self.rate_function._params['KS_IP'] = K * P_mw

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
        subst = ('xb', 'su', 'aa', 'fa', 'va', 'bu', 'pro', 'PHA')
        for s in subst:
            f_tot = sum([stoichio[k] for k in self._stoichio_params[:-7] \
                         if k.endswith(s)])
            if f_tot != 1:
                raise ValueError(f"the sum of 'f_()_{s}' values must equal 1")
