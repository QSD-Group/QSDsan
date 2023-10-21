# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

from thermosteam.utils import chemicals_user
from thermosteam import settings
from chemicals.elements import molecular_weight as get_mw
from qsdsan import Component, Components, Process, Processes, CompiledProcesses
import numpy as np, qsdsan.processes as pc, qsdsan as qs
from qsdsan.utils import ospath, data_path
from qsdsan.processes._adm1 import (
    R,
    create_adm1_cmps, 
    ADM1,
    mass2mol_conversion,
    T_correction_factor,
    substr_inhibit,
    non_compet_inhibit,
    Hill_inhibit
    )
# from scipy.optimize import brenth
# from warnings import warn


__all__ = ('create_madm1_cmps', 'ModifiedADM1')

_path = ospath.join(data_path, 'process_data/_madm1.tsv')

#%% components
# C_mw = get_mw({'C':1})
N_mw = get_mw({'N':1})
P_mw = get_mw({'P':1})
S_mw = get_mw({'S':1})
Fe_mw = get_mw({'Fe':1})
O_mw = get_mw({'O':1})

def create_madm1_cmps(set_thermo=True, ASF_L=0.31, ASF_H=1.2):
    '''
    Create a set of components for the modified ADM1.

    Parameters
    ----------
    set_thermo : bool, optional
        Whether to set thermo with the returned set of components. The default is True.
    ASF_L : float, optional
        Active site factor for X_HFO_L [mol P sites/mol Fe]. The default is 0.31.
    ASF_H : float, optional
        Active site factor for X_HFO_H [mol P sites/mol Fe]. The default is 1.2.

    Returns
    -------
    cmps_madm1 : class:`CompiledComponents`

    '''
    
    # Components from the original ADM1
    # *********************************
    _cmps = create_adm1_cmps(False)
    S_aa = _cmps.S_aa
    X_pr = _cmps.X_pr
    S_aa.i_C = X_pr.i_C = 0.36890
    S_aa.i_N = X_pr.i_N = 0.11065
    S_aa.i_P = X_pr.i_P = 0.
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
    
    ion_properties = dict(
        particle_size='Soluble',
        degradability='Undegradable',
        organic=False)
    S_K = Component.from_chemical('S_K', chemical='K+', description='Potassium ion', 
                                  measured_as='K', **ion_properties)
    S_Mg = Component.from_chemical('S_Mg', chemical='Mg2+', description='Magnesium ion',
                                  measured_as='Mg',**ion_properties)
    S_SO4 = Component.from_chemical('S_SO4', chemical='SO4-2', description='Sulfate',
                                  measured_as='S', **ion_properties)
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
                                  description='Elemental sulfur',
                                  measured_as='COD',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
    S_Fe3 = Component.from_chemical('S_Fe3', chemical='Fe3+', description='Iron (III)',
                                  measured_as='Fe',**ion_properties)
    S_Fe2 = Component.from_chemical('S_Fe2', chemical='Fe2+', description='Iron (II)',
                                  measured_as='Fe',**ion_properties)
    S_Fe2.i_COD = 0.5*O_mw/Fe_mw
    S_Fe2.measured_as = 'COD'
    
    # Multiple mineral precipitation
    # ******************************
    mineral_properties = dict(
        particle_size='Particulate',
        degradability='Undegradable',
        organic=False)
    
    X_HFO_H = Component('X_HFO_H', formula='FeO(OH)',
                        description='Hydrous ferric oxide with high number of active sites',
                        measured_as='Fe',**mineral_properties)
    X_HFO_L = X_HFO_H.copy('X_HFO_L')
    X_HFO_L.description = 'Hydrous ferric oxide with low number of active sites'
    
    X_HFO_old = X_HFO_H.copy('X_HFO_old')
    X_HFO_old.description = 'Inactive hydrous ferric oxide'
    
    X_HFO_HP = Component('X_HFO_HP', formula=f'FeO(OH)P{ASF_H}',
                         description='X_HFO_H with phosphorus-bounded adsorption sites',
                         measured_as='Fe', **mineral_properties)
    X_HFO_HP_old = X_HFO_HP.copy('X_HFO_HP_old')
    X_HFO_HP_old.description = 'Old ' + X_HFO_HP.description
    
    X_HFO_LP = Component('X_HFO_LP', formula=f'FeO(OH)P{ASF_L}',
                         description='X_HFO_L with phosphorus-bounded adsorption sites',
                         measured_as='Fe', **mineral_properties)
    X_HFO_LP_old = X_HFO_LP.copy('X_HFO_LP_old')
    X_HFO_LP_old.description = 'Old ' + X_HFO_LP.description
    
    X_CCM = Component.from_chemical('X_CCM', chemical='calcite', description='Calcite', **mineral_properties)
    X_ACC = Component.from_chemical('X_ACC', chemical='aragonite', description='Aragonite', **mineral_properties)
    X_ACP = Component.from_chemical('X_ACP', chemical='Ca3(PO4)2', description='Amorphous calcium phosphate', **mineral_properties)
    X_HAP = Component.from_chemical('X_HAP', chemical='hydroxylapatite', description='Hydroxylapatite', **mineral_properties)
    X_DCPD = Component.from_chemical('X_DCPD', chemical='CaHPO4', description='Dicalcium phosphate', **mineral_properties)
    X_OCP = Component('X_OCP', formula='Ca4HP3O12', description='Octacalcium phosphate', **mineral_properties)
    X_struv = Component.from_chemical('X_struv', chemical='MgNH4PO4', description='Struvite', **mineral_properties)
    X_newb = Component.from_chemical('X_newb', chemical='MgHPO4', description='Newberyite', **mineral_properties)
    X_magn = Component.from_chemical('X_magn', chemical='MgCO3', description='Magnesite', **mineral_properties)
    X_kstruv = Component('X_kstruv', formula='MgKPO4', description='K-struvite', **mineral_properties)
    X_FeS = Component.from_chemical('X_FeS', chemical='FeS', description='Iron sulfide', **mineral_properties)
    X_Fe3PO42 = Component('X_Fe3PO42', formula='Fe3(PO4)2', description='Ferrous phosphate', **mineral_properties)
    X_AlPO4 = Component.from_chemical('X_AlPO4', chemical='AlPO4', description='Aluminum phosphate', **mineral_properties)

    S_Ca = Component.from_chemical('S_Ca', chemical='Ca2+', description='Calsium ion',
                                   measured_as='Ca', **ion_properties)
    S_Al = Component.from_chemical('S_Al', chemical='Al3+', description='Aluminum ion',
                                   measured_as='Al', **ion_properties)
    S_Na = Component.from_chemical('S_Na', chemical='Na+', description='Sodium ion',
                                   measured_as='Na', **ion_properties)
    S_Cl = Component.from_chemical('S_Cl', chemical='Cl-', description='Chloride',
                                   measured_as='Cl', **ion_properties)
    
    cmps_madm1 = Components([_cmps.S_su, S_aa, S_fa, _cmps.S_va, S_bu, 
                             S_pro, S_ac, _cmps.S_h2, _cmps.S_ch4, 
                             _cmps.S_IC, _cmps.S_IN, S_IP, S_I, 
                             X_ch, X_pr, X_li, *adm1_biomass, X_I,
                             X_PHA, asm_cmps.X_PP, X_PAO, S_K, S_Mg, 
                             S_SO4, S_IS, X_hSRB, X_aSRB, X_pSRB, X_c4SRB,
                             S_S0, S_Fe3, S_Fe2, X_HFO_H, X_HFO_L, X_HFO_old,
                             X_HFO_HP, X_HFO_LP, X_HFO_HP_old, X_HFO_LP_old,
                             S_Ca, S_Al, X_CCM, X_ACC, X_ACP, X_HAP, X_DCPD,
                             X_OCP, X_struv, X_newb, X_magn, X_kstruv, X_FeS,
                             X_Fe3PO42, X_AlPO4, 
                             S_Na, S_Cl, _cmps.H2O])
    cmps_madm1.default_compile()
    
    if set_thermo: qs.set_thermo(cmps_madm1)
    return cmps_madm1

#%% rate functions

# https://wiki.dynamita.com/en/biokinetic_process_models#chemical-phosphorus-removal-with-metal-salts-addition-iron-or-aluminium

# =============================================================================
# state_variable_indices = {
#     'S_su': 0, 'S_aa': 1, 'S_fa': 2, 'S_va': 3, 'S_bu': 4, 'S_pro': 5, 'S_ac': 6, 'S_h2': 7,
#     'S_ch4': 8, 'S_IC': 9, 'S_IN': 10, 'S_IP': 11, 'S_I': 12,
#     'X_ch': 13, 'X_pr': 14, 'X_li': 15, 
#     'X_su': 16, 'X_aa': 17, 'X_fa': 18, 'X_c4': 19, 'X_pro': 20, 'X_ac': 21, 'X_h2': 22, 'X_I': 23,
#     'X_PHA': 24, 'X_PP': 25, 'X_PAO': 26, 'S_K': 27, 'S_Mg': 28,
#     'S_SO4': 29, 'S_IS': 30, 'X_hSRB': 31, 'X_aSRB': 32, 'X_pSRB': 33, 'X_c4SRB': 34,
#     'S_S0': 35, 'S_Fe3': 36, 'S_Fe2': 37,
#     'X_HFO_H': 38, 'X_HFO_L': 39, 'X_HFO_old': 40, 'X_HFO_HP': 41, 'X_HFO_LP': 42, 'X_HFO_HP_old': 43, 'X_HFO_LP_old': 44,
#     'S_Ca': 45, 'S_Al': 46,
#     'X_CCM': 47, 'X_ACC': 48, 'X_ACP': 49, 'X_HAP': 50, 'X_DCPD': 51, 'X_OCP': 52,
#     'X_struv': 53, 'X_newb': 54, 'X_magn': 55, 'X_kstruv': 56, 
#     'X_FeS': 57, 'X_Fe3PO42': 58,
#     'X_AlPO4': 59,
#     'S_Na': 60, 'S_Cl': 61, 'H2O': 62
#     }
# =============================================================================

def calc_pH():
    pass

def calc_biogas():
    pass

def pcm():
    pass

def saturation_index():
    pass

rhos = np.zeros(38+8+13+4) # 38 biological + 8 chemical P removal by HFO + 13 MMP + 4 gas transfer
Cs = np.empty(38+8)
sum_stoichios = np.array([2, 2, 5, 9, 3, 8, 3, 3, 2, 3, 2, 2])

def rhos_madm1(state_arr, params, T_op):
    ks = params['rate_constants']
    Ks = params['half_sat_coeffs']
    K_PP = params['K_PP']
    K_so4 = params['K_so4']
    cmps = params['components']
    # n = len(cmps)
    pH_LLs, pH_ULs = params['pH_limits']
    KS_IN = params['KS_IN']
    KS_IP = params['KS_IP']
    KI_nh3 = params['KI_nh3']
    KIs_h2 = params['KIs_h2']
    KIs_h2s = params['KIs_h2s']
    KHb = params['K_H_base']
    Kab = params['Ka_base']
    KH_dH = params['K_H_dH']
    Ka_dH = params['Ka_dH']
    kLa = params['kLa']
    k_cryst = params['k_cryst']
    n_cryst = params['n_cryst']
    Kspb = params['Ksp_base']
    Ksp_dH = params['Ksp_dH']
    T_base = params['T_base']
    
    Cs[:7] = state_arr[13:20]                   # original ADM1 processes
    Cs[7:11] = state_arr[19:23]
    Cs[11:18] = state_arr[16:23]
    Cs[18:23] = X_PAO = state_arr[26]           # P extension processes
    Cs[23:25] = X_PP, X_PHA = state_arr[[25,24]]
    Cs[25:27] = state_arr[31]                   # S extension processes
    Cs[27:29] = state_arr[32]
    Cs[29:31] = state_arr[33]
    Cs[31:34] = state_arr[34]
    Cs[34:36] = Cs[36:38] = Cs[38:40] = Cs[40:42] = state_arr[38:40]   # Fe extension processes + HFO module
    Cs[42:44] = Cs[44:46] = state_arr[41:43]
    
    rhos[:46] = ks * Cs
    primary_substrates = state_arr[:8]
    
    rhos[3:11] *= substr_inhibit(primary_substrates, Ks[:8])
    c4 = primary_substrates[[3,4]]
    if sum(c4) > 0: rhos[[6,7]] *= c4/sum(c4)
    
    vfas = primary_substrates[3:7]
    rhos[18:22] *= substr_inhibit(vfas, Ks[8])
    if sum(vfas) > 0: rhos[18:22] *= vfas/sum(vfas)
    if X_PAO > 0: rhos[18:22] *= substr_inhibit(X_PP/X_PAO, K_PP)
    
    srb_subs = np.flip(primary_substrates[3:])
    S_SO4, S_IS = state_arr[29:31]
    rhos[[25,27,29,31,32]] *= substr_inhibit(srb_subs, Ks[9:13]) * substr_inhibit(S_SO4, K_so4)
    if sum(srb_subs[-2:]) > 0: rhos[[31,32]] *= srb_subs[-2:]/sum(srb_subs[-2:])
    
    #!!! why divide by 16 or 64?
    S_h2 = primary_substrates[-1]
    rhos[34:36] *= S_h2 / 16 
    rhos[36:38] *= S_IS / 64
    
    KPbind, KPdiss = Ks[-2:]
    S_IP = state_arr[11]
    rhos[40:42] *= substr_inhibit(S_IP, KPbind)
    rhos[44:46] *= non_compet_inhibit(S_IP, KPdiss)
    
    # inhibition factors
    # ******************
    unit_conversion = mass2mol_conversion(cmps)    
    if T_op == T_base:
        Ka = Kab
        KH = KHb / unit_conversion[[7,8,9,30]]
        Ksp = Kspb
    else:
        T_temp = params.pop('T_op', None)
        if T_op == T_temp:
            params['T_op'] = T_op
            Ka = params['Ka']
            KH = params['KH']
            Ksp = params['Ksp']
        else:
            params['T_op'] = T_op
            Ka = params['Ka'] = Kab * T_correction_factor(T_base, T_op, Ka_dH)
            KH = params['KH'] = KHb * T_correction_factor(T_base, T_op, KH_dH) / unit_conversion[[7,8,9,30]]
            Ksp = params['Ksp'] = Kspb * T_correction_factor(T_base, T_op, Ksp_dH)
            
    S_IN, S_IP = state_arr[[10,11]]
    I_nutrients = substr_inhibit(S_IN, KS_IN) * substr_inhibit(S_IP, KS_IP)
    rhos[3:11] *= I_nutrients
    rhos[[25,27,29,31,32]] *= I_nutrients
    
# =============================================================================
#     !!! place holder for PCM (speciation)
# =============================================================================
    pH, nh3, co2, acts = pcm(state_arr, params)
    Is_pH = Hill_inhibit(10**(-pH), pH_ULs, pH_LLs)
    rhos[3:9] *= Is_pH[0]
    rhos[9:11] *= Is_pH[1:3]
    rhos[[25,27]] *= Is_pH[3:5]
    rhos[[29,31,32]] *= Is_pH[-1]
    
    Is_h2 = non_compet_inhibit(S_h2, KIs_h2)
    rhos[5:9] *= Is_h2
    Inh3 = non_compet_inhibit(nh3, KI_nh3)
    rhos[9] *= Inh3
    
    Z_h2s = calc_biogas() # should be a function of pH, like co2 and nh3
    Is_h2s = non_compet_inhibit(Z_h2s, KIs_h2s)
    rhos[6:11] *= Is_h2s[:5]
    rhos[[25,27,29,31,32]] *= Is_h2s[5:]
    
    # multiple mineral precipitation
    # ******************************
    SIs = np.maximum(1.0, saturation_index(acts, Ksp))  # should be an array
    rhos[46:59] = k_cryst * state_arr[47:60] * (SIs**(1/sum_stoichios) - 1)**n_cryst
    
    # gas transfer
    # ************
    biogas_S = state_arr[[7,8,9,30]].copy()
    biogas_S[2] = co2 / unit_conversion[9]
    biogas_S[3] = Z_h2s / unit_conversion[30]
    biogas_p = R * T_op * state_arr[63:67]
    rhos[-4:] = kLa * (biogas_S - KH * biogas_p)
    
    return rhos

#%% modified ADM1 class
_load_components = settings.get_default_chemicals

def fun(q_aging_H=450.0, q_aging_L=0.1, q_Pcoprec=360, q_Pbinding=0.3, q_diss_H=36.0, q_diss_L=36.0,
        K_Pbind=37.2, K_Pdiss=0.93):
    '''
    

    Parameters
    ----------

    Returns
    -------
    None.

    '''
    pass    
    
@chemicals_user
class ModifiedADM1(CompiledProcesses):
    """
    Modified Anaerobic Digestion Model no.1 [1]_, [2]_, [3]_

    Parameters
    ----------
    f_ch_xb : float, optional
        Fraction of carbohydrates as biomass decay product. The default is 0.275.
    f_pr_xb : flaot, optional
        Fraction of proteins as biomass decay product. The default is 0.275.
    f_li_xb : float, optional
        Fraction of lipids as biomass decay product. The default is 0.35.
    f_xI_xb : float, optional
        Fraction of inert particulates as biomass decay product. The default is 0.1.
    f_va_pha : float, optional
        Fraction of valerate as PHA lysis product. The default is 0.1.
    f_bu_pha : float, optional
        Fraction of butyrate as PHA lysis product. The default is 0.1.
    f_pro_pha : float, optional
        Fraction of propionate as PHA lysis product. The default is 0.4.
    Y_PO4 : float, optional
        Poly-phosphorus (PP) required for PHA storage [kg P/kg COD]. The default is 0.4.
    Y_hSRB : float, optional
        Sulfide-reducing biomass (SRB) yield of hydrogen uptake [kg COD/kg COD]. 
        The default is 0.05.
    Y_aSRB : float, optional
        SRB yield of acetate uptake [kg COD/kg COD]. The default is 0.05.
    Y_pSRB : float, optional
        SRB yield of propionate uptake [kg COD/kg COD]. The default is 0.04.
    Y_c4SRB : float, optional
        SRB yield of butyrate or valerate uptake [kg COD/kg COD]. 
        The default is 0.06.
    q_pha : float, optional
        Maximum specific rate constant for PHA storage by phosphorus-accumulating
        organisms (PAOs) [d^(-1)]. The default is 3.0.
    b_pao : float, optional
        PAO lysis rate constant [d^(-1)]. The default is 0.2.
    b_pp : float, optional
        PP lysis rate constant [d^(-1)]. The default is 0.2.
    b_pha : float, optional
        PHA lysis rate constant [d^(-1)]. The default is 0.2.
    K_A : float, optional
        Substrate half saturation coefficient for PHA storage [kg COD/m3]. 
        The default is 4e-3.
    K_PP : float, optional
        PP half saturation coefficient for PHA storage [kg P (X_PP)/kg COD (X_PHA)]. 
        The default is 0.01.
    k_hSRB : float, optional
        Maximum specific growth rate constant of hydrogen-uptaking SRB [d^(-1)]. 
        The default is 41.125.
    k_aSRB : float, optional
        Maximum specific growth rate constant of acetate-uptaking SRB [d^(-1)]. 
        The default is 10..
    k_pSRB : float, optional
        Maximum specific growth rate constant of propionate-uptaking SRB [d^(-1)]. 
        The default is 16.25.
    k_c4SRB : float, optional
        Maximum specific growth rate constant of butyrate- or valerate-uptaking 
        SRB [d^(-1)]. The default is 23.
    b_hSRB : float, optional
        Hydrogen-uptaking SRB decay rate constant [d^(-1)]. The default is 0.02.
    b_aSRB : float, optional
        Acetate-uptaking SRB decay rate constant [d^(-1)]. The default is 0.02.
    b_pSRB : float, optional
        Propionate-uptaking SRB decay rate constant [d^(-1)]. The default is 0.02.
    b_c4SRB : float, optional
        Butyrate- or valerate-uptaking SRB decay rate constant [d^(-1)]. 
        The default is 0.02.
    K_hSRB : float, optional
        Substrate half saturation coefficient of hydrogen uptake by SRB 
        [kg COD/m3]. The default is 5.96e-6.
    K_aSRB : float, optional
        Substrate half saturation coefficient of acetate uptake by SRB 
        [kg COD/m3]. The default is 0.176.
    K_pSRB : float, optional
        Substrate half saturation coefficient of propionate uptake by SRB 
        [kg COD/m3]. The default is 0.088.
    K_c4SRB : float, optional
        Substrate half saturation coefficient of butyrate or valerate uptake by  
        SRB [kg COD/m3]. The default is 0.1739.
    K_so4_hSRB : float, optional
        Sulfate half saturation coefficient of SRB uptaking hydrogen [kg S/m3]. 
        The default is 3.335e-3.
    K_so4_aSRB : float, optional
        Sulfate half saturation coefficient of SRB uptaking acetate [kg S/m3]. 
        The default is 6.413e-3.
    K_so4_pSRB : float, optional
        Sulfate half saturation coefficient of SRB uptaking propionate [kg S/m3]. 
        The default is 6.413e-3.
    K_so4_c4SRB : float, optional
        Sulfate half saturation coefficient of SRB uptaking butyrate or valerate  
        [kg S/m3]. The default is 6.413e-3.
    k_Fe3t2_h2 : float, optional
        Fe(3+) reduction rate constant [m3∙kg^(-1) Fe(III)∙d^(-1)] using hydrogen
        as electron donor. The default is 1.79e7.
    k_Fe3t2_is : float, optional
        Fe(3+) reduction rate constant [m3∙kg^(-1) Fe(III)∙d^(-1)] using sulfide
        as electron donor. The default is 1.79e7.
    KS_IP : float, optional
        Inorganic phosphorus (nutrient) inhibition coefficient for soluble 
        substrate uptake [M]. The default is 2e-5.
    q_aging_H : float, optional
        Aging rate constant of X_HFO_H and X_HFO_HP [d^(-1)]. The default is 450.0.
    q_aging_L : float, optional
        Aging rate constant of X_HFO_L and X_HFO_LP [d^(-1)]. The default is 0.1.
    q_Pcoprec : float, optional
        Rate constant of P binding and coprecipitation on X_HFO_H [d^(-1)]. 
        The default is 360.
    q_Pbinding : float, optional
        Rate constant of P binding on X_HFO_L [d^(-1)]. The default is 0.3.
    q_diss_H : float, optional
        Dissolution rate constant of X_HFO_HP [d^(-1)]. The default is 36.0.
    q_diss_L : float, optional
        Dissolution rate constant of X_HFO_HP [d^(-1)]. The default is 36.0.
    K_Pbind : float, optional
        S_IP half saturation coefficient for binding with X_HFO_H or X_HFO_L
        [kg P/m3]. The default is 37.2, i.e., 1.20 kmol P/m3.
    K_Pdiss : float, optional
        S_IP half inhibition coefficient for dissolution of X_HFO_HP or X_HFO_LP
        [kg P/m3]. The default is 0.93, i.e., 0.03 kmol P/m3.
    KI_h2s_c4 : float, optional
        H2S half inhibition coefficient for butyrate or valerate uptake 
        [kg COD/m3]. The default is 0.481.
    KI_h2s_pro : float, optional
        H2S half inhibition coefficient for propionate uptake [kg COD/m3]. 
        The default is 0.481.
    KI_h2s_ac : float, optional
        H2S half inhibition coefficient for acetate uptake [kg COD/m3]. 
        The default is 0.460.
    KI_h2s_h2 : float, optional
        H2S half inhibition coefficient for hydrogen uptake [kg COD/m3]. 
        The default is 0.400.
    KI_h2s_c4SRB : float, optional
        H2S half inhibition coefficient for butyrate or valerate uptake by SRB
        [kg COD/m3]. The default is 0.520.
    KI_h2s_pSRB : float, optional
        H2S half inhibition coefficient for propionate uptake by SRB [kg COD/m3]. 
        The default is 0.520.
    KI_h2s_aSRB : float, optional
        H2S half inhibition coefficient for acetate uptake by SRB [kg COD/m3]. 
        The default is 0.499.
    KI_h2s_hSRB : float, optional
        H2S half inhibition coefficient for hydrogen uptake by SRB [kg COD/m3]. 
        The default is 0.499.        
    pH_limits_aa_SRB : 2-tuple, optional
        Lower and upper limits of pH inhibition for acetogenosis by SRB, 
        unitless. The default is (6,7).
    pH_limits_ac_SRB : 2-tuple, optional
        Lower and upper limits of pH inhibition for acetate uptake by SRB, 
        unitless. The default is (6,7).
    pH_limits_h2_SRB : 2-tuple, optional
        Lower and upper limits of pH inhibition for hydrogen uptake by SRB, 
        unitless. The default is (5,6).
    k_cryst : iterable[float], optional
        Mineral precipitation rate constants [h^(-1)], following the order of 
        `ModifiedADM1._precipitates`. The default is 
        [0.35, 1e-3, 3.0, 1e-3, 2.0, 0.76, 5.0, 1e-3, 1e-3, 1e-3, 1e2, 1e-3, 1e-3].
    n_cryst : iterable[int], optional
        The effect orders of mineral precipitation reactions [unitless], following 
        the order of `ModifiedADM1._precipitates`. The default is 
        [2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2].
    
    
    Examples
    --------
    ...

    References
    ----------
    .. [1] Flores-Alsina, X., Solon, K., Kazadi Mbamba, C., Tait, S., 
        Gernaey, K. V., Jeppsson, U., Batstone, D. J. (2016). 
        Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions 
        for dynamic simulations of anaerobic digestion processes. 
        Water Research, 95, 370–382. https://doi.org/10.1016/J.WATRES.2016.03.012
    .. [2] Solon, K., Flores-Alsina, X., Kazadi Mbamba, C., Ikumi, D., 
        Volcke, E. I. P., Vaneeckhaute, C., Ekama, G., Vanrolleghem, P. A., 
        Batstone, D. J., Gernaey, K. v., Jeppsson, U. (2017). Plant-wide 
        modelling of phosphorus transformations in wastewater treatment systems: 
        Impacts of control and operational strategies. Water Research, 
        113, 97–110. https://doi.org/10.1016/J.WATRES.2017.02.007
    .. [3] Hauduc, H., Takács, I., Smith, S., Szabo, A., Murthy, S., Daigger, G. T., 
        Spérandio, M. (2015). A dynamic physicochemical model for chemical phosphorus 
        removal. Water Research, 73, 157–170. https://doi.org/10.1016/J.WATRES.2014.12.053
    
    See Also
    --------
    `qsdsan.processes.ADM1 <https://qsdsan.readthedocs.io/en/latest/api/processes/ADM1.html>`_  

    """
        
    _cmp_dependent_stoichio = ('K_XPP', 'Mg_XPP', 
                               'MW_S0', 'MW_IS', 
                               'i_mass_S0', 'i_mass_IS', 'i_mass_Fe2')
    _stoichio_params = (*ADM1._stoichio_params[5:],
                        'f_ch_xb', 'f_pr_xb', 'f_li_xb', 'f_xI_xb', 'f_sI_xb',
                        'f_va_pha',	'f_bu_pha',	'f_pro_pha', 'f_ac_pha',
                        'f_is_pro', 'f_is_bu', 'f_is_va',
                        'Y_PO4', 'Y_hSRB', 'Y_aSRB', 'Y_pSRB', 'Y_c4SRB',
                        *_cmp_dependent_stoichio
                        )
    _kinetic_params = ('rate_constants', 'half_sat_coeffs', 'K_PP', 'K_so4', 
                       'pH_limits', 'KS_IN', 'KS_IP', 'KI_nh3', 'KIs_h2', 'KIs_h2s'
                       'Ka_base', 'Ka_dH', 'K_H_base', 'K_H_dH', 'kLa', 
                       'k_cryst', 'n_cryst', 'Ksp_base', 'Ksp_dH',
                       'T_base', 'components', 
                       # 'root'
                       )
    _acid_base_pairs = ADM1._acid_base_pairs
    _biogas_IDs = (*ADM1._biogas_IDs, 'S_IS')
    _biomass_IDs = (*ADM1._biomass_IDs, 'X_PAO', 'X_hSRB', 'X_aSRB', 'X_pSRB', 'X_c4SRB')
    _precipitates = ('X_CCM', 'X_ACC', 'X_ACP', 'X_HAP', 'X_DCPD', 'X_OCP',
                    'X_struv', 'X_newb', 'X_magn', 'X_kstruv', 
                    'X_FeS', 'X_Fe3PO42', 'X_AlPO4')
    _T_base = 298.15
    _K_H_base = [7.8e-4, 1.4e-3, 3.5e-2, 0.105]    # biogas species Henry's Law constant [M/bar]
    _K_H_dH = [-4180, -14240, -19410, -19180]      # Heat of reaction of liquid-gas transfer of biogas species [J/mol]
    
    _pKsp_base = [8.48, 8.3, 28.92, 44.333, 18.995, 47.08, 
                  13.6, 18.175, 7.46, 11.5508, 
                  2.95, 37.76, 18.2]
    _Ksp_dH = [8000, -12000, 54000, 0, 31000, 0, 
               -22600, -22600, -20000, -22600, 
               -11000, 5060, 0]
    
    def __new__(cls, components=None, path=None, 
                f_ch_xb=0.275, f_pr_xb=0.275, f_li_xb=0.35, f_xI_xb=0.1,
                f_fa_li=0.95, f_bu_su=0.1328, f_pro_su=0.2691, f_ac_su=0.4076,
                f_va_aa=0.23, f_bu_aa=0.26, f_pro_aa=0.05, f_ac_aa=0.4,
                f_ac_fa=0.7, f_pro_va=0.54, f_ac_va=0.31, f_ac_bu=0.8, f_ac_pro=0.57,
                Y_su=0.1, Y_aa=0.08, Y_fa=0.06, Y_c4=0.06, Y_pro=0.04, Y_ac=0.05, Y_h2=0.06,
                f_va_pha=0.1, f_bu_pha=0.1, f_pro_pha=0.4,
                Y_PO4=0.4, Y_hSRB=0.05, Y_aSRB=0.05, Y_pSRB=0.04, Y_c4SRB=0.06,                
                q_ch_hyd=10, q_pr_hyd=10, q_li_hyd=10,
                k_su=30, k_aa=50, k_fa=6, k_c4=20, k_pro=13, k_ac=8, k_h2=35,
                K_su=0.5, K_aa=0.3, K_fa=0.4, K_c4=0.2, K_pro=0.1, K_ac=0.15, K_h2=7e-6,
                b_su=0.02, b_aa=0.02, b_fa=0.02, b_c4=0.02, b_pro=0.02, b_ac=0.02, b_h2=0.02,
                q_pha=3.0, b_pao=0.2, b_pp=0.2, b_pha=0.2, K_A=4e-3, K_PP=0.01, 
                k_hSRB=41.125, k_aSRB=10., k_pSRB=16.25, k_c4SRB=23, 
                b_hSRB=0.02, b_aSRB=0.02, b_pSRB=0.02, b_c4SRB=0.02,
                K_hSRB=5.96e-6, K_aSRB=0.176, K_pSRB=0.088, K_c4SRB=0.1739,
                K_so4_hSRB=1.04e-4*S_mw, K_so4_aSRB=2e-4*S_mw, K_so4_pSRB=2e-4*S_mw, K_so4_c4SRB=2e-4*S_mw,
                k_Fe3t2_h2=1e9/Fe_mw, k_Fe3t2_is=1e9/Fe_mw,
                q_aging_H=450.0, q_aging_L=0.1, q_Pcoprec=360, q_Pbinding=0.3, q_diss_H=36.0, q_diss_L=36.0,
                K_Pbind=37.2, K_Pdiss=0.93, # 1.20 and 0.03 in MATLAB, assuming in kmol-P/m3 ?
                KI_h2_fa=5e-6, KI_h2_c4=1e-5, KI_h2_pro=3.5e-6, KI_nh3=1.8e-3, KS_IN=1e-4, KS_IP=2e-5,
                KI_h2s_c4=0.481, KI_h2s_pro=0.481, KI_h2s_ac=0.460, KI_h2s_h2=0.400,
                KI_h2s_c4SRB=0.520, KI_h2s_pSRB=0.520, KI_h2s_aSRB=0.499, KI_h2s_hSRB=0.499,
                pH_limits_aa=(4,5.5), pH_limits_ac=(6,7), pH_limits_h2=(5,6),
                pH_limits_aa_SRB=(6,7), pH_limits_ac_SRB=(6,7), pH_limits_h2_SRB=(5,6),
                kLa=200, pKa_base=[14, 9.25, 6.35, 4.76, 4.88, 4.82, 4.86],
                Ka_dH=[55900, 51965, 7646, 0, 0, 0, 0],
                k_cryst=[0.35, 1e-3, 3.0, 1e-3, 2.0, 0.76, 5.0, 1e-3, 1e-3, 1e-3, 1e2, 1e-3, 1e-3],
                n_cryst=[2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2],
                **kwargs):
        
        cmps = _load_components(components)

        if not path: path = _path
        self = Processes.load_from_file(path,
                                        components=cmps,
                                        conserved_for=('C', 'N', 'P'),
                                        parameters=cls._stoichio_params,
                                        compile=False)
        
        for i in ('fast_P_binding', 'slow_P_sorption', 'dissolution_HFO_HP', 'dissolution_HFO_LP'):
            p = getattr(self, i)
            p.ref_component = 'S_IP'
        
        precipitation = []
        for i in cls._precipitates[:-3]:
            new_p = Process('precipitation_%s' % i.lstrip('X_'),
                            reaction='[?]S_IC + [?]S_IN + [?]S_IP + [?]S_K + [?]S_Mg + [?]S_Ca -> %s' % i,
                            ref_component=i,
                            conserved_for=('C', 'N', 'P', 'K', 'Mg', 'Ca'),
                            parameters=())
            precipitation.append(new_p)
        
        i_mass_IS = cmps.S_IS.i_mass
        i_mass_Fe2 = cmps.S_Fe2.i_mass
        FeS_mw = cmps.X_FeS.chem_MW
        new_p = Process('precipitation_FeS',
                        reaction={'S_Fe2': -Fe_mw/FeS_mw/i_mass_Fe2,
                                  'S_IS': -S_mw/FeS_mw/i_mass_IS,
                                  'X_FeS': 1},
                        ref_component='X_FeS',
                        conserved_for=())
        precipitation.append(new_p)
        
        Fe3PO42_mw = cmps.X_Fe3PO42.chem_MW
        new_p = Process('precipitation_Fe3PO42',
                        reaction={'S_Fe2': -3*Fe_mw/Fe3PO42_mw/i_mass_Fe2,
                                  'S_IP': '?',
                                  'X_Fe3PO42': 1},
                        ref_component='X_Fe3PO42',
                        conserved_for=('P',))
        precipitation.append(new_p)
        
        AlPO4_mw = cmps.X_AlPO4.chem_MW
        Al_mw = cmps.S_Al.chem_MW
        new_p = Process('precipitation_AlPO4',
                        reaction={'S_Al': -Al_mw/AlPO4_mw,
                                  'S_IP': '?',
                                  'X_AlPO4': 1},
                        ref_component='X_AlPO4',
                        conserved_for=('P',))        
        precipitation.append(new_p)

        self.extend(precipitation)
        
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
                         f_ch_xb, f_pr_xb, f_li_xb, f_xI_xb, round(1.0-f_ch_xb-f_pr_xb-f_li_xb-f_xI_xb, 4),
                         f_va_pha, f_bu_pha, f_pro_pha, 1-f_va_pha-f_bu_pha-f_pro_pha,
                         1-f_ac_pro, 1-f_ac_bu, 1-f_pro_va-f_ac_va,
                         Y_PO4, Y_hSRB, Y_aSRB, Y_pSRB, Y_c4SRB,
                         cmps.X_PP.i_K, cmps.X_PP.i_Mg,
                         cmps.S_S0.chem_MW, cmps.S_IS.chem_MW, 
                         cmps.S_S0.i_mass, i_mass_IS, i_mass_Fe2)
        
        pH_limits = np.array([pH_limits_aa, pH_limits_ac, pH_limits_h2, 
                              pH_limits_h2_SRB, pH_limits_ac_SRB, pH_limits_aa_SRB]).T

        ks = np.array((q_ch_hyd, q_pr_hyd, q_li_hyd,
                       k_su, k_aa, k_fa, k_c4, k_c4, k_pro, k_ac, k_h2,
                       b_su, b_aa, b_fa, b_c4, b_pro, b_ac, b_h2,               # original ADM1
                       q_pha, q_pha, q_pha, q_pha, b_pao, b_pp, b_pha,          # P extension
                       k_hSRB, b_hSRB, k_aSRB, b_aSRB, k_pSRB, b_pSRB, k_c4SRB, k_c4SRB, b_c4SRB, # S extension
                       k_Fe3t2_h2, k_Fe3t2_h2, k_Fe3t2_is, k_Fe3t2_is,          # Fe extension
                       q_aging_H, q_aging_L, q_Pcoprec, q_Pbinding,             # HFO module
                       q_aging_H, q_aging_L, q_diss_H, q_diss_L))         
        
        Ks = np.array((K_su, K_aa, K_fa, K_c4, K_c4, K_pro, K_ac, K_h2,         # original ADM1
                       K_A,                                                     # P extension
                       K_hSRB, K_aSRB, K_pSRB, K_c4SRB,                         # S extension  
                       K_Pbind, K_Pdiss))                                       # HFO module                             
        K_so4 = np.array((K_so4_hSRB, K_so4_aSRB, K_so4_pSRB, K_so4_c4SRB))
        
        KIs_h2 = np.array((KI_h2_fa, KI_h2_c4, KI_h2_c4, KI_h2_pro))
        KIs_h2s = np.array((KI_h2s_c4, KI_h2s_c4, KI_h2s_pro, KI_h2s_ac, KI_h2s_h2,
                            KI_h2s_hSRB, KI_h2s_aSRB, KI_h2s_pSRB, KI_h2s_c4SRB, KI_h2s_c4SRB))
        K_H_base = np.array(cls._K_H_base)
        K_H_dH = np.array(cls._K_H_dH)
        Ka_base = np.array([10**(-pKa) for pKa in pKa_base])
        Ka_dH = np.array(Ka_dH)
        k_cryst = np.array(k_cryst) * 24    # converted to d^(-1)
        n_cryst = np.array(n_cryst)
        Ksp_base = np.array([10**(-pK) for pK in cls.pKsp_base])
        Ksp_dH = np.array(cls.Ksp_dH)
        # root = TempState()
        dct = self.__dict__
        dct.update(kwargs)

        dct['_parameters'] = dict(zip(cls._stoichio_params, stoichio_vals))
        self.set_rate_function(rhos_madm1)
        self.rate_function._params = dict(zip(cls._kinetic_params,
                                              [ks, Ks, K_PP, K_so4, 
                                                pH_limits, KS_IN*N_mw, KS_IP*P_mw,
                                                KI_nh3, KIs_h2, KIs_h2s, 
                                                Ka_base, Ka_dH, K_H_base, K_H_dH, kLa, 
                                                k_cryst, n_cryst, Ksp_base, Ksp_dH,
                                                cls.T_base, self._components, 
                                                # root,
                                                ]))
        return self