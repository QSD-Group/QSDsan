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
import qsdsan as qs
import numpy as np
# from scipy.optimize import fsolve
from math import exp
# from flexsolve import

# __all__ = ('load_adm1_cmps', 'ADM1')

# _path = os.path.join(data_path, 'process_data/_adm1.tsv')
# _path_cmps = os.path.join(data_path, '_adm1_cmps.pckl')
# _load_components = settings.get_default_chemicals

cmps_all = qs.Components.load_default(default_compile=False)

# varies
X_c = cmps_all.X_OHO.copy('X_c')
X_c.description = 'Composite'
X_c.i_N = 0.002 * 14


X_ch = qs.Component.from_chemical('X_ch', chemical='glycogen', Tc=1011.4, # glucan
                                  description='Carbohydrates',
                                  measured_as='COD',
                                  particle_size='Particulate',
                                  degradability='Slowly',
                                  organic=True)

# varies
X_pr = cmps_all.X_B_Subst.copy('X_pr')
X_pr.i_N = 0.007 * 14

X_li = qs.Component.from_chemical('X_li', chemical='tripalmitin',
                                  description='Lipids',
                                  measured_as='COD',
                                  particle_size='Particulate',
                                  degradability='Slowly',
                                  organic=True)

# both varies
X_I = cmps_all.X_U_Inf.copy('X_I')
S_I = cmps_all.S_U_Inf.copy('S_I')
X_I.i_N = S_I.i_N = 0.002 * 14

S_su = qs.Component.from_chemical('S_su', chemical='glucose',
                                  description='Monosaccharides',
                                  measured_as='COD',
                                  particle_size='Soluble',
                                  degradability='Readily',
                                  organic=True)

# varies
# S_aa = qs.Component.from_chemical('S_aa', chemical='valine',
#                                   description='Amino acids',
#                                   measured_as='COD',
#                                   particle_size='Soluble',
#                                   degradability='Readily',
#                                   organic=True)
S_aa = cmps_all.S_F.copy('S_aa')
S_aa.i_N = 0.007*14
S_aa.i_P = 0
S_aa.i_C = 0.313

S_fa = qs.Component.from_chemical('S_fa', chemical='palmitate',
                                  description='Total long-chain fatty acids',
                                  measured_as='COD',
                                  particle_size='Soluble',
                                  degradability='Readily',
                                  organic=True)

S_va = qs.Component.from_chemical('S_va', chemical='valerate',
                                  description='Total valerate',
                                  measured_as='COD',
                                  particle_size='Soluble',
                                  degradability='Readily',
                                  organic=True)

S_bu = qs.Component.from_chemical('S_bu', chemical='butyrate',
                                  description='Total butyrate',
                                  measured_as='COD',
                                  particle_size='Soluble',
                                  degradability='Readily',
                                  organic=True)

S_pro = cmps_all.S_Prop.copy('S_pro')
S_ac = cmps_all.S_Ac.copy('S_Ac')
S_h2 = cmps_all.S_H2.copy('S_h2')
S_ch4 = cmps_all.S_CH4.copy('S_ch4')

S_IC = qs.Component.from_chemical('S_IC', chemical='CO2',
                                  description='Inorganic carbon',
                                  particle_size='Dissolved gas',
                                  degradability='Undegradable',
                                  organic=False)

S_IN = qs.Component.from_chemical('S_IN', chemical='NH3',
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

for bio in (X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2): bio.formula = 'C5H7O2N'

S_cat = cmps_all.S_CAT.copy('S_cat')
S_an = cmps_all.S_AN.copy('S_an')

cmps_adm1 = qs.Components([S_su, S_aa, S_fa, S_va, S_bu, S_pro, S_ac, S_h2,
                           S_ch4, S_IC, S_IN, S_I, X_c, X_ch, X_pr, X_li,
                           X_su, X_aa, X_fa, X_c4, X_pro, X_ac, X_h2, X_I,
                           S_cat, S_an])
# cmps_adm1.default_compile()


pairs = ('H+', 'OH-', 'NH4+', 'NH3', 'CO2', 'HCO3-', 'HAc', 'Ac-',
         'HPr', 'Pr-', 'HBu', 'Bu-', 'HVa', 'Va-')

def acid_base_rxn(mols, state_arr, pKas):
    h, oh, nh4, nh3, co2, hco3, hac, ac, hpr, pr, hbu, bu, hva, va = mols
    S_cat, S_an, S_IN, S_IC, S_ac, S_pro, S_bu, S_va = state_arr  # in M
    Kw, Ka_nh, Ka_co2, Ka_ac, Ka_pr, Ka_bu, Ka_va = 10 ** pKas
    rhs = np.zeros_like(mols)
    rhs[0] = h + nh4 + S_cat - oh - hco3 - ac - pr - bu - va - S_an
    rhs[1] = S_IN - nh3 - nh4
    rhs[2] = S_IC - co2 - hco3
    rhs[3] = S_ac - ac - hac
    rhs[4] = S_pro - pr - hpr
    rhs[5] = S_bu - bu - hbu
    rhs[6] = S_va - va - hva
    rhs[7] = h*oh - Kw
    rhs[8] = h*nh3 - Ka_nh*nh4
    rhs[9] = h*hco3 - Ka_co2*co2
    rhs[10] = h*ac - Ka_ac*hac
    rhs[11] = h*pr - Ka_pr*hpr
    rhs[12] = h*bu - Ka_bu*hbu
    rhs[13] = h*va - Ka_va*hva
    return rhs

def jac(mols, state_arr, pKas):
    h, oh, nh4, nh3, co2, hco3, hac, ac, hpr, pr, hbu, bu, hva, va = mols
    S_cat, S_an, S_IN, S_IC, S_ac, S_pro, S_bu, S_va = state_arr  # in M
    Kw, Ka_nh, Ka_co2, Ka_ac, Ka_pr, Ka_bu, Ka_va = 10 ** pKas
    n = len(mols)
    fprime = np.zeros((n,n))
    fprime[0,[0,2]] = 1
    fprime[0,[1,5,7,9,11,13]] = -1
    fprime[1,2:4] = -1
    fprime[2,4:6] = -1
    fprime[3,6:8] = -1
    fprime[4,8:10] = -1
    fprime[5,10:12] = -1
    fprime[6,-2:] = -1
    fprime[7,[0,1]] = [oh, h]
    fprime[8,[2,3]] = [-Ka_nh, h]
    fprime[9,[4,5]] = [-Ka_co2, h]
    fprime[10,[6,7]] = [-Ka_ac, h]
    fprime[11,[8,9]] = [-Ka_pr, h]
    fprime[12,[10,11]] = [-Ka_bu, h]
    fprime[13,[12,13]] = [-Ka_va, h]
    return fprime

def low_pH_inhibit(pH, ul, ll):
    if pH >= ul: return 1
    else: return exp(-3 * ((pH-ul)/(ul-ll))**2)

def pH_inhibit(pH, ul, ll):
    return (1+2*10**(0.5*(ll-ul)))/(1+10**(pH-ul)+10**(ll-pH))

class ADM1(qs.Processes):

    _stoichio_params = ('f_ch_xc', 'f_pr_xc', 'f_li_xc', 'f_xI_xc',
                        'f_fa_li', 'f_bu_su', 'f_pro_su', 'f_ac_su',
                        'f_va_aa', 'f_bu_aa', 'f_pro_aa', 'f_ac_aa',
                        'f_ac_fa', 'f_pro_va', 'f_ac_va', 'f_ac_bu', 'f_ac_pro',
                        'Y_su', 'Y_aa', 'Y_fa', 'Y_c4', 'Y_pro', 'Y_ac', 'Y_h2')
    _kinetic_params = ('q_dis', 'q_ch_hyd', 'q_pr_hyd', 'q_li_hyd',
                       'k_su', 'k_aa', 'k_fa', 'k_c4', 'k_pro', 'k_ac', 'k_h2',
                       'K_su', 'K_aa', 'K_fa', 'K_c4', 'K_pro', 'K_ac', 'K_h2',
                       'b_su', 'b_aa', 'b_fa', 'b_c4', 'b_pro', 'b_ac', 'b_h2')

    def __new__(cls, components=None, path=None, N_xc=2e-3, N_I=2e-3, N_aa=7e-3,
                f_ch_xc=0.2, f_pr_xc=0.2, f_li_xc=0.25, f_xI_xc=0.25,
                f_fa_li=0.95, f_bu_su=0.13, f_pro_su=0.27, f_ac_su=0.41,
                f_va_aa=0.23, f_bu_aa=0.26, f_pro_aa=0.05, f_ac_aa=0.4,
                f_ac_fa=0.7, f_pro_va=0.54, f_ac_va=0.31, f_ac_bu=0.8, f_ac_pro=0.57,
                Y_su=0.1, Y_aa=0.08, Y_fa=0.06, Y_c4=0.06, Y_pro=0.04, Y_ac=0.05, Y_h2=0.06,
                q_dis=0.4, q_ch_hyd=0.25, q_pr_hyd=0.2, q_li_hyd=0.1,
                k_su=30, k_aa=50, k_fa=6, k_c4=20, k_pro=13, k_ac=8, k_h2=35,
                K_su=0.5, K_aa=0.3, K_fa=0.4, K_c4=0.3, K_pro=0.3, K_ac=0.15, K_h2=2.5e-5,
                b_su=0.02, b_aa=0.02, b_fa=0.02, b_c4=0.02, b_pro=0.02, b_ac=0.02, b_h2=0.02,
                KI_h2_fa=5e-6, KI_h2_c4=1e-5, KI_h2_pro=3.5e-6, KI_nh3=1.8e-3,
                **kwargs):
        pass