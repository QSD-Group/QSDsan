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
from chemicals.elements import molecular_weight as get_mw
from qsdsan import Component, Components, Process, Processes, CompiledProcesses
import numpy as np
from qsdsan.utils import ospath, data_path, save_pickle, load_pickled_cmps
from scipy.optimize import brenth

from warnings import warn

__all__ = ('create_adm1_cmps','load_adm1_cmps', 'ADM1')

_path = ospath.join(data_path, 'process_data/_adm1.tsv')
_path_cmps = ospath.join(data_path, '_adm1_cmps.pckl')
_load_components = settings.get_default_chemicals

#%%
# =============================================================================
# ADM1-specific components
# =============================================================================

C_mw = get_mw({'C':1})
N_mw = get_mw({'N':1})

def create_adm1_cmps(pickle=False):
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
    cmps_adm1.default_compile()
    
    if pickle:
        save_pickle(cmps_adm1, _path_cmps)
    return cmps_adm1

# cmps = create_adm1_cmps(False)
# create_adm1_cmps(True)

def load_adm1_cmps(pickle=None):
    return load_pickled_cmps(create_adm1_cmps, _path_cmps, pickle)


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

def calc_Kas(pKas, T_base, T_op, theta):
    pKas = np.asarray(pKas)
    return 10**(-pKas) * T_correction_factor(T_base, T_op, theta)

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
    
rhos = np.zeros(22) # 22 kinetic processes
Cs = np.empty(19)

def rhos_adm1(state_arr, params):
    ks = params['rate_constants']
    Ks = params['half_sat_coeffs']
    cmps = params['components']
    # n = len(cmps)
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
    # h_0 = params['root'].data
    
    # Cs_ids = cmps.indices(['X_c', 'X_ch', 'X_pr', 'X_li', 'X_su', 'X_aa', 
    #                        'X_fa', 'X_c4', 'X_c4', 'X_pro', 'X_ac', 'X_h2',
    #                        'X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro', 'X_ac', 'X_h2'])
    # Cs = state_arr[Cs_ids]
    Cs[:8] = state_arr[12:20]
    Cs[8:12] = state_arr[19:23]
    Cs[12:] = state_arr[16:23]
    # substrates_ids = cmps.indices(['S_su', 'S_aa', 'S_fa', 'S_va',
    #                                'S_bu', 'S_pro', 'S_ac', 'S_h2'])
    # substrates = state_arr[substrates_ids]
    substrates = state_arr[:8]
    # S_va, S_bu, S_h2, S_IN = state_arr[cmps.indices(['S_va', 'S_bu', 'S_h2', 'S_IN'])]
    # S_va, S_bu, S_h2, S_ch4, S_IC, S_IN = state_arr[[3,4,7,8,9,10]]
    S_va, S_bu, S_h2, S_IN = state_arr[[3,4,7,10]]
    unit_conversion = mass2mol_conversion(cmps)
    cmps_in_M = state_arr[:27] * unit_conversion
    weak_acids = cmps_in_M[[24, 25, 10, 9, 6, 5, 4, 3]]
    
    T_op = state_arr[-1]
    biogas_S = state_arr[7:10].copy()
    biogas_p = R * T_op * state_arr[27:30]
    Kas = Kab * T_correction_factor(T_base, T_op, Ka_dH)
    KH = KHb * T_correction_factor(T_base, T_op, KH_dH) / unit_conversion[7:10]
    
    rhos[:-3] = ks * Cs
    rhos[4:12] *= substr_inhibit(substrates, Ks)
    if S_va > 0: rhos[7] *= 1/(1+S_bu/S_va)
    if S_bu > 0: rhos[8] *= 1/(1+S_va/S_bu)
    
    # params['root'].data = h = newton(acid_base_rxn, h_0, fprime=fprime_abr, 
    #                                   args=(weak_acids, Kas), 
    #                                   tol=1e-12, maxiter=100)
    h = brenth(acid_base_rxn, 1e-14, 1.0, 
               args=(weak_acids, Kas), 
               xtol=1e-12, maxiter=100)
    # pH = 7.4655
    # params['root'].data = h = 10**(-pH)
    #!!! nh3 != IN, CO2 != IC
    nh3 = Kas[1] * weak_acids[2] / (Kas[1] + h)
    co2 = weak_acids[3] - Kas[2] * weak_acids[3] / (Kas[2] + h)
    biogas_S[-1] = co2 / unit_conversion[9]
    
    # rhos[4:12] *= pH_inhibit(pH, pH_ULs, pH_LLs) * substr_inhibit(S_IN, KS_IN)
    rhos[4:12] *= Hill_inhibit(h, pH_ULs, pH_LLs) * substr_inhibit(S_IN, KS_IN)
    rhos[6:10] *= non_compet_inhibit(S_h2, KIs_h2)
    rhos[10] *= non_compet_inhibit(nh3, KI_nh3)
    rhos[-3:] = kLa * (biogas_S - KH * biogas_p)
    
    return rhos

#%%        
# =============================================================================
# ADM1 class
# =============================================================================
class TempState:
    def __init__(self, length=0):
        if length > 0:
            self.data = np.zeros(length)
        else:
            self.data = 0.

@chemicals_user
class ADM1(CompiledProcesses):
    
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
                # Ka_theta=[0.076, 0.070, 0.010, 0, 0, 0, 0],
                Ka_dH=[55900, 51965, 7646, 0, 0, 0, 0],
                kLa=200, K_H_base=[7.8e-4, 1.4e-3, 3.5e-2], 
                # K_H_theta=[-5.66e-3, -1.929e-2, -2.629e-2],
                # !!! assume R unit is bar/M/K
                K_H_dH=[-4180, -14240, -19410],
                **kwargs):
        
        cmps = _load_components(components)
        cmps.X_c.i_N = N_xc * N_mw
        cmps.X_I.i_N = cmps.S_I.i_N = N_I * N_mw
        cmps.S_aa.i_N = cmps.X_pr.i_N = N_aa * N_mw
        
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
        
        stoichio_vals = (f_ch_xc, f_pr_xc, f_li_xc, f_xI_xc, 1-f_ch_xc-f_pr_xc-f_li_xc-f_xI_xc,
                         f_fa_li, f_bu_su, f_pro_su, f_ac_su, 1-f_bu_su-f_pro_su-f_ac_su,
                         f_va_aa, f_bu_aa, f_pro_aa, f_ac_aa, 1-f_va_aa-f_bu_aa-f_pro_aa-f_ac_aa,
                         f_ac_fa, 1-f_ac_fa, f_pro_va, f_ac_va, 1-f_pro_va-f_ac_va,
                         f_ac_bu, 1-f_ac_bu, f_ac_pro, 1-f_ac_pro,
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
        root.data = 10**(-7.4655)
        dct = self.__dict__
        dct.update(kwargs)
        
        self.set_rate_function(rhos_adm1)
        dct['_parameters'] = dict(zip(cls._stoichio_params, stoichio_vals))
        self.rate_function._params = dict(zip(cls._kinetic_params, 
                                              [ks, Ks, pH_ULs, pH_LLs, KS_IN*N_mw, 
                                               KI_nh3, KIs_h2, Ka_base, Ka_dH,
                                               K_H_base, K_H_dH, kLa, 
                                               T_base, self._components, root]))
        
        return self
    
    def set_pKas(self, pKas):
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
        i = self._find_index(process)
        self.rate_function._params['ks'][i] = k
    
    def set_half_sat_K(self, K, process):
        i = self._find_index(process)
        self.rate_function._params['Ks'][i-4] = K
    
    def set_pH_inhibit_bounds(self, process, lower=None, upper=None):
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
        '''Set values to stoichiometric parameters.'''
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
        stoichio = self.parameters   
        subst = ('xc', 'su', 'aa', 'fa', 'va', 'bu', 'pro')
        for s in subst:
            f_tot = sum([stoichio[k] for k in self._stoichio_params[:-7] \
                         if k.endswith(s)])
            if f_tot != 1: 
                raise ValueError(f"the sum of 'f_()_{s}' values must equal 1")      