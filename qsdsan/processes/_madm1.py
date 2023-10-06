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

def create_madm1_cmps(set_thermo=True):
    
    # Components from the original ADM1
    # *********************************
    _cmps = pc.create_adm1_cmps(False)
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
    S_K = Component.from_chemical('S_K', chemical='K+',
                                  description='Potassium ion',
                                  measured_as='K',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
    S_Mg = Component.from_chemical('S_Mg', chemical='Mg2+',
                                  description='Magnesium ion',
                                  measured_as='Mg',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
    S_SO4 = Component.from_chemical('S_SO4', chemical='SO4-2',
                                  description='Sulfate',
                                  measured_as='S',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
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
                                  description='Sulfur',
                                  measured_as='COD',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
    S_Fe3 = Component.from_chemical('S_Fe3', chemical='Fe3+',
                                  description='Iron (III)',
                                  measured_as='Fe',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
    S_Fe2 = Component.from_chemical('S_Fe2', chemical='Fe2+',
                                  description='Iron (II)',
                                  measured_as='Fe',
                                  particle_size='Soluble',
                                  degradability='Undegradable',
                                  organic=False)
    S_Fe2.i_COD = 0.5*O_mw/Fe_mw
    S_Fe2.measured_as = 'COD'
    
    cmps_madm1 = Components([_cmps.S_su, S_aa, S_fa, _cmps.S_va, S_bu, 
                             S_pro, S_ac, _cmps.S_h2, _cmps.S_ch4, 
                             _cmps.S_IC, _cmps.S_IN, S_IP, S_I, 
                             X_ch, X_pr, X_li, *adm1_biomass, X_I,
                             X_PHA, asm_cmps.X_PP, X_PAO, S_K, S_Mg, 
                             S_SO4, S_IS, X_hSRB, X_aSRB, X_pSRB, X_c4SRB,
                             S_S0, S_Fe3, S_Fe2,
                             _cmps.H2O])
    cmps_madm1.default_compile()
    
    if set_thermo: qs.set_thermo(cmps_madm1)
    return cmps_madm1

#%% rate functions

def rhos_madm1():
    pass


#%% modified ADM1 class
_load_components = settings.get_default_chemicals

@chemicals_user
class ModifiedADM1(CompiledProcesses):
    """
    Modified Anaerobic Digestion Model no.1 [1]_

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
    k_Fe3t2 : float, optional
        Fe(3+) reduction rate constant [m3∙kg^(-1) Fe(III)∙d^(-1)]. 
        The default is 1.79e7.
    KS_IP : float, optional
        Inorganic phosphorus (nutrient) inhibition coefficient for soluble 
        substrate uptake [M]. The default is 2e-5.
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

    Examples
    --------
    ...

    References
    ----------
    .. [1] Flores-Alsina, X., Solon, K., Kazadi Mbamba, C., Tait, S., 
        Gernaey, K. V., Jeppsson, U., ; Batstone, D. J. (2016). 
        Modelling phosphorus (P), sulfur (S) and iron (Fe) interactions 
        for dynamic simulations of anaerobic digestion processes. 
        Water Research, 95, 370–382. https://doi.org/10.1016/J.WATRES.2016.03.012
            
    See Also
    --------
    `qsdsan.processes.ADM1 <https://qsdsan.readthedocs.io/en/latest/api/processes/ADM1.html>`_  

    """
        
    _cmp_dependent_stoichio = ('K_XPP', 'Mg_XPP', 
                               'MW_S0', 'MW_IS', 
                               'i_mass_S0', 'i_mass_IS', 'i_mass_Fe2')
    _stoichio_params = (*pc.ADM1._stoichio_params[5:],
                        'f_ch_xb', 'f_pr_xb', 'f_li_xb', 'f_xI_xb', 'f_sI_xb',
                        'f_va_pha',	'f_bu_pha',	'f_pro_pha', 'f_ac_pha',
                        'f_is_pro', 'f_is_bu', 'f_is_va',
                        'Y_PO4', 'Y_hSRB', 'Y_aSRB', 'Y_pSRB', 'Y_c4SRB',
                        *_cmp_dependent_stoichio
                        )
    _kinetic_params = ('rate_constants', 'half_sat_coeffs', 'K_PP', 'K_so4', 
                       'pH_limits', 'KS_IN', 'KI_nh3', 'KIs_h2',
                       'Ka_base', 'Ka_dH', 'K_H_base', 'K_H_dH', 'kLa',
                       'T_base', 'components', 'root',
                       )
    _acid_base_pairs = pc.ADM1._acid_base_pairs
    _biogas_IDs = (*pc.ADM1._biogas_IDs, 'S_IS')
    _biomass_IDs = (*pc.ADM1._biomass_IDs, 'X_PAO', 'X_hSRB', 'X_aSRB', 'X_pSRB', 'X_c4SRB')
    _T_base = 298.15
    _K_H_base = [7.8e-4, 1.4e-3, 3.5e-2, 0.105]    # biogas species Henry's Law constant [M/bar]
    _K_H_dH = [-4180, -14240, -19410, -19180]      # Heat of reaction of liquid-gas transfer of biogas species [J/mol]
    
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
                k_Fe3t2=1e9/Fe_mw,
                KI_h2_fa=5e-6, KI_h2_c4=1e-5, KI_h2_pro=3.5e-6, KI_nh3=1.8e-3, KS_IN=1e-4, KS_IP=2e-5,
                KI_h2s_c4=0.481, KI_h2s_pro=0.481, KI_h2s_ac=0.460, KI_h2s_h2=0.400,
                KI_h2s_c4SRB=0.520, KI_h2s_pSRB=0.520, KI_h2s_aSRB=0.499, KI_h2s_hSRB=0.499,
                pH_limits_aa=(4,5.5), pH_limits_ac=(6,7), pH_limits_h2=(5,6),
                pH_limits_aa_SRB=(6,7), pH_limits_ac_SRB=(6,7), pH_limits_h2_SRB=(5,6),
                kLa=200, pKa_base=[14, 9.25, 6.35, 4.76, 4.88, 4.82, 4.86],
                Ka_dH=[55900, 51965, 7646, 0, 0, 0, 0],**kwargs):
        
        cmps = _load_components(components)

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
                         f_ch_xb, f_pr_xb, f_li_xb, f_xI_xb, round(1.0-f_ch_xb-f_pr_xb-f_li_xb-f_xI_xb, 4),
                         f_va_pha, f_bu_pha, f_pro_pha, 1-f_va_pha-f_bu_pha-f_pro_pha,
                         1-f_ac_pro, 1-f_ac_bu, 1-f_pro_va-f_ac_va,
                         Y_PO4, Y_hSRB, Y_aSRB, Y_pSRB, Y_c4SRB,
                         cmps.X_PP.i_K, cmps.X_PP.i_Mg,
                         cmps.S_S0.chem_MW, cmps.S_IS.chem_MW, 
                         cmps.S_S0.i_mass, cmps.S_IS.i_mass, cmps.S_Fe2.i_mass)
        
        pH_limits = np.array([pH_limits_aa, pH_limits_ac, pH_limits_h2, 
                              pH_limits_h2_SRB, pH_limits_ac_SRB, pH_limits_aa_SRB]).T

        ks = np.array((q_ch_hyd, q_pr_hyd, q_li_hyd,
                       k_su, k_aa, k_fa, k_c4, k_c4, k_pro, k_ac, k_h2,
                       b_su, b_aa, b_fa, b_c4, b_pro, b_ac, b_h2,               # original ADM1
                       q_pha, q_pha, q_pha, q_pha, b_pao, b_pp, b_pha,          # P extension
                       k_hSRB, b_hSRB, k_aSRB, b_aSRB, k_pSRB, b_pSRB, k_c4SRB, k_c4SRB, b_c4SRB, # S extension
                       k_Fe3t2, k_Fe3t2))                                       # Fe extension
        
        Ks = np.array((K_su, K_aa, K_fa, K_c4, K_c4, K_pro, K_ac, K_h2,         # original ADM1
                       K_A,                                                     # P extension
                       K_hSRB, K_aSRB, K_pSRB, K_c4SRB))                        # S extension                       
        K_so4 = np.array((K_so4_hSRB, K_so4_aSRB, K_so4_pSRB, K_so4_c4SRB))
        
        KIs_h2 = np.array((KI_h2_fa, KI_h2_c4, KI_h2_c4, KI_h2_pro))
        KIs_h2s = np.array((KI_h2s_c4, KI_h2s_c4, KI_h2s_pro, KI_h2s_ac, KI_h2s_h2,
                            KI_h2s_hSRB, KI_h2s_aSRB, KI_h2s_pSRB, KI_h2s_c4SRB, KI_h2s_c4SRB))
        K_H_base = np.array(cls._K_H_base)
        K_H_dH = np.array(cls._K_H_dH)
        Ka_base = np.array([10**(-pKa) for pKa in pKa_base])
        Ka_dH = np.array(Ka_dH)
        # root = TempState()
        dct = self.__dict__
        dct.update(kwargs)

        self.set_rate_function(rhos_madm1)
        dct['_parameters'] = dict(zip(cls._stoichio_params, stoichio_vals))
        self.rate_function._params = dict(zip(cls._kinetic_params,
                                              [ks, Ks, K_PP, K_so4, 
                                               pH_limits, KS_IN*N_mw, KS_IP*P_mw,
                                               KI_nh3, KIs_h2, KIs_h2s, 
                                               Ka_base, Ka_dH, K_H_base, K_H_dH, kLa,
                                               cls.T_base, self._components, 
                                               # root,
                                               ]))
        return self