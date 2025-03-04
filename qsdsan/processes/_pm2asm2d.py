# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''
# from thermosteam.utils import chemicals_user
from thermosteam import settings
from qsdsan import Component, Components, Process, Processes, CompiledProcesses
from qsdsan.utils import ospath, data_path
import numpy as np

__all__ = ('create_pm2asm2d_cmps', 'PM2ASM2d')

_path = ospath.join(data_path, 'process_data/_pm2asm2d_1.tsv')
_path_2 = ospath.join(data_path, 'process_data/_pm2asm2d_2.tsv')

# _load_components = settings.get_default_chemicals

#%%
# =============================================================================
# PM2ASM2d-specific components
# =============================================================================

def create_pm2asm2d_cmps(set_thermo=True):
    cmps = Components.load_default()

    # X_CHL (g Chl/m^3)
    X_CHL = Component(ID = 'X_CHL',
                      formula = 'C55H72MgN4O5',
                      description = 'Chlorophyll content of cells',
                      particle_size = 'Particulate',
                      degradability = 'Slowly',
                      organic = True)

    # X_ALG (g COD/m^3)
    X_ALG = cmps.X_OHO.copy('X_ALG')
    X_ALG.description = 'Concentration of carbon-accumulating mixotrophic organisms'
    X_ALG.formula = 'CH1.8O0.5N0.2P0.018'
    X_ALG.f_BOD5_COD = X_ALG.f_uBOD_COD = None
    X_ALG.f_Vmass_Totmass = 0.89

    # X_CH (g COD/m^3)
    X_CH = cmps.X_GAO_Gly.copy('X_CH')
    X_CH.description = 'Concentration of stored carbohydrates'
    X_CH.formula = 'CH2O'
    X_CH.f_BOD5_COD = X_CH.f_uBOD_COD = None

    # X_LI (g COD/m^3)
    X_LI = cmps.X_GAO_Gly.copy('X_LI')
    X_LI.description = 'Concentration of stored lipids'
    X_LI.formula = 'CH1.92O0.118'
    X_LI.f_BOD5_COD = X_LI.f_uBOD_COD = None

    # S_CO2 (g CO2/m^3)
    S_CO2 = Component.from_chemical(ID = 'S_CO2',
                                    chemical = 'CO2',
                                    description = 'Soluble carbon dioxide',
                                    particle_size = 'Soluble',
                                    degradability = 'Undegradable',
                                    organic = False)

    # S_A (g COD/m^3)
    S_A = cmps.S_Ac.copy('S_A')
    S_A.description = 'Concentration of extracellular dissolved organic carbon (acetate)'

    # S_F (g COD/m^3)
    S_F = Component.from_chemical(ID = 'S_F',
                                  chemical = 'glucose',
                                  description = 'Concentration of extracellular dissolved organic carbon (glucose)',
                                  measured_as = 'COD',
                                  particle_size = 'Soluble',
                                  degradability = 'Readily',
                                  organic = True)

    # S_O2 (g O2/m^3)
    S_O2 = cmps.S_O2.copy('S_O2')
    S_O2.description = ('Concentration of dissolved oxygen')

    # S_NH (g N/m^3)
    S_NH = cmps.S_NH4.copy('S_NH')
    S_NH.description = ('Concentration of dissolved ammonium')

    # S_NO (g N/m^3)
    S_NO = cmps.S_NO3.copy('S_NO')
    S_NO.description = ('Concentration of dissolved nitrate/nitrite')

    # S_P (g P/m^3)
    S_P = cmps.S_PO4.copy('S_P')
    S_P.description = ('Concentration of dissolved phosphorus')

    # X_N_ALG (g N/m^3)
    X_N_ALG = cmps.X_B_Subst.copy('X_N_ALG')
    X_N_ALG.description = 'Concentration of algal cell-associated nitrogen'
    X_N_ALG.measured_as = 'N'
    X_N_ALG.i_C = X_N_ALG.i_P = X_N_ALG.i_COD = X_N_ALG.f_BOD5_COD = X_N_ALG.f_uBOD_COD = X_N_ALG.f_Vmass_Totmass = 0
    X_N_ALG.i_mass = 1

    # X_P_ALG (g P/m^3)
    X_P_ALG = cmps.X_B_Subst.copy('X_P_ALG')
    X_P_ALG.description = 'Concentration of algal cell-associated phosphorus'
    X_P_ALG.measured_as = 'P'
    X_P_ALG.i_C = X_P_ALG.i_N = X_P_ALG.i_COD = X_P_ALG.f_BOD5_COD = X_P_ALG.f_uBOD_COD = X_P_ALG.f_Vmass_Totmass = 0
    X_P_ALG.i_mass = 1

    '''added from asm2d'''
    # S_N2 (g N/m^3)
    S_N2 = cmps.S_N2.copy('S_N2')
    S_N2.description = ('Concentration of dinitrogen')

    # S_ALK (g C/m^3)
    S_ALK = cmps.S_CO3.copy('S_ALK')                      # measured as g C, not as mole HCO3-
    S_ALK.description = ('Concentration of alkalinity')

    # S_I (g COD/m^3)
    S_I = cmps.S_U_E.copy('S_I')
    S_I.description = ('Concentration of inert soluble organic material')

    # X_I (g COD/m^3)
    X_I = cmps.X_U_OHO_E.copy('X_I')
    X_I.description = ('Concentration of inert particulate organic material')

    # X_S (g COD/m^3)
    X_S = cmps.X_B_Subst.copy('X_S')
    X_S.description = ('Concentration of slowly biodegradable substrates')

    # X_H (g COD/m^3)
    X_H = cmps.X_OHO.copy('X_H')
    X_H.description = ('Concentration of heterotrophic organisms (including denitrifer)')

    # X_AUT (g COD/m^3)
    X_AUT = cmps.X_AOO.copy('X_AUT')
    X_AUT.description = ('Concentration of nitrifying organisms')

    S_I.i_N = 0.01
    # S_F.i_N = 0.03
    X_I.i_N = 0.02
    X_S.i_N = 0.04
    X_H.i_N = X_AUT.i_N = 0.07

    S_I.i_P = 0.00
    # S_F.i_P = 0.01
    X_I.i_P = 0.01
    X_S.i_P = 0.01
    X_H.i_P = X_AUT.i_P = 0.02

    X_I.i_mass = 0.75
    X_S.i_mass = 0.75
    X_H.i_mass = X_AUT.i_mass = 0.9

    cmps_pm2asm2d = Components([X_CHL, X_ALG, X_CH, X_LI, S_CO2, S_A, S_F,
                           S_O2, S_NH, S_NO, S_P, X_N_ALG, X_P_ALG,
                           S_N2, S_ALK, S_I, X_I, X_S, X_H, X_AUT, cmps.H2O])

    cmps_pm2asm2d.default_compile()

    if set_thermo: settings.set_thermo(cmps_pm2asm2d)
    return cmps_pm2asm2d

# create_pm2asm2d_cmps()

#%%
# =============================================================================
# kinetic rate functions
# =============================================================================

# Calculation of ratio
def ratio(numerator, denominator, minimum, maximum):
    return min(max(minimum, numerator / denominator), maximum)

# Calculation of 'I_0' (for initial sensitivity analysis using calculated I_0)
def calc_irrad(t):
    '''
    :param t: time [days]
    :return: I_0, calculated irradiance [uE/m^2/s]

    -Assumes 14 hours of daylight
    '''
    daylight_hours = 14.0  # hours
    start_time = (12.0 - daylight_hours / 2) / 24.0  # days
    end_time = (12.0 + daylight_hours / 2) / 24.0  # days
    if t-np.floor(t) < start_time or t-np.floor(t) > end_time:
        return 0
    else:
        return 400.0 * (np.sin(2 * np.pi * (((t - np.floor(t)) - 5 / 24) / (14 / 24)) - np.pi / 2) + 1) / 2

# Calculation of 'I' from 'I_0' (Beer-Lambert)
def attenuation(light, X_TSS, a_c, b_reactor):
    '''
    :param light: I_0, calculated irradiance from 'calc_irrad' method (for sensitivity analysis) or
                       photosynthetically active radiation (PAR) imported from input excel file (for calibration & validation) [uE/m^2/s]
    :param X_TSS: total biomass concentration (X_ALG + X_CH + X_LI) * i_mass [g TSS/m^3]
    :param a_c: PAR absorption coefficient on a TSS (total suspended solids) basis [m^2/g TSS]
    :parma b_reactor: thickness of reactor along light path [m]
    :return: I, depth-averaged irradiance [uE/m^2/s]
    '''
    if X_TSS > 0:
        i_avg = (light * (1 - np.exp(-a_c * X_TSS * b_reactor))) / (a_c * X_TSS * b_reactor)
        return min(i_avg, light)
    else:
        return light

# Calculation of 'f_I' from 'I' (Eilers & Peeters)
def irrad_response(i_avg, X_CHL, X_carbon, I_n, I_opt):
    '''
    :param i_avg: I, depth-averaged irradiance (calculated from 'attenuation' method) [uE/m^2/s]
    :param X_CHL: chlorophyll content of cells [g Chl/m^3]
    :param X_carbon: carbon content of cells (X_ALG + X_CH + X_LI) * i_C [g C/m^3]
    :param I_n: maximum incident PAR irradiance (“irradiance at noon”) [uE/m^2/s]
    :param I_opt: optimal irradiance [uE/m^2/s]
    :return: f_I, irradiance response function [unitless]
    '''
    if X_carbon > 0:
        f_I = i_avg / (i_avg + I_n * (0.25 - (5 * X_CHL/X_carbon)) * ((i_avg ** 2 / I_opt ** 2) - (2 * i_avg / I_opt) + 1))
        return min(1, max(0, f_I))
    else:
        return 0

# Droop model
def droop(quota, subsistence_quota, exponent):
    '''
    :param quota: Q_N or Q_P [g N or g P/g COD]
    :param subsistence_quota: Q_N_min or Q_P_min [g N or g P/g COD]
    :param exponent: exponent to allow for more rapid transitions from growth to storage (see Guest et al., 2013) [unitless]
    :return: rate [unitless]
    '''
    return 1 - (subsistence_quota / quota) ** exponent

# Monod model
def monod(substrate, half_sat_const, exponent):
    '''
    :param substrate: S_NH, S_NO or S_P [g N or g P/m^3]
    :param half_sat_const: K_N or K_P [g N or g P/m^3]
    :param exponent: exponent to allow for more rapid transitions from growth to storage (see Guest et al., 2013) [unitless]
    :return: rate [unitless]
    '''
    return (substrate / (half_sat_const + substrate)) ** exponent

# Temperature model (Arrhenius)
def temperature(temp, arr_a, arr_e):
    '''
    :param temp: temperature (will be imported from input excel file) [K]
    :param arr_a: arrhenius constant (A) (Goldman et al., 1974) [unitless]
    :param arr_e: arrhenius exponential constant (E/R) (Goldman et al., 1974) [K]
    :return: temperature component of overall growth equation [unitless]
    '''
    return arr_a * np.exp(-arr_e / temp)  # Used equation from Goldman et al., 1974

# Photoadaptation (_p1)
def photoadaptation(i_avg, X_CHL, X_carbon, I_n, k_gamma):
    '''
    :param i_avg: I, depth-averaged irradiance (calculated from 'attenuation' method) [uE/m^2/s]
    :param X_CHL: chlorophyll content of cells [g Chl/m^3]
    :param X_carbon: carbon content of cells (X_ALG + X_CH + X_LI) * i_C [g C/m^3]
    :param I_n: maximum incident PAR irradiance (“irradiance at noon”) [uE/m^2/s]
    :param k_gamma: photoadaptation coefficient [unitless]
    :return: photoadaptation rate [g Chl/m^3/d]
    '''
    if X_carbon > 0:
        return 24 * ((0.2 * i_avg / I_n) / (k_gamma + (i_avg / I_n))) *\
               (0.01 + 0.03 * ((np.log(i_avg / I_n + 0.005)) / (np.log(0.01))) - X_CHL/X_carbon) * X_carbon
    else: return 0

# Nutrients uptake (_p2, _p3, _p4, _p5, _p6)
def nutrient_uptake(X_ALG, quota, substrate, uptake_rate, half_sat_const, maximum_quota, subsistence_quota):
    '''
    :param X_ALG: algae biomass concentration (i.e., no storage products) [g COD/m^3]
    :param quota: Q_N or Q_P [g N or g P/g COD]
    :param substrate: S_NH, S_NO or S_P [g N or g P/m^3]
    :param uptake_rate: V_NH, V_NO or V_P [g N or g P/g COD/d]
    :param half_sat_const: K_N or K_P [g N or g P/m^3]
    :param maximum_quota: Q_N_max or Q_P_max [g N or g P/g COD]
    :param subsistence_quota: Q_N_min or Q_P_min [g N or g P/g COD]
    :return: nutrient uptake rate [g N or g P/m^3/d]
    '''
    return uptake_rate * monod(substrate, half_sat_const, 1) * X_ALG * \
           ((maximum_quota - quota) / (maximum_quota - subsistence_quota)) ** 0.01

# Maximum total photoautotrophic or heterotrophic-acetate or heterotrophic-glucose growth rate (_p7, _p10, _p11, _p15, _p18, _p19, _p23, _p26, _p27)
def max_total_growth(X_ALG, mu_max, f_np, f_temp):
    '''
    :param X_ALG: algae biomass concentration (i.e., no storage products) [g COD/m^3]
    :param mu_max: maximum specific growth rate [d^(-1)]
    :param f_np: inhibition factor by nitrogen or phosphorus (between 0 and 1) [unitless]
    :param f_temp: temperature correction factor (between 0 and 1) [unitless]
    :return: maximum total growth rate for a particular mechanism,
             without considering carbon source or light inhibition (= product of shared terms in growth-related equations) [g COD/m^3/d]
    '''
    return mu_max * f_np * X_ALG * f_temp

# Split the total growth rate between three processes (_p7, _p10, _p11, _p15, _p18, _p19, _p23, _p26, _p27)
def growth_split(f_I, f_CH, f_LI, rho, Y_CH, Y_LI, K_STO):
    '''
    :param f_I: irradiance response function (calculated from 'irrad_response' method) [unitless]
    :param f_CH: ratio of stored carbohydrates to cells (X_CH / X_ALG) [g COD/g COD]
    :param f_LI: ratio of stored lipids to cells (X_LI / X_ALG) [g COD/g COD]
    :param rho: carbohydrate relative preference factor (calibrated in Guest et al., 2013) [unitless]
    :param Y_CH: yield of storage carbohydrates (as polyglucose, PG), Y_CH_PHO, Y_CH_NR_HET_ACE, or Y_CH_NR_HET_GLU [g COD/g COD]
    :param Y_LI: yield of storage lipids (as triacylglycerol, TAG), Y_LI_PHO, Y_LI_NR_HET_ACE, or Y_LI_NR_HET_GLU [g COD/g COD]
    :param K_STO: half-saturation constant for stored organic carbon (calibrated in Guest et al., 2013) [g COD/g COD]
    :return: splits the total growth rate between three processes,
             growth, growth on stored carbohydrates, and growth on stored lipid (= process-specific terms) [unitless]
    '''
    numerators = np.asarray([K_STO * (1 - f_I), rho * f_CH, f_LI * Y_CH / Y_LI])
    return numerators/(sum(numerators))

# Part of storage equations (_p8, _p9, _p16, _p17, _p24, _p25)
def storage_saturation(f, f_max, beta):
    '''
    :param f: f_CH or f_LI [g COD/g COD]
    :param f_max: f_CH_max or f_LI_max [g COD/g COD]
    :param beta: beta_1 or beta_2 [unitless]
    :return: part of storage equations [unitless]
    '''
    return 1 - (f / f_max) ** beta

# Maximum total photoautotrophic or heterotrophic-acetate or heterotrophic-glucose maintenance rate (_p12, _p13, _p14, _p20, _p21, _p22, _p28, _p29, _p30)
def max_total_maintenance(X_ALG, m_ATP):
    '''
    :param X_ALG: algae biomass concentration (i.e., no storage products) [g COD/m^3]
    :param m_ATP: specific maintenance rate [g ATP/g COD/d]
    :return: maximum total maintenance rate for a particular mechanism
             (= product of shared terms in maintenance-related equations) [g COD/m^3/d]
    '''
    return m_ATP * X_ALG

# Split the total maintenance rate between three processes (_p12, _p13, _p14, _p20, _p21, _p22, _p28, _p29, _p30)
def maintenance_split(f_CH, f_LI, rho, Y_CH, Y_LI, Y_X_ALG, Y_ATP, K_STO):
    '''
    :param f_CH: ratio of stored carbohydrates to cells (X_CH / X_ALG) [g COD/g COD]
    :param f_LI: ratio of stored lipids to cells (X_LI / X_ALG) [g COD/g COD]
    :param rho: carbohydrate relative preference factor (calibrated in Guest et al., 2013) [unitless]
    :param Y_CH: yield of storage carbohydrates (as polyglucose, PG), Y_CH_PHO, Y_CH_NR_HET_ACE, or Y_CH_NR_HET_GLU [g COD/g COD]
    :param Y_LI: yield of storage lipids (as triacylglycerol, TAG), Y_LI_PHO, Y_LI_NR_HET_ACE, or Y_LI_NR_HET_GLU [g COD/g COD]
    :param Y_X_ALG: yield of carbon-accumulating phototrophic organisms, Y_X_ALG_PHO, Y_X_ALG_HET_ACE, or Y_X_ALG_HET_GLU [g COD/g COD]
    :param Y_ATP: yield of ATP, Y_ATP_PHO, Y_ATP_HET_ACE, or Y_ATP_HET_GLU [g ATP/g COD]
    :param K_STO: half-saturation constant for stored organic carbon (calibrated in Guest et al., 2013) [g COD/g COD]
    :return: splits the total maintenance rate between three processes,
             stored carbohydrate degradation, stored lipid degradation, and endogenous respiration (= process-specific terms) [unitless]
    '''
    numerators = np.asarray([rho * f_CH, f_LI * Y_CH / Y_LI, K_STO])
    yield_ratios = np.asarray([Y_CH, Y_LI, Y_X_ALG]) / Y_ATP
    return numerators/(sum(numerators)) * yield_ratios

# Storage of carbohydrate/lipid (_p8, _p9, _p16, _p17, _p24, _p25)
def storage(X_ALG, f_np, response, saturation, storage_rate):
    '''
    :param X_ALG: algae biomass concentration (i.e., no storage products) [g COD/m^3]
    :param f_np: inhibition factor by nitrogen or phosphorus (between 0 and 1) [unitless]
    :param response: f_I (irradiance response function, calculated from 'irrad_response' method), acetate_response (monod(S_A, K_A, 1)), or glucose_response (monod(S_F, K_F, 1)) [unitless]
    :param saturation: 1 - (f / f_max) ** beta (calculated from 'storage_saturation' method) [unitless]
    :param storage_rate: q_CH or q_LI [g COD/g COD/d]
    :return: storage rate [g COD/m^3/d]
    '''
    return storage_rate * saturation * (1 - f_np) * response * X_ALG

'''added from asm2d'''

# Hydrolysis (_p31, _p32, _p33)
def hydrolysis(X_S, X_H, K_h, K_X):
    '''
    :param X_S: concentration of slowly biodegradable substrates
    :param X_H: concentration of heterotrophic organisms (including denitrifer)
    :param K_h: hydrolysis rate constant
    :param K_X: slowly biodegradable substrate half saturation coefficient for hydrolysis
    :return: shared parts of hydrolysis equations
    '''
    return K_h * (X_S/X_H) / (K_X + X_S/X_H) * X_H

# Growth in ASM2d (_p34, _p35, _p36, _p37, _p40)
def growth_asm2d(S_NH, S_P, S_ALK, mu, X, K_NH4, K_P, K_ALK):
    '''
    :param S_NH: concentration of dissolved ammonium
    :param S_P: concentration of dissolved phosphorus
    :param S_ALK: concentration of alkalinity
    :param mu: maximum specific growth rate (mu_H or mu_AUT)
    :param X: concentration of biomass (X_H or X_AUT)
    :param K_NH4: ammonium (nutrient) half saturation coefficient (K_NH4_H or K_NH4_AUT)
    :param K_P: phosphorus (nutrient) half saturation coefficient (K_P_H or K_P_AUT)
    :param K_ALK: alkalinity half saturation coefficient  (K_ALK_H or K_ALK_AUT)
    :return: shared parts of growth-related equations
    '''
    return mu * S_NH/(K_NH4+S_NH) * S_P/(K_P+S_P) * S_ALK/(K_ALK + S_ALK) * X

def rhos_pm2asm2d(state_arr, params):

    # extract values of state variables
    c_arr = state_arr[:21]
    temp = state_arr[22]
    light = state_arr[23]     # imported from input file assumed

    # Q = state_arr[21]                             # Flow rate
    # t = state_arr[22]                             # time

    X_CHL, X_ALG, X_CH, X_LI, S_CO2, S_A, S_F, S_O2, S_NH, S_NO, S_P, X_N_ALG, X_P_ALG, S_N2, S_ALK, S_I, X_I, X_S, X_H, X_AUT, H2O = c_arr

    # extract values of parameters
    cmps = params['cmps']
    a_c = params['a_c']
    I_n = params['I_n']
    arr_a = params['arr_a']
    arr_e = params['arr_e']
    beta_1 = params['beta_1']
    beta_2 = params['beta_2']
    b_reactor = params['b_reactor']
    I_opt = params['I_opt']
    k_gamma = params['k_gamma']
    K_N = params['K_N']
    K_P = params['K_P']
    K_A = params['K_A']
    K_F = params['K_F']
    rho = params['rho']
    K_STO = params['K_STO']
    f_CH_max = params['f_CH_max']
    f_LI_max = params['f_LI_max']
    m_ATP = params['m_ATP']
    mu_max = params['mu_max']
    q_CH = params['q_CH']
    q_LI = params['q_LI']
    Q_N_max = params['Q_N_max']
    Q_N_min = params['Q_N_min']
    Q_P_max = params['Q_P_max']
    Q_P_min = params['Q_P_min']
    V_NH = params['V_NH']
    V_NO = params['V_NO']
    V_P = params['V_P']
    exponent = params['exponent']
    Y_ATP_PHO = params['Y_ATP_PHO']
    Y_CH_PHO = params['Y_CH_PHO']
    Y_LI_PHO = params['Y_LI_PHO']
    Y_X_ALG_PHO = params['Y_X_ALG_PHO']
    Y_ATP_HET_ACE = params['Y_ATP_HET_ACE']
    Y_CH_NR_HET_ACE = params['Y_CH_NR_HET_ACE']
    Y_LI_NR_HET_ACE = params['Y_LI_NR_HET_ACE']
    Y_X_ALG_HET_ACE = params['Y_X_ALG_HET_ACE']
    Y_ATP_HET_GLU = params['Y_ATP_HET_GLU']
    Y_CH_NR_HET_GLU = params['Y_CH_NR_HET_GLU']
    Y_LI_NR_HET_GLU = params['Y_LI_NR_HET_GLU']
    Y_X_ALG_HET_GLU = params['Y_X_ALG_HET_GLU']
    n_dark = params['n_dark']

    '''added from asm2d'''
    # f_SI = params['f_SI']
    # Y_H = params['Y_H']
    # f_XI_H = params['f_XI_H']
    # Y_A = params['Y_A']
    # f_XI_AUT = params['f_XI_AUT']
    K_h = params['K_h']
    eta_NO3 = params['eta_NO3']
    eta_fe = params['eta_fe']
    K_O2 = params['K_O2']
    K_NO3 = params['K_NO3']
    K_X = params['K_X']
    mu_H = params['mu_H']
    q_fe = params['q_fe']
    eta_NO3_H = params['eta_NO3_H']
    b_H = params['b_H']
    K_O2_H = params['K_O2_H']
    K_F_H = params['K_F_H']                 # K_F overlaps with PM2 -> change into K_F_H
    K_fe = params['K_fe']
    K_A_H = params['K_A_H']
    K_NO3_H = params['K_NO3_H']
    K_NH4_H = params['K_NH4_H']
    K_P_H = params['K_P_H']
    K_ALK_H = params['K_ALK_H']
    mu_AUT = params['mu_AUT']
    b_AUT = params['b_AUT']
    K_O2_AUT = params['K_O2_AUT']
    K_NH4_AUT = params['K_NH4_AUT']
    K_ALK_AUT = params['K_ALK_AUT']
    K_P_AUT = params['K_P_AUT']

# intermediate variables
    f_CH = ratio(X_CH, X_ALG, 0, f_CH_max)
    f_LI = ratio(X_LI, X_ALG, 0, f_LI_max)

    # Q_N = ratio(X_N_ALG, X_ALG, Q_N_min, Q_N_max)
    # Q_P = ratio(X_P_ALG, X_ALG, Q_P_min, Q_P_max)

    alg_iN, alg_iP = cmps.X_ALG.i_N, cmps.X_ALG.i_P
    Q_N = ratio(X_N_ALG+X_ALG*alg_iN, X_ALG, Q_N_min, Q_N_max)
    Q_P = ratio(X_P_ALG+X_ALG*alg_iP, X_ALG, Q_P_min, Q_P_max)

    idx = cmps.indices(['X_ALG', 'X_CH', 'X_LI'])
    X_bio = np.array([X_ALG, X_CH, X_LI])
    X_TSS = sum(X_bio * cmps.i_mass[idx])
    X_carbon = sum(X_bio * cmps.i_C[idx])

    i_avg = attenuation(light, X_TSS, a_c, b_reactor)
    f_I = irrad_response(i_avg, X_CHL, X_carbon, I_n, I_opt)
    dark_response = max(f_I, n_dark)
    acetate_response = monod(S_A, K_A, 1)
    glucose_response = monod(S_F, K_F, 1)

    f_np = min(droop(Q_N, Q_N_min, exponent), droop(Q_P, Q_P_min, exponent))
    f_temp = temperature(temp, arr_a, arr_e)

    f_sat_CH = storage_saturation(f_CH, f_CH_max, beta_1)
    f_sat_LI = storage_saturation(f_LI, f_LI_max, beta_2)

    max_total_growth_rho = max_total_growth(X_ALG, mu_max, f_np, f_temp)
    max_maintenance_rho = max_total_maintenance(X_ALG, m_ATP)
    # light = calc_irrad(t)

    # calculate kinetic rate values
    rhos = np.empty(41)

    rhos[0] = photoadaptation(i_avg, X_CHL, X_carbon, I_n, k_gamma)

    rhos[1] = nutrient_uptake(X_ALG, Q_N, S_NH, V_NH, K_N, Q_N_max, Q_N_min)
    rhos[[2,3,4]] = nutrient_uptake(X_ALG, Q_N, S_NO, V_NO, K_N, Q_N_max, Q_N_min) * (K_N/(K_N + S_NH))

    rhos[5] = nutrient_uptake(X_ALG, Q_P, S_P, V_P, K_P, Q_P_max, Q_P_min)

    rhos[[6,9,10]] = max_total_growth_rho \
        * growth_split(f_I, f_CH, f_LI, rho, Y_CH_PHO, Y_LI_PHO, K_STO)
    rhos[6] *= f_I
    rhos[[9,10]] *= dark_response

    rhos[[14,17,18]] = max_total_growth_rho \
        * acetate_response \
        * growth_split(f_I, f_CH, f_LI, rho, Y_CH_NR_HET_ACE, Y_LI_NR_HET_ACE, K_STO)

    rhos[[22,25,26]] = max_total_growth_rho \
        * glucose_response \
        * growth_split(f_I, f_CH, f_LI, rho, Y_CH_NR_HET_GLU, Y_LI_NR_HET_GLU, K_STO)

    rhos[[11,12,13]] = max_maintenance_rho \
        * maintenance_split(f_CH, f_LI, rho, Y_CH_PHO, Y_LI_PHO, Y_X_ALG_PHO, Y_ATP_PHO, K_STO)

    rhos[[19,20,21]] = max_maintenance_rho \
        * maintenance_split(f_CH, f_LI, rho, Y_CH_NR_HET_ACE, Y_LI_NR_HET_ACE, Y_X_ALG_HET_ACE, Y_ATP_HET_ACE, K_STO)

    rhos[[27,28,29]] = max_maintenance_rho \
        * maintenance_split(f_CH, f_LI, rho, Y_CH_NR_HET_GLU, Y_LI_NR_HET_GLU, Y_X_ALG_HET_GLU, Y_ATP_HET_GLU, K_STO)

    rhos[7] = storage(X_ALG, f_np, f_I, f_sat_CH, q_CH)
    rhos[8] = storage(X_ALG, f_np, f_I, f_sat_LI, q_LI) * (f_CH / f_CH_max)
    rhos[15] = storage(X_ALG, f_np, acetate_response, f_sat_CH, q_CH)
    rhos[16] = storage(X_ALG, f_np, acetate_response, f_sat_LI, q_LI) * (f_CH / f_CH_max)
    rhos[23] = storage(X_ALG, f_np, glucose_response, f_sat_CH, q_CH)
    rhos[24] = storage(X_ALG, f_np, glucose_response, f_sat_LI, q_LI) * (f_CH / f_CH_max)

    rhos[30] = hydrolysis(X_S, X_H, K_h, K_X) * monod(S_O2, K_O2, 1)
    rhos[31] = hydrolysis(X_S, X_H, K_h, K_X) * eta_NO3 * monod(K_O2, S_O2, 1) * monod(S_NO, K_NO3, 1)
    rhos[32] = hydrolysis(X_S, X_H, K_h, K_X) * eta_fe * monod(K_O2, S_O2, 1) * monod(K_NO3, S_NO, 1)

    rhos[33] = growth_asm2d(S_NH, S_P, S_ALK, mu_H, X_H, K_NH4_H, K_P_H, K_ALK_H) * monod(S_O2, K_O2_H, 1) * monod(S_F, K_F_H, 1) * monod(S_F, S_A, 1)
    rhos[34] = growth_asm2d(S_NH, S_P, S_ALK, mu_H, X_H, K_NH4_H, K_P_H, K_ALK_H) * monod(S_O2, K_O2_H, 1) * monod(S_A, K_A_H, 1) * monod(S_A, S_F, 1)
    rhos[35] = growth_asm2d(S_NH, S_P, S_ALK, mu_H, X_H, K_NH4_H, K_P_H, K_ALK_H) * eta_NO3_H * monod(K_O2_H, S_O2, 1) * monod(S_NO, K_NO3_H, 1) * monod(S_F, K_F_H, 1) * monod(S_F, S_A, 1)
    rhos[36] = growth_asm2d(S_NH, S_P, S_ALK, mu_H, X_H, K_NH4_H, K_P_H, K_ALK_H) * eta_NO3_H * monod(K_O2_H, S_O2, 1) * monod(S_NO, K_NO3_H, 1) * monod(S_A, K_A_H, 1) * monod(S_A, S_F, 1)
    rhos[37] = q_fe * monod(K_O2_H, S_O2, 1) * monod(K_NO3_H, S_NO, 1) * monod(S_F, K_fe, 1) * monod(S_ALK, K_ALK_H, 1) * X_H
    rhos[38] = b_H * X_H

    rhos[39] = growth_asm2d(S_NH, S_P, S_ALK, mu_AUT, X_AUT, K_NH4_AUT, K_P_AUT, K_ALK_AUT) * monod(S_O2, K_O2_AUT, 1)
    rhos[40] = b_AUT * X_AUT

    return rhos

#%%
# =============================================================================
# PM2ASM2d class
# =============================================================================

class PM2ASM2d(CompiledProcesses):
    '''
    Parameters
    ----------
    components: class:`CompiledComponents`, optional
              Components corresponding to each entry in the stoichiometry array,
              defaults to thermosteam.settings.chemicals.
    a_c : float, optional
              PAR absorption coefficient on a TSS (total suspended solids) basis, in [m^2/g TSS].
              The default is 0.049.
    I_n : float, optional
              Maximum incident PAR irradiance (“irradiance at noon”), in [uE/m^2/s].
              The default is 250.
    arr_a : float, optional
              Arrhenius constant (A), in [unitless].
              The default is 1.8 * 10**10.
    arr_e : float, optional
              Arrhenius exponential constant (E/R), in [K].
              The default is 6842.
    beta_1 : float, optional
              Power coefficient for carbohydrate storage inhibition, in [unitless].
              The default is 2.90.
    beta_2 : float, optional
              Power coefficient for lipid storage inhibition, in [unitless].
              The default is 3.50.
    b_reactor : float, optional
              Thickness of reactor along light path, in [m].
              The default is 0.03.
    I_opt : float, optional
              Optimal irradiance, in [uE/m^2/s].
              The default is 300.
    k_gamma : float, optional
              Photoadaptation coefficient, in [unitless].
              The default is 0.00001.
    K_N : float, optional
              Nitrogen half-saturation constant, in [g N/m^3].
              The default is 0.1.
    K_P : float, optional
              Phosphorus half-saturation constant, in [g P/m^3].
              The default is 1.0.
    K_A : float, optional
              Organic carbon half-saturation constant (acetate) (Wagner, 2016), in [g COD/m^3].
              The default is 6.3.
    K_F : float, optional
              Organic carbon half-saturation constant (glucose); assumes K_A = K_F, in [g COD/m^3].
              The default is 6.3.
    rho : float, optional
              Carbohydrate relative preference factor (calibrated in Guest et al., 2013), in [unitless].
              The default is 1.186.
    K_STO : float, optional
              Half-saturation constant for stored organic carbon (calibrated in Guest et al., 2013), in [g COD/g COD].
              The default is 1.566.
    f_CH_max : float, optional
              Maximum achievable ratio of stored carbohydrates to functional cells, in [g COD/g COD].
              The default is 0.819.
    f_LI_max : float, optional
              Maximum achievable ratio of stored lipids to functional cells, in [g COD/g COD].
              The default is 3.249.
    m_ATP : float, optional
              Specific maintenance rate, in [g ATP/g COD/d].
              The default is 15.835.
    mu_max : float, optional
              Maximum specific growth rate, in [d^(-1)].
              The default is 1.969.
    q_CH : float, optional
              Maximum specific carbohydrate storage rate, in [g COD/g COD/d].
              The default is 0.594.
    q_LI : float, optional
              Maximum specific lipid storage rate, in [g COD/g COD/d].
              The default is 0.910.
    Q_N_max : float, optional
              Maximum nitrogen quota, in [g N/g COD].
              The default is 0.417.
    Q_N_min : float, optional
              Nitrogen subsistence quota, in [g N/g COD].
              The default is 0.082.
    Q_P_max : float, optional
              Maximum phosphorus quota, in [g P/g COD].
              The default is 0.092.
    Q_P_min : float, optional
              Phosphorus subsistence quota; assumes N:P ratio of 5:1, in [g P/g COD].
              The default is 0.0163.
    V_NH : float, optional
              Maximum specific ammonium uptake rate (calibrated in Guest et al., 2013), in [g N/g COD/d].
              The default is 0.254.
    V_NO : float, optional
              Maximum specific nitrate uptake rate (calibrated in Guest et al., 2013), in [g N/g COD/d].
              The default is 0.254.
    V_P : float, optional
              Maximum specific phosphorus uptake rate (calibrated in Guest et al., 2013), in [g P/g COD/d].
              The default is 0.016.
    exponent : float, optional
              Exponent to allow for more rapid transitions from growth to storage (see Guest et al., 2013), in [unitless]
              The default is 4.
    Y_ATP_PHO : float, optional
              Yield of ATP on CO2 fixed to G3P, in [g ATP/g CO2].
              The default is 55.073.
    Y_CH_PHO : float, optional
              Yield of storage carbohydrate (as polyglucose, PG) on CO2 fixed to G3P, in [g COD/g CO2].
              The default is 0.754.
    Y_LI_PHO : float, optional
              Yield of storage lipids (as triacylglycerol, TAG) on CO2 fixed to G3P, in [g COD/g CO2].
              The default is 0.901.
    Y_X_ALG_PHO : float, optional
              Yield of carbon-accumulating phototrophic organisms on CO2 fixed to G3P, in [g COD/g CO2].
              The default is 0.450.
    Y_ATP_HET_ACE : float, optional
              Yield of ATP on acetate fixed to acetyl-CoA, in [g ATP/g COD].
              The default is 39.623.
    Y_CH_NR_HET_ACE : float, optional
              Yield of storage carbohydrates (as polyglucose, PG) on acetate fixed to acetyl-CoA under nutrient-replete condition, in [g COD/g COD].
              The default is 0.625.
    Y_CH_ND_HET_ACE : float, optional
              Yield of storage carbohydrates (as polyglucose, PG) on acetate fixed to acetyl-CoA under nutrient-deplete condition, in [g COD/g COD].
              The default is 0.600.
    Y_LI_NR_HET_ACE : float, optional
              Yield of storage lipids (as triacylglycerol, TAG) on acetate fixed to acetyl-CoA under nutrient-replete condition, in [g COD/g COD].
              The default is 1.105.
    Y_LI_ND_HET_ACE : float, optional
              Yield of storage lipids (as triacylglycerol, TAG) on acetate fixed to acetyl-CoA under nutrient-deplete condition, in [g COD/g COD].
              The default is 0.713.
    Y_X_ALG_HET_ACE : float, optional
              Yield of carbon-accumulating phototrophic organisms on acetate fixed to acetyl-CoA, in [g COD/g COD].
              The default is 0.216.
    Y_ATP_HET_GLU : float, optional
              Yield of ATP on glucose fixed to G6P, in [g ATP/g COD].
              The default is 58.114.
    Y_CH_NR_HET_GLU : float, optional
              Yield of storage carbohydrates (as polyglucose, PG) on glucose fixed to G6P under nutrient-replete condition, in [g COD/g COD].
              The default is 0.917.
    Y_CH_ND_HET_GLU : float, optional
              Yield of storage carbohydrates (as polyglucose, PG) on glucose fixed to G6P under nutrient-deplete condition, in [g COD/g COD].
              The default is 0.880.
    Y_LI_NR_HET_GLU : float, optional
              Yield of storage lipids (as triacylglycerol, TAG) on glucose fixed to G6P under nutrient-replete condition, in [g COD/g COD].
              The default is 1.620.
    Y_LI_ND_HET_GLU : float, optional
              Yield of storage lipids (as triacylglycerol, TAG) on glucose fixed to G6P under nutrient-deplete condition, in [g COD/g COD].
              The default is 1.046.
    Y_X_ALG_HET_GLU : float, optional
              Yield of carbon-accumulating phototrophic organisms on glucose fixed to G6P, in [g COD/g COD].
              The default is 0.317.
    n_dark: float, optional
              Dark growth reduction factor, in [unitless]
              The default is 0.7.              
    f_SI : float, optional
              Production of soluble inerts in hydrolysis, in [g COD/g COD]. 
              The default is 0.0.
    Y_H : float, optional
              Heterotrophic yield coefficient, in[g COD/g COD]. 
              The default is 0.625.
    f_XI_H : float, optional
              Fraction of inert COD generated in heterotrophic biomass lysis, in [g COD/g COD]. 
              The default is 0.1.
    Y_A : float, optional
              Autotrophic yield, in [g COD/g N]. 
              The default is 0.24.
    f_XI_AUT : float, optional
              Fraction of inert COD generated in autotrophic biomass lysis, in [g COD/g COD]. 
              The default is 0.1.    
    K_h : float, optional
              Hydrolysis rate constant, in [d^(-1)]. 
              The default is 3.0.
    eta_NO3 : float, optional
              Reduction factor for anoxic hydrolysis, dimensionless. 
              The default is 0.6.
    eta_fe : float, optional
              Anaerobic hydrolysis reduction factor, dimensionless. 
              The default is 0.4.
    K_O2 : float, optional
              O2 half saturation coefficient for hydrolysis, in [g O2/m^3]. 
              The default is 0.2.
    K_NO3 : float, optional
              Nitrate half saturation coefficient for hydrolysis, in [g N/m^3].
              The default is 0.5.
    K_X : float, optional
              Slowly biodegradable substrate half saturation coefficient for hydrolysis, in [g COD/g COD]. 
              The default is 0.1.
    mu_H : float, optional
              Heterotrophic maximum specific growth rate, in [d^(-1)]. 
              The default is 6.0.
    q_fe : float, optional
              Fermentation maximum rate, in [d^(-1)]. 
              The default is 3.0.
    eta_NO3_H : float, optional
              Reduction factor for anoxic heterotrophic growth, dimensionless.
              The default is 0.8.
    b_H : float, optional
              Lysis and decay rate constant, in [d^(-1)]. 
              The default is 0.4.
    K_O2_H : float, optional
              O2 half saturation coefficient for heterotrophs, in [g O2/m^3].
              The default is 0.2.
    K_F_H : float, optional    
              Fermentable substrate half saturation coefficient for heterotrophic growth (K_F in ASM2d), in [g COD/m^3]. 
              The default is 4.0.
    K_fe : float, optional
              Fermentable substrate half saturation coefficient for fermentation, in [g COD/m^3]. 
              The default is 4.0.
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
              Phosphorus (nutrient) half saturation coefficient for heterotrophs, in [g P/m^3]. 
              The default is 0.01.
    K_ALK_H : float, optional
              Alkalinity half saturation coefficient for heterotrophs, in [mol(HCO3-)/m^3]. (user input unit, converted as C)
              The default is 0.1.          
    mu_AUT : float, optional
              Autotrophic maximum specific growth rate, in [d^(-1)]. 
              The default is 1.0.
    b_AUT : float, optional
              Autotrophic decay rate, in [d^(-1)]. 
              The default is 0.15.
    K_O2_AUT : float, optional
              O2 half saturation coefficient for autotrophs, in [g O2/m^3].
              The default is 0.5.
    K_NH4_AUT : float, optional
              Ammonium (nutrient) half saturation coefficient for autotrophs, in [g N/m^3].
              The default is 1.0.
    K_ALK_AUT : float, optional
              Alkalinity half saturation coefficient for autotrophs, in [mol(HCO3-)/m^3]. (user input unit, converted as C)
              The default is 0.5.
    K_P_AUT : float, optional
              Phosphorus (nutrient) half saturation coefficient for autotrophs, in [g P/m^3].
              The default is 0.01.
    path : str, optional
              Alternative file path for the Petersen matrix.
              The default is None.

    Examples
    --------
    >>> from qsdsan import processes as pc
    >>> cmps = pc.create_pm2asm2d_cmps()
    >>> pm2asm2d = pc.PM2ASM2d()
    >>> pm2asm2d.show()
    PM2ASM2d([photoadaptation, ammonium_uptake, nitrate_uptake_pho, nitrate_uptake_ace, nitrate_uptake_glu, phosphorus_uptake,
         growth_pho, carbohydrate_storage_pho, lipid_storage_pho, carbohydrate_growth_pho, lipid_growth_pho,
         carbohydrate_maintenance_pho, lipid_maintenance_pho, endogenous_respiration_pho,
         growth_ace, carbohydrate_storage_ace, lipid_storage_ace, carbohydrate_growth_ace, lipid_growth_ace,
         carbohydrate_maintenance_ace, lipid_maintenance_ace, endogenous_respiration_ace,
         growth_glu, carbohydrate_storage_glu, lipid_storage_glu, carbohydrate_growth_glu, lipid_growth_glu,
         carbohydrate_maintenance_glu, lipid_maintenance_glu, endogenous_respiration_glu,
         aero_hydrolysis, anox_hydrolysis, anae_hydrolysis,
         hetero_growth_S_F, hetero_growth_S_A, denitri_S_F, denitri_S_A, ferment, hetero_lysis,
         auto_aero_growth, auto_lysis])
    
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

    _shared_params = ('Y_CH_PHO', 'Y_LI_PHO', 'Y_X_ALG_PHO',
               'Y_CH_NR_HET_ACE', 'Y_LI_NR_HET_ACE', 'Y_X_ALG_HET_ACE',
               'Y_CH_NR_HET_GLU', 'Y_LI_NR_HET_GLU', 'Y_X_ALG_HET_GLU')

    _stoichio_params = ('Y_CH_ND_HET_ACE', 'Y_LI_ND_HET_ACE', 'Y_CH_ND_HET_GLU', 'Y_LI_ND_HET_GLU',
                        'f_SI', 'Y_H', 'f_XI_H', 'Y_A', 'f_XI_AUT',
                        *_shared_params)

    _kinetic_params = ('a_c', 'I_n', 'arr_a', 'arr_e', 'beta_1', 'beta_2', 'b_reactor', 'I_opt', 'k_gamma',
                       'K_N', 'K_P', 'K_A', 'K_F', 'rho', 'K_STO', 'f_CH_max', 'f_LI_max', 'm_ATP', 'mu_max',
                       'q_CH', 'q_LI', 'Q_N_max', 'Q_N_min', 'Q_P_max', 'Q_P_min', 'V_NH', 'V_NO', 'V_P', 'exponent',
                       'Y_ATP_PHO', 'Y_ATP_HET_ACE', 'Y_ATP_HET_GLU', *_shared_params, 'n_dark', 'cmps',
                       'K_h', 'eta_NO3', 'eta_fe', 'K_O2', 'K_NO3', 'K_X', 'mu_H', 'q_fe', 'eta_NO3_H', 
                       'b_H', 'K_O2_H', 'K_F_H', 'K_fe', 'K_A_H', 'K_NO3_H', 'K_NH4_H', 'K_P_H', 
                       'K_ALK_H', 'mu_AUT', 'b_AUT', 'K_O2_AUT', 'K_NH4_AUT', 'K_ALK_AUT', 'K_P_AUT')

    def __new__(cls, components=None,
                a_c=0.049, I_n=250, arr_a=1.8e10, arr_e=6842, beta_1=2.90, beta_2=3.50, b_reactor=0.03, I_opt=300, k_gamma=1e-5,
                K_N=0.1, K_P=1.0, K_A=6.3, K_F=6.3, rho=1.186, K_STO=1.566,
                f_CH_max=0.819, f_LI_max=3.249, m_ATP=15.835, mu_max=1.969, q_CH=0.594, q_LI=0.910,
                Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163, V_NH=0.254, V_NO=0.254, V_P=0.016, exponent=4,
                Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
                Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
                Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
                Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
                Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317, n_dark=0.7,
                f_SI=0.0, Y_H=0.625, f_XI_H=0.1, Y_A=0.24, f_XI_AUT=0.1,
                K_h=3.0, eta_NO3=0.6, eta_fe=0.4, K_O2=0.2, K_NO3=0.5, K_X=0.1,
                mu_H=6.0, q_fe=3.0, eta_NO3_H=0.8, b_H=0.4, K_O2_H=0.2, K_F_H=4.0,
                K_fe=4.0, K_A_H=4.0, K_NO3_H=0.5, K_NH4_H=0.05, K_P_H=0.01, K_ALK_H=0.1,
                mu_AUT=1.0, b_AUT=0.15, K_O2_AUT=0.5, K_NH4_AUT=1.0, K_ALK_AUT=0.5, K_P_AUT=0.01,
                path=None, **kwargs):

        if not path: path = _path
        
        self = Processes.load_from_file(path,
                                        components=components,
                                        conserved_for=('COD', 'C', 'N', 'P'),
                                        parameters=cls._stoichio_params,
                                        compile=False)
        
        asm2d_processes = Processes.load_from_file(_path_2,
                                        components=components,
                                        conserved_for=('COD', 'N', 'P', 'charge'),
                                        parameters=cls._stoichio_params,
                                        compile=False)
        self.extend(asm2d_processes)

        if path == _path:
            _p3 = Process('nitrate_uptake_pho',
                           'S_NO -> [?]S_O2 + X_N_ALG',
                           components=components,
                           ref_component='X_N_ALG',
                           conserved_for=('COD', 'C'))

            _p4 = Process('nitrate_uptake_ace',
                           'S_NO + [?]S_A -> [?]S_CO2 + X_N_ALG',
                           components=components,
                           ref_component='X_N_ALG',
                           conserved_for=('COD', 'C'))

            _p5 = Process('nitrate_uptake_glu',
                           'S_NO + [?]S_F -> [?]S_CO2 + X_N_ALG',
                           components=components,
                           ref_component='X_N_ALG',
                           conserved_for=('COD', 'C'))

            self.insert(2, _p3)
            self.insert(3, _p4)
            self.insert(4, _p5)

        self.compile(to_class=cls)

        self.set_rate_function(rhos_pm2asm2d)
        shared_values = (Y_CH_PHO, Y_LI_PHO, Y_X_ALG_PHO,
                         Y_CH_NR_HET_ACE, Y_LI_NR_HET_ACE, Y_X_ALG_HET_ACE,
                         Y_CH_NR_HET_GLU, Y_LI_NR_HET_GLU, Y_X_ALG_HET_GLU)
        stoichio_values = (Y_CH_ND_HET_ACE, Y_LI_ND_HET_ACE, Y_CH_ND_HET_GLU, Y_LI_ND_HET_GLU,
                           f_SI, Y_H, f_XI_H, Y_A, f_XI_AUT,
                           *shared_values)       
        Q_N_min = max(self.Th_Q_N_min, Q_N_min)
        Q_P_min = max(self.Th_Q_P_min, Q_P_min)
        kinetic_values = (a_c, I_n, arr_a, arr_e, beta_1, beta_2, b_reactor, I_opt, k_gamma,
                          K_N, K_P, K_A, K_F, rho, K_STO, f_CH_max, f_LI_max, m_ATP, mu_max,
                          q_CH, q_LI, Q_N_max, Q_N_min, Q_P_max, Q_P_min, V_NH, V_NO, V_P, exponent,
                          Y_ATP_PHO, Y_ATP_HET_ACE, Y_ATP_HET_GLU,
                          *shared_values, n_dark, self._components,
                          K_h, eta_NO3, eta_fe, K_O2, K_NO3, K_X, mu_H, q_fe, eta_NO3_H, 
                          b_H, K_O2_H, K_F_H, K_fe, K_A_H, K_NO3_H, K_NH4_H, K_P_H, 
                          K_ALK_H*12, mu_AUT, b_AUT, K_O2_AUT, K_NH4_AUT, K_ALK_AUT*12, K_P_AUT,
                          )

        dct = self.__dict__
        dct.update(kwargs)
        dct['_parameters'] = dict(zip(cls._stoichio_params, stoichio_values))
        self.rate_function._params = dict(zip(cls._kinetic_params, kinetic_values))

        return self

    def set_parameters(self, **parameters):
        '''Set values to stoichiometric and/or kinetic parameters.'''
        stoichio_only = {k:v for k,v in parameters.items() if k in self._stoichio_params}
        self._parameters.update(stoichio_only)
        if self._stoichio_lambdified is not None:
            self.__dict__['_stoichio_lambdified'] = None
        if 'Q_N_min' in parameters.keys():
            if parameters['Q_N_min'] < self.Th_Q_N_min:
                raise ValueError(f'Value for Q_N_min must not be less than the '
                                 f'theoretical minimum {self.Th_Q_N_min}')
        if 'Q_P_min' in parameters.keys():
            if parameters['Q_P_min'] < self.Th_Q_P_min:
                raise ValueError(f'Value for Q_P_min must not be less than the '
                                 f'theoretical minimum {self.Th_Q_P_min}')
        self.rate_function.set_param(**parameters)

    @property
    def Th_Q_N_min(self):
        return abs(self.stoichiometry.loc['growth_pho', 'X_N_ALG'])*1.001

    @property
    def Th_Q_P_min(self):
        return abs(self.stoichiometry.loc['growth_pho', 'X_P_ALG'])*1.001