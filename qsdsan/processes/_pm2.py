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

from thermosteam import settings
from qsdsan import Component, Components, Process, Processes, CompiledProcesses
from qsdsan.utils import ospath, data_path
import numpy as np

__all__ = ('create_pm2_cmps', 'PM2')

_path = ospath.join(data_path, 'process_data/_pm2.tsv')

#%%
# =============================================================================
# PM2-specific components
# =============================================================================

def create_pm2_cmps(set_thermo=True):
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

    # X_PG (g COD/m^3)
    X_PG = cmps.X_GAO_Gly.copy('X_PG')
    X_PG.description = 'Concentration of stored carbohydrates'
    X_PG.formula = 'CH2O'
    X_PG.f_BOD5_COD = X_PG.f_uBOD_COD = None

    # X_TAG (g COD/m^3)
    X_TAG = cmps.X_GAO_Gly.copy('X_TAG')
    X_TAG.description = 'Concentration of stored lipids'
    X_TAG.formula = 'CH1.92O0.118'
    X_TAG.f_BOD5_COD = X_TAG.f_uBOD_COD = None

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

    # S_G (g COD/m^3)
    S_G = Component.from_chemical(ID = 'S_G',
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
    X_N_ALG.description = 'Concentration of stored nitrogen in microalgal cell'
    X_N_ALG.measured_as = 'N'
    X_N_ALG.i_C = X_N_ALG.i_P = X_N_ALG.i_COD = X_N_ALG.f_BOD5_COD = X_N_ALG.f_uBOD_COD = X_N_ALG.f_Vmass_Totmass = 0
    X_N_ALG.i_mass = 1

    # X_P_ALG (g P/m^3)
    X_P_ALG = cmps.X_B_Subst.copy('X_P_ALG')
    X_P_ALG.description = 'Concentration of stored phosphorus in microalgal cell'
    X_P_ALG.measured_as = 'P'
    X_P_ALG.i_C = X_P_ALG.i_N = X_P_ALG.i_COD = X_P_ALG.f_BOD5_COD = X_P_ALG.f_uBOD_COD = X_P_ALG.f_Vmass_Totmass = 0
    X_P_ALG.i_mass = 1

    cmps_pm2 = Components([X_CHL, X_ALG, X_PG, X_TAG, S_CO2, S_A, S_G,
                           S_O2, S_NH, S_NO, S_P, X_N_ALG, X_P_ALG, cmps.H2O])

    cmps_pm2.default_compile()

    if set_thermo: settings.set_thermo(cmps_pm2)
    return cmps_pm2

# create_pm2_cmps()

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
    :param X_TSS: total biomass concentration [g TSS/m^3]
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
    :param X_carbon: carbon content of cells [g C/m^3]
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
    :param X_carbon: carbon content of cells [g C/m^3]
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
    :param f_CH: ratio of stored carbohydrates to cells (X_PG / X_ALG) [g COD/g COD]
    :param f_LI: ratio of stored lipids to cells (X_TAG / X_ALG) [g COD/g COD]
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
    :param f_CH: ratio of stored carbohydrates to cells (X_PG / X_ALG) [g COD/g COD]
    :param f_LI: ratio of stored lipids to cells (X_TAG / X_ALG) [g COD/g COD]
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
    :param response: f_I (irradiance response function, calculated from 'irrad_response' method), acetate_response (monod(S_A, K_A, 1)), or glucose_response (monod(S_G, K_G, 1)) [unitless]
    :param saturation: 1 - (f / f_max) ** beta (calculated from 'storage_saturation' method) [unitless]
    :param storage_rate: q_CH or q_LI [g COD/g COD/d]
    :return: storage rate [g COD/m^3/d]
    '''
    return storage_rate * saturation * (1 - f_np) * response * X_ALG

def rhos_pm2(state_arr, params):

    # extract values of state variables
    c_arr = state_arr[:14]
    temp = state_arr[15]
    light = state_arr[16]     # imported from input file assumed

    # Q = state_arr[14]                             # Flow rate
    # t = state_arr[15]                             # time

    # light = calc_irrad(t)                         # when to use calculated light (I_0)

    X_CHL, X_ALG, X_PG, X_TAG, S_CO2, S_A, S_G, S_O2, S_NH, S_NO, S_P, X_N_ALG, X_P_ALG, H2O = c_arr

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
    K_G = params['K_G']
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

# intermediate variables
    f_CH = ratio(X_PG, X_ALG, 0, f_CH_max)
    f_LI = ratio(X_TAG, X_ALG, 0, f_LI_max)

    # Q_N = ratio(X_N_ALG, X_ALG, Q_N_min, Q_N_max)
    # Q_P = ratio(X_P_ALG, X_ALG, Q_P_min, Q_P_max)

    alg_iN, alg_iP = cmps.X_ALG.i_N, cmps.X_ALG.i_P
    Q_N = ratio(X_N_ALG+X_ALG*alg_iN, X_ALG, Q_N_min, Q_N_max)
    Q_P = ratio(X_P_ALG+X_ALG*alg_iP, X_ALG, Q_P_min, Q_P_max)

    idx = cmps.indices(['X_ALG', 'X_PG', 'X_TAG'])
    X_bio = np.array([X_ALG, X_PG, X_TAG])
    X_TSS = sum(X_bio * cmps.i_mass[idx])
    X_carbon = sum(X_bio * cmps.i_C[idx])

    i_avg = attenuation(light, X_TSS, a_c, b_reactor)
    f_I = irrad_response(i_avg, X_CHL, X_carbon, I_n, I_opt)
    dark_response = max(f_I, n_dark)
    acetate_response = monod(S_A, K_A, 1)
    glucose_response = monod(S_G, K_G, 1)

    f_np = min(droop(Q_N, Q_N_min, exponent), droop(Q_P, Q_P_min, exponent))
    f_temp = temperature(temp, arr_a, arr_e)

    f_sat_CH = storage_saturation(f_CH, f_CH_max, beta_1)
    f_sat_LI = storage_saturation(f_LI, f_LI_max, beta_2)

    max_total_growth_rho = max_total_growth(X_ALG, mu_max, f_np, f_temp)
    max_maintenance_rho = max_total_maintenance(X_ALG, m_ATP)
    # light = calc_irrad(t)

    # calculate kinetic rate values
    rhos = np.empty(30)

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

    return rhos

#%%
# =============================================================================
# PM2 class
# =============================================================================

class PM2(CompiledProcesses):
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
    K_G : float, optional
              Organic carbon half-saturation constant (glucose); assumes K_A = K_G, in [g COD/m^3].
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
    path : str, optional
              Alternative file path for the Petersen matrix.
              The default is None.

    Examples
    --------
    >>> from qsdsan import processes as pc
    >>> cmps = pc.create_pm2_cmps()
    >>> pm2 = pc.PM2()
    >>> pm2.show()
    PM2([photoadaptation, ammonium_uptake, nitrate_uptake_pho, nitrate_uptake_ace, nitrate_uptake_glu, phosphorus_uptake,
         growth_pho, carbohydrate_storage_pho, lipid_storage_pho, carbohydrate_growth_pho, lipid_growth_pho,
         carbohydrate_maintenance_pho, lipid_maintenance_pho, endogenous_respiration_pho,
         growth_ace, carbohydrate_storage_ace, lipid_storage_ace, carbohydrate_growth_ace, lipid_growth_ace,
         carbohydrate_maintenance_ace, lipid_maintenance_ace, endogenous_respiration_ace,
         growth_glu, carbohydrate_storage_glu, lipid_storage_glu, carbohydrate_growth_glu, lipid_growth_glu,
         carbohydrate_maintenance_glu, lipid_maintenance_glu, endogenous_respiration_glu])

    >>> # Evaluate the rate of reaction at initial condition
    >>> import numpy as np
    >>> init_cond = {
    ...     'X_CHL':2.81,
    ...     'X_ALG':561.57,
    ...     'X_PG':13.74,
    ...     'X_TAG':62.22,
    ...     'S_CO2':30.0,
    ...     'S_A':5.0,
    ...     'S_G':5.0,
    ...     'S_O2':20.36,
    ...     'S_NH':25,
    ...     'S_NO':9.30,
    ...     'S_P':0.383,
    ...     'X_N_ALG':3.62,
    ...     'X_P_ALG':12.60,
    ...     }
    >>> state_arr = np.append(cmps.kwarray(init_cond), [1000, 298, 112.6])   # flowrate, temperature, & irradiance
    >>> pm2.rate_function(state_arr)
    array([  4.437, 142.045,   0.562,   0.562,   0.562,   2.48 , 304.319,
           193.505,   8.856,  24.737,  79.042,   2.093,   7.992,  67.419,
           186.095, 110.903,   5.076,  15.127,  32.669,   2.455,   9.375,
            45.795, 186.075, 110.903,   5.076,  15.126,  32.691,   2.456,
             9.378,  45.822])

    >>> pm2.set_parameters(I_opt = 200) # Change optimal irradiance
    >>> pm2.rate_function(state_arr)
    array([  4.437, 142.045,   0.562,   0.562,   0.562,   2.48 , 299.409,
           209.95 ,   9.609,  34.175, 109.197,   2.093,   7.992,  67.419,
           171.898, 110.903,   5.076,  19.621,  42.373,   2.455,   9.375,
            45.795, 171.874, 110.903,   5.076,  19.618,  42.4  ,   2.456,
             9.378,  45.822])
    '''

    _shared_params = ('Y_CH_PHO', 'Y_LI_PHO', 'Y_X_ALG_PHO',
               'Y_CH_NR_HET_ACE', 'Y_LI_NR_HET_ACE', 'Y_X_ALG_HET_ACE',
               'Y_CH_NR_HET_GLU', 'Y_LI_NR_HET_GLU', 'Y_X_ALG_HET_GLU')

    _stoichio_params = ('Y_CH_ND_HET_ACE', 'Y_LI_ND_HET_ACE', 'Y_CH_ND_HET_GLU', 'Y_LI_ND_HET_GLU',
                        *_shared_params)

    _kinetic_params = ('a_c', 'I_n', 'arr_a', 'arr_e', 'beta_1', 'beta_2', 'b_reactor', 'I_opt', 'k_gamma',
                       'K_N', 'K_P', 'K_A', 'K_G', 'rho', 'K_STO', 'f_CH_max', 'f_LI_max', 'm_ATP', 'mu_max',
                       'q_CH', 'q_LI', 'Q_N_max', 'Q_N_min', 'Q_P_max', 'Q_P_min', 'V_NH', 'V_NO', 'V_P', 'exponent',
                       'Y_ATP_PHO', 'Y_ATP_HET_ACE', 'Y_ATP_HET_GLU', *_shared_params, 'n_dark', 'cmps')

    def __new__(cls, components=None,
                a_c=0.049, I_n=250, arr_a=1.8e10, arr_e=6842, beta_1=2.90, beta_2=3.50, b_reactor=0.03, I_opt=300, k_gamma=1e-5,
                K_N=0.1, K_P=1.0, K_A=6.3, K_G=6.3, rho=1.186, K_STO=1.566,
                f_CH_max=0.819, f_LI_max=3.249, m_ATP=15.835, mu_max=1.969, q_CH=0.594, q_LI=0.910,
                Q_N_max=0.417, Q_N_min=0.082, Q_P_max=0.092, Q_P_min=0.0163, V_NH=0.254, V_NO=0.254, V_P=0.016, exponent=4,
                Y_ATP_PHO=55.073, Y_CH_PHO=0.754, Y_LI_PHO=0.901, Y_X_ALG_PHO=0.450,
                Y_ATP_HET_ACE=39.623, Y_CH_NR_HET_ACE=0.625, Y_CH_ND_HET_ACE=0.600,
                Y_LI_NR_HET_ACE=1.105, Y_LI_ND_HET_ACE=0.713, Y_X_ALG_HET_ACE=0.216,
                Y_ATP_HET_GLU=58.114, Y_CH_NR_HET_GLU=0.917, Y_CH_ND_HET_GLU=0.880,
                Y_LI_NR_HET_GLU=1.620, Y_LI_ND_HET_GLU=1.046, Y_X_ALG_HET_GLU=0.317, n_dark=0.7,
                path=None, **kwargs):

        if not path: path = _path
        self = Processes.load_from_file(path,
                                        components=components,
                                        conserved_for=('COD', 'C', 'N', 'P'),
                                        parameters=cls._stoichio_params,
                                        compile=False)

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
                           'S_NO + [?]S_G -> [?]S_CO2 + X_N_ALG',
                           components=components,
                           ref_component='X_N_ALG',
                           conserved_for=('COD', 'C'))

            self.insert(2, _p3)
            self.insert(3, _p4)
            self.insert(4, _p5)

        self.compile(to_class=cls)

        self.set_rate_function(rhos_pm2)
        shared_values = (Y_CH_PHO, Y_LI_PHO, Y_X_ALG_PHO,
                         Y_CH_NR_HET_ACE, Y_LI_NR_HET_ACE, Y_X_ALG_HET_ACE,
                         Y_CH_NR_HET_GLU, Y_LI_NR_HET_GLU, Y_X_ALG_HET_GLU)
        stoichio_values = (Y_CH_ND_HET_ACE, Y_LI_ND_HET_ACE, Y_CH_ND_HET_GLU, Y_LI_ND_HET_GLU,
                           *shared_values)
        Q_N_min = max(self.Th_Q_N_min, Q_N_min)
        Q_P_min = max(self.Th_Q_P_min, Q_P_min)
        kinetic_values = (a_c, I_n, arr_a, arr_e, beta_1, beta_2, b_reactor, I_opt, k_gamma,
                          K_N, K_P, K_A, K_G, rho, K_STO, f_CH_max, f_LI_max, m_ATP, mu_max,
                          q_CH, q_LI, Q_N_max, Q_N_min, Q_P_max, Q_P_min, V_NH, V_NO, V_P, exponent,
                          Y_ATP_PHO, Y_ATP_HET_ACE, Y_ATP_HET_GLU,
                          *shared_values, n_dark, self._components)

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
        self.rate_function.set_params(**parameters)

    @property
    def Th_Q_N_min(self):
        return abs(self.stoichiometry.loc['growth_pho', 'X_N_ALG'])*1.001

    @property
    def Th_Q_P_min(self):
        return abs(self.stoichiometry.loc['growth_pho', 'X_P_ALG'])*1.001