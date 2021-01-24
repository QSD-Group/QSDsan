#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 08:20:30 2020

@author: lewisrowles

References:
(1) 	Guidelines for Using Activated Sludge Models; Rieger, L., Ed.; Scientific and technical report; IWA Publ: London, 2013.
(2) 	Henze, M.; Gujer, W.; Mino, T.; Matsuo, T.; Wentzel, M. C.; Marais, G. v. R.; Van Loosdrecht, M. C. M. Activated Sludge Model No.2d, ASM2D. Water Sci. Technol. 1999, 39 (1), 165–182. https://doi.org/10.2166/wst.1999.0036.

"""


import numpy as np
from scipy.integrate import odeint
import sympy as sp
import matplotlib.pyplot as plt
from sympy.parsing.sympy_parser import parse_expr
from datetime import datetime

start_time = datetime.now()


#conversion of parameters 

#SVFA = SAc + SProp
# XU = XU_Inf + XU_OHO_E + XU_PAO_E	
# XANO = XAOO + XNOO
# XPAO_PP = XPAO_PP_Lo + XPAO_PP_Hi
# XS = XB_Subst + CB_Subst + CB_BAP + CB_UAP
# XPHA = XPAO_PHA + XPAO_Gly

# define influent concentrations
SN20 = 15.0
SNH40 = 16.0
SNO30 = 0.0001
XMeOH0 = 0.0001
SPO40 = 3.6
XPAO_PP0 =0.0
XMeP0 = 0.0
XOHO0 = 30.0
XU0 = 25.0
SVFA0 = 20
SAlk0 = 5.0
XPAO0 = 0.0001
XS0 = 125.0
XTSS0 = 180.0
SF0 = 30.0
SU0 = 30.0
XANO0 = 0.0001
XPHA0 = 0.0001


# define influent concentrations
s_n20 = SN20
s_nh40 = SNH40
s_no30 = SNO30
x_meoh0 = XMeOH0
s_po40 = SPO40 
x_pp0 = XPAO_PP0
x_mep0 = XMeP0
x_h0 = XOHO0
x_i0 = XU0
s_a0 = SVFA0
s_alk0 = SAlk0
x_pao0 = XPAO0 
x_s0 = XS0 
x_tss0 = XTSS0 
s_f0 = SF0
s_i0 = SU0
x_aut0 = XANO0
x_pha0 = XPHA0

C0 = [s_n20, s_nh40, s_no30,  x_meoh0, s_po40, x_pp0, x_mep0, x_h0, x_i0, s_a0, s_alk0, x_pao0, x_s0, x_tss0, s_f0, s_i0, x_aut0, x_pha0]

# set flowrate, volume, and oxygen concentration to be supplied
flowrate = 10000
V_aerobic = 5000
O2_conc = 2.0
s_o2 = O2_conc
calculate_materials = 'yes'

# table 9 ASM2d typical conversion factors for conservation equation
i_nsi = 0.01
i_nsf = 0.03
i_nxi = 0.02
i_nxs = 0.04
i_nbm = 0.07
i_psi = 0.00
i_psf = 0.01
i_pxi = 0.01
i_pxs = 0.01
i_pbm = 0.02
i_tssxi = 0.75
i_tssxs = 0.75
i_tssbm = 0.9

# table 9 ASM2d typical stoichiometric parameters 
f_si = 0
y_h = 0.625
f_xi = 0.1
y_pao = 0.625
y_po4 = 0.40
y_pha = 0.2
y_a = 0.24

#table 10 typical values for kinetic parameters of ASM2d at 20 degC
k_h = 3.00
n_no3 = 0.60
n_fe = 0.40
k_o2 = 0.20
k_no3 = 0.50
k_x = 0.10
u_h = 6.00
q_fe = 3.00
b_h = 0.40
k_f = 4
k_fe = 4
k_a = 4
k_nh4 = 0.05
k_p = 0.01
k_alk = 0.10
q_pha = 3.00
q_pp = 1.50
u_pao = 1.00
b_pao = 0.20
b_pp = 0.2
b_pha = 0.2
k_ps = 0.20
k_pp = 0.01
k_max = 0.34
k_ipp = 0.02
k_pha = 0.01
u_aut = 1.00
b_aut = 0.15
k_pre = 1.00
k_red = 0.60


#from ASM2d Tables 1 and 2 input as matrices
#stoichiometric coefficients can be computed using Conservation equation 2 and calling on these tables

#table 1 conversation factors i_ci

conversion_factors = {
    's_o2': {
        'cod': -1,
        'n': 0,
        'p': 0,
        'charge': 0,
        'mass': 0
    },
    's_f': {
        'cod': 1,
        'n': i_nsf,
        'p': i_psf,
        'charge': 0,
        'mass': 0
    },
    's_a': {
        'cod': 1,
        'n': 0,
        'p': 0,
        'charge': -1/64,
        'mass': 0
    },
    's_nh4': {
        'cod': 0,
        'n': 1,
        'p': 0,
        'charge': 1/14,
        'mass': 0
    },
    's_no3': {
        'cod': -64/14,
        'n': 1,
        'p': 0,
        'charge': -1/14,
        'mass': 0
    },
    's_po4': {
        'cod': 0,
        'n': 0,
        'p': 1,
        'charge': -1.5/31,
        'mass': 0
    },
    's_i': {
        'cod': 1,
        'n': i_nsi,
        'p': i_psi,
        'charge': 0,
        'mass': 0
    },
    's_alk': {
        'cod': 0,
        'n': 0,
        'p': 0,
        'charge': -1,
        'mass': 0
    },
    's_n2': {
        'cod': -24/14,
        'n': 1,
        'p': 0,
        'charge': 0,
        'mass': 0
    },
    'x_i': {
        'cod': 1,
        'n': i_nxi,
        'p': i_pxi,
        'charge': 0,
        'mass': i_tssxi
    },
    'x_s': {
        'cod': 1,
        'n': i_nxs,
        'p': i_pxs,
        'charge': 0,
        'mass': i_tssxs
    },
    'x_h': {
        'cod': 1,
        'n': i_nbm,
        'p': i_pbm,
        'charge': 0,
        'mass': i_tssbm
    },
    'x_pao': {
        'cod': 1,
        'n': i_nbm,
        'p': i_pbm,
        'charge': 0,
        'mass': i_tssbm
    },
    'x_pp': {
        'cod': 0,
        'n': 0,
        'p': 1,
        'charge': -1/31,
        'mass': 3.23
    },
    'x_pha': {
        'cod': 1,
        'n': 0,
        'p': 0,
        'charge': 0,
        'mass': 0.6
    },
    'x_aut': {
        'cod': 1,
        'n': i_nbm,
        'p': i_nbm,
        'charge': 0,
        'mass': i_tssbm
    },
    'x_tss': {
        'cod': 0,
        'n': 0,
        'p': 0,
        'charge': 0,
        'mass': -1
    },
    'x_meoh': {
        'cod': 0,
        'n': 0,
        'p': 0,
        'charge': 0,
        'mass': 1
    },
    'x_mep': {
        'cod': 0,
        'n': 0,
        'p': 0.205,
        'charge': 0,
        'mass': 1
    }
}



v_1nh4 = -((f_si)*i_nsi + (-1)*i_nxs + (1 - f_si)*i_nsf)
v_1po4 = -((f_si)*i_psi + (-1)*i_pxs + (1 - f_si)*i_psf) 
v_1alk = ((v_1po4)*-0.04838709677419355 + (v_1nh4)*0.07142857142857142)
v_1tss = ((-1)*i_tssxs)
v_2nh4 = -((f_si)*i_nsi + (-1)*i_nxs + (1 - f_si)*i_nsf)
v_2po4 = -((f_si)*i_psi + (-1)*i_pxs + (1 - f_si)*i_psf)
v_2alk = ((v_2po4)*-0.04838709677419355 + (v_2nh4)*0.07142857142857142)
v_2tss = ((-1)*i_tssxs)
v_3nh4 = -((f_si)*i_nsi + (-1)*i_nxs + (1 - f_si)*i_nsf)
v_3po4 = -((f_si)*i_psi + (-1)*i_pxs + (1 - f_si)*i_psf)
v_3alk = ((v_3po4)*-0.04838709677419355 + (v_3nh4)*0.07142857142857142)
v_3tss = ((-1)*i_tssxs)
v_12no3 = 0
v_13o2 = ((1)*1 + (-1/y_h)*1)
v_14no3 = 0
v_15po4 = -((1 - f_xi)*i_pxs + (-1)*i_pbm + (f_xi)*i_pxi)
v_18nh4 = -((1)*i_nbm + (1/y_a)*1)
v_19nh4 = -((-1)*i_nbm + (1 - f_xi)*i_nxs + (f_xi)*i_nxi)
v_19po4 = -((-1)*i_nbm + (1 - f_xi)*i_pxs + (f_xi)*i_pxi)
v_20alk = ((-1)*-0.04838709677419355)
v_21alk = -((1)*-0.04838709677419355)






#table 2 to 6 stocihiometry
stoichiometry = {
    'aerobic_hydrolysis' : {
        's_f' : (1 - f_si),
        's_nh4' : v_1nh4,
        's_po4' : v_1po4,
        's_i' : f_si,
        's_alk' : v_1alk,
        'x_s' : -1,
        'x_tss' : v_1tss
    },
    'anoxic_hydrolysis' : {
        's_f' : (1 - f_si),
        's_nh4' : v_2nh4,
        's_po4' : v_2po4,
        's_i' : f_si,
        's_alk' : v_2alk,
        'x_s' : -1,
        'x_tss' : v_2tss
    },
    'anaerobic_hydrolysis' : {
        's_f' : (1 - f_si),
        's_nh4' : v_3nh4,
        's_po4' : v_3po4,
        's_i' : f_si,
        's_alk' : v_3alk,
        'x_s' : -1,
        'x_tss' : v_3tss
    },
    'aerobic_growth_on_s_f' : {
        's_o2' : (1 - (1 / y_h)),
        's_f' : -(1 / y_h),
        's_a' : 0,
        's_no3' : 0,
        's_n2' : 0,
        'x_i' : 0,
        'x_s' : 0,
        'x_h' : -1
    },
    'aerobic_growth_on_s_a' : {
        's_o2' : (1 - (1 / y_h)),
        's_f' : 0,
        's_a' : - (1 / y_h),
        's_no3' : 0,
        's_n2' : 0,
        'x_i' : 0,
        'x_s' : 0,
        'x_h' : 1
    },
    'anoxic_growth_on_s_f' : {
        's_o2' : 0,
        's_f' : - (1 / y_h),
        's_a' : 0,
        's_no3' : (-(1 - y_h) / (2.86 * y_h)),
        's_n2' : ((1 - y_h) / (2.86 * y_h)),
        'x_i' : 0,
        'x_s' : 0,
        'x_h' : 1
    },
    'denitrification' : {
        's_o2' : 0,
        's_f' : 0,
        's_a' : - (1 / y_h),
        's_no3' : (-(1 - y_h) / (2.86 * y_h)),
        's_n2' : ((1 - y_h) / (2.86 * y_h)),
        'x_i' : 0,
        'x_s' : 0,
        'x_h' : 1
    },
    'fermentation' : {
        's_o2' : 0,
        's_f' : -1,
        's_a' : 1,
        's_no3' : 0,
        's_n2' : 0,
        'x_i' : 0,
        'x_s' : 0,
        'x_h' : 0
    },
    'lysis' : {
        's_o2' : 0,
        's_f' : 0,
        's_a' : 0,
        's_no3' : 0,
        's_n2' : 0,
        'x_i' : f_xi,
        'x_s' : 1 - f_xi,
        'x_h' : -1
    },
    'storage_of_x_pha' : {
        's_o2' : 0,
        's_a' : -1,
        's_n2' : 0,
        's_no3' : 0,
        's_po4' : y_po4,
        'x_i' : 0,
        'x_s' : 0,
        'x_pao' : 0,
        'x_pp' : -y_po4,
        'x_pha' : 0,
    },
    'aerobic_storage_of_x_pp' : {
        's_o2' : -y_pha,
        's_a' : 0,
        's_n2' : 0,
        's_no3' : 0,
        's_po4' : -1,
        'x_i' : 0,
        'x_s' : 0,
        'x_pao' : 0,
        'x_pp' : 1,
        'x_pha' : -y_pha,
    },
    'anoxic_storage_of_x_pp' : {
        's_o2' : 0,
        's_a' : 0,
        's_n2' : -v_12no3,
        's_no3' : v_12no3,
        's_po4' : -1,
        'x_i' : 0,
        'x_s' : 0,
        'x_pao' : 0,
        'x_pp' : 1,
        'x_pha' : -y_pha,
    },
    'aerobic_growth_of_x_pao' : {
        's_o2' : v_13o2,
        's_a' : 0,
        's_n2' : 0,
        's_no3' : 0,
        's_po4' : -i_pbm,
        'x_i' : 0,
        'x_s' : 0,
        'x_pao' : 1,
        'x_pp' : 0,
        'x_pha' : (-1 / y_h),
    },
    'anoxic_growth_of_x_pao' : {
        's_o2' : 0,
        's_a' : 0,
        's_n2' : -v_14no3,
        's_no3' : v_14no3,
        's_po4' : -i_pbm,
        'x_i' : 0,
        'x_s' : 0,
        'x_pao' : 1,
        'x_pp' : 0,
        'x_pha' : (-1 / y_h),
    },
    'lysis_of_x_pao' : {
        's_o2' : 0,
        's_a' : 0,
        's_n2' : 0,
        's_no3' : 0,
        's_po4' : v_15po4,
        'x_i' : f_xi,
        'x_s' : 1 - f_xi,
        'x_pao' : -1,
        'x_pp' : 0,
        'x_pha' : 0,
    },
    'lysis_of_x_pp' : {
        's_o2' : 0,
        's_a' : 0,
        's_n2' : 0,
        's_no3' : 0,
        's_po4' : 1,
        'x_i' : 0,
        'x_s' : 0,
        'x_pao' : 0,
        'x_pp' : -1,
        'x_pha' : 0,
    },
    'lysis_of_x_pha' : {
        's_o2' : 0,
        's_a' : 1,
        's_n2' : 0,
        's_no3' : 0,
        's_po4' : 0,
        'x_i' : 0,
        'x_s' : 0,
        'x_pao' : 0,
        'x_pp' : 0,
        'x_pha' : -1,
    },
    'aerobic_growth_of_x_aut' : {
        's_o2' : (-(4.57 - y_a) / y_a),
        's_nh4' : v_18nh4,
        's_no3' : (1 / y_a),
        's_po4' : -i_pbm,
        'x_i' : 0,
        'x_s' : 0,
        'x_aut' : 1,
    },
    'lysis_of_x_aut' : {
        's_o2' : 0,
        's_nh4' : v_19nh4,
        's_no3' : 0,
        's_po4' : v_19po4,
        'x_i' : f_xi,
        'x_s' : 1 - f_xi,
        'x_aut' : -1,
    },
    'precipitation' : {
        's_po4' : -1,
        's_alk' : v_20alk,
        'x_meoh' : -3.45,
        'x_mep' : 4.87,
        'x_tss' : 1.42,
    },
    'redissolution' : {
    's_po4': 1,
    's_alk': v_21alk,
    'x_meoh': 3.45,
    'x_mep': -4.87,
    'x_tss': -1.42,
    }
}

s_o2, x_s, x_h, s_no3, s_f, s_a, s_nh4, s_po4, s_alk, x_pp, x_pao, x_pha,  x_aut, x_meoh, x_mep = sp.symbols('s_o2, x_s, x_h, s_no3, s_f, s_a, s_nh4, s_po4, s_alk, x_pp, x_pao, x_pha,  x_aut, x_meoh, x_mep')

# #table 7 process rate equations
process_rate_equations = {
    'aerobic_hydrolysis' :
        (k_h * (s_o2 / (k_o2 + s_o2)) * ((x_s / x_h) / (k_x + (x_s / x_h))) * x_h),
    'anoxic_hydrolysis' :
        (k_h * n_no3 * (k_o2 / (k_o2 + s_o2)) * (s_no3 / (k_no3 + s_no3)) * ((x_s / x_h) / (k_x + (x_s / x_h))) * x_h),
    'anaerobic_hydrolysis' :
        (k_h * n_fe * (k_o2 / (k_o2 + s_o2)) * (k_no3 / (k_no3 + s_no3)) * ((x_s / x_h) / (k_x + (x_s / x_h))) * x_h),
    'aerobic_growth_on_s_f' :
        (u_h * (s_o2 / (k_o2 + s_o2)) * (s_f / (k_f + s_f)) * (s_f / (s_f + s_a)) * (s_nh4 / (k_nh4 + s_nh4)) * (s_po4 / (k_p * s_po4)) * (s_alk / (k_alk + s_alk)) * x_h),
    'aerobic_growth_on_S_a' :
        (u_h * (s_o2 / (k_o2 + s_o2)) * (s_a / (k_f + s_f)) * (s_a / (s_f + s_a)) * (s_nh4 / (k_nh4 + s_nh4)) * (s_po4 / (k_p + s_po4)) * (s_alk / (k_alk + s_alk)) * x_h),
    'anoxic_growth_on_s_f' :
        (u_h * n_no3 * (k_o2 / (k_o2 + s_o2)) * (k_no3 / (k_no3 + s_no3)) * (s_f / (k_f + s_f)) * (s_f / (s_f + s_a)) * (s_nh4 / (k_nh4 + s_nh4)) * (s_po4 / (k_p + s_po4)) * (s_alk / (k_alk + s_alk)) * x_h),
    'denitrification' :
        (u_h * n_no3 * (k_o2 / (k_o2 + s_o2)) * (k_no3 / (k_no3 + s_no3)) * (s_a / (k_a + s_a)) * (s_a / (s_f + s_a)) * (s_nh4 / (k_nh4 + s_nh4)) * (s_po4 / (k_p + s_po4)) * (s_alk / (k_alk + s_alk)) * x_h),
    'fermentation' :
        (q_fe * (k_o2 / (k_o2 + s_o2)) * (k_no3 / (k_no3 + s_no3)) * (s_f / (k_f + s_f)) * (s_alk / (k_alk + s_alk)) * x_h),
    'lysis' :
        (b_h * x_h),
    'storage_of_x_pha' :
        (q_pha * (s_a / (k_a + s_a)) * (s_alk / (k_alk + s_alk)) * ((x_pp / x_pao) / (k_pp + (x_pp / x_pao))) * x_pao),
    'aerobic_storage_of_x_pp' :
        (q_pp * (s_o2 / (k_o2 + s_o2)) * (s_po4 / (k_ps + s_po4)) * (s_alk / (k_alk + s_alk)) * ((x_pha / x_pao) / (k_pha + (x_pha / x_pao))) * ((k_max - (x_pp / x_pao)) / (k_pp + k_max - (x_pp / x_pao))) * x_pao),
    'anoxic_storage_of_x_pp' :
        ((q_pp * (s_o2 / (k_o2 + s_o2)) * (s_po4 / (k_ps + s_po4)) * (s_alk / (k_alk + s_alk)) * ((x_pha / x_pao) / (k_pha + (x_pha / x_pao))) * ((k_max - (x_pp / x_pao)) / (k_pp + k_max - (x_pp / x_pao))) * x_pao) * n_no3 * (k_o2 / s_o2) * (s_no3 / (k_no3 + s_no3))),
    'aerobic_growth_of_x_pao' :
        (u_pao * (s_o2 / (k_o2 + s_o2)) * (s_nh4 / (k_nh4 + s_nh4)) * (s_alk / (k_alk + s_alk)) * ((x_pha / x_pao) / (k_pha + (x_pha / x_pao))) * x_pao),
    'anoxic_growth_of_x_pao' :
        (u_pao * (s_o2 / (k_o2 + s_o2)) * (s_nh4 / (k_nh4 + s_nh4)) * (s_alk / (k_alk + s_alk)) * ((x_pha / x_pao) / (k_pha + (x_pha / x_pao))) * x_pao * n_no3 * (k_o2 / s_no3) * (s_no3 / (k_no3 + s_no3))),
    'lysis_of_x_pao' :
        (b_pao * x_pao * (s_alk / (k_alk + s_alk))),
    'lysis_of_x_pp' :
        (b_pp * x_pp * (s_alk / (k_alk + s_alk))),
    'lysis_of_x_pha' :
        (b_pha * x_pha * (s_alk / (k_alk + s_alk))),
    'aerobic_growth_of_x_aut' :
        (u_aut * (s_o2 / (k_o2 + s_o2)) * (s_nh4 / (k_nh4 + s_nh4)) * (s_po4 / (k_p * s_po4)) * (s_alk / (k_alk + s_alk)) * x_aut),
    'lysis_of_x_aut' :
        (b_aut * x_aut),
    'precipitation' :
        (k_pre * s_po4 * x_meoh),
    'redissolution' :
        (k_red * x_mep * (s_alk / (k_alk + s_alk)))
}


#function that will return conservation equation for a component (i) in process (j) which is equal to zero
def conservation_equation_solver(process, component):
    call_process = stoichiometry[process]
    call_component =  dict((k, v[component]) for k, v in conversion_factors.items())
    components_from_process = set(call_component.keys()).intersection(set(call_process))
    v_process_component = []
    for overlap in components_from_process:
        a = call_component[overlap]
        b = call_process[overlap]
        v_process_component.append('(' +str(b) + ')' + '*' +  str(a))
    v_process_component_string = ' + '.join(v_process_component)
    return v_process_component_string

#function that return rate of production for the specified component (i) over all processes(j)
def rate_of_production(component):
    component_stoichiometry = []
    component_process_rate_equation = []
    for process_num, process_info in stoichiometry.items():
        for process_component in process_info:
            if component == process_component:
                component_stoichiometry.append(process_num)
    for process_num, process_rate_equation in process_rate_equations.items():
        if process_num in component_stoichiometry:
            a = process_rate_equation
            b = stoichiometry[process_num][component]
            component_process_rate_equation.append('(' +str(a) + ')' + '*' +  str(b))
            component_process_rate_equation_string = ' + '.join(component_process_rate_equation)
    return component_process_rate_equation_string


components = ['s_n2', 's_nh4', 's_no3', 'x_meoh', 's_po4', 'x_pp', 'x_mep', 'x_h', 'x_i', 's_a', 's_alk', 'x_pao', 'x_s', 'x_tss', 's_f', 's_i', 'x_aut', 'x_pha']
components_0 = {'s_n2' : s_n20, 's_nh4' : s_nh40, 's_no3' : s_no30, 'x_meoh' : x_meoh0, 's_po4' : s_po40, 'x_pp' : x_pp0, 'x_mep' : x_mep0, 'x_h' : x_h0, 'x_i' : x_i0, 's_a' : s_a0, 's_alk' : s_alk0, 'x_pao' : x_pao0, 'x_s' : x_s0, 'x_tss' : x_tss0, 's_f' : s_f0, 's_i' : s_i0, 'x_aut' : x_aut0, 'x_pha' : x_pha0}
args_components = sp.symbols('s_n2  s_nh4, s_no3  x_meoh s_po4 x_pp x_mep x_h x_i s_a s_alk x_pao x_s x_tss s_f s_i x_aut x_pha')

def reactor(reactor_type, C):
    equations = {}
    s_n2 , s_nh4, s_no3,  x_meoh, s_po4, x_pp, x_mep, x_h, x_i, s_a, s_alk, x_pao, x_s, x_tss, s_f, s_i, x_aut, x_pha = C
    if reactor_type == 'batch':
        for component in components:
            equations.update({component: rate_of_production(component)})

    if reactor_type == 'cstr':
        for component in components:
            volume = V_aerobic
            equations.update({component: (str((flowrate/volume)) + '*' + str(components_0.get(component)) + '+' + str((flowrate/volume)) + '*' + component + '+' + str(rate_of_production(component)))})
    # write function to replace s_o2 with value
    #equations = [equations.replace(key, value) for key  ]
    for key, value in equations.items():
        value = value.replace('s_o2', str(O2_conc))
        equations[key] = value 
    ds_n2 = equations['s_n2']
    ds_nh4 = equations['s_nh4']
    ds_no3 = equations[ 's_no3']
    dx_meoh = equations['x_meoh']
    ds_po4 = equations['s_po4']
    dx_pp = equations['x_pp']
    dx_mep = equations['x_mep']
    dx_h = equations['x_h']
    dx_i = equations['x_i']
    ds_a = equations['s_a']
    ds_alk = equations['s_alk']
    dx_pao = equations['x_pao']
    dx_s = equations['x_s']
    dx_tss = equations['x_tss']
    ds_f = equations['s_f']
    ds_i = equations['s_i']
    dx_aut = equations['x_aut']
    dx_pha = equations['x_pha']
    equations = [ds_n2, ds_nh4, ds_no3, dx_meoh, ds_po4, dx_pp, dx_mep, dx_h, dx_i, ds_a, ds_alk, dx_pao, dx_s, dx_tss, ds_f, ds_i, dx_aut, dx_pha]
    syms = [parse_expr(item) for item in equations]
    funcs = [sp.lambdify(args_components, f) for f in syms]
    return funcs

funcs_batch = reactor('batch', C0)
funcs_cstr = reactor('cstr', C0)

def lambdified_odes_batch(C, t):
    s_n2 , s_nh4, s_no3,  x_meoh, s_po4, x_pp, x_mep, x_h, x_i, s_a, s_alk, x_pao, x_s, x_tss, s_f, s_i, x_aut, x_pha = C
    xdot = [f(s_n2 , s_nh4, s_no3,  x_meoh, s_po4, x_pp, x_mep, x_h, x_i, s_a, s_alk, x_pao, x_s, x_tss, s_f, s_i, x_aut, x_pha) for f in funcs_batch]
    return xdot
   
def lambdified_odes_cstr(C, t):
    s_n2 , s_nh4, s_no3,  x_meoh, s_po4, x_pp, x_mep, x_h, x_i, s_a, s_alk, x_pao, x_s, x_tss, s_f, s_i, x_aut, x_pha = C
    xdot = [f(s_n2 , s_nh4, s_no3,  x_meoh, s_po4, x_pp, x_mep, x_h, x_i, s_a, s_alk, x_pao, x_s, x_tss, s_f, s_i, x_aut, x_pha) for f in funcs_cstr]
    return xdot

def run_reactor(reactor_type, time_start, time_stop, time_points, C0, export_results):
    t = np.linspace(time_start, time_stop, time_points)
    if reactor_type == 'batch':
        C = odeint(lambdified_odes_batch, C0, t)        
    if reactor_type == 'cstr':
        C = odeint(lambdified_odes_cstr, C0, t)        
    s_n2 = C[:, 0]
    s_nh4 = C[:, 1]
    s_no3 = C[:, 2]
    x_meoh = C[:, 3]
    s_po4 = C[:, 4]
    x_pp = C[:, 5]
    x_mep = C[:, 6]
    x_h = C[:, 7]
    x_i = C[:, 8]
    s_a = C[:, 9]
    s_alk = C[:, 10]
    x_pao = C[:, 11]
    x_s = C[:, 12]
    x_tss = C[:, 13]
    s_f = C[:, 14]
    s_i = C[:, 15]
    x_aut = C[:, 16]
    x_pha = C[:, 17]
    data = np.column_stack((t, C))
    plt.figure(1)
    plt.subplot(2, 1, 1)
    plt.plot(t, s_nh4, 'k.-', label='s_nh4', linewidth=2.0)
    plt.plot(t, x_aut, 'b--', label='x_aut', linewidth=2.0)
    plt.plot(t, s_no3, 'y-', label='s_no3', linewidth=2.0)
    plt.plot(t, x_h, 'm--', label='x_h', linewidth=2.0)
    #plt.plot(t, x_s, 'c-', label='x_s', linewidth=2.0)
    plt.legend(loc='best')
    plt.xlabel('Time(d)')
    plt.yscale("log")
    plt.ylabel('Conc (g/m^3)')
    plt.show()
    if export_results =='yes':
        np.savetxt("output_ASM.csv", data, delimiter=",", header='t, s_n2 , s_nh4, s_no3,  x_meoh, s_po4, x_pp, x_mep, x_h, x_i, s_a, s_alk, x_pao, x_s, x_tss, s_f, s_i, x_aut, x_pha')


run_reactor('batch', 0, 2, 1000, C0, 'no') 

print(datetime.now() - start_time)  



#%% concrete required for tanks
# reactor dimensions and concrete needs
if calculate_materials == 'yes':
    depth_reactor = 5   # m
    freeboard_reactor = 0.3     # m 
    WD_ratio_reactor = 1.5
    t_reactor = 0.3     # m

    eff_depth_reactor = depth_reactor - freeboard_reactor   # filled depth
    width_reactor = depth_reactor * WD_ratio_reactor
    length_aerobic = V_aerobic / (width_reactor * eff_depth_reactor)
    V_aerobic_concrete_floor = t_reactor * (width_reactor*length_aerobic)
    if V_aerobic == 0:
        V_aerobic_concrete_wall = 0
    else:
        V_aerobic_concrete_wall = t_reactor * (2*width_reactor*depth_reactor + 2*width_reactor*length_aerobic)
    
    #add concrete calcs for anoxic tank
    
    
    #add concrete calcs for clarifer
    
    
    # total floor and wall concrete for activated sludge reactors
    V_concrete_reactor_floor = V_aerobic_concrete_floor
    V_concrete_reactor_wall =  V_aerobic_concrete_wall
    
    
    
#excavation & land (cost only)
# necessary land area based on size of reactors, clarifiers, thickeners, digesters, plus a little extra
# assume some fraction of reactors, clarifiers, etc. is underground and must be excavated




#%% aeration

    
    # oxygen transfer efficiency       
    # eq 5-70 (Tchobanoglous et al., 2014)
    DO_concentration = O2_conc
    DO_sat = 9.092  #mg/L at 20 C
    DO_sat_20 = 9.092  #mg/L at 20 C
    mid_depth_correction = 0.4
    diffuser_depth = 4.5
    beta = 0.95
    pressure_correction = 1
    SOTE = 0.38
    T = 20
    T_correction = 1.024
    T_standard = 20
    alpha = 0.6
    fouling_factor = 0.8
    O2_consumed = O2_conc * 0.5 * flowrate
    O2_MW = 32
    air_O2_fraction = 0.21
    air_MW = 28.97
    universal_gas_constant = 8.314
    compressor_eff = 0.80
    n = 0.283
    
    
    tau = DO_sat / DO_sat_20
    DO_sat_dif = DO_sat_20*(1 + mid_depth_correction*(diffuser_depth/10.33))
    OTE = SOTE * (((tau*beta*pressure_correction*DO_sat_dif - DO_concentration)/DO_sat_dif) 
                  * T_correction**(T - T_standard) * alpha * fouling_factor)
    
    # air pumping requirement
    O2_daily_consumption = O2_consumed * V_aerobic / 1000   # kg/d
    airflow = ((((O2_daily_consumption / OTE) / O2_MW) / air_O2_fraction) * air_MW)/24/60/60    # kg/s
    
    # power requirement of blower
    # eq 5-77 (Tchobanoglous et al., 2014)
    P_blower = ((airflow*universal_gas_constant*(T+273.15)) / (air_MW*n*compressor_eff)) * (((10.33+diffuser_depth)/10.33)**n - 1) # kW
    
    # daily energy needs of blower
    E_blower = P_blower * 24    # kWh/d
    
    # calculated as a point of interest: energy per kg of O2 delivered
    aeration_energy = E_blower / O2_daily_consumption # kWh/kg O2
    
    # blower costing/LCA impacts
    air_density = (101325 * air_MW) / (universal_gas_constant * (T + 273.15)) / 1000    # kg/m3
    airflow_scfm = airflow / air_density * 35.3147 * 60 # ft3/min
























