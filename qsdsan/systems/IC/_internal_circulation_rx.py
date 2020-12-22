#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 10:24:11 2020

@author: yalinli_cabbi
"""



# %%

import numpy as np
import pandas as pd
import flexsolve as fs
from warnings import warn
from biosteam.exceptions import DesignError

param_info = (
    ['Qi', 'influent liquid flow', '[m3/hr]'],
    ['Qe', 'effluent liquid flow', '[m3/hr]'],
    ['Qw', 'waste sludge flow', '[m3/hr]'],
    ['Qg', 'gas flow', '[m3/hr]'],
    ['OLR', 'organic loading rate', '[kg/m3/hr]'],
    ['SS', 'soluble substrate conc.', '[kg/m3 as COD]'],
    ['rSS', 'SS reaction rate', '[kg/m3/hr as COD]'],
    ['X', 'biomass conc.', '[kg/m3]'],
    ['rX', 'biomass growth rate', '[kg/m3/hr]'],
    ['V', 'reactor volume', '[m3]'],
    ['A', 'reactor area', '[m2]'],
    ['H', 'reactor height', '[m]'],
    ['D', 'reactor diameter', '[m]'],
    ['HtoD', 'reactor aspect ratio', ''],
    ['v', 'upflow velocity', '[m/hr]'],
    ['HRT', 'hydraulic retention time', '[hr]'],
    ['SRT', 'solid retention time', '[hr]'],
    ['COD_rm', 'COD removal', '']
    )

indices = pd.MultiIndex.from_tuples(
    param_info, names=('Parameter', 'Description', 'Unit'))




all_rx = dict.fromkeys([i[0] for i in param_info], 0)
b_rx = all_rx.copy()
t_rx = all_rx.copy()


# Upstream inputs
all_rx['Qi'] = 325 # influent volumetric flow rate, [m3/hr]
all_rx['SS'] = 5040/1e3 # influent soluble substrate conc., [kg/m3 as COD]
all_rx['X'] = 230/1e3 # influent biomass conc., [kg/m3]


# Assumptions based on literature
mu_max = 0.01 # maximum specific growth rate, [/hr]
Y = 0.07 # biomass yield, [kg biomass/kg substrate as COD]
F_xb = 0.0032 # biomass transfer ratio from bottom to top, 0-1 (ideal to no retention)
F_xt = 0.0281 # biomass transfer ratio from the top rx to effluent
lit_assumptions = (mu_max, Y, F_xb, F_xt)

def update_design():
    all_rx['Qe'] = all_rx['Qi'] - all_rx['Qw']
    for i in ('Qi', 'Qw', 'Qe'):
        b_rx[i] = all_rx[i]
    t_rx['Qi'] = t_rx['Qe'] = all_rx['Qe']
    
    all_rx['SRT'] = (b_rx['X']*b_rx['V']+t_rx['X']*t_rx['V']) / \
        (b_rx['Qw']*b_rx['V']+F_xt*t_rx['X']*t_rx['Qe'])
    b_rx['SRT'] = b_rx['X']*b_rx['V']/(b_rx['X']*b_rx['Qw']+F_xb*b_rx['X']*b_rx['Qe'])
    t_rx['SRT'] = t_rx['X']*t_rx['V']/(t_rx['X']*t_rx['Qw']+F_xt*t_rx['X']*t_rx['Qe'])
    
    for i in (all_rx, b_rx, t_rx):
        if i is all_rx:
            i['A'] = i['Qe']/i['v'] #!!! Qe or Qi?
            i['D'] = (i['A']/(np.pi/4))**(1/2)
            i['rX'] = mu_max * i['X']
            i['rSS'] = - i['rX']/Y
        else:
            i['OLR'] = i['Qi']*i['SS']/i['V']
            i['A'] = all_rx['A']
            i['D'] = all_rx['D']
            i['v'] = all_rx['v']
            i['COD_rm'] = 1 - (i['SS']*i['Qe'])/(all_rx['SS']*all_rx['Qi'])
        i['H'] = i['V']/i['A']
        i['HtoD'] = i['H']/i['D']
        i['HRT'] = i['V']/i['Qi']
        



def check_constraints(max_OLR=None, min_COD_rm=None):
    update_design()
    max_OLR = max_OLR or 35/24
    min_COD_rm = min_COD_rm or 0
    assert 0<=all_rx['OLR']<=max_OLR
    assert 0<=t_rx['SS']<=b_rx['SS']<=all_rx['SS']
    assert 0<=t_rx['X']<=b_rx['X']
    assert min_COD_rm<=b_rx['COD_rm']<=all_rx['COD_rm']<= 1


def get_cod_rm(mu_max=0.01, Y=0.07, F_xb=0.0032, F_xt=0.0281,
               v=3, # v cannot be solved based on other parameters
               OLR_all=None, V_b_to_V_t=None, waste_ratio=None):
    
    # Overall
    all_rx['OLR'] = OLR_all
    all_rx['v'] = v
    all_rx['Qw'] = waste_ratio*all_rx['Qi']
    all_rx['Qe'] = all_rx['Qi'] - all_rx['Qw']
    all_rx['V'] = all_rx['Qi']*all_rx['SS']/OLR_all
    
    # Bottom rx
    b_rx['V'] = all_rx['V']/(1+1/V_b_to_V_t)
    # Biomass mass balance of the bottom rx
    b_rx['X'] = (all_rx['Qi']*all_rx['X'])/(F_xb*all_rx['Qe']+all_rx['Qw']-mu_max*b_rx['V'])
    b_rx['rX'] = mu_max * b_rx['X']
    b_rx['rSS'] = - b_rx['rX']/Y
    # Substrate mass balance of the bottom rx
    b_rx['SS'] = all_rx['SS'] + b_rx['V']*b_rx['rSS']/all_rx['Qi']
    
    # Top rx
    t_rx['V'] = b_rx['V']/V_b_to_V_t
    # Biomass mass balance of the top rx
    t_rx['X'] = (F_xb*all_rx['Qe']*b_rx['X'])/(F_xt*all_rx['Qe']-mu_max*t_rx['V'])
    t_rx['rX'] = mu_max * t_rx['X']
    t_rx['rSS'] = - t_rx['rX']/Y
    # Substrate mass balance of the top rx
    t_rx['SS'] = b_rx['SS'] + t_rx['V']*t_rx['rSS']/all_rx['Qe']
    
    all_rx['COD_rm'] = 1 - (t_rx['SS']*all_rx['Qe'])/(all_rx['SS']*all_rx['Qi'])
    
    return all_rx['COD_rm']



def rm_at_OLR(OLR_all, v, V_b_to_V_t, waste_ratio, COD_rm):
    rm = get_cod_rm(v=v, OLR_all=OLR_all, V_b_to_V_t=V_b_to_V_t, waste_ratio=waste_ratio)
    return rm-COD_rm

def rm_at_V_ratio(V_b_to_V_t, v, OLR_all, waste_ratio, COD_rm):
    rm = get_cod_rm(v=v, OLR_all=OLR_all, V_b_to_V_t=V_b_to_V_t, waste_ratio=waste_ratio)
    return rm-COD_rm

def rm_at_waste_ratio(waste_ratio, v, OLR_all, V_b_to_V_t, COD_rm):
    rm = get_cod_rm(v=v, OLR_all=OLR_all, V_b_to_V_t=V_b_to_V_t, waste_ratio=waste_ratio)
    return rm-COD_rm


def solve_param_from_rm(param='OLR_all', bounds=(), COD_rm=0.9, v=3,
                        OLR_all=None, V_b_to_V_t=None, waste_ratio=None):
    if not bounds:
        lb = ub = 0.
    else:
        lb, ub  = bounds
    if param == 'OLR_all':
        f = rm_at_OLR
        args = (v, V_b_to_V_t, waste_ratio, COD_rm)
        if OLR_all:
            warn(f'Solving for OLR_all, the given value {OLR_all} is ignored.')    
        if lb == 0:
            lb = 1e-4
        if ub == 0:
            ub = 35/24
    elif param == 'V_b_to_V_t':
        f = rm_at_V_ratio
        args = (v, OLR_all, waste_ratio, COD_rm)
        if V_b_to_V_t:
            warn(f'Solving for V_b_to_V_t, the given value {V_b_to_V_t} is ignored.')
        if lb == 0:
            lb = 1e-4
        if ub == 0:
            ub = 100
    elif param == 'waste_ratio':
        f = rm_at_waste_ratio
        args = (v, OLR_all, V_b_to_V_t, COD_rm)
        if waste_ratio:
            warn(f'Solving for waste_ratio, the given value {waste_ratio} is ignored.')    
        if ub == 0 or ub >= 1:
            ub = 1-1e-4
    else:
        raise ValueError("param can only be in ('OLR_all', 'V_b_to_V_t', "
                         f"'waste_ratio'), not {param}).")
    lb0, ub0 = lb, ub
    step = (ub-lb)/100
    # Boundary matters, so try to solve it both ways
    while lb <= ub:
        val = fs.IQ_interpolation(f=f, x0=lb, x1=ub,
                                  xtol=1e-4, ytol=1e-5, maxiter=50,
                                  args=args, checkbounds=False)
        diff = f(val, *args)
        try:
            assert abs(diff)<1e-4
            check_constraints()
            return val
        except:
            lb += step
            continue
    lb, ub = lb0, ub0
    while lb <= ub:
        val = fs.IQ_interpolation(f=f, x0=lb, x1=ub,
                                  xtol=1e-4, ytol=1e-5, maxiter=50,
                                  args=args, checkbounds=False)
        diff = f(val, *args)
        try:
            assert abs(diff)<1e-4
            check_constraints()
            return val
        except:
            ub -= step
            continue
    raise DesignError('Did not find a feasible design for the given specifications.')



def get_design(index):
    update_design()
    df = pd.DataFrame(data={'Overall': (i for i in all_rx.values()),
                            'Bottom': (i for i in b_rx.values()),
                            'Top': (i for i in t_rx.values()),},
                      index=index)
    return df


# Can get COD removal with given OLR_all, V_b_to_V_t, v, and waste_ratio
COD_rm = get_cod_rm(mu_max=0.01, Y=0.07, F_xb=0.0032, F_xt=0.0281,
                OLR_all=1.084, V_b_to_V_t=1.5, v=3, waste_ratio=0.05)
print('\n------------Get COD removal------------')
print(f'COD removal is {COD_rm:.1%}.')

# Or can solve for any of OLR_all, V_b_to_V_t, and waste_ratio with target COD removal
nOLR_all = solve_param_from_rm(param='OLR_all', bounds=(0, 35/24), OLR_all=None,
                               v=3, V_b_to_V_t=1, waste_ratio=0.05, COD_rm=0.92)
print('\n------------Solve for overall organic loading rate ------------')
print(f'\nOLR_all is {nOLR_all:.3f}.')
design = get_design(indices)
print(design)

V_b_to_V_t = solve_param_from_rm(param='V_b_to_V_t', bounds=(0, 100), V_b_to_V_t=None,
                                 v=3, OLR_all=1.084, waste_ratio=0.05, COD_rm=0.92)
print('\n------------Solve for bottom-to-top rx volume ratio ------------')
print(f'\nV_b_to_V_t is {V_b_to_V_t:.1f}.')
design = get_design(indices)
print(design)

waste_ratio = solve_param_from_rm(param='waste_ratio', bounds=(0, 1), waste_ratio=None,
                                  OLR_all=1.084, V_b_to_V_t=1.5, v=3, COD_rm=0.92)
print('\n------------Solve for waste sludge ratio ------------')
print(f'\nwaste_ratio is {waste_ratio:.1%}.')
design = get_design(indices)
print(design)

















