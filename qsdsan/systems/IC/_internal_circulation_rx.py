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
import matplotlib.pyplot as plt
import biosteam as bst
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
    
    b_rx['OLR'] = b_rx['Qi']*all_rx['SS']/b_rx['V']    
    t_rx['OLR'] = t_rx['Qi']*b_rx['SS']/t_rx['V']
    
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
    assert 0<=t_rx['OLR']<=b_rx['OLR']
    assert 0<=t_rx['SS']<=b_rx['SS']<=all_rx['SS']
    assert 0<=t_rx['X']<=b_rx['X']
    assert min_COD_rm<=b_rx['COD_rm']<=all_rx['COD_rm']<= 1


def get_cod_rm(mu_max=0.01, Y=0.07, F_xb=0.0032, F_xt=0.0281,
               v=3, # v cannot be solved based on other parameters
               OLR_all=None, V_ratio=None, waste_ratio=None):
    
    # Overall
    all_rx['OLR'] = OLR_all
    all_rx['v'] = v
    all_rx['Qw'] = waste_ratio*all_rx['Qi']
    all_rx['Qe'] = all_rx['Qi'] - all_rx['Qw']
    all_rx['V'] = all_rx['Qi']*all_rx['SS']/OLR_all
    
    # Bottom rx
    b_rx['V'] = all_rx['V']/(1+1/V_ratio)
    # Biomass mass balance of the bottom rx
    b_rx['X'] = (all_rx['Qi']*all_rx['X'])/(F_xb*all_rx['Qe']+all_rx['Qw']-mu_max*b_rx['V'])
    b_rx['rX'] = mu_max * b_rx['X']
    b_rx['rSS'] = - b_rx['rX']/Y
    # Substrate mass balance of the bottom rx
    b_rx['SS'] = all_rx['SS'] + b_rx['V']*b_rx['rSS']/all_rx['Qi']
    
    # Top rx
    t_rx['V'] = b_rx['V']/V_ratio
    # Biomass mass balance of the top rx
    t_rx['X'] = (F_xb*all_rx['Qe']*b_rx['X'])/(F_xt*all_rx['Qe']-mu_max*t_rx['V'])
    t_rx['rX'] = mu_max * t_rx['X']
    t_rx['rSS'] = - t_rx['rX']/Y
    # Substrate mass balance of the top rx
    t_rx['SS'] = b_rx['SS'] + t_rx['V']*t_rx['rSS']/all_rx['Qe']
    
    all_rx['COD_rm'] = 1 - (t_rx['SS']*all_rx['Qe'])/(all_rx['SS']*all_rx['Qi'])
    return all_rx['COD_rm']



def rm_at_OLR(OLR_all, v, V_ratio, waste_ratio, COD_rm):
    rm = get_cod_rm(v=v, OLR_all=OLR_all, V_ratio=V_ratio, waste_ratio=waste_ratio)
    return rm-COD_rm

def rm_at_V_ratio(V_ratio, v, OLR_all, waste_ratio, COD_rm):
    rm = get_cod_rm(v=v, OLR_all=OLR_all, V_ratio=V_ratio, waste_ratio=waste_ratio)
    return rm-COD_rm

def rm_at_waste_ratio(waste_ratio, v, OLR_all, V_ratio, COD_rm):
    rm = get_cod_rm(v=v, OLR_all=OLR_all, V_ratio=V_ratio, waste_ratio=waste_ratio)
    return rm-COD_rm


def solve_param_from_rm(param='OLR_all', bounds=(), COD_rm=0.9, v=3,
                        OLR_all=None, V_ratio=None, waste_ratio=None):
    if not bounds:
        lb = ub = 0.
    else:
        lb, ub  = bounds
    if param == 'OLR_all':
        f = rm_at_OLR
        args = (v, V_ratio, waste_ratio, COD_rm)
        if OLR_all:
            warn(f'Solving for OLR_all, the given value {OLR_all} is ignored.')    
        if lb == 0:
            lb = 1e-4
        if ub == 0:
            ub = 35/24
    elif param == 'V_ratio':
        f = rm_at_V_ratio
        args = (v, OLR_all, waste_ratio, COD_rm)
        if V_ratio:
            warn(f'Solving for V_ratio, the given value {V_ratio} is ignored.')
        if lb == 0:
            lb = 1e-4
        if ub == 0:
            ub = 100
    elif param == 'waste_ratio':
        f = rm_at_waste_ratio
        args = (v, OLR_all, V_ratio, COD_rm)
        if waste_ratio:
            warn(f'Solving for waste_ratio, the given value {waste_ratio} is ignored.')    
        if ub == 0 or ub >= 1:
            ub = 1-1e-4
    else:
        raise ValueError("param can only be in ('OLR_all', 'V_ratio', "
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



def get_design():
    update_design()
    df = pd.DataFrame(data={'Overall': (i for i in all_rx.values()),
                            'Bottom': (i for i in b_rx.values()),
                            'Top': (i for i in t_rx.values()),},
                      index=indices)
    return df

# %%

def get_demo_results():
    # Can get COD removal with given OLR_all, V_ratio, v, and waste_ratio
    COD_rm = get_cod_rm(mu_max=0.01, Y=0.07, F_xb=0.0032, F_xt=0.0281,
                        OLR_all=1.084, V_ratio=1.5, v=3, waste_ratio=0.05)
    print('\n------------Get COD removal------------')
    print(f'COD removal is {COD_rm:.1%}.')
    
    # Or can solve for any of OLR_all, V_ratio, and waste_ratio with target COD removal
    print('\n------------Solve for overall organic loading rate ------------')
    nOLR_all = solve_param_from_rm(param='OLR_all', bounds=(0, 35/24), OLR_all=None,
                                   v=3, V_ratio=1.5, waste_ratio=0.05, COD_rm=0.918)
    print(f'\nOLR_all is {nOLR_all:.3f}.')
    design = get_design()
    print(design)
    
    print('\n------------Solve for bottom-to-top rx volume ratio ------------')
    V_ratio = solve_param_from_rm(param='V_ratio', bounds=(0, 10), V_ratio=None,
                                  v=3, OLR_all=1.084, waste_ratio=0.05, COD_rm=0.918)
    print(f'\nV_ratio is {V_ratio:.1f}.')
    design = get_design()
    print(design)
    
    print('\n------------Solve for waste sludge ratio ------------')
    waste_ratio = solve_param_from_rm(param='waste_ratio', bounds=(0, 1), waste_ratio=None,
                                      OLR_all=1.084, V_ratio=1.5, v=3, COD_rm=0.918)
    print(f'\nwaste_ratio is {waste_ratio:.1%}.')
    design = get_design()
    print(design)


get_demo_results()


# %%

OLRs = np.arange(start=0.5, stop=35/24, step=0.01)
waste_ratios = np.arange(start=0, stop=1, step=0.01)
V_ratios = np.arange(start=0.5, stop=3.1, step=0.5)

results = {
    'OLR': [],
    'waste_ratio': [],
    'COD_rm': [],
    'Vtot': []
    }

def evalute_IC_design(OLRs, waste_ratios, V_ratios, excel_name=None,
                      dropna=True):
    dfs = []
    for r in V_ratios:
        for val in results.values():
            val.clear()
        for i in OLRs:
            for j in waste_ratios:
                COD_rm = get_cod_rm(OLR_all=i, waste_ratio=j, V_ratio=r)
                results['OLR'].append(i)
                results['waste_ratio'].append(j)
                results['Vtot'].append(all_rx['V'])
                try:
                    check_constraints()
                    results['COD_rm'].append(COD_rm)
                except:
                    results['COD_rm'].append(np.nan)
        df = pd.DataFrame.from_dict(results)
        dfs.append(df)
    if excel_name:
        writer = pd.ExcelWriter(path='IC_design.xlsx')
        for df in dfs:
            df.dropna().to_excel(writer, sheet_name=f'V_ratio={r}')
        writer.save()
    return dfs

dfs = evalute_IC_design(OLRs, waste_ratios, V_ratios, excel_name='COD_rm.xlsx')


def plot_COD_rm(data, x_axis, y_axis, nplts, ncols, fig_name=None):
    X, Y = np.meshgrid(x_axis, y_axis)
    nrows = int(np.ceil((nplts/ncols)))
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, sharex=True, sharey=True,
                             constrained_layout=True)
    num = 0
    for m in range(nrows):
        for n in range(ncols):
            df = data[num]['COD_rm'].to_numpy().reshape(X.shape) 
            cs = axes[m, n].contourf(X, Y, df, levels=np.arange(0, 1.1, 0.1))
            if m == nrows-1:
                axes[m, n].set_xlabel('Waste Ratio [%]')
            if n == 0:            
                axes[m, n].set_ylabel('OLR [kg/m3/hr]')
            axes[m, n].set_title(f'V_ratio={V_ratios[num]}')
            num += 1
    fig.colorbar(cs)
    if fig_name:
        fig.savefig(fig_name, dpi=300)
    return fig, axes

fig, axes = plot_COD_rm(dfs, waste_ratios, OLRs, nplts=len(V_ratios), ncols=3,
                        fig_name='COD_rm.png')

# %%

# =============================================================================
# Methane calculation
# =============================================================================

from biorefineries import cornstover as cs
cs.cornstover_sys.simulate()
chems = cs.chemicals

'''
Steps:
    Set all gases and inorganics to 0
    Set tiers 1 and 2 (make a list/dict of them)
    In uncertainty analysis, set %BD to 70%, 80%, 90% or similar for theoretical ones


Formula to calculate theoretical methane production being:
    CvHwOxNySz + (v-w/4-x/2+3y/4+z/2)H2O -> \
        (v/2+w/8-x/4-3y/8-z/4)CH4 + (v/2-w/8+x/4+3y/8+z/4)CO2 + yNH3 + zH2S

'''

#TODO: Make these into Reaction objects, similar to combustion
def get_yield_dct(chemicals):
    '''
    Get dicts of theoretical gas product yields mol as product/mol chemical.
    '''
    yield_dct = {
        'CH4': {},
        'CO2': {},
        'NH3': {},
        'H2S': {}
        }
    for chem in chemicals:
        if not chem.formula:
            for dct in yield_dct.values():
               dct[chem.ID] = 0
        else:
            d = dict.fromkeys(('C', 'H', 'O', 'N', 'S'), 0)
            d.update(chem.atoms)
            if d['C'] == 0:
                yield_dct['CH4'][chem.ID] = yield_dct['CO2'][chem.ID] = 0
            else:
                yield_dct['CH4'][chem.ID] = \
                    d['C']/2+d['H']/8-d['O']/4-d['N']*3/8-d['S']/4
                yield_dct['CO2'][chem.ID] = \
                    d['C']/2-d['H']/8+d['O']/4+d['N']*3/8+d['S']/4
                yield_dct['NH3'][chem.ID] = d['N']
                yield_dct['H2S'][chem.ID] = d['S']
    # CSL stream is modeled as 50% water, 25% protein, and 25% lactic acid
    for dct in yield_dct.values():
       dct['CH4'] = dct['Protein']/4 + dct['LacticAcid']/4
    return yield_dct

yield_dct = get_yield_dct(chems)

BD_dct = dict.fromkeys((i.ID for i in cs.chemicals), 1)

# Modify BMP and CO2 based on experimental data
#!!! Does the biodegradability apply for CO2 as well?
BD_dct['AceticAcid'] = 0.87
BD_dct['LacticAcid'] = 0.85
BD_dct['HMF'] = 0.85
BD_dct['Glucose'] = 0.87
BD_dct['Xylose'] = 0.8
BD_dct['Arabinose'] = 0.04
BD_dct['Lignin'] = BD_dct['SolubleLignin'] = 0
BD_dct['Tar'] = 0


def get_yield_array(BD_dct, yield_dct):
    BD_array = chems.kwarray(BD_dct)
    yield_array = {
        'CH4': chems.kwarray(yield_dct['CH4']) * BD_array,
        'CO2': chems.kwarray(yield_dct['CO2']) * BD_array,
        'NH3': chems.kwarray(yield_dct['NH3']) * BD_array,
        'H2S': chems.kwarray(yield_dct['H2S']) * BD_array
        }
    return yield_array

yield_array = get_yield_array(BD_dct, yield_dct)


AD_inf = cs.R601.ins[0]

get_CH4_mol = lambda stream: (stream.mol*yield_array['CH4']).sum()
# Calculated is 334 kmol/hr, cs.R601.outs[0].imol['CH4'] is around 328 kmol/hr
CH4_mol = get_CH4_mol(AD_inf)
# COD in kg/m3/hr
#!!! Is H2S counted in COD calculation?
get_COD = lambda stream: get_CH4_mol(stream)*chems.CH4.MW/stream.F_vol

def get_gas_volume(stream, yield_array):
    vol = 0
    for ID, array in yield_array.items():
        # AD temperature around 35°C
        V = getattr(chems, ID).V(35+298.15, 101325) # in m3/mol
        vol += (stream.mol*array).sum()*V*1e3 # stream.mol in kmol/hr
    return vol

# Calculated is 18385 m3/hr, cs.R601.outs[0].F_vol is around 18560 m3/hr
AD_gas_vol = get_gas_volume(AD_inf, yield_array)
































