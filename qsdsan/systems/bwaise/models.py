#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
Copyright (C) 2020, Quantitative Sustainable Design Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

# %%

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from chaospy import distributions as shape
from biosteam import PowerUtility
from biosteam.evaluation import Model, Metric
from qsdsan import currency, ImpactItem
from qsdsan.utils.loading import load_data, data_path
from qsdsan.utils.setters import Setter, DictAttrSetter
from qsdsan.systems import bwaise as bw


# %%

# =============================================================================
# Functions for batch-making metrics and -setting parameters
# =============================================================================

systems = bw.systems
sys_dct = systems.sys_dct
price_dct = systems.price_dct
GWP_dct = systems.GWP_dct
GWP = systems.GWP
get_summarizing_fuctions = systems.get_summarizing_fuctions
func = get_summarizing_fuctions()

def make_metrics(system):
    sys_ID = system.ID
    tea = sys_dct['TEA'][sys_ID]
    lca = sys_dct['LCA'][sys_ID]
    ppl = sys_dct['ppl'][sys_ID]
    unit = f'{currency}/cap/yr'
    cat = 'TEA results'
    metrics = [
        Metric('Annual cost', lambda: func['get_annual_cost'](tea, ppl), unit, cat),
        Metric('Annual CAPEX', lambda: func['get_annual_CAPEX'](tea, ppl), unit, cat),
        Metric('Annual OPEX', lambda: func['get_annual_OPEX'](tea, ppl), unit, cat),
        ]
    unit = f'{GWP.unit}/cap/yr'
    cat = 'LCA results'
    metrics.extend([
        Metric('Net emission', lambda: func['get_annual_GWP'](lca, ppl), unit, cat),
        Metric('Construction', lambda: func['get_constr_GWP'](lca, ppl), unit, cat),
        Metric('Transportation', lambda: func['get_trans_GWP'](lca, ppl), unit, cat),
        Metric('Direct emission', lambda: func['get_direct_emission_GWP'](lca, ppl), unit, cat),
        Metric('Offset', lambda: func['get_offset_GWP'](lca, ppl), unit, cat),
        Metric('Other', lambda: func['get_other_GWP'](lca, ppl), unit, cat),
        ])
    for i in ('COD', 'N', 'P', 'K'):
        cat = f'{i} recovery'
        metrics.extend([
            Metric(f'Liquid {i}', lambda: func[f'get_liq_{i}_recovery'](system, i), '', cat),
            Metric(f'Solid {i}', lambda: func[f'get_sol_{i}_recovery'](system, i), '', cat),
            Metric(f'Gas {i}', lambda: func[f'get_gas_{i}_recovery'](system, i), '', cat),
            Metric(f'Total {i}', lambda: func[f'get_tot_{i}_recovery'](system, i), '', cat)
            ])
    return metrics


def batch_setting_unit_para(df, model, unit, exclude=()):
    for para in df.index:
        if para in exclude: continue
        b = getattr(unit, para)
        lower = float(df.loc[para]['low'])
        upper = float(df.loc[para]['high'])
        dist = df.loc[para]['distribution']
        if dist == 'uniform':
            D = shape.Uniform(lower=lower, upper=upper)
        elif dist == 'triangular':
            D = shape.Triangle(lower=lower, midpoint=b, upper=upper)
        elif dist == 'constant': continue
        else:
            raise ValueError(f'Distribution {dist} not recognized for unit {unit}.')     
        model.parameter(name=para,
                        setter=Setter(unit, para),
                        element=unit, kind='coupled', units=df.loc[para]['unit'],
                        baseline=b, distribution=D)


# %%

# =============================================================================
# Scenario A (sysA)
# =============================================================================

sysA = systems.sysA
sysA.simulate()
metricsA = make_metrics(sysA)
modelA = Model(sysA, metricsA)
paramA = modelA.parameter

# Assumptions
A2 = systems.A2
b = systems.tau_deg
D = shape.Uniform(lower=1, upper=3)
@paramA(name='Full degradation time', element=A2, kind='coupled', units='yr',
        baseline=b, distribution=D)
def set_tau_deg(i):
    systems.tau_deg = i

b = systems.log_deg
D = shape.Uniform(lower=2, upper=4)
@paramA(name='Log degradation', element=A2, kind='coupled', units='',
        baseline=b, distribution=D)
def set_log_deg(i):
    systems.log = i

b = systems.get_exchange_rate()
D = shape.Triangle(lower=3600, midpoint=b, upper=3900)
@paramA(name='Exchange rate', element='TEA', kind='isolated', units='UGX/USD',
        baseline=b, distribution=D)
def set_exchange_rate(i):
    systems.exchange_rate = i

b = systems.get_discount_rate()
D = shape.Uniform(lower=0.03, upper=0.06)
@paramA(name='Discount rate', element='TEA', kind='isolated', units='fraction',
        baseline=b, distribution=D)
def set_discount_rate(i):
    systems.discount_rate = i

b = GWP_dct['CH4']
D = shape.Uniform(lower=28, upper=34)
@paramA(name='CH4 CF', element='LCA', kind='isolated', units='kg CO2-eq/kg CH4',
        baseline=b, distribution=D)
def set_CH4_CF(i):
    GWP_dct['CH4'] = i

b = GWP_dct['N2O']
D = shape.Uniform(lower=265, upper=298)
@paramA(name='N2O CF', element='LCA', kind='isolated', units='kg CO2-eq/kg N2O',
        baseline=b, distribution=D)
def set_N2O_CF(i):
    GWP_dct['N2O'] = i

b = price_dct['Electricity']
D = shape.Triangle(lower=0.08, midpoint=b, upper=0.21)
@paramA(name='Electricity price', element='TEA', kind='isolated',
        units='$/kWh', baseline=b, distribution=D)
def set_electricity_price(i):
    PowerUtility.price = i

b = GWP_dct['Electricity']
D = shape.Uniform(lower=0.106, upper=0.121)
@paramA(name='Electricity CF', element='LCA', kind='isolated',
        units='kg CO2-eq/kWh', baseline=b, distribution=D)
def set_electricity_item(i):
    GWP_dct['Electricity'] = i


# Excretion
A1 = systems.A1
su_data_path = data_path + 'sanunit_data/'
path = su_data_path + '_excretion.csv'
data = load_data(path)
batch_setting_unit_para(data, modelA, A1)
    
# PitLatrine
A2 = systems.A2
path = su_data_path + '_toilet.csv'
toilet_data = load_data(path)
path = su_data_path + '_pit_latrine.csv'
pit_latrine_data = load_data(path)
data = pd.concat((toilet_data, pit_latrine_data))
batch_setting_unit_para(data, modelA, A2, exclude=('MCF_decay', 'N2O_EF_decay'))

#!!! Need to add other ones if to consider other decision variables considered
b = A2.MCF_decay
D = shape.Triangle(lower=0.4, midpoint=b, upper=0.6)
@paramA(name='MCF_decay', element=A2, kind='coupled',
        units='fraction of anaerobic conversion of degraded COD',
        baseline=b, distribution=D)
def set_MCF_decay(i):
    A2._MCF_decay['communal_above_water'] = i
    
b = A2.N2O_EF_decay
D = shape.Triangle(lower=0, midpoint=b, upper=0.001)
@paramA(name='N2O_EF_decay', element=A2, kind='coupled',
        units='fraction of N emitted as N2O',
        baseline=b, distribution=D)
def set_N2O_EF_decay(i):
    A2._N2O_EF_decay['communal_above_water'] = i


# Uncertainties with Trucking were set through ImpactItem

# SedimentationTank
A5 = systems.A5
path = su_data_path + '_sedimentation_tank.csv'
data = load_data(path)
batch_setting_unit_para(data, modelA, A5)

# AnaerobicLagoon
A6 = systems.A6
path = su_data_path + '_anaerobic_lagoon.csv'
data = load_data(path)
batch_setting_unit_para(data, modelA, A6)

# FacultativeLagoon
A7 = systems.A7
path = su_data_path + '_facultative_lagoon.csv'
data = load_data(path)
batch_setting_unit_para(data, modelA, A7)

# DryingBed
A8 = systems.A8
path = su_data_path + '_drying_bed.csv'
data = load_data(path)
batch_setting_unit_para(data, modelA, A8, exclude=('sol_frac', 'bed_H'))

b = A8.sol_frac
D = shape.Uniform(lower=0.3, upper=0.4)
@paramA(name='sol_frac', element=A8, kind='coupled', units='fraction',
        baseline=b, distribution=D)
def set_sol_frac(frac):
    A8._sol_frac['unplanted'] = frac

b = A8.bed_H['covered']
D = shape.Uniform(lower=0.45, upper=0.75)
@paramA(name='non_storage_bed_H', element=A8, kind='coupled', units='m',
        baseline=b, distribution=D)
def set_non_storage_bed_H(H):
    A8.bed_H['covered'] = H
    A8.bed_H['uncovered'] = H

b = A8.bed_H['storage']
D = shape.Uniform(lower=1.2, upper=1.8)
@paramA(name='storage_bed_H', element=A8, kind='coupled', units='m',
        baseline=b, distribution=D)
def set_storage_bed_H(H):
    A8.bed_H['storage'] = H


# CropApplication
A9 = systems.A9
D = shape.Uniform(lower=0, upper=0.1)
@paramA(name='NH3 application loss', element=A9, kind='coupled',
        units='fraction of applied', baseline=0.02, distribution=D)
def set_app_NH3_loss(ratio):
    A9.loss_ratio['NH3'] = ratio
    
# Mg, Ca, C actually not affecting results
for key in A9.loss_ratio.keys():
    if key == 'NH3': continue
    D = shape.Uniform(lower=0, upper=0.05)
    def set_loss_ratio(ratio):
        A9.loss_ratio[key] = ratio
    modelA.parameter(name=key+' applicaiton loss',
                      setter=set_loss_ratio,
                      element=A9,
                      kind='coupled',
                      units='fraction of applied',
                      baseline=0.02, distribution=D)

# ImpactItems
path = data_path + '_impact_item.xlsx'
data = load_data(path, sheet='GWP')
items = ImpactItem._items

for para in data.index:
    item = ImpactItem._items[para]
    b = item.CFs['GlobalWarming']
    lower = float(data.loc[para]['low'])
    upper = float(data.loc[para]['high'])
    dist = data.loc[para]['distribution']
    if dist == 'uniform':
        D = shape.Uniform(lower=lower, upper=upper)
    elif dist == 'triangular':
        D = shape.Triangle(lower=lower, midpoint=b, upper=upper)
    elif dist == 'constant': continue
    else:
        raise ValueError(f'Distribution {dist} not recognized.')
    modelA.parameter(name=para,
                      setter=DictAttrSetter(item, 'CFs', 'GlobalWarming'),
                      element='LCA', kind='isolated',
                      units=f'kg CO2-eq/{item.functional_unit}',
                      baseline=b, distribution=D)


# %%

# =============================================================================
# Functions to run simulation and generate plots
# =============================================================================

try: result_dct
except:
    result_dct = {
        'sysA': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        'sysB': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        'sysC': dict.fromkeys(('parameters', 'data', 'percentiles', 'spearman')),
        }

def run_uncertainty(model, N_sample=1000, seed=None,
                    percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                    return_dct=False):
    if seed:
        np.random.seed(seed)
    if model is not modelA:
        return f'{model} has not been added.'
    samples = model.sample(N_sample, 'L')
    model.load_samples(samples)
    model.evaluate()
    dct = result_dct[model._system.ID]
    # Data organization
    index_p = len(model.get_parameters())
    dct['parameters'] = model.table.iloc[:, :index_p].copy()
    dct['data'] = model.table.iloc[:, index_p:].copy()
    dct['percentiles'] = dct['data'].quantile(q=percentiles)
    # Spearman's rank correlation
    spearman_metrics = [metricsA[i] for i in (0, 3, 12, 16, 20, 24)]
    spearman_results = model.spearman(modelA.get_parameters(), spearman_metrics)
    spearman_results.columns = pd.Index([i.name_with_units for i in spearman_metrics])
    dct['spearman'] = spearman_results
    if return_dct:
        return dct

def save_uncertainty_results(model, path=None):
    if not path:
        import os
        path = os.path.dirname(os.path.realpath(__file__))
        path += f'/results/model{model._system.ID[-1]}.xlsx'
        del os
    dct = result_dct[model._system.ID]
    try: len(dct['parameters'])
    except: raise ValueError('No cached result, run model first.')
    with pd.ExcelWriter(path) as writer:
        dct['parameters'].to_excel(writer, sheet_name='Parameters')
        dct['data'].to_excel(writer, sheet_name='Uncertainty results')
        dct['percentiles'].to_excel(writer, sheet_name='Percentiles')
        dct['spearman'].to_excel(writer, sheet_name='Spearman')
        model.table.to_excel(writer, sheet_name='Raw data')


# %%

# =============================================================================
# Functions to make quick box plots
# =============================================================================

from biosteam.utils import colors
light_color = colors.brown_tint.RGBn
dark_color = colors.brown_shade.RGBn

def plot_series_bp(title, dfs, light_color, dark_color, xlabels, ylabel, ylim):
    fig, ax1 = plt.subplots(figsize=(len(dfs), 6))
    fig.canvas.set_window_title(title)
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    for i, df in enumerate(dfs):      
        ax1.boxplot(x=df, positions=(i,), patch_artist=True,
                    widths=0.8, whis=[5, 95],
                    boxprops={'facecolor':light_color,
                              'edgecolor':dark_color},
                    medianprops={'color':dark_color,
                                  'linewidth':1.5},
                    flierprops = {'marker':'D',
                                  'markerfacecolor': light_color,
                                  'markeredgecolor': dark_color,
                                  'markersize':6})
    
    ax1.set_xticklabels(xlabels, rotation=45, fontsize=12)
    ax1.set_ylabel(ylabel, fontsize=12)
    ax1.set_ylim(ylim)
    return fig


def plot_cost_emission(model):
    sys_ID = model._system.ID
    data = result_dct[sys_ID]['data']
    cost_df = data.iloc[:, 0]
    emission_df = data.iloc[:, 3:4]
    dfs = (cost_df, emission_df)
    plot_series_bp(f'{sys_ID} Uncertainty Results',
                   dfs, colors.brown_tint.RGBn, colors.brown_shade.RGBn,
                   ('Net cost', 'Net GWP'),
                   f'Costs and Emissions [{currency} or {GWP.unit}/cap/yr]',
                   (0, 100))

def plot_recovery(model, resource):
    sys_ID = model._system.ID
    data = result_dct[sys_ID]['data']
    resources = ('COD', 'N', 'P', 'K')
    try:
        index = (resources.index(resource)+2)*4+1
    except:
        return f'resource can only be "COD", "N", "P", or "K", not "{resource}".'
    dfs = []
    for i in range(4):
        dfs.append(data.iloc[:, index+i:index+i+1])
    plot_series_bp(f'{sys_ID} Uncertainty Results',
                   dfs, colors.brown_tint.RGBn, colors.brown_shade.RGBn,
                   ('Liquid', 'Solid', 'Gas', 'Total'),
                   f'{resource} Recovery [%]',
                   (0, 1))




























