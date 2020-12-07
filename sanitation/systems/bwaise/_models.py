#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Sanitation Explorer: Sustainable design of non-sewered sanitation technologies
Copyright (C) 2020, Sanitation Explorer Development Group

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the UIUC open-source license. Please refer to 
https://github.com/QSD-for-WaSH/sanitation/blob/master/LICENSE.txt
for license details.

Ref:
    [1] Trimmer et al., Navigating Multidimensional Social–Ecological System
        Trade-Offs across Sanitation Alternatives in an Urban Informal Settlement.
        Environ. Sci. Technol. 2020, 54 (19), 12641–12653.
        https://doi.org/10.1021/acs.est.0c03296.

'''

# %%

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from chaospy import distributions as shape
from biosteam import plots
from biosteam.evaluation import Model, Metric
from sanitation import ImpactItem
from sanitation.utils.loading import load_data, data_path
from bwaise import systems

np.random.seed(3221)


# %%

# =============================================================================
# Scenario A
# =============================================================================

sysA = systems.sysA
teaA = systems.teaA
lcaA = systems.lcaA

# Metrics
metrics = [
    Metric('Net cost without CAPEX', systems.get_AOC_cap, '$/cap/yr'),
    Metric('Net cost with CAPEX', systems.get_EAC_cap, '$/cap/yr'),
    Metric('Net emission', systems.get_GWP, 'kg CO2-eq/cap/yr'),
    Metric('Total N recovery', systems.get_N_recovery, ''),
    Metric('Liquid N recovery', systems.get_liq_N_recovery, ''),
    Metric('Solid N recovery', systems.get_sol_N_recovery, ''),
    Metric('Total P recovery', systems.get_P_recovery, ''),
    Metric('Liquid P recovery', systems.get_liq_P_recovery, ''),
    Metric('Solid P recovery', systems.get_sol_P_recovery, ''),
    Metric('Total K recovery', systems.get_K_recovery, ''),
    Metric('Liquid K recovery', systems.get_liq_K_recovery, ''),
    Metric('Solid K recovery', systems.get_sol_K_recovery, ''),
    Metric('Total COD recovery', systems.get_COD_recovery, ''),
    Metric('Liquid COD recovery', systems.get_liq_COD_recovery, ''),
    Metric('Solid COD recovery', systems.get_sol_COD_recovery, '')
    ]

modelA = Model(sysA, metrics)
metric_dict = {i.name: i for i in modelA.metrics}

# Add parameters
param = modelA.parameter

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
                        setter=lambda i: setattr(unit, para, i),
                        element=unit, kind='coupled', units=df.loc[para]['unit'],
                        baseline=b, distribution=D)

# Excretion
A1 = systems.A1
path = data_path + 'unit_data/_excretion.csv'
data = load_data(path)
batch_setting_unit_para(data, modelA, A1)
    
# PitLatrine
A2 = systems.A2
path = data_path + 'unit_data/_toilet.csv'
toilet_data = load_data(path)
path = data_path + 'unit_data/_pit_latrine.csv'
pit_latrine_data = load_data(path)
data = pd.concat((toilet_data, pit_latrine_data))
batch_setting_unit_para(data, modelA, A2, exclude=('MCF_decay', 'N2O_EF_decay'))

#!!! Need to add other ones if other situations considered
b = A2.MCF_decay
D = shape.Triangle(lower=0.4, midpoint=b, upper=0.6)
@param(name='MCF_decay', element=A2, kind='coupled',
        units='fraction of anaerobic conversion of degraded COD',
        baseline=b, distribution=D)
def set_MCF_decay(i):
    A2._MCF_decay['communal_above_water'] = i
    
b = A2.N2O_EF_decay
D = shape.Triangle(lower=0, midpoint=b, upper=0.001)
@param(name='N2O_EF_decay', element=A2, kind='coupled',
        units='fraction of N emitted as N2O',
        baseline=b, distribution=D)
def set_N2O_EF_decay(i):
    A2._N2O_EF_decay['communal_above_water'] = i


# Trucking - #!!! set through ImpactItem

# SedimentationTank
A4 = systems.A4
path = data_path + 'unit_data/_sedimentation_tank.csv'
data = load_data(path)
batch_setting_unit_para(data, modelA, A4)

# AnaerobicLagoon
A5 = systems.A5
path = data_path + 'unit_data/_anaerobic_lagoon.csv'
data = load_data(path)
batch_setting_unit_para(data, modelA, A5)

# FacultativeLagoon
A6 = systems.A6
path = data_path + 'unit_data/_facultative_lagoon.csv'
data = load_data(path)
batch_setting_unit_para(data, modelA, A6)

# DryingBed
A7 = systems.A7
path = data_path + 'unit_data/_drying_bed.csv'
data = load_data(path)
batch_setting_unit_para(data, modelA, A7, exclude=('sol_frac', 'bed_H'))

b = A7.sol_frac
D = shape.Uniform(lower=0.3, upper=0.4)
@param(name='sol_frac', element=A7, kind='coupled', units='fraction',
        baseline=b, distribution=D)
def set_sol_frac(frac):
    A7._sol_frac['unplanted'] = frac

b = A7.bed_H['covered']
D = shape.Uniform(lower=0.45, upper=0.75)
@param(name='non_storage_bed_H', element=A7, kind='coupled', units='m',
        baseline=b, distribution=D)
def set_non_storage_bed_H(H):
    A7.bed_H['covered'] = H
    A7.bed_H['uncovered'] = H

b = A7.bed_H['storage']
D = shape.Uniform(lower=1.2, upper=1.8)
@param(name='storage_bed_H', element=A7, kind='coupled', units='m',
        baseline=b, distribution=D)
def set_storage_bed_H(H):
    A7.bed_H['storage'] = H


# CropApplication
A8 = systems.A8
D = shape.Uniform(lower=0, upper=0.1)
@param(name='NH3 application loss', element=A8, kind='coupled',
        units='fraction of applied', baseline=0.02, distribution=D)
def set_app_NH3_loss(ratio):
    A8.loss_ratio['NH3'] = ratio
    
# Mg, Ca, C actually not affecting results
for element in A8.loss_ratio.keys():
    if element == 'NH3': continue
    D = shape.Uniform(lower=0, upper=0.05)
    def set_loss_ratio(ratio):
        A8.loss_ratio[element] = ratio
    modelA.parameter(name=element+' applicaiton loss',
                     setter=set_loss_ratio,
                     element=A8, kind='coupled',
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
    def set_CF(CF):
        item.CFs['GlobalWarming'] = CF
    modelA.parameter(name=para,
                      setter=set_CF,
                      element='Biorefinery', kind='isolated',
                      units=f'kg CO2-eq/{item.functional_unit}',
                      baseline=b, distribution=D)



#!!! Still need to add a bunch of global parameters


# %%

# =============================================================================
# Run simulation and generate plots
# =============================================================================

N_sample = 100
samples = modelA.sample(N_sample, 'L')
modelA.load_samples(samples)
modelA.evaluate()

modelA.table.to_excel('results/analysis_sysA.xlsx')

def plot_montecarlo(metrics, N=100):
    """
    Create a matplotlib boxplot and return plot objects.
    
    Examples
    --------
    >>> # plot_montecarlo('Excess electricity', 100) -> creates plot
    {'whiskers': [<matplotlib.lines.Line2D at 0x27ac8469b88>,
      <matplotlib.lines.Line2D at 0x27ab70a9c08>],
     'caps': [<matplotlib.lines.Line2D at 0x27ac852ed88>,
      <matplotlib.lines.Line2D at 0x27ac8551f48>],
     'boxes': [<matplotlib.patches.PathPatch at 0x27ac84692c8>],
     'medians': [<matplotlib.lines.Line2D at 0x27ac855b2c8>],
     'fliers': [<matplotlib.lines.Line2D at 0x27ac8551fc8>],
     'means': []}
    
    """
    samples = modelA.sample(N, 'L')
    modelA.load_samples(samples)
    modelA.evaluate()
    for metric in metrics:
        df = modelA.table[[metric.index]]
        bx = plots.plot_montecarlo(df, transpose=True)
        plt.xticks([], [])
        plt.ylabel(metric.name_with_units)
    return bx


# def plot_montecarlo(metric, N=100):
#     """
#     Create a matplotlib boxplot and return plot objects.
    
#     Examples
#     --------
#     >>> # plot_montecarlo('Excess electricity', 100) -> creates plot
#     {'whiskers': [<matplotlib.lines.Line2D at 0x27ac8469b88>,
#       <matplotlib.lines.Line2D at 0x27ab70a9c08>],
#      'caps': [<matplotlib.lines.Line2D at 0x27ac852ed88>,
#       <matplotlib.lines.Line2D at 0x27ac8551f48>],
#      'boxes': [<matplotlib.patches.PathPatch at 0x27ac84692c8>],
#      'medians': [<matplotlib.lines.Line2D at 0x27ac855b2c8>],
#      'fliers': [<matplotlib.lines.Line2D at 0x27ac8551fc8>],
#      'means': []}
    
#     """
#     if isinstance(metric, str): metric = metric_dict[metric]
#     samples = modelA.sample(N, 'L')
#     modelA.load_samples(samples)
#     modelA.evaluate()
#     df = modelA.table[[metric.index]]
#     bx = plots.plot_montecarlo(df, transpose=True)
#     plt.xticks([], [])
#     plt.ylabel(metric.name_with_units)
#     return bx





















