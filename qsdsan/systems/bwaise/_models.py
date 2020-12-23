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
from biosteam import plots, PowerUtility
from biosteam.evaluation import Model, Metric

from qsdsan import ImpactItem
from qsdsan.utils.loading import load_data, data_path
from qsdsan.utils.setter import Setter, DictAttrSetter
from bwaise import systems

np.random.seed(3221)
percentiles = [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1]


# %%

# =============================================================================
# Scenario A
# =============================================================================

sysA = systems.sysA
teaA = systems.teaA
lcaA = systems.lcaA

# Metrics
metrics = [
    Metric('Net cost without CAPEX', systems.get_AOC_cap, '$/cap/yr', 'TEA Results'),
    Metric('Net cost with CAPEX', systems.get_EAC_cap, '$/cap/yr', 'TEA Results'),
    Metric('Net emission', systems.get_GWP, 'kg CO2-eq/cap/yr', 'LCA Results'),
    Metric('Total N recovery', systems.get_N_recovery, '', 'N recovery'),
    Metric('Liquid N recovery', systems.get_liq_N_recovery, '', 'N recovery'),
    Metric('Solid N recovery', systems.get_sol_N_recovery, '', 'N recovery'),
    Metric('Total P recovery', systems.get_P_recovery, '', 'P recovery'),
    Metric('Liquid P recovery', systems.get_liq_P_recovery, '', 'P recovery'),
    Metric('Solid P recovery', systems.get_sol_P_recovery, '', 'P recovery'),
    Metric('Total K recovery', systems.get_K_recovery, '', 'K recovery'),
    Metric('Liquid K recovery', systems.get_liq_K_recovery, '', 'K recovery'),
    Metric('Solid K recovery', systems.get_sol_K_recovery, '', 'K recovery'),
    Metric('Total COD recovery', systems.get_COD_recovery, '', 'COD recovery'),
    Metric('Liquid COD recovery', systems.get_liq_COD_recovery, '', 'COD recovery'),
    Metric('Solid COD recovery', systems.get_sol_COD_recovery, '', 'COD recovery')
    ]

modelA = Model(sysA, metrics)
metric_dict = {i.name: i for i in modelA.metrics}

# Add parameters
param = modelA.parameter


# Assumptions
b = systems.tau_deg
D = shape.Uniform(lower=1, upper=3)
@param(name='Full degradation time', element=systems.A2, kind='coupled', units='yr',
        baseline=b, distribution=D)
def set_tau_deg(i):
    systems.tau_deg = i

b = systems.log_deg
D = shape.Uniform(lower=2, upper=4)
@param(name='Log degradation', element=systems.A2, kind='coupled', units='',
        baseline=b, distribution=D)
def set_log_deg(i):
    systems.log = i

b = systems.get_exchange_rate()
D = shape.Triangle(lower=3600, midpoint=b, upper=3900)
@param(name='Exchange rate', element='TEA', kind='isolated', units='UGX/USD',
        baseline=b, distribution=D)
def set_exchange_rate(i):
    systems.exchange_rate = i
    systems.A3.fee=systems.truck_cost['TankerTruck1']

b = teaA.discount_rate
D = shape.Uniform(lower=0.03, upper=0.06)
@param(name='Discount rate', element='TEA', kind='isolated', units='fraction',
        baseline=b, distribution=D)
def set_discount_rate(i):
    teaA.discount_rate = i

b = systems.CH4_item.CFs['GlobalWarming']
D = shape.Uniform(lower=28, upper=34)
@param(name='CH4 CF', element=systems.fugitive_CH4, kind='isolated',
        units='kg CO2-eq/kg CH4', baseline=b, distribution=D)
def set_CH4_CF(i):
    systems.CH4_item.CFs['GlobalWarming'] = i

b = systems.N2O_item.CFs['GlobalWarming']
D = shape.Uniform(lower=265, upper=298)
@param(name='N2O CF', element=systems.fugitive_N2O, kind='isolated',
        units='kg CO2-eq/kg N2O', baseline=b, distribution=D)
def set_N2O_CF(i):
    systems.N2O_item.CFs['GlobalWarming'] = i

b = PowerUtility.price
D = shape.Triangle(lower=0.08, midpoint=b, upper=0.21)
@param(name='Electricity price', element='TEA', kind='isolated',
        units='$/kWh', baseline=b, distribution=D)
def set_electricity_price(i):
    PowerUtility.price = i

b = systems.e_item.CFs['GlobalWarming']
D = shape.Uniform(lower=0.106, upper=0.121)
@param(name='Electricity CF', element='LCA', kind='isolated',
        units='kg CO2-eq/kWh', baseline=b, distribution=D)
def set_electricity_item(i):
    systems.e_item.CFs['GlobalWarming'] = i

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

# A2 = systems.A2
# b = A2.P_leaching
# D = shape.Uniform(lower=0, upper=0.37)
# @param(name='P leaching', element=systems.A2, kind='coupled', units='',
#         baseline=b, distribution=D)
# def set_P_leaching(i):
#     A2.P_leaching = i




#!!! Need to add other ones if to consider other decision variables considered
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


# Uncertainties with Trucking were set through ImpactItem

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
    modelA.parameter(name=para,
                     setter=DictAttrSetter(item, 'CFs', 'GlobalWarming'),
                     element='LCA', kind='isolated',
                     units=f'kg CO2-eq/{item.functional_unit}',
                     baseline=b, distribution=D)





# %%

# =============================================================================
# Run simulation and generate plots
# =============================================================================

N_sample = 100
samples = modelA.sample(N_sample, 'L')
modelA.load_samples(samples)
modelA.evaluate()

# Data organization
index_p = len(modelA.get_parameters())
p_df = modelA.table.iloc[:, :index_p].copy()
results = modelA.table.iloc[:, index_p:].copy()
results_percentiles = results.quantile(q=percentiles)

# Spearman's rank correlation
spearman_results = modelA.spearman(modelA.get_parameters(), metrics)
spearman_results.columns = pd.Index([i.name_with_units for i in metrics])

with pd.ExcelWriter('results/analysis_sysA.xlsx') as writer:
    p_df.to_excel(writer, sheet_name='Parameters')
    results.to_excel(writer, sheet_name='Uncertainty results')
    results_percentiles.to_excel(writer, sheet_name='Percentiles')
    # spearman_results.to_excel(writer, sheet_name='Spearman')
    modelA.table.to_excel(writer, sheet_name='Raw data')


writer.save()




from biosteam.utils import colors
light_color = colors.brown_tint.RGBn
dark_color = colors.brown_shade.RGBn


def plot_series_bp(title, df, start, end, light_color, dark_color, labels, ylim):
    fig, ax1 = plt.subplots(figsize=(10, 6))
    fig.canvas.set_window_title('title')
    fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
    for i in range(end-start):
        data = df.iloc[:, start+i:start+i+1]
        data = data.transpose()
        ax1.boxplot(x=data, positions=(i,), patch_artist=True,
                    widths=0.8, whis=[5, 95],
                    boxprops={'facecolor':light_color,
                              'edgecolor':dark_color},
                    medianprops={'color':dark_color,
                                  'linewidth':1.5},
                    flierprops = {'marker':'D',
                                  'markerfacecolor': light_color,
                                  'markeredgecolor': dark_color,
                                  'markersize':6})
    
    ax1.set_xticklabels(labels, rotation=45, fontsize=12)
    ax1.set_ylim(ylim)
    return fig


plot_series_bp('Costs and Emissions',
               results, 0, 3, colors.brown_tint.RGBn, colors.brown_shade.RGBn,
               ('Cost w CAPEX', 'Cost w/o CAPEX', 'GWP'), (0, 100))


plot_series_bp('Nitrogen Recovery',
               results, 3, 6, colors.yellow_tint.RGBn, colors.yellow_shade.RGBn,
               ('Total N recovery', 'Liq N recovery', 'Sol N recovery'), (0, 1))

plot_series_bp('Phosphorus Recovery',
               results, 6, 9, colors.purple_tint.RGBn, colors.purple_shade.RGBn,
               ('Total P recovery', 'Liq P recovery', 'Sol P recovery'), (0, 1))


plot_series_bp('Potassium Recovery',
               results, 3, 6, colors.CABBI_blue_light.RGBn, colors.CABBI_teal_green.RGBn,
               ('Total K recovery', 'Liq K recovery', 'Sol K recovery'), (0, 1))


























