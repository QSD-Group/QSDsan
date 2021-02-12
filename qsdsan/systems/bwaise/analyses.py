#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt
for license details.
'''

# %%

import pandas as pd
import numpy as np
from qsdsan.stats import define_inputs, generate_samples, \
    get_correlation, morris_analysis, plot_morris_results, sobol_analysis
from qsdsan.systems import bwaise as bw

modelA = bw.models.modelA



imp_metrics = [i for i in modelA.metrics if 'Net' in i.name or 'Total' in i.name]


# %%

# =============================================================================
# Pearson and Spearman
# =============================================================================

# # This is the new, recommended method, but not seems to be widely used
# from numpy.random import MT19937, RandomState, SeedSequence
# rs = RandomState(MT19937(SeedSequence(3221)))
# rs.seed(3221)
# rs.random.sample(5)

#!!! PAUSED, USE FUNCTIONS TO FILTER PARAMETERS TO ONLY 0.2 SPEARMAN

np.random.seed(3221)
# For Monte Carlo uncertainty analysis
mc_samples = modelA.sample(N=100, rule='L')
modelA.load_samples(mc_samples)
modelA.evaluate()

pearson_r, pearson_p = get_correlation(modelA, kind='Pearson',
                                       nan_policy='raise')
spearman_rho, spearman_p = get_correlation(modelA, kind='Spearman',
                                           nan_policy='raise')
cost_r = pearson_r[('TEA results', 'Annual cost [USD/cap/yr]')]
emission_r = pearson_r[('LCA results', 'Net emission [kg CO2-eq/cap/yr]')]
cost_rho = spearman_rho[( 'TEA results', 'Annual cost [USD/cap/yr]')]
emission_rho = spearman_rho[('LCA results', 'Net emission [kg CO2-eq/cap/yr]')]

def filter_df(dfs, threshhold):
    new_dfs = []
    for df in dfs:
        new_df = pd.concat((df[df>threshhold], df[df<-threshhold]))
        new_dfs.append(new_df)
    return new_dfs

cost_r, emission_r, cost_rho, emission_rho = \
    filter_df((cost_r, emission_r, cost_rho, emission_rho), 0.2)


# %%

# =============================================================================
# Morris one at a time
# =============================================================================

inputs = define_inputs(modelA)
morris_samples = generate_samples(inputs, kind='Morris', N=100, seed=3221)
modelA.load_samples(morris_samples)

modelA.evaluate()

morris_dct = morris_analysis(modelA, morris_samples, inputs, 
                             nan_policy='fill_mean', seed=3221,
                             print_to_console=True, file='Morris.xlsx')
morris_plot = plot_morris_results(morris_dct, metric=modelA.metrics[0])


# %%

# =============================================================================
# Sobol
# =============================================================================

saltelli_samples = generate_samples(inputs, kind='Saltelli', N=100,
                                    calc_second_order=True)
modelA.load_samples(saltelli_samples)

modelA.evaluate()

sobol_dct = sobol_analysis(modelA, saltelli_samples, inputs,
                           calc_second_order=True, conf_level=0.95,
                           print_to_console=True, nan_policy='fill_mean',
                           file='Sobol.xlsx', seed=3221)























