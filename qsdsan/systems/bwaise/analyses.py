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

import warnings
warnings.filterwarnings(action='ignore')
# # The warning is because of some outdated packages
# warnings.filterwarnings(action='ignore', message='`should_run_async`')

import pandas as pd
from qsdsan.stats import define_inputs, generate_samples, \
    get_correlation, morris_analysis, plot_morris_results, sobol_analysis
from qsdsan.utils.decorators import time_printer
from qsdsan.systems import bwaise as bw

m = bw.models
modelA = m.modelA

# Net cost, net GWP, and total COD/N/P/K recovery
key_metrics = [i for i in modelA.metrics if 'Net' in i.name or 'Total' in i.name]

get_param_dct = lambda model: {p.name_with_units:p for p in model.get_parameters()}

def filter_parameters(model, df, threshold):
    new_df = pd.concat((df[df>=threshold], df[df<=-threshold]))
    filtered = new_df.dropna(how='all')
    param_dct = get_param_dct(model)
    parameters = set(param_dct[i[1]] for i in filtered.index)
    return parameters

@time_printer
def evaluate(model, samples, print_time=False):
    model.load_samples(samples)
    model.evaluate()


# %%

'''
All analyses below are run for sysA/modelA, simply change "A" to "B" or "C"
for alternative systems.
'''

# =============================================================================
# Pearson and Spearman
# =============================================================================

# # This is the new, recommended method for setting seed, but not seems to be widely used
# from numpy.random import MT19937, RandomState, SeedSequence
# rs = RandomState(MT19937(SeedSequence(3221)))
# rs.seed(3221)
# rs.random.sample(5)
        

import os
result_path = os.path.dirname(os.path.realpath(__file__)) + '/results/'
del os

modelA_dct = m.run_uncertainty(modelA, seed=3221, N_sample=100, rule='L',
                               percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                               print_time=True)

pearson_rA, pearson_pA = get_correlation(modelA, input_y=key_metrics,
                                         kind='Pearson',
                                         nan_policy='raise',
                                         # file=result_path+'PearsonA.xlsx'
                                         )

spearman_rhoA, spearman_pA = get_correlation(modelA, kind='Spearman',
                                             input_y=key_metrics,
                                             nan_policy='raise',
                                             # file=result_path+'SpearmanA.xlsx'
                                             )

key_paramsA = filter_parameters(modelA, spearman_rhoA, 0.5)


# %%

# =============================================================================
# Morris one-at-a-time
# =============================================================================

modelA.set_parameters(key_paramsA)
inputs = define_inputs(modelA)
morris_samples = generate_samples(inputs, kind='Morris', N=5, seed=3221)

evaluate(modelA, morris_samples, print_time=True)

morris_dctA = morris_analysis(modelA, morris_samples, inputs,
                              metrics=key_metrics,
                              nan_policy='fill_mean', seed=3221,
                              print_to_console=True,
                              # file=result_path+'MorrisA.xlsx'
                              )

figs = []
for metric in key_metrics:
    fig = plot_morris_results(morris_dctA, metric=metric)
    fig.suptitle(metric.name)
    figs.append(fig)
    


# %%

# =============================================================================
# Sobol
# =============================================================================

saltelli_samples = generate_samples(inputs, kind='Saltelli', N=10,
                                    calc_second_order=True)

evaluate(modelA, saltelli_samples, print_time=True)

sobol_dctA = sobol_analysis(modelA, saltelli_samples, inputs,
                            metrics=key_metrics,
                            calc_second_order=True, conf_level=0.95,
                            print_to_console=True,
                            nan_policy='fill_mean',
                            # file=result_path+'SobolA.xlsx',
                            seed=3221)























