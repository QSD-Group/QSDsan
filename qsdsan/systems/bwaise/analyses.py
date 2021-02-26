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

import pandas as pd
from qsdsan import stats as s
from qsdsan.utils.decorators import time_printer
from qsdsan.systems import bwaise as bw

m = bw.models
modelA = m.modelA
modelB = m.modelB
modelC = m.modelC

import os
result_path = os.path.dirname(os.path.realpath(__file__)) + '/results/'
figure_path = os.path.dirname(os.path.realpath(__file__)) + '/figures/'
del os

# Net cost, net GWP, and total COD/N/P/K recovery
key_metrics = [i for i in modelA.metrics if 'Net' in i.name or 'Total' in i.name]

get_param_dct = lambda model: {p.name_with_units:p for p in model.get_parameters()}

def filter_parameters(model, df, threshold):
    new_df = pd.concat((df[df>=threshold], df[df<=-threshold]))
    filtered = new_df.dropna(how='all')
    param_dct = get_param_dct(model)
    parameters = set(param_dct[i[1]] for i in filtered.index)
    return list(parameters)

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

def run_correlation(model, N, metrics=key_metrics, threshold=0.5,
                    auto_filter_parameters=True, file=''):
    if file == 'default':
        file=f'{result_path}Spearman{model._system.ID[-1]}.xlsx'
    corr_dct = m.run_uncertainty(model, seed=3221, N=N, rule='L',
                                 percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                                 print_time=True)

    spearman_rho, spearman_p = s.get_correlation(model, kind='Spearman',
                                                 input_y=metrics,
                                                 nan_policy='raise',
                                                 file=file)
    all_params = model.get_parameters()
    if auto_filter_parameters:
        key_params = filter_parameters(model, spearman_rho, threshold)
        if len(key_params) > 5:
            model.set_parameters(key_params)
    
    return (corr_dct, all_params)


# %%

# =============================================================================
# Morris one-at-a-time
# =============================================================================

def run_plot_morris(model, N, test_convergence=False, metrics=key_metrics,
                    file_prefix=''):
    inputs = s.define_inputs(model)
    
    suffix = model._system.ID[-1] if file_prefix=='default' else ''
    if file_prefix=='default':
        suffix = model._system.ID[-1]
        dct_file = f'{result_path}Morris{suffix}.xlsx'
        fig_file = f'{result_path}Morris{suffix}.png'
    else:
        dct_file = f'{file_prefix}.xlsx' if file_prefix else ''
        fig_file = f'{file_prefix}.png' if file_prefix else ''
    
    if not test_convergence:
        morris_samples = s.generate_samples(inputs, kind='Morris', N=N, seed=3221)        
        evaluate(model, morris_samples, print_time=True)

        dct = s.morris_analysis(model, inputs, metrics=metrics,
                                nan_policy='fill_mean', seed=3221,
                                file=dct_file)
        fig, ax = s.plot_morris_results(dct, metric=metrics[0], file=fig_file)
   
    else:
        dct = s.morris_till_convergence(model, inputs, metrics=metrics,
                                        N_max=N, print_time=True,
                                        file=dct_file)
        
        fig, ax = s.plot_morris_convergence(dct, metrics[0], plot_rank=True,
                                            file=fig_file)

    return (dct, fig, ax)
    

# %%

# =============================================================================
# Sobol
# =============================================================================

def run_plot_sobol(model, N, metrics=key_metrics, file_prefix=''):
    inputs = s.define_inputs(model)
    saltelli_samples = s.generate_samples(inputs, kind='Saltelli', N=N,
                                          calc_second_order=True)

    evaluate(model, saltelli_samples, print_time=True)

    if file_prefix=='default':
        suffix = model._system.ID[-1]
        dct_file = f'{result_path}Sobol{suffix}.xlsx'
        fig_file = f'{result_path}Sobol{suffix}.png'
    else:
        dct_file = f'{file_prefix}.xlsx' if file_prefix else ''
        fig_file = f'{file_prefix}.png' if file_prefix else ''

    sobol_dct = s.sobol_analysis(model, inputs,
                                 metrics=metrics,
                                 calc_second_order=True, conf_level=0.95,
                                 nan_policy='fill_mean', seed=3221,
                                 file=dct_file)
    
    fig, ax = s.plot_sobol_results(sobol_dct, metric=metrics[0],
                                   error_bar=True, annotate_heatmap=False,
                                   file=fig_file)

    return (sobol_dct, fig, ax)


# %%

# =============================================================================
# Below are testing codes, please leave as commented
# =============================================================================

# uncertainty_outs = run_correlation(modelA, N=100)

# morris_dct, fig, ax = run_plot_morris(modelA, 10, test_convergence=False)

# morris_dct_conv, fig, ax = run_plot_morris(modelA, 10, test_convergence=True)

# sobol_dct, fig, ax = run_plot_sobol(modelA, 10, file_prefix='')


# fig, ax = s.plot_morris_results(morris_dct, key_metrics[0])

# fig, ax = s.plot_morris_convergence(morris_dct_conv, key_metrics[0], plot_rank=True)


# fig, ax = s.plot_sobol_results(sobol_dct, metric=key_metrics[0], plot_in_diagonal='ST')









