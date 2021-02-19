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

'''
TODO:
    1. Add FAST, eFAST, and RBD-FAST, and MCF
    2. Potentially add other sampling techniques for Sobol
'''


# %%

__all__ = ('get_correlation', 'define_inputs', 'generate_samples',
           'morris_analysis', 'plot_morris_results', 'sobol_analysis',)

import numpy as np
import pandas as pd
import biosteam as bst
from scipy.stats import pearsonr, spearmanr
from SALib.sample import (morris as morris_sampling, saltelli)
from SALib.analyze import morris, sobol
from matplotlib import pyplot as plt
from .utils.decorators import time_printer

isinstance = isinstance
getattr = getattr
var_indices = bst.evaluation._model.var_indices
indices_to_multiindex = bst.evaluation._model.indices_to_multiindex

def _update_input(input_val, default_val):
    if input_val is None:
        return default_val
    else:
        try:
            iter(input_val)
            return input_val if not isinstance(input_val, str) else (input_val,)
        except:
            return (input_val,)


def _update_nan(df, nan_policy, legit=('propagate', 'raise', 'omit')):
    if not nan_policy in legit:
        raise ValueError(f'nan_policy can only be in {legit}, not "{nan_policy}".')
    if nan_policy == 'propagate':
        return 'nan'
    elif nan_policy == 'raise':
        raise ValueError('"NaN" values in inputs, cannot run analysis.')
    elif nan_policy == 'omit':
        return df.dropna()
    elif nan_policy == 'fill_mean':
        return df.fillna(df.dropna().mean())
    # Shouldn't get to this step
    else:
        return df
    

# %%

def get_correlation(model, input_x=None, input_y=None,
                    kind='Pearson', nan_policy='propagate', file=''):
    '''
    Get Pearson's r between two inputs using ``scipy``.
    
    Parameters
    ----------
    model : :class:`biosteam.Model`
        Uncertainty model with defined paramters and metrics.
    input_x : :class:`biosteam.Parameter` or :class:`biosteam.Metric`
        First set of input, can be single values or iteral,
        will be defaulted to all model parameters if not provided.
    input_x : :class:`biosteam.Parameter` or :class:`biosteam.Metric`
        Second set of input, can be single values or iteral,
        will be defaulted to all model parameters if not provided.
    kind : str
        Can be "Pearson" for Pearson's r or "Spearman" for Spearman's rho
    nan_policy : str
        - "propagate": returns nan.
        - "raise": raise an error.
        - "omit": drop the pair from analysis.
    file : str
        If provided, the results will be saved as an Excel file.

    Returns
    -------
    Two :class:`pandas.DataFrame` containing Pearson'r or Spearman's rho and p-values.
    
    See Also
    --------
    `scipy.stats.pearsonr <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html>`_
    `scipy.stats.spearmanr <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html>`_
    
    '''
    table = model.table
    input_x = _update_input(input_x, model.get_parameters())
    input_y = _update_input(input_y, model.metrics)
    x_indices = var_indices(input_x)
    x_data = [table[i] for i in x_indices]
    y_indices = var_indices(input_y)
    y_data = [table[i] for i in y_indices]
    df_index = indices_to_multiindex(x_indices, ('Element', 'Input x'))
    df_column = indices_to_multiindex(y_indices, ('Element', 'Input y'))
    rs, ps = [], []
    for x in x_data:
        rs.append([])
        ps.append([])
        for y in y_data:
            df = pd.concat((x, y), axis=1)
            if True in df.isna().any().values:
                df = _update_nan(df, nan_policy)
            if isinstance(df, str):
                r, p = (np.nan, np.nan)
            else:
                if kind.capitalize() == 'Pearson':
                    r, p = pearsonr(df.iloc[:,0], df.iloc[:,1])
                    sheet_name = 'r'
                elif kind.capitalize() == 'Spearman':
                    r, p = spearmanr(df.iloc[:,0], df.iloc[:,1])
                    sheet_name = 'rho'
                else:
                    raise ValueError('kind can only be "Pearson" or "Spearman", ' \
                                      f'not "{kind}".')
            rs[-1].append(r)
            ps[-1].append(p)
    r_df = pd.DataFrame(rs, index=df_index, columns=df_column)
    p_df = pd.DataFrame(ps, index=df_index, columns=df_column)
    if file:
        with pd.ExcelWriter(file) as writer:
            r_df.to_excel(writer, sheet_name=sheet_name)
            p_df.to_excel(writer, sheet_name='p-value')
    return r_df, p_df

    
# %%

# =============================================================================
# Used for all SALib modules
# =============================================================================

def define_inputs(model):
    '''
    Define the model inputs (referred to as "problem") to be used for sampling by ``SALib``.
    
    Parameters
    ----------
    model : :class:`biosteam.Model`
        Uncertainty model with defined paramters and metrics.

    Returns
    -------
    inputs : dict
        A dict containing model inputs for sampling by ``SALib``.

    See Also
    --------
    `SALib Basics <https://salib.readthedocs.io/en/latest/basics.html#an-example>`_

    '''
    params = model.get_parameters()
    problem = {
        'num_vars': len(params),
        'names': [i.name for i in params],
        'bounds': [i.bounds if i.bounds
                   else (i.distribution.lower[0], i.distribution.upper[0])
                   for i in params]
        }
    return problem

def generate_samples(inputs, kind, N, seed=None, **kwargs):
    '''
    Generate samples for sensitivity analysis using ``SALib``.
    
    Parameters
    ----------
    model : :class:`biosteam.Model`
        Uncertainty model with defined paramters and metrics.
    inputs : dict
        A dict generated by :func:`~.sensitivity.define_inputs` to be used for ``SALib``,
        keys should include "num_vars", "names", and "bounds".
    kind : str
        Can be "Morris" (for Morris analysis) or "Saltelli" (for Sobol analysis).
    N : int
        The number of trajectories (Morris) or samples.
    seed : int
        Seed to generate a random number.
    
    Returns
    -------
    samples: array
        Samples to be used for the indicated sensitivies analyses.
    
    See Also
    --------
    `SALib.sample.morris <https://salib.readthedocs.io/en/latest/api.html?highlight=morris#method-of-morris>`_
    `SALib.sample.saltelli <https://salib.readthedocs.io/en/latest/api/SALib.sample.html?highlight=saltelli#module-SALib.sample.saltelli>`_
    '''
    if kind.capitalize() == 'Morris':
        return morris_sampling.sample(inputs, N=N, seed=seed, **kwargs)
    elif kind.capitalize() == 'Saltelli':
        return saltelli.sample(inputs, N=N, seed=seed, **kwargs)
    else:
        raise ValueError('kind can only be "Morris" or "Saltelli", ' \
                         f'not "{kind}".')


# %%

# =============================================================================
# Morris
# =============================================================================

@time_printer
def morris_analysis(model, samples, inputs, metrics=None, nan_policy='propagate',
                    conf_level=0.95, print_to_console=False,
                    print_time=False, file='', **kwargs):
    '''
    Run Morris sensitivity analysis using ``SALib``.
    
    Parameters
    ----------
    model : :class:`biosteam.Model`
        Uncertainty model with defined paramters and metrics.
    samples : :class:`numpy.array`
        Samples for Morris analysis.
    inputs : dict
        A dict generated by :func:`~.sensitivity.define_inputs` to be used for ``SALib``,
        keys should include "num_vars", "names", and "bounds".
    metrics : :class:`biosteam.Metric`
        Metrics to be included for Morris analysis, must be a subset of
        (i.e., included in the `metrics` attribute of the given model).
    nan_policy : str
        - "propagate": returns nan.
        - "raise": raise an error.
        - "fill_mean": fill nan with mean of the results.
    conf_level : float
        Confidence level of results.
    print_to_console : bool
        Whether to show results in the console.
    print_time : bool
        Whether to show simulation time in the console. 
    file : str
        If provided, the results will be saved as an Excel file.
    
    Returns
    -------
    morris_dct : dict
        A dict of Morris analysis results.
    
    See Also
    --------
    `SALib.analyze.morris <https://salib.readthedocs.io/en/latest/api.html?highlight=SALib.analyze.morris.analyze#method-of-morris>`_
    
    '''
    morris_dct = {}
    table = model.table.copy()
    table = _update_nan(table, nan_policy, legit=('propagate', 'raise', 'fill_mean'))
    if isinstance(table, str):
        table = model.table.copy()
    param_val = table.iloc[:, :len(model.get_parameters())]
    metrics = _update_input(metrics, model.metrics)
    metric_val = pd.concat([table[metric.index] for metric in metrics], axis=1)
    for metric in metrics:
        results = metric_val[metric.index]
        si = morris.analyze(inputs, param_val.to_numpy(), results.to_numpy(),
                            conf_level=conf_level, print_to_console=print_to_console,
                            **kwargs)
        morris_dct[metric.name] = si.to_df()
    if file:
        writer = pd.ExcelWriter(file)
        for name, si_df in morris_dct.items():
            si_df.to_excel(writer, sheet_name=name)
        writer.save()
    return morris_dct

def plot_morris_results(morris_dct, metric, x_axis='mu_star',
                        k1=0.1, k2=0.5, k3=1, label_kind='number',
                        single_color=False, file=''):
    '''
    Visualize the results from Morris One-at-A-Time analysis.
    
    Parameters
    ----------
    morris_dct : dict
        Results dict generated by :func:`morris_analysis`
    metric : :class:`biosteam.Metric`
        The metric of interest for the plot.
    x_axis : str
        X-axis parameter, should be either "mu_star" (the commonly used one) or "mu".
    k1 : float
        The slope to differentiate monotonic (above the line)
        and linear (below the line).
    k2 : float
        The slope to differentiate almost monotonic (above the line)
        and monotonic (below the line).
    k3 : float
        The slope to differentiate non-linear and/or non-monotonic (above the line)
        and almost monotonic (below the line).        
    label_kind : str
        How to label the points, can be either "number" (use index number of the result table)
        of "name" (use index name of the result table).
    single_color : bool
        Whether to use a single color or not in plotting the points.
    file : str
        If provided, the generated figure will be saved as a png file.
        
    Returns
    -------
    figure : :class:`matplotlib.figure.Figure`
        The generated figure.
    axis : :class:`matplotlib.axes._subplots.AxesSubplot`
        The generated figure axis.
    '''
    df = morris_dct[metric.name]
    x_data = getattr(df, x_axis)
    y_data = df.sigma
    num = len(x_data)
    if label_kind == 'number':
        labels = range(num)
    elif label_kind == 'name':
        labels = df.index.values
    else:
        raise ValueError(f'label_kind can only be "number" or "name", not "{label_kind}".')
    plt.scatter(x_data, y_data, c=range(len(x_data)) if not single_color else None)
    for x, y, label in zip(x_data, y_data, labels):
        plt.annotate(label, (x, y), xytext=(10, 10), textcoords='offset points',
                      ha='center')
    xlim = plt.gca().get_xlim()
    ylim = plt.gca().get_ylim()
    x_range = np.arange(-1, np.ceil(xlim[1])+1)
    plt.plot(x_range, k1*x_range, color='grey', linestyle='dashed')
    plt.plot(x_range, k2*x_range, color='grey', linestyle='dashed')
    plt.plot(x_range, k3*x_range, color='grey', linestyle='dashed')
    if x_axis == 'mu_star':
        xlim = (0, xlim[1])
        ylim = (0, ylim[1])
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel(r'$\mu^*$')
    plt.ylabel(r'$\sigma$')
    fig = plt.gcf()
    plt.close()
    if file:
        fig.savefig(file, dpi=300)
    return fig
    


# %%

# =============================================================================
# Sobol
# =============================================================================

@time_printer
def sobol_analysis(model, samples, inputs, metrics=None, nan_policy='propagate',
                   calc_second_order=True, conf_level=0.95, print_to_console=False,
                   print_time=False, file='', **kwargs):
    '''
    Run Sobol sensitivity analysis using ``SALib``.
    
    Parameters
    ----------
    model : :class:`biosteam.Model`
        Uncertainty model with defined paramters and metrics.
    samples : :class:`numpy.array`
        Samples for Sobol analysis.
    inputs : dict
        A dict generated by :func:`~.sensitivity.define_inputs` to be used for ``SALib``,
        keys should include "num_vars", "names", and "bounds".
    metrics : :class:`biosteam.Metric`
        Metrics to be included for Sobol analysis, must be a subset of
        (i.e., included in the `metrics` attribute of the given model).
    nan_policy : str
        - "propagate": returns nan.
        - "raise": raise an error.
        - "fill_mean": fill nan with mean of the results.
    calc_second_order : bool
        Whether to calculate second-order interaction effects.
    conf_level : float
        Confidence level of results.
    print_to_console : bool
        Whether to show results in the console.
    print_time : bool
        Whether to show simulation time in the console. 
    file : str
        If provided, the results will be saved as an Excel file.

    Returns
    -------
    si_dct : dict
        A dict of Sobol analysis results.    

    See Also
    --------
    `SALib.analyze.sobol <https://salib.readthedocs.io/en/latest/api.html#sobol-sensitivity-analysis>`_
    
    '''
    sobol_dct = {}
    metrics = _update_input(metrics, model.metrics)
    df = pd.concat([model.table[metric.index] for metric in metrics], axis=1)
    results = _update_nan(df, nan_policy, legit=('propagate', 'raise', 'fill_mean'))
    if isinstance(results, str):
        results = df
    for metric in metrics:
        result = results[metric.index]
        si = sobol.analyze(inputs, result.to_numpy(),
                           calc_second_order=calc_second_order,
                           conf_level=conf_level, print_to_console=print_to_console,
                           **kwargs)
        sobol_dct[metric.name] = dict(zip(('ST', 'S1', 'S2'), si.to_df()))
    if file:
        writer = pd.ExcelWriter(file)
        for name, si_df in sobol_dct.items():
            n_row = 0
            for df in si_df:
                df.to_excel(writer, sheet_name=name, startrow=n_row)
                n_row += len(df.index) + 2 + len(df.columns.names)
        writer.save()
    return sobol_dct









