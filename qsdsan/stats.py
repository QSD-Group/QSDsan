#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>
    
With contributions from:
    Yoel Rene Cortés-Peña <yoelcortes@gmail.com>
    Joy Cheung <joycheung1994@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''


# %%

__all__ = ('get_correlations', 'define_inputs', 'generate_samples',
           'morris_analysis', 'morris_till_convergence', 'sobol_analysis',
           'plot_uncertainties', 'plot_correlations',
           'plot_morris_results', 'plot_morris_convergence', 'plot_sobol_results')

from collections.abc import Iterable
import numpy as np
import pandas as pd
import seaborn as sns
import biosteam as bst
from warnings import warn
from matplotlib import pyplot as plt
from SALib.sample import (
    morris as morris_sampler,
    fast_sampler,
    latin as rbd_sampler,
    saltelli as sobol_sampler)
from SALib.analyze import morris, fast, rbd_fast, sobol
from SALib.plotting import morris as sa_plt_morris
from biosteam.plots import plot_spearman
from .utils.decorators import time_printer

isinstance = isinstance
getattr = getattr
var_indices = bst.evaluation._model.var_indices
indices_to_multiindex = bst.evaluation._model.indices_to_multiindex


# %%

# =============================================================================
# Utility functions
# =============================================================================

def _update_input(input_val, default_val):
    if input_val is None:
        return default_val
    else:
        try:
            iter(input_val)
            if len(input_val)==0: # empty iterable
                return default_val
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
    
def _update_df_names(df, columns=True, index=True):
    new_df = df.copy()

    if columns and not new_df.columns.empty:
        try:
            iter(new_df.columns)
            new_df.columns = [i[-1].split(' [')[0] for i in new_df.columns]
        except: pass

    if index and not new_df.index.empty:
        try:
            iter(new_df.index)
            new_df.index = [i[-1].split(' [')[0] for i in new_df.index]
        except: pass        

    return new_df
    

def _save_fig_return(fig, ax, file, close_fig):
    if file:
        fig.savefig(file, dpi=300)
    
    if close_fig:
        plt.close()
    
    return fig, ax


# %%

# =============================================================================
# Correlations
# =============================================================================

def get_correlations(model, input_x=None, input_y=None,
                     kind='Spearman', nan_policy='propagate', file='',
                     **kwargs):
    '''
    Get correlation coefficients between two inputs using ``scipy``.
    
    Parameters
    ----------
    model : :class:`biosteam.Model`
        Uncertainty model with defined paramters and metrics.
    input_x : :class:`biosteam.Parameter` or :class:`biosteam.Metric`
        First set of input, can be single values or an iterable,
        will be defaulted to all model parameters if not provided.
    input_y : :class:`biosteam.Parameter` or :class:`biosteam.Metric`
        Second set of input, can be single values or an iterable,
        will be defaulted to all model parameters if not provided.
    kind : str
        Can be "Spearman" for Spearman's rho, "Pearson" for Pearson's r,
        "Kendall" for Kendall's tau, or "KS" for Kolmogorov–Smirnov's D.
    nan_policy : str
        - "propagate": returns nan.
        - "raise": raise an error.
        - "omit": drop the pair from analysis.
    file : str
        If provided, the results will be saved as an Excel file.
    kwargs
        Other kwargs that will be passed to ``scipy``.

    Returns
    -------
    Two :class:`pandas.DataFrame` containing the test statistics and p-values.
    
    See Also
    --------    
    :func:`scipy.stats.spearmanr`

    :func:`scipy.stats.pearsonr`
    
    :func:`scipy.stats.kendalltau`
    
    :func:`scipy.stats.kstest`
    
    '''

    if input_x and not isinstance(input_x, Iterable): input_x = (input_x,)
    if input_y and not isinstance(input_y, Iterable): input_y = (input_y,)
    if nan_policy not in ('propagate', 'raise', 'omit'):
        raise ValueError(f'nan_policy can only be in ("omit", "propagate", "raise"), ' \
                         f'not "{nan_policy}".')

    name = kind.lower()
    if name == 'spearman':
        correlation = model.spearman_r
        sheet_name = 'rho'
    elif name == 'pearson':
        correlation = model.pearson_r
        sheet_name = 'r'
    elif name == 'kendall':
        correlation = model.kendall_tau
        sheet_name = 'tau'
    elif name == 'ks':
        correlation = model.kolmogorov_smirnov_d
        sheet_name = 'D'
    else:
        raise ValueError('kind can only be "Spearman", "Pearson", ' 
                        f'"Kendall", or "KS", not "{kind}".')
    r_df, p_df = dfs = correlation(input_x, input_y, nan_policy+' nan', **kwargs)
    for df in dfs:
        df.index.names = ('Element', 'Input x')
        df.columns.names = ('Element', 'Input y')
        
    if file:
        with pd.ExcelWriter(file) as writer:
            r_df.to_excel(writer, sheet_name=sheet_name)
            p_df.to_excel(writer, sheet_name='p-value')
    return r_df, p_df

    
# %%

# =============================================================================
# SALib modules
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
    return model.problem()


# TODO: incorporate sampling from different distributions
def generate_samples(inputs, kind, N, seed=None, **kwargs):
    '''
    Generate samples for sensitivity analysis using ``SALib``.
    
    Parameters
    ----------
    model : :class:`biosteam.Model`
        Uncertainty model with defined paramters and metrics.
    inputs : dict
        A dict generated by :func:`qsdsan.sensitivity.define_inputs` to be used for ``SALib``,
        keys should include "num_vars", "names", and "bounds".
    kind : str
        Can be "FAST" for Fourier amplitude sensitivity test (FAST),
        "RBD" for random balance design FAST (latin hypercube sampling [LHS]),
        "Morris" for Morris One-at-A-Time (OAT),
        or "Sobol" for Sobol analysis (Saltelli sampling).
    N : int
        The number of samples or trajectories (Morris).
    seed : int
        Seed to generate random samples.
    kwargs
        Other kwargs that will be passed to ``SALib``.
    
    Returns
    -------
    samples: array
        Samples to be used for the indicated sensitivies analyses.
    
    See Also
    --------
    ``SALib`` `API <https://salib.readthedocs.io/en/latest/api.html>`_
    
    '''
    lower = kind.lower()
    if lower == 'morris':
        return morris_sampler.sample(inputs, N=N, seed=seed, **kwargs)
    elif lower in ('fast', 'efast'):
        return fast_sampler.sample(inputs, N=N, seed=seed, **kwargs)
    elif lower == 'rbd':
        return rbd_sampler.sample(inputs, N=N, seed=seed, **kwargs)
    elif lower == 'sobol':
        return sobol_sampler.sample(inputs, N=N, seed=seed, **kwargs)
    else:
        raise ValueError('kind can only be "FAST", "RBD", "Morris", or "Sobol", ' \
                         f'not "{kind}".')

# =============================================================================
# Morris
# =============================================================================

@time_printer
def morris_analysis(model, inputs, metrics=None, nan_policy='propagate',
                    conf_level=0.95, print_to_console=False,
                    print_time=False, file='', **kwargs):
    '''
    Run Morris sensitivity analysis using ``SALib``.
    
    Parameters
    ----------
    model : :class:`biosteam.Model`
        Uncertainty model with defined paramters and metrics.
    inputs : dict
        A dict generated by :func:`qsdsan.sensitivity.define_inputs` to be used for ``SALib``,
        keys should include "num_vars", "names", and "bounds".
    metrics : :class:`biosteam.Metric`
        Metrics to be included for Morris analysis, must be a subset of
        the metrics of the model to be analyzed.
        (i.e., included in the `metrics` attribute of the given model).
        If None is provided, all metrics in the model will be included.
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
    kwargs
        Other kwargs that will be passed to ``SALib``.
    
    Returns
    -------
    morris_dct : dict
        A dict of Morris analysis results.
    
    See Also
    --------
    ``SALib`` `API <https://salib.readthedocs.io/en/latest/api.html>`_
    
    '''
    morris_dct = {}
    model = model.copy()
    metrics = _update_input(metrics, model.metrics)
    model.metrics = metrics
    table = model.table.astype('float64')
    table = _update_nan(table, nan_policy, legit=('propagate', 'raise', 'fill_mean'))
    
    if isinstance(table, str):
        table = model.table.astype('float64')
    
    param_val = table.iloc[:, :len(model.get_parameters())]
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


@time_printer
def morris_till_convergence(model, inputs, metrics=None,
                            N_max=20, seed=None, threshold=0.1,
                            nan_policy='propagate',
                            conf_level=0.95, print_to_console=False,
                            print_time=False, file='', **kwargs):
    '''
    Run Morris analysis from N=2 to N=N_max until the results converge
    (i.e., mu_star_conf/mu_star_max < threshold for all parameters,
     where as mu_star_max is the maximum :math:`{\mu^*}` value for a certain metric,
     and this should be satisfied for all metrics).
    
    Parameters
    ----------
    model : :class:`biosteam.Model`
        Uncertainty model with defined paramters and metrics.
    inputs : dict
        A dict generated by :func:`qsdsan.sensitivity.define_inputs` to be used for ``SALib``,
        keys should include "num_vars", "names", and "bounds".
    metrics : :class:`biosteam.Metric`
        Metrics to be included for Morris analysis, must be a subset of
        the metrics of the model to be analyzed.
        (i.e., included in the `metrics` attribute of the given model).
        If None is provided, all metrics in the model will be included.
    N_max : int
        Maximum number of trajectories to be considered.
    seed : int
        Seed to generate random samples.
    threshold : float
        Threshold for the convergence.
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
    kwargs
        Other kwargs that will be passed to ``SALib``.

    See Also
    --------    
    :func:`qsdsan.stats.morris_analysis`
    
    '''
    num_levels = kwargs['num_levels'] if 'num_levels' in kwargs.keys() else 4
    kwargs = {i:kwargs[i] for i in kwargs.keys() if i!='num_levels'}
    samples = generate_samples(inputs=inputs, kind='Morris', N=N_max,
                               seed=seed, num_levels=num_levels)
    model.load_samples(samples)
    
    param_num = len(model.get_parameters())
    cum_model = model.copy()
    cum_model.load_samples(samples[0: 2*(param_num+1)])
    cum_model.evaluate()
    cum_dct = dict(mu_star={}, mu_star_conf={})
    metrics = _update_input(metrics, model.metrics)
    temp_dct = morris_analysis(model=cum_model, inputs=inputs, metrics=metrics,
                               nan_policy=nan_policy, conf_level=conf_level,
                               print_to_console=print_to_console, **kwargs)
    for m in metrics:
        for idx in ('mu_star', 'mu_star_conf'):
            data0 = getattr(temp_dct[m.name], idx)
            df = pd.DataFrame(columns=data0.index, index=(2,))
            df.index.name = idx
            df.loc[2] = data0.copy()
            cum_dct[idx][m.name] = df
    
    for n in range(2, N_max):
        temp_model = model.copy()
        temp_model.load_samples(samples[n*(param_num+1): (n+1)*(param_num+1)])
        temp_model.evaluate()
        cum_model.table = pd.concat((cum_model.table, temp_model.table))
        
        temp_dct = morris_analysis(model=cum_model, inputs=inputs, metrics=metrics,
                                   nan_policy=nan_policy, conf_level=conf_level,
                                   print_to_console=print_to_console, **kwargs)
        all_converged = True
        for m in metrics:
            mu_star = temp_dct[m.name].mu_star
            mu_star_conf = temp_dct[m.name].mu_star_conf
            cum_dct['mu_star'][m.name].loc[n+1] = mu_star
            cum_dct['mu_star_conf'][m.name].loc[n+1] = mu_star_conf
            
            converged = False if (mu_star_conf/mu_star.max()>threshold).any() else True
            all_converged = all_converged & converged

        if all_converged:
            print(f'mu_star converges at {n} trajectories.')
            break
        elif n == N_max-1:
            print(f'mu_star has not converged with {n} trajectories.')
    
    if file:
        writer = pd.ExcelWriter(file)
        for m in metrics:
            cum_dct['mu_star'][m.name].to_excel(writer, sheet_name=m.name)
            cum_dct['mu_star_conf'][m.name].to_excel(
                writer, sheet_name=m.name, startrow=N_max)
        writer.save()
    
    return cum_dct


# =============================================================================
# (e)FAST and RBD-FAST
# =============================================================================

@time_printer
def fast_analysis(model, inputs, kind, metrics=None, nan_policy='propagate',
                  conf_level=0.95, print_to_console=False, print_time=False,
                  file='', **kwargs):
    '''
    Run Fourier amplitude sensitivity test (Saltelli's extended FAST) or 
    random balance design (RBD) FAST using ``SALib``.
    
    Parameters
    ----------
    model : :class:`biosteam.Model`
        Uncertainty model with defined paramters and metrics.
    inputs : dict
        A dict generated by :func:`qsdsan.sensitivity.define_inputs` to be used for ``SALib``,
        keys should include "num_vars", "names", and "bounds".
    kind : str
        Either "FAST" or "RBD".
    metrics : :class:`biosteam.Metric`
        Metrics to be included for Morris analysis, must be a subset of
        the metrics of the model to be analyzed.
        (i.e., included in the `metrics` attribute of the given model).
        If None is provided, all metrics in the model will be included.
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
    kwargs
        Other kwargs that will be passed to ``SALib``.
    
    Returns
    -------
    fast_dct : dict
        A dict of FAST analysis results.
    
    See Also
    --------
    ``SALib`` `API <https://salib.readthedocs.io/en/latest/api.html>`_
    
    '''
    fast_dct = {}
    model = model.copy()
    metrics = _update_input(metrics, model.metrics)
    model.metrics = metrics
    table = model.table.astype('float64')
    table = _update_nan(table, nan_policy, legit=('propagate', 'raise', 'fill_mean'))
    
    if isinstance(table, str):
        table = model.table.astype('float64')
        
    param_val = table.iloc[:, :len(model.get_parameters())]
    metric_val = pd.concat([table[metric.index] for metric in metrics], axis=1)
    
    if kind.lower() in ('fast', 'efast'):
        for metric in metrics:
            results = metric_val[metric.index]
            si = fast.analyze(inputs, results.to_numpy(), conf_level=conf_level,
                              print_to_console=print_to_console, **kwargs)
            fast_dct[metric.name] = si.to_df()
    
    elif kind.lower() == 'rbd':
        for metric in metrics:
            results = metric_val[metric.index]
            si = rbd_fast.analyze(inputs, param_val.to_numpy(), results.to_numpy(),
                                  conf_level=conf_level, print_to_console=print_to_console,
                                  **kwargs)
            fast_dct[metric.name] = si.to_df()
    
    else:
        raise ValueError(f'kind can only be "FAST" or "RBD", not "{kind}".')
    
    if file:
        writer = pd.ExcelWriter(file)
        for name, si_df in fast_dct.items():
            si_df.to_excel(writer, sheet_name=name)
        writer.save()
    
    return fast_dct


# =============================================================================
# Sobol
# =============================================================================

@time_printer
def sobol_analysis(model, inputs, metrics=None, nan_policy='propagate',
                   calc_second_order=True, conf_level=0.95, print_to_console=False,
                   print_time=False, file='', **kwargs):
    '''
    Run Sobol sensitivity analysis using ``SALib``.
    
    Parameters
    ----------
    model : :class:`biosteam.Model`
        Uncertainty model with defined paramters and metrics.
    inputs : dict
        A dict generated by :func:`qsdsan.sensitivity.define_inputs` to be used for ``SALib``,
        keys should include "num_vars", "names", and "bounds".
    metrics : :class:`biosteam.Metric`
        Metrics to be included for Morris analysis, must be a subset of
        the metrics of the model to be analyzed.
        (i.e., included in the `metrics` attribute of the given model).
        If None is provided, all metrics in the model will be included.
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
    kwargs
        Other kwargs that will be passed to ``SALib``.

    Returns
    -------
    si_dct : dict
        A dict of Sobol analysis results.    

    See Also
    --------
    ``SALib`` `API <https://salib.readthedocs.io/en/latest/api.html>`_
    
    '''
    sobol_dct = {}
    metrics = _update_input(metrics, model.metrics)
    model = model.copy()
    model.metrics = metrics
    table = model.table.astype('float64')
    
    df = pd.concat([table[metric.index] for metric in metrics], axis=1)
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


# %%

# =============================================================================
# Plot uncertainty analysis results
# =============================================================================

def plot_uncertainties(model, x_axis=(), y_axis=(), kind='box',
                       file='', close_fig=True, center_kws={}, margin_kws={}):
    '''
    Visualize uncertainty analysis results as one of the following depending on inputs:
    
    +---------------------------------+------------------------------------------+
    | input                           | returned plot                            |
    +----------+----------+-----------+-------------+-----------+----------------+    
    | x_axis   | y_axis   | kind      | orientation | center    | margin         |
    +==========+==========+===========+=============+===========+================+
    | single   | None     | box       | 1D          | box       | N/A            |
    | or       |          +-----------+ horizontal  +-----------+                +
    | iterable |          | hist      |             | histogram |                |
    |          |          +-----------+             +-----------+                +
    |          |          | kde       |             | kde       |                |
    +----------+----------+-----------+-------------+-----------+----------------+
    | None     | single   | box       | 1D          | box       | N/A            |
    |          | or       +-----------+ vertical    +-----------+                +
    |          | iterable | hist      |             | histogram |                |
    |          |          +-----------+             +-----------+                +
    |          |          | kde       |             | kde       |                |
    +----------+----------+-----------+-------------+-----------+----------------+
    | single   | single   | hist-box  | 2D          | histogram | box            |
    |          |          +-----------+             +           +----------------+
    |          |          | hist-kde  |             |           | kernel density |
    |          |          +-----------+             +           +----------------+
    |          |          | hist-hist |             |           | histogram      |
    |          |          +-----------+             +-----------+----------------+
    |          |          | kde-box   |             | kde       | box            |
    |          |          +-----------+             +           +----------------+
    |          |          | kde-hist  |             |           | histogram      |
    |          |          +-----------+             +           +----------------+ 
    |          |          | kde-kde   |             |           | kernel density |
    +----------+----------+-----------+-------------+-----------+----------------+ 

    .. note::
        When both x_axis and y_axis are not None (i.e., the figure is 2D),
        for clarity reasons, this function requires len(x_axis) and len(y_axis) to be both 1.

        If wanted, it is possible to use colors (``hue`` in ``seaborn``) to
        differentiated the different parameters or metrics in either x_axis or
        y_axis. To achieve this, please refer to the documentation of ``seaborn``.
        
        Similarly, there are other potential combinations of the center and
        margin plots that are supported by ``seaborn`` but not included here.

    Parameters
    ----------
    model : :class:`biosteam.Model`
        The model with uncertainty analysis (in <:class:`Model`.table>) results for plotting.
    x_axis : :class:`biosteam.Parameter`, :class:`biosteam.Metric` or iterable
        What to plot on the x-axis, can be parameters or metrics of the model,
        default to all model metrics included in the model result table
        if neither x nor y is provided.
    y_axis : :class:`biosteam.Parameter`, :class:`biosteam.Metric` or iterable
        What to plot on the y-axis, can be parameters or metrics of the model,
        default to None.
    kind : str
        What kind of plot to be returned, refer to the summary table for valid inputs.
    file : str
        If provided, the generated figure will be saved as a png file.
    close_fig : bool
        Whether to close the figure
        (if not close, new figure will be overlaid on the current figure).
    center_kws : dict
        Will be passed to ``seaborn`` for the center plot.
    margin_kws : dict
        Will be passed to ``seaborn`` for the margin plots.
        
    Returns
    -------
    figure : :class:`matplotlib.figure.Figure`
        The generated figure.
    axis : :class:`matplotlib.axes._subplots.AxesSubplot` or iterable
        The generated figure axis (or axes for 2D figure).
        
    See Also
    --------    
    :func:`seaborn.jointplot` `docs <https://seaborn.pydata.org/generated/seaborn.jointplot.html>`_
    
    '''

    kind_lower = kind.lower()
    table = model.table.astype('float64')
    df = _update_df_names(table)
    
    twoD = False
    x_df = y_df = None
    sns_df = pd.DataFrame(columns=('x_group', 'x_data', 'y_group', 'y_data'))
    
    dfs = []
    if not y_axis: # no data provided or only x, 1D, horizontal
        x_axis = _update_input(x_axis, model.metrics)
        x_df = df[[i.name for i in x_axis]]
        
        for x in x_df.columns:
            temp_df = sns_df.copy()
            temp_df['x_data'] = x_df[x]
            temp_df['x_group'] = x
            dfs.append(temp_df)

        sns_df = pd.concat(dfs)

    elif not x_axis: # only y, 1D, vertical
        y_df = df[[i.name for i in y_axis]]
    
        for y in y_df.columns:
            temp_df = sns_df.copy()
            temp_df['y_data'] = y_df[y]
            temp_df['y_group'] = y
            dfs.append(temp_df)
        
        sns_df = pd.concat(dfs)
    
    else: # x and y, 2D
        len_x = len_y = 1
        try:
            len_x = len(x_axis)
        except:
            try:
                len_y = len(y_axis)
            except:
                pass
        if not (len_x==1 and len_y==1):
            raise ValueError('Only single input allowed when both x_axis and y_axis are provided.')
        
        twoD = True
        x_name, y_name = x_axis.name, y_axis.name
        sns_df = df[[x_name, y_name]]
        
    sns.set_theme(style='ticks')
    if not twoD: # 1D plot
        if kind_lower == 'box':
            center_kws.setdefault('dodge', False)
            
            if x_df is not None: # horizontal box plot
                ax = sns.boxplot(data=sns_df, x='x_group', y='x_data', hue='x_group',
                                 **center_kws)
                xlabel = ''
                ylabel = 'Values'
            else:
                ax = sns.boxplot(data=sns_df, x='y_data', y='y_group', hue='y_group',
                                 **center_kws)
                xlabel = 'Values'
                ylabel = ''

            ax.set(xlabel=xlabel, ylabel=ylabel)
                
        elif kind_lower == 'hist':
            if x_df is not None: # horizontal hist plot
                ax = sns.histplot(data=sns_df, x='x_data', hue='x_group', **center_kws)
                ax.set_xlabel('Values')
            
            else: # vertical hist plot
                ax = sns.histplot(data=sns_df, y='y_data', hue='y_group', **center_kws)
                ax.set_ylabel('Values')
                
        elif kind_lower == 'kde':
            if x_df is not None: # horizontal kde plot
                ax = sns.kdeplot(data=sns_df, x='x_data', hue='x_group', **center_kws)
            else: # vertical kde plot
                ax = sns.kdeplot(data=sns_df, y='y_data', hue='y_group', **center_kws)
        
        else:
            raise ValueError('kind can only be "box", "kde", "hist", or "hist-kde", ' \
                             f'for 1D plot, not "{kind}".')
        
        ax.get_legend().set_title('')
        return _save_fig_return(ax.figure, ax, file, close_fig)
        
    else: # 2D plot
        g = sns.JointGrid(data=sns_df, x=x_name, y=y_name)

        kind_split = kind_lower.split('-')
        if kind_split[0] in ('hist', 'kde'):
            func = getattr(sns, f'{kind_split[0]}plot')
            g.plot_joint(func, **center_kws)
        else:
            raise ValueError(f'The provided kind "{kind}" is not valid for 2D plot.')

        if kind_split[1] in ('box', 'hist', 'kde'):
            func = getattr(sns, f'{kind_split[1]}plot')
            g.plot_marginals(func, **margin_kws)
        else:
            raise ValueError(f'The provided kind "{kind}" is not valid for 2D plot.')
        
        return _save_fig_return(g.fig, g.fig.axes, file, close_fig)


# =============================================================================
# Plot correlations
# =============================================================================

def _plot_corr_tornado(corr_df, top):
    fig, ax = plot_spearman(corr_df.iloc[:,0], top=top)
    
    ax.set_xlabel(corr_df.columns[0])
    
    for ax in fig.axes:
        for key in ax.spines.keys():
            ax.spines[key].set(color='k', linewidth=0.5, visible=True)
            ax.grid(False)
    
    return fig, ax


def _plot_corr_bubble(corr_df, ratio, **kwargs):
    sns.set_theme(style="whitegrid")

    margin_x = kwargs['margin_x'] if 'margin_x' in kwargs.keys() else 0.1/ratio
    margin_y = kwargs['margin_y'] if 'margin_y' in kwargs.keys() else 0.1
    kwargs = {i: kwargs[i] for i in kwargs.keys() if 'margin' not in i}
    
    keys = ('height', 'palette', 'hue_norm', 'sizes', 'size_norm', 'edgecolor')
    values = (5+ratio, 'vlag', (-1, 1), (0, 1000), (0, 2), '0.5')
    
    for num, k in enumerate(keys):
        kwargs.setdefault(keys[num], values[num])
    
    g = sns.relplot(data=corr_df, x='metric', y='parameter',
                    hue='correlation', size='size',  **kwargs)
    
    g.set(xlabel='', ylabel='', aspect=1)
    g.ax.margins(x=margin_x, y=margin_y)
    
    for label in g.ax.get_xticklabels():
        label.set_rotation(90)
    
    for artist in g.legend.legendHandles:
        artist.set_edgecolor('0.5')
    
    for key in g.ax.spines.keys():
        g.ax.spines[key].set(color='k', linewidth=0.5, visible=True)
    
    g.ax.grid(True, which='major', color='k',linestyle='--', linewidth=0.7)
    g.tight_layout()
    
    return g
    

def plot_correlations(result_df, parameters=(), metrics=(), top=None,
                      file='', close_fig=True, **kwargs):
    '''
    Visualize the correlations between model parameters and metric results
    as tornado (single metric) or bubble plots (multiple metrics).
    
    Parameters
    ----------
    result_df : :class:`pandas.DataFrame`
        Result table generated by :func:`get_correlations` containing
        correlation indices.
    parameters : :class:`biosteam.Parameter`
        Metric(s) of interest for the plot, will be default to all parameters
        included in ``corr_dct`` if not provided.
    metrics : :class:`biosteam.Metric`
        Metric(s) of interest for the plot, will be default to all metrics
        included in ``corr_dct`` if not provided.
    top : int
        Plot the top X parameters with the highest absolute correlation indices,
        this is only applicable for the case of just one metric.
    file : str
        If provided, the generated figure will be saved as a png file.
    close_fig : bool
        Whether to close the figure
        (if not close, new figure will be overlaid on the current figure).
    kwargs
        Other kwargs that will be passed to :func:`seaborn.relplot`.
        
    Returns
    -------
    figure : :class:`matplotlib.figure.Figure`
        The generated figure.
    axis : :class:`matplotlib.axes._subplots.AxesSubplot`
        The generated figure axis.
        
    See Also
    --------
    :func:`biosteam.plots.plot_spearman`
    
    :func:`seaborn.relplot` and `scatter heat map <https://seaborn.pydata.org/examples/heat_scatter.html>`_
    
    '''
    
    df = _update_df_names(result_df)

    param_names = _update_input(parameters, df.index)
    param_names = param_names if isinstance(param_names[0], str) \
                              else [p.name for p in param_names]

    metric_names = _update_input(metrics, df.columns)    
    metric_names = metric_names if isinstance(metric_names[0], str) \
                                else [m.name for m in metric_names]
    
    df = df[metric_names].loc[param_names]
    
    if len(param_names)*len(metric_names) == 0:
        raise ValueError('No correlation data for plotting.')
    
    elif len(metric_names) == 1: # one metric, tornado plot
        fig, ax = _plot_corr_tornado(df, top)
        return _save_fig_return(fig, ax, file, close_fig)
    
    else: # multiple metrics, bubble plot
        corr_df = pd.DataFrame()
        
        for m in df.columns:
            temp_df = pd.DataFrame(columns=('parameter', 'metric', 'correlation', 'size'))
            temp_df['parameter'] = df.index
            temp_df['metric'] = [m]*temp_df.shape[0]
            temp_df['correlation'] = df[m].values
            temp_df['size'] = np.abs(df[m].values)
            corr_df = pd.concat((corr_df, temp_df))
        
        g = _plot_corr_bubble(corr_df, len(metric_names)/len(param_names), **kwargs)
        
        return _save_fig_return(g.fig, g.ax, file, close_fig)
    


# =============================================================================
# Plot Morris analysis results
# =============================================================================

def plot_morris_results(morris_dct, metric, kind='scatter',
                        x_axis='mu_star',
                        k1=0.1, k2=0.5, k3=1, label_kind='number',
                        file='', close_fig=True, **kwargs):
    '''
    Visualize the results from Morris One-at-A-Time analysis as either scatter
    or bar plots.
    In scatter plots, the x values are :math:`{\mu^*}` and the y values are :math:`{\sigma}`.
    In bar plots, bar length indicate the :math:`{\mu^*}` values with error bars
    representing confidence intervals of the analysis.
    
    Parameters
    ----------
    morris_dct : dict
        Results dict generated by :func:`morris_analysis`.
    metric : :class:`biosteam.Metric`
        The metric of interest for the plot.
    kind : str
        Either "scatter" (:math:`{\sigma}`) vs. :math:`{\mu^*}` or "bar" (:math:`{\mu^*}` with confidence interval) plot.
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
    file : str
        If provided, the generated figure will be saved as a png file.
    close_fig : bool
        Whether to close the figure
        (if not close, new figure will be overlaid on the current figure).
    kwargs
        Other kwargs that will be passed to :func:`morris.horizontal_bar_plot` in ``SALib.plotting``.
        
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
    
    ax = plt.subplot()
    sns.set_theme(style='ticks')
    
    if kind == 'scatter':
        ax.scatter(x_data, y_data, color='k')
        for x, y, label in zip(x_data, y_data, labels):
            ax.annotate(label, (x, y), xytext=(10, 10), textcoords='offset points',
                        ha='center')
        x_range = np.arange(-1, np.ceil(ax.get_xlim()[1])+1)
        line1, = ax.plot(x_range, k1*x_range, color='black', linestyle='-.')
        line2, = ax.plot(x_range, k2*x_range, color='black', linestyle='--')
        line3, = ax.plot(x_range, k3*x_range, color='black', linestyle='-')
        ax.legend((line3, line2, line1), (r'$\sigma/\mu^*$'+f'={k3}',
                                          r'$\sigma/\mu^*$'+f'={k2}',
                                          r'$\sigma/\mu^*$'+f'={k1}'),
                  loc='best')
        if x_axis == 'mu_star':
            ax.set_xlim(0,)
        ax.set_ylim(0,)
        ax.set_xlabel(r'$\mu^*$')
        ax.set_ylabel(r'$\sigma$')
        fig = ax.figure
    
    elif kind == 'bar':
        if x_axis == 'mu':
            raise ValueError('Bar plot can only be made for mu_star, not mu.')
        df = morris_dct[metric.name]
        df['names'] = df.index
        fig = sa_plt_morris.horizontal_bar_plot(ax, df, opts=kwargs)

    # for ax in fig.axes:
    for key in ax.spines.keys():
        ax.spines[key].set(color='k', linewidth=0.5, visible=True)
        ax.grid(False)

    return _save_fig_return(fig, ax, file, close_fig)


def plot_morris_convergence(result_dct, metric, parameters=(), plot_rank=False,
                            kind='line', show_error=True, palette='pastel', file='', close_fig=True):
    '''
    Plot the evolution of :math:`{\mu^*}` or its rank with the number of trajectories.
    
    Parameters
    ----------
    result_dct : dict
        Result dictionary generated from :func:`qsdsan.stats.morris_till_convergence`
    metric : :class:`biosteam.Metric`
        The metric of interest for the plot.
    parameters : :class:`biosteam.Parameter`
        Single or a iterable of model parameters whose :math:`{\mu^*}` will be
        included in the plot.
        Will be set to all parameters in retult_dct will be used if not provided.
    plot_rank : bool
        If True, will plot rank of :math:`{\mu^*}` instead of value.
        
        .. note::
            If plot_rank is True, error bars or bands will not be included.
    kind : str
        Can be either 'line' or 'scatter'.
    show_error : bool
        Whether to include the confidence interval in the plot,
        will be bars for scatter plot and bands for line plot.
    palette : string, list, dict, or :class:`matplotlib.colors.Colormap`
        Will be passed on to :func:`seaborn.color_palette`.
    file : str
        If provided, the generated figure will be saved as a png file.
    close_fig : bool
        Whether to close the figure
        (if not close, new figure will be overlaid on the current figure).
    
    Returns
    -------
    figure : :class:`matplotlib.figure.Figure`
        The generated figure.
    axis : :class:`matplotlib.axes._subplots.AxesSubplot`
        The generated figure axis.
    
    '''
    ax = plt.subplot()
    df = result_dct['mu_star'][metric.name].copy().astype('float64')
    conf_df = result_dct['mu_star_conf'][metric.name].copy().astype('float64')

    param_names = _update_input(parameters, df.columns)    
    param_names = param_names if isinstance(param_names[0], str) \
                              else [p.name for p in param_names]
    
    if plot_rank:
        df = df.rank(axis=1)
        ylabel = f'Rank for {metric.name.lower()}'
        loc = 'lower left'
    else:
        ylabel = f'$\mu^*$ for {metric.name.lower()}'
        loc = 'best'
    
    palette = sns.color_palette('deep', n_colors=len(param_names))
    sns.set_theme(style='ticks', palette=palette)
    
    for n, param in enumerate(param_names):
        if kind == 'line':
            ax.plot(df.index, df[param], color=palette[n], linewidth=1.5, label=param)
            if not plot_rank and show_error:
                ax.fill_between(df.index, df[param]-conf_df[param], df[param]+conf_df[param],
                                color=palette[n], linewidth=0, alpha=0.2)
        elif kind == 'scatter':
            ax.scatter(df.index, df[param], color=palette[n], label=param)
            if not plot_rank and show_error:
                ax.errorbar(df.index, df[param], conf_df[param],
                            color=palette[n], alpha=0.5)
        else:
            raise ValueError(f'kind can only be "line" or "scatter", not "{kind}".')
    
    ax.legend(loc=loc)
    ax.set(xlabel='Number of trajectories', ylabel=ylabel, ylim=(0, ax.get_ylim()[1]))
    
    return _save_fig_return(ax.figure, ax, file, close_fig)
    

# =============================================================================
# Plot variance breakdown
# =============================================================================

def _plot_bar(kind, df, error, ax=None):
    ax = ax if ax else plt.subplot()

    sns.set_theme(style='white')
    if 'ST' in kind:
        sns.set_color_codes('pastel')
        sns.barplot(x=df.ST, y=df.index, data=df,
                    ax=ax, label='Total', color='b')
        if error:
            ax.errorbar(df.ST, df.index, xerr=df.ST_conf, fmt='none', ecolor='b')
    
    if 'S1' in kind:
        sns.set_color_codes('muted')
        sns.barplot(x=df.S1, y=df.index, data=df,
                    ax=ax, label='Main', color='b')
        if error:
            ax.errorbar(df.S1, df.index, xerr=df.S1_conf, fmt='none', ecolor='b')        
    
    ax.set_xlabel('Variance')
    ax.legend(ncol=2, loc='lower right', frameon=True)
    ax.set_ylim(df.shape[0]-0.5, -0.5)

    return ax

def _plot_heatmap(hmap_df, ax=None, annot=False, diagonal='', sts1_df=None,
                  default_cbar=True):
    ax = ax if ax else plt.subplot()
    ax_cbar = ax.figure.add_axes([0.03, 0.3, 0.02, 0.4]) if not default_cbar else None

    if diagonal:
        np.fill_diagonal(hmap_df.values, getattr(sts1_df, diagonal))
        hmap_df = hmap_df.astype('float64')
        k = -1
        title = 'Total/Interaction Effects' if diagonal=='ST' else 'Main/Interaction Effects'
    else:
        hmap_df = hmap_df.fillna(0)
        k = 0
        title = 'Interaction Effects'
    
    mask = np.tril(np.ones_like(hmap_df, dtype=bool), k)
    
    sns.set_theme(style='white')
    cmap = sns.diverging_palette(230, 20, as_cmap=True)
    sns.heatmap(hmap_df,
                mask=mask,
                ax=ax, cmap=cmap, center=0, linewidths=.5,
                annot=annot, cbar_kws={'shrink': 0.5}, cbar_ax=ax_cbar)
    ax.set_title(title)
    
    return ax


def plot_fast_results(result_dct, metric, parameters=(), error_bar=True,
                      file='', close_fig=True):
    '''
    Visualize the results from FAST or RBD-FAST analysis as a bar plot.
    
    Parameters
    ----------
    result_dct : dict
        Result dictionary generated from :func:`qsdsan.stats.fast_analysis`
    metric : :class:`biosteam.Metric`
        The metric of interest for the plot.
    parameters : :class:`biosteam.Parameter`
        Single or a iterable of model parameters whose :math:`{\mu^*}` will be
        included in the plot.
        Will be set to all parameters in retult_dct will be used if not provided.
    error_bar : bool
        Whether to include the confidence interval as error bars in the plot.
    file : str
        If provided, the generated figure will be saved as a png file.
    close_fig : bool
        Whether to close the figure
        (if not close, new figure will be overlaid on the current figure).
        
    Returns
    -------
    figure : :class:`matplotlib.figure.Figure`
        The generated figure.
    axis : :class:`matplotlib.axes._subplots.AxesSubplot`
        The generated figure axis.
        
    '''
    
    param_names = _update_input(parameters, result_dct[metric.name]['S1'].index)    
    param_names = param_names if isinstance(param_names[0], str) \
                              else [p.name for p in param_names]
    
    df = result_dct[metric.name]
    kind = 'STS1' if 'ST' in df.columns else 'S1'
    
    ax = _plot_bar(kind, df, error_bar)
    
    return _save_fig_return(ax.figure, ax, file, close_fig)


def plot_sobol_results(result_dct, metric, parameters=(), kind='all',
                       annotate_heatmap=False, plot_in_diagonal='',
                       error_bar=True, file='', close_fig=True):
    '''
    Visualize the results from Sobol analysis as a bar plot and/or heat map.
    Total (:math:`S_{Ti}`) and main (:math:`S_{1i}`) effects can be drawn in
    the bar plot or diagnoal of the heat map;
    second-order interaction effects between two parameters (:math:`S_{2ij}`)
    are shown in the heat map.

    Parameters
    ----------
    result_dct : dict
        Result dictionary generated from :func:`qsdsan.stats.sobol_analysis`
    metric : :class:`biosteam.Metric`
        The metric of interest for the plot.
    parameters : :class:`biosteam.Parameter`
        Single or a iterable of model parameters whose :math:`{\mu^*}` will be
        included in the plot.
        Will be set to all parameters in retult_dct will be used if not provided.
    kind : str
        Which sensitivity index or indices to plot:
            
        +------+---------------------------------------------------------------+
        | kind | returned plot                                                 |
        +======+===============================================================+
        | ST   | total effects (bar)                                           |
        +------+---------------------------------------------------------------+
        | S1   | main effects (bar)                                            |
        +------+---------------------------------------------------------------+
        | S2   | interaction effects (heat map)                                |
        +------+---------------------------------------------------------------+
        | STS1 | total and main effects (bar)                                  |
        +------+---------------------------------------------------------------+
        | STS2 | total and interaction effects (heat map or bar and heat map)  |
        +------+---------------------------------------------------------------+
        | S1S2 | main and interaction effects (heat map or bar and heat map)   |
        +------+---------------------------------------------------------------+
        | all  | all effects (bar and heat map)                                |
        +------+---------------------------------------------------------------+
    annotate_heatmap : bool
        Whether to annotate the index values in the heat map.
    plot_in_diagonal : str
        Plot total or main effects in the diagonal of the interaction heat map,
        can be "ST", "S1", or "".
        This is applicable when kind is "STS2", "S1S2", or "all".
    error_bar : bool
        Whether to include the confidence interval as error bars in the plot,
        this is only applicable for the bar plot.
    file : str
        If provided, the generated figure will be saved as a png file.
    close_fig : bool
        Whether to close the figure
        (if not close, new figure will be overlaid on the current figure).
    
    Returns
    -------
    figure : :class:`matplotlib.figure.Figure`
        The generated figure.
    axis : :class:`matplotlib.axes._subplots.AxesSubplot`
        The generated figure axis.
        If generating bar plot and heat map, a tuple of two axes will be returned
        for the respective plot.
    
    '''
    kind_upper = kind.upper()
    if kind_upper=='ALL' or set(kind_upper)==set('STS1S2'):
        kind_upper = 'STS1S2'
    elif not set(kind_upper).union(set('STS1S2'))==set('STS1S2'):
        raise ValueError(f'The plot kind of "{kind}" is invalid.')
    
    ax_sts1 = ax_s2 = None
    
    param_names = _update_input(parameters, result_dct[metric.name]['ST'].index)    
    param_names = param_names if isinstance(param_names[0], str) \
                              else [p.name for p in param_names]
    
    st_df = result_dct[metric.name]['ST'].loc[[p for p in param_names]]
    s1_df = result_dct[metric.name]['S1'].loc[[p for p in param_names]]
    sts1_df = pd.concat((st_df, s1_df), axis=1).sort_values('ST', ascending=False)

    if kind_upper in ('ST', 'S1', 'STS1', 'S1ST'): # no S2, bar plot only
        ax_sts1 = _plot_bar(kind_upper, sts1_df, error_bar)
        return _save_fig_return(ax_sts1.figure, ax_sts1, file, close_fig)
    
    else: # has S2, need heat map
        s2_df = result_dct[metric.name]['S2']
        hmap_df = pd.DataFrame(columns=sts1_df.index, index=sts1_df.index)
        for (p1, p2) in s2_df.index:
            if not (p1 in hmap_df.index and p2 in hmap_df.index):
                continue
            hmap_df[p1][p2] = hmap_df[p2][p1] = s2_df.S2[(p1, p2)]

        if kind_upper == 'S2': # only S2, only heat map
            ax_s2 = _plot_heatmap(hmap_df, annot=annotate_heatmap)
            return _save_fig_return(ax_s2.figure, ax_s2, file, close_fig)
            
        else: # has S2, need heat map
            not_s2 = ''.join(i for i in kind_upper.split('S2')) # 'ST', 'S1', or 'STS1'

            plot_in_diagonal = plot_in_diagonal.upper()
            if not_s2 and not_s2 == plot_in_diagonal: # ST or S1 in heat map, only heat map
                ax_s2 = _plot_heatmap(hmap_df, ax=ax_s2, annot=annotate_heatmap,
                                      diagonal=plot_in_diagonal,
                                      sts1_df=sts1_df, default_cbar=True)
                return _save_fig_return(ax_s2.figure, ax_s2, file, close_fig)

            # ST or S1 in bar, need bar and heat map
            elif not_s2 != 'STS1':
                warn(f'The plot_in_diagonal value of "{plot_in_diagonal}" ' \
                     f'is invalid for kind "{kind}" and is ignored.')
                plot_in_diagonal = ''

            if plot_in_diagonal and plot_in_diagonal.upper() not in ('ST', 'S1'):
                raise ValueError('plot_in_diagonal must be "ST", "S1", or "", '\
                                 f'not "{plot_in_diagonal}".')
            
            fig, (ax_s2, ax_sts1) = plt.subplots(1, 2, figsize=(8, 5))
            bar = not_s2.replace(plot_in_diagonal, '') # 'ST', 'S1', or 'STS1'

            ax_sts1 = _plot_bar(bar, sts1_df, error_bar, ax=ax_sts1)
            ax_sts1.yaxis.set_visible(False)
            if bar == 'ST':
                ax_sts1_title = 'Total Effects'
            elif bar == 'S1':
                ax_sts1_title = 'Main Effects'
            else:
                ax_sts1_title = 'Total/Main Effects'
            ax_sts1.set_title(ax_sts1_title)
            
            ax_s2 = _plot_heatmap(hmap_df, ax=ax_s2, annot=annotate_heatmap,
                                        diagonal=plot_in_diagonal,
                                        sts1_df=sts1_df, default_cbar=False)
            
            labels = [i if len(i)<=15 else i[0:15]+'...' for i in hmap_df.index]
            ax_s2.yaxis.set_label_position('right')
            ax_s2.yaxis.tick_right()
            ax_s2.set_yticklabels(labels, rotation=0, ha='center',
                                  position=(1.2, 0.0))
            
            ax_s2.tick_params(length=0)
            xlabels = labels.copy()
            xlabels[0] = '' if not plot_in_diagonal else xlabels[0]
            ax_s2.set_xticklabels(xlabels)

            plt.tight_layout()
            plt.subplots_adjust(wspace=0.4, top=0.85)
            fig.suptitle(f'Variance breakdown for {metric.name.lower()}')

        return _save_fig_return(fig, (ax_sts1, ax_s2), file, close_fig)



