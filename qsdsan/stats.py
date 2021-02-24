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
           'morris_analysis', 'plot_morris_results',
           'morris_till_convergence', 'plot_morris_convergence',
           'sobol_analysis', 'plot_sobol_results',)

import numpy as np
import pandas as pd
import biosteam as bst
from scipy.stats import pearsonr, spearmanr, kendalltau
from SALib.sample import (morris as morris_sampling, saltelli)
from SALib.analyze import morris, sobol
from SALib.plotting import morris as sa_plt_morris
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
                    kind='Pearson', nan_policy='propagate', file='',
                    **kwargs):
    '''
    Get correlation coefficients between two inputs using ``scipy``.
    
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
        Can be "Pearson" for Pearson's r, "Spearman" for Spearman's rho,
        or "Kendall" for Kendall's tau.
    nan_policy : str
        - "propagate": returns nan.
        - "raise": raise an error.
        - "omit": drop the pair from analysis.
    file : str
        If provided, the results will be saved as an Excel file.

    Returns
    -------
    Two :class:`pandas.DataFrame` containing Pearson'r, Spearman's rho, or Kendall's tau
    and p-values.
    
    See Also
    --------
    `scipy.stats.pearsonr <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html>`_
    
    `scipy.stats.spearmanr <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html>`_
    
    `scipy.stats.kendalltau <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kendalltau.html>`_
    
    '''
    table = model.table.astype('float64')
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
                    r, p = pearsonr(df.iloc[:,0], df.iloc[:,1], **kwargs)
                    sheet_name = 'r'
                elif kind.capitalize() == 'Spearman':
                    r, p = spearmanr(df.iloc[:,0], df.iloc[:,1], **kwargs)
                    sheet_name = 'rho'
                elif kind.capitalize() == 'Kendall':
                    r, p = kendalltau(df.iloc[:,0], df.iloc[:,1], **kwargs)
                    sheet_name = 'tau'
                else:
                    raise ValueError('kind can only be "Pearson", "Spearman", ' \
                                      f'or "Kendall", not "{kind}".')
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
    `SALib basics <https://salib.readthedocs.io/en/latest/basics.html#an-example>`_

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
        A dict generated by :func:`qsdsan.sensitivity.define_inputs` to be used for ``SALib``,
        keys should include "num_vars", "names", and "bounds".
    kind : str
        Can be "Morris" (for Morris analysis) or "Saltelli" (for Sobol analysis).
    N : int
        The number of trajectories (Morris) or samples.
    seed : int
        Seed to generate random samples.
    
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
    
    Returns
    -------
    morris_dct : dict
        A dict of Morris analysis results.
    
    See Also
    --------
    `SALib.analyze.morris <https://salib.readthedocs.io/en/latest/api.html?highlight=SALib.analyze.morris.analyze#method-of-morris>`_
    
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


def plot_morris_results(morris_dct, metric,
                        axis=None, kind='scatter',
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
        Results dict generated by :func:`morris_analysis`
    metric : :class:`biosteam.Metric`
        The metric of interest for the plot.
    axis : :class:`matplotlib.axes._subplots.AxesSubplot`
        Axis for the figure, a new axis will be created if not provided.
    kind : str
        Either "scatter" (:math:`{\mu^*}` vs. :math:`{\sigma}`) or "bar" (:math:`{\mu^*}` with confidence interval) plot.
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
        Whether to close the figure (if not close, new figure will be overlaid on the current figure).
        
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
    ax = axis if axis else plt.subplot()
    if kind == 'scatter':
        ax.scatter(x_data, y_data)
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

    if file:
        fig.savefig(file, dpi=300)
    if close_fig:
        plt.close()
    return fig, ax


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

    See Also
    --------
    :func:`qsdsan.stats.generate_samples`
    
    :func:`qsdsan.stats.morris_analysis`
    
    '''
    num_levels = kwargs['num_levels'] if 'num_levels' in kwargs.keys() else 4
    kwargs = {i:kwargs[i] for i in kwargs.keys() if i!='num_levels'}
    samples = generate_samples(inputs=inputs, kind='Morris', N=N_max, seed=seed, num_levels=num_levels)
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
    
    for n in range(3, N_max):
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
            cum_dct['mu_star'][m.name].loc[n] = mu_star
            cum_dct['mu_star_conf'][m.name].loc[n] = mu_star_conf
            
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


def plot_morris_convergence(result_dct, metric, axis=None, parameters=(),
                            error_bar=True, file='', close_fig=True):
    '''
    Plot the evolution of :math:`{\mu^*}` with the number of trajectories.
    
    Parameters
    ----------
    result_dct : dict
        Result dictionary generated from :func:`qsdsan.stats.morris_till_convergence`
    metric : :class:`biosteam.Metric`
        The metric of interest for the plot.
    axis : :class:`matplotlib.axes._subplots.AxesSubplot`
        Axis for the figure, a new axis will be created if not provided.
    parameters : :class:`biosteam.Parameter`
        Single or a collection of model parameters whose :math:`{\mu^*}` will be
        included in the plot.
        Will be set to all parameters in retult_dct will be used if not provided.
    error_bar : bool
        Whether to include the confidence interval as error bars in the plot.
    file : str
        If provided, the generated figure will be saved as a png file.
    close_fig : bool
        Whether to close the figure (if not close, new figure will be overlaid on the current figure).
    
    Returns
    -------
    figure : :class:`matplotlib.figure.Figure`
        The generated figure.
    axis : :class:`matplotlib.axes._subplots.AxesSubplot`
        The generated figure axis.
    
    '''
    ax = axis if axis else plt.subplot()
    df = result_dct['mu_star'][metric.name]
    conf_df = result_dct['mu_star_conf'][metric.name]

    try: iter(parameters)
    except: parameters = (parameters,)
    if (parameters and not isinstance(parameters[0], str)):
        param_names = [p.name for p in parameters]
    else:
        param_names = df.columns

    for param in param_names:
        scatter = ax.scatter(df.index, df[param])
        if error_bar:
            ax.errorbar(df.index, df[param], conf_df[param])
    legend = ax.legend(scatter.legend_elements()[0], labels=param_names, loc='best')
    ax.add_artist(legend)
    ax.set_xlabel('Number of trajectories')
    ax.set_ylabel(f'$\mu^*$ for {metric.name.lower()}')
    
    fig = ax.figure
    if file:
        fig.savefig(file, dpi=300)
    if close_fig:
        plt.close()
    return fig, ax
    

# %%

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



def plot_sobol_results():
    '''NOT YET IMPLEMENTED'''





