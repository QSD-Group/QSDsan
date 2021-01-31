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

# =============================================================================
# Sobol analysis with Saltelli sampling
# =============================================================================

__all__ = ('define_saltelli_problem', 'sobol_analysis')

import pandas as pd
# #!!! Can potentially change the sampling method
# from SALib.sample.saltelli import sample
from SALib.analyze.sobol import analyze

def define_saltelli_problem(model):    
    '''
    Define the problem to be used in Saltelli sampling for ``SALib``.
    
    Parameters
    ----------
    model : :class:`biosteam.Model`
        Uncertainty model with defined metrics and parameters.    

    See Also
    --------
    `SALib.sample.saltelli.sample <https://salib.readthedocs.io/en/latest/api.html#sobol-sensitivity-analysis>`_

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

def sobol_analysis(model, samples, problem, metrics, fill_nan_with_mean=False,
                   calc_second_order=True, conf_level=0.95, print_to_console=False,
                   file='', **sobol_kwargs):
    '''
    Run Sobol sensitivity analysis using ``SALib``.
    
    Parameters
    ----------
    model : :class:`biosteam.Model`
        Uncertainty model with defined metrics and parameters.
    samples : :class:`numpy.array`
        Samples for Sobol analysis.
    problem : dict
        Keys including "num_vars", "names", and "bounds".
    metrics : iterable
        Metrics to be included for Sobol analysis, must be a subset of
        (i.e., included in the `metrics` attribute of the given model).
    fill_nan_with_mean : bool
        Whether to use the mean value to replace the "NaN" values in the metric results.
    calc_second_order : bool
        Whether to calculate second-order interaction effects.
    conf_level : float
        Confidence level of results.
    print_to_console : float
        Whether to show results in the console.
    file : str
        If provided, the results will be saved as an Excel file.

    Returns
    -------
    si_dct : dict
        A dict of Sobol analysis results.    

    See Also
    --------
    `SALib.analyze.sobol.analyze <https://salib.readthedocs.io/en/latest/api.html#sobol-sensitivity-analysis>`_
    
    '''
    si_dct = {}
    for metric in metrics:
        results = model.table[metric.index]
        if results.isna().values.any():
            if fill_nan_with_mean:
                results.fillna(results.mean(), inplace=True)
            else:
                raise ValueError('"NaN" values in metrics, cannot run analysis.')
        si = analyze(problem, results.to_numpy(),
                     calc_second_order=calc_second_order,
                     conf_level=conf_level, print_to_console=print_to_console,
                     **sobol_kwargs)
        si_dct[metric.name] = si
    if file:
        writer = pd.ExcelWriter(file)
        for name, si in si_dct.items():
            n_row = 0
            for df in si.to_df():
                df.to_excel(writer, sheet_name=name, startrow=n_row)
                n_row += len(df.index) + 2 + len(df.columns.names)
        writer.save()
    return si_dct









