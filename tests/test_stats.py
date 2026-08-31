#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt
for license details.
'''

import matplotlib
matplotlib.use('Agg')

import pytest
from qsdsan import stats as s
from qsdsan.utils import create_example_model

__all__ = (
    'test_correlations',
    'test_define_inputs_and_sampling',
    'test_morris_analysis',
    'test_morris_till_convergence',
    'test_fast_analysis',
    'test_sobol_analysis',
    'test_plot_uncertainties',
    'test_plot_correlations',
    'test_plot_morris_results',
    'test_plot_fast_results',
    'test_plot_sobol_results',
    'test_plot_api_harmonization',
)


# =============================================================================
# Module-scoped fixtures — each analysis runs exactly once
# =============================================================================

@pytest.fixture(scope='module')
def evaluated_model():
    return create_example_model(evaluate=True, N=100, rule='L', seed=554)


@pytest.fixture(scope='module')
def inputs(evaluated_model):
    return s.define_inputs(evaluated_model)


@pytest.fixture(scope='module')
def morris_analyzed(evaluated_model, inputs):
    samples = s.generate_samples(inputs, kind='Morris', N=10, seed=554)
    m = evaluated_model.copy()
    m.load_samples(samples)
    m.evaluate()
    dct = s.morris_analysis(m, inputs, nan_policy='fill_mean')
    return m, dct


@pytest.fixture(scope='module')
def morris_convergence(evaluated_model, inputs):
    # morris_till_convergence mutates the model it receives, so pass a copy
    dct = s.morris_till_convergence(evaluated_model.copy(), inputs, seed=554, N_max=5)
    return dct


@pytest.fixture(scope='module')
def fast_analyzed(evaluated_model, inputs):
    samples = s.generate_samples(inputs, kind='FAST', N=100, seed=554)
    m = evaluated_model.copy()
    m.load_samples(samples)
    m.evaluate()
    return m, s.fast_analysis(m, inputs, kind='FAST', nan_policy='fill_mean')


@pytest.fixture(scope='module')
def rbd_analyzed(evaluated_model, inputs):
    samples = s.generate_samples(inputs, kind='RBD', N=100, seed=554)
    m = evaluated_model.copy()
    m.load_samples(samples)
    m.evaluate()
    return m, s.fast_analysis(m, inputs, kind='RBD', nan_policy='fill_mean')


@pytest.fixture(scope='module')
def sobol_analyzed(evaluated_model, inputs):
    samples = s.generate_samples(inputs, kind='Sobol', N=10, calc_second_order=True)
    m = evaluated_model.copy()
    m.load_samples(samples)
    m.evaluate()
    return m, s.sobol_analysis(m, inputs, calc_second_order=True, nan_policy='fill_mean')


# =============================================================================
# Correlations
# =============================================================================

def test_correlations(evaluated_model):
    model = evaluated_model
    p, m = model.parameters, model.metrics

    r_sp, _ = s.get_correlations(model, kind='Spearman')
    assert r_sp.shape[0] >= 1

    s.get_correlations(model, kind='Pearson')
    s.get_correlations(model, kind='Kendall')
    s.get_correlations(model, input_x=p[0], input_y=m[-1], kind='KS', thresholds=[0.1])

    r_sub, _ = s.get_correlations(model, input_x=p[0], input_y=m[-1], kind='Spearman')
    assert r_sub.shape[0] >= 1

    with pytest.raises(ValueError, match='nan_policy'):
        s.get_correlations(model, nan_policy='drop')
    with pytest.raises(ValueError, match='kind'):
        s.get_correlations(model, kind='InvalidKind')


# =============================================================================
# SALib inputs and sampling
# =============================================================================

def test_define_inputs_and_sampling(inputs, evaluated_model):
    assert 'num_vars' in inputs
    assert inputs['num_vars'] == len(evaluated_model.parameters)

    morris_samples = s.generate_samples(inputs, kind='Morris', N=5, seed=554)
    assert morris_samples.shape[1] == inputs['num_vars']

    fast_samples = s.generate_samples(inputs, kind='FAST', N=100, seed=554)
    assert fast_samples.shape[1] == inputs['num_vars']

    rbd_samples = s.generate_samples(inputs, kind='RBD', N=50, seed=554)
    assert rbd_samples.shape[1] == inputs['num_vars']

    sobol_samples = s.generate_samples(inputs, kind='Sobol', N=10, calc_second_order=True)
    assert sobol_samples.shape[1] == inputs['num_vars']

    with pytest.raises(ValueError, match='seed'):
        s.generate_samples(inputs, kind='Sobol', N=10, seed=554)

    with pytest.raises(ValueError, match='kind'):
        s.generate_samples(inputs, kind='Unknown', N=10)


# =============================================================================
# Morris
# =============================================================================

def test_morris_analysis(morris_analyzed, evaluated_model):
    m, dct = morris_analyzed
    assert set(dct.keys()) == {metric.name for metric in m.metrics}

    m0 = m.metrics[0]
    dct_single = s.morris_analysis(m, s.define_inputs(evaluated_model),
                                   metrics=m0, nan_policy='fill_mean')
    assert m0.name in dct_single


def test_morris_till_convergence(morris_convergence):
    assert 'mu_star' in morris_convergence
    assert 'mu_star_conf' in morris_convergence


# =============================================================================
# FAST / RBD
# =============================================================================

def test_fast_analysis(fast_analyzed, rbd_analyzed, evaluated_model):
    m_fast, dct_fast = fast_analyzed
    assert set(dct_fast.keys()) == {metric.name for metric in m_fast.metrics}

    m_rbd, dct_rbd = rbd_analyzed
    assert set(dct_rbd.keys()) == {metric.name for metric in m_rbd.metrics}

    with pytest.raises(ValueError, match='kind'):
        s.fast_analysis(m_fast, s.define_inputs(evaluated_model), kind='Invalid')


# =============================================================================
# Sobol
# =============================================================================

def test_sobol_analysis(sobol_analyzed):
    m, dct = sobol_analyzed
    assert set(dct.keys()) == {metric.name for metric in m.metrics}
    for metric_name in dct:
        assert 'ST' in dct[metric_name]
        assert 'S1' in dct[metric_name]


# =============================================================================
# Plots — uncertainties
# =============================================================================

def test_plot_uncertainties(evaluated_model):
    m = evaluated_model.metrics

    fig, ax = s.plot_uncertainties(evaluated_model, x_axis=m[0], kind='box', close_fig=True)
    assert fig is not None
    s.plot_uncertainties(evaluated_model, y_axis=m[1], kind='hist', close_fig=True)
    s.plot_uncertainties(evaluated_model, x_axis=m[2], kind='kde', close_fig=True)
    s.plot_uncertainties(evaluated_model, x_axis=m[0], y_axis=m[1],
                         kind='hist-hist', close_fig=True)
    s.plot_uncertainties(evaluated_model, x_axis=m[0], y_axis=m[1],
                         kind='kde-kde', close_fig=True)


# =============================================================================
# Plots — correlations
# =============================================================================

def test_plot_correlations(evaluated_model):
    m = evaluated_model.metrics
    r_df, _ = s.get_correlations(evaluated_model, kind='Spearman')

    fig, ax = s.plot_correlations(r_df, metrics=m[-2], close_fig=True)
    assert fig is not None
    s.plot_correlations(r_df, close_fig=True)


# =============================================================================
# Plots — Morris
# =============================================================================

def test_plot_morris_results(morris_analyzed, morris_convergence, evaluated_model):
    _, dct = morris_analyzed
    m0 = evaluated_model.metrics[0]

    fig, ax = s.plot_morris_results(dct, metric=m0, close_fig=True)
    assert fig is not None

    fig, ax = s.plot_morris_convergence(
        morris_convergence, metric=evaluated_model.metrics[-2],
        plot_rank=False, close_fig=True)
    s.plot_morris_convergence(
        morris_convergence, metric=evaluated_model.metrics[-2],
        plot_rank=True, close_fig=True)


# =============================================================================
# Plots — FAST
# =============================================================================

def test_plot_fast_results(fast_analyzed, evaluated_model):
    _, dct = fast_analyzed
    fig, ax = s.plot_fast_results(dct, metric=evaluated_model.metrics[-3], close_fig=True)
    assert fig is not None


# =============================================================================
# Plots — Sobol
# =============================================================================

def test_plot_sobol_results(sobol_analyzed, evaluated_model):
    _, dct = sobol_analyzed
    m = evaluated_model.metrics

    result = s.plot_sobol_results(dct, metric=m[-1], kind='STS1', close_fig=True)
    assert result is not None
    s.plot_sobol_results(dct, metric=m[-1], kind='STS2',
                         plot_in_diagonal='ST', close_fig=True)
    s.plot_sobol_results(dct, metric=m[0], kind='all', close_fig=True)


# =============================================================================
# Plot API harmonization — ax= is honored where applicable, **plot_kws reaches
# the underlying call. Smoke tests; numeric output is not asserted.
# =============================================================================

def _to_rgba(color):
    from matplotlib.colors import to_rgba
    return to_rgba(color)


def test_plot_api_harmonization(evaluated_model, morris_analyzed,
                                morris_convergence, fast_analyzed, sobol_analyzed):
    """All six plot_* helpers accept ax= and **plot_kws uniformly."""
    import matplotlib.pyplot as plt

    m = evaluated_model.metrics
    red = _to_rgba('red')

    # plot_uncertainties: 1D honors ax= and plot_kws reaches sns.histplot
    fig, host_ax = plt.subplots()
    fig_u, ax_u = s.plot_uncertainties(
        evaluated_model, x_axis=m[0], kind='hist', ax=host_ax,
        color='red', close_fig=False,
    )
    assert ax_u is host_ax
    # seaborn applies alpha to bar facecolor; compare only the RGB triple.
    red_rgb = red[:3]
    assert any(_to_rgba(p.get_facecolor())[:3] == red_rgb for p in ax_u.patches)
    plt.close(fig_u)

    # plot_uncertainties: 2D ignores ax= (jointplot builds its own grid) but
    # accepts plot_kws and merges it into center_kws.
    _, throwaway = plt.subplots()
    fig2, _ = s.plot_uncertainties(
        evaluated_model, x_axis=m[0], y_axis=m[1], kind='kde-kde',
        ax=throwaway, color='red', close_fig=False,
    )
    assert fig2 is not throwaway.figure  # built its own figure
    plt.close(fig2)
    plt.close(throwaway.figure)

    # plot_correlations: ax= is accepted (currently ignored) and **plot_kws
    # is forwarded. Smoke check: tornado path runs with a forwarded kwarg.
    r_df, _ = s.get_correlations(evaluated_model, kind='Spearman')
    fig_c, _ = s.plot_correlations(r_df, metrics=m[-2], ax=None, close_fig=False)
    assert fig_c is not None
    plt.close(fig_c)

    # plot_morris_results: ax= honored, plot_kws reaches ax.scatter
    # Use x_axis='mu' to exercise the pure scatter path (mu_star uses errorbar
    # which accepts a different kwarg set).
    _, morris_dct = morris_analyzed
    fig, host_ax = plt.subplots()
    fig_m, ax_m = s.plot_morris_results(
        morris_dct, metric=m[0], ax=host_ax, kind='scatter', x_axis='mu',
        s=80, close_fig=False,
    )
    assert ax_m is host_ax
    sizes = ax_m.collections[0].get_sizes()
    assert len(sizes) > 0 and sizes[0] == 80
    plt.close(fig_m)

    # plot_morris_results: bar branch — plot_kws goes via opts= to SALib;
    # smoke test that it does not raise with an empty plot_kws dict.
    fig, host_ax = plt.subplots()
    s.plot_morris_results(
        morris_dct, metric=m[0], ax=host_ax, kind='bar', close_fig=False,
    )
    plt.close(host_ax.figure)

    # plot_morris_convergence: ax= honored, plot_kws reaches ax.plot
    fig, host_ax = plt.subplots()
    fig_mc, ax_mc = s.plot_morris_convergence(
        morris_convergence, metric=m[-2], ax=host_ax,
        linestyle=':', close_fig=False,
    )
    assert ax_mc is host_ax
    assert any(line.get_linestyle() == ':' for line in ax_mc.get_lines())
    plt.close(fig_mc)

    # plot_fast_results: ax= honored, plot_kws reaches sns.barplot
    _, fast_dct = fast_analyzed
    fig, host_ax = plt.subplots()
    fig_f, ax_f = s.plot_fast_results(
        fast_dct, metric=m[-3], ax=host_ax, alpha=0.4, close_fig=False,
    )
    assert ax_f is host_ax
    assert any(abs(p.get_alpha() - 0.4) < 1e-9
               for p in ax_f.patches if p.get_alpha() is not None)
    plt.close(fig_f)

    # plot_sobol_results: ax= honored for bar-only kinds, plot_kws reaches
    # sns.barplot.
    _, sobol_dct = sobol_analyzed
    fig, host_ax = plt.subplots()
    fig_s, ax_s = s.plot_sobol_results(
        sobol_dct, metric=m[-1], kind='STS1', ax=host_ax,
        alpha=0.3, close_fig=False,
    )
    assert ax_s is host_ax
    assert any(abs(p.get_alpha() - 0.3) < 1e-9
               for p in ax_s.patches if p.get_alpha() is not None)
    plt.close(fig_s)


if __name__ == '__main__':
    m = create_example_model(evaluate=True, N=100, rule='L', seed=554)
    inp = s.define_inputs(m)

    morris_samples = s.generate_samples(inp, kind='Morris', N=10, seed=554)
    m_morris = m.copy(); m_morris.load_samples(morris_samples); m_morris.evaluate()
    morris_dct = s.morris_analysis(m_morris, inp, nan_policy='fill_mean')

    morris_conv = s.morris_till_convergence(m.copy(), inp, seed=554, N_max=5)

    fast_samples = s.generate_samples(inp, kind='FAST', N=100, seed=554)
    m_fast = m.copy(); m_fast.load_samples(fast_samples); m_fast.evaluate()
    fast_dct = s.fast_analysis(m_fast, inp, kind='FAST', nan_policy='fill_mean')

    rbd_samples = s.generate_samples(inp, kind='RBD', N=100, seed=554)
    m_rbd = m.copy(); m_rbd.load_samples(rbd_samples); m_rbd.evaluate()
    rbd_dct = s.fast_analysis(m_rbd, inp, kind='RBD', nan_policy='fill_mean')

    sobol_samples = s.generate_samples(inp, kind='Sobol', N=10, calc_second_order=True)
    m_sobol = m.copy(); m_sobol.load_samples(sobol_samples); m_sobol.evaluate()
    sobol_dct = s.sobol_analysis(m_sobol, inp, calc_second_order=True, nan_policy='fill_mean')

    test_correlations(m)
    test_define_inputs_and_sampling(inp, m)
    test_morris_analysis((m_morris, morris_dct), m)
    test_morris_till_convergence(morris_conv)
    test_fast_analysis((m_fast, fast_dct), (m_rbd, rbd_dct), m)
    test_sobol_analysis((m_sobol, sobol_dct))
    test_plot_uncertainties(m)
    test_plot_correlations(m)
    test_plot_morris_results((m_morris, morris_dct), morris_conv, m)
    test_plot_fast_results((m_fast, fast_dct), m)
    test_plot_sobol_results((m_sobol, sobol_dct), m)
