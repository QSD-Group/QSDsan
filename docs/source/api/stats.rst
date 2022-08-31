stats
=====

``QSDsan`` has many functions uncertainty analysis visualization and sensitivity analysis/visualization. You can find the complete script of following code snippets from `stats_demo.py <https://github.com/QSD-Group/EXPOsan/blob/main/exposan/bwaise/stats_demo.py>`_ in the ``bwaise`` module of Exposan.

Uncertainties
-------------
.. automethod:: qsdsan.stats.plot_uncertainties

Examples
^^^^^^^^

Box plot
********
.. code:: python

	import pandas as pd
	from qsdsan import stats as s
	from exposan import bwaise as bw

	m = bw.models
	modelA = bw.modelA
		
	# Total COD/N/P/K recovery and net cost/GWP
	modelA.metrics = key_metrics = bw.get_key_metrics(
	    modelA, alt_names={'Annual net cost': 'Cost',
	                       'Net emission GlobalWarming': 'GWP'})
		
	seed = 3221 # set numpy seed for sample reproducibility
		
	# Run Monte Carlo uncertainty analysis and get Spearman rank correlations,
	# here we use a small sample size for demonstrative purpose
	m.run_uncertainty(modelA, N=100, seed=seed, rule='L',
	                  percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1))
		
	# Pass a path to `file` or use `fig.savefig` if want to save the figure,
	# the `file` kwarg exists for pretty much all of the plotting functions
	fig, ax = s.plot_uncertainties(modelA,
	                               x_axis=key_metrics[:-2], # only recoveries
	                               kind='box', file='')
	   
	# Trim figure
	fig.subplots_adjust(bottom=0.25)
	for label in ax.get_xticklabels():
	    label.set_rotation(45)

.. figure:: ../images/stats/plot_uncer_box.png
   :width: 50%


Histogram plot
**************
.. code:: python
	
	# Kernel density curve can be added to the histogram,
	# with a log scale, we can have all metric results in the same plot
	fig, ax = s.plot_uncertainties(modelA, y_axis=key_metrics, kind='hist',
	                               center_kws={'kde':True, 'log_scale': 10})

.. figure:: ../images/stats/plot_uncer_hist.png
   :width: 60%


.. code:: python
	
	# We can also have 2D histogram plot
	fig, axes = s.plot_uncertainties(modelA,
	                                 x_axis=key_metrics[-2], # cost
	                                 y_axis=key_metrics[-1], # GWP
	                                 kind='hist-box')

.. figure:: ../images/stats/plot_uncer_hist-box.png
   :width: 50%


Kernel density plots
********************
.. code:: python
	
	# Similar to histogram plots, kernel density plots can be 1D
	fig, ax = s.plot_uncertainties(modelA, x_axis=key_metrics, kind='kde',
	                               center_kws={'fill': True, 'log_scale': 2})                        

.. figure:: ../images/stats/plot_uncer_kde.png
   :width: 60%


.. code:: python
	
	# Or 2D with different kinds of margins
	fig, axes = s.plot_uncertainties(modelA, x_axis=key_metrics[-2],
	                                 y_axis=key_metrics[-1], kind='kde-kde',
	                                 center_kws={'fill': True})

.. figure:: ../images/stats/plot_uncer_kde-kde.png
   :width: 50%


.. code:: python
	
	fig, axes = s.plot_uncertainties(modelA, x_axis=key_metrics[-2],
	                                 y_axis=key_metrics[-1], kind='kde-hist',
	                                 center_kws={'fill': True},
	                                 margin_kws={'kde': True, 'fill': False})

.. figure:: ../images/stats/plot_uncer_kde-hist.png
   :width: 50%


Correlations
------------
.. automethod:: qsdsan.stats.get_correlations
.. automethod:: qsdsan.stats.plot_correlations


Examples
^^^^^^^^

Bar plot for single metric
**************************
.. code:: python

	spearman_rho, spearman_p = s.get_correlations(
	    modelA, kind='Spearman', nan_policy='raise',
	    file='') # pass a path to `file` if you want to save the results as an Excel

	# Filter out parameters that only meet a certain threshold
	def filter_parameters(model, df, threshold):
	    new_df = pd.concat((df[df>=threshold], df[df<=-threshold]))
	    filtered = new_df.dropna(how='all')
	    param_dct = {p.name_with_units:p for p in model.get_parameters()}
	    parameters = set(param_dct[i[1]] for i in filtered.index)
	    return list(parameters)

	# Only want parameters with Spearman's rho >= 0.4 or <= -0.4
	modelA.parameters = key_parameters = \
	    filter_parameters(modelA, spearman_rho, threshold=0.4)

	fig, ax = s.plot_correlations(spearman_rho, parameters=key_parameters,
		                          metrics=key_metrics[-2])

	fig.subplots_adjust(left=0.25)


.. figure:: ../images/stats/plot_corr_bar.png
   :width: 60%


Bubble plot for multiple metrics
********************************
.. code:: python

	fig, ax = s.plot_correlations(
	    spearman_rho, parameters=key_parameters, metrics=key_metrics)


.. figure:: ../images/stats/plot_corr_bubble.png
   :width: 80%


Input and sample preparation
----------------------------
.. automethod:: qsdsan.stats.define_inputs
.. automethod:: qsdsan.stats.generate_samples


Morris
------
.. automethod:: qsdsan.stats.morris_analysis
.. automethod:: qsdsan.stats.morris_till_convergence
.. automethod:: qsdsan.stats.plot_morris_results
.. automethod:: qsdsan.stats.plot_morris_convergence

Examples
^^^^^^^^

:math:`\sigma` vs. :math:`\mu^*`
********************************
.. code:: python
	
	# Run Morris analysis without testing the convergence,
	# here we use a small sample size for demonstrative purpose
	inputs = s.define_inputs(modelA)
	morris_samples = s.generate_samples(inputs, kind='Morris', N=10, seed=seed)

	evaluate = bw.evaluate
	evaluate(modelA, morris_samples)

	dct = s.morris_analysis(modelA, inputs, metrics=key_metrics, seed=seed,
	                        nan_policy='fill_mean')

	# Unfortunately the auto-labelling is not good when you have close points,
	# so you'll have to do some manual manipulation
	fig, ax = s.plot_morris_results(dct, key_metrics[-2])

	fig.subplots_adjust(bottom=0.3)


.. figure:: ../images/stats/plot_morris.png
   :width: 60%


Line plot with error bands for evolutionary of :math:`\mu^*`
************************************************************
.. code:: python
	
	# Test if mu_star can converge within 100 trajectories
	# (spoiler: it cannot because we already sort of selected the key parameters,
	# and you will get a message prompt)
	dct = s.morris_till_convergence(modelA, inputs, metrics=key_metrics, seed=seed,
	                                N_max=100)

	# Look at mu_star values for two parameters with regard to cost
	fig, ax = s.plot_morris_convergence(dct,
	                                    parameters=key_parameters[:2],
	                                    metric=key_metrics[-2], plot_rank=False)


.. figure:: ../images/stats/plot_morris_conv.png
   :width: 80%


Line plot for evolutionary of :math:`\mu^*` rank
************************************************
.. code:: python
	
	# Look at ranks of mu_star values for all parameters with regard to cost
	fig, ax = s.plot_morris_convergence(dct, parameters=key_parameters,
	                                    metric=key_metrics[-2], plot_rank=True)


.. figure:: ../images/stats/plot_morris_conv_rank.png
   :width: 80%


FAST
----
.. automethod:: qsdsan.stats.fast_analysis
.. automethod:: qsdsan.stats.plot_fast_results

Examples
^^^^^^^^

Bar plot for (e)FAST
********************
.. code:: python
	
	# Total and main effects from FAST analysis,
	# here we use a small sample size for demonstrative purpose
	fast_samples = s.generate_samples(inputs, kind='FAST', N=100, M=4, seed=seed)

	evaluate(modelA, fast_samples)

	dct = s.fast_analysis(modelA, inputs, kind='FAST', metrics=key_metrics,
	                      M=4, seed=seed, nan_policy='fill_mean')

	fig, ax = s.plot_fast_results(dct, metric=key_metrics[-2])

	fig.subplots_adjust(left=0.4)


.. figure:: ../images/stats/plot_fast.png
   :width: 60%


Bar plot for RBD-FAST
*********************
.. code:: python
	
	# Main effects from RBD-FAST analysis,
	# here we use a small sample size for demonstrative purpose
	fast_samples = s.generate_samples(inputs, kind='RBD', N=100, seed=seed)

	evaluate(modelA, fast_samples)

	dct = s.fast_analysis(modelA, inputs, kind='RBD', metrics=key_metrics,
	                      seed=seed, nan_policy='fill_mean')

	fig, ax = s.plot_fast_results(dct, metric=key_metrics[-2])

	fig.subplots_adjust(left=0.4)


.. figure:: ../images/stats/plot_rbd.png
   :width: 60%


Sobol
-----
.. automethod:: qsdsan.stats.sobol_analysis
.. automethod:: qsdsan.stats.plot_sobol_results

Examples
^^^^^^^^

Bar plot for total and main effects
***********************************
.. code:: python
	
	# Run Sobol analysis, here we use a small sample size for demonstrative purpose
	sobol_samples = s.generate_samples(inputs, kind='Sobol', N=10,
	                                   calc_second_order=True)

	evaluate(modelA, sobol_samples)

	dct = s.sobol_analysis(modelA, inputs, metrics=key_metrics, seed=seed,
	                       calc_second_order=True, conf_level=0.95,
	                       nan_policy='fill_mean')

	fig, ax = s.plot_sobol_results(dct, metric=key_metrics[-1], kind='STS1')

	fig.subplots_adjust(left=0.4, top=0.95)


.. figure:: ../images/stats/plot_sobol_sts1.png
   :width: 60%


Heat map for total and second-order effects
*******************************************
.. code:: python
	
	fig, ax = s.plot_sobol_results(dct, metric=key_metrics[-1], kind='STS2',
	                               plot_in_diagonal='ST')

	for label in ax.get_xticklabels():
	    label.set_rotation(45)

	fig.subplots_adjust(left=0.4, bottom=0.4)


.. figure:: ../images/stats/plot_sobol_sts2.png
   :width: 80%


Bar plot and heat map for total, main, and second-order effects
***************************************************************
.. code:: python
	
	fig, ax = s.plot_sobol_results(dct, metric=key_metrics[-1], kind='all')


.. figure:: ../images/stats/plot_sobol_all.png