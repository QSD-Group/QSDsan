Change Log
==========

This document records notable changes to `~ <https://github.com/QSD-Group/QSDsan>`_. We aim to follow `Semantic Versioning <https://semver.org/>`_.


Ongoing
-------
- Added an :class:`Equipment` class for design and costing of unit equipment.
- For the ``stats`` module:

	- More statistical tests:

		- :func:`qsdsan.stats.fast_analysis` for (extended) Fourier amplitude sensitivity test (FAST) and random balance design (RBD) FAST.
		- :func:`qsdsan.stats.morris_till_convergence` to run Morris analysis until the results converge.
		- Added Kendall's tau and Kolmogorov–Smirnov test to :func:`qsdsan.stats.get_correlations`.
	
	- Plotting functions to visualize all test results:

		- :func:`qsdsan.stats.plot_uncertainties` for results from uncertainty analysis.
		- :func:`qsdsan.stats.plot_correlations` for results from :func:`qsdsan.stats.get_correlation`.
		- Bar plot option for :func:`qsdsan.stats.plot_morris_results`.
		- :func:`qsdsan.stats.plot_morris_convergence` to plot :math:`{\mu^*}` against the number of trajectories.
		- :func:`qsdsan.stats.plot_fast_results` for results from FAST and/or RBD-FAST analyses.
		- :func:`qsdsan.stats.plot_sobol_results` for results from Sobol analysis.


`0.1.0`_ (2021-02-14)
---------------------
- Added a ``stats`` module including:

	- Pearson and Spearman correlations: :func:`qsdsan.stats.get_correlations`.
	- Morris One-at-A-Time (OAT) screening method: :func:`qsdsan.stats.morris_analysis`.

		- Also added a function for plotting: :func:`qsdsan.stats.plot_morris_results`.

	- Sobol sensitivity analysis: :func:`qsdsan.stats.sobol_analysis`.

- Added all uncertainty parameters for all of the scenarios in the bwaise system, also added demonstrative Morris and Sobol analysis.
- :func:`LCA.get_normalized_impacts` was replaced by :func:`qsdsan.LCA.get_allocated_impacts` for :class:`~.LCA` to enable flexible allocation options.
- Reformatted all documents, added instructions on documentation.
- Added brief instructions on contributing and code of conduct.
- Updated UML diagram.


`0.0.3`_ (2021-01-10)
---------------------
- More flexible setting of :class:`~.ImpactItem` for :class:`~.WasteStream`.
- Add status badge to README.rst
- Add CHANGELOG.rst
- Tutorial updates:

	- New:
		- :class:`~.TEA` and :class:`~.LCA`
	- Updated:
		-  :class:`~.Component` and :class:`~.WasteStream`
		-  :class:`~.SanUnit` and :class:`~.System`


`0.0.2`_ (2021-01-07)
---------------------
- Added the all three sanitation scenarios as described in `Trimmer et al.`_, including uncertainty/sensitivity analyses with tutorial.
- Inclusion of GPX models for estimation of :class:`~.WasteStream` properties.
- Live documentation for the `stable package`_ and `beta version`_.
- New classes:

    - All units in `Trimmer et al.`_
    - Added descriptors (``qsdsan.utils.descriptors``) and decorators (``qsdsan.utils.checkers``) to check user-input values.
    - :class:`~.utils.setters.AttrSetter`, :class:`~.utils.setters.DictAttrSetter`, and :class:`~.utils.getters.FuncGetter` for batch-setting of uncertainty analysis parameters.

- Added :func:`save_report` function to :class:`~.LCA` for report exporting.


`0.0.1`_ (2020-12-23)
---------------------
- First public release.


.. Other links
.. _stable package: https://qsdsan.readthedocs.io/en/latest/
.. _beta version: https://qsdsan-beta.readthedocs.io/en/latest/
.. _Trimmer et al.: https://doi.org/10.1021/acs.est.0c03296

.. Commit links
.. _0.1.0: https://github.com/yalinli2/QSDsan/commit/a3164b257d95889305aa94186bb072ad3d7b5f77
.. _0.0.3: https://github.com/QSD-Group/QSDsan/commit/e20222caccc58d9ee414ca08d8ec55f3a44ffca7
.. _0.0.2: https://github.com/QSD-Group/QSDsan/commit/84653f5979fbcd76a80ffb6b22ffec1c5ca2a084
.. _0.0.1: https://github.com/yalinli2/QSDsan/commit/f95e6172780cfe24ab68cd27ba19837e010b3d99

