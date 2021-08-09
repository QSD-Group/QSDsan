Change Log
==========

This document records notable changes to `QSDsan <https://github.com/QSD-Group/QSDsan>`_. We aim to follow `Semantic Versioning <https://semver.org/>`_.


`0.3.0`_
--------
- Now LCA data can be imported from external databases using the newly made `BW2QSD <https://github.com/QSD-Group/BW2QSD>`_ package.
- New subclasses of :class:`~.SanUnit`:

	- :class:`~.sanunits.Clarifier`
	- :class:`~.sanunits.CSTR`

	- :class:`~.sanunits.ElectrochemicalCell` using the following :class:`~.Equipment`:

		- :class:`~.equipments.Column`
		- :class:`~.equipments.Electrode`
		- :class:`~.equipments.Machine`
		- :class:`~.equipments.Membrane`

- New subclasses of :class:`~.Process`:

	- :class:`~.processes.DiffusedAeration`
	- :class:`~.processes.ASM1`
	- :class:`~.processes.ASM2d`

- Updated :class:`~.SanUnit` so that it can be initialized with any of :class:`thermosteam.Stream`, :class:`~.SanStream`, or :class:`~.WasteStream`.

	- These three classes can now be mixed.

- Added :class:`~.SanStream` for non-waste streams (e.g., gases).
- Updated the ``add_OPEX`` attribute of :class:`~.SanUnit` and ``system_add_OPEX`` attribute of :class:`~.SimpleTEA` so that they take :class:`dict` as the default to allow display of multiple additional operating expenses.
- Split the ``systems`` module into an individual package `EXPOsan <https://github.com/QSD-Group/exposan>`_.
- Now using :class:`thermosteam.utils.Registry` to manage :class:`~.ImpactIndicator` and :class:`~.ImpactItem`.
- Added `AppVeyor CI <https://ci.appveyor.com/project/yalinli2/qsdsan>`_.
- Renamed the ``master`` branch to ``main``.


`0.2.0`_
--------
- Added :class:`~.Process`, :class:`~.Processes`, and :class:`~.CompiledProcesses` classes for stoichiometric process and its kinetics.
- Added an :class:`~.Equipment` class for design and costing of unit equipment.
- For the ``stats`` module:

	- More statistical tests:

		- :func:`qsdsan.stats.fast_analysis` for (extended) Fourier amplitude sensitivity test (FAST) and random balance design (RBD) FAST.
		- :func:`qsdsan.stats.morris_till_convergence` to run Morris analysis until the results converge.
		- Added Kendall's tau and Kolmogorov–Smirnov test to :func:`qsdsan.stats.get_correlations`.
	
	- Plotting functions to visualize all test results:

		- :func:`qsdsan.stats.plot_uncertainties` fpr results from uncertainty analysis as different 1D or 2D plots.
		- :func:`qsdsan.stats.plot_correlations` for results from :func:`qsdsan.stats.get_correlation`.
		- Bar plot option for :func:`qsdsan.stats.plot_morris_results`.
		- :func:`qsdsan.stats.plot_morris_convergence` to plot :math:`{\mu^*}` against the number of trajectories.
		- :func:`qsdsan.stats.plot_fast_results` for results from FAST and/or RBD-FAST analyses.
		- :func:`qsdsan.stats.plot_sobol_results` for results from Sobol analysis.

- Changed all .csv data files to .tsv so that they can be viewed on GitHub.
- Added more clear guidelines on `contribution <https://qsdsan.readthedocs.io/en/latest/CONTRIBUTING.html>`_ and a `author list <https://qsdsan.readthedocs.io/en/latest/AUTHORS.html>`_ in the document.


`0.1.0`_
--------
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


`0.0.3`_
--------
- More flexible setting of :class:`~.ImpactItem` for :class:`~.WasteStream`.
- Add status badge to README.rst
- Add CHANGELOG.rst
- Tutorial updates:

	- New:
		- :class:`~.TEA` and :class:`~.LCA`
	- Updated:
		-  :class:`~.Component` and :class:`~.WasteStream`
		-  :class:`~.SanUnit` and :class:`~.System`


`0.0.2`_
--------
- Added the all three sanitation scenarios as described in `Trimmer et al.`_, including uncertainty/sensitivity analyses with tutorial.
- Inclusion of GPX models for estimation of :class:`~.WasteStream` properties.
- Live documentation for the `latest`_ and `beta`_ version.
- New classes:

    - All units in `Trimmer et al.`_
    - Added descriptors (``qsdsan.utils.descriptors``) and decorators (``qsdsan.utils.checkers``) to check user-input values.
    - :class:`~.utils.setters.AttrSetter`, :class:`~.utils.setters.DictAttrSetter`, and :class:`~.utils.getters.FuncGetter` for batch-setting of uncertainty analysis parameters.

- Added :func:`save_report` function to :class:`~.LCA` for report exporting.


`0.0.1`_
--------
- First public release.


.. Other links
.. _latest: https://qsdsan.readthedocs.io/en/latest/
.. _beta: https://qsdsan.readthedocs.io/en/beta/
.. _Trimmer et al.: https://doi.org/10.1021/acs.est.0c03296

.. Commit links
.. _0.3.0: https://github.com/QSD-Group/QSDsan/commit/3c19aebd5503433120217228c3388533cee4bd30
.. _0.2.0: https://github.com/QSD-Group/QSDsan/commit/286943eb206ebd89f58e50b9fdd1bed486e894ae
.. _0.1.0: https://github.com/QSD-Group/QSDsan/commit/1c3d11d9f72421c8b5dbdf6b537775ca35ec65c0
.. _0.0.3: https://github.com/QSD-Group/QSDsan/commit/e20222caccc58d9ee414ca08d8ec55f3a44ffca7
.. _0.0.2: https://github.com/QSD-Group/QSDsan/commit/84653f5979fbcd76a80ffb6b22ffec1c5ca2a084
.. _0.0.1: https://github.com/QSD-Group/QSDsan/commit/f95e6172780cfe24ab68cd27ba19837e010b3d99

