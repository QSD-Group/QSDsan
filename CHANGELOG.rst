Change Log
==========

This document records notable changes to `QSDsan <https://github.com/QSD-Group/QSDsan>`_. We aim to follow `Semantic Versioning <https://semver.org/>`_.


`1.5.3`_
--------
- :func:`qsdsan.default` now takes ``utilities``, ``CEPCI``, and ``flowsheet`` arguments (all default ``True``, mirroring ``biosteam.default``) so users can choose what to reset; ``flowsheet`` defaults to ``True`` (vs. BioSTEAM's ``False``) to preserve prior behavior.

- Constructing a :class:`~.Component` by name/identifier now produces the same fully-initialized component as :meth:`~.Component.from_chemical`, including phase-overridden species (which previously kept an enthalpy model referenced to the wrong phase, giving a wrong reference enthalpy or a ``simulate()`` crash). ``search_ID`` now takes only a lookup name and ``chemical`` only a pre-built ``thermosteam.Chemical``; ``from_chemical`` is a thin wrapper. One behavior change: ``from_chemical(chemical='<name>', formula=<override>)`` now reports the override formula's molecular weight.

- Added a ``formula_override`` argument (default ``False``) to :class:`~.Component`/:meth:`~.Component.from_chemical` that gates overriding a resolved chemical's formula with a different atomic composition (a re-spelling of the same atoms or a blank custom component is not gated). EXPOsan components that legitimately override now pass ``formula_override=True`` (``Struvite``, ``FilterBag``, ``Zeolite``).

- Added runnable docstring ``Examples`` to many built-in unit operations (sanitation, static, and wastewater-treatment units), backed by new :func:`qsdsan.utils.create_example_sanitation_components` and :func:`qsdsan.utils.create_example_wwt_components` helpers; all run under ``pytest --doctest-modules``. Fixed the :class:`~.unit_operations.Incinerator` example and added BioSTEAM API ``See Also`` links to the ``bst``-namespace units.

- Fixed equipment and construction lifetimes so they reach ``equipment_lifetime``, the attribute TEA and LCA actually read (they were previously written to ``_default_equipment_lifetime`` and silently dropped, so replacement costs/impacts were never charged). This shifts TEA/LCA results for systems with sub-project-lifetime equipment, within existing EXPOsan tolerances. ``lifetime`` is now the recommended name for a unit's equipment lifetime.

- Fixed :class:`~.TEA` net present value and replacement-cost calculations for an ``AgileSystem`` (per-unit capital objects expose ``equipment_lifetime``, not ``lifetime``).

- Fixed :meth:`~.TEA.get_unit_annualized_equipment_cost`: a shadowed loop variable mis-computed the annualized cost when a unit declared per-equipment lifetimes; each item is now annualized over its own lifetime.

- Fixed :class:`~.unit_operations.DynamicInfluent` for non-cyclic input data: the padding branch used the removed ``DataFrame.append`` (a silent no-op before pandas 2.0, an ``AttributeError`` after); it now uses ``pd.concat``. Default cyclic input is unchanged.

- Added a public :attr:`~.Process.dynamic_parameters` accessor on :class:`~.Process`/:class:`~.CompiledProcesses` (mirrors :attr:`~.Process.parameters`), so users no longer reach into the private ``_dyn_params``.

- :meth:`~.Processes.load_from_file` now supports per-process conservation rules: ``conserved_for`` is keyword-only and required, accepting a uniform tuple, a per-process dict, or ``None`` to defer to a new ``conserved_for`` data-file column. Backward incompatibility: callers that omitted ``conserved_for`` must now pass it explicitly.

- Fixed the SanUnit-mixin plumbing on eight ``bst``-namespace wrappers (distillation columns, ``Flash``, ``IsothermalCompressor``, ``ProcessWaterCenter``, ``HeatExchangerNetwork``) where ``SanUnit.__init__`` never ran, so LCA/add-on attributes and kwargs (``lifetime``, ``add_OPEX``, etc.) were silently dropped. The mixin-install logic was extracted into a shared ``self._init_sanunit_addons(...)``.

- Harmonized the six ``qsdsan.stats`` plotting helpers: every ``plot_*`` now accepts ``ax=None`` and ``**plot_kws`` forwarded to the underlying plotting call. ``plot_correlations``'s ``**kwargs`` was renamed to ``**plot_kws`` (named callers unaffected), and a latent ``ax=``-dropping bug in ``plot_sobol_results`` was fixed.

- Made ``qsdsan`` a more self-contained namespace so users rarely need to import biosteam/thermosteam directly. Newly re-exported: ``settings``, ``preferences``, ``stream_utility_prices``, ``Thermo``, ``UtilityAgent``, ``Facility``, ``AgileSystem``, ``get_OSBL``, ``MissingStream``, and ``report``; ``qsdsan.unit_operations.bst`` surfaces more un-customized BioSTEAM units; ``qsdsan.utils`` re-exports ``rho_to_V``, ``V_to_rho``, ``@cost``, ``var_columns``/``var_indices``. A new :ref:`BioSTEAM API <biosteam_api>` page lists the full surface.

- ``qsdsan.utils`` now wraps BioSTEAM's ``@cost`` decorator to also accept ``CEPCI`` as an alias for ``CE`` (BioSTEAM's ``CE`` keyword still works).

- Fixed ``qsdsan.utils.indices.ChemPPI_by_year``: the 2021/2022 entries had been mistakenly populated with CEPCI values and have been removed, so the series ends at 2020 until authentic values are added.

- Fixed compatibility with newer Thermosteam (``>=0.53.4``) where ``Chemical.MW`` is read-only: the ``MW = 1.`` placeholder for formula-less components now goes through ``qsdsan._compat.set_chemical_MW``. This unblocked ``Components.load_default`` and nearly every system build.

- Unified report generation across :class:`~.System`, :class:`~.TEA`, and :class:`~.LCA`: ``save_report`` on any of the three now produces the same workbook (with an ``LCA`` sheet when the system has an :class:`~.LCA`). :meth:`~.LCA.save_report` delegates to it and its default filename changed to ``{system.ID}_report.xlsx``; use :meth:`~.LCA.get_impact_table` for the LCA tables alone.

- Added a ``time_frame`` argument (``'lifetime'``/``'all'`` default, ``'yr'``, ``'month'``, ``'week'``, ``'day'``, ``'hr'``) to the :class:`~.LCA` results methods. The existing ``annual`` flag is kept as a backward-compatible alias for ``time_frame='yr'``.

- Added :meth:`~.LCA.get_normalized_impacts` (impacts per kg/m³/MJ of reference streams, the LCA counterpart to ``TEA.solve_price``) and :meth:`~.LCA.get_allocated_impact_table` (tabular counterpart of :meth:`~.LCA.get_allocated_impacts`, available as an opt-in report sheet).

- Several :class:`~.LCA` fixes: :meth:`~.LCA.get_unit_impacts` no longer double-counts stream impacts; :meth:`~.LCA.get_allocated_impacts` accepts a function for ``allocate_by``; :meth:`~.LCA.get_impact_table` writes the ``Sum`` row via ``.loc`` (fixing blank totals under pandas ≥ 2.x) and handles empty ``'Stream'``/``'Other'`` categories under pandas 3.0; :meth:`~.LCA.add_other_item` now names the missing ID in its error.

- Added :func:`qsdsan.utils.create_example_treatment_systems` (the aerobic/anaerobic wastewater systems shared by the TEA and LCA tutorials); the LCA tutorial now uses it instead of the ``bwaise`` EXPOsan system.

- In EXPOsan, restored the documented ``exposan.bsm1.cmps``/``components``/``asm`` module attributes (their assignment in ``bsm1.load()`` was commented out, so they stayed ``None`` after a fresh ``load()``), and replaced deprecated per-test ``clear_lca_registries()`` calls with an autouse ``conftest.py`` fixture.

- Substantially expanded and reorganized the tutorials and FAQ (new tutorials and sections on operational flexibility, recycle convergence, specifications, auxiliary units, cost/environmental trade-offs, and the dynamic process-modeling tutorials; a multi-page FAQ; dark-mode rendering fixes).


`1.5.2`_
--------
- Fixed packaging: ``qsdsan/units_of_measure.txt`` (the pint unit-definition file loaded at import) was not declared in ``package-data``, so non-editable installs (wheels) omitted it and ``import qsdsan`` raised ``FileNotFoundError``. It is now included in the distributed package.

- Fixed :meth:`~.SanStream.copy_like`: when called with ``copy_price=True`` or ``copy_impact_item=True``, the price/impact item were copied in the wrong direction (overwriting the source stream and leaving the target unchanged). They are now correctly copied from the source into the target.

- Added doctest examples to :meth:`~.SanStream.copy`, :meth:`~.SanStream.copy_like`, and :meth:`~.SanStream.copy_flow` documenting what each method copies (flows, temperature/pressure, price, and impact item).

- Added ``.. warning::`` notes to :attr:`~.WasteStream.pH` and :attr:`~.WasteStream.SAlk` clarifying that these are not calculated from stream composition (no acid-base model yet) and should be treated as user-provided inputs.

- Fixed :class:`~.sanunits.HXutility`: its ``_units`` dictionary was inadvertently set to ``None`` (it was assigned the return value of ``dict.update``, which is always ``None``), which broke :meth:`results` with an ``AttributeError`` and mutated BioSTEAM's shared ``_units``. The unit-of-measure entries are now built with a dict merge.

- Re-exported the Thermosteam reaction classes (``Reaction``, ``ReactionItem``, ``ReactionSet``, ``ParallelReaction``, ``SeriesReaction``, ``ReactionSystem``, and the ``Rxn``/``RxnI``/``RxnS``/``PRxn``/``SRxn``/``RxnSys`` aliases) from the top-level ``qsdsan`` namespace, so they can be imported with ``from qsdsan import Reaction`` instead of reaching into BioSTEAM/Thermosteam.

- Added ``qsdsan.CEPCI``, a settable view of the global Chemical Engineering Plant Cost Index (which BioSTEAM abbreviates as ``CE``), so the costing index can be read and set (e.g., ``qsdsan.CEPCI = qsdsan.CEPCI_by_year[2023]``) without importing biosteam.

- Fixed :class:`~.TEA`: its ``CEPCI`` argument defaulted to ``bst.CE``, which Python binds once at import time (freezing it at 567.5). Every ``TEA`` created without an explicit ``CEPCI`` therefore reset the global cost index, silently overriding a deliberately set ``qsdsan.CEPCI``/``bst.CE``. ``CEPCI`` now defaults to ``None`` (the current index is left untouched), and a provided ``CEPCI`` is applied *before* simulation so it actually affects costing.

- :attr:`~.TEA.CEPCI_by_year` now returns ``qsdsan.CEPCI_by_year`` (was BioSTEAM's table) for consistency, and ``qsdsan.CEPCI_by_year`` now merges in BioSTEAM's ``design_tools.CEPCI_by_year`` years not already present (e.g., pre-1990), with qsdsan's more precise values taking precedence on overlapping years. QSDsan's built-in hydroprocessing/hydrothermal units now source their ``@cost`` reference indices from ``qsdsan.CEPCI_by_year`` as well (a sub-0.1% cost change from the more precise values).

- In EXPOsan, systems that set the global cost index (``htl`` and ``saf``) now do so via ``qsdsan.CEPCI = ...`` instead of ``bst.CE = ...``, for consistency with the new ``qsdsan.CEPCI`` handle (functionally identical, since ``qsdsan.CEPCI`` is a live view of ``bst.CE``).

- Renamed the cost-index dicts in ``qsdsan.utils.indices`` and the keys of ``qsdsan.utils.tea_indices`` to the ``*_by_year`` convention (``'CEPCI_by_year'``, ``'ChemPPI_by_year'``, ``'labor_by_year'``, ``'PCEPI_by_year'``). Update any code that used the old short keys, e.g., ``tea_indices['CEPCI']`` becomes ``tea_indices['CEPCI_by_year']``.

- Tutorial updates ongoing.


`1.5.1`_
--------
- :class:`~.Component` now validates the ``particle_size``, ``degradability``, and ``organic`` arguments at creation, for both the constructor and :meth:`~.Component.from_chemical`. Invalid values (e.g., a misspelled ``particle_size``) that were previously accepted silently now raise a ``ValueError``.

- The ``f_BOD5_COD``, ``f_uBOD_COD``, and ``f_Vmass_Totmass`` fractions are now range-checked to ``[0, 1]`` (this check was previously unreachable).

- The :class:`~.Component` constructor now accepts a ``chemical`` keyword to build a component directly from an existing ``thermosteam.Chemical``. :meth:`~.Component.from_chemical` is now a thin wrapper around it, with unchanged behavior.

- Documented the ``ignore_inaccurate_molar_weight`` and ``adjust_MW_to_measured_as`` options of :meth:`~.Components.compile`, and substantially revised the topical tutorials.


`1.5.0`_
--------
- All LCA registry types (:class:`~.ImpactIndicator`, :class:`~.ImpactItem`, :class:`~.Construction`, :class:`~.Transportation`) are now isolated per flowsheet. Switching between systems via :meth:`~.SanMainFlowsheet.set_flowsheet` atomically swaps all four registries, so no manual :func:`~.utils.clear_lca_registries` calls are needed between systems. :func:`~.utils.clear_lca_registries` is deprecated.

- Added :attr:`SanUnit._construction_specs <.SanUnit._construction_specs>` — a class-level tuple of dicts for declaring default construction materials. Specs are resolved lazily by :class:`~.LCA` at creation time, so :class:`~.ImpactItem` objects do not need to exist when the unit is instantiated.

- Reorganized unit operations and process models into clearer namespaces.

  - ``qsdsan.sanunits`` is renamed to :mod:`qsdsan.unit_operations` and reorganized into three behavior-based sub-namespaces:

    - :mod:`qsdsan.unit_operations.bst` — BioSTEAM-inherited unit operations (mixers, splitters, pumps, heat exchangers, distillation columns, tanks, etc.)
    - :mod:`qsdsan.unit_operations.static` — steady-state QSDsan unit operations (sanitation fixtures, treatment beds, clarifiers, sludge handling, hydrothermal/hydroprocessing units, etc.)
    - :mod:`qsdsan.unit_operations.dynamic` — unit operations with explicit dynamic-state behavior (bioreactors, dynamic influent, junctions, membrane bioreactors, etc.)

  - ``qsdsan.processes`` is renamed to :mod:`qsdsan.process_models`.
  - All existing imports via ``qsdsan.sanunits`` and ``qsdsan.processes`` remain valid for backward compatibility.

- Restructured API documentation to mirror the new package layout, with dedicated pages for each sub-namespace.

- Bug fixes in :func:`~.process_models.rhos_asm2d`:

  - Added ``if X_MeOH > 0`` and ``if X_MeP > 0`` guards to the precipitation and redissolution reactions, consistent with the existing guards for ``X_H``, ``X_PAO``, and ``X_AUT``. Without these guards, the BDF solver's polynomial extrapolation could produce tiny positive floating-point values for ``X_MeOH``, causing spurious ``X_MeP`` accumulation over long simulations.
  - Added a ``S_F + S_A > 0`` guard to the heterotrophic growth substrate-partitioning terms. When the BDF Newton iterations drive both fermentable substrate (``S_F``) and acetate (``S_A``) to zero simultaneously, the partition fractions ``S_F/(S_F+S_A)`` and ``S_A/(S_F+S_A)`` produce ``0/0`` and raise a ``FloatingPointError``; zero substrate correctly implies zero growth.

- Multiple bug fixes and improvements to :class:`~.unit_operations.PolishingFilter`:

  - Moved O2 deficit calculation before the effluent split so dissolved oxygen is correctly accounted for in all outlet streams.
  - Added an aerobic-only guard for air injection: air is now only added when ``self.has_pump`` is ``False`` and the unit operates in aerobic mode.
  - Fixed a missing ``_freeboard`` attribute that caused ``AttributeError`` on initialization.
  - Restored the ``_design_anaerobic`` method that had been inadvertently removed.
  - Corrected the slab concrete volume formula.
  - Added a ``biomass_ID`` parameter to allow users to specify which component tracks active biomass.
  - Renamed internal attributes ``gas``/``soluble``/``solid`` to ``gases``/``solubles``/``solids`` for consistency with the rest of the codebase.
  - Fixed the condition in ``get_digestion_rxns`` and corrected the argument order in ``_refresh_rxns``.
  - Fixed a ``SanStream.degassing`` ``AttributeError`` that occurred when the polishing filter effluent was a plain ``SanStream``.
  - Fixed an O2 double-counting error in ``air_out`` that produced a ~0.6% mass-balance error.

- Fixed :class:`~.unit_operations.HydraulicDelay`: added missing ``_update_state`` and ``_update_dstate`` methods so the unit correctly propagates state during dynamic simulation.

- Added a deprecation warning to :class:`~.SimpleTEA` to guide users toward the updated TEA interface.

- Lazy-imported optional heavy dependencies (``SALib``, ``seaborn``, ``chaospy``, ``sympy``) so that ``import qsdsan`` no longer pays the startup cost of those libraries unless they are actually used.

- Added a GitHub Actions release workflow that automatically publishes to PyPI and creates a GitHub release when a ``v*.*.*`` tag is pushed.


`1.4.0`_
--------
- A lot of the updates have been focused on the dynamic simulation, now the open-loop Benchmark Simulation Model No. 2 (`BSM2 <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm2>`_) configuration has been implemented with new process models and unit operation including

	- :class:`qsdsan.processes.ADM1p`
	- :class:`qsdsan.processes.ADM1_p_extension`
	- :class:`qsdsan.processes.ModifiedADM1`
	- :class:`qsdsan.processes.mASM2d`
	- :class:`qsdsan.sanunits.IdealClarifier`
	- :class:`qsdsan.sanunits.PrimaryClarifier`
	- :class:`qsdsan.sanunits.PrimaryClarifierBSM2`
	- :class:`qsdsan.sanunits.GasExtractionMembrane`
	- :class:`qsdsan.sanunits.Thickener`
	- :class:`qsdsan.sanunits.Centrifuge`
	- :class:`qsdsan.sanunits.Incinerator`
	- :class:`qsdsan.sanunits.BatchExperiment`
	- :class:`qsdsan.sanunits.PFR`
	- :class:`qsdsan.sanunits.BeltThickener`
	- :class:`qsdsan.sanunits.SludgeCentrifuge`
	- :class:`qsdsan.sanunits.SludgeThickener`

- New publications

	- Feng et al., *Environmental Science & Technology*, on the sustainability of `hydrothermal liquefaction (HTL) <https://doi.org/10.1021/acs.est.3c07394>`_ for resource recovery from a range of wet organic wastes.


`1.3.0`_
--------
- Enhance and use QSDsan's capacity for dynamic simulation for emerging technologies and benchmark configurations (see EXPOsan METAB and PM2 (on the algae branch, still under development) modules).
- New publications

	- The paper introducing `DMsan <https://doi.org/10.1021/acsenvironau.2c00067>`_, the package developed for decision-making of sanitation and resource recovery technologies, is published in *ACS Environmental Au*!
	- QSDsan was used to evaluate the sustainability of the `NEWgenerator <https://doi.org/10.1021/acsenvironau.3c00001>`_ system as in this paper on *ACS Environmental Au*!

- New modules

	- :class:`qsdsan.processes.KineticReaction`

- ``QSDsan`` now has a `website <https://qsdsan.com/>`_ to host all of the resources!
- ``QSDsan``'s `documentation <https://qsdsan.readthedocs.io/en/latest/index.html>`_ is getting a new look!
- Add new units to enable dynamic simulation of systems with multiple process models. Check out :class:`qsdsan.sanunits.Junction`, :class:`qsdsan.sanunits.ADMtoASM`, :class:`qsdsan.sanunits.ASMtoADM` and their use in the `interface system demo <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/interface>`_.
- In `online testing <https://github.com/QSD-Group/QSDsan/actions>`_, we dropped the test for Python 3.8 and added Python 3.10. The main developing environment for QSDsan is 3.9.


`1.2.0`_
--------
- The `QSDsan paper <https://www.doi.org/10.1039/d2ew00455k>`_ is accepted by *Environmental Science: Water Research & Technology*!
- The first paper using QSDsan for the design of sanitation is accepted by *ACS Environmental Au*! Read the `Biogenic Refinery <https://pubs.acs.org/doi/10.1021/acsenvironau.2c00022>`_ paper and check out the system module in ``QSDsan``/``EXPOsan``.
- Added multiple systems (including their unit operations), check out the details on the `Developed System <https://qsdsan.readthedocs.io/en/latest/Developed_Systems.html>`_ page!

	- Biogenic Refinery
	- Eco-San
	- Reclaimer

- Added the anaerobic digestion model no. 1 (ADM1) process model and the unit :class:`qsdsan.sanunits.AnaerobicCSTR`, the corresponding `system <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/adm>`_ can be found in EXPOsan.
- Other new unit operations:

	- Encapsulation Bioreactors:

		- :class:`qsdsan.sanunits.CH4E`
		- :class:`qsdsan.sanunits.H2E`


`1.1.0`_
--------
- Fully tested dynamic simulation capacity, refer to the `BSM1 system <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm1>`_ in EXPOsan for an example implementation.
- Added many new :class:`qsdsan.SanUnit` and reorganized package/documentation structure, new unit operations include:

	- :class:`qsdsan.sanunits.AnMBR`
	- :class:`qsdsan.sanunits.CHP`
	- :class:`qsdsan.sanunits.InternalCirculationRx`
	- :class:`qsdsan.sanunits.SludgeHandling`

		- :class:`qsdsan.sanunits.BeltThickener`
		- :class:`qsdsan.sanunits.SludgeCentrifuge`

	- :class:`qsdsan.sanunits.PolishingFilter`
	- :class:`qsdsan.sanunits.WWTpump`

- Continue to enhance documentation (e.g., :class:`qsdsan.Process`, `qsdsan.stats`, util functions).


`1.0.0`_
--------
Official release of ``QSDsan`` v1.0.0!

- Added system-wise dynamic simulation capacity. To use the dynamic simulation function, a unit needs to have several supporting methods to initialize its state and compile ordinary differential equations (ODEs), refer to the units included in the BSM1 system below for usage, documentation and tutorial will be coming soon!
- Developed the `benchmark simulation system no.1 (BSM1) model on EXPOsan <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm1>`_ with comparison against the MATLAB/Simulink model developed by the International Water Association (IWA) Task Group on Benchmarking of Control Strategies. See the `README <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm1>`_ for details
- Significantly expanded the tutorials with demo videos on `YouTube <https://www.youtube.com/playlist?list=PL-tj_uM0mIdFv72MAULnWjS6lx_cCyi2N>`_. Now tutorials cover all non-dynamic major classes (tutorials on dynamic classes will be included in the next major release).


`0.3.0`_
--------
- Now LCA data can be imported from external databases using the newly made `BW2QSD <https://github.com/QSD-Group/BW2QSD>`_ package.
- New subclasses of :class:`qsdsan.SanUnit`:

	- :class:`qsdsan.sanunits.Clarifier`
	- :class:`qsdsan.sanunits.CSTR`

	- :class:`qsdsan.sanunits.ElectrochemicalCell` using the following :class:`qsdsan.Equipment`:

		- :class:`qsdsan.equipments.Column`
		- :class:`qsdsan.equipments.Electrode`
		- :class:`qsdsan.equipments.Machine`
		- :class:`qsdsan.equipments.Membrane`

- New subclasses of :class:`qsdsan.Process`:

	- :class:`qsdsan.processes.DiffusedAeration`
	- :class:`qsdsan.processes.ASM1`
	- :class:`qsdsan.processes.ASM2d`

- Updated :class:`qsdsan.SanUnit` so that it can be initialized with any of :class:`thermosteam.Stream`, :class:`qsdsan.SanStream`, or :class:`qsdsan.WasteStream`.

	- These three classes can now be mixed.

- Added :class:`qsdsan.SanStream` for non-waste streams (e.g., gases).
- Updated the ``add_OPEX`` attribute of :class:`qsdsan.SanUnit` and ``system_add_OPEX`` attribute of :class:`qsdsan.SimpleTEA` so that they take :class:`dict` as the default to allow display of multiple additional operating expenses.
- Split the ``systems`` module into an individual package `EXPOsan`_.
- Now using :class:`thermosteam.utils.Registry` to manage :class:`qsdsan.ImpactIndicator` and :class:`qsdsan.ImpactItem`.
- Added `AppVeyor CI <https://ci.appveyor.com/project/yalinli2/qsdsan>`_.
- Renamed the ``master`` branch to ``main``.


`0.2.0`_
--------
- Added :class:`qsdsan.Process`, :class:`qsdsan.Processes`, and :class:`qsdsan.CompiledProcesses` classes for stoichiometric process and its kinetics.
- Added an :class:`qsdsan.Equipment` class for design and costing of unit equipment.
- For the ``stats`` module:

	- More statistical tests:

		- :func:`qsdsan.stats.fast_analysis` for (extended) Fourier amplitude sensitivity test (FAST) and random balance design (RBD) FAST.
		- :func:`qsdsan.stats.morris_till_convergence` to run Morris analysis until the results converge.
		- Added Kendall's tau and Kolmogorov–Smirnov test to :func:`qsdsan.stats.get_correlations`.

	- Plotting functions to visualize all test results:

		- :func:`qsdsan.stats.plot_uncertainties` for results from uncertainty analysis as different 1D or 2D plots.
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
- :func:`LCA.get_normalized_impacts` was replaced by :func:`qsdsan.LCA.get_allocated_impacts` for :class:`qsdsan.LCA` to enable flexible allocation options.
- Reformatted all documents, added instructions on documentation.
- Added brief instructions on contributing and code of conduct.
- Updated UML diagram.


`0.0.3`_
--------
- More flexible setting of :class:`qsdsan.ImpactItem` for :class:`qsdsan.WasteStream`.
- Add status badge to README.rst
- Add CHANGELOG.rst
- Tutorial updates:

	- New:
		- :class:`qsdsan.TEA` and :class:`qsdsan.LCA`
	- Updated:
		-  :class:`qsdsan.Component` and :class:`qsdsan.WasteStream`
		-  :class:`qsdsan.SanUnit` and :class:`qsdsan.System`


`0.0.2`_
--------
- Added the all three sanitation scenarios as described in `Trimmer et al.`_, including uncertainty/sensitivity analyses with tutorial.
- Inclusion of GPX models for estimation of :class:`qsdsan.WasteStream` properties.
- Live documentation for the `latest`_ and `beta`_ version.
- New classes:

    - All units in `Trimmer et al.`_
    - Added descriptors (``qsdsan.utils.descriptors``) and decorators (``qsdsan.utils.checkers``) to check user-input values.
    - :class:`qsdsan.utils.setters.AttrSetter`, :class:`qsdsan.utils.setters.DictAttrSetter`, and :class:`qsdsan.utils.getters.FuncGetter` for batch-setting of uncertainty analysis parameters.

- Added :func:`save_report` function to :class:`qsdsan.LCA` for report exporting.


`0.0.1`_
--------
- First public release.


.. Other links
.. _latest: https://qsdsan.readthedocs.io/en/latest
.. _beta: https://qsdsan.readthedocs.io/en/beta
.. _EXPOsan:  https://github.com/QSD-Group/exposan
.. _Trimmer et al.: https://doi.org/10.1021/acs.est.0c03296

.. Commit links
.. _1.5.3: https://github.com/QSD-Group/QSDsan/releases/tag/v1.5.3
.. _1.5.2: https://github.com/QSD-Group/QSDsan/releases/tag/v1.5.2
.. _1.5.1: https://github.com/QSD-Group/QSDsan/releases/tag/v1.5.1
.. _1.5.0: https://github.com/QSD-Group/QSDsan/releases/tag/v1.5.0
.. _1.4.0: https://github.com/QSD-Group/QSDsan/releases/tag/v1.4.0
.. _1.3.0: https://github.com/QSD-Group/QSDsan/releases/tag/v1.3.0
.. _1.2.0: https://github.com/QSD-Group/QSDsan/releases/tag/v1.2.0
.. _1.1.0: https://github.com/QSD-Group/QSDsan/releases/tag/v1.1.0
.. _1.0.0: https://github.com/QSD-Group/QSDsan/releases/tag/v1.0.0
.. _0.3.0: https://github.com/QSD-Group/QSDsan/releases/tag/v0.3.0
.. _0.2.0: https://github.com/QSD-Group/QSDsan/commit/286943eb206ebd89f58e50b9fdd1bed486e894ae
.. _0.1.0: https://github.com/QSD-Group/QSDsan/commit/1c3d11d9f72421c8b5dbdf6b537775ca35ec65c0
.. _0.0.3: https://github.com/QSD-Group/QSDsan/commit/e20222caccc58d9ee414ca08d8ec55f3a44ffca7
.. _0.0.2: https://github.com/QSD-Group/QSDsan/commit/84653f5979fbcd76a80ffb6b22ffec1c5ca2a084
.. _0.0.1: https://github.com/QSD-Group/QSDsan/commit/f95e6172780cfe24ab68cd27ba19837e010b3d99
