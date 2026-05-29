Change Log
==========

This document records notable changes to `QSDsan <https://github.com/QSD-Group/QSDsan>`_. We aim to follow `Semantic Versioning <https://semver.org/>`_.


`1.5.3`_
--------
- Added ``tutorials/14_Modeling_Notes_and_Pitfalls.ipynb`` covering common modeling surprises across streams/components, unit design, TEA/LCA, and behavior inherited from BioSTEAM. Extended tutorial 5 §1.1 with ``simulate()`` / ``_summary`` call-graph mechanics.

- Restructured the FAQ into a four-page ``faq/`` subdirectory (``errors``, ``tips``, ``styling``, ``ai_assisted_coding``), absorbing the previously standalone AI-assisted coding page from the tutorials section. The old ``FAQ.html`` bookmark redirects to the new index.

- Fixed :class:`~.unit_operations.DynamicInfluent` for non-cyclic input data: when the first and last rows of the data file differed, the constructor entered a branch that called ``df.append(y_end, ignore_index=True)`` to pad the time series with a phantom final point. The line was always a silent no-op (``DataFrame.append`` returned a new frame rather than modifying in place; the result was discarded), and starting with pandas 2.0 ``append`` was removed entirely, so the same call now raises ``AttributeError: 'DataFrame' object has no attribute 'append'`` and the constructor fails outright. The branch now uses ``pd.concat`` and reassigns ``self._data`` so the phantom point is actually appended; ``self._t_end`` and the interpolants pick it up. Default-file behavior is unchanged because the shipped ``_inf_dry_2006.tsv`` is cyclic (first row equals last row) and the branch is skipped.

- In EXPOsan, restored the documented ``exposan.bsm1.cmps`` / ``components`` / ``asm`` module attributes (the canonical handles used by the dynamic-simulation tutorial). ``bsm1.load()`` declared the three names ``global`` but the matching assignment was commented out, so after a fresh ``load()`` they stayed ``None`` and ``qs.set_thermo(bsm1.cmps)`` raised ``TypeError: 'NoneType' object is not iterable``. The assignments (``cmps = components = O1.components``, ``asm = O1.suspended_growth_model``) are now active.

- Added a public :attr:`~.Process.dynamic_parameters` accessor on :class:`~.Process` and :class:`~.CompiledProcesses` that returns the dictionary of attached :class:`~.DynamicParameter` objects (keyed by symbol; empty when all parameters are static). This mirrors the existing :attr:`~.Process.parameters` property, so users no longer have to reach into the private ``_dyn_params`` attribute to inspect what dynamic parameters a process has. Added test coverage in ``tests/test_process.py``.

- :meth:`~.Processes.load_from_file` now supports per-process conservation rules through three composable forms. Previously every process in a batch-built model shared one ``conserved_for`` tuple, forcing users to fall back to per-process manual construction whenever (e.g.) decay shouldn't be required to conserve carbon. ``conserved_for`` is now keyword-only and required (no default value) and accepts: a **tuple** (uniform — applied to every process, unchanged from before), a **dict** keyed by process ID (per-process rules; IDs absent from the dict fall back to ``('COD', 'N', 'P', 'charge')``), or **None** (defers to a new optional ``conserved_for`` column in the data file, which carries comma-separated material names per row; empty cells mean no enforcement). The file column is consumed during parsing and does not appear in the stoichiometry matrix. Resolution order, highest first: kwarg-dict entry for the listed process ID → kwarg-tuple (uniform) → file column → default tuple. Backward incompatibility: existing callers that omitted ``conserved_for`` (relying on the previous tuple default) will now raise ``TypeError`` and must pass an explicit value. The tutorial's example datafile is now split into ``_bkm.csv`` (no column, demonstrates the tuple-then-dict workflow) and ``_bkm_with_conserved.csv`` (column-bearing, demonstrates ``conserved_for=None``). Added test coverage in ``tests/test_process.py``.

- Fixed the SanUnit-mixin plumbing on eight ``bst``-namespace wrappers (:class:`~.unit_operations.BinaryDistillation`, :class:`~.unit_operations.ShortcutColumn`, :class:`~.unit_operations.MESHDistillation`, :class:`~.unit_operations.AdiabaticMultiStageVLEColumn`, :class:`~.unit_operations.Flash`, :class:`~.unit_operations.IsothermalCompressor`, :class:`~.unit_operations.ProcessWaterCenter`, :class:`~.unit_operations.HeatExchangerNetwork`) whose MRO put the BioSTEAM parent ahead of ``SanUnit`` (or, for ``HeatExchangerNetwork``, explicitly overwrote ``__init__``). ``SanUnit.__init__`` never ran, so the LCA/add-on attributes (``construction``, ``transportation``, ``equipment``, ``add_OPEX``, ``uptime_ratio``, ``lifetime``, ``include_construction``) silently did not exist on instances of these classes, and kwargs like ``lifetime=`` or ``add_OPEX=`` were dropped. Each affected wrapper now defines an ``__init__`` that calls the BioSTEAM parent and then ``self._init_sanunit_addons(...)`` (a new method extracted from ``SanUnit.__init__`` so the mixin-install logic has a single source of truth). Added a three-layer test for every ``bst`` class in ``tests/test_bst_units.py`` (smoke + parity-where-feasible + add-on persistence through ``simulate()``) so any future mixin regression fails CI rather than silently propagating.

- API: harmonized the six plotting helpers in ``qsdsan.stats`` so users can place a plot on a chosen axis and forward styling kwargs uniformly. Every ``plot_*`` now accepts ``ax=None`` and ``**plot_kws`` (forwarded to the primary underlying plotting call, i.e., the ``seaborn`` function or the matplotlib ``scatter``/``bar``/``errorbar``). On ``plot_uncertainties``, the existing ``center_kws``/``margin_kws`` are kept because the 2D ``seaborn.JointGrid`` has two axes; ``plot_kws`` merges into ``center_kws`` (``center_kws`` wins on conflict). On ``plot_uncertainties`` 2D plots and on ``plot_correlations``, ``ax=`` is documented as ignored because the underlying calls build their own figure. ``plot_correlations``'s old ``**kwargs`` was renamed to ``**plot_kws``; existing callers that pass named kwargs (every caller surveyed in EXPOsan and the tutorials) are unaffected. Also fixed a latent bug in the bar-only branch of ``plot_sobol_results`` that was dropping the user-supplied ``ax=``.

- Made ``qsdsan`` a more self-contained namespace so users typically don't  need to ``import biosteam``/``import thermosteam`` directly. Newly re-exported at the top level: ``settings``, ``preferences``, ``stream_utility_prices``, ``Thermo``, ``UtilityAgent``, ``Facility``, ``AgileSystem``, ``get_OSBL``, ``MissingStream``, and the ``report`` submodule. ``qsdsan.unit_operations.bst`` now also surfaces the BioSTEAM units QSDsan does not customize (``IsenthalpicValve``, ``Stripper``, ``MolecularSieve``, ``BatchBioreactor``, ``VacuumSystem``, ``Boiler``, ``BoilerTurbogenerator``, ``ChilledWaterPackage``, ``CoolingTower``, ``SolidsCentrifuge``), and ``qsdsan.utils`` now re-exports ``rho_to_V``, ``V_to_rho``, the ``@cost`` decorator, and ``var_columns``/``var_indices``. A new :ref:`Public API <public_api>` documentation page lists the full surface, and ``tests/test_public_api.py`` asserts each re-export still matches its BioSTEAM/Thermosteam source (so an upstream rename fails loudly in CI rather than reaching a user).

- Fixed :meth:`~.TEA.get_unit_annualized_equipment_cost` (used by :attr:`~.TEA.annualized_equipment_cost`): when a unit declared per-equipment lifetimes through an ``equipment_lifetime`` dict, the loop variable shadowed the cost accumulator, so the annualized equipment cost was computed incorrectly (the running sum was dropped and the bare-module factors mis-applied). It now iterates each unit's ``installed_costs`` with a separate variable, annualizing every item over its own lifetime (falling back to the unit or TEA lifetime). Added ``tests/test_tea.py`` covering this case along with the other ``TEA`` cost metrics.

- Fixed ``qsdsan.utils.indices.ChemPPI_by_year``: the 2021 and 2022 entries had been mistakenly populated with CEPCI values (708.8, 816.0). They have been removed, so the chemical-PPI series now ends at 2020 until authentic values are added.

- ``qsdsan.utils`` now thinly wraps BioSTEAM's ``@cost`` decorator so it also accepts ``CEPCI`` as an alias for the reference cost-index argument ``CE`` (e.g., ``@cost(..., CEPCI=522)``), for consistency with ``qsdsan.CEPCI``/``qsdsan.CEPCI_by_year``. BioSTEAM's ``CE`` keyword continues to work unchanged.

- Fixed compatibility with newer Thermosteam releases that make ``Chemical.MW`` a read-only property: :class:`~.Component` creation assigned ``self.MW = 1.`` as a placeholder for components without a formula or molecular weight, which raised ``AttributeError: cannot set molecular weight``. The placeholder is now set via ``qsdsan._compat.set_chemical_MW``, which uses the public ``MW`` setter when it is writable and falls back to Thermosteam's internal constant-reset helper when it is read-only (so it works on both old and new Thermosteam, ``>=0.53.4``). This unblocked ``Components.load_default`` and therefore nearly every system build.

- Unified report generation across :class:`~.System`, :class:`~.TEA`, and :class:`~.LCA`: calling ``save_report`` on any of the three now produces the same Excel workbook (system design, costs, and utilities, plus the QSDsan LCA tables on an ``LCA`` sheet whenever the system has an :class:`~.LCA`). This works by wrapping BioSTEAM's ``System.save_report`` at import (``qsdsan.System`` is BioSTEAM's ``System``, and ``TEA.save_report`` already delegates to it); systems without an LCA behave exactly as before. :meth:`~.LCA.save_report` now delegates to the system report, and its default filename changed from ``{system.ID}_lca.xlsx`` to ``{system.ID}_report.xlsx``; for the LCA tables alone, use :meth:`~.LCA.get_impact_table`. Added ``tests/test_lca.py`` covering the unified output. In EXPOsan, ``bwaise``'s ``save_reports`` no longer writes a separate LCA file.

- Added ``qsdsan.utils.create_example_treatment_systems``, which builds the aerobic and anaerobic wastewater treatment systems shared by the TEA and LCA tutorials (a compact, realistic substrate for techno-economic and life cycle analyses). The LCA tutorial now uses it instead of the ``bwaise`` EXPOsan system.

- Added a ``time_frame`` argument to the :class:`~.LCA` results methods (:meth:`~.LCA.get_total_impacts`, :meth:`~.LCA.get_unit_impacts`, :meth:`~.LCA.get_impact_table`, :meth:`~.LCA.get_allocated_impacts`/:meth:`~.LCA.get_allocated_impact_table`, and ``save_report``). It normalizes the results to a chosen time frame: ``'lifetime'`` (or ``'all'``, the default), ``'yr'`` (equivalent to ``annual=True``), ``'month'``, ``'week'``, ``'day'``, or ``'hr'``. The existing ``annual`` flag is kept as a backward-compatible alias (``annual=True`` ≡ ``time_frame='yr'``); ``time_frame`` takes precedence when both are given.

- Added :meth:`~.LCA.get_normalized_impacts`, which expresses impacts per functional unit (per kg, per m³, or per MJ) of one or more reference streams, the LCA counterpart to ``TEA.solve_price``. By default the total impacts are divided by the streams' combined throughput; pass ``allocate_by`` to normalize the impacts allocated to those streams instead.

- Fixed :meth:`~.LCA.get_unit_impacts`: stream impacts were added twice (the accumulator was initialized to the stream-impact dict and then the stream impacts were added again), overstating the result. Each category is now counted once.

- Fixed :meth:`~.LCA.get_allocated_impacts` (and the new :meth:`~.LCA.get_allocated_impact_table`): passing a function for ``allocate_by`` (a documented option) raised ``TypeError`` because ``callable`` was tested *after* ``iter(allocate_by)``, which fails on a non-iterable function. ``callable`` is now checked first, so a function returning the allocation ratios works; an invalid ``allocate_by`` now raises a clear ``ValueError``.

- Added :meth:`~.LCA.get_allocated_impact_table`, the tabular counterpart of :meth:`~.LCA.get_allocated_impacts` (impacts allocated to two or more streams, one row per stream, indicator columns, plus an ``'Allocation factor'`` column). The unified report can include it as an opt-in ``'LCA allocation'`` sheet by passing ``lca_allocate_streams`` (and optionally ``lca_allocate_by``) to ``save_report``. The shared allocation-ratio logic was factored into a private helper used by both methods.

- Fixed :meth:`~.LCA.get_impact_table`: the per-category ``Sum`` row was written via pandas chained assignment, which silently no-ops under Copy-on-Write (pandas ≥ 2.x), so the total row came out blank and ``ChainedAssignmentError`` warnings were emitted. The row totals are now assigned with ``.loc``.

- Fixed :meth:`~.LCA.get_impact_table` for empty ``'Stream'`` and ``'Other'`` categories: building the total row on an empty table assigned a float into a column pandas 3.0 had inferred as the string dtype, raising ``TypeError: Invalid value ... for dtype 'str'``. These categories now return a ``'No ...-related impacts.'`` message when empty (matching the existing ``'Construction'``/``'Transportation'`` behavior), which is also skipped when writing the report.

- Fixed :meth:`~.LCA.add_other_item`: when given a string ID with no matching :class:`~.ImpactItem`, the raised ``ValueError`` reported ``None`` (the failed-lookup result) instead of the requested ID. It now names the missing ID.

- Testing: added a ``conftest.py`` autouse fixture that resets the LCA registries and auto-ID ticket counters before each doctest, so the LCA doctests no longer leak state into one another (e.g., an indicator alias or auto-generated ID carrying over) when the modules run together.

- In EXPOsan, replaced the per-test ``clear_lca_registries()`` calls (deprecated) with an autouse ``conftest.py`` fixture that resets the LCA registries before each test.

- Documentation: expanded the tutorials with new sections on defining component groups (``Components.define_group``), unit specifications (``add_specification``), and inferring a :class:`~.System` from a list of units (``System.from_units``), plus notes on flowsheet retrieval, recycle convergence, and exporting results.

- Documentation: fixed dark-mode rendering of DataFrame tables and ``stderr`` (warning) output, standardized ``.diagram`` usage with a cross-reference between the System and Dynamic Simulation tutorials, and repaired a broken hyperlink.

- Documentation: concluded the topical-tutorial polish pass. Tutorial 12 (Anaerobic Digestion Model No. 1) gained a new §1 with hand-authored light/dark SVG pairs of the ADM1 reaction network and DAE structure (under ``docs/source/images/adm1/``), cross-tutorial connections to tutorials 10 and 11, and a data-grounded §3 discussion of the t = 10 d start-up transient (acidified plateau, biomass split into growing vs. declining groups). Tutorial 10 (Process) also gained a short Petersen-matrix paragraph that tutorial 12 references. Tutorial 13 (Process Modeling 101) was repositioned as the flowsheet-scale companion to tutorials 10-12: §1 was rewritten, cross-tutorial connections were added at the top of each §2 subsection (Component → tutorial 2, WasteStream → tutorial 3, Process → tutorial 10, SanUnit → tutorials 4/5/11, System → tutorials 6/11), §3 gained prose on the integrator-plus-recycle-convergence loop and on interpreting the steady-state output (SRT ≈ 10 d at the design point; nitrification / denitrification signatures; KLa-driven dissolved oxygen), and the prerequisites and learning objectives were updated to match. Tutorial 13's seven embedded image attachments were extracted to ``assets/tutorial_13/`` as snake_case files, clearing a Sphinx ``image.not_readable`` warning that the URL-encoded ``%20`` paths produced. Two ``Run ?su.X`` runtime-help prompts in tutorial 13 were replaced with documentation links into the autodoc-generated API pages (with the ``?su.X`` form retained as a parenthetical fallback).


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
