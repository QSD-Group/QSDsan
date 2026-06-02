.. _biosteam_api:

BioSTEAM API in QSDsan
======================

``QSDsan`` is built on top of `BioSTEAM <https://github.com/BioSTEAMDevelopmentGroup/biosteam>`_
and BioSTEAM's thermodynamic library `Thermosteam <https://github.com/BioSTEAMDevelopmentGroup/thermosteam>`_. QSDsan re-exports the parts of those libraries that a typical user needs.

Everything below is reachable through ``import qsdsan as qs``.

Streams, chemicals, and thermo configuration
--------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Name
     - Notes
   * - ``qs.Chemical``
     - Thermosteam ``Chemical`` (see also ``qs.Component``).
   * - ``qs.Chemicals``
     - Thermosteam ``Chemicals`` (see also ``qs.Components``).
   * - ``qs.CompiledChemicals``
     - Compiled chemical package.
   * - ``qs.Stream``
     - BioSTEAM ``Stream`` (see also ``qs.SanStream`` / ``qs.WasteStream``).
   * - ``qs.MultiStream``
     - Multiphase stream.
   * - ``qs.MissingStream``
     - Placeholder used during unit initialization.
   * - ``qs.set_thermo``
     - Set the thermodynamic property package.
   * - ``qs.get_thermo``
     - Get the active thermo package.
   * - ``qs.get_components``
     - Get the active components.
   * - ``qs.settings``
     - Process-wide settings object (e.g. ``qs.settings.thermo.ideal()``).
   * - ``qs.preferences``
     - Display/diagram defaults (dark mode, stream labels, graphviz format, etc.).
   * - ``qs.Thermo``
     - Build a custom thermo package.

Units, systems, and TEA
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Name
     - Notes
   * - ``qs.Unit``
     - Base BioSTEAM unit (see also ``qs.SanUnit``).
   * - ``qs.System``
     - Process system (``qsdsan`` subclass).
   * - ``qs.AgileSystem``
     - System evaluated across multiple operating scenarios.
   * - ``qs.Facility``
     - Base class for facility units.
   * - ``qs.get_OSBL``
     - Collect outside-battery-limit (facility) units of a system.
   * - ``qs.TEA``
     - Techno-economic analysis (``qsdsan`` subclass of BioSTEAM ``TEA``).
   * - ``qs.Model``
     - Uncertainty/sensitivity model.
   * - ``qs.Metric`` / ``qs.Parameter``
     - Model metrics and parameters.
   * - ``qs.Scope``
     - Dynamic-simulation scope.
   * - ``qs.report``
     - Reporting helpers (e.g. ``qs.report.voc_table``).

Utilities and costing globals
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Name
     - Notes
   * - ``qs.HeatUtility``
     - Heating/cooling agents (e.g. ``qs.HeatUtility.get_agent('cooling_water')``).
   * - ``qs.PowerUtility``
     - Electricity; set the price via ``qs.PowerUtility.price = ...``.
   * - ``qs.UtilityAgent``
     - Define a custom heating/cooling agent.
   * - ``qs.default``
     - Reset utilities, the CEPCI, and/or the active flowsheet to ``qsdsan``'s
       defaults; choose what to reset with ``utilities`` / ``CEPCI`` / ``flowsheet``
       (e.g. ``qs.default(flowsheet=False)``).
   * - ``qs.default_utilities``
     - Reset utilities to their defaults.
   * - ``qs.stream_utility_prices``
     - Mapping of stream-based utility prices.
   * - ``qs.CEPCI``
     - Chemical Engineering Plant Cost Index; a live, settable view of BioSTEAM's
       ``CE`` (e.g. ``qs.CEPCI = qs.CEPCI_by_year[2023]``).
   * - ``qs.CEPCI_by_year``
     - Lookup table of CEPCI values by year.

.. note::

   ``qs.CEPCI`` is the only cost/utility global that needs a special accessor: it
   is a settable property that writes through to BioSTEAM's ``CE``. The others
   above are the same objects/classes as in BioSTEAM, so setting an attribute on
   them (for example ``qs.PowerUtility.price`` or ``qs.settings.thermo``) updates
   the shared, process-wide state directly.

Reactions
---------

The full reaction API is re-exported: ``qs.Reaction`` (``qs.Rxn``),
``qs.ParallelReaction`` (``qs.PRxn``), ``qs.SeriesReaction`` (``qs.SRxn``),
``qs.ReactionItem`` (``qs.RxnI``), ``qs.ReactionSet`` (``qs.RxnS``), and
``qs.ReactionSystem`` (``qs.RxnSys``).

BioSTEAM-inherited units
------------------------

``qsdsan`` wraps many BioSTEAM units (mixers, splitters, pumps, tanks, heat
exchangers, distillation, flash) under :mod:`qsdsan.unit_operations.bst`. Several
BioSTEAM units that ``qsdsan`` does not customize are also surfaced there as raw
re-exports so systems can be built without ``import biosteam``::

    from qsdsan.unit_operations import bst as su

    su.IsenthalpicValve, su.Stripper, su.MolecularSieve, su.BatchBioreactor,
    su.VacuumSystem, su.Boiler, su.BoilerTurbogenerator, su.ChilledWaterPackage,
    su.CoolingTower, su.SolidsCentrifuge

See the :ref:`unit operations <unit_operations>` reference for the full list.

Helper functions in :mod:`qsdsan.utils`
---------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Name
     - Notes
   * - ``qs.utils.rho_to_V``
     - Convert density to molar volume.
   * - ``qs.utils.V_to_rho``
     - Convert molar volume to density.
   * - ``qs.utils.cost``
     - The ``@cost`` decorator for adding cost/design to a unit.
   * - ``qs.utils.var_columns``
     - Column labels for a ``Model`` parameter/metric DataFrame.
   * - ``qs.utils.var_indices``
     - Index labels for ``Model`` parameters/metrics.

What is intentionally not re-exported
-------------------------------------

A small set of lower-level, unit-authoring tools are deliberately left in
``biosteam`` (import them directly when you need them):

* ``biosteam.units.design_tools`` (e.g. ``PressureVessel``, ``flash_vessel_design``,
  ``size_batch``, ``material_densities_lb_per_ft3``)
* ``biosteam.exceptions`` (e.g. ``DesignError``, ``DesignWarning``, ``bounds_warning``)

These are needed only when building custom unit operations, and they fall outside
the common surface ``qsdsan`` aims to cover.
