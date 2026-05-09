.. _unit_operations:

Unit Operations
===============

This page lists the unit operation models available in ``QSDsan``. Use the behavior-based namespaces for the clearest public import paths, and use the detailed API pages below when you need implementation-level documentation for a specific class or module.

The catalog keeps existing documentation paths stable:

- Uppercase page names usually document one major class, such as ``Excretion`` or ``Flash``.
- Lowercase page names usually document a module containing multiple related classes or helper functions, such as ``tank`` or ``heat_exchanging``.
- System-level collection pages group units commonly used together in a specific sanitation system design.

List of Units
-------------
.. csv-table::
   :file: _index.csv
   :header-rows: 1

#. Dynamic units support dynamic simulation.
#. Abstract units do not have design and cost algorithms, they are for simulation purpose only.


Behavior-Based Namespaces
-------------------------
.. toctree::
   :maxdepth: 1

   bst
   static
   dynamic


System-Level Unit Collections
-----------------------------
.. toctree::
   :maxdepth: 1

   biogenic_refinery
   eco_san
   reclaimer


Single-Class Pages
------------------
.. toctree::
   :maxdepth: 1

   ActivatedSludgeProcess
   CropApplication
   DynamicInfluent
   ElectrochemicalCell
   Excretion
   Flash
   InternalCirculationRx
   Lagoon
   MembraneDistillation
   MembraneGasExtraction
   PolishingFilter
   Reactor
   Screening
   Sedimentation
   SepticTank
   SludgePasteurization
   Trucking


Module Pages
------------
.. toctree::
   :maxdepth: 1

   abstract
   anaerobic_reactor
   clarifier
   combustion
   compressor
   distillation
   facilities
   heat_exchanging
   hydroprocessing
   hydrothermal
   junction
   membrane_bioreactor
   non_reactive
   pumping
   sludge_thickening
   sludge_treatment
   suspended_growth_bioreactor
   tank
   toilet
   treatment_bed
