.. _unit_operations:

Unit Operations
===============

This page lists the unit operation modules available in ``QSDsan``. Unit operations are organized into three primary categories:

- ``BST`` units inherit from BioSTEAM unit operations and add QSDsan behavior while remaining steady-state.
- ``QSDsan static`` units are QSDsan-native unit operations without dynamic state equations.
- ``QSDsan dynamic`` units are QSDsan-native unit operations with dynamic state, derivative state, ODE, or algebraic-equation behavior.

``Abstract`` is a unit property, not a separate category. Abstract units are simulation helpers without design, costing, or environmental impact algorithms, and they may appear in any of the three primary categories.

The catalog keeps existing documentation paths stable:

- Uppercase page names usually document one major class, such as ``Excretion`` or ``Flash``.
- Lowercase page names usually document a module containing multiple related classes or helper functions, such as ``tank`` or ``heat_exchanging``.
- System-level collection pages group units commonly used together in a specific sanitation system design.

List of Unit Operations
-----------------------
.. raw:: html

   <div class="unit-operation-filters" aria-label="Unit operation table filters">
     <label for="unit-operation-category-filter">Category</label>
     <select id="unit-operation-category-filter">
       <option value="">All</option>
       <option value="BST">BST</option>
       <option value="QSDsan static">QSDsan static</option>
       <option value="QSDsan dynamic">QSDsan dynamic</option>
     </select>
     <label for="unit-operation-abstract-filter">Abstract</label>
     <select id="unit-operation-abstract-filter">
       <option value="">All</option>
       <option value="Yes">Yes</option>
       <option value="No">No</option>
     </select>
     <span id="unit-operation-filter-count" class="unit-operation-filter-count"></span>
   </div>

.. csv-table::
   :file: _index.csv
   :header-rows: 1
   :class: unit-operation-table

#. ``Abstract`` indicates that design, costing, and environmental impact algorithms are not implemented.
#. ``API page`` links to the detailed documentation page. Some units appear both on a category page and on their implementation module page.


.. toctree::
   :maxdepth: 1
   :hidden:

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
