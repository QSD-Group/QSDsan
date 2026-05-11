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

   bst/index
   static/index
   dynamic/index


System-Level Unit Collections
-----------------------------
.. toctree::
   :maxdepth: 1

   static/biogenic_refinery
   static/eco_san
   static/reclaimer


BST Units
---------
.. toctree::
   :maxdepth: 1

   bst/Flash
   bst/abstract
   bst/compressor
   bst/distillation
   bst/facilities
   bst/heat_exchanging
   bst/pumping
   bst/tank


Static Units
------------
.. toctree::
   :maxdepth: 1

   static/ActivatedSludgeProcess
   static/abstract
   static/anaerobic_reactor
   static/clarifier
   static/combustion
   static/CropApplication
   static/ElectrochemicalCell
   static/Excretion
   static/hydroprocessing
   static/hydrothermal
   static/InternalCirculationRx
   static/Lagoon
   static/MembraneDistillation
   static/MembraneGasExtraction
   static/non_reactive
   static/PolishingFilter
   static/pumping
   static/Reactor
   static/Screening
   static/Sedimentation
   static/SepticTank
   static/sludge_thickening
   static/sludge_treatment
   static/SludgePasteurization
   static/toilet
   static/treatment_bed
   static/Trucking


Dynamic Units
-------------
.. toctree::
   :maxdepth: 1

   dynamic/abstract
   dynamic/anaerobic_reactor
   dynamic/DynamicInfluent
   dynamic/junction
   dynamic/membrane_bioreactor
   dynamic/pumping
   dynamic/suspended_growth_bioreactor
