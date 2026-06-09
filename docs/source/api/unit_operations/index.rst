.. _unit_operations:

Unit Operations
===============

This page lists the unit operation modules available in ``QSDsan``. Unit operations are organized into three primary categories:

- ``bst`` units inherit from BioSTEAM unit operations and add QSDsan behavior while remaining steady-state. They are used the same way as native QSDsan units; for parameters and behavior specific to the BioSTEAM parent class, see the `BioSTEAM documentation <https://biosteam.readthedocs.io/en/latest/>`_.
- ``QSDsan static`` units are QSDsan-native unit operations without dynamic state equations.
- ``QSDsan dynamic`` units are QSDsan-native unit operations with dynamic state, derivative state, ODE, or algebraic-equation behavior.

.. toctree::
   :maxdepth: 1

   bst/index
   static/index
   dynamic/index


List of Unit Operations
-----------------------

#. ``Category`` indicates the unit operation category (BST, QSDsan static, or QSDsan dynamic).
#. ``Abstract`` indicates that design, costing, and environmental impact algorithms are not implemented.
#. ``API page`` links to the detailed documentation page. Some units appear both on a category page and on their implementation module page.

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
   :file: index.csv
   :header-rows: 1
   :class: unit-operation-table


