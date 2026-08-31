.. _unit_operations:

Unit Operations
===============

This page lists the unit operation modules available in ``QSDsan``. Unit operations are organized into three primary categories:

- ``bst`` units are thin wrappers that make a real BioSTEAM unit operation usable in QSDsan (same core BioSTEAM simulate/design/cost logic, plus QSDsan stream compatibility). They are used the same way as native QSDsan units; for parameters and behavior specific to the BioSTEAM parent class, see the `BioSTEAM documentation <https://biosteam.readthedocs.io/en/latest/>`_.
- ``static`` units have no working dynamic-state contract -- running them with ``isdynamic=True`` would error.
- ``dynamic`` units can run in dynamic mode -- ``isdynamic=True`` works, whether that capability is native to the unit or inherited from a BioSTEAM-wrapper base class (e.g. ``Mixer``, ``Splitter``, ``Pump``).

.. toctree::
   :maxdepth: 1

   bst/index
   static/index
   dynamic/index


List of Unit Operations
-----------------------

#. ``Category`` indicates the unit operation category (bst, static, or dynamic).
#. ``Abstract`` indicates that design, costing, and environmental impact algorithms are not implemented.
#. ``API page`` links to the detailed documentation page. Some units appear both on a category page and on their implementation module page.

.. raw:: html

   <div class="unit-operation-filters" aria-label="Unit operation table filters">
     <label for="unit-operation-category-filter">Category</label>
     <select id="unit-operation-category-filter">
       <option value="">All</option>
       <option value="bst">bst</option>
       <option value="static">static</option>
       <option value="dynamic">dynamic</option>
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


