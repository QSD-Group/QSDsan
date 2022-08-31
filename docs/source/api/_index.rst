.. _api:

API
===

The UML below shows the structure and core Python classes implemented in ``QSDsan``. Each class is represented by a box containing the class name (bold, top part of the box) with select data (middle part of the box) and method (end with parentheses, bottom part of the box) attributes. The ``Component`` and ``WasteStream`` classes in blue are inherited from the ``Chemical`` and ``Stream`` classes in `Thermosteam <https://github.com/BioSTEAMDevelopmentGroup/thermosteam>`_ with the addition of wastewater-related attributes. The `Process` class in red enables dynamic simulation of Component objects’ transformation during kinetic processes (e.g., degradation of substrates). The ``SanUnit``, ``System``, ``TEA``, and ``Model`` classes in yellow are inherited from `BioSTEAM <https://github.com/BioSTEAMDevelopmentGroup/biosteam>`_ with added capacities for dynamic simulation and handling of construction inventories. Green boxes including ``mpactItem``, ``ImpactIndicator``, and ``LCA`` are implemented in ``QSDsan`` to enable LCA functionalities. 

.. figure:: https://lucid.app/publicSegments/view/c8de361f-7292-47e3-8870-d6f668691728/image.png
   :width: 100%

   Simplified unified modeling language (UML) diagram of ``QSDsan``


Module Links
------------

.. toctree::
   :maxdepth: 1

   Component
   Construction
   Equipment
   ImpactIndicator
   ImpactItem
   LCA
   Process
   SanUnit
   SimpleTEA
   streams
   Transportation
   equipments/_index
   processes/_index
   sanunits/_index
   stats
   utils/_index
