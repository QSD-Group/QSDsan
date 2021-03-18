QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
====================================================================================
.. figure:: ./images/various_configurations.png
   
   Collage of sanitation units inclucded in QSDsan


What is ``QSDsan``?
-------------------

``QSDsan`` is a package for the quantitative sustainable design of sanitation and resource recovery systems leveraging the structure and modules developed in ``BioSTEAM`` [1]_. As an open-source and impact-driven platform, QSDsan aims to identify configuration combinations, systematically probe interdependencies across technologies, and identify key sensitivities to contextual assumptions through the use of quantitative sustainable design methods (techno-economic analysis and life cycle assessment and under uncertainty). 

All systems developed with ``QSDsan`` will be included in another repository in the future.


Getting Started
---------------
First install the package at `PyPI <https://pypi.org/>`_. If you use pip:

.. code:: bash

    pip install qsdsan


Follow the tutorial to get started!

.. toctree::
   :maxdepth: 1
   :caption: Tutorial

   tutorials/Component_and_WasteStream
   tutorials/SanUnit_and_System
   tutorials/TEA_and_LCA


How does ``QSDsan`` work?
-------------------------
.. figure:: https://lucid.app/publicSegments/view/c8de361f-7292-47e3-8870-d6f668691728/image.png

   Simplified unified modeling language (UML) diagram of ``QSDsan``

``QSDsan`` follows the structure of `BioSTEAM <https://github.com/BioSTEAMDevelopmentGroup/biosteam>`_, a fast and flexible package for the design, simulation, and techno-economic analysis of biorefineries under uncertainty, but ``QSDsan`` is enhanced with features geared toward quantitative sustainable design of sanitation systems.

The above Unified Modeling Language (UML) diagram of the package shows the relationship between ``QSDsan`` and its dependencies `BioSTEAM <https://github.com/BioSTEAMDevelopmentGroup/biosteam>`_ and `Thermosteam <https://github.com/BioSTEAMDevelopmentGroup/thermosteam>`_.

In particular, ``QSDsan`` introduces:

- :class:`~.Component`, a subclass of :class:`thermosteam.Chemical`, instance of this class does not necessarily corresponds to a specific chemical, but represents commonly used/modeled component such as biodegradable colloidal substrate.
- :class:`~.WasteStream`, a sublcass of :class:`thermosteam.Stream`, instance of this class has additional composite properties such as chemical oxygen demand (COD) that are widely used in sanitation systems.
- :class:`~.Process`, a new class that describes a certain biological, chemical, or physical process in a unit operation, it is similar in concept with :class:`thermosteam.Reaction`, but has unique features and utilities.


.. toctree::
   :maxdepth: 1
   :caption: API

   Component
   Components
   Construction
   Equipment
   ImpactIndicator
   ImpactItem
   LCA
   Process
   SanUnit
   SimpleTEA
   Transportation
   WasteStream
   sanunits/sanunits
   stats


About the developers
--------------------
Development and maintenance of the package is supported by the Quantitative Sustainable Design Group led by members of the `Guest Group <http://engineeringforsustainability.com/>`_ at the `University of Illinois Urbana-Champaign (UIUC) <https://illinois.edu/>`_. Core contributors are listed below, please refer to the `author page <https://qsdsan-beta.readthedocs.io/en/latest/for_developers/AUTHORS.html>`_ for the full list of authors.

**Lead developers:**
   - `Yalin Li <zoe.yalin.li@gmail.com>`_
   - `Joy Cheung <joycheung1994@gmail.com>`_
   - `Stetson Rowles <lsr@illinois.edu>`_

**Project conception & funding support:**
   - `Jeremy Guest <jsguest@illinois.edu>`_

**Special acknowledgement:**
   - Yoel Cortés-Peña for helping many of the ``QSDsan`` members get started on Python and package development.


Join the community
------------------
We would like to build an open and welcoming community, you can always post issues on our `GitHub homepage <https://github.com/QSD-Group/QSDsan/issues>`_ or contact any of the Quantitative Sustainable Design Group memebers. We are always excited to have new members in our team.

If you would like to contribute, please follow our contribution guide, thank you for making ``QSDsan`` better!

.. toctree::
   :maxdepth: 1
   :caption: Contributions

   AUTHORS
   CODE_OF_CONDUCT
   CONTRIBUTING
   for_developers/tutorial_template


``QSDsan`` is and will stay open source under University of Illinois/NCSA Open Source License. Any third-party packages copied from ``QSDsan`` must be strictly open-source (not copy-left nor open-access). Please refer to `LICENSE <https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt>`_ and `CONTRIBUTION <https://qsdsan.readthedocs.io/en/latest/CONTRIBUTING.html>`_ for details.


More resources
--------------
.. toctree::
   :maxdepth: 1
   :caption: FAQ

   FAQ.rst

Additionally, to get the full value of ``QSDsan``, we highly recommend reading through the documents of these packages:

- `biosteam docs <https://biosteam.readthedocs.io/en/latest/index.html>`_
- `thermosteam docs <https://thermosteam.readthedocs.io/en/latest/index.html>`_
- `chemicals docs <https://chemicals.readthedocs.io/en/latest/>`_ (the thermodynamic property package for thermosteam)


.. toctree::
   :maxdepth: 1
   :caption: What's new

   CHANGELOG


References
----------
.. [1] Cortés-Peña, Y.; Kumar, D.; Singh, V.; Guest, J. S. BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020. https://doi.org/10.1021/acssuschemeng.9b07040.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
