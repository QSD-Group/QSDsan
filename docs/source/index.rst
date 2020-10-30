SanPlorer: Exploring Sustainable Sanitation Technologies
========================================================
.. figure:: ./images/various_configurations.png
:caption: blah...

What is SanPlorer?
------------------
SanPlorer (Sanitation Explorer) is a package for the sustainable design of non-sewered sanitation technologies that are affordable and appropriate in low-income, resource-limited settings leveraging the structure and modules developed in BioSTEAM [1]_. As an open-source and impact-driven platform, SanPlorer aims to identify configuration combinations, systematically probe interdependencies across technologies, and identify key sensitivities to contextual assumptions through the use of quantitative sustainable design methods (techno-economic analysis and life cycle assessment under uncertainty).

# Here also mention SanExpo


Getting Started
---------------
Want to try out SanPlorer? Follow the tutorial to get started!

.. toctree::
   :maxdepth: 1
   :caption: Tutorial:

   tutorial/Installation
   tutorial/Make_a_simple_system
   tutorial/Define_your_own_Component
   tutorial/Create_a_new_Unit
   tutorial/TEA_and_LCA


How does SanPlorer work?
------------------------
.. figure:: ./images/sanitation_UML.png

SanPlorer follows the structure of `BioSTEAM <https://github.com/BioSTEAMDevelopmentGroup/biosteam>`_, a fast and flexible package for the design, simulation, and techno-economic analysis of biorefineries under uncertainty, but SanPlorer is enhanced with features geared toward quantitative sustainable design of sanitation systems.

The above Unified Modeling Language (UML) diagram of the package shows the relationship between SanPlorer and its dependencies `BioSTEAM <https://github.com/BioSTEAMDevelopmentGroup/biosteam>`_, `thermosteam <https://github.com/BioSTEAMDevelopmentGroup/thermosteam>`_, and `chemicals <https://github.com/CalebBell/chemicals>`_.

In particular, SanPlorer introduces:
 - ``Component``, a subclass of ``Chemical`` in thermosteam, instance of this class does not necessarily corresponds to a specific chemical, but represents commonly used/modeled component such as biodegradable colloidal substrate.

 - ``WasteStream``, a sublcass of ``Stream`` in thermosteam, instance of this class has additional composite properties such as chemical oxygen demand (COD) that are widely used in sanitation systems.

 - ``Process``, a new class that describes a certain biological, chemical, or physical process in a unit operation, it has some similarities with the ``reaction`` class in thermosteam, but has unique features and utilities.


.. toctree::
   :maxdepth: 1
   :caption: Inside SanPlorer:

   Component
   Components
   WasteStream
   SanUnit
   units


More resources
--------------
To get the full value of SanPlorer, we highly recommend reading through the documents of these packages:
 - `biosteam <https://biosteam.readthedocs.io/en/latest/index.html>`_
 - `thermosteam <https://thermosteam.readthedocs.io/en/latest/index.html>`_
 - `chemicals <https://chemicals.readthedocs.io/en/latest/>`_


About the authors
-----------------
Development and maintenance of the package is supported by the SanPlorer Development Group led by members of the `Guest Group <http://engineeringforsustainability.com/>`_ at the `University of Illinois Urbana-Champaign (UIUC) <https://illinois.edu/>`_. Core contributors include:

**Code development**
 - `Yalin Li <zoe.yalin.li@gmail.com>`_
 - Joy Cheung

**Unit design**
 - `Yalin Li <zoe.yalin.li@gmail.com>`_
 - Joy Cheung
 - Stetson Rowles

 **Project conception & funding support**
 - `Jeremy Guest <jsguest@illinois.edu>`_


Join the community
------------------
We would like to build an open and welcoming community, you can always post issues on our GitHub homepage or contact any of the SanPlorer Development Group memebers. If you would like to contribute, please follow our contribution guide:

.. toctree::
   :maxdepth: 1
   :caption: How to contribute?:

   Contributing
   Code_of_Conduct
   License


References
----------
.. [1] Cortes-Peña, Y.; Kumar, D.; Singh, V.; Guest, J. S. BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020. https://doi.org/10.1021/acssuschemeng.9b07040.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
