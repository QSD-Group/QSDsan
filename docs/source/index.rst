QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
====================================================================================

.. figure:: ./_static/logo_light_full.png
   :width: 500
   :align: center

|

What is ``QSDsan``?
-------------------
``QSDsan`` is an open-source, community-led platform for the quantitative sustainable design (QSD) of sanitation and resource recovery systems. It is one of a series of platforms that are being developed for the execution of QSD - a methodology for the research, design, and deployment of technologies and inform decision-making. [1]_ It leverages the structure and modules developed in the `BioSTEAM <https://github.com/BioSTEAMDevelopmentGroup/biosteam>`_ platform [2]_ with additional functions tailored to sanitation processes.

As an open-source and impact-driven platform, QSDsan aims to identify configuration combinations, systematically probe interdependencies across technologies, and identify key sensitivities to contextual assumptions through the use of quantitative sustainable design methods (techno-economic analysis and life cycle assessment and under uncertainty). 

All systems developed with ``QSDsan`` are included in the package `EXPOsan <https://github.com/QSD-Group/EXPOsan>`_ - exposition of sanitation and resource recovery systems.

Additionally, another package, `DMsan <https://github.com/QSD-Group/DMsan>`_ (decision-making for sanitation and resource recovery systems), is being developed for decision-making among multiple dimensions of sustainability with consideration of location-specific contextual parameters.


To capitalize, or not to capitalize?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
That is a good question. We prefer to use the capitalized version (e.g., ``QSDsan`` instead of ``qsdsan``) because:

- It refers to the platform, not just the core package (i.e., it includes the entire ecosystem that supports the core package).
- We style the name to convey the name's meaning (e.g., the "QSD" part stands for "quantitative sustainable design").

But names of the actual packages are all in lower cases per `PEP-8 <https://www.python.org/dev/peps/pep-0008/#package-and-module-names>`_:

   *Modules should have short, all-lowercase names. Underscores can be used in the module name if it improves readability. Python packages should also have short, all-lowercase names, although the use of underscores is discouraged.*


Installation
------------
The easiest way is through ``pip``, in command-line interface (Anaconda prompt, terminal):

.. code::

    pip install qsdsan

If you need to update:

.. code::

    pip install -U qsdsan

Or for a specific version (replace X.X.X with the version number):

.. code::

    pip install qsdsan==X.X.X

If you want to install the latest GitHub version at the `main branch <https://github.com/qsd-group/qsdsan>`_:

.. code::

    pip install git+https://github.com/QSD-Group/QSDsan.git


You can also download the package from `PyPI <https://pypi.org/project/qsdsan/>`_.

Note that development of this package is currently under initial stage with limited backward compatibility, please feel free to `submit an issue <https://github.com/QSD-Group/QSDsan/issues>`_ for any questions regarding package upgrading.

If you are a developer and want to contribute to ``QSDsan``, please follow the steps in the `Contributing to QSDsan <https://qsdsan.readthedocs.io/en/latest/CONTRIBUTING.html>`_ section of the documentation to clone the repository.


Getting Started
---------------
To help you familiarize with ``QSDsan``, we include two types of tutorials - topical tutorials that cover on one or multiple classes of ``QSDsan``, or design tutorials that show you how to translate a practical design problem into coding in QSDsan. Follow them to get started!

All tutorials are written using Jupyter Notebook, you can run your own Jupyter environment, or you can visit `this page <https://mybinder.org/v2/gh/QSD-Group/QSDsan/main?filepath=%2Fdocs%2Fsource%2Ftutorials>`_ to run the Jupyter environment in your browser.

For each of these tutorials and examples, we are also recording videos where one of the QSD group members will go through the material step-by-step. We are gradually releasing these videos on our `YouTube channel <https://www.youtube.com/channel/UC8fyVeo9xf10KeuZ_4vC_GA>`_ so subscribe to receive updates!

.. toctree::
   :maxdepth: 1
   :caption: Topical Tutorials

   tutorials/0_Quick_Overview
   tutorials/1_Helpful_Basics
   tutorials/2_Component
   tutorials/3_WasteStream
   tutorials/4_SanUnit_basic
   tutorials/5_SanUnit_advanced
   tutorials/6_System
   tutorials/7_TEA
   tutorials/8_LCA
   tutorials/9_Uncertainty_and_Sensitivity_Analyses

..
   .. toctree::
   :maxdepth: 1
   :caption: Design Tutorials

   tutorials/12_Chlorination


How does ``QSDsan`` work?
-------------------------
.. figure:: https://lucid.app/publicSegments/view/c8de361f-7292-47e3-8870-d6f668691728/image.png

   Simplified unified modeling language (UML) diagram of ``QSDsan``

The core package ``QSDsan`` follows the structure of `biosteam <https://github.com/BioSTEAMDevelopmentGroup/biosteam>`_, a rapid and agile package for the design, simulation, and techno-economic analysis of biorefineries under uncertainty, but ``QSDsan`` is enhanced with features geared toward quantitative sustainable design of sanitation systems.

Work is also in progress to develop the package ``DMsan`` and connect ``QSDsan`` to ``DMsan`` for decision-making under different contexts.

The above Unified Modeling Language (UML) diagram of the package shows the relationship between ``QSDsan`` and its dependencies.

In particular, ``QSDsan`` introduces:

- :class:`~.Component`, a subclass of :class:`thermosteam.Chemical`, instance of this class does not necessarily corresponds to a specific chemical, but represents commonly used/modeled component such as biodegradable colloidal substrate.
- :class:`~.SanStream` and :class:`~.WasteStream`. :class:`~.SanStream` is a sublcass of :class:`thermosteam.Stream` that has an additional attribute ``impact_item`` for life cycle assessment. :class:`~.WasteStream` is a subclass of :class:`~.SanStream` with additional composite properties such as chemical oxygen demand (COD) that are widely used in sanitation systems.
- :class:`~.Process`, a new class that describes a certain biological, chemical, or physical process in a unit operation, it is similar in concept with :class:`thermosteam.Reaction`, but has unique features and utilities.


.. toctree::
   :maxdepth: 1
   :caption: API

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
   processes/processes
   sanunits/_index
   stats
   utils/_index


About the developers
--------------------
Development and maintenance of the platform is supported by the Quantitative Sustainable Design Group led by members of the `Guest Group <http://engineeringforsustainability.com/>`_ at the `University of Illinois Urbana-Champaign (UIUC) <https://illinois.edu/>`_, as well as `other developers <https://github.com/QSD-Group/QSDsan/graphs/contributors>`_ that have contributed to the repository.


.. toctree::
   :maxdepth: 1
   :caption: Core Developers

   authors/Yalin_Li
   authors/Joy_Zhang
   authors/Tori_Morgan
   authors/Hannah_Lohman
   authors/Stetson_Rowles


Roles
^^^^^
**Lead developers:**
   - `Yalin Li`_ (current maintainer)
   - `Joy Zhang`_


**Tutorials and videos:**
   - `Yalin Li`_
   - `Joy Zhang`_
   - `Tori Morgan <https://qsdsan.readthedocs.io/en/beta/authors/Tori_Morgan.html>`_
   - `Hannah Lohman <https://qsdsan.readthedocs.io/en/beta/authors/Hannah_Lohman.html>`_


**Module development:** `developers <https://github.com/QSD-Group/QSDsan/graphs/contributors>`_ that have contributed to the repository.


**Funding support:**
   - `Jeremy Guest <mailto:jsguest@illinois.edu>`_
   - `Ro Cusick <mailto:rcusick@illinois.edu>`_
   - `William Tarpeh <mailto:wtarpeh@stanford.edu>`_


**Special acknowledgement:**
   - Yoel Cortés-Peña for helping many of the ``QSDsan`` members get started on Python and package development.


Join the community
------------------
We would like to build an open and welcoming community, you can always post issues on our `GitHub homepage <https://github.com/QSD-Group/QSDsan/issues>`_ or contact any of the Quantitative Sustainable Design Group memebers. We are always excited to have new members in our team.

If you would like to contribute, please follow our contribution guide, thank you for making ``QSDsan`` better!

.. toctree::
   :maxdepth: 1
   :caption: Contributions

   CODE_OF_CONDUCT
   CONTRIBUTING
   for_developers/Tutorial_Template


``QSDsan`` is and will stay open source under University of Illinois/NCSA Open Source License. Any third-party packages copied from ``QSDsan`` must be strictly open-source (not copy-left nor open-access). Please refer to `LICENSE <https://github.com/QSD-Group/QSDsan/blob/master/LICENSE.txt>`_ and `CONTRIBUTION <https://qsdsan.readthedocs.io/en/latest/CONTRIBUTING.html>`_ for details.

If you would like to receive news related to the QSDsan platform, you can subscribe to email updates using `this form <https://groups.webservices.illinois.edu/subscribe/154591>`_ (don't worry, you will be able to unsubscribe :)). Thank you in advance for your interest!


.. toctree::
   :caption: Events
   :maxdepth: 1

   Events.rst

We will keep the calendar up-to-date as we organize more events (office hours, workshops, etc.), click on the events in the calendar to see the details (including meeting links).

.. raw:: html

    <embed>
        <iframe src="https://calendar.google.com/calendar/embed?src=ep1au561lj8knfumpcd2a7ml08%40group.calendar.google.com&ctz=America%2FChicago" style="border: 0" width="100%" height="600" scrolling="no"></iframe>
    </embed>

|

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
.. [1] Li, Y.; Trimmer, J. T.; Hand, S.; Zhang, X.; Chambers, K. G.; Lohman, H. A. C.; Shi, R.; Byrne, D. M.; Cook, S. M.; Guest, J. S. Quantitative Sustainable Design (QSD): A Methodology for the Prioritization of Research, Development, and Deployment of Technologies. Submitted to Environmental Science & Technology, 2022.

.. [2] Cortés-Peña, Y.; Kumar, D.; Singh, V.; Guest, J. S. BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040.


.. Links
.. _Yalin Li: https://qsdsan.readthedocs.io/en/beta/authors/Yalin_Li.html
.. _Joy Zhang: https://qsdsan.readthedocs.io/en/beta/authors/Joy_Zhang.html