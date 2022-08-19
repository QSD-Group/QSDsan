QSDsan: Quantitative Sustainable Design for sanitation and resource recovery systems
====================================================================================

.. figure:: ./_static/logo_light_full.png
   :width: 50%
   :align: center

|

What is ``QSDsan``?
-------------------
``QSDsan`` is an open-source, community-led platform for the quantitative sustainable design (QSD) [1]_ of sanitation and resource recovery systems [2]_. It leverages existing platforms such as `BioSTEAM <https://biosteam.readthedocs.io>`_ [3]_  with enhanced features tailored to sanitation an resource recovery technologies. Through the integration with `DMsan <https://github.com/QSD-Group/DMsan>`_ (decision-making for sanitation and resource recovery systems), this platform can be used to guide the research, development, and deployment (RD&D) of early-stage technologies considering location-specific parameters and stakeholder priorities.


All systems developed with ``QSDsan`` are included in the package `EXPOsan <https://github.com/QSD-Group/EXPOsan>`_ - exposition of sanitation and resource recovery systems. Refer to the `Developed Systems <https://qsdsan.readthedocs.io/en/latest/Developed_Systems.html>`_ page for a list of the developed systems and the related literature.

.. toctree::
   :maxdepth: 1
   :caption: List of Systems

   Developed_Systems


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


Under the Hood
--------------
.. figure:: https://lucid.app/publicSegments/view/c8de361f-7292-47e3-8870-d6f668691728/image.png
   :width: 100%

   Simplified unified modeling language (UML) diagram of ``QSDsan``

|

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
   processes/_index
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
   - `Tori Morgan <https://qsdsan.readthedocs.io/en/latest/authors/Tori_Morgan.html>`_
   - `Hannah Lohman <https://qsdsan.readthedocs.io/en/latest/authors/Hannah_Lohman.html>`_


**Module development:** `developers <https://github.com/QSD-Group/QSDsan/graphs/contributors>`_ that have contributed to the repository.


**Funding support:**
   - `Jeremy Guest <mailto:jsguest@illinois.edu>`_
   - `Roland Cusick <mailto:rcusick@illinois.edu>`_
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
.. [1] Li, Y.; Trimmer, J.T.; Hand, S.; Zhang, X.; Chambers, K.G.; Lohman, H.A.C.; Shi, R.; Byrne, D.M.; Cook, S.M.; Guest, J.S. Quantitative Sustainable Design (QSD) for the Prioritization of Research, Development, and Deployment of Technologies: A Tutorial and Review, ChemRxiv, DOI:10.26434/chemrxiv-2022-rjqbn.

.. [2] Li, Y.; Zhang, X.; Morgan, V.L.; Lohman, H.A.C.; Rowles, L.S.; Mittal, S.; Kogler, A.; Cusick, R.D.; Tarpeh, W.A.; Guest, J.S. QSDsan: An integrated platform for quantitative sustainable design of sanitation and resource recovery systems. arXiv: `2203.06243 <https://arxiv.org/abs/2203.06243>`_ [cs.CY], March 7, 2022.

.. [3] Cortés-Peña, Y.; Kumar, D.; Singh, V.; Guest, J.S. BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040.


.. Links
.. _Yalin Li: https://qsdsan.readthedocs.io/en/latest/authors/Yalin_Li.html
.. _Joy Zhang: https://qsdsan.readthedocs.io/en/latest/authors/Joy_Zhang.html