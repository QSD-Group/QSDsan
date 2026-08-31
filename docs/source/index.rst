QSDsan: Quantitative Sustainable Design for Sanitation and Resource Recovery Systems
====================================================================================

What is ``QSDsan``?
-------------------
``QSDsan`` is an open-source, community-led platform for the quantitative sustainable design (QSD) [1]_ of sanitation and resource recovery systems [2]_. Built in Python, it integrates process modeling, system simulation, techno-economic analysis (TEA), and life cycle assessment (LCA) to support transparent, reproducible, and comprehensive evaluation of emerging technologies. Leveraging BioSTEAM [3]_, the platform provides modular and extensible tools to compare treatment configurations, evaluate resource recovery opportunities, and assess energy, cost, and sustainability tradeoffs. The platform has a range of built-in :ref:`unit operations <unit_operations>` and :ref:`systems <systems>`, including commonly used treatment processes and biokinetic models. These capabilities support the research, development, and deployment (RD&D) of early-stage technologies for advancing sustainable water management and circular resource recovery.

.. grid:: 2 4 4 4

   .. grid-item-card::  Tutorials
      :text-align: center
      :link: tutorials
      :link-type: ref
      :class-title: nav-card-title

      .. image:: images/icons/tutorials_icon_light.png
         :height: 100
         :align: center
         :class: only-light

      .. image:: images/icons/tutorials_icon_dark.png
         :height: 100
         :align: center
         :class: only-dark
          
   .. grid-item-card::  API
      :text-align: center
      :link: api
      :link-type: ref
      :class-title: nav-card-title

      .. image:: images/icons/api_icon_light.png
         :height: 100
         :align: center
         :class: only-light

      .. image:: images/icons/api_icon_dark.png
         :height: 100
         :align: center
         :class: only-dark

   .. grid-item-card::  Systems & Publications
      :text-align: center
      :link: systems
      :link-type: ref
      :class-title: nav-card-title

      .. image:: images/icons/systems_icon_light.png
         :height: 100
         :align: center
         :class: only-light

      .. image:: images/icons/systems_icon_dark.png
         :height: 100
         :align: center
         :class: only-dark

   .. grid-item-card::  Learning
      :text-align: center
      :link: learning
      :link-type: ref
      :class-title: nav-card-title

      .. image:: images/icons/learning_icon_light.png
         :height: 100
         :align: center
         :class: only-light

      .. image:: images/icons/learning_icon_dark.png
         :height: 100
         :align: center
         :class: only-dark


Installation
------------
``QSDsan`` is currently tested against Python 3.12, the minimum required to install it. The easiest way to install ``QSDsan`` is through ``pip`` in a command-line interface (e.g., a terminal or Anaconda Prompt):

.. code::

    pip install qsdsan

.. note::

    ``pip`` is the standard installer and works inside any environment, including a ``conda`` environment if you already use Anaconda. For a faster, drop-in alternative you can use `uv <https://docs.astral.sh/uv/>`_ (``uv pip install qsdsan``); it produces the same result. Installing into a dedicated virtual environment is recommended; the :doc:`Contributing guidelines <CONTRIBUTING>` give a step-by-step setup.

To upgrade an existing installation:

.. code::

    pip install -U qsdsan

To install a specific version, replace ``X.X.X`` with the version number:

.. code::

    pip install qsdsan==X.X.X

To install the latest GitHub version from the `main branch <https://github.com/QSD-Group/QSDsan>`_:

.. code::

    pip install git+https://github.com/QSD-Group/QSDsan.git

To install from another fork and/or branch, replace ``<USERNAME_OF_THE_FORK>`` and ``<BRANCH_NAME>``:

.. code::

    pip install git+https://github.com/<USERNAME_OF_THE_FORK>/QSDsan.git@<BRANCH_NAME>

You can also download the package from `PyPI <https://pypi.org/project/qsdsan/>`_.

For diagram generation, ``QSDsan`` uses Graphviz. If diagrams fail to render, install Graphviz following the `official Graphviz download instructions <https://graphviz.org/download/>`_ and see :ref:`graphviz-installation` for a quick check.


Join the Community
------------------
We would like to build an open and welcoming community, you can always post issues on our `GitHub homepage <https://github.com/QSD-Group/QSDsan/issues>`_ or contact any of the Quantitative Sustainable Design Group members. We are always excited to have new members in our team.

You can receive updates by joining `QSDsan's Google Group <https://groups.google.com/g/qsdsan>`_. To be added, log in to your Google account and click the ``Ask to join group`` button. If you run into any problems or would like to use a non-Gmail account, please `contact us <mailto:quantitative.sustainable.design@gmail.com>`_.

If you would like to contribute, please follow our `Contributing Guidelines`_ and the `Contributor Covenant <https://www.contributor-covenant.org/>`_, thank you for making ``QSDsan`` better!


``QSDsan`` is and will stay open source under University of Illinois/NCSA Open Source License. Any third-party packages copied from ``QSDsan`` must be strictly open-source (not copy-left nor open-access). Please refer to the `license <https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt>`_ page for details.


.. Hidden TOCs for navigation bar
.. toctree::
   :maxdepth: 1
   :hidden:

   tutorials/index


.. toctree::
   :maxdepth: 1
   :hidden:

   api/index


.. toctree::
   :maxdepth: 1
   :hidden:

   systems/index


.. toctree::
   :maxdepth: 1
   :hidden:

   learning/index


.. toctree::
   :maxdepth: 2
   :hidden:

   faq/index


.. toctree::
   :maxdepth: 1
   :hidden:

   app/index


.. toctree::
   :maxdepth: 1
   :hidden:

   CONTRIBUTING


.. toctree::
   :maxdepth: 1
   :hidden:

   CHANGELOG


References
----------
.. [1] Li, Y.; Trimmer, J.T.; Hand, S.; Zhang, X.; Chambers, K.G.; Lohman, H.A.C.; Shi, R.; Byrne, D.M.; Cook, S.M.; Guest, J.S. Quantitative Sustainable Design (QSD): A Methodology for the Prioritization of Research, Development, and Deployment of Technologies. (Tutorial Review) Environ. Sci.: Water Res. Technol. 2022, 8 (11), 2439–2465. https://doi.org/10.1039/D2EW00431C.

.. [2] Li, Y.; Zhang, X.; Morgan, V.L.; Lohman, H.A.C.; Rowles, L.S.; Mittal, S.; Kogler, A.; Cusick, R.D.; Tarpeh, W.A.; Guest, J.S. QSDsan: An integrated platform for quantitative sustainable design of sanitation and resource recovery systems. Environ. Sci.: Water Res. Technol. 2022, 8 (10), 2289-2303. https://doi.org/10.1039/d2ew00455k.

.. [3] Cortes-Peña, Y.; Kumar, D.; Singh, V.; Guest, J.S. BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040.

.. Links
.. _Contributing Guidelines: CONTRIBUTING.html#contributing-guidelines
