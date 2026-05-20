====================================================================================
QSDsan: Quantitative Sustainable Design for Sanitation and Resource Recovery Systems
====================================================================================

.. License
.. image:: https://img.shields.io/pypi/l/qsdsan?color=blue&logo=UIUC&style=flat
   :target: https://github.com/QSD-Group/QSDsan/blob/main/LICENSE.txt

.. Tested Python version
.. image:: https://img.shields.io/pypi/pyversions/qsdsan?style=flat
   :target: https://pypi.python.org/pypi/qsdsan

.. PyPI version
.. image:: https://img.shields.io/pypi/v/qsdsan?style=flat&color=blue
   :target: https://pypi.org/project/qsdsan

.. Zenodo release DOI
.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.20256569.svg
   :target: https://doi.org/10.5281/zenodo.20256569

.. Paper DOI
.. image:: https://img.shields.io/badge/qsdsan--paper-10.1039%2Fd2ew00455k-blue?style=flat
   :target: https://doi.org/10.1039/d2ew00455k

.. Documentation build
.. image:: https://readthedocs.org/projects/qsdsan/badge/?version=latest
   :target: https://qsdsan.readthedocs.io/en/latest

.. GitHub test and coverage of the main branch
.. image:: https://github.com/QSD-Group/QSDsan/actions/workflows/build-coverage.yml/badge.svg?branch=main
   :target: https://github.com/QSD-Group/QSDsan/actions/workflows/build-coverage.yml

.. Codecov
.. image:: https://codecov.io/gh/QSD-Group/QSDsan/branch/main/graph/badge.svg?token=Z1CASBXEOE
   :target: https://codecov.io/gh/QSD-Group/QSDsan

.. Binder launch of tutorials
.. image:: ./docs/source/images/custom_binder_logo.svg
   :target: https://mybinder.org/v2/gh/QSD-Group/QSDsan-env/main?urlpath=git-pull%3Frepo%3Dhttps%253A%252F%252Fgithub.com%252FQSD-group%252FQSDsan%26urlpath%3Dlab%252Ftree%252FQSDsan%252Fdocs%252Fsource%252Ftutorials%26branch%3Dmain

.. .. Email subscription form
.. .. image:: https://img.shields.io/badge/news-subscribe-F3A93C?style=flat&logo=rss
..    :target: https://groups.webservices.illinois.edu/subscribe/154591

.. YouTube video
.. image:: https://img.shields.io/endpoint?color=%23ff0000&label=YouTube%20 @qsd-group&url=https%3A%2F%2Fyoutube-channel-badge-blond.vercel.app%2Fapi%2Fvideos
   :target: https://www.youtube.com/@qsd-group

|

.. contents::

|

What is ``QSDsan``?
-------------------
``QSDsan`` is an open-source, community-led platform for the quantitative sustainable design (QSD) [1]_ of sanitation and resource recovery systems [2]_. Built in Python, it integrates process modeling, system simulation, techno-economic analysis (TEA), and life cycle assessment (LCA) to support transparent, reproducible, and comprehensive evaluation of emerging technologies. Leveraging BioSTEAM [3]_, the platform provides modular and extensible tools to compare treatment configurations, evaluate resource recovery opportunities, and assess energy, cost, and sustainability tradeoffs. These capabilities support the research, development, and deployment (RD&D) of early-stage technologies for advancing sustainable water management and circular resource recovery.

All systems developed with ``QSDsan`` are included in the package `EXPOsan <https://github.com/QSD-Group/EXPOsan>`_ - exposition of sanitation and resource recovery systems.


Installation
------------
``QSDsan`` requires Python 3.12 or newer. The easiest way to install ``QSDsan`` is through ``pip`` in a command-line interface (e.g., terminal, PowerShell, etc.):

.. code::

    pip install qsdsan

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

To get the git version (use the ``depth`` flag to choose how many commit histories you want to clone):

.. code:: bash

    git clone https://github.com/QSD-Group/QSDsan.git --depth=1

Then navigate into the repository (``cd QSDsan``) and install in editable mode with development dependencies:

.. code:: bash

    pip install -e ".[dev]"


.. note::

   Using the ``depth`` flag will only clone the main branch by default. If you need other branches, add the ``--no-single-branch`` flag:

   .. code:: bash

       git clone https://github.com/<YOUR_USERNAME>/QSDsan.git --depth=1 --no-single-branch


For diagram generation, ``QSDsan`` uses Graphviz. If diagrams fail to render, install Graphviz following the `official Graphviz download instructions <https://graphviz.org/download/>`_ and see the `FAQ <docs/source/FAQ.rst#graphviz-installation>`_ for a quick check.


Documentation
-------------
You can find tutorials and documents at:

   https://qsdsan.readthedocs.io

All tutorials are written using Jupyter Notebook, you can run your own Jupyter environment, or you can click the ``launch binder`` badge on the top to launch the environment in your browser.

For each of these tutorials, we are also recording videos where one of the QSD group members will go through the tutorial step-by-step. We are gradually releasing these videos on our `YouTube channel <https://www.youtube.com/channel/UC8fyVeo9xf10KeuZ_4vC_GA>`_ so subscribe to receive updates!


Authors and Contributing
------------------------
``QSDsan`` and its related packages are developed by the Quantitative Sustainable Design Group and the broader community. `Yalin Li <https://github.com/yalinli2>`_ is the currently maintainer.  See `commit history <https://github.com/QSD-Group/QSDsan/graphs/contributors>`_ for contributors who have contributed to the repository. 

If you want to contribute to ``QSDsan``, please refer to the `Contributing Guidelines <https://qsdsan.readthedocs.io/en/latest/CONTRIBUTING.html>`_ section of the documentation for instructions and guidelines.


License Information
-------------------
Please refer to the ``LICENSE.txt`` for information on the terms & conditions for usage of this software, and a DISCLAIMER OF ALL WARRANTIES.


References
----------
.. [1] Li, Y.; Zhang, X.; Morgan, V.L.; Lohman, H.A.C.; Rowles, L.S.; Mittal, S.; Kogler, A.; Cusick, R.D.; Tarpeh, W.A.; Guest, J.S. QSDsan: An integrated platform for quantitative sustainable design of sanitation and resource recovery systems. Environ. Sci.: Water Res. Technol. 2022, 8 (10), 2289-2303. https://doi.org/10.1039/d2ew00455k.

.. [2] Li, Y.; Trimmer, J.T.; Hand, S.; Zhang, X.; Chambers, K.G.; Lohman, H.A.C.; Shi, R.; Byrne, D.M.; Cook, S.M.; Guest, J.S. Quantitative Sustainable Design (QSD): A Methodology for the Prioritization of Research, Development, and Deployment of Technologies. (Tutorial Review) Environ. Sci.: Water Res. Technol. 2022, 8 (11), 2439–2465. https://doi.org/10.1039/D2EW00431C.

.. [3] Cortés-Peña, Y.; Kumar, D.; Singh, V.; Guest, J.S. BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040.

.. create custom binder badge: https://mybinder.readthedocs.io/en/latest/howto/badges.html